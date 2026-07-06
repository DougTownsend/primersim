"""primersim — Python interface to the C++ PCR equilibrium simulator.

The heavy lifting (thermodynamics + parallel numerics) lives in the
compiled ``primersim._core`` extension. This module owns the *experiment
description*: address tokens, binom(N, 2) enumeration, touchdown temperature
profiles, scalar->vector broadcasting, output shaping, and validation.

Example
-------
    from primersim import Simulator

    sim = Simulator("primers.csv", mv=100, dv=1.5, dntp=0.2, num_cpu=8)
    sim.add_address(1, "2r")          # primer 1 fwd, RC of primer 2 reverse
    sim.add_address(3, 4)
    sim.set_address_conc(3e-15 / 2)   # per-address initial template conc
    sim.set_primer_conc(250e-9)
    sim.set_temperature(55.0)         # constant; pass a list for touchdown
    res = sim.simulate_all(cycles=30)
    print(res.ratios)                 # one amplification ratio per address
    print(res.concentrations)         # N x N final F+R totals
"""

from __future__ import annotations

import itertools
import warnings
from typing import Iterable, List, Optional, Sequence, Tuple, Union

from . import _core

__all__ = [
    "Simulator",
    "AddressResult",
    "BatchResult",
    "parse_token",
]

Token = Union[int, str]


def parse_token(tok: Token) -> Tuple[int, bool]:
    """Parse a primer-slot token into ``(primer_index, is_reverse_complement)``.

    An ``int`` (or a bare digit string) names the primer's literal sequence.
    A ``"<idx>r"`` token (e.g. ``"2r"``) names that primer's reverse
    complement — the same convention as ``pairings.csv`` (``"0r,4"``).
    """
    if isinstance(tok, bool):  # bool is an int subclass; reject to avoid surprises
        raise TypeError(f"invalid primer token: {tok!r}")
    if isinstance(tok, int):
        if tok < 0:
            raise ValueError(f"primer index must be non-negative: {tok}")
        return tok, False
    if isinstance(tok, str):
        s = tok.strip()
        rc = False
        if s and s[-1] in ("r", "R"):
            rc = True
            s = s[:-1]
        if not s or not s.isdigit():
            raise ValueError(f"invalid primer token: {tok!r}")
        return int(s), rc
    raise TypeError(f"primer token must be int or str, got {type(tok).__name__}")


class AddressResult:
    """Result of :meth:`Simulator.simulate_address`.

    Attributes
    ----------
    ratio : float
        target concentration / sum of off-target concentrations.
    concentrations : list[float]
        Final total DNA (F-strand + R-strand) of every address present in
        this run — one float per address.
    """

    __slots__ = ("ratio", "concentrations")

    def __init__(self, ratio: float, concentrations: List[float]):
        self.ratio = ratio
        self.concentrations = concentrations

    def __repr__(self) -> str:
        return f"AddressResult(ratio={self.ratio!r}, concentrations=[{len(self.concentrations)} values])"


class BatchResult:
    """Result of :meth:`Simulator.simulate_all`.

    Attributes
    ----------
    ratios : list[float]
        One amplification ratio per address.
    concentrations : list[list[float]]
        N x N matrix. Row ``a`` is address ``a``'s simulation; entry
        ``[a][i]`` is the final F+R total of address ``i`` in that run.
    """

    __slots__ = ("ratios", "concentrations")

    def __init__(self, ratios: List[float], concentrations: List[List[float]]):
        self.ratios = ratios
        self.concentrations = concentrations

    def __repr__(self) -> str:
        n = len(self.ratios)
        return f"BatchResult(ratios=[{n} values], concentrations={n}x{n})"


class Simulator:
    """Stateful wrapper around the compiled PCR simulator.

    Parameters
    ----------
    primers_file : str
        Path to a primer pool file (one ACGT sequence per line).
    mv, dv, dntp : float
        Monovalent cation, divalent cation, and dNTP concentrations. These
        are baked into the dH/dS cache at construction; changing them later
        via :meth:`set_cations` rebuilds the (expensive) cache.
    num_cpu : int
        Default worker-thread count for :meth:`simulate_all`.
    """

    def __init__(self, primers_file: str, mv: float, dv: float, dntp: float,
                 num_cpu: int):
        self._engine = _core.Engine(primers_file, float(mv), float(dv),
                                    float(dntp), int(num_cpu))
        self._num_cpu = int(num_cpu)
        self._mv, self._dv, self._dntp = float(mv), float(dv), float(dntp)

        # Experiment description, filled via the setters below.
        self._addresses: List[Tuple[int, int, bool, bool]] = []
        self._address_conc: Optional[Union[float, List[float]]] = None
        self._primer_f: Optional[float] = None
        self._primer_r: Optional[float] = None
        self._temperature: Optional[Union[float, List[float]]] = None

    # ---- introspection -------------------------------------------------

    @property
    def num_primers(self) -> int:
        return self._engine.num_primers()

    @property
    def num_addresses(self) -> int:
        return len(self._addresses)

    @property
    def addresses(self) -> List[Tuple[int, int, bool, bool]]:
        """Current addresses as ``(f_idx, r_idx, f_rc, r_rc)`` tuples."""
        return list(self._addresses)

    def address_tokens(self) -> List[Tuple[str, str]]:
        """Addresses rendered back to ``("<idx>[r]", "<idx>[r]")`` tokens.

        Round-trips against the ``pairings.csv`` / ``write_pairings`` format.
        """
        def tok(idx: int, rc: bool) -> str:
            return f"{idx}r" if rc else f"{idx}"
        return [(tok(f, frc), tok(r, rrc)) for (f, r, frc, rrc) in self._addresses]

    # ---- address configuration ----------------------------------------

    def add_address(self, f: Token, r: Token) -> None:
        """Append one address = (forward slot, reverse slot).

        Each slot is a primer index (literal sequence) or a ``"<idx>r"``
        token (that primer's reverse complement). e.g. ``add_address(1, "2r")``.
        """
        f_idx, f_rc = parse_token(f)
        r_idx, r_rc = parse_token(r)
        self._check_index(f_idx)
        self._check_index(r_idx)
        self._addresses.append((f_idx, r_idx, f_rc, r_rc))

    def set_addresses(self, addresses: Iterable[Sequence[Token]]) -> None:
        """Bulk-replace all addresses. Each element is an ``(f, r)`` pair of
        tokens (see :meth:`add_address`)."""
        new: List[Tuple[int, int, bool, bool]] = []
        for pair in addresses:
            if len(pair) != 2:
                raise ValueError(f"each address must be an (f, r) pair, got {pair!r}")
            f_idx, f_rc = parse_token(pair[0])
            r_idx, r_rc = parse_token(pair[1])
            self._check_index(f_idx)
            self._check_index(r_idx)
            new.append((f_idx, r_idx, f_rc, r_rc))
        self._addresses = new

    def enumerate_pairs(self, indices: Optional[Sequence[int]] = None) -> None:
        """Set addresses to every binom(N, 2) unordered pair of primers.

        Enumeration is done here in Python (``itertools.combinations``), not
        in C++. ``indices`` defaults to all primers; both slots use each
        primer's literal (forward) sequence.
        """
        if indices is None:
            indices = range(self._engine.num_primers())
        self.set_addresses(list(itertools.combinations(indices, 2)))

    # ---- concentrations / cations / temperature ------------------------

    def set_address_conc(self, conc: Union[float, Sequence[float]]) -> None:
        """Initial template concentration seeded into each address.

        A scalar applies the same value to every address; a sequence gives a
        per-address value (length must equal the address count at simulate
        time).
        """
        if _is_scalar(conc):
            self._address_conc = float(conc)
        else:
            self._address_conc = [float(x) for x in conc]

    def set_primer_conc(self, x: Optional[float] = None, *,
                        f: Optional[float] = None,
                        r: Optional[float] = None) -> None:
        """Set forward/reverse primer concentrations.

        ``set_primer_conc(x)`` sets both F and R to ``x``.
        ``set_primer_conc(f=..., r=...)`` sets them independently — F and R
        are distinct, independently-depleted species in the engine.
        """
        if x is not None:
            if f is not None or r is not None:
                raise ValueError("pass either a positional scalar or f=/r=, not both")
            self._primer_f = self._primer_r = float(x)
        else:
            if f is None or r is None:
                raise ValueError("set_primer_conc needs either x, or both f= and r=")
            self._primer_f = float(f)
            self._primer_r = float(r)

    def set_cations(self, mv: float, dv: float, dntp: float) -> None:
        """Change cations and **rebuild the dH/dS cache** (expensive)."""
        warnings.warn(
            "set_cations rebuilds the dimer cache (full pairwise thal sweep) — "
            "this is expensive; avoid calling it in a loop.",
            stacklevel=2,
        )
        self._engine.set_cations(float(mv), float(dv), float(dntp))
        self._mv, self._dv, self._dntp = float(mv), float(dv), float(dntp)

    def set_temperature(self, temp: Union[float, Sequence[float]]) -> None:
        """Set the anneal temperature (deg C).

        A scalar is a constant-temperature run; a sequence is a per-cycle
        profile (touchdown PCR) whose length must equal ``cycles`` at
        simulate time.
        """
        if _is_scalar(temp):
            self._temperature = float(temp)
        else:
            self._temperature = [float(t) for t in temp]

    # ---- simulation ----------------------------------------------------

    def simulate_address(self, addr: int, cycles: int) -> AddressResult:
        """Simulate a single address by index; see :class:`AddressResult`."""
        self._require_ready()
        if not (0 <= addr < len(self._addresses)):
            raise IndexError(f"address index {addr} out of range [0, {len(self._addresses)})")
        profile = self._build_profile(cycles)
        conc = self._build_address_conc()
        ratio, concs = self._engine.simulate_address(
            self._addresses, int(addr), profile, conc,
            self._primer_f, self._primer_r)
        return AddressResult(ratio, list(concs))

    def simulate_all(self, cycles: int,
                     threads: Optional[int] = None) -> BatchResult:
        """Simulate every address; see :class:`BatchResult`.

        ``threads`` defaults to the constructor's ``num_cpu``.
        """
        self._require_ready()
        profile = self._build_profile(cycles)
        conc = self._build_address_conc()
        nthreads = self._num_cpu if threads is None else int(threads)
        if nthreads < 1:
            raise ValueError("threads must be >= 1")
        ratios, mat = self._engine.simulate_all(
            self._addresses, profile, conc,
            self._primer_f, self._primer_r, nthreads)
        return BatchResult(list(ratios), [list(row) for row in mat])

    # ---- helpers -------------------------------------------------------

    def _check_index(self, idx: int) -> None:
        n = self._engine.num_primers()
        if not (0 <= idx < n):
            raise IndexError(f"primer index {idx} out of range [0, {n})")

    def _require_ready(self) -> None:
        if not self._addresses:
            raise RuntimeError("no addresses set; call add_address / set_addresses first")
        if self._address_conc is None:
            raise RuntimeError("address concentration not set; call set_address_conc")
        if self._primer_f is None or self._primer_r is None:
            raise RuntimeError("primer concentration not set; call set_primer_conc")
        if self._temperature is None:
            raise RuntimeError("temperature not set; call set_temperature")

    def _build_profile(self, cycles: int) -> List[float]:
        cycles = int(cycles)
        if cycles < 1:
            raise ValueError("cycles must be >= 1")
        temp = self._temperature
        if _is_scalar(temp):
            return [float(temp)] * cycles
        if len(temp) != cycles:
            raise ValueError(
                f"temperature profile has {len(temp)} entries but cycles={cycles}; "
                "a per-cycle profile must match the cycle count")
        return list(temp)

    def _build_address_conc(self) -> List[float]:
        n = len(self._addresses)
        conc = self._address_conc
        if _is_scalar(conc):
            return [float(conc)] * n
        if len(conc) != n:
            raise ValueError(
                f"address_conc has {len(conc)} entries but there are {n} addresses")
        return list(conc)


def _is_scalar(x) -> bool:
    return isinstance(x, (int, float)) and not isinstance(x, bool)
