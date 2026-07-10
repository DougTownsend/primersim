// pybind11 binding for primersim. Exposes a thin `Engine` class over
// Primeanneal: it owns the primer pool + the (primers, cations)-keyed
// dimer cache, and runs the parallel numerics. Everything the biologist
// touches — address tokens, binom enumeration, touchdown profiles,
// scalar→vector broadcasting, output shaping, validation — lives in the
// Python `Simulator` wrapper (primersim/__init__.py).
//
// The compute methods release the GIL (py::gil_scoped_release): the C++
// thread pool touches no Python objects, so this gives real parallelism
// without the multi-MB cache duplication that multiprocessing would cost.

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <tuple>
#include <vector>
#include <string>
#include <stdexcept>

#include "eq.hpp"

namespace py = pybind11;
using primersim::Primeanneal;
using primersim::pairing;

// One address from Python, as a (f_idx, r_idx, f_rc, r_rc) tuple. Python
// has already resolved "<idx>r" tokens into indices + rc flags.
using AddrTuple = std::tuple<int, int, bool, bool>;

namespace {

std::vector<pairing> to_pairings(const std::vector<AddrTuple> &addrs) {
    std::vector<pairing> out;
    out.reserve(addrs.size());
    for (const auto &a : addrs) {
        pairing p;
        p.f_idx = std::get<0>(a);
        p.r_idx = std::get<1>(a);
        p.f_rc  = std::get<2>(a);
        p.r_rc  = std::get<3>(a);
        out.push_back(p);
    }
    return out;
}

// Engine: a Primeanneal plus the salt the cache was last built with.
// Storing the salt here (rather than accepting it per-call) keeps the
// dimer cache and the live hairpin recomputation in sim_pcr consistent —
// the one correctness constraint the design calls out.
class Engine {
public:
    Engine(const std::string &primers_file, double mv, double dv,
           double dntp, unsigned int num_cpu)
        : mv_(mv), dv_(dv), dntp_(dntp) {
        pa_.num_cpu = num_cpu > 0 ? num_cpu : 1;
        pa_.read_primer_pool(primers_file.c_str());
        if (pa_.primer_pool.empty())
            throw std::runtime_error("no primers read from '" + primers_file + "'");
        rebuild_cache();
    }

    size_t num_primers() const { return pa_.primer_pool.size(); }

    std::string primer_seq(size_t i) const {
        if (i >= pa_.primer_pool.size())
            throw std::out_of_range("primer index out of range");
        return std::string(pa_.primer_pool[i].seq);
    }

    unsigned int num_cpu() const { return pa_.num_cpu; }

    // Constant dH (cal/mol) added to the 3'-anchored (productive) bindings
    // only — models extra stability from polymerase binding the primer's
    // 3' junction. Applied per-simulate at K computation; no cache rebuild,
    // so this is cheap to sweep. Negative = more stable.
    void set_3p_dh_offset(double dh) { pa_.dh_3p_offset = dh; }
    double get_3p_dh_offset() const { return pa_.dh_3p_offset; }

    // --- Manual thermodynamic-cache control -------------------------------
    // Override the dH/dS the engine uses, so primersim can be driven by an
    // external predictor (ViennaRNA / oxDNA / ...) instead of thal. The cache
    // is keyed by SLOT = 2*primer_index + kind, where kind 0 = the primer's
    // literal sequence and kind 1 = its reverse complement; num_slots() =
    // 2*num_primers(). Two caches:
    //   total (thal_any)  — symmetric interaction        [get/set_cache_total]
    //   3'  (thal_end1)   — slot `a`'s 3' end anchored    [get/set_cache_3p]
    // dh in cal/mol, ds in cal/mol/K; the engine forms K(T) = exp(-(dh - T*ds)
    // / (R*T)) each cycle. Construction and set_cations() repopulate from thal,
    // so apply manual values AFTER those. Entries left unset keep their thal
    // values (partial override is safe); note ds must not accidentally be left
    // at 0 with dh=0 for a real interaction (that is a phantom K=1 duplex).
    size_t num_slots() const { return 2 * pa_.primer_pool.size(); }

    void set_cache_total(int a, int b, double dh, double ds) {
        auto &e = pa_.dimer_cache[total_idx(a, b)];
        e.dh = (ThalReal)dh; e.ds = (ThalReal)ds;
    }
    std::tuple<double, double> get_cache_total(int a, int b) const {
        const auto &e = pa_.dimer_cache[total_idx(a, b)];
        return std::make_tuple((double)e.dh, (double)e.ds);
    }
    void set_cache_3p(int a, int b, double dh, double ds) {
        auto &e = pa_.dimer_3p_cache[tp_idx(a, b)];
        e.dh = (ThalReal)dh; e.ds = (ThalReal)ds;
    }
    std::tuple<double, double> get_cache_3p(int a, int b) const {
        const auto &e = pa_.dimer_3p_cache[tp_idx(a, b)];
        return std::make_tuple((double)e.dh, (double)e.ds);
    }

    // Rebuild the dH/dS cache under new cations. Expensive (runs the full
    // pairwise thal sweep over num_cpu threads); the salt correction is
    // baked into the stored entropy, so this is mandatory when cations
    // change. Temperature does NOT go here — it is applied per-cycle.
    void set_cations(double mv, double dv, double dntp) {
        mv_ = mv; dv_ = dv; dntp_ = dntp;
        py::gil_scoped_release nogil;
        rebuild_cache();
    }

    // Single-address simulation (item 6). Returns (ratio, final_concs)
    // where final_concs[i] is address i's final F+R total in this run.
    std::tuple<double, std::vector<double>>
    simulate_address(const std::vector<AddrTuple> &addrs, unsigned int addr,
                     const std::vector<double> &temp_profile,
                     const std::vector<double> &address_conc,
                     double primer_f, double primer_r) {
        std::vector<pairing> pairings = to_pairings(addrs);
        if (addr >= pairings.size())
            throw std::out_of_range("address index out of range");
        if (address_conc.size() != pairings.size())
            throw std::invalid_argument("address_conc length must equal number of addresses");
        double ratio;
        std::vector<double> concs;
        {
            py::gil_scoped_release nogil;
            ratio = pa_.sim_pcr(pairings, nullptr, addr,
                                (unsigned int)temp_profile.size(), temp_profile,
                                address_conc, primer_f, primer_r,
                                mv_, dv_, dntp_, /*log_cycles=*/false, &concs);
        }
        return std::make_tuple(ratio, std::move(concs));
    }

    // Batch simulation (item 7). Returns (ratios, NxN conc matrix).
    std::tuple<std::vector<double>, std::vector<std::vector<double>>>
    simulate_all(const std::vector<AddrTuple> &addrs,
                 const std::vector<double> &temp_profile,
                 const std::vector<double> &address_conc,
                 double primer_f, double primer_r, unsigned int threads) {
        std::vector<pairing> pairings = to_pairings(addrs);
        if (address_conc.size() != pairings.size())
            throw std::invalid_argument("address_conc length must equal number of addresses");
        std::vector<double> ratios;
        std::vector<std::vector<double>> mat;
        {
            py::gil_scoped_release nogil;
            pa_.sim_all(pairings, (unsigned int)temp_profile.size(), temp_profile,
                        address_conc, primer_f, primer_r, mv_, dv_, dntp_,
                        threads, &ratios, &mat);
        }
        return std::make_tuple(std::move(ratios), std::move(mat));
    }

private:
    void rebuild_cache() {
        pa_.populate_dimer_cache(mv_, dv_, dntp_, /*temp_c=*/55.0);
    }

    void check_slot(int s) const {
        if (s < 0 || (size_t)s >= 2 * pa_.primer_pool.size())
            throw std::out_of_range("cache slot out of range [0, 2*num_primers)");
    }
    // Canonical (symmetric) index into the triangular total cache.
    size_t total_idx(int a, int b) const {
        check_slot(a); check_slot(b);
        if (a > b) { int t = a; a = b; b = t; }
        size_t M = 2 * pa_.primer_pool.size();
        return (size_t)a * (2 * M - a - 1) / 2 + (size_t)b;
    }
    // Row-major index into the square 3'-anchored cache (a = anchored primer).
    size_t tp_idx(int a, int b) const {
        check_slot(a); check_slot(b);
        size_t M = 2 * pa_.primer_pool.size();
        return (size_t)a * M + (size_t)b;
    }

    Primeanneal pa_;
    double mv_, dv_, dntp_;
};

}  // namespace

PYBIND11_MODULE(_core, m) {
    m.doc() = "C++ core for primersim (PCR equilibrium simulator).";

    py::class_<Engine>(m, "Engine")
        .def(py::init<const std::string &, double, double, double, unsigned int>(),
             py::arg("primers_file"), py::arg("mv"), py::arg("dv"),
             py::arg("dntp"), py::arg("num_cpu"))
        .def("num_primers", &Engine::num_primers)
        .def("primer_seq", &Engine::primer_seq, py::arg("i"))
        .def("num_cpu", &Engine::num_cpu)
        .def("set_cations", &Engine::set_cations,
             py::arg("mv"), py::arg("dv"), py::arg("dntp"))
        .def("set_3p_dh_offset", &Engine::set_3p_dh_offset, py::arg("dh"))
        .def("get_3p_dh_offset", &Engine::get_3p_dh_offset)
        .def("num_slots", &Engine::num_slots)
        .def("set_cache_total", &Engine::set_cache_total,
             py::arg("a"), py::arg("b"), py::arg("dh"), py::arg("ds"))
        .def("get_cache_total", &Engine::get_cache_total, py::arg("a"), py::arg("b"))
        .def("set_cache_3p", &Engine::set_cache_3p,
             py::arg("a"), py::arg("b"), py::arg("dh"), py::arg("ds"))
        .def("get_cache_3p", &Engine::get_cache_3p, py::arg("a"), py::arg("b"))
        .def("simulate_address", &Engine::simulate_address,
             py::arg("addresses"), py::arg("addr"), py::arg("temp_profile"),
             py::arg("address_conc"), py::arg("primer_f"), py::arg("primer_r"))
        .def("simulate_all", &Engine::simulate_all,
             py::arg("addresses"), py::arg("temp_profile"),
             py::arg("address_conc"), py::arg("primer_f"),
             py::arg("primer_r"), py::arg("threads"));
}
