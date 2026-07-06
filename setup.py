"""Build config for the primersim._core C++ extension.

Pybind11Extension is imperative, so the extension lives here while
metadata lives in pyproject.toml. Compiles the four existing translation
units plus src/pybind.cpp into primersim._core.

Portability note: -mavx / -march=native are opt-in (set PRIMERSIM_NATIVE=1)
so plain `pip install .` produces a portable module; a native build gets
the vectorized codegen the standalone test_eq uses.
"""
import os

from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup

extra_compile_args = ["-O3"]
if os.environ.get("PRIMERSIM_NATIVE"):
    extra_compile_args.append("-march=native")

ext_modules = [
    Pybind11Extension(
        "primersim._core",
        sources=[
            "src/thal.cpp",
            "src/eq.cpp",
            "src/sim.cpp",
            "src/address_eval.cpp",
            "src/pybind.cpp",
        ],
        include_dirs=["include"],
        extra_compile_args=extra_compile_args,
        libraries=["pthread"],
        cxx_std=17,
    ),
]

setup(ext_modules=ext_modules, cmdclass={"build_ext": build_ext})
