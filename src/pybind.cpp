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
        .def("simulate_address", &Engine::simulate_address,
             py::arg("addresses"), py::arg("addr"), py::arg("temp_profile"),
             py::arg("address_conc"), py::arg("primer_f"), py::arg("primer_r"))
        .def("simulate_all", &Engine::simulate_all,
             py::arg("addresses"), py::arg("temp_profile"),
             py::arg("address_conc"), py::arg("primer_f"),
             py::arg("primer_r"), py::arg("threads"));
}
