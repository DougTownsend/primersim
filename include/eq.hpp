#include <vector>
#include <string>
#include "thal.h"
#include <mutex>
#include <cstdlib>
#include <cstdint>
#include <cmath>
#include <limits>

#ifndef EQ_HPP
#define EQ_HPP

namespace primersim{
    // Switchable underlying float type. Default: 64-bit IEEE double
    // (SSE2-friendly on Linux x86-64). Pass -DPSIM_USE_LONG_DOUBLE to
    // promote to 80-bit x87 long double — the 19-digit margin is below
    // the 9-digit print precision so it's not normally needed.
#ifdef PSIM_USE_LONG_DOUBLE
    typedef long double Real;
#else
    typedef double Real;
#endif

    // Concentration species (used to index EQ::c and EQ::c0):
    //   F:  Forward primer       FF: Forward primer dimer
    //   R:  Reverse primer       RR: Reverse primer dimer
    //   X:  Nonspec binding site FR: Forward-Reverse primer dimer
    //   FH: Forward primer hairpin
    //   RH: Reverse primer hairpin
    //   FX: Forward bound to nonspec
    //   RX: Reverse bound to nonspec
    enum class Species : int {
        F=0, R=1, X=2, FH=3, RH=4, FF=5, RR=6, FR=7, FX=8, RX=9
    };

    // Equilibrium constants (used to index EQ::k).
    enum class Rate : int {
        K_FH=0, K_RH=1, K_FF=2, K_RR=3, K_FR=4, K_FX=5, K_RX=6
    };

    // Aliases so the existing `c[F]`, `k[K_FH]` syntax works without
    // qualification. Indexing c with a Rate (e.g. `c[K_FF]`) is a
    // compile error — IndexedArray below only matches one Idx type.
    inline constexpr Species F  = Species::F;
    inline constexpr Species R  = Species::R;
    inline constexpr Species X  = Species::X;
    inline constexpr Species FH = Species::FH;
    inline constexpr Species RH = Species::RH;
    inline constexpr Species FF = Species::FF;
    inline constexpr Species RR = Species::RR;
    inline constexpr Species FR = Species::FR;
    inline constexpr Species FX = Species::FX;
    inline constexpr Species RX = Species::RX;

    inline constexpr Rate K_FH = Rate::K_FH;
    inline constexpr Rate K_RH = Rate::K_RH;
    inline constexpr Rate K_FF = Rate::K_FF;
    inline constexpr Rate K_RR = Rate::K_RR;
    inline constexpr Rate K_FR = Rate::K_FR;
    inline constexpr Rate K_FX = Rate::K_FX;
    inline constexpr Rate K_RX = Rate::K_RX;

    // Real array typed by an index enum. operator[] only accepts the
    // matching Idx type, so c (Species-indexed) and k (Rate-indexed)
    // can't be confused. Use .at(int) for raw integer access (only
    // print_state needs this for its dump loops).
    template<typename Idx, std::size_t N>
    struct IndexedArray {
        Real data[N];
        Real& operator[](Idx i)             { return data[static_cast<int>(i)]; }
        const Real& operator[](Idx i) const { return data[static_cast<int>(i)]; }
        Real& at(int i)             { return data[i]; }
        const Real& at(int i) const { return data[i]; }
        static constexpr std::size_t size() { return N; }
    };

    const int addr_end = 0;
    const int primerf = 1;
    const int primerfrc = 2;
    const int primerr = 3;
    const int primerrrc = 4;

    const int dh = 0;
    const int ds = 1;

    struct dh_ds{
        double dh;
        double ds;
    };

    // A single primer in the flat primer pool. seq is the literal
    // sequence (read from primers.csv); rc is its reverse complement.
    // Both are owned char* allocations, freed by ~Primeanneal.
    struct primer_seq {
        char *seq;
        char *rc;
    };

    // An (F-primer, R-primer) tuple naming one PCR pair. *_idx points
    // into Primeanneal::primer_pool; *_rc=true means use that primer's
    // reverse complement (primer_pool[idx].rc) instead of its canonical
    // sequence (primer_pool[idx].seq). On disk this is encoded as
    // "<idx>[r],<idx>[r]" — e.g. "0r,4" pairs the RC of primer 0 with
    // the seq of primer 4. A "pairing" — the full assignment of F/R
    // partners across an experiment — is the vector
    // Primeanneal::pairings.
    struct pairing {
        int  f_idx;
        int  r_idx;
        bool f_rc;
        bool r_rc;
    };

    // Cached calc_dimer result (one entry per (addr1, kind1, addr2, kind2)
    // tuple). Width depends on ThalReal — float halves it.
    struct dh_ds_cache_entry {
        ThalReal dh;
        ThalReal ds;
    };

    class address_k_conc{
        public:
            //fstrand[0][0] = 5' F -> F Strand Data -> RRC 3'
            //rstrand[0][0] = 5' R -> R strand data -> FRC 3'
            Real fstrand[5][5];
            Real rstrand[5][5];

            Real fstrand_change[5][5];
            Real rstrand_change[5][5];

            Real total_f_conc;
            Real total_r_conc;
            Real last_f_conc;
            Real last_r_conc;

            //first index:  0=primer_f, 1=primer_r
            //second index: 0=addr_f, 1=addr_frc, 2=addr_r, 3=addr_rrc
            // tmp_dhds: total interaction (thal_any). tmp_dhds_3p: only
            // alignments where the primer's 3' end is at the duplex end
            // (thal_end1 with primer as seq1) — the productive /
            // amplifying portion. dh/ds values come from thal_results
            // (double) but are stored as ThalReal here so they can
            // shrink to 4 bytes/value in -DTHAL_USE_FLOAT builds.
            ThalReal tmp_dhds[2][4][2];
            ThalReal tmp_dhds_3p[2][4][2];

            // Specific (primer × address-region) interactions used by
            // update_strand_concs. _3p variants are the
            // 3'-end-anchored (productive) portion only.
            ThalReal dhds_primer_f_addr_f[2];
            ThalReal dhds_primer_f_addr_frc[2];
            ThalReal dhds_primer_f_addr_r[2];
            ThalReal dhds_primer_f_addr_rrc[2];
            ThalReal dhds_primer_r_addr_f[2];
            ThalReal dhds_primer_r_addr_frc[2];
            ThalReal dhds_primer_r_addr_r[2];
            ThalReal dhds_primer_r_addr_rrc[2];
            ThalReal dhds_primer_f_addr_frc_3p[2];
            ThalReal dhds_primer_r_addr_frc_3p[2];
            ThalReal dhds_primer_f_addr_rrc_3p[2];
            ThalReal dhds_primer_r_addr_rrc_3p[2];

            // Per-cycle K caches (populated once per address per cycle
            // in sim_pcr).
            //   tmp_dhds_K        = K(thal_any)             — total
            //   tmp_dhds_K_3p     = K(thal_end1)            — amplifying
            //   tmp_dhds_K_nonamp = K(any − end1, dh/ds sub) — non-amp;
            //                       set to 0 when total == 3p (i.e. all
            //                       binding is already 3'-end-anchored)
            //                       so we don't model a phantom 0-ΔG
            //                       (K=1) interaction.
            Real tmp_dhds_K[2][4];
            Real tmp_dhds_K_3p[2][4];
            Real tmp_dhds_K_nonamp[2][4];

            // Per-cycle cache for the dhds_primer_*_addr_{frc,rrc}
            // fields read by update_strand_concs. Only the _3p (amp)
            // version drives extension; the original total K is no
            // longer needed in update_strand_concs.
            Real dhds_K_primer_f_addr_frc_3p;
            Real dhds_K_primer_r_addr_frc_3p;
            Real dhds_K_primer_f_addr_rrc_3p;
            Real dhds_K_primer_r_addr_rrc_3p;
    };

    class EQ{
        public:
            // True after the first solve_eq has populated c[F]/c[R].
            // Subsequent calls warm-start from the prior solution instead
            // of resetting to c0/2.
            bool warm_start_valid = false;
            // last_val is only ever last_val[F] / last_val[R]; size 2 is
            // sufficient since F=0, R=1.
            IndexedArray<Species, 2>  last_val;
            IndexedArray<Species, 10> c;
            IndexedArray<Species, 3>  c0;  // only F, R, X are used
            IndexedArray<Rate,    7>  k;
            Real spec_total;
            Real nonspec_total;
            Real spec_fwd_amp;
            Real spec_rev_amp;
            Real nonspec_exp_amp;
            Real nonspec_lin_amp;
            Real best_spec_exp_amp;
            Real best_spec_lin_amp;
            Real best_nonspec_exp_amp;
            Real best_nonspec_lin_amp;
            Real last_nonspec_frc_total;
            Real last_nonspec_rrc_total;
            std::vector<address_k_conc> address_k_conc_vec;

            void calc_cx();
            void calc_bound_concs();
            void solve_eq();
            void print_state(std::string out_filename, std::string s);
    };

    class Primeanneal{
        public:
            struct primer_info{
                char *f;
                char *rc;
                double spec_dg;
                double nonspec_dg_f;
                double nonspec_dg_rc;
                double nonspec_dg_f_f;
                double nonspec_dg_f_rc;
                double nonspec_dg_rc_f;
                double nonspec_dg_rc_rc;
            };
            struct address{
                char *f;
                char *r;
                char *f_rc;
                char *r_rc;
                double temp_c;

                char *get_seq(int i){
                    switch(i){
                        case 0: return f;
                        case 1: return f_rc;
                        case 2: return r;
                        case 3: return r_rc;
                        default: return nullptr;
                    }
                }
            };
            std::vector<Primeanneal::primer_info> primers;     // legacy, used by address_eval.cpp
            std::vector<Primeanneal::address> addresses;       // legacy, used by address_eval.cpp

            // Flat pool of primer sequences (read from primers.csv) and the
            // current F/R pairing assignment (read from pairings.csv). The
            // pairing is what sim_pcr operates on; the primer_pool is what
            // the dimer cache is keyed by.
            std::vector<primer_seq> primer_pool;
            std::vector<pairing> pairings;

            // Symmetric / triangular cache of calc_dimer (thal_any)
            // results, keyed by primer (not pairing) so it survives
            // pairing reshuffles. M = 2*primer_pool.size() sequence
            // slots; each primer contributes seq=kind 0 and rc=kind 1.
            // Canonical pair (a, b) with a <= b at:
            //   idx(a, b) = a*(2*M - a - 1)/2 + b
            // Populated by populate_dimer_cache after read_primer_pool.
            std::vector<dh_ds_cache_entry> dimer_cache;
            // Square cache of thal_end1 results — primer-3'-end-anchored
            // (productive) alignment energies. seq1 is index `a`, seq2
            // is index `b`; thal_end1 anchors seq1's 3' end so the
            // result is asymmetric → no canonical-order swap. Layout:
            //   idx(a, b) = a*M + b
            // Same M as dimer_cache. Populated alongside dimer_cache.
            std::vector<dh_ds_cache_entry> dimer_3p_cache;
            unsigned int address_index;
            std::mutex addr_mtx;
            std::mutex outfile_mtx;
            unsigned int num_cpu;

        public:
            Primeanneal();
            ~Primeanneal();
            void reverse_comp(const char *fwd, char *rc);
            double dg_to_eq_const(double dg, double temp_c);
            double eq_const_to_dg(double eq_const, double temp_c);
            void dg_to_eq_const_mpfr(Real &ret, double dg, double temp_c);
            double eq_const_to_dg_mpfr(Real &tmp, Real &eq_const, double temp_c);
            thal_results calc_dimer(char *seq1, char *seq2, thal_args a);
            thal_results calc_hairpin(char *seq1, thal_args a);
            void read_primers_individual(std::string filename);
            void swap_primer_rc(Primeanneal::primer_info *p);
            bool primer_compare(const Primeanneal::primer_info& a, const Primeanneal::primer_info& b);
            void assign_addresses(const char *out_filename, double temp_c, double mv_conc, double dv_conc, double dntp_conc);
            void assign_addresses_nosort(const char *out_filename);
            void read_addresses(const char *filename, bool with_temp_c);
            // Read a flat primer-per-line CSV into primer_pool. Frees any
            // previously-loaded primers. Each entry gets its reverse
            // complement computed eagerly.
            void read_primer_pool(const char *filename);
            // Read an "f_idx,r_idx" CSV into pairings (one F/R partnership
            // per row, indices reference primer_pool). Replaces any
            // previous pairings; primer_pool is untouched, so the dimer
            // cache survives this call.
            void read_pairings(const char *filename);
            // Dump the current pairings to an "f_idx,r_idx" CSV. Format
            // round-trips through read_pairings. Used by pairing sweeps
            // to save the trial's assignment alongside its histogram.
            void write_pairings(const char *filename);
            // Precompute every canonical thal pair (primer, kind, primer,
            // kind) under a fixed buffer condition. Uses num_cpu threads.
            // Call after read_primer_pool. Output is pairing-agnostic —
            // re-shuffling pairings does NOT invalidate this cache.
            void populate_dimer_cache(double mv_conc, double dv_conc, double dntp_conc, double temp_c);
            void evaluate_addresses(const char *in_filename, const char *out_filename, double dna_conc, double primer_conc, double mv_conc, double dv_conc, double dntp_conc);
            void eval_thread(const char *out_filename, double dna_conc, double primer_conc, double mv_conc, double dv_conc, double dntp_conc);
            // pairings_in names the F/R partner indices for this
            // simulation; primer_pool and dimer_cache are shared
            // (read-only) state. addr indexes into pairings_in. Taking
            // pairings as a parameter — rather than reading
            // this->pairings — lets multiple trials run sim_pcr in
            // parallel against different pairing assignments.
            //
            // log_cycles=true writes the 12-column per-cycle CSV trace
            // to out_filename (append mode). When false (e.g. pairing
            // sweeps that only need the final ratio), out_filename is
            // ignored and no file I/O happens — saves ~12 fopen/fclose
            // cycles per call.
            double sim_pcr(const std::vector<pairing> &pairings_in, const char *out_filename, unsigned int addr, unsigned int pcr_cycles, const std::vector<double> &temp_c_profile, double dna_conc, double primer_f_conc, double primer_r_conc, double mv_conc, double dv_conc, double dntp_conc, bool log_cycles);
            void update_strand_concs(EQ &eq, int i, int addr, const std::vector<double> &temp_c_profile, int cycle, int end5, int end3);
            void calc_strand_bindings(EQ &eq, const std::vector<double> &temp_c_profile, int i, int cycle, int addr, int end5, int end3,
                                       Real &nonspec_total, Real &sum_f_weighted, Real &sum_r_weighted);
            double dhds_to_eq_const(ThalReal dhds[2], double temp_c);
            // Worker for eval_addresses: pulls pair indices off
            // address_index and stores each sim_pcr ratio into
            // (*ratios)[i]. No file I/O happens here — eval_addresses
            // writes results after threads join.
            void eval_addresses_thread(std::vector<double> *ratios, unsigned int pcr_cycles, const std::vector<double> &temp_c_profile, double dna_conc, double primer_f_conc, double primer_r_conc, double mv_conc, double dv_conc, double dntp_conc);
            // Run sim_pcr for every pairing, multithreaded. Writes a
            // pair_idx,ratio CSV to out_filename (NULL skips). If
            // hist_filename is non-NULL, also writes a log10-spaced
            // histogram of the ratios — one row per bin.
            void eval_addresses(const char *out_filename, const char *hist_filename, unsigned int pcr_cycles, const std::vector<double> &temp_c_profile, double dna_conc, double primer_f_conc, double primer_r_conc, double mv_conc, double dv_conc, double dntp_conc);
            // Multithreaded sweep over n_trials random pairing
            // assignments. Each trial draws a fresh full-pool random
            // pairing (no F/R role constraint, with random RC of each
            // primer), then runs sim_pcr per pair × per temperature
            // (constant temp for all pcr_cycles).
            //
            // Output is a sequence of rolling CSV files in out_dir:
            //   sweep.NNNN.csv  (NNNN starts at 0000, rolled at ~500 MB)
            // Each row is one (trial, pair) result with columns:
            //   trial_id, f_seq, r_seq, <T1>C, <T2>C, ...
            // A single mutex serializes the writes so all 50 rows of a
            // trial land contiguously in one file (never split across
            // file boundaries). Concurrent trials may still be
            // interleaved by trial_id within one file — group by
            // trial_id during analysis.
            //
            // Termination:
            //   • n_trials > 0: stop after that many trials.
            //   • n_trials = 0: loop indefinitely until external
            //     SIGTERM (e.g. LSF -W cap).
            //   • max_data_gb > 0: also stop once total bytes written
            //     exceeds the budget. 0 = unlimited size.
            //
            // Each worker draws its own seed from std::random_device
            // (/dev/urandom on Linux), so back-to-back processes and
            // concurrent workers always get independent sequences.
            void sweep_pairings(const char *out_dir, int n_trials, double max_data_gb, const std::vector<int> &temperatures_c, unsigned int pcr_cycles, double dna_conc, double primer_f_conc, double primer_r_conc, double mv_conc, double dv_conc, double dntp_conc);
            void shuffle_addresses(void);
    };
}
#endif
