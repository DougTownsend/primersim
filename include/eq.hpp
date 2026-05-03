#include <vector>
#include <string>
#include "thal.h"
#include <mutex>
#include <cstdlib>
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
            double tmp_dhds[2][4][2];

            double dhds_primer_f_addr_f[2];
            double dhds_primer_f_addr_frc[2];
            double dhds_primer_f_addr_r[2];
            double dhds_primer_f_addr_rrc[2];
            double dhds_primer_r_addr_f[2];
            double dhds_primer_r_addr_frc[2];
            double dhds_primer_r_addr_r[2];
            double dhds_primer_r_addr_rrc[2];

            // Per-cycle cache of dhds_to_eq_const(tmp_dhds[..], temp_c).
            // Populated once per (address, cycle) in sim_pcr; read by
            // calc_strand_bindings (25 reads per address per cycle).
            Real tmp_dhds_K[2][4];

            // Per-cycle cache for the four dhds_primer_*_addr_{frc,rrc}
            // fields read by update_strand_concs.
            Real dhds_K_primer_f_addr_frc;
            Real dhds_K_primer_r_addr_frc;
            Real dhds_K_primer_f_addr_rrc;
            Real dhds_K_primer_r_addr_rrc;
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
            std::vector<Primeanneal::primer_info> primers;
            std::vector<Primeanneal::address> addresses;
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
            void evaluate_addresses(const char *in_filename, const char *out_filename, double dna_conc, double primer_conc, double mv_conc, double dv_conc, double dntp_conc);
            void eval_thread(const char *out_filename, double dna_conc, double primer_conc, double mv_conc, double dv_conc, double dntp_conc);
            double sim_pcr(const char *out_filename, unsigned int addr, unsigned int pcr_cycles, const std::vector<double> &temp_c_profile, double dna_conc, double primer_f_conc, double primer_r_conc, double mv_conc, double dv_conc, double dntp_conc);
            void update_strand_concs(EQ &eq, int i, int addr, const std::vector<double> &temp_c_profile, int cycle, int end5, int end3);
            void calc_strand_bindings(EQ &eq, const std::vector<double> &temp_c_profile, int i, int cycle, int addr, int end5, int end3,
                                       Real &nonspec_total, Real &sum_f_weighted, Real &sum_r_weighted);
            double dhds_to_eq_const(double dhds[2], double temp_c);
            void eval_addresses_thread(const char *out_filename, unsigned int pcr_cycles, const std::vector<double> &temp_c_profile, double dna_conc, double primer_f_conc, double primer_r_conc, double mv_conc, double dv_conc, double dntp_conc);
            void eval_addresses(const char *out_filename, unsigned int pcr_cycles, const std::vector<double> &temp_c_profile, double dna_conc, double primer_f_conc, double primer_r_conc, double mv_conc, double dv_conc, double dntp_conc);
            void shuffle_addresses(void);
    };
}
#endif
