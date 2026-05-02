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

    /*
    Concentration array indices (active set after A/B/Y removal):
    0, F:  Forward primer
    1, R:  Reverse primer
    2, X:  Nonspec binding site
    3, FH: Forward primer hairpin
    4, RH: Reverse primer hairpin
    5, FF: Forward primer dimer
    6, RR: Reverse primer dimer
    7, FR: Forward-Reverse primer dimer
    8, FX: Forward bound to nonspec
    9, RX: Reverse bound to nonspec

    Equilibrium constant indices:
    0, K_FH:   Forward hairpin
    1, K_RH:   Reverse hairpin
    2, K_FF:   Forward dimer
    3, K_RR:   Reverse dimer
    4, K_FR:   Forward-Reverse dimer
    5, K_FX:   Forward to nonspec
    6, K_RX:   Reverse to nonspec
    */
    const int F = 0;
    const int R = 1;
    const int X = 2;
    const int FH = 3;
    const int RH = 4;
    const int FF = 5;
    const int RR = 6;
    const int FR = 7;
    const int FX = 8;
    const int RX = 9;

    const int K_FH = 0;
    const int K_RH = 1;
    const int K_FF = 2;
    const int K_RR = 3;
    const int K_FR = 4;
    const int K_FX = 5;
    const int K_RX = 6;

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
            Real last_val[2];
            Real tmp[4];
            Real c[19];
            Real c0[6];
            Real k[13];
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
            void calc_strand_bindings(EQ &eq, const std::vector<double> &temp_c_profile, int i, int cycle, int addr, int end5, int end3);
            double dhds_to_eq_const(double dhds[2], double temp_c);
            void eval_addresses_thread(const char *out_filename, unsigned int pcr_cycles, const std::vector<double> &temp_c_profile, double dna_conc, double primer_f_conc, double primer_r_conc, double mv_conc, double dv_conc, double dntp_conc);
            void eval_addresses(const char *out_filename, unsigned int pcr_cycles, const std::vector<double> &temp_c_profile, double dna_conc, double primer_f_conc, double primer_r_conc, double mv_conc, double dv_conc, double dntp_conc);
            void shuffle_addresses(void);
    };
}
#endif
