#include <vector>
#include <string>
#include "thal.h"
#include <mutex>
#include <cstdlib>

#ifndef EQ_HPP
#define EQ_HPP

#ifndef FLOAT_PREC
#define FLOAT_PREC (1024)
#endif

namespace primersim{
    /*
    Initial concentrations only include indecies 0-5 inclusive
    Final concentrations include all indecies
    Concentration array indecies:
    0, F:  Forward primer
    1, R:  Reverse primer
    2, A:  Forward binding site
    3, B:  Reverse binding site
    4, X:  Nonspec binding site 1
    5, Y:  Nonspec bidning site 2
    6, FH: Forward primer hairpin
    7, RH: Reverse primer hairpin
    8, FF: Forward primer dimer
    9, RR: Reverse primer dimer
    10, FR: Forward-Reverse primer dimer
    11, FA: Forward specific binding
    12, FB: Forward bound to reverse binding site
    13, RB: Reverse specific binding
    14, RA: Reverse bound to forward binding site
    15, FX: Forward bound to nonspec 1
    16, FY: Forward bound to nonspec 2
    17, RX: Reverse bound to nonspec 1
    18, RY: Reverse bound to nonspec 2

    Equilibrium constant array indecies:
    0, K_FH:   Forward hairpin
    1, K_RH:   Reverse hairpin
    2, K_FF:   Forward dimer
    3, K_RR:   Reverse dimer
    4, K_FR:   Forward-Reverse dimer
    5, K_FA:   Forward specific
    6, K_FB:   Forward to reverse binding site
    7, K_RB:   Reverse specific
    8, K_RA:   Reverse to forward binding site
    9, K_FX:   Forward to nonspec 1
    10, K_FY:   Forward to nonspec 2
    11, K_RX:   Reverse to nonspec 1
    12, K_RY:   Reverse to nonspec 2
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
    /*
    #define F   0
    #define R   1
    #define A   2
    #define B   3
    #define X   4
    #define Y   5
    #define FH  6
    #define RH  7
    #define FF  8
    #define RR  9
    #define FR  10
    #define FA  11
    #define FB  12
    #define RB  13
    #define RA  14
    #define FX  15
    #define FY  16
    #define RX  17
    #define RY  18
    */
    const int K_FH = 0;
    const int K_RH = 1;
    const int K_FF = 2;
    const int K_RR = 3;
    const int K_FR = 4;
    const int K_FX = 5;
    const int K_RX = 6;
    /*
    #define K_FH    0
    #define K_RH    1
    #define K_FF    2
    #define K_RR    3
    #define K_FR    4
    #define K_FA    5
    #define K_FB    6
    #define K_RB    7
    #define K_RA    8
    #define K_FX    9
    #define K_FY    10
    #define K_RX    11
    #define K_RY    12
    */

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

    class Psim_f{
        private:
        class Temps{
            public:
            std::vector<Psim_f> vec;
            uint32_t idx;

            Temps(){
                vec.resize(100); //Should be more than enough
                idx = 0;
            }

            Psim_f& get(){
                idx++;
                idx %= vec.size();
                return vec[idx];
            }
        };

        thread_local inline static Temps tmps;
        public:
        
        mpfr_t val;
        Psim_f(){
            mpfr_init2(val, FLOAT_PREC);
        }
        ~Psim_f(){
            mpfr_clear(val);
        }
        void mul(const Psim_f &o1, const Psim_f &o2){
            mpfr_mul(val, o1.val, o2.val, MPFR_RNDN);
        }
        Psim_f& operator*(const Psim_f& o){
            Psim_f &ret = tmps.get();
            mpfr_mul(ret.val, val, o.val, MPFR_RNDN);
            return ret;
        }
        void mul_d(const Psim_f &o1, double o2){
            mpfr_mul_d(val, o1.val, o2, MPFR_RNDN);
        }
        Psim_f& operator*(const double& o){
            Psim_f &ret = tmps.get();
            mpfr_mul_d(ret.val, val, o, MPFR_RNDN);
            return ret;
        }
        void add(const Psim_f &o1, const Psim_f &o2){
            mpfr_add(val, o1.val, o2.val, MPFR_RNDN);
        }
        Psim_f& operator+(const Psim_f& o){
            Psim_f &ret = tmps.get();
            mpfr_add(ret.val, val, o.val, MPFR_RNDN);
            return ret;
        }
        void add_d(const Psim_f &o1, double o2){
            mpfr_add_d(val, o1.val, o2, MPFR_RNDN);
        }
        Psim_f& operator+(const double& o){
            Psim_f &ret = tmps.get();
            mpfr_add_d(ret.val, val, o, MPFR_RNDN);
            return ret;
        }
        void sub(const Psim_f &o1, const Psim_f &o2){
            mpfr_sub(val, o1.val, o2.val, MPFR_RNDN);
        }
        Psim_f& operator-(const Psim_f& o){
            Psim_f &ret = tmps.get();
            mpfr_sub(ret.val, val, o.val, MPFR_RNDN);
            return ret;
        }
        void sub_d(const Psim_f &o1, double o2){
            mpfr_sub_d(val, o1.val, o2, MPFR_RNDN);
        }
        Psim_f& operator-(const double& o){
            Psim_f &ret = tmps.get();
            mpfr_sub_d(ret.val, val, o, MPFR_RNDN);
            return ret;
        }
        void div(const Psim_f &o1, const Psim_f &o2){
            mpfr_div(val, o1.val, o2.val, MPFR_RNDN);
        }
        Psim_f& operator/(const Psim_f& o){
            Psim_f &ret = tmps.get();
            mpfr_div(ret.val, val, o.val, MPFR_RNDN);
            return ret;
        }
        void div_d(const Psim_f &o1, double o2){
            mpfr_div_d(val, o1.val, o2, MPFR_RNDN);
        }
        Psim_f& operator/(const double& o){
            Psim_f &ret = tmps.get();
            mpfr_div_d(ret.val, val, o, MPFR_RNDN);
            return ret;
        }

        void log(const Psim_f &o){
            mpfr_log(val, o.val, MPFR_RNDN);
        }
        void exp(const Psim_f &o){
            mpfr_exp(val, o.val, MPFR_RNDN);
        }

        double get_d(){
            return mpfr_get_d(val, MPFR_RNDN);
        }
        Psim_f& operator=(double o){
            mpfr_set_d(val, o, MPFR_RNDN);
            return *this;
        }
        void set_d(double o){
            mpfr_set_d(val, o, MPFR_RNDN);
        }
        Psim_f& operator=(Psim_f &o){
            mpfr_set(val, o.val, MPFR_RNDN);
            return *this;
        }
        void set(Psim_f &o){
            mpfr_set(val, o.val, MPFR_RNDN);
        }

        Psim_f& operator+=(double o){
            *this = *this + o;
            return *this;
        }
        Psim_f& operator+=(Psim_f &o){
            *this = *this + o;
            return *this;
        }
        Psim_f& operator-=(double o){
            *this = *this - o;
            return *this;
        }
        Psim_f& operator-=(Psim_f &o){
            *this = *this - o;
            return *this;
        }

        int cmp_d(double o){
            return mpfr_cmp_d(val, o);
        }
        int cmp(Psim_f &o){
            return mpfr_cmp(val, o.val);
        }

        int min(Psim_f &o1, Psim_f &o2){
            return mpfr_min(val, o1.val, o2.val, MPFR_RNDN);
        }
        int max(Psim_f &o1, Psim_f &o2){
            return mpfr_max(val, o1.val, o2.val, MPFR_RNDN);
        }
    };
    
    class address_k_conc{
        public:
            /*
            mpfr_t conc_f;
            mpfr_t conc_r;
            mpfr_t conc_frc;
            mpfr_t conc_rrc;
            */

            //strands with all possible 5' and 3' ends
            //addr_end = original end of the strand
            //fstrand[0][0] = 5' F -> F Strand Data -> RRC 3'
            //rstrand[0][0] = 5' R -> R strand data -> FRC 3'
            //primerf = forward primer sequence
            //primerfrc = forward primer RC sequence
            //primerr = reverse primer sequence
            //primerrrc = reverse primer RC sequence
            Psim_f fstrand[5][5];
            Psim_f rstrand[5][5];

            Psim_f fstrand_change[5][5];
            Psim_f rstrand_change[5][5];

            Psim_f total_f_conc;
            Psim_f total_r_conc;
            Psim_f last_f_conc;
            Psim_f last_r_conc;

            //first index:
            //  primer_f = 0
            //  primer_r = 1

            //second index:
            //  addr_f = 0
            //  addr_frc = 1
            //  addr_r = 2
            //  addr_rrc = 3

            dh_ds dhds[2][4];

            double tmp_dhds[2][4][2];

            double dhds_primer_f_addr_f[2];
            double dhds_primer_f_addr_frc[2];
            double dhds_primer_f_addr_r[2];
            double dhds_primer_f_addr_rrc[2];
            double dhds_primer_r_addr_f[2];
            double dhds_primer_r_addr_frc[2];
            double dhds_primer_r_addr_r[2];
            double dhds_primer_r_addr_rrc[2];
    };

    class EQ{
        public:
            Psim_f last_val[2];
            Psim_f bounds[3];
            Psim_f solutions[3];
            Psim_f tmp[4];
            Psim_f c[19];
            Psim_f c0[6];
            Psim_f k[13];
            Psim_f spec_total;
            Psim_f nonspec_total;
            Psim_f spec_fwd_amp;
            Psim_f spec_rev_amp;
            Psim_f nonspec_exp_amp;
            Psim_f nonspec_lin_amp;
            Psim_f best_spec_exp_amp;
            Psim_f best_spec_lin_amp;
            Psim_f best_nonspec_exp_amp;
            Psim_f best_nonspec_lin_amp;
            Psim_f saved_frc_conc;
            Psim_f saved_rrc_conc;
            Psim_f avg_nonspec_amp;
            Psim_f last_nonspec_frc_total;
            Psim_f last_nonspec_rrc_total;
            std::vector<address_k_conc> address_k_conc_vec;

            void calc_cx();
            void calc_cf(Psim_f &ret);
            void calc_cr(Psim_f &ret);
            void calc_bound_concs();
            void solve_cf_cr(int F_or_R);
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
            //For multithreading address evaluation
            unsigned int address_index;
            std::mutex addr_mtx;
            std::mutex outfile_mtx;
            std::mutex mpfr_mtx;
            pthread_barrier_t bar;
            unsigned int num_cpu;


        public:
            Primeanneal();
            ~Primeanneal();
            void read_primers_individual(std::string filename);
            void swap_primer_rc(Primeanneal::primer_info *p);
            void reverse_comp(const char *fwd, char *rc);
            double dg_to_eq_const(double dg, double temp_c);
            double eq_const_to_dg(double eq_const, double temp_c);
            void dg_to_eq_const_mpfr(Psim_f &ret, double dg, double temp_c);
            double eq_const_to_dg_mpfr(Psim_f &tmp, Psim_f &eq_const, double temp_c);
            void assign_addresses(const char *out_filename, double temp_c, double mv_conc, double dv_conc, double dntp_conc);
            void assign_addresses_nosort(const char *out_filename);
            //void address_access(Primeanneal::address, double temp_c, double mv_conc, double dv_conc, double dntp_conc);
            //void address_access(Primeanneal::address addr, double *f_nsrc_dh, double *f_nsrc_ds, double *r_nsrc_dh, double *r_nsrc_ds,  double *f_nsf_dh, double *f_nsf_ds, double *r_nsf_dh, double *r_nsf_ds, double temp_c, double mv_conc, double dv_conc, double dntp_conc);
            thal_results calc_dimer(char *seq1, char *seq2, thal_args a);
            thal_results calc_hairpin(char *seq1, thal_args a);
            bool primer_compare(const Primeanneal::primer_info& a, const Primeanneal::primer_info& b);
            void read_addresses(const char *filename, bool with_temp_c);
            void evaluate_addresses(const char * in_filename, const char *out_filename, double dna_conc, double primer_conc, double mv_conc, double dv_conc, double dntp_conc);
            void eval_thread(const char * out_filename, double dna_conc, double primer_conc, double mv_conc, double dv_conc, double dntp_conc);
            //void sim_pcr_thread(unsigned int tid, unsigned int num_cpu, const char *out_filename, unsigned int addr, unsigned int pcr_cycles, double mv_conc, double dv_conc, double dntp_conc);
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

