#include <stdio.h>
#include <math.h>
#include <vector>
#include <string.h>
#include <stdlib.h>
#include <thread>

#include "thal.h"
#include "eq.hpp"

namespace primersim{

    Primeanneal::Primeanneal(){
        num_cpu = 1;
    }

    Primeanneal::~Primeanneal(){
        for(auto &p : primers){
            delete[] p.f;
            delete[] p.rc;
        }
        primers.clear();
        for(auto &a : addresses){
            delete[] a.f;
            delete[] a.r;
            delete[] a.f_rc;
            delete[] a.r_rc;
        }
        addresses.clear();
    }

    double Primeanneal::dg_to_eq_const(double dg, double temp_c){
        double temp_k = temp_c + 273.15;
        const double ideal_gas_constant = 1.98720425864083;
        return exp(dg / (-ideal_gas_constant*temp_k));
    }

    double Primeanneal::dhds_to_eq_const(ThalReal dhds[2], double temp_c){
        // dhds entries promote to double for the math; sim_pcr's chain
        // keeps full double precision regardless of ThalReal.
        double dg = (double)dhds[dh] - (temp_c+273.15) * (double)dhds[ds];
        return dg_to_eq_const(dg, temp_c);
    }

    void Primeanneal::reverse_comp(const char *fwd, char *rc){
        int end = strlen(fwd)-1;
        for (int i = 0; i <= end; i++){
            if (fwd[i] == 'A')
                rc[end-i] = 'T';
            else if (fwd[i] == 'T')
                rc[end-i] = 'A';
            else if (fwd[i] == 'G')
                rc[end-i] = 'C';
            else if (fwd[i] == 'C')
                rc[end-i] = 'G';
        }
        rc[end+1] = '\0';
    }

    thal_results Primeanneal::calc_dimer(char *seq1, char *seq2, thal_args a){
        unsigned char *seq1_uchar = (unsigned char *)seq1;
        unsigned char *seq2_uchar = (unsigned char *)seq2;
        thal_results o;
        thal(seq1_uchar, seq2_uchar, &a, THL_FAST, &o);
        return o;
    }

    thal_results Primeanneal::calc_hairpin(char *seq1, thal_args a){
        unsigned char *seq1_uchar = (unsigned char *)seq1;
        a.type = thal_hairpin;
        thal_results o;
        thal(seq1_uchar, seq1_uchar, &a, THL_FAST, &o);
        return o;
    }

    void Primeanneal::populate_dimer_cache(double mv_conc, double dv_conc, double dntp_conc, double temp_c){
        thal_args ta;
        set_thal_default_args(&ta);
        ta.mv = mv_conc;
        ta.dv = dv_conc;
        ta.dntp = dntp_conc;
        ta.temp = 273.15 + temp_c;

        // Symmetric / triangular cache: thal(seq1, seq2) is approximately
        // symmetric in its arguments, so we store only one of (a, b) and
        // (b, a). With M = 4*N "sequence slots" (each of N addresses has
        // 4 kinds: f, f_rc, r, r_rc), the cache holds the M*(M+1)/2
        // canonical pairs (a, b) with a <= b. Index formula:
        //   idx(a, b) = a*(2*M - a - 1)/2 + b
        // Lookup callers swap to a <= b before indexing — see
        // cached_dimer in sim_pcr.
        const size_t N = addresses.size();
        const size_t M = 4 * N;
        dimer_cache.assign(M * (M + 1) / 2, {0, 0});

        // Cyclic distribution by `a` for load balance: each thread handles
        // a's where a % nthreads == t. Triangular-row work decreases as
        // a grows (row a has M-a entries), so a contiguous slice would
        // give thread 0 most of the work and thread N-1 almost none.
        const unsigned int nthreads = num_cpu > 0 ? num_cpu : 1;
        std::vector<std::thread> threads;
        threads.reserve(nthreads);
        for (unsigned int t = 0; t < nthreads; t++) {
            threads.emplace_back([this, t, nthreads, M, ta]() {
                thal_args local_ta = ta;
                for (size_t a = t; a < M; a += nthreads) {
                    size_t i = a / 4;
                    int   ki = (int)(a % 4);
                    size_t row_offset = a * (2*M - a - 1) / 2;
                    for (size_t b = a; b < M; b++) {
                        size_t j = b / 4;
                        int   kj = (int)(b % 4);
                        thal_results o = this->calc_dimer(
                            this->addresses[i].get_seq(ki),
                            this->addresses[j].get_seq(kj),
                            local_ta);
                        auto &p = this->dimer_cache[row_offset + b];
                        p.dh = (ThalReal)o.dh;
                        p.ds = (ThalReal)o.ds;
                    }
                }
            });
        }
        for (auto &th : threads) th.join();
    }

    void Primeanneal::read_addresses(const char *filename, bool with_temp_c = false){
        // Free any previously-allocated address strings before clearing
        // the vector — otherwise repeated calls leak.
        for(auto &a : addresses){
            delete[] a.f;
            delete[] a.r;
            delete[] a.f_rc;
            delete[] a.r_rc;
        }
        addresses.clear();
        FILE *infile = fopen(filename, "r");
        address a;
        char tmp_f[100];
        char tmp_r[100];
        while(true){
            if(with_temp_c){
                if(fscanf(infile, "%[ACGT],%[ACGT],%lf,%*[^\n]\n", tmp_f, tmp_r, &(a.temp_c)) == EOF)
                    break;
            } else {
                if(fscanf(infile, "%[ACGT],%[ACGT]%*[^\n]\n", tmp_f, tmp_r) == EOF)
                    break;
            }
            a.f    = new char[strlen(tmp_f)+1];
            a.r    = new char[strlen(tmp_r)+1];
            a.f_rc = new char[strlen(tmp_f)+1];
            a.r_rc = new char[strlen(tmp_r)+1];
            strcpy(a.f, tmp_f);
            strcpy(a.r, tmp_r);
            reverse_comp(a.f, a.f_rc);
            reverse_comp(a.r, a.r_rc);
            addresses.push_back(a);
        }
    }

    void Primeanneal::update_strand_concs(EQ &eq, int i, int addr, const std::vector<double> &temp_c_profile, int cycle, int end5, int end3){
        int end5_rc = 0;
        switch (end5){
            case 0: end5_rc = 0; break; //addr_end
            case 1: end5_rc = 2; break; //1 = primer f, 2 = primer frc
            case 2: end5_rc = 1; break;
            case 3: end5_rc = 4; break; //3 = primer r, 4 = primer rrc
            case 4: end5_rc = 3; break;
            default: break;
        }

        int bind_addr = i;
        if (end3 != 0) bind_addr = addr;

        //R strand bindings
        eq.c0[X] = eq.address_k_conc_vec[i].rstrand[end5][end3];
        eq.k[K_FX] = eq.address_k_conc_vec[bind_addr].dhds_K_primer_f_addr_frc;
        eq.k[K_RX] = eq.address_k_conc_vec[bind_addr].dhds_K_primer_r_addr_frc;
        eq.calc_cx();
        eq.calc_bound_concs();
        //Remove bound primer from primer concentration
        eq.c0[F] -= eq.c[FX];
        eq.c0[R] -= eq.c[RX];
        //Add to concentration of elongated strand
        eq.address_k_conc_vec[i].fstrand_change[primerf][end5_rc] += eq.c[FX];
        eq.address_k_conc_vec[i].fstrand_change[primerr][end5_rc] += eq.c[RX];

        //F Strand Bindings
        eq.c0[X] = eq.address_k_conc_vec[i].fstrand[end5][end3];
        eq.k[K_FX] = eq.address_k_conc_vec[bind_addr].dhds_K_primer_f_addr_rrc;
        eq.k[K_RX] = eq.address_k_conc_vec[bind_addr].dhds_K_primer_r_addr_rrc;
        eq.calc_cx();
        eq.calc_bound_concs();
        //Remove bound primer from primer concentration
        eq.c0[F] -= eq.c[FX];
        eq.c0[R] -= eq.c[RX];
        //Add to concentration of elongated strand
        eq.address_k_conc_vec[i].rstrand_change[primerf][end5_rc] += eq.c[FX];
        eq.address_k_conc_vec[i].rstrand_change[primerr][end5_rc] += eq.c[RX];
    }

    // Called 25 times per address per cycle from sim_pcr's inner loop.
    // The three Real& parameters are accumulators owned by the caller:
    //   nonspec_total   — sum of nonspecific strand concentrations across
    //                     all (end5, end3) pairs and addresses
    //   sum_f_weighted  — sum(strand_conc * F-primer-binding K), used by
    //                     the caller to derive the average k_FX
    //   sum_r_weighted  — same but for the R primer / k_RX
    void Primeanneal::calc_strand_bindings(EQ &eq, const std::vector<double> &temp_c_profile, int i, int cycle, int addr, int end5, int end3,
                                            Real &nonspec_total, Real &sum_f_weighted, Real &sum_r_weighted){
        (void)temp_c_profile; (void)cycle;  // K values now come from the cache
        int bind_addr5 = (end5 == 0) ? i : addr;
        int bind_addr3 = (end3 == 0) ? i : addr;

        int binding_end5 = 0;
        switch(end5){
            case 0: binding_end5 = 2; break; //binding to R of target addr
            case 1: binding_end5 = 0; break; //binding to primer F
            case 2: binding_end5 = 1; break; //binding to primer FRC
            case 3: binding_end5 = 2; break; //binding to primer R
            case 4: binding_end5 = 3; break; //binding to primer RRC
            default: break;
        }
        int binding_end3 = 0;
        switch(end3){
            case 0: binding_end3 = 1; break;
            case 1: binding_end3 = 0; break;
            case 2: binding_end3 = 1; break;
            case 3: binding_end3 = 2; break;
            case 4: binding_end3 = 3; break;
            default: break;
        }

        // K values come from the per-cycle cache populated in sim_pcr.
        const auto &K_i = eq.address_k_conc_vec[i].tmp_dhds_K;

        // Reverse strand bindings
        // 5' R -> R Strand -> FRC 3'
        nonspec_total += eq.address_k_conc_vec[i].rstrand[end5][end3];
        sum_f_weighted += eq.address_k_conc_vec[bind_addr5].rstrand[end5][end3] * K_i[0][binding_end5];
        sum_r_weighted += eq.address_k_conc_vec[bind_addr5].rstrand[end5][end3] * K_i[1][binding_end5];

        // 3' Binding
        nonspec_total += eq.address_k_conc_vec[i].rstrand[end5][end3];
        sum_f_weighted += eq.address_k_conc_vec[bind_addr3].rstrand[end5][end3] * K_i[0][binding_end3];
        sum_r_weighted += eq.address_k_conc_vec[bind_addr3].rstrand[end5][end3] * K_i[1][binding_end3];


        if (end5 == 0) binding_end5 = 0; //Address F
        if (end3 == 0) binding_end3 = 3; //Address RRC
        // 5' F -> F Strand -> RRC 3'
        nonspec_total += eq.address_k_conc_vec[i].fstrand[end5][end3];
        sum_f_weighted += eq.address_k_conc_vec[bind_addr5].fstrand[end5][end3] * K_i[0][binding_end5];
        sum_r_weighted += eq.address_k_conc_vec[bind_addr5].fstrand[end5][end3] * K_i[1][binding_end5];

        // 3' Binding
        nonspec_total += eq.address_k_conc_vec[i].fstrand[end5][end3];
        sum_f_weighted += eq.address_k_conc_vec[bind_addr3].fstrand[end5][end3] * K_i[0][binding_end3];
        sum_r_weighted += eq.address_k_conc_vec[bind_addr3].fstrand[end5][end3] * K_i[1][binding_end3];
    }

    double Primeanneal::sim_pcr(const char *out_filename, unsigned int addr, unsigned int pcr_cycles, const std::vector<double> &temp_c_profile, double dna_conc, double primer_f_conc, double primer_r_conc, double mv_conc, double dv_conc, double dntp_conc){
        EQ eq;
        thal_args ta;
        thal_results o;
        set_thal_default_args(&ta);
        ta.mv = mv_conc;
        ta.dv = dv_conc;
        ta.dntp = dntp_conc;
        ta.temp = 273.15 + 55;
        eq.address_k_conc_vec.resize(addresses.size());
        eq.last_nonspec_frc_total = 0.0;
        eq.last_nonspec_rrc_total = 0.0;

        for(unsigned int i = 0; i < addresses.size(); i++){
            eq.address_k_conc_vec[i].fstrand[0][0] = (Real)dna_conc / addresses.size();
            eq.address_k_conc_vec[i].rstrand[0][0] = (Real)dna_conc / addresses.size();
            for(int ii = 0; ii < 5; ii++){
                for(int jj = 0; jj < 5; jj++){
                    if (!jj && !ii)
                        continue;
                    eq.address_k_conc_vec[i].fstrand[ii][jj] = 0.0;
                    eq.address_k_conc_vec[i].rstrand[ii][jj] = 0.0;
                }
            }

            // Read pre-computed dh/ds values from dimer_cache instead of
            // calling calc_dimer. The cache is symmetric / triangular —
            // it stores only canonical pairs (a <= b in the M = 4*N
            // sequence-slot space). Swap if needed so we hit the same
            // entry regardless of argument order.
            const size_t N_addr = addresses.size();
            const size_t M_total = 4 * N_addr;
            auto cached = [this, M_total](size_t a_addr, int ka, size_t b_addr, int kb) -> const dh_ds_cache_entry& {
                size_t a = a_addr * 4 + ka;
                size_t b = b_addr * 4 + kb;
                if (a > b) std::swap(a, b);
                return this->dimer_cache[a * (2*M_total - a - 1) / 2 + b];
            };

            for (int ii = 0; ii < 2; ii++){
                for (int jj = 0; jj < 4; jj++){
                    const auto &p = cached(addr, ii, i, jj);
                    eq.address_k_conc_vec[i].tmp_dhds[ii][jj][dh] = p.dh;
                    eq.address_k_conc_vec[i].tmp_dhds[ii][jj][ds] = p.ds;
                }
            }

            // Map: addresses[].get_seq() returns f at 0, f_rc at 1, r at 2, r_rc at 3.
            { const auto &p = cached(addr, 0, i, 0); eq.address_k_conc_vec[i].dhds_primer_f_addr_f[dh]   = p.dh; eq.address_k_conc_vec[i].dhds_primer_f_addr_f[ds]   = p.ds; }
            { const auto &p = cached(addr, 0, i, 1); eq.address_k_conc_vec[i].dhds_primer_f_addr_frc[dh] = p.dh; eq.address_k_conc_vec[i].dhds_primer_f_addr_frc[ds] = p.ds; }
            { const auto &p = cached(addr, 0, i, 2); eq.address_k_conc_vec[i].dhds_primer_f_addr_r[dh]   = p.dh; eq.address_k_conc_vec[i].dhds_primer_f_addr_r[ds]   = p.ds; }
            { const auto &p = cached(addr, 0, i, 3); eq.address_k_conc_vec[i].dhds_primer_f_addr_rrc[dh] = p.dh; eq.address_k_conc_vec[i].dhds_primer_f_addr_rrc[ds] = p.ds; }
            { const auto &p = cached(addr, 2, i, 0); eq.address_k_conc_vec[i].dhds_primer_r_addr_f[dh]   = p.dh; eq.address_k_conc_vec[i].dhds_primer_r_addr_f[ds]   = p.ds; }
            { const auto &p = cached(addr, 2, i, 1); eq.address_k_conc_vec[i].dhds_primer_r_addr_frc[dh] = p.dh; eq.address_k_conc_vec[i].dhds_primer_r_addr_frc[ds] = p.ds; }
            { const auto &p = cached(addr, 2, i, 2); eq.address_k_conc_vec[i].dhds_primer_r_addr_r[dh]   = p.dh; eq.address_k_conc_vec[i].dhds_primer_r_addr_r[ds]   = p.ds; }
            { const auto &p = cached(addr, 2, i, 3); eq.address_k_conc_vec[i].dhds_primer_r_addr_rrc[dh] = p.dh; eq.address_k_conc_vec[i].dhds_primer_r_addr_rrc[ds] = p.ds; }
        }

        ThalReal dhds_f_hairpin[2];
        o = calc_hairpin(addresses[addr].f, ta);
        dhds_f_hairpin[dh] = o.dh;
        dhds_f_hairpin[ds] = o.ds;

        ThalReal dhds_r_hairpin[2];
        o = calc_hairpin(addresses[addr].r, ta);
        dhds_r_hairpin[dh] = o.dh;
        dhds_r_hairpin[ds] = o.ds;

        eq.c0[F] = primer_f_conc;
        eq.c0[R] = primer_r_conc;

        if(out_filename != NULL){
            // Each cycle row written by sim_pcr has 12 comma-separated columns.
            // Cycle 0 fills columns 3, 4, 10, 12 with literal zeros (no growth
            // measurement is meaningful before any cycle has run); subsequent
            // cycles fill them with computed ratios. Layout:
            //
            //   1  cycle              (00..N, %02u)
            //   2  temp_c             (the cycle's temperature in C)
            //   3  f_growth_ratio     (target_F_total / target_F_last - 1)
            //   4  r_growth_ratio     (target_R_total / target_R_last - 1)
            //   5  c0[F]              (free F primer remaining)
            //   6  c0[R]              (free R primer remaining)
            //   7  total_f_conc       (target's F-strand total)
            //   8  total_r_conc       (target's R-strand total)
            //   9  total_nonspec_f    (sum of F-strand totals across non-target addrs)
            //  10  nonspec_f_growth   (total_nonspec_f / last_nonspec_f - 1)
            //  11  total_nonspec_r    (sum of R-strand totals across non-target addrs)
            //  12  nonspec_r_growth   (total_nonspec_r / last_nonspec_r - 1)
            FILE *outfile = fopen(out_filename, "a");
            fprintf(outfile, "%02u,%lf,0.000000000,0.000000000,%.9Le,%.9Le,%.9Le,%.9Le", 0, temp_c_profile[0], (long double)eq.c0[F], (long double)eq.c0[R], (long double)eq.address_k_conc_vec[addr].fstrand[addr_end][addr_end], (long double)eq.address_k_conc_vec[addr].rstrand[addr_end][addr_end]);
            // Cycle-0 dump: sum the initial fstrand/rstrand totals across
            // all non-target addresses for the nonspec columns.
            Real init_sum_fstrand = 0.0;
            Real init_sum_rstrand = 0.0;
            for(unsigned int i = 0; i < addresses.size(); i++){
                if (i == addr)
                    continue;
                init_sum_fstrand += eq.address_k_conc_vec[i].fstrand[addr_end][addr_end];
                init_sum_rstrand += eq.address_k_conc_vec[i].rstrand[addr_end][addr_end];
            }
            fprintf(outfile, ",%.9Le,0.000000000,%.9Le,0.000000000\n", (long double)init_sum_fstrand, (long double)init_sum_rstrand);
            fclose(outfile);
        }
        for(unsigned int cycle = 1; cycle <= pcr_cycles; cycle++){
            eq.k[K_FF] = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_f, temp_c_profile[cycle-1]);
            eq.k[K_FR] = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_r, temp_c_profile[cycle-1]);
            eq.k[K_RR] = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_r, temp_c_profile[cycle-1]);
            eq.k[K_FH] = dhds_to_eq_const(dhds_f_hairpin, temp_c_profile[cycle-1]);
            eq.k[K_RH] = dhds_to_eq_const(dhds_r_hairpin, temp_c_profile[cycle-1]);

            // Populate per-address K caches for this cycle. Each value
            // would otherwise be recomputed 25 times by calc_strand_bindings
            // and update_strand_concs.
            for(unsigned int i = 0; i < addresses.size(); i++){
                auto &kc = eq.address_k_conc_vec[i];
                for (int p = 0; p < 2; p++)
                    for (int e = 0; e < 4; e++)
                        kc.tmp_dhds_K[p][e] = dhds_to_eq_const(kc.tmp_dhds[p][e], temp_c_profile[cycle-1]);
                kc.dhds_K_primer_f_addr_frc = dhds_to_eq_const(kc.dhds_primer_f_addr_frc, temp_c_profile[cycle-1]);
                kc.dhds_K_primer_r_addr_frc = dhds_to_eq_const(kc.dhds_primer_r_addr_frc, temp_c_profile[cycle-1]);
                kc.dhds_K_primer_f_addr_rrc = dhds_to_eq_const(kc.dhds_primer_f_addr_rrc, temp_c_profile[cycle-1]);
                kc.dhds_K_primer_r_addr_rrc = dhds_to_eq_const(kc.dhds_primer_r_addr_rrc, temp_c_profile[cycle-1]);
            }

            // Per-cycle accumulators populated by calc_strand_bindings:
            //   nonspec_total   — total nonspecific strand concentration
            //   sum_f_weighted  — sum(strand_conc * F-primer K), used to
            //                     derive the average k_FX
            //   sum_r_weighted  — same for R primer / k_RX
            Real nonspec_total  = 0.0;
            Real sum_f_weighted = 0.0;
            Real sum_r_weighted = 0.0;

            for(unsigned int i = 0; i < addresses.size(); i++){
                for(int end5 = 0; end5 < 5; end5++){
                    for (int end3 = 0; end3 < 5; end3++){
                        calc_strand_bindings(eq, temp_c_profile, i, cycle, addr, end5, end3,
                                              nonspec_total, sum_f_weighted, sum_r_weighted);
                    }
                }
            }
            //Set k[K_FX] and k[K_RX] to average eq constant
            eq.k[K_FX] = sum_f_weighted / nonspec_total;
            eq.k[K_RX] = sum_r_weighted / nonspec_total;

            eq.c0[X] = nonspec_total;

            eq.solve_eq();

            //Solve individual nonspecific concentrations and update c0 for next cycle
            for(unsigned int i = 0; i < addresses.size(); i++){
                eq.address_k_conc_vec[i].last_f_conc = 0.0;
                eq.address_k_conc_vec[i].last_r_conc = 0.0;

                for(int ii = 0; ii < 5; ii++){
                    for(int jj = 0; jj < 5; jj++){
                        eq.address_k_conc_vec[i].last_f_conc += eq.address_k_conc_vec[i].fstrand[ii][jj];
                        eq.address_k_conc_vec[i].last_r_conc += eq.address_k_conc_vec[i].rstrand[ii][jj];

                        eq.address_k_conc_vec[i].fstrand_change[ii][jj] = 0.0;
                        eq.address_k_conc_vec[i].rstrand_change[ii][jj] = 0.0;
                    }
                }

                for (int end5 = 0; end5 < 5; end5++){
                    for (int end3 = 0 ; end3 < 5; end3++){
                        update_strand_concs(eq, i, addr, temp_c_profile, cycle, end5, end3);
                    }
                }

                //Add changes and sum all f and r strands
                eq.address_k_conc_vec[i].total_f_conc = 0.0;
                eq.address_k_conc_vec[i].total_r_conc = 0.0;
                for(int ii = 0; ii < 5; ii++){
                    for(int jj = 0; jj < 5; jj++){
                        eq.address_k_conc_vec[i].fstrand[ii][jj] += eq.address_k_conc_vec[i].fstrand_change[ii][jj];
                        eq.address_k_conc_vec[i].rstrand[ii][jj] += eq.address_k_conc_vec[i].rstrand_change[ii][jj];
                        eq.address_k_conc_vec[i].total_f_conc += eq.address_k_conc_vec[i].fstrand[ii][jj];
                        eq.address_k_conc_vec[i].total_r_conc += eq.address_k_conc_vec[i].rstrand[ii][jj];
                    }
                }
            }

            // Per-cycle target growth ratios for the active address's F
            // and R strands (printed in the cycle's first output row).
            Real f_growth_ratio = eq.address_k_conc_vec[addr].total_f_conc / eq.address_k_conc_vec[addr].last_f_conc - 1.0;
            Real r_growth_ratio = eq.address_k_conc_vec[addr].total_r_conc / eq.address_k_conc_vec[addr].last_r_conc - 1.0;

            eq.spec_total = eq.address_k_conc_vec[addr].total_f_conc + eq.address_k_conc_vec[addr].total_r_conc;

            if(out_filename != NULL){
                // Columns 1..8 of the cycle row. See the cycle-0 fprintf
                // above for the full 12-column layout.
                FILE *outfile = fopen(out_filename, "a");
                fprintf(outfile, "%02u,%lf,%.9Lf,%.9Lf,%.9Le,%.9Le,%.9Le,%.9Le", cycle, temp_c_profile[cycle-1], (long double)f_growth_ratio, (long double)r_growth_ratio, (long double)eq.c0[F], (long double)eq.c0[R], (long double)eq.address_k_conc_vec[addr].total_f_conc, (long double)eq.address_k_conc_vec[addr].total_r_conc);
                fclose(outfile);
            }

            // Nonspec totals across all non-target addresses, for the
            // continuation row of the cycle's output.
            Real total_nonspec_f = 0.0;
            Real last_nonspec_f  = 0.0;
            Real total_nonspec_r = 0.0;
            Real last_nonspec_r  = 0.0;
            eq.nonspec_total = 0.0;
            for(unsigned int i = 0; i < addresses.size(); i++){
                if (i == addr)
                    continue;
                total_nonspec_f += eq.address_k_conc_vec[i].total_f_conc;
                last_nonspec_f  += eq.address_k_conc_vec[i].last_f_conc;
                total_nonspec_r += eq.address_k_conc_vec[i].total_r_conc;
                last_nonspec_r  += eq.address_k_conc_vec[i].last_r_conc;
                eq.nonspec_total += eq.address_k_conc_vec[i].total_f_conc;
                eq.nonspec_total += eq.address_k_conc_vec[i].total_r_conc;
            }
            Real nonspec_f_growth = total_nonspec_f / last_nonspec_f - 1.0;
            Real nonspec_r_growth = total_nonspec_r / last_nonspec_r - 1.0;
            eq.last_nonspec_frc_total = total_nonspec_f;
            eq.last_nonspec_rrc_total = nonspec_f_growth;
            if(out_filename != NULL){
                // Columns 9..12 of the cycle row, completing the line. See
                // the cycle-0 fprintf for the full 12-column layout.
                FILE *outfile = fopen(out_filename, "a");
                fprintf(outfile, ",%.9Le,%.9Lf,%.9Le,%.9Lf\n", (long double)total_nonspec_f, (long double)nonspec_f_growth, (long double)total_nonspec_r, (long double)nonspec_r_growth);
                fclose(outfile);
            }
        }

        return (double)(eq.spec_total / eq.nonspec_total);
    }


    void Primeanneal::eval_addresses_thread(const char *out_filename, unsigned int pcr_cycles, const std::vector<double> &temp_c_profile, double dna_conc, double primer_f_conc, double primer_r_conc, double mv_conc, double dv_conc, double dntp_conc){
        while(true){
            unsigned int i;
            addr_mtx.lock();
            i = address_index;
            address_index++;
            addr_mtx.unlock();
            if(i >= addresses.size())
                break;
            double ratio = sim_pcr(NULL, i, temp_c_profile.size(), temp_c_profile, dna_conc, primer_f_conc, primer_r_conc, mv_conc, dv_conc, dntp_conc);
            outfile_mtx.lock();
            FILE *outfile = fopen(out_filename, "a");
            fprintf(outfile, "%u,%e\n", i, ratio);
            fclose(outfile);
            outfile_mtx.unlock();
        }
    }

    void Primeanneal::eval_addresses(const char *out_filename, unsigned int pcr_cycles, const std::vector<double> &temp_c_profile, double dna_conc, double primer_f_conc, double primer_r_conc, double mv_conc, double dv_conc, double dntp_conc){
        address_index = 0;
        std::vector<std::thread> threads;
        for(unsigned int i = 0; i < num_cpu; i++){
            threads.push_back(std::thread(&Primeanneal::eval_addresses_thread, this, out_filename, temp_c_profile.size(), temp_c_profile,  dna_conc, primer_f_conc, primer_r_conc, mv_conc, dv_conc, dntp_conc));
        }

        for (auto &t : threads)
            t.join();

    }
}
