#include <stdio.h>
#include <math.h>
#include <vector>
#include <string.h>
#include <stdlib.h>
#include <mpfr.h>
#include <thread>

#include "thal.h"
#include "eq.hpp"

namespace primersim{

    Primeanneal::Primeanneal(){
        num_cpu = 1;
    }

    Primeanneal::~Primeanneal(){
        for(auto &p : primers){
            free(p.f);
            free(p.rc);
        }
        primers.clear();
    }

    double Primeanneal::dg_to_eq_const(double dg, double temp_c){
        double temp_k = temp_c + 273.15;
        const double ideal_gas_constant = 1.98720425864083;
        return exp(dg / (-ideal_gas_constant*temp_k));
    }

    double Primeanneal::dhds_to_eq_const(double dhds[2], double temp_c){
        double dg = dhds[dh] - (temp_c+273.15) * dhds[ds];
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

    void Primeanneal::read_addresses(const char *filename, bool with_temp_c = false){
        addresses.clear(); //TODO: Fix potential memory leak here.
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
            a.f = (char *)malloc(strlen(tmp_f)+1);
            a.r = (char *)malloc(strlen(tmp_f)+1);
            a.f_rc = (char *)malloc(strlen(tmp_f)+1);
            a.r_rc = (char *)malloc(strlen(tmp_f)+1);
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
        eq.k[K_FX] = dhds_to_eq_const(eq.address_k_conc_vec[bind_addr].dhds_primer_f_addr_frc, temp_c_profile[cycle-1]);
        eq.k[K_RX] = dhds_to_eq_const(eq.address_k_conc_vec[bind_addr].dhds_primer_r_addr_frc, temp_c_profile[cycle-1]);
        eq.calc_cx();
        eq.calc_bound_concs();
        //Remove bound primer from primer concentration
        eq.c0[F] -= eq.c[FX];
        eq.c0[R] -= eq.c[RX];
        //Add to concentration of elongated strand
        eq.address_k_conc_vec[i].fstrand_change[primerf][end5_rc] += eq.c[FX];
        eq.address_k_conc_vec[i].fstrand_change[primerr][end5_rc] += eq.c[RX];

        //F Strand Bindings
        //[0][0]
        eq.c0[X] = eq.address_k_conc_vec[i].fstrand[end5][end3];
        eq.k[K_FX] = dhds_to_eq_const(eq.address_k_conc_vec[bind_addr].dhds_primer_f_addr_rrc, temp_c_profile[cycle-1]);
        eq.k[K_RX] = dhds_to_eq_const(eq.address_k_conc_vec[bind_addr].dhds_primer_r_addr_rrc, temp_c_profile[cycle-1]);
        eq.calc_cx();
        eq.calc_bound_concs();
        //Remove bound primer from primer concnetration
        eq.c0[F] -= eq.c[FX];
        eq.c0[R] -= eq.c[RX];
        //Add to concentration of elongated strand
        eq.address_k_conc_vec[i].rstrand_change[primerf][end5_rc] += eq.c[FX];
        eq.address_k_conc_vec[i].rstrand_change[primerr][end5_rc] += eq.c[RX];
    }

    void Primeanneal::calc_strand_bindings(EQ &eq, const std::vector<double> &temp_c_profile, int i, int cycle, int addr, int end5, int end3){
        //will be args
        //int end5 = 0;
        //int end3 = 0;
        double eq_const;
        //dhds indicies
        int bind_addr5;
        int bind_addr3;
        if (end5 == 0) bind_addr5 = i;
        else bind_addr5 = addr;
        if (end3 == 0) bind_addr3 = i;
        else bind_addr3 = addr;
        //primer f to 5' end
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
            case 0: binding_end3 = 1; break; //binding to R of target addr
            case 1: binding_end3 = 0; break; //binding to primer F
            case 2: binding_end3 = 1; break; //binding to primer FRC
            case 3: binding_end3 = 2; break; //binding to primer R
            case 4: binding_end3 = 3; break; //binding to primer RRC
            default: break;
        }

        // calc c[F] and c[R]
        // eq.tmp[0] will hold the sum of the nonspecific concentrations
        // eq.tmp[1] will hold the sum of (nonspec_conc * nonspec_k_f)
        // eq.tmp[2] will hold the sum of (nonspec_conc * nonspec_k_r)
        // eq.tmp[1] / eq.tmp[0] is average nonspec k for forward primer
        // eq.tmp[2] / eq.tmp[0] is average nonspec k for reverse primer
        // eq/tmp[3] will hold partial results
        // Reverse strand bindings:

        // 5' R -> R Strand -> FRC 3' [0][0]
        // 5' binding
        eq.tmp[0] += eq.address_k_conc_vec[i].rstrand[end5][end3];

        eq_const = dhds_to_eq_const(eq.address_k_conc_vec[i].tmp_dhds[0][binding_end5], temp_c_profile[cycle - 1]);
        eq.tmp[1] += eq.address_k_conc_vec[bind_addr5].rstrand[end5][end3] * eq_const;

        eq_const = dhds_to_eq_const(eq.address_k_conc_vec[i].tmp_dhds[1][binding_end5], temp_c_profile[cycle - 1]);
        eq.tmp[2] += eq.address_k_conc_vec[bind_addr5].rstrand[end5][end3] * eq_const;

        // 3' Binding
        eq.tmp[0] += eq.address_k_conc_vec[i].rstrand[end5][end3];

        eq_const = dhds_to_eq_const(eq.address_k_conc_vec[i].tmp_dhds[0][binding_end3], temp_c_profile[cycle - 1]);
        eq.tmp[1] += eq.address_k_conc_vec[bind_addr3].rstrand[end5][end3] * eq_const;

        eq_const = dhds_to_eq_const(eq.address_k_conc_vec[i].tmp_dhds[1][binding_end3], temp_c_profile[cycle - 1]);
        eq.tmp[2] += eq.address_k_conc_vec[bind_addr3].rstrand[end5][end3] * eq_const;



        if (end5 == 0) binding_end5 = 0; //Address F
        if (end3 == 0) binding_end3 = 3; //Address RRC
        // 5' F -> F Strand -> RRC 3' [0][0]
        //5' binding
        eq.tmp[0] += eq.address_k_conc_vec[i].fstrand[end5][end3];

        eq_const = dhds_to_eq_const(eq.address_k_conc_vec[i].tmp_dhds[0][binding_end5], temp_c_profile[cycle-1]);
        eq.tmp[1] += eq.address_k_conc_vec[bind_addr5].fstrand[end5][end3] * eq_const;

        eq_const = dhds_to_eq_const(eq.address_k_conc_vec[i].tmp_dhds[1][binding_end5], temp_c_profile[cycle-1]);
        eq.tmp[2] += eq.address_k_conc_vec[bind_addr5].fstrand[end5][end3] * eq_const;

        //3' Binding
        eq.tmp[0] += eq.address_k_conc_vec[i].fstrand[end5][end3];

        eq_const = dhds_to_eq_const(eq.address_k_conc_vec[i].tmp_dhds[0][binding_end3], temp_c_profile[cycle-1]);
        eq.tmp[1] += eq.address_k_conc_vec[bind_addr3].fstrand[end5][end3] * eq_const;

        eq_const = dhds_to_eq_const(eq.address_k_conc_vec[i].tmp_dhds[1][binding_end3], temp_c_profile[cycle-1]);
        eq.tmp[2] += eq.address_k_conc_vec[bind_addr3].fstrand[end5][end3] * eq_const;
    }

    double Primeanneal::sim_pcr(const char *out_filename, unsigned int addr, unsigned int pcr_cycles, const std::vector<double> &temp_c_profile, double dna_conc, double primer_f_conc, double primer_r_conc, double mv_conc, double dv_conc, double dntp_conc){
        //addresses.clear();
        //read_addresses(in_filename, true);
        EQ eq;
        thal_args ta;
        thal_results o;
        set_thal_default_args(&ta);
        ta.mv = mv_conc;
        ta.dv = dv_conc;
        ta.dntp = dntp_conc;
        //ta.temp = 273.15 + addresses[addr].temp_c;
        ta.temp = 273.15 + 55;
        //double temp_c = addresses[addr].temp_c;
        eq.address_k_conc_vec.resize(addresses.size());
        eq.last_nonspec_frc_total = 0.;
        eq.last_nonspec_rrc_total = 0.;

        for(unsigned int i = 0; i < addresses.size(); i++){
            eq.address_k_conc_vec[i].fstrand[0][0] = dna_conc / addresses.size();
            eq.address_k_conc_vec[i].rstrand[0][0] = dna_conc / addresses.size();
            for(int ii = 0; ii < 5; ii++){
                for(int jj = 0; jj < 5; jj++){
                    if (!jj && !ii)
                        continue;
                    eq.address_k_conc_vec[i].fstrand[ii][jj] = 0.0;
                    eq.address_k_conc_vec[i].rstrand[ii][jj] = 0.0;
                }
            }

            for (int ii = 0; ii < 2; ii++){
                for (int jj = 0; jj < 4; jj++){
                    o = calc_dimer(addresses[addr].get_seq(ii), addresses[i].get_seq(jj), ta);
                    eq.address_k_conc_vec[i].tmp_dhds[ii][jj][dh] = o.dh;
                    eq.address_k_conc_vec[i].tmp_dhds[ii][jj][ds] = o.ds;
                }
            }

            o = calc_dimer(addresses[addr].f, addresses[i].f, ta);
            eq.address_k_conc_vec[i].dhds_primer_f_addr_f[dh] = o.dh;
            eq.address_k_conc_vec[i].dhds_primer_f_addr_f[ds] = o.ds;
            //eq.address_k_conc_vec[i].k_primer_f_addr_f.set_d(dg_to_eq_const(o.dg, temp_c));
            o = calc_dimer(addresses[addr].f, addresses[i].f_rc, ta);
            eq.address_k_conc_vec[i].dhds_primer_f_addr_frc[dh] = o.dh;
            eq.address_k_conc_vec[i].dhds_primer_f_addr_frc[ds] = o.ds;
            //eq.address_k_conc_vec[i].k_primer_f_addr_frc.set_d(dg_to_eq_const(o.dg, temp_c));
            o = calc_dimer(addresses[addr].f, addresses[i].r, ta);
            eq.address_k_conc_vec[i].dhds_primer_f_addr_r[dh] = o.dh;
            eq.address_k_conc_vec[i].dhds_primer_f_addr_r[ds] = o.ds;
            //eq.address_k_conc_vec[i].k_primer_f_addr_r.set_d(dg_to_eq_const(o.dg, temp_c));
            o = calc_dimer(addresses[addr].f, addresses[i].r_rc, ta);
            eq.address_k_conc_vec[i].dhds_primer_f_addr_rrc[dh] = o.dh;
            eq.address_k_conc_vec[i].dhds_primer_f_addr_rrc[ds] = o.ds;
            //eq.address_k_conc_vec[i].k_primer_f_addr_rrc.set_d(dg_to_eq_const(o.dg, temp_c));
            o = calc_dimer(addresses[addr].r, addresses[i].f, ta);
            eq.address_k_conc_vec[i].dhds_primer_r_addr_f[dh] = o.dh;
            eq.address_k_conc_vec[i].dhds_primer_r_addr_f[ds] = o.ds;
            //eq.address_k_conc_vec[i].k_primer_r_addr_f.set_d(dg_to_eq_const(o.dg, temp_c));
            o = calc_dimer(addresses[addr].r, addresses[i].f_rc, ta);
            eq.address_k_conc_vec[i].dhds_primer_r_addr_frc[dh] = o.dh;
            eq.address_k_conc_vec[i].dhds_primer_r_addr_frc[ds] = o.ds;
            //eq.address_k_conc_vec[i].k_primer_r_addr_frc.set_d(dg_to_eq_const(o.dg, temp_c));
            o = calc_dimer(addresses[addr].r, addresses[i].r, ta);
            eq.address_k_conc_vec[i].dhds_primer_r_addr_r[dh] = o.dh;
            eq.address_k_conc_vec[i].dhds_primer_r_addr_r[ds] = o.ds;
            //eq.address_k_conc_vec[i].k_primer_r_addr_r.set_d(dg_to_eq_const(o.dg, temp_c));
            o = calc_dimer(addresses[addr].r, addresses[i].r_rc, ta);
            eq.address_k_conc_vec[i].dhds_primer_r_addr_rrc[dh] = o.dh;
            eq.address_k_conc_vec[i].dhds_primer_r_addr_rrc[ds] = o.ds;
            //eq.address_k_conc_vec[i].k_primer_r_addr_rrc.set_d(dg_to_eq_const(o.dg, temp_c));
        }

        double dhds_f_hairpin[2];
        o = calc_hairpin(addresses[addr].f, ta);
        dhds_f_hairpin[dh] = o.dh;
        dhds_f_hairpin[ds] = o.ds;

        double dhds_r_hairpin[2];
        o = calc_hairpin(addresses[addr].r, ta);
        dhds_r_hairpin[dh] = o.dh;
        dhds_r_hairpin[ds] = o.ds;

        eq.c0[F] = primer_f_conc;
        eq.c0[R] = primer_r_conc;

        if(out_filename != NULL){
            FILE *outfile = fopen(out_filename, "w");
            //cycle,f_primer,r_primer,spec_f,spec_r,nonspec_f,nonspec_r
            mpfr_fprintf(outfile, "%02u,%lf,0.000000000,0.000000000,%.9Re,%.9Re,%.9Re,%.9Re", 0, temp_c_profile[0], eq.c0[F].val, eq.c0[R].val, eq.address_k_conc_vec[addr].fstrand[addr_end][addr_end].val, eq.address_k_conc_vec[addr].rstrand[addr_end][addr_end].val);
            eq.tmp[0] = 0.;
            eq.tmp[1] = 0.;
            for(unsigned int i = 0; i < addresses.size(); i++){
                if (i == addr)
                    continue;
                eq.tmp[0] = eq.tmp[0] + eq.address_k_conc_vec[i].fstrand[addr_end][addr_end];
                eq.tmp[1] = eq.tmp[1] + eq.address_k_conc_vec[i].rstrand[addr_end][addr_end];
            }
            mpfr_fprintf(outfile, ",%.9Re,0.000000000,%.9Re,0.000000000\n", eq.tmp[0].val, eq.tmp[1].val);
            fclose(outfile);
        }
        for(unsigned int cycle = 1; cycle <= pcr_cycles; cycle++){
            //calc c[F] and c[R]
            //eq.tmp[0] will hold the sum of the nonspecific concentrations
            //eq.tmp[1] will hold the sum of (nonspec_conc * nonspec_k_f)
            //eq.tmp[2] will hold the sum of (nonspec_conc * nonspec_k_r)
            //eq.tmp[1] / eq.tmp[0] is average nonspec k for forward primer
            //eq.tmp[2] / eq.tmp[0] is average nonspec k for reverse primer
            //eq/tmp[3] will hold partial results

            //Set the eq constant of all specific bindings, primer dimers, and primer hairpins
            eq.k[K_FF] = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_f, temp_c_profile[cycle-1]);
            eq.k[K_FR] = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_r, temp_c_profile[cycle-1]);
            eq.k[K_RR] = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_r, temp_c_profile[cycle-1]);
            eq.k[K_FH] = dhds_to_eq_const(dhds_f_hairpin, temp_c_profile[cycle-1]);
            eq.k[K_RH] = dhds_to_eq_const(dhds_r_hairpin, temp_c_profile[cycle-1]);

            eq.tmp[0] = 0.;//Total concentration
            eq.tmp[1] = 0.;//sum(f primer eq const * partial concentration)
            eq.tmp[2] = 0.;//sum(r primer eq const * partial concentration)

            for(unsigned int i = 0; i < addresses.size(); i++){
                for(int end5 = 0; end5 < 5; end5++){
                    for (int end3 = 0; end3 < 5; end3++){
                        calc_strand_bindings(eq, temp_c_profile, i, cycle, addr, end5, end3);
                    }
                }
            }
            //Set k[K_FX] and k[K_RX] to average eq constant
            eq.k[K_FX].div(eq.tmp[1], eq.tmp[0]);
            eq.k[K_RX].div(eq.tmp[2], eq.tmp[0]);
            //All k values now set

            //Set initial concentrations
            //c0[F] and c0[R] are set by previous iteration or initialized before the loop
            //eq.tmp[0] still holds the total nonspecific concentration
            eq.c0[X].set(eq.tmp[0]);
            //All initial concs now set

            eq.solve_eq();

            //Now solve individual nonspecific concentrations and update c0 for next cycle
            for(unsigned int i = 0; i < addresses.size(); i++){
                eq.address_k_conc_vec[i].last_f_conc.set_d(0.0);
                eq.address_k_conc_vec[i].last_r_conc.set_d(0.0);

                for(int ii = 0; ii < 5; ii++){
                    for(int jj = 0; jj < 5; jj++){
                        eq.address_k_conc_vec[i].last_f_conc.add(eq.address_k_conc_vec[i].last_f_conc, eq.address_k_conc_vec[i].fstrand[ii][jj]);
                        eq.address_k_conc_vec[i].last_r_conc.add(eq.address_k_conc_vec[i].last_r_conc, eq.address_k_conc_vec[i].rstrand[ii][jj]);

                        eq.address_k_conc_vec[i].fstrand_change[ii][jj].set_d(0.0);
                        eq.address_k_conc_vec[i].rstrand_change[ii][jj].set_d(0.0);
                    }
                }

                for (int end5 = 0; end5 < 5; end5++){
                    for (int end3 = 0 ; end3 < 5; end3++){
                        update_strand_concs(eq, i, addr, temp_c_profile, cycle, end5, end3);
                    }
                }

                //Add changes to strand concentrations and sum all f and r strands
                eq.address_k_conc_vec[i].total_f_conc.set_d(0.0);
                eq.address_k_conc_vec[i].total_r_conc.set_d(0.0);
                for(int ii = 0; ii < 5; ii++){
                    for(int jj = 0; jj < 5; jj++){
                        eq.address_k_conc_vec[i].fstrand[ii][jj].add(eq.address_k_conc_vec[i].fstrand[ii][jj], eq.address_k_conc_vec[i].fstrand_change[ii][jj]);
                        eq.address_k_conc_vec[i].rstrand[ii][jj].add(eq.address_k_conc_vec[i].rstrand[ii][jj], eq.address_k_conc_vec[i].rstrand_change[ii][jj]);
                        eq.address_k_conc_vec[i].total_f_conc.add(eq.address_k_conc_vec[i].total_f_conc, eq.address_k_conc_vec[i].fstrand[ii][jj]);
                        eq.address_k_conc_vec[i].total_r_conc.add(eq.address_k_conc_vec[i].total_r_conc, eq.address_k_conc_vec[i].rstrand[ii][jj]);
                    }
                }

            }

            eq.tmp[0].div(eq.address_k_conc_vec[addr].total_f_conc, eq.address_k_conc_vec[addr].last_f_conc);
            eq.tmp[0].sub_d(eq.tmp[0], 1.0);

            eq.tmp[1].div(eq.address_k_conc_vec[addr].total_r_conc, eq.address_k_conc_vec[addr].last_r_conc);
            eq.tmp[1].sub_d(eq.tmp[1], 1.0);

            eq.spec_total.add(eq.address_k_conc_vec[addr].total_f_conc, eq.address_k_conc_vec[addr].total_r_conc);

            if(out_filename != NULL){
                FILE *outfile = fopen(out_filename, "a");
                //cycle,f_primer,r_primer,spec_f,spec_r,nonspec_f,nonspec_r
                mpfr_fprintf(outfile, "%02u,%lf,%.9Rf,%.9Rf,%.9Re,%.9Re,%.9Re,%.9Re", cycle, temp_c_profile[cycle-1], eq.tmp[0].val, eq.tmp[1].val, eq.c0[F].val, eq.c0[R].val, eq.address_k_conc_vec[addr].total_f_conc.val, eq.address_k_conc_vec[addr].total_r_conc.val);
                fclose(outfile);
            }
            eq.tmp[0].set_d(0.);
            eq.tmp[1].set_d(0.);
            eq.tmp[2].set_d(0.);
            eq.tmp[3].set_d(0.);
            eq.nonspec_total.set_d(0.0);
            for(unsigned int i = 0; i < addresses.size(); i++){
                if (i == addr)
                    continue;
                eq.tmp[0].add(eq.tmp[0], eq.address_k_conc_vec[i].total_f_conc);
                eq.tmp[1].add(eq.tmp[1], eq.address_k_conc_vec[i].last_f_conc);
                eq.tmp[2].add(eq.tmp[2], eq.address_k_conc_vec[i].total_r_conc);
                eq.tmp[3].add(eq.tmp[3], eq.address_k_conc_vec[i].last_r_conc);
                eq.nonspec_total.add(eq.nonspec_total, eq.address_k_conc_vec[i].total_f_conc);
                eq.nonspec_total.add(eq.nonspec_total, eq.address_k_conc_vec[i].total_r_conc);
            }
            eq.tmp[1].div(eq.tmp[0], eq.tmp[1]);
            eq.tmp[1].sub_d(eq.tmp[1], 1.);
            eq.tmp[3].div(eq.tmp[2], eq.tmp[3]);
            eq.tmp[3].sub_d(eq.tmp[3], 1.);
            eq.last_nonspec_frc_total.set(eq.tmp[0]);
            eq.last_nonspec_rrc_total.set(eq.tmp[1]);
            if(out_filename != NULL){
                FILE *outfile = fopen(out_filename, "a");
                mpfr_fprintf(outfile, ",%.9Re,%.9Rf,%.9Re,%.9Rf\n", eq.tmp[0].val, eq.tmp[1].val, eq.tmp[2].val, eq.tmp[3].val);
                fclose(outfile);
            }
        }

        eq.tmp[0].div(eq.spec_total, eq.nonspec_total);
        return eq.tmp[0].get_d();
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
