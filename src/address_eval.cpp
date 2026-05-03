#include <stdio.h>
#include <math.h>
#include <vector>
#include <string.h>
#include <string>
#include <stdlib.h>
#include <cmath>
#include <algorithm>
#include <thread>

#include "thal.h"
#include "eq.hpp"

namespace primersim{

    double Primeanneal::eq_const_to_dg(double eq_const, double temp_c){
        double temp_k = temp_c + 273.15;
        const double ideal_gas_constant = 1.98720425864083;
        return log(eq_const) * (-ideal_gas_constant*temp_k);
    }

    void Primeanneal::dg_to_eq_const_mpfr(Real &ret, double dg, double temp_c){
        double temp_k = temp_c + 273.15;
        ret = std::exp((Real)(dg / (-1.98720425864083 * temp_k)));
    }

    double Primeanneal::eq_const_to_dg_mpfr(Real &tmp, Real &eq_const, double temp_c){
        double temp_k = temp_c + 273.15;
        tmp = std::log(eq_const) * (Real)(-1.98720425864083 * temp_k);
        return (double)tmp;
    }

    void Primeanneal::read_primers_individual(std::string filename){
        FILE *infile = fopen(filename.c_str(), "r");
        char buffer[50];
        // Free any previously-allocated primer strings before clearing
        // the vector — otherwise repeated calls leak.
        for(auto &p : primers){
            delete[] p.f;
            delete[] p.rc;
        }
        primers.clear();
        while(fscanf(infile, "%[ACGT]%*[^\n]\n", buffer) != EOF){
            //printf("%s\n", buffer);
            int len = strlen(buffer)+1;
            if(len != 21)
                printf("%s\n", buffer);
            Primeanneal::primer_info primer_pair;
            primer_pair.f = new char[len];
            primer_pair.rc = new char[len];
            strcpy(primer_pair.f, buffer);
            if(strlen(primer_pair.f) != 20)
                printf("pp:%s\n", primer_pair.f);
            reverse_comp(primer_pair.f, primer_pair.rc);
            primers.push_back(primer_pair);
            if(strlen(primers[primers.size()-1].f) != 20)
                printf("%s\n", primers[primers.size()-1].f);
            if(strlen(primers[primers.size()-1].rc) != 20)
                printf("%s\n", primers[primers.size()-1].rc);
        }
        fclose(infile);
    }

    void Primeanneal::swap_primer_rc(Primeanneal::primer_info *p){
        char *tmp = p->f;
        p->f = p->rc;
        p->rc = tmp;
    }

    bool Primeanneal::primer_compare(const Primeanneal::primer_info& a, const Primeanneal::primer_info& b){
        return (a.nonspec_dg_f_rc - a.spec_dg > b.nonspec_dg_f_rc - b.spec_dg);
    }

    void Primeanneal::assign_addresses(const char *out_filename, double temp_c, double mv_conc, double dv_conc, double dntp_conc){
        thal_args a;
        set_thal_default_args(&a);
        a.mv = mv_conc;
        a.dv = dv_conc;
        a.dntp = dntp_conc;
        a.temp = temp_c + 273.15;
        int half_num_nonspec = primers.size() - 1;
        bool swap = true;
        int iterations = 0;
        while(swap){
            swap = false;
            int num_swaps = 0;
            int count = 0;
            for (auto &p : primers){
                printf("%d\n", count++);
                if(strlen(p.f) != 20)
                    printf("%s\n", p.f);
                if(strlen(p.rc) != 20)
                    printf("%s\n", p.rc);
                p.spec_dg = calc_dimer(p.f, p.rc, a).dg;
                double f_nonspec_k = 0.0;
                double rc_nonspec_k = 0.0;
                double f_f_nonspec_k = 0.;
                double f_rc_nonspec_k = 0.;
                double rc_f_nonspec_k = 0.;
                double rc_rc_nonspec_k = 0.;
                for(auto &p2 : primers){
                    if(!strcmp(p.f, p2.f))
                        continue;
                    //Fwd-Fwd
                    double nonspec_dg;// = calc_dimer(p.f, p2.f, a).dg;
                    //f_nonspec_k += dg_to_eq_const(nonspec_dg, temp_c) / num_nonspec;
                    //f_f_nonspec_k += dg_to_eq_const(nonspec_dg, temp_c) / half_num_nonspec;
                    //Fwd-RC
                    nonspec_dg = calc_dimer(p.f, p2.rc, a).dg;
                    //f_nonspec_k += dg_to_eq_const(nonspec_dg, temp_c) / num_nonspec;
                    f_rc_nonspec_k += dg_to_eq_const(nonspec_dg, temp_c) / half_num_nonspec;
                    //RC-Fwd
                    //nonspec_dg = calc_dimer(p.rc, p2.f, a).dg;
                    //rc_nonspec_k += dg_to_eq_const(nonspec_dg, temp_c) / num_nonspec;
                    //rc_f_nonspec_k += dg_to_eq_const(nonspec_dg, temp_c) / half_num_nonspec;
                    //RC-RC
                    nonspec_dg = calc_dimer(p.rc, p2.rc, a).dg;
                    //rc_nonspec_k += dg_to_eq_const(nonspec_dg, temp_c) / num_nonspec;
                    rc_rc_nonspec_k += dg_to_eq_const(nonspec_dg, temp_c) / half_num_nonspec;
                }
                p.nonspec_dg_f = eq_const_to_dg(f_nonspec_k, temp_c);
                p.nonspec_dg_rc = eq_const_to_dg(rc_nonspec_k, temp_c);
                p.nonspec_dg_f_f = eq_const_to_dg(f_f_nonspec_k, temp_c);
                p.nonspec_dg_f_rc = eq_const_to_dg(f_rc_nonspec_k, temp_c);
                p.nonspec_dg_rc_f = eq_const_to_dg(rc_f_nonspec_k, temp_c);
                p.nonspec_dg_rc_rc = eq_const_to_dg(rc_rc_nonspec_k, temp_c);
                if(p.nonspec_dg_f_rc - p.spec_dg < p.nonspec_dg_rc_rc - p.spec_dg){
                    swap_primer_rc(&p);
                    swap = true;
                    num_swaps++;
                }
            }
            iterations++;
            printf("%d: %d\n", iterations, num_swaps);
        }
        std::sort(primers.begin(), primers.end(), std::bind(&Primeanneal::primer_compare, this, std::placeholders::_1, std::placeholders::_2));
        FILE *outfile = fopen("addresses.csv", "w+");
        int num_addresses = primers.size() / 2;
        for (int i = 0; i < num_addresses; i++){
            fprintf(outfile, "%s,%s,%lf,%lf,%lf,%lf,%lf,%lf\n", primers[i].f, primers[i+num_addresses].f,primers[i].spec_dg,
                                primers[i+num_addresses].spec_dg, primers[i].nonspec_dg_f_rc, primers[i+num_addresses].nonspec_dg_f_rc,
                                primers[i].nonspec_dg_f_rc - primers[i].spec_dg,
                                primers[i+num_addresses].nonspec_dg_f_rc - primers[i+num_addresses].spec_dg);
        }
        fclose(outfile);
    }

    void Primeanneal::assign_addresses_nosort(const char *out_filename){
        FILE *outfile = fopen(out_filename, "w+");
        int num_addresses = primers.size() / 2;
        for (int i = 0; i < num_addresses; i++){
            fprintf(outfile, "%s,%s,%lf,%lf,%lf,%lf,%lf,%lf\n", primers[i].f, primers[i+num_addresses].f,primers[i].spec_dg,
                                primers[i+num_addresses].spec_dg, primers[i].nonspec_dg_f_rc, primers[i+num_addresses].nonspec_dg_f_rc,
                                primers[i].nonspec_dg_f_rc - primers[i].spec_dg,
                                primers[i+num_addresses].nonspec_dg_f_rc - primers[i+num_addresses].spec_dg);
        }
        fclose(outfile);
    }

    void Primeanneal::eval_thread(const char * out_filename, double dna_conc, double primer_conc, double mv_conc, double dv_conc, double dntp_conc){
        EQ eq;
        eq.c0[F] = primer_conc;
        eq.c0[R] = primer_conc;
        thal_results o;
        thal_args ta;
        set_thal_default_args(&ta);
        ta.mv = mv_conc;
        ta.dv = dv_conc;
        ta.dntp = dntp_conc;
        ThalReal f_hp_dh;
        ThalReal f_hp_ds;
        ThalReal r_hp_dh;
        ThalReal r_hp_ds;
        struct dh_ds_info{
            ThalReal f_f_dh;
            ThalReal f_f_ds;
            ThalReal f_r_dh;
            ThalReal f_r_ds;
            ThalReal f_frc_dh;
            ThalReal f_frc_ds;
            ThalReal f_rrc_dh;
            ThalReal f_rrc_ds;
            ThalReal r_f_dh;
            ThalReal r_f_ds;
            ThalReal r_r_dh;
            ThalReal r_r_ds;
            ThalReal r_frc_dh;
            ThalReal r_frc_ds;
            ThalReal r_rrc_dh;
            ThalReal r_rrc_ds;
        };

        std::vector<struct dh_ds_info> dh_ds;
        dh_ds.resize(addresses.size());
        while(true){
            unsigned int i = 0;
            addr_mtx.lock();
            i = address_index;
            address_index++;
            addr_mtx.unlock();
            if (i >= addresses.size()){
                break;
            }

            double best_temp = 0.0;
            bool first = true;
            o = calc_hairpin(addresses[i].f, ta);
            f_hp_dh = o.dh;
            f_hp_ds = o.ds;
            o = calc_hairpin(addresses[i].r, ta);
            r_hp_dh = o.dh;
            r_hp_ds = o.ds;
            for (unsigned int j = 0; j < addresses.size(); j++){
                struct dh_ds_info s;
                o = calc_dimer(addresses[i].f, addresses[j].f, ta);
                s.f_f_dh = o.dh;
                s.f_f_ds = o.ds;
                o = calc_dimer(addresses[i].f, addresses[j].r, ta);
                s.f_r_dh = o.dh;
                s.f_r_ds = o.ds;
                o = calc_dimer(addresses[i].f, addresses[j].f_rc, ta);
                s.f_frc_dh = o.dh;
                s.f_frc_ds = o.ds;
                o = calc_dimer(addresses[i].f, addresses[j].r_rc, ta);
                s.f_rrc_dh = o.dh;
                s.f_rrc_ds = o.ds;
                o = calc_dimer(addresses[i].r, addresses[j].f, ta);
                s.r_f_dh = o.dh;
                s.r_f_ds = o.ds;
                o = calc_dimer(addresses[i].r, addresses[j].r, ta);
                s.r_r_dh = o.dh;
                s.r_r_ds = o.ds;
                o = calc_dimer(addresses[i].r, addresses[j].f_rc, ta);
                s.r_frc_dh = o.dh;
                s.r_frc_ds = o.ds;
                o = calc_dimer(addresses[i].r, addresses[j].r_rc, ta);
                s.r_rrc_dh = o.dh;
                s.r_rrc_ds = o.ds;
                dh_ds[j] = s;
            }

            eq.best_spec_exp_amp = 0.0;
            eq.best_spec_lin_amp = 0.0;
            eq.best_nonspec_exp_amp = 0.0;
            eq.best_nonspec_lin_amp = 0.0;

            for(double temp_c = 40.; temp_c < 80.01; temp_c += 1.){
                eq.c0[X] = dna_conc * ((4.0 * (addresses.size() - 1.0)) + 2.0); //4 nonspecific binding sites per double strand DNA, plus two 3' binding sites on target DNA

                eq.spec_fwd_amp = 0.0;
                eq.spec_rev_amp = 0.0;
                eq.nonspec_exp_amp = 0.0;
                eq.nonspec_lin_amp = 0.0;

                eq.k[K_FF] = dg_to_eq_const(dh_ds[i].f_f_dh - (dh_ds[i].f_f_ds * (temp_c + 273.15)), temp_c);
                eq.k[K_RR] = dg_to_eq_const(dh_ds[i].r_r_dh - (dh_ds[i].r_r_ds * (temp_c + 273.15)), temp_c);
                eq.k[K_FR] = dg_to_eq_const(dh_ds[i].f_r_dh - (dh_ds[i].f_r_ds * (temp_c + 273.15)), temp_c);

                eq.k[K_FH] = dg_to_eq_const(f_hp_dh - (f_hp_ds * (temp_c + 273.15)), temp_c);
                eq.k[K_RH] = dg_to_eq_const(r_hp_dh - (r_hp_ds * (temp_c + 273.15)), temp_c);

                eq.k[K_FX] = 0.;
                eq.k[K_RX] = 0.;
                int count = 0;
                for(unsigned int j = 0; j < addresses.size(); j++){
                    eq.k[K_FX] += dg_to_eq_const(dh_ds[j].f_f_dh - (dh_ds[j].f_f_ds * (temp_c + 273.15)), temp_c);
                    eq.k[K_RX] += dg_to_eq_const(dh_ds[j].r_f_dh - (dh_ds[j].r_f_ds * (temp_c + 273.15)), temp_c);
                    count++;

                    eq.k[K_FX] += dg_to_eq_const(dh_ds[j].f_frc_dh - (dh_ds[j].f_frc_ds * (temp_c + 273.15)), temp_c);
                    eq.k[K_RX] += dg_to_eq_const(dh_ds[j].r_frc_dh - (dh_ds[j].r_frc_ds * (temp_c + 273.15)), temp_c);
                    count++;

                    eq.k[K_FX] += dg_to_eq_const(dh_ds[j].f_r_dh - (dh_ds[j].f_r_ds * (temp_c + 273.15)), temp_c);
                    eq.k[K_RX] += dg_to_eq_const(dh_ds[j].r_r_dh - (dh_ds[j].r_r_ds * (temp_c + 273.15)), temp_c);
                    count++;

                    eq.k[K_FX] += dg_to_eq_const(dh_ds[j].f_rrc_dh - (dh_ds[j].f_rrc_ds * (temp_c + 273.15)), temp_c);
                    eq.k[K_RX] += dg_to_eq_const(dh_ds[j].r_rrc_dh - (dh_ds[j].r_rrc_ds * (temp_c + 273.15)), temp_c);
                    count++;
                }
                eq.k[K_FX] = eq.k[K_FX] / count;
                eq.k[K_RX] = eq.k[K_RX] / count;
                eq.solve_eq();
                //c[F] and c[R] are now solved, all nonspec concs can now be solved individually

                //Set c[X] to concs of individual strands
                eq.c0[X] = dna_conc / addresses.size();
                for (unsigned int j = 0; j < addresses.size(); j++){
                    if (i == j)
                        continue;
                    eq.k[K_FX] = dg_to_eq_const(dh_ds[j].f_frc_dh - (dh_ds[j].f_frc_ds * (temp_c + 273.15)), temp_c);
                    eq.k[K_RX] = dg_to_eq_const(dh_ds[j].r_frc_dh - (dh_ds[j].r_frc_ds * (temp_c + 273.15)), temp_c);
                    eq.calc_cx();
                    eq.calc_bound_concs();

                    // bound_total — total primer bound to this nonspec address
                    // Note: the original code mixed this with an
                    // uninitialized tmp[1] in the min/max split (see
                    // docs/address_suggestions.md #1). Replacing tmp[1]
                    // with a 0-initialized local takes the degenerate-
                    // but-defined branch noted there: min collapses to 0,
                    // so all bound mass goes to nonspec_lin_amp.
                    Real bound_total = eq.c[FX] + eq.c[RX];
                    Real other_pool = 0.0;  // was eq.tmp[1], never written in this loop
                    Real bound_min = std::fmin(bound_total, other_pool);
                    eq.nonspec_exp_amp += bound_min;
                    Real bound_max = std::fmax(bound_total, other_pool);
                    eq.nonspec_lin_amp += (bound_max - bound_min);
                }
                // cur_ratio  = spec_fwd_amp / nonspec_exp_amp at this temp_c
                // best_ratio = the running best across temp_c iterations
                Real cur_ratio  = eq.spec_fwd_amp / eq.nonspec_exp_amp;
                Real best_ratio = eq.best_spec_exp_amp / eq.best_nonspec_exp_amp;
                if(first || (cur_ratio > best_ratio)){
                    first = false;
                    best_temp = temp_c;
                    eq.best_spec_exp_amp = eq.spec_fwd_amp;
                    eq.best_spec_lin_amp = eq.spec_rev_amp;
                    eq.best_nonspec_exp_amp = eq.nonspec_exp_amp;
                    eq.best_nonspec_lin_amp = eq.nonspec_lin_amp;
                }
            }
            // Final spec/nonspec ratios for the output row.
            Real best_exp_ratio = eq.best_spec_exp_amp / eq.best_nonspec_exp_amp;
            Real best_lin_ratio = eq.best_spec_lin_amp / eq.best_nonspec_lin_amp;
            outfile_mtx.lock();
            FILE *outfile = fopen(out_filename, "a");
            fprintf(outfile, "%s,%s,%lf,%u,%.9Le,%.9Le,%.9Le,%.9Le,%.9Le,%.9Le\n", addresses[i].f, addresses[i].r, best_temp, i,
                    (long double)best_exp_ratio, (long double)best_lin_ratio,
                    (long double)eq.best_spec_exp_amp, (long double)eq.best_nonspec_exp_amp,
                    (long double)eq.best_spec_lin_amp, (long double)eq.best_nonspec_lin_amp);
            fclose(outfile);
            outfile_mtx.unlock();
        }
    }

    void Primeanneal::evaluate_addresses(const char *in_filename, const char * out_filename, double dna_conc, double primer_conc, double mv_conc, double dv_conc, double dntp_conc){
        if(addresses.size() == 0)
            read_addresses(in_filename, false);
        //Clear the outfile
        FILE *outfile = fopen(out_filename, "w");
        fclose(outfile);
        std::vector<std::thread> threads;
        address_index = 0;
        for(unsigned int i = 0; i < num_cpu; i++){
            threads.push_back(std::thread(&Primeanneal::eval_thread, this, out_filename, dna_conc, primer_conc, mv_conc, dv_conc, dntp_conc));
        }

        for (auto &t : threads)
            t.join();
    }

    void Primeanneal::shuffle_addresses(void){
        int dest_idx;
        int dest_rc;
        int dest_fr;
        for(unsigned int i = 0; i < addresses.size(); i++){
            char *tmp;

            dest_idx = std::rand() % addresses.size();
            dest_rc = std::rand() % 2;
            dest_fr = std::rand() % 2;

            if(!dest_rc && !dest_fr){
                tmp = addresses[dest_idx].f;
                addresses[dest_idx].f = addresses[i].f;
                addresses[i].f = tmp;
                tmp = addresses[dest_idx].f_rc;
                addresses[dest_idx].f_rc = addresses[i].f_rc;
                addresses[i].f_rc = tmp;
            }
            if(!dest_rc && dest_fr){
                tmp = addresses[dest_idx].r;
                addresses[dest_idx].r = addresses[i].f;
                addresses[i].f = tmp;
                tmp = addresses[dest_idx].r_rc;
                addresses[dest_idx].r_rc = addresses[i].f_rc;
                addresses[i].f_rc = tmp;
            }
            if(dest_rc && !dest_fr){
                tmp = addresses[dest_idx].f_rc;
                addresses[dest_idx].f_rc = addresses[i].f;
                addresses[i].f = tmp;
                tmp = addresses[dest_idx].f;
                addresses[dest_idx].f = addresses[i].f_rc;
                addresses[i].f_rc = tmp;
            }
            if(dest_rc && dest_fr){
                tmp = addresses[dest_idx].r;
                addresses[dest_idx].r = addresses[i].f;
                addresses[i].f = tmp;
                tmp = addresses[dest_idx].r_rc;
                addresses[dest_idx].r_rc = addresses[i].f_rc;
                addresses[i].f_rc = tmp;
            }

            dest_idx = std::rand() % addresses.size();
            dest_rc = std::rand() % 2;
            dest_fr = std::rand() % 2;

            if(!dest_rc && !dest_fr){
                tmp = addresses[dest_idx].f;
                addresses[dest_idx].f = addresses[i].r;
                addresses[i].r = tmp;
                tmp = addresses[dest_idx].f_rc;
                addresses[dest_idx].f_rc = addresses[i].r_rc;
                addresses[i].r_rc = tmp;
            }
            if(!dest_rc && dest_fr){
                tmp = addresses[dest_idx].r;
                addresses[dest_idx].r = addresses[i].r;
                addresses[i].r = tmp;
                tmp = addresses[dest_idx].r_rc;
                addresses[dest_idx].r_rc = addresses[i].r_rc;
                addresses[i].r_rc = tmp;
            }
            if(dest_rc && !dest_fr){
                tmp = addresses[dest_idx].f_rc;
                addresses[dest_idx].f_rc = addresses[i].r;
                addresses[i].r = tmp;
                tmp = addresses[dest_idx].f;
                addresses[dest_idx].f = addresses[i].r_rc;
                addresses[i].r_rc = tmp;
            }
            if(dest_rc && dest_fr){
                tmp = addresses[dest_idx].r;
                addresses[dest_idx].r = addresses[i].r;
                addresses[i].r = tmp;
                tmp = addresses[dest_idx].r_rc;
                addresses[dest_idx].r_rc = addresses[i].r_rc;
                addresses[i].r_rc = tmp;
            }
        }
    }
}
