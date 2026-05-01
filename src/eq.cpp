#include <stdio.h>
#include <math.h>
#include <vector>
#include <string.h>
#include <string>
#include <stdlib.h>
#include <mpfr.h>
#include <algorithm>
#include <thread>

#include "thal.h"
#include "eq.hpp"

namespace primersim{

    void EQ::print_state(std::string out_filename, std::string s){
        int i;
        FILE *outfile = fopen(out_filename.c_str(), "a");
        fprintf(outfile, "%s,", s.c_str());
        for(i = 0; i < 2; i++)
            mpfr_fprintf(outfile, "last_val[%d],%.9Re,", i, last_val[i].val);
        for(i = 0; i < 4; i++)
            mpfr_fprintf(outfile, "tmp[%d],%.9Re,",i,tmp[i].val);
        for(i = 0; i < 19; i++)
            mpfr_fprintf(outfile, "c[%d],%.9Re,",i,c[i].val);
        for(i = 0; i < 6; i++)
            mpfr_fprintf(outfile, "c0[%d],%.9Re,",i,c0[i].val);
        for(i = 0; i < 13; i++)
            mpfr_fprintf(outfile, "k[%d],%.9Re,",i,k[i].val);
        mpfr_fprintf(outfile, "spec_exp,%.9Re,",spec_fwd_amp.val);
        mpfr_fprintf(outfile, "spec_lin,%.9Re,",spec_rev_amp.val);
        mpfr_fprintf(outfile, "nonspec_exp,%.9Re,",nonspec_exp_amp.val);
        mpfr_fprintf(outfile, "nonspec_lin%.9Re,",nonspec_lin_amp.val);
        mpfr_fprintf(outfile, "best_spec_exp,%.9Re,",best_spec_exp_amp.val);
        mpfr_fprintf(outfile, "best_spec_lin,%.9Re,",best_spec_lin_amp.val);
        mpfr_fprintf(outfile, "best_nonspec_exp,%.9Re,",best_nonspec_exp_amp.val);
        mpfr_fprintf(outfile, "best_nonspec_lin,%.9Re\n",best_nonspec_lin_amp.val);
        fclose(outfile);
    }

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

    double Primeanneal::eq_const_to_dg(double eq_const, double temp_c){
        double temp_k = temp_c + 273.15;
        const double ideal_gas_constant = 1.98720425864083;
        return log(eq_const) * (-ideal_gas_constant*temp_k);
    }

    void Primeanneal::dg_to_eq_const_mpfr(Psim_f &ret, double dg, double temp_c){
        double temp_k = temp_c + 273.15;
        ret.set_d(dg / (-1.98720425864083 * temp_k));
        ret.exp(ret);
    }

    
    double Primeanneal::dhds_to_eq_const(double dhds[2], double temp_c){
        double dg = dhds[dh] - (temp_c+273.15) * dhds[ds];
        return dg_to_eq_const(dg, temp_c);
    }

    double Primeanneal::eq_const_to_dg_mpfr(Psim_f &tmp, Psim_f &eq_const, double temp_c){
        double temp_k = temp_c + 273.15;
        tmp.log(eq_const);
        tmp.mul_d(tmp, -1.98720425864083 * temp_k);
        return tmp.get_d();
    }

    //c[F], c[R], c0[:], and k[:] must be defined
    void EQ::calc_cx(){
        //c[X] = c0[X] / (1 + k[K_FX]*c[F] + k[K_RX]*c[R])
        tmp[2].mul(k[K_FX], c[F]);
        tmp[3].mul(k[K_RX], c[R]);
        tmp[2].add(tmp[2], tmp[3]);
        tmp[2].add_d(tmp[2], 1.0);
        c[X].div(c0[X], tmp[2]);
    }

    void EQ::calc_bound_concs(){
        c[FH] = k[K_FH] * c[F];
        c[RH] = k[K_RH] * c[R];
        c[FF] = k[K_FF] * c[F] * c[F];
        c[RR] = k[K_RR] * c[R] * c[R];
        c[FR] = k[K_FR] * c[F] * c[R];
        c[FX] = k[K_FX] * c[F] * c[X];
        c[RX] = k[K_RX] * c[R] * c[X];
    }

    // Solves the coupled equilibrium for c[F], c[R], and c[X].
    // Inner Newton solver handles the F/R 2x2 system with c[X] held constant;
    // outer loop refreshes c[X] = c0[X] / (1 + k_FX*c[F] + k_RX*c[R]) until
    // c[F] and c[R] stop changing.
    void EQ::solve_eq(){
        constexpr int max_iter_outer = FLOAT_PREC * 2;
        constexpr int max_iter = 15;

        Psim_f B1, B2;
        Psim_f f, r;
        Psim_f F1, F2;
        Psim_f J11, J12, J21, J22;
        Psim_f det, df, dr;
        Psim_f f_new, r_new;
        Psim_f tmp1, tmp2;
        Psim_f prev_res;

        c[F] = c0[F] / 2.0;
        c[R] = c0[R] / 2.0;
        calc_cx();
        last_val[F] = c[F];
        last_val[R] = c[R];

        for (int outer = 0; outer < max_iter_outer; ++outer) {
            mpfr_set_inf(prev_res.val, 1);

            B1 = k[K_RX]*c[X] + 1.0;
            B2 = k[K_FX]*c[X] + k[K_FH] + 1.0;

            // Decoupled initial guesses (exact when k_FR = K_RH = 0)
            mpfr_sqrt(tmp1.val, (B2*B2 + k[K_FF]*c0[F]*8.0).val, MPFR_RNDN);
            f = (B2*(-1.0) + tmp1) / (k[K_FF]*4.0);
            mpfr_sqrt(tmp1.val, (B1*B1 + k[K_RR]*c0[R]*8.0).val, MPFR_RNDN);
            r = (B1*(-1.0) + tmp1) / (k[K_RR]*4.0);

            for (int iter = 0; iter < max_iter; ++iter) {
                F1 = k[K_RR]*r*r*2.0 + (B1 + k[K_FR]*f)*r + k[K_RH]*f - c0[R];
                F2 = k[K_FF]*f*f*2.0 + (B2 + k[K_FR]*r)*f - c0[F];

                // Stop when residual stops decreasing — that's the precision
                // floor of mpfr arithmetic; further iterations can't improve it.
                mpfr_abs(tmp1.val, F1.val, MPFR_RNDN);
                mpfr_abs(tmp2.val, F2.val, MPFR_RNDN);
                mpfr_add(tmp1.val, tmp1.val, tmp2.val, MPFR_RNDN);
                if (tmp1.cmp(prev_res) >= 0) break;
                mpfr_set(prev_res.val, tmp1.val, MPFR_RNDN);

                J11 = k[K_RH] + k[K_FR]*r;
                J12 = k[K_RR]*r*4.0 + B1 + k[K_FR]*f;
                J21 = k[K_FF]*f*4.0 + B2 + k[K_FR]*r;
                J22 = k[K_FR]*f;

                det = J11*J22 - J12*J21;
                df  = (J12*F2 - J22*F1) / det;
                dr  = (J21*F1 - J11*F2) / det;

                f_new = f + df;
                r_new = r + dr;

                if      (f_new.cmp_d(0.0) < 0) f_new = f*0.5;
                else if (f_new.cmp(c0[F]) > 0) f_new = (f + c0[F])*0.5;
                if      (r_new.cmp_d(0.0) < 0) r_new = r*0.5;
                else if (r_new.cmp(c0[R]) > 0) r_new = (r + c0[R])*0.5;

                f = f_new;
                r = r_new;
            }

            c[F] = f;
            c[R] = r;
            calc_cx();

            if (!c[F].cmp(last_val[F]) && !c[R].cmp(last_val[R]))
                break;
            last_val[F] = c[F];
            last_val[R] = c[R];
        }

        calc_bound_concs();
    }

    void Primeanneal::read_primers_individual(std::string filename){
        FILE *infile = fopen(filename.c_str(), "r");
        char buffer[50];
        primers.clear();
        while(fscanf(infile, "%[ACGT]%*[^\n]\n", buffer) != EOF){
            //printf("%s\n", buffer);
            int len = strlen(buffer)+1;
            if(len != 21)
                printf("%s\n", buffer);
            Primeanneal::primer_info primer_pair;
            primer_pair.f = (char *) malloc(sizeof(char) * len);
            primer_pair.rc = (char *) malloc(sizeof(char) * len);
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

    void Primeanneal::swap_primer_rc(Primeanneal::primer_info *p){
        char *tmp = p->f;
        p->f = p->rc;
        p->rc = tmp;
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

    void Primeanneal::eval_thread(const char * out_filename, double dna_conc, double primer_conc, double mv_conc, double dv_conc, double dntp_conc){
        EQ eq;
        eq.c0[F] = primer_conc;
        eq.c0[R] = primer_conc;
        //eq.c0[A].set_d(dna_conc / addresses.size());
        //eq.c0[B].set_d(dna_conc / addresses.size());
        thal_results o;
        thal_args ta;
        set_thal_default_args(&ta);
        ta.mv = mv_conc;
        ta.dv = dv_conc;
        ta.dntp = dntp_conc;
        double f_hp_dh;
        double f_hp_ds;
        double r_hp_dh;
        double r_hp_ds;
        struct dh_ds_info{
            double f_f_dh;
            double f_f_ds;
            double f_r_dh;
            double f_r_ds;
            double f_frc_dh;
            double f_frc_ds;
            double f_rrc_dh;
            double f_rrc_ds;
            double r_f_dh;
            double r_f_ds;
            double r_r_dh;
            double r_r_ds;
            double r_frc_dh;
            double r_frc_ds;
            double r_rrc_dh;
            double r_rrc_ds;
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
                    eq.k[K_FX] = eq.k[K_FX] + dg_to_eq_const(dh_ds[j].f_f_dh - (dh_ds[j].f_f_ds * (temp_c + 273.15)), temp_c);
                    eq.k[K_RX] = eq.k[K_RX] + dg_to_eq_const(dh_ds[j].r_f_dh - (dh_ds[j].r_f_ds * (temp_c + 273.15)), temp_c);
                    count++;

                    eq.k[K_FX] = eq.k[K_FX] + dg_to_eq_const(dh_ds[j].f_frc_dh - (dh_ds[j].f_frc_ds * (temp_c + 273.15)), temp_c);
                    eq.k[K_RX] = eq.k[K_RX] + dg_to_eq_const(dh_ds[j].r_frc_dh - (dh_ds[j].r_frc_ds * (temp_c + 273.15)), temp_c);
                    count++;

                    eq.k[K_FX] = eq.k[K_FX] + dg_to_eq_const(dh_ds[j].f_r_dh - (dh_ds[j].f_r_ds * (temp_c + 273.15)), temp_c);
                    eq.k[K_RX] = eq.k[K_RX] + dg_to_eq_const(dh_ds[j].r_r_dh - (dh_ds[j].r_r_ds * (temp_c + 273.15)), temp_c);
                    count++;

                    eq.k[K_FX] = eq.k[K_FX] + dg_to_eq_const(dh_ds[j].f_rrc_dh - (dh_ds[j].f_rrc_ds * (temp_c + 273.15)), temp_c);
                    eq.k[K_RX] = eq.k[K_RX] + dg_to_eq_const(dh_ds[j].r_rrc_dh - (dh_ds[j].r_rrc_ds * (temp_c + 273.15)), temp_c);
                    count++;

                }
                eq.k[K_FX] = eq.k[K_FX] / count;
                eq.k[K_RX] = eq.k[K_RX] / count;
                eq.solve_eq();
                //c[F] and c[R] are now solved, all nonspec concs can now be solved individually

                //Set c[X] to concs of individual strands
                eq.c0[X].set_d(dna_conc / addresses.size());
                for (unsigned int j = 0; j < addresses.size(); j++){
                    if (i == j)
                        continue;
                    eq.k[K_FX].set_d(dg_to_eq_const(dh_ds[j].f_frc_dh - (dh_ds[j].f_frc_ds * (temp_c + 273.15)), temp_c));
                    eq.k[K_RX].set_d(dg_to_eq_const(dh_ds[j].r_frc_dh - (dh_ds[j].r_frc_ds * (temp_c + 273.15)), temp_c));
                    eq.calc_cx();
                    eq.calc_bound_concs();

                    //nonspec_exp_amp += min(fwd bound, rev bound)
                    eq.tmp[0].add(eq.c[FX], eq.c[RX]);
                    eq.tmp[2].min(eq.tmp[0], eq.tmp[1]);
                    eq.nonspec_exp_amp.add(eq.nonspec_exp_amp, eq.tmp[2]);

                    //nonspec_lin_amp += max(fwd bound, rev bound) - min(fwd bound, rev bound)
                    eq.tmp[0].max(eq.tmp[0], eq.tmp[1]);
                    eq.tmp[2].sub(eq.tmp[0], eq.tmp[2]);
                    eq.nonspec_lin_amp.add(eq.nonspec_lin_amp, eq.tmp[2]);
                }
                eq.tmp[0].div(eq.spec_fwd_amp, eq.nonspec_exp_amp);
                eq.tmp[1].div(eq.best_spec_exp_amp, eq.best_nonspec_exp_amp);
                if(first || (eq.tmp[0].cmp(eq.tmp[1]) > 0)){
                    first = false;
                    best_temp = temp_c;
                    eq.best_spec_exp_amp.set(eq.spec_fwd_amp);
                    eq.best_spec_lin_amp.set(eq.spec_rev_amp);
                    eq.best_nonspec_exp_amp.set(eq.nonspec_exp_amp);
                    eq.best_nonspec_lin_amp.set(eq.nonspec_lin_amp);
                }
            }
            eq.tmp[0].div(eq.best_spec_exp_amp, eq.best_nonspec_exp_amp);
            eq.tmp[1].div(eq.best_spec_lin_amp, eq.best_nonspec_lin_amp);
            outfile_mtx.lock();
            FILE *outfile = fopen(out_filename, "a");
            mpfr_fprintf(outfile, "%s,%s,%lf,%u,%.9Re,%.9Re,%.9Re,%.9Re,%.9Re,%.9Re\n",addresses[i].f, addresses[i].r, best_temp, i, eq.tmp[0].val, eq.tmp[1].val, eq.best_spec_exp_amp.val, eq.best_nonspec_exp_amp.val, eq.best_spec_lin_amp.val, eq.best_nonspec_lin_amp.val);
            fclose(outfile);
            outfile_mtx.unlock();
        }
    }

    void Primeanneal::evaluate_addresses(const char *in_filename, const char * out_filename, double dna_conc, double primer_conc, double mv_conc, double dv_conc, double dntp_conc){
        if(addresses.size() == 0)
            read_addresses(in_filename, dna_conc);
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

