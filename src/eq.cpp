#include <stdio.h>
#include <string>
#include <mpfr.h>

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

}
