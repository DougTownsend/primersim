#include <stdio.h>
#include <string>
#include <cmath>

#include "eq.hpp"

namespace primersim{

    void EQ::print_state(std::string out_filename, std::string s){
        FILE *outfile = fopen(out_filename.c_str(), "a");
        fprintf(outfile, "%s,", s.c_str());
        for(int i = 0; i < 2; i++)
            fprintf(outfile, "last_val[%d],%.9Le,", i, (long double)last_val[i]);
        for(int i = 0; i < 10; i++)
            fprintf(outfile, "c[%d],%.9Le,", i, (long double)c[i]);
        for(int i = 0; i < 3; i++)
            fprintf(outfile, "c0[%d],%.9Le,", i, (long double)c0[i]);
        for(int i = 0; i < 7; i++)
            fprintf(outfile, "k[%d],%.9Le,", i, (long double)k[i]);
        fprintf(outfile, "\n");
        fclose(outfile);
    }

    void EQ::calc_cx(){
        // c[X] = c0[X] / (1 + k_FX*c[F] + k_RX*c[R])
        c[X] = c0[X] / (1.0 + k[K_FX]*c[F] + k[K_RX]*c[R]);
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
        constexpr int max_iter_outer = 64;
        constexpr int max_iter = 15;

        Real B1, B2;
        Real f, r;
        Real F1, F2;
        Real J11, J12, J21, J22;
        Real det, df, dr;
        Real f_new, r_new;
        Real residual, prev_res;

        // First call: cold-start in the middle of the box.
        // Subsequent calls: keep c[F], c[R] from the prior solve as a
        // warm start. c0/k change slowly between calls so the prior
        // equilibrium is usually within a couple of Newton iters of the
        // new one. calc_cx is always called to refresh c[X] under
        // current c0[X], k_FX, k_RX.
        if (!warm_start_valid) {
            c[F] = c0[F] / 2.0;
            c[R] = c0[R] / 2.0;
            warm_start_valid = true;
        }
        calc_cx();
        last_val[F] = c[F];
        last_val[R] = c[R];

        for (int outer = 0; outer < max_iter_outer; ++outer) {
            prev_res = std::numeric_limits<Real>::infinity();

            B1 = k[K_RX]*c[X] + k[K_RH] + 1.0;
            B2 = k[K_FX]*c[X] + k[K_FH] + 1.0;

            // Decoupled initial guesses (exact when k_FR = 0)
            f = (-B2 + std::sqrt(B2*B2 + 8.0*k[K_FF]*c0[F])) / (4.0*k[K_FF]);
            r = (-B1 + std::sqrt(B1*B1 + 8.0*k[K_RR]*c0[R])) / (4.0*k[K_RR]);

            for (int iter = 0; iter < max_iter; ++iter) {
                F1 = 2.0*k[K_RR]*r*r + (B1 + k[K_FR]*f)*r - c0[R];
                F2 = 2.0*k[K_FF]*f*f + (B2 + k[K_FR]*r)*f - c0[F];

                // Stop when residual stops decreasing — that's the precision
                // floor of Real arithmetic.
                residual = std::fabs(F1) + std::fabs(F2);
                if (residual >= prev_res) break;
                prev_res = residual;

                J11 = k[K_FR]*r;
                J12 = 4.0*k[K_RR]*r + B1 + k[K_FR]*f;
                J21 = 4.0*k[K_FF]*f + B2 + k[K_FR]*r;
                J22 = k[K_FR]*f;

                det = J11*J22 - J12*J21;
                df  = (J12*F2 - J22*F1) / det;
                dr  = (J21*F1 - J11*F2) / det;

                f_new = f + df;
                r_new = r + dr;

                if      (f_new < 0.0)    f_new = f * 0.5;
                else if (f_new > c0[F])   f_new = (f + c0[F]) * 0.5;
                if      (r_new < 0.0)    r_new = r * 0.5;
                else if (r_new > c0[R])   r_new = (r + c0[R]) * 0.5;

                f = f_new;
                r = r_new;
            }

            c[F] = f;
            c[R] = r;
            calc_cx();

            // Stop when c[F] and c[R] both stop moving by more than a few
            // ulps. 16*epsilon gives a relative tolerance of ~3.6e-15 for
            // double, ~1.7e-18 for long double — well below the ~5e-10
            // precision floor of the %.9 output format.
            constexpr Real conv_eps = 16 * std::numeric_limits<Real>::epsilon();
            if (std::fabs(c[F] - last_val[F]) <= conv_eps * std::fabs(c[F]) &&
                std::fabs(c[R] - last_val[R]) <= conv_eps * std::fabs(c[R]))
                break;
            last_val[F] = c[F];
            last_val[R] = c[R];
        }

        calc_bound_concs();
    }

}
