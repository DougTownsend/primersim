#include <stdio.h>
#include <iostream>
#include <random>
#include <string>
#include <cstdlib>
#include "eq.hpp"

void random_address_assignment(primersim::Primeanneal &pa, int count);

// strip dirname + ".csv" suffix to get a basename usable in filenames.
//   "primers_100p.csv"          → "primers_100p"
//   "/path/to/primers_1958.csv" → "primers_1958"
static std::string basename_no_ext(const std::string &path) {
    std::string s = path;
    auto slash = s.find_last_of('/');
    if (slash != std::string::npos) s = s.substr(slash + 1);
    auto dot = s.find_last_of('.');
    if (dot != std::string::npos) s = s.substr(0, dot);
    return s;
}

int main(int argc, char **argv){
    // Usage:
    //   ./test_eq                              regression, num_cpu=24
    //   ./test_eq <num_cpu>                    regression with custom thread cap
    //   ./test_eq <num_cpu> sweep <n_trials> [out_dir]
    //
    // num_cpu is the max worker thread count. The HPC scheduler does
    // not pin thread counts — the binary self-caps to whatever was
    // requested via -n in the bsub script (e.g. $LSB_DJOB_NUMPROC).
    char regression_fname[] = "regression_new.csv";
    primersim::Primeanneal pa;
    pa.num_cpu = (argc > 1) ? std::atoi(argv[1]) : 24;
    if (pa.num_cpu < 1) {
        std::cerr << "usage: " << argv[0]
                  << " <num_cpu> [sweep <n_trials> [out_dir]]\n";
        return 1;
    }
    double dna_conc = 3e-15;
    double primer_conc = 250e-9;
    double mv = 100;
    double dv = 1.5;
    double dntp = 0.2;
    // Input file paths are env-var configurable so different primer
    // sets can be swapped without renaming. Defaults match the legacy
    // hardcoded names. The basename (without dir and .csv) is also
    // used to tag sweep output filenames so it's clear which primer
    // set produced which sweep.
    const char *primers_path  = std::getenv("PRIMERS_FILE");
    if (!primers_path)  primers_path  = "primers.csv";
    const char *pairings_path = std::getenv("PAIRINGS_FILE");
    if (!pairings_path) pairings_path = "pairings.csv";
    std::string primers_basename = basename_no_ext(primers_path);

    pa.read_primer_pool(primers_path);
    pa.read_pairings(pairings_path);
    std::cout << pa.primer_pool.size() << " primers (" << primers_path << "), "
              << pa.pairings.size() << " pairings (" << pairings_path << "), "
              << "num_cpu=" << pa.num_cpu << "\n";
    // Precompute every (primer, kind, primer, kind) thal pair once.
    // The cache is pairing-agnostic: shuffling pairings (e.g. for
    // exploring assignments) does NOT invalidate it. Uses pa.num_cpu
    // threads internally — the dominant work in this benchmark.
    pa.populate_dimer_cache(mv, dv, dntp, /*temp_c=*/55.0);

    if (argc > 2 && std::string(argv[2]) == "sweep") {
        // CLI: ./test_eq <num_cpu> sweep <n_trials> [out_dir] [max_data_gb]
        // n_trials = 0 means loop until externally killed (e.g. LSF -W).
        // max_data_gb = 0 (default) means unlimited; otherwise the
        // sweep stops once the rolling output files exceed the budget.
        int n_trials = (argc > 3) ? std::atoi(argv[3]) : 10;
        const char *out_dir = (argc > 4) ? argv[4] : "sweep_out";
        double max_gb = (argc > 5) ? std::atof(argv[5]) : 0.0;
        std::vector<int> temps_c;
        for (int t = 50; t <= 70; t++) temps_c.push_back(t);
        std::cout << "sweep: "
                  << (n_trials == 0 ? "infinite" : std::to_string(n_trials))
                  << " trials × " << temps_c.size() << " temps -> " << out_dir
                  << " (max " << (max_gb > 0 ? std::to_string(max_gb) + " GB" : "unlimited")
                  << ")\n";
        pa.sweep_pairings(out_dir, primers_basename.c_str(), n_trials, max_gb,
                          temps_c, /*pcr_cycles=*/30,
                          dna_conc, primer_conc, primer_conc, mv, dv, dntp);
        return 0;
    }

    FILE *f_tmp = fopen(regression_fname, "w"); //clear the regression file
    fclose(f_tmp);
    std::vector<double> temps;
    double temp = 55;
    for(int i = 0; i < 30; i++){
        temps.push_back(temp);
    }
    for (size_t i = 0; i < pa.pairings.size(); i++){
        FILE *regfile = fopen(regression_fname, "a");
        fprintf(regfile, "\nSimulating address %lu\n", i);
        fclose(regfile);
        pa.sim_pcr(pa.pairings, "regression_new.csv", i, temps.size(), temps, dna_conc, primer_conc, primer_conc, mv, dv, dntp, /*log_cycles=*/true);
    }
    //printf("%03u: %e\n", i, ratio);
    /*
    for (double power = -20; power < -5.5; power += 1.){
        char filename[80];
        sprintf(filename, "sort_eval_new_250e%d.csv", (int) power);
        primer_conc = 250 *pow(10, power);
        pa.sim_pcr(pa.pairings, "sort_eval_new.csv", filename, 0, 30, dna_conc, primer_conc, primer_conc, mv, dv, dntp, true);
    }
    */
    return 0;
}

void random_address_assignment(primersim::Primeanneal &pa, int count){
    double dna_conc = 3e-15;
    double primer_conc = 250e-9;
    double mv = 100;
    double dv = 1.5;
    double dntp = 0.2;
    std::vector<double> temps;
    double temp = 55;
    for(int i = 0; i < 30; i++){
        temps.push_back(temp);
    }
    printf("reading addresses\n");
    fflush(stdout);
    pa.read_addresses("addresses.csv", false);
    printf("done\n");
    fflush(stdout);
    for(int i = 0; i < count; i++){
        char out_filename[100];
        //pa.shuffle_addresses();
        sprintf(out_filename, "rand_addresses_%04d.csv", i);
        FILE *f = fopen(out_filename, "w");
        for(auto a : pa.addresses){
            fprintf(f, "%s,%s\n", a.f, a.r);
        }
        fclose(f);
        printf("%s\n", pa.addresses[0].f);
        fflush(stdout);
        for(double temp_c = 54; temp_c < 60.1; temp_c += 1){
            for(unsigned int j = 0; j < temps.size(); j++){
                temps[j] = temp_c;
            }
            if(temp_c < 54.5){
                double touchdown_temp = 60.0;
                for(unsigned int j = 0; j < temps.size(); j++){
                    temps[j] = touchdown_temp;
                    if(touchdown_temp > 55.5)
                        touchdown_temp -= 1.0;
                }
            }
            sprintf(out_filename, "rand_addr_eval_%04d_%dC.csv", i, (int) temp_c);
            pa.eval_addresses(out_filename, NULL, temps.size(), temps, dna_conc, primer_conc, primer_conc, mv, dv, dntp);
            FILE *infile = fopen(out_filename, "r");
            double min_amp = INFINITY;
            double amp;
            while(fscanf(infile, "%*d,%le\n", &amp) != EOF){
                if (amp < min_amp)
                    min_amp = amp;
            }
            fclose(infile);
            sprintf(out_filename, "shuffled_min_amps_%dC.csv", (int) temp_c);
            FILE *outfile = fopen(out_filename, "a");
            fprintf(outfile, "%e\n", min_amp);
            fclose(outfile);
        }
    }

}