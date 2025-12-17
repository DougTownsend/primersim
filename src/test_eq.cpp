#include <stdio.h>
#include <iostream>
#include <mpfr.h>
#include "eq.hpp"

void random_address_assignment(primersim::Primeanneal &pa, int count);

int main(void){
    primersim::Primeanneal pa;
    pa.num_cpu = 32;
    //double temp_c = 60;
    double dna_conc = 3e-15;
    double primer_conc = 250e-9;
    double mv = 100;
    double dv = 1.5;
    double dntp = 0.2;

    //pa.read_primers_individual("./codes.csv");
    //pa.assign_addresses("addresses_new.csv", 60, mv, dv, dntp);
    //pa.evaluate_addresses("addresses.csv", "sort_eval.csv", dna_conc, primer_conc, mv, dv, dntp);
    std::vector<double> temps;
    double temp = 55;
    for(int i = 0; i < 30; i++){
        temps.push_back(temp);
        /*
        temp -= 1.;
        if (temp < 40)
            temp = 40;
            */
    }
    //pa.read_primers_individual("codes.csv");
    //pa.assign_addresses_nosort("addresses.csv");
    pa.read_addresses("addresses.csv", false);
    std::cout << pa.addresses.size() << "\n";
    //random_address_assignment(pa, 1);

    pa.sim_pcr("test_bad3.csv", 0, temps.size(), temps, dna_conc, primer_conc, primer_conc, mv, dv, dntp);
    //printf("%03u: %e\n", i, ratio);
    /*
    for (double power = -20; power < -5.5; power += 1.){
        char filename[80];
        sprintf(filename, "sort_eval_new_250e%d.csv", (int) power);
        primer_conc = 250 *pow(10, power);
        pa.sim_pcr("sort_eval_new.csv", filename, 0, 30, dna_conc, primer_conc, primer_conc, mv, dv, dntp);
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
            pa.eval_addresses(out_filename, temps.size(), temps, dna_conc, primer_conc, primer_conc, mv, dv, dntp);
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