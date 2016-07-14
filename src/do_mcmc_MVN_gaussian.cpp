// *********************************************************************
// Copyright (C) 2016 James T. MacDonald
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
// *********************************************************************

#include "mcmc_multispecies_utils.h"
#include "model.h"

#include <iostream>
#include <fstream>


int main(int argc, char **argv){

    uint64_t seed = 315875;


	bmeg_functor::para_type p = bmeg_functor::para_type();

	const int mat_GFP_spec_index = bmeg_functor::mat_GFP_spec_index;
	const int MRNA_spec_index = bmeg_functor::MRNA_spec_index;
	uint_vec species_vec(0);
	species_vec.push_back(mat_GFP_spec_index);
    species_vec.push_back(MRNA_spec_index);


    const double mRNA_scale_nM = 1.0/1000.0;

    double min_abs_GFP_stddev = 0.1;
    double min_abs_mRNA_stddev = 10.0 * mRNA_scale_nM;

    string prior_sigma_file = "uninf";
    string prior_mean_file = "uninf";
    bool using_mvn_prior = false;



    seed = lexical_cast<uint64_t>(argv[1]);
    string sigma_file = string(argv[2]);
    cout << "# seed:\t" << seed << endl;

	sim_condition_map conds_map;
	double_vec times;

	string GFP_list_file = string(argv[3]);

    string mRNA_list_file = string(argv[4]);

	string param_file = string(argv[5]);


	if (argc >= 7){
        cout << "# reading min_abs_GFP_stddev from command line: " << argv[6] << endl;
	    min_abs_GFP_stddev = lexical_cast<double>(argv[6]);
	}
	if (argc >= 8){
	    cout << "# reading min_abs_mRNA_stddev from command line: " << argv[7] << endl;
        min_abs_mRNA_stddev = lexical_cast<double>(argv[7]) * mRNA_scale_nM;
	}
	if (argc >= 10){
	    cout << "# reading MVN prior from files: mean: " << argv[8] << " sigma: " << argv[9] << endl;
        using_mvn_prior = true;
        prior_mean_file = argv[8];
        prior_sigma_file = argv[9];
	}






	cout << "# parsing GFP list file" << endl;
	parse_MVN_list_file(GFP_list_file, conds_map, mat_GFP_spec_index, times, 1.0, min_abs_GFP_stddev*min_abs_GFP_stddev);


    cout << "# parsing mRNA list file" << endl;
    parse_MVN_list_file(mRNA_list_file, conds_map, MRNA_spec_index, times, mRNA_scale_nM, min_abs_mRNA_stddev*min_abs_mRNA_stddev);


    cout << "# number of time steps: " << times.size() << endl;

    rand_gen.seed(seed);
	// burn in random num generator
	for (int ii = 0; ii < 1000; ii++){
		rand_gen.randInt();
	}


	cout << "# Parsing start parameters" << endl;
	parse_initial_params(param_file, p);


    cout << "# Loading proposal distribution Sigma" << endl;

	move_parameters_functor move_functor = move_parameters_functor();
	move_functor.loadSigma(sigma_file, seed+1);


	unsigned long steps = 100000000;
	unsigned long burn_in = 1000; //5000
	unsigned long output_freq = 10;  //10

	if (using_mvn_prior == false){
        double ll = get_LL_condition_set(bmeg_functor(p), init_conds_functor<bmeg_functor>()(bmeg_functor(p)), species_vec, times, conds_map);
        cout << "# Initial LL: " << ll << endl;
        do_mcmc(bmeg_functor(p), init_conds_functor<bmeg_functor>(), species_vec, times, conds_map, steps, burn_in, output_freq,
            uninf_prior_functor<bmeg_functor::para_type>(), move_functor);
    }
    else {
        cout << "# parsing MVN prior files" << endl;
        Eigen::Matrix<double, 9, 1> prior_mean = readMatrix<double, 9, 1>(prior_mean_file);
        Eigen::Matrix<double, 9, 9> prior_sigma = readMatrix<double, 9, 9>(prior_sigma_file);
        multigauss_prior_functor<bmeg_functor::para_type, 9> prior_functor;
        uint_vec mapping(9);

        mapping[0] = 1; // tl_rate,
        mapping[1] = 2; // tx_rate,
        mapping[2] = 3; // rna_deg_rate,
        mapping[3] = 8; // mat_rate
        mapping[4] = 9; // init_mrna_cap,
        mapping[5] = 10; // mrna_cap_deg_rate,
        mapping[6] = 11; //km_mrna_cap,
        mapping[7] = 12; //km_dna,
        mapping[8] = 17; //km_mrna,
        prior_functor.set_covar(prior_sigma);
        prior_functor.set_mean(prior_mean);
        prior_functor.set_mapping(mapping);

        double ll = get_LL_condition_set(bmeg_functor(p), init_conds_functor<bmeg_functor>()(bmeg_functor(p)), species_vec, times, conds_map);
        cout << "# Initial LL: " << ll << endl;
        cout << "# Initial prior: " << prior_functor(p) << endl;
        do_mcmc(bmeg_functor(p), init_conds_functor<bmeg_functor>(), species_vec, times, conds_map, steps, burn_in, output_freq,
            prior_functor, move_functor);


    }

}




