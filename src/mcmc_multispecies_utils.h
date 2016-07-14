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

#ifndef MCMC_H_
#define MCMC_H_


#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <cstdio>
#include <sstream>

#include<boost/algorithm/string.hpp>
#include<boost/lexical_cast.hpp>
#include<boost/algorithm/string/split.hpp>

#include <boost/array.hpp>

#include <boost/numeric/odeint.hpp>

#include <boost/math/distributions/normal.hpp>
#include  <cmath>

#include "MersenneTwister.h"

#include "multivariate_normal.h"

using namespace std;
using namespace boost::numeric::odeint;
using boost::math::normal_distribution;
using boost::math::normal;
using boost::math::pdf;
using std::log;
using std::exp;

using std::cout;
using std::endl;
using std::cerr;
using std::ios;
using std::ifstream;
using std::ofstream;
using std::istream;
using std::strcpy;
using std::strcat;
using std::strcmp;
using std::pow;
using std::sqrt;
using std::asin;
using std::string;

using boost::lexical_cast;
using boost::split;
using boost::is_any_of;
using boost::trim;



typedef std::vector<double> double_vec;
typedef std::vector<string> string_vector;
typedef std::vector<unsigned int> uint_vec;
typedef std::vector<std::vector<double> > double_vec_vec;
typedef std::map< unsigned int, std::vector<double> > double_vec_map;
typedef std::map< unsigned int, MultivariateNormalPDF > MultivariateNormalPDF_map;


MTRand rand_gen;




template<class state_type_>
class observer_whole_traj_functor {
public:
    double_vec& store_times;
    vector<state_type_>& store_traj;

public:
    observer_whole_traj_functor(double_vec& store_times_, vector<state_type_>& store_traj_) : store_times(store_times_), store_traj(store_traj_){
    }

    void operator()(const state_type_ &x , const double t){
        store_traj.push_back(x);
        store_times.push_back(t);
    }

};




template<class state_type_>
void write_traj( const state_type_ &x , const double t )
{
    cout << t << '\t' << x[0] << '\t' << x[1] << '\t' << x[2]  << "\t" << x[3] << endl;
}



template <typename state_type_>
double get_LL(vector<state_type_>& sim_traj, const uint_vec species,
              MultivariateNormalPDF_map& mvnpdf_map){


        double sum_ll = 0;

        for (unsigned int iter = 0; iter != species.size(); iter++){
                const unsigned int this_species = species[iter];



        if (mvnpdf_map.count(this_species) > 0 )
        {
            Eigen::Matrix<double,Eigen::Dynamic,1> sim_traj_mat(sim_traj.size());
            for (unsigned long i = 0; i < sim_traj.size(); i++)
            {
                sim_traj_mat(i) = sim_traj[i][this_species];
            }
            const double this_ll = mvnpdf_map[this_species].getLogPdf(sim_traj_mat);

            sum_ll += this_ll;
        }



        }

        return sum_ll;
}


template <typename system_>
vector<typename system_::state_type> run_trajectory(system_ sys, typename system_::state_type init_conds, double_vec times){

        typedef runge_kutta_fehlberg78< typename system_::state_type > error_stepper_type;
        typedef controlled_runge_kutta< error_stepper_type > controlled_stepper_type;
        controlled_stepper_type controlled_stepper;

        double_vec store_times;
        vector< typename system_::state_type >  store_traj;

        observer_whole_traj_functor< typename system_::state_type > obs = observer_whole_traj_functor< typename system_::state_type >(store_times, store_traj);
        integrate_times( controlled_stepper , sys, init_conds , times.begin() , times.end() , 0.1 , obs );

	return store_traj;
}

template <typename system_>
double get_loglikelihood(system_ sys, typename system_::state_type init_conds, const uint_vec species,
                         double_vec times, MultivariateNormalPDF_map& mvnpdf_map ){

	vector<typename system_::state_type> store_traj = run_trajectory(sys, init_conds, times);

	return get_LL(store_traj, species,  mvnpdf_map);

}



double get_mean(const double_vec& vec){
	double sum = 0;
	double count = 0;
	for (double_vec::const_iterator it = vec.begin(); it != vec.end(); it++){
		sum += *it;
		count += 1;
	}
	return sum/count;
}


class sim_condition{

public:
	// note all conc units in uM
	double xylr_conc;
	double xylose_conc;
	double dna_conc;

	double_vec_map exp_traj;
	double_vec_map exp_traj_stddev;

	MultivariateNormalPDF_map mvn_pdf_map;

	sim_condition(double xylr_conc_, double xylose_conc_, double dna_conc_){
		exp_traj = double_vec_map();
		exp_traj_stddev = double_vec_map();
		xylr_conc = xylr_conc_;
		xylose_conc = xylose_conc_;
		dna_conc = dna_conc_;
	}
	sim_condition(){
		exp_traj = double_vec_map();
		exp_traj_stddev = double_vec_map();
	}

	void print_summary() const{
		cout << "xylr_conc:" << xylr_conc << "\t"
		     << "xylose_conc:" << xylose_conc << "\t"
		     << "dna_conc:" << dna_conc << "\t"
		     << "exp_traj.size():" << exp_traj.size() << "\t"
		     //<< "get_mean(exp_traj):" << get_mean(exp_traj) << "\t"
		     << "exp_traj_stddev.size():" << exp_traj_stddev.size() << "\t";
		     //<< "get_mean(exp_traj_stddev):" << get_mean(exp_traj_stddev) << "\t";
//		     << endl;

	}

};

typedef std::vector<sim_condition> sim_condition_vec;
typedef std::map<string, sim_condition> sim_condition_map;






string make_new_key(double dna_conc, double xylose_conc, double xylr_conc){
	std::stringstream skey(std::stringstream::in | std::stringstream::out);
	skey << "d"  << setfill('0') << setw(5) << (dna_conc*1000000) << "_x" << setfill('0') << setw(5) << (xylose_conc) << "_r" << setfill('0') << setw(5) << (xylr_conc*1000);
	return skey.str();
}

bool readMeans(string filename, double_vec& times, double_vec& means){
    times.clear();
    means.clear();
    ifstream input(filename.c_str(), ios::in);
    string_vector SplitVec;
    string lineStr;
    unsigned long length, lineNum = 0;

    while ( !input.eof() )
    {
        getline(input, lineStr);
        string resStr;
        lineNum++;

        length = lineStr.length();


        if (length > 1)
        {
            split( SplitVec, lineStr, is_any_of("\t") );
            const double time = lexical_cast<double>(SplitVec[0]);
            const double val = lexical_cast<double>(SplitVec[1]);
            times.push_back(time);
            means.push_back(val);
        }
    }

    return true;
}

bool parse_MVN_list_file(string filename, sim_condition_map& conds, const unsigned int species, double_vec& times, const double scale = 1.0, const double min_variance = 0.0)
{
   //const bool use_new_key = true;
    //times.clear();
    ifstream input(filename.c_str(), ios::in);
    string_vector SplitVec;
    string lineStr;
    unsigned long length, lineNum = 0;

    while ( !input.eof() )
    {
        getline(input, lineStr);
        string resStr;
        lineNum++;

        length = lineStr.length();


        if (length > 0)
        {
            split( SplitVec, lineStr, is_any_of("\t") );

            // parsing trajectories

            double dna_conc = lexical_cast<double>(SplitVec[0])  /1000.0;
            double xylose_conc = lexical_cast<double>(SplitVec[1]) * 1000.0;
            double xylr_conc = lexical_cast<double>(SplitVec[2]);
            string sigma_filename = SplitVec[3];
            string means_filename = SplitVec[4];
            string key = make_new_key(dna_conc, xylose_conc, xylr_conc);


            if (conds.find(key) == conds.end())
            {
                // create new element
                cout << "# creating new condition: " << key  << endl;
                conds[key] = sim_condition(xylr_conc, xylose_conc, dna_conc);
            }

            cout << "# reading LL Sigma: " << sigma_filename << endl;
            Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> cov_mat = (scale*scale) * readDynamicMatrix(sigma_filename);
            for (int i = 0; i < cov_mat.rows(); i++){
                if (cov_mat(i,i) < min_variance){
                    cov_mat(i,i) = min_variance;
                }
            }
            double_vec means_vec;
            cout << "# reading LL Means: " << means_filename << endl;
            readMeans(means_filename, times, means_vec);
            Eigen::Matrix<double,Eigen::Dynamic,1> means_mat(means_vec.size(), 1);
            for (unsigned int i = 0; i< means_vec.size(); i++){
                means_mat(i) = scale * means_vec[i];
            }

            if (means_mat.rows() != cov_mat.rows()){
                cerr << "ERROR: Dimensions of means and covs do not match: " << means_mat.rows() << "\t" << cov_mat.rows() << endl;
                return false;
            }

            conds[key].mvn_pdf_map[species].setMean(means_mat);
            conds[key].mvn_pdf_map[species].setCovar(cov_mat);


        }
    }

    return true;
}



template <typename para_type_>
bool parse_initial_params(string filename, para_type_& p){
        ifstream input(filename.c_str(), ios::in);
        string_vector SplitVec;
        string lineStr;
        unsigned long length, lineNum = -1;

        while ( !input.eof() ) {
                getline(input, lineStr);
		trim(lineStr);
                string resStr;
                lineNum++;

                length = lineStr.length();


                if (length > 0) {
                        split( SplitVec, lineStr, is_any_of("\t ") );
                        p[lineNum]= lexical_cast<double>(SplitVec[0]);
                }
        }
    return true;
}



template <typename system_>
double get_LL_condition_set(system_ sys, typename system_::state_type init_conds, const uint_vec species, double_vec& times, sim_condition_map& conds ){
	sim_condition_map::iterator iter;
	double ll_sum = 0;
	for (iter = conds.begin(); iter != conds.end(); iter++){
		sys.p[13] = iter->second.dna_conc; // dna
		sys.p[14] = iter->second.xylr_conc; // xylr
		sys.p[15] = iter->second.xylose_conc; // xylose
		ll_sum += get_loglikelihood(sys, init_conds, species, times, iter->second.mvn_pdf_map);
	}
	return ll_sum;
}


double get_gaussian_sample(double value, double stddev, double low, double high){
	bool ok = false;
	double new_value = value;
	while (ok == false){
		new_value = rand_gen.randNorm(value, stddev);
		if (new_value >= low && new_value <= high) ok = true;
	}
	return new_value;
}



// uninformative prior function
template <typename para_type_>
class uninf_prior_functor{
public:
    double operator()(const para_type_&  p){
        return 0.0;
    }
};

// multivariate gaussian prior functor
template <typename para_type_, unsigned int DIM>
class multigauss_prior_functor{
private:
    uint_vec mapping;
    Eigen::Matrix<double,DIM,1> mapped_vec;
    MultivariateNormalPDF pdf;

public:


    void set_mapping(uint_vec mapping_){
        if (DIM != mapping_.size()){
            cout << "ERROR: multigauss_prior_functor.set_mapping()  - dimensions don't match" << endl;
        }
        assert(DIM == mapping_.size());
        mapping = mapping_;
    }

    void set_covar(Eigen::Matrix<double,DIM,DIM> covar_){
        pdf.setCovar(covar_);
    }

    void set_mean(Eigen::Matrix<double,DIM,1> mean_){
        pdf.setMean(mean_);
    }

    double operator()(const para_type_&  p){
        for (unsigned int ii = 0; ii < DIM; ii++){
            mapped_vec(ii) = p[mapping[ii]];
        }
        const double logpdf = pdf.getLogPdf(mapped_vec);
        return logpdf;
    }
};




template<typename system_, typename init_conds_functor_, typename prior_, typename mover_>
int do_mcmc(system_ sys, init_conds_functor_ init_functor, uint_vec species, double_vec& times, sim_condition_map& conds_map,
            unsigned long steps, unsigned long burn_in, unsigned long output_freq, prior_ log_prior, mover_ mover){


	double start_ll = get_LL_condition_set(sys, init_functor(sys), species, times, conds_map) + log_prior(sys.p);
	double acc_ll = start_ll; //last accepted ll
	double this_ll = acc_ll;
	typename system_::para_type acc_p = sys.p; // last accepted p
	unsigned long step_num = 0;
	unsigned long num_trials = 0;
	const unsigned long output_acc_rate_freq = output_freq * 100;

	while (step_num < steps){
		sys.p = mover(acc_p);
		this_ll = get_LL_condition_set(sys, init_functor(sys), species, times, conds_map) + log_prior(sys.p);
		bool accept = false;
		if (!std::isinf(this_ll)){
			if (this_ll > acc_ll){
				accept = true;
			}
			else {
				double p_acc = exp(this_ll - acc_ll);
				double rand = rand_gen.rand();
				if (rand < p_acc){
					accept = true;
				}
				else {
				}
			}
		}
		if (accept){
			if (step_num >= burn_in){
				if (step_num%output_freq == 0){
					cout << "ACCEPTED\tstep_num:\t" << step_num << "\tLL:\t" << this_ll;

					for (int ii = 0 ; ii < sys.p.size(); ii++){
						cout << "\t" << sys.p[ii];
					}
					cout << endl;
				}
			}



			acc_ll = this_ll;
			acc_p = sys.p;



			step_num++;
		}
		num_trials++;
		if (num_trials % output_acc_rate_freq == 0){
			const double acc_rate = static_cast<double>(step_num)/ static_cast<double>(num_trials);
			cout << "#trials:\t" << num_trials << "\tstep_num:\t" << step_num << "\tAcceptance_rate:\t" << acc_rate << endl;
		}
	}
	cout << "# MCMC finished\tnum_trials:\t" << num_trials << "\t\tsteps\t" << steps << endl;
	return true;
}


template <typename para_type_>
vector<para_type_> parse_parameters_list(string filename)
{
    ifstream input(filename.c_str(), ios::in);
    string_vector SplitVec;
    string lineStr;
    unsigned long length, lineNum = -1;

    vector<para_type_> rtn_vec;

    while ( !input.eof() )
    {
        getline(input, lineStr);
        trim(lineStr);
        string resStr;
        lineNum++;

        length = lineStr.length();

        para_type_ p;
        if (length > 0)
        {
            split( SplitVec, lineStr, is_any_of("\t ") );
            if ( SplitVec[0].substr(0,1).compare("#") != 0 )
            {
                para_type_ p;
                for (int i = 5; i < (5 + p.size()); i++)
                {
                    double val = lexical_cast<double>(SplitVec[i]);
                    p[i-5]= val;
                }
                rtn_vec.push_back(p);
            }
        }
    }
    return rtn_vec;
}

template<typename system_, typename init_conds_functor_>
std::map<string, vector<typename system_::state_type> > get_traj_condition_set(system_ sys, init_conds_functor_ init_functor, double_vec& times, sim_condition_map& conds ){

	std::map<string, vector<typename system_::state_type> > rtn_vec;

        sim_condition_map::iterator iter;
        double ll_sum = 0;
        for (iter = conds.begin(); iter != conds.end(); iter++){
                sys.p[13] = iter->second.dna_conc; // dna
                sys.p[14] = iter->second.xylr_conc; // xylr
                sys.p[15] = iter->second.xylose_conc; // xylose
                vector<typename system_::state_type> traj = run_trajectory(sys, init_functor(sys), times);
		rtn_vec[iter->first] = traj;
        }
        return rtn_vec;
}

template <class state_type_, class para_type_>
void print_summary(para_type_  p, double_vec& times, sim_condition_map& conds, const double min_stddev, string prefix){
        sim_condition_map::const_iterator iter;
        for (iter = conds.begin(); iter != conds.end(); iter++){
		double val = get_loglikelihood<state_type_>(p, times, iter->second.exp_traj, iter->second.exp_traj_stddev, min_stddev);
                cout << prefix << "\t" << iter->first << "\t";
                iter->second.print_summary();
                cout << val << endl;
        }
}



#endif /* MCMC_H_ */



