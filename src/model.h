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

#ifndef BMEG_MODEL_H_
#define BMEG_MODEL_H_

#include "eigenmvn.h"

class bmeg_functor {

private:

    double tl_MM_rxn_new(double param_0, double modif_0, double param_1) 	//tl_MM_rxn_new
        {return  param_0*(modif_0/param_1/(1.00000000000000000+modif_0/param_1));}
    double tx_ind_rxn_new(double param_0, double param_1, double param_2, double param_3, double sub_0, double param_4, double param_5, double param_6, double param_7, double param_8, double param_9, double param_10) 	//tx_ind_rxn_new
        {
		// param_1 is DNA concentration
		const double hill_function = (1.00000000000000000/(1.00000000000000000+pow((param_5*(1.00000000000000000/(1.00000000000000000+pow((param_6/param_7),param_8)))/param_9),param_10)));

		return  param_0*(pow(((param_1 * hill_function)/param_2),param_3)*(sub_0/param_4)/(pow((1.00000000000000000+(param_1 * hill_function)/param_2),param_3)*(1.00000000000000000+sub_0/param_4)));

}



public:

    typedef boost::array< double , 19 > para_type;
    typedef boost::array< double , 4 > state_type;
    para_type p;

    bmeg_functor(){
        p = para_type();
    }

    bmeg_functor(para_type p_init){
        p = para_type();
        p = p_init;
    }

    static const int mat_GFP_spec_index = 3;
    static const int MRNA_spec_index = 1;

    void operator()( const state_type &x , state_type &dxdt , double t )
    {

        dxdt[0] = tl_MM_rxn_new(p[1], x[1], p[17])*p[0]-(p[8] * x[0]) *p[0];	//
        dxdt[1] = -p[3]*pow(max(x[1],0.0), p[18])*p[0] + tx_ind_rxn_new(p[2], p[13], p[12], p[16], x[2], p[11], p[14], p[15], p[5], p[7], p[4], p[6])*p[0];	//
        dxdt[2] = -tx_ind_rxn_new(p[2], p[13], p[12], p[16], x[2], p[11], p[14], p[15], p[5], p[7], p[4], p[6])*p[0]-(p[10] * x[2]) *p[0];	//
        dxdt[3] = (p[8] * x[0]) *p[0];	//


    }
};




class move_parameters_functor{
public:
        bmeg_functor::para_type stddev, low, high;
        Eigen::Matrix<double,Eigen::Dynamic,1> mean ;
        Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> covar;
        Eigen::EigenMultivariateNormal<double> norm_cholesk;

        move_parameters_functor() : mean(Eigen::Matrix<double, 11, 1>::Zero()), covar(Eigen::Matrix<double, 11, 11>::Zero()), norm_cholesk(mean, covar){
                const double c = 10.0;//0.5;
                stddev[0] =  c*0;               low[0]=1;               high[0]=1;       //compartment 'TestTube'
                stddev[1] =  c*0.0001;            low[1]=0.005;               high[1]=20; //17.6586; //global quantity 'tl_rate'                              7       ACF 0 around 400 lag - prob OK
                stddev[2] =  c*0.000001;        low[2]=0.0005;          high[2]=5.0;//0.00144542;      //global quantity 'tx_rate'                    8       ACF 0 around 250 lag - prob OK
                stddev[3] =  c*0.0001;          low[3]=0.01;            high[3]=0.4;//0.137127;        //global quantity 'rna_deg_rate'                 9       ACF 0 around 150 lag -
                stddev[4] =  c*0;          low[4]=0;               high[4]=1;//0.0750182;       //global quantity 'kd'                             10      ACF 0 around  lag        #reduce these
                stddev[5] =  c*0;             low[5]=0;               high[5]=500;//18.3364; //global quantity 'kx'                                   11      ACF 0 around  lag       #reduce these
                stddev[6] =  c*0;               low[6]=0.1;             high[6]=5;//4      //global quantity 'n'                                        12      N/A
                stddev[7] =  c*0;               low[7]=0.1;             high[7]=2;//1      //global quantity 'm'                                        13      N/A
                stddev[8] =  c*0.00004;          low[8]=0.01;            high[8]=0.2;//0.174011;        //global quantity 'mat_rate'                     14      ACF 0 around 100 lag - prob OK
                stddev[9] =  c*0.001;            low[9]=0.1;             high[9]=20.0;//1.32792; //global quantity 'init_mrna_cap'                       15      ACF 0 around 400 lag - prob OK
                stddev[10] = c*0.0001;         low[10]=0.001;          high[10]=0.25;//0.0191179;      //global quantity 'mrna_cap_deg_rate'           16      ACF 0 around 250 lag
                stddev[11] = c*0.0001;          low[11]=0.0001;              high[11]=5.0;//0.183292;       //global quantity 'km_mrna_cap'                  17      ACF 0 around 1000 lag - prob OK     # a problem parameter
                stddev[12] = c*0.00001;          low[12]=0.0001;         high[12]=0.1;//0.0340584;      //global quantity 'km_dna'                       18      ACF 0 around 400 lag - prob OK
                stddev[13] = c*0;               low[13]=0;              high[13]=1000; //global quantity 'dna_conc_d2_5'
                stddev[14] = c*0;               low[14]=0;              high[14]=1000;    //global quantity 'xylr_conc_r0_0'
                stddev[15] = c*0;               low[15]=0;              high[15]=1000;   //global quantity 'xylose_conc_x0_0'
                stddev[16] = c*0.0;          low[16]=1.0;              high[16]=1;   //global quantity 'd'                                               22      ACF  0 around 20 lag
                stddev[17] = c*0.0005;          low[17]=0.0;              high[17]=5;   //global quantity 'km_mrna'                                         23      ACF  0 around 500 lag
                stddev[18] = c*0.001;          low[18]=1.0;              high[18]=2;   //global quantity 'rna_deg_order'					24
        }

        bool check_bounds(bmeg_functor::para_type  p){
            for  (unsigned int ii = 0; ii < p.size(); ii++){
                if (p[ii] < low[ii] ){
                    //cout << "low\t" << ii << "\t" << p[ii] << endl;
                    return false;
                }
                else if (p[ii] > high[ii]){
                    //cout << "high\t" << ii << "\t" << p[ii] << endl;
                    return false;
                }
            }
            //cout << "OK" << endl;
            return true;
        }

        bmeg_functor::para_type operator()(bmeg_functor::para_type  orig_p){



            bmeg_functor::para_type  p = orig_p;
            bool bounds_ok = false;

            while (bounds_ok == false){
                p = move_MultiGauss(orig_p);
                bounds_ok = check_bounds(p);
            }


            return p;



        }

        bmeg_functor::para_type move_all(bmeg_functor::para_type  p){
                for (unsigned int ii = 0; ii < p.size(); ii++){
                        if (stddev[ii] != 0 ){
                                p[ii] = get_gaussian_sample(p[ii], stddev[ii], low[ii], high[ii]);
                        }
                }
                return p;
        }

        bmeg_functor::para_type move_MultiGauss(bmeg_functor::para_type  p){


            const double scale = 0.4;//0.1;

            Eigen::Matrix<double,Eigen::Dynamic,1> move_ = this->norm_cholesk.samples(1);
            //p[] = scale * move_(,0);
            p[1] += scale * move_(0,0);
            p[2] += scale * move_(1,0);
            p[3] += scale * move_(2,0);
            p[8] += scale * move_(3,0);
            p[9] += scale * move_(4,0);
            p[10] += scale * move_(5,0);
            p[11] += scale * move_(6,0);
            p[12] += scale * move_(7,0);
            p[17] += scale * move_(8,0);

            return p;
        }

        bmeg_functor::para_type move_one(bmeg_functor::para_type  p){

            const double scale = 1;

            bool found_suitable =false;
            int para = 0; //rand_gen.randInt(p.size());
            while (found_suitable == false){
                para = rand_gen.randInt(p.size()-1);
                if (stddev[para] != 0 ){
                    found_suitable = true;
                }
            }

            p[para] = get_gaussian_sample(p[para], scale * stddev[para], low[para], high[para]);

            return p;
        }

        void loadSigma(string filename, const uint64_t seed){

            ifstream input(filename.c_str(), ios::in);

            const unsigned int dim = 9;

            covar = Eigen::Matrix<double, dim, dim>();
            mean = Eigen::Matrix<double, dim, 1>::Zero();


            string_vector SplitVec;
            string lineStr;
            unsigned long length, lineNum = 0;

            while ( !input.eof() ) {
                getline(input, lineStr);
                string resStr;
                //lineNum++;

                length = lineStr.length();


                if (length > 0) {
                        split( SplitVec, lineStr, is_any_of("\t") );

                        for (unsigned int ii = 0; ii < SplitVec.size(); ii++){
                            covar(ii, lineNum) =   ((2.4*2.4)/double(dim)) * lexical_cast<double>(SplitVec[ii]);
                        }
                        lineNum++;
                }
            }

            norm_cholesk = Eigen::EigenMultivariateNormal<double>(mean,covar,true, seed);

        }
};


template <typename system_>
class init_conds_functor{
public:

	typename system_::state_type operator()(system_ sys){
		typename system_::state_type init_conds = { 0.0 , 0.0 , sys.p[9], 0.0 };
                return init_conds;
        }



};


#endif /* BMEG_MODEL_H_ */
