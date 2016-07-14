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

#ifndef MULTIVARIATE_NORMAL_H_INCLUDED
#define MULTIVARIATE_NORMAL_H_INCLUDED

#ifndef M_PI
#define M_PI REAL(3.1415926535897932384626433832795029)
#endif

#include <eigen3/Eigen/Dense>
#include <iostream>
#include <fstream>
#include<boost/algorithm/string.hpp>
#include<boost/lexical_cast.hpp>
#include<boost/algorithm/string/split.hpp>
#include <vector>
#include  <cmath>

#define MAXBUFSIZE  ((int) 1e6)


typedef std::vector<std::string> string_vector;



Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> readDynamicMatrix(std::string filename)
{
    int cols = 0, rows = 0;
    double buff[MAXBUFSIZE];

    // Read numbers from file into buffer.
    std::ifstream infile;
    infile.open(filename);
    while (! infile.eof())
    {
        std::string line;
        getline(infile, line);

        int temp_cols = 0;
        std::stringstream stream(line);
        while(! stream.eof())
            stream >> buff[cols*rows+temp_cols++];

        if (temp_cols == 0)
            continue;

        if (cols == 0)
            cols = temp_cols;

        rows++;
    }

    infile.close();

    rows--;

    // Populate matrix with numbers.
    Eigen::MatrixXd result(rows,cols);
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            result(i,j) = buff[ cols*i+j ];

    return result;
}


        template <typename SCALAR,unsigned int DIM_ROW, unsigned int DIM_COL>
        Eigen::Matrix<double, DIM_ROW, DIM_COL> readMatrix(std::string filename){

            std::ifstream input(filename.c_str(), std::ios::in);

            //const unsigned int dim = 11;

            Eigen::Matrix<SCALAR, DIM_ROW, DIM_COL> mat = Eigen::Matrix<SCALAR, DIM_ROW, DIM_COL>::Zero();
            //mean = Eigen::Matrix<double, dim, 1>::Zero();


            string_vector SplitVec;
            std::string lineStr;
            unsigned long length, lineNum = 0;

            while ( !input.eof() ) {
                getline(input, lineStr);
                std::string resStr;
                //lineNum++;

                length = lineStr.length();


                if (length > 0) {
                        boost::split( SplitVec, lineStr, boost::is_any_of("\t") );

                        for (unsigned int ii = 0; ii < SplitVec.size(); ii++){
                            //cout << lineNum << "\t" << ii << "\t" << SplitVec[ii] << endl;
                            mat(lineNum, ii) =   boost::lexical_cast<SCALAR>(SplitVec[ii]);
                        }
                        lineNum++;
                }
            }

            return mat;
        }

class MultivariateNormalPDF{
private:

    Eigen::Matrix<double,Eigen::Dynamic,1> mean ;
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> covar;
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> covar_inv;

    double det, logdet;
    double normalising_term, log_normalising_term;

public:




    MultivariateNormalPDF(){
    }

    Eigen::Matrix<double,Eigen::Dynamic,1> get_mean(){
        return mean;
    }

    void setCovar(Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> covar_, const bool use_factorization = true){
        covar = covar_;

        if (use_factorization == false){
            covar_inv = covar.inverse();
            det = covar.determinant();
            normalising_term = 1.0 / std::sqrt(std::pow(M_PI, covar.rows()) * det);
            log_normalising_term = std::log(normalising_term);
        }
        else {

            Eigen::FullPivHouseholderQR< Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > qr(covar);
            covar_inv = qr.inverse();
            logdet = qr.logAbsDeterminant();
            det = std::exp(logdet);
            log_normalising_term = 0.5 * (-std::log(std::pow(M_PI, covar.rows())) - logdet);
            normalising_term = std::exp(log_normalising_term);
        }
    }

    void setMean(Eigen::Matrix<double,Eigen::Dynamic,1> mean_){
        mean = mean_;
    }

    double  getPdf(Eigen::Matrix<double,Eigen::Dynamic,1> val){

        return std::exp(this->getLogPdf(val)); //prob;
    }


    double  getLogPdf(Eigen::Matrix<double,Eigen::Dynamic,1> val){
        Eigen::Matrix<double,Eigen::Dynamic,1>  val_cent = val - mean;

        const double exp_term =  (-0.5 * (val_cent.transpose()* covar_inv * val_cent))(0,0);

        const double log_prob = log_normalising_term + exp_term;

        return log_prob;
    }

    double  getLogPdfNoNorm(Eigen::Matrix<double,Eigen::Dynamic,1> val){
        Eigen::Matrix<double,Eigen::Dynamic,1>  val_cent = val - mean;
        const double exp_term =  (-0.5 * (val_cent.transpose()* covar_inv * val_cent))(0,0);
        const double log_prob = exp_term;

        return log_prob;
    }

};


#endif // MULTIVARIATE_NORMAL_H_INCLUDED
