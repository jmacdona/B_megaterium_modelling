
Please go to our new repository for a more general and better version of this software - [ODE\_MCMC\_tools](https://github.com/jmacdona/ODE_MCMC_tools)

# Bayesian statistical inference of Bacillus megaterium cell-free transcription and translation using Markov chain Monte Carlo (MCMC)

## Requirements

The software has currently been tested with Linux (RHEL 7) and Mac OS X with MacPorts. The software requires the Eigen3 and Boost C++ libraries (including the development header files).

## Compiling the C++ code

To compile the software go into the src/ directory and run the command:
```
./make.com
```

## Running MCMC on the example data

To run MCMC on the example data, go into the example\_data/ directory and run the command
```
../src/do_mcmc_MVN_gaussian 3453245  Sigma.dat pKMBm5-MGapt_GFP_simple_diag_list_of_inputs.dat pKMBm5-MGapt_mRNA_simple_diag_list_of_inputs.dat start_params1.dat  0.1 20
```

Where 3453245 is the random seed, Sigma.dat is the multivariate normal (MVN) proposal distribution covariance matrix, pKMBm5-MGapt\_GFP\_simple\_diag\_list\_of\_inputs.dat and pKMBm5-MGapt\_mRNA\_simple\_diag\_list\_of\_inputs.dat are files listing the experimental conditions and corresponding experimental data files (in the form of a vector of means and a covariance matrix), start\_params1.dat is a vector describing the initial start parameters, and 0.1 and 20 are minimum standard deviation values permitted for each experimental data point.

## License

This software is distributed under the GNU GPL license, version 3.

(C) James T. MacDonald, 2016. 
Imperial College London.





