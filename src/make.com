#!/usr/bin/sh

g++ -std=c++11 -Wall -O3 do_mcmc_MVN_gaussian.cpp -o do_mcmc_MVN_gaussian -I/opt/local/include/



