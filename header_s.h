/*
 *  header_s.h
 *  Porj_EffectiveDiffusionNEW
 *
 *  Created by Marino Raffaele on 19/10/15.
 *  Copyright 2015 Raffaele Marino. All rights reserved.
 *
 */

#ifndef HEADER_S_H
#define HEADER_S_H
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iterator>
#include <algorithm>
#include <valarray>
#include <numeric>
#include <complex>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <stdio.h>
#include <unistd.h>
#include <time.h>
#include <random>
#include <armadillo>


#define MAX_SIZE 1000000 



using namespace std;
using namespace arma;
enum{CONF=100000}; //MAX1024*1000000 with -mcmodel=medium
enum{PASSI=1000};
const double PI=3.14159265;
//const size_t  MAX_SIZE_ARRAY=CONF*PASSI;
const double k_B=+1.3806E-2; //femtoNewton*Micrometr*Kelvin{-1}
const double dt=0.001;
const double deltat=0.1;
const double the_time=1000.;
const double myzero=1E-15;

#endif

