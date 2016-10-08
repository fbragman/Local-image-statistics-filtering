#ifndef _Stats_H
#define _Stats_H

//  Statistical analysis functions

#include <iostream>
#include <string>
#include "nifti1_io.h"

using namespace std;

//  given a 1D array of size length, arrayMedian returns the median value of the array
double arrayMedian(double *array,int length);

//  given a 1D array of size length, arrayMax returns the maximum value
double arrayMax(double *array,int length);

//  given a 1D array of size length, arrayMin returns the minimum value
double arrayMin(double *array,int length);

//  given a 1D array of size length, arrayStd returns the standard deviation of the value
double arrayStd(double *array,int length);


#endif

