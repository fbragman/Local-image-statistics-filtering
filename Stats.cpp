#include "Image.h"
#include <iostream>
#include <string>
#include <cmath>

using namespace std;

double arrayMedian(double *array,int length)
{
    double arrayMedian = 0.0;
    int idx = 0;
    
    //  sort the array
    sort(array,array+length);

    //  find median value
    if (length % 2 == 0)    // even
    {
        idx = length/2;
        arrayMedian = (array[idx-1] + array[idx])/2.0;
    }
    else                    // odd
    {
        idx = (length-1)/2;        
        arrayMedian = (array[idx]);
    }
return arrayMedian;
}

double arrayMax(double *array,int length)
{
    double arrayMax = 0.0;
    
    //  sort the array
    sort(array,array+length);
    
    arrayMax = array[length-1];
    return arrayMax;
}

double arrayMin(double *array, int length)
{
    double arrayMin = 0.0;
    
    //  sort the array
    sort(array,array+length);
    
    arrayMin = array[0];
    
    return arrayMin;
}

double arrayStd(double *array,int length)
{
    double arrayStd = 0.0;
    double average = 0.0;
    double tmp = 0.0;
    length = (double)(length);
    
    //  average calculation
    for (int i = 0; i<length; ++i)
    {
        average += array[i];
    }
    
    average = average / length;
    //  variance
    for (int i = 0; i<length; ++i)
    {
        tmp += pow((array[i] - average),2);
    }
    
    //  standard deviation
    arrayStd = pow((tmp/length),0.5);

    return arrayStd;
}

 
