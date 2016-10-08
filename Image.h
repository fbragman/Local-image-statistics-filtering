#ifndef _Image_H
#define _Image_H

//  Class Description
//  Multi-dimensional image data based on nifticlib libraries

#include <iostream>
#include <string>
#include "nifti1_io.h"

using namespace std;

class Image {
	private:
    
    int size[3];            //  Grid dimension
    double spacing[3];      //  Voxel size [mm]
    double origin[3];       //  Physical location of voxel centre at [0,0,0] [mm]
    int noOfComponents;     //  Number of components per voxel: 1 = Scalar n = Vector
    double *data;           //  Image data stored as a 1D array


    int sub2ind(const int i,const int j,const int k,const int nCo) const; //  using [i,j,k] coordinates, map to 1D index of 3D volumetric data
	
	public:

	//  constructor
    Image (string filename);
    
    //  copy constructor
    Image (const Image& rhs);
    
	//  destructor
    ~Image();

	//  accessors
    int *getSize();
    double *getOrigin();
    double *getRes();
    int getNoOfComponents();
    double getVoxel(const int i, const int j, const int k, const int nCo) const; // using grid index
    double *physLoc(const int i, const int j, const int k);
    double *getContLoc(const double x, const double y, const double z);
    
    //  modifier
    void setVoxC(const int i, const int j, const int k, const int nCo,double setMag);
    
    //  trilinear interpolation
    //  step 1: find lattice point coordinates surrounding point C to be interpolated
    //  step 2: linear interpolation to find coo,co1,c1o,c11
    //  step 3: linear interpolation to find co,c1
    //  step 4: linear interpolation to find c
    
    double triInterpolate(const int x, const int y, const int z,const int nCo); // i = iterator for voxel components
    int    phys2idx(const int x, const int y, const int z);
    
    //  test function
    int getIdx(const int i,const int j, const int k);
        
    //  print
    void print() const;
    
    //  operator << overloading
    friend ostream& operator<< (ostream& out,const Image& nim);
    
    //  operator + overloading
    Image& operator+= (const Image& rhs);
    friend Image operator+ (const Image& rhs,const Image& lhs);
    
    //  operator = overloading
    Image& operator= (const Image& rhs);
    
    //  saving image
    void saveVol(const string& filename, const string& outputFilename);
};

#endif

