//
//  VolumeQuery.cpp
//  CW3
//
//  Created by Felix Bragman on 18/12/2012.
//  Copyright (c) 2012 Felix Bragman. All rights reserved.
//

#include <iostream>
#include <string>
#include <sstream>
#include "Image.h"
#include "Stats.h"
#include <cmath>
#include <stdlib.h>
#include <iomanip>

using namespace std;

int main(int argc, const char * argv[])
{
    if (argc < 6)
    {
		cout << "Error - " << "command linke syntax - ";
        cout << "VolumeQuery nifti_filename x y z size_cube" << endl;
		return 1;
    }
    
    //  command line inputs
    string filename = argv[1];
    
    string outputFilename = argv[2];
    
    string argument = argv[3];
    double x;
    istringstream(argument) >> x;
    
    argument = argv[4];
    double y;
    istringstream(argument) >> y;
    
    argument = argv[5];
    double z;
    istringstream(argument) >> z;
    
    argument = argv[6];
    int nhood;
    istringstream(argument) >> nhood;
    
    if (nhood%2 == 0 | nhood < 0)
    {
        cout << "size of neighbourhood cube must be positive odd voxel in width" << endl;
        return 1;
    }
        
    //  load image
    Image *nim = new Image(filename);
    
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    //  image properties
    
    //  obtain image properties
    int *size = nim->getSize();             // grid size
    double *spacing = nim->getRes();        // voxel size (mm)
    int nCo = nim->getNoOfComponents();     // number of components per voxel
    double *ori = nim->getOrigin();         // origin voxel centre location (mm)
    
    //  kernel properties
    int kernelN = nhood*nhood*nhood;        // number of voxels in user defined grid
       
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    //  local neighbourhood
    
    //  set up on component 1 then use index of values to perform statistical analysis on all components
    //  note:   input size is in voxels
    
    //  get continuous voxel index at centre of user-defined cube/grid
    double *centreIdx = nim->getContLoc(x,y,z);
    
    //  voxel centre of boundary voxels in user-defined grid
    //  i.e.    if boundPhys[0] = 10 and spacing[0] = 2, physical boundary xo = 9
    //  boundPhys[0] = xo, boundPhys[1] = x1, boundPhys[2] = yo, boundPhys[3] = y1, boundPhys[4] = zo, boundPhys[5] = z1,
    double *boundPhys = new double[6];
    boundPhys[0] = x - ((nhood-1)/2)*spacing[0];
    boundPhys[1] = x + ((nhood-1)/2)*spacing[0];
    boundPhys[2] = y - ((nhood-1)/2)*spacing[1];
    boundPhys[3] = y + ((nhood-1)/2)*spacing[1];
    boundPhys[4] = z - ((nhood-1)/2)*spacing[2];
    boundPhys[5] = z + ((nhood-1)/2)*spacing[2];
    
    //  continuous index of user-defined grid voxel centre boundaries
    double *boundCont = new double[6];
    boundCont[0] = (boundPhys[0] - ori[0])/spacing[0];
    boundCont[1] = (boundPhys[1] - ori[0])/spacing[0];
    boundCont[2] = (boundPhys[2] - ori[1])/spacing[1];
    boundCont[3] = (boundPhys[3] - ori[1])/spacing[1];
    boundCont[4] = (boundPhys[4] - ori[2])/spacing[2];
    boundCont[5] = (boundPhys[5] - ori[2])/spacing[2];

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    //  looping over all components
    
    for (int component = 0; component < nCo; ++component)
    {
        double maxVal = 0.0;
        double minVal = 0.0;
        double medVal = 0.0;
        double stdVal = 0.0;
        double voxInt = nim->triInterpolate(x,y,z,component); // triInterpolate deals with voxels outside of image
        
        //  convert physical boundaries to continuous indices
        
        //  algorithm
        //  starting from physical location of the cube (C) boundaries: iterate over physical location of voxel centres of grid C
        //  i,j,k are physical locations
        //  convert i,j,k to continuous indices and then image indices: check if image index completely enclosed by physical boundaries
        
        int total = 0;
                        
        double *dummy = new double[3];
        double *cInt = new double[kernelN]; // array for all voxels in local neighbourhood

        int voxelI = 0;
        for (double i = boundPhys[0]; i <= boundPhys[1]; i += spacing[0])
        {
        int voxelJ = 0;
            for (double j = boundPhys[2]; j <= boundPhys[3]; j += spacing[1])
            {
            int voxelK = 0;
                for (double k = boundPhys[4]; k <= boundPhys[5]; k += spacing[2])
                {
                    dummy = nim->getContLoc(i,j,k);
                                        
                    //  calculate boundary of voxel
                    //  if voxel boundary is within C: get intensity
                    //  if not: trilinear interpolation from the physical location of the voxel not enclosed
                                    
                    //  find nearest voxel in image grid to voxel centre boundary in user defined grid
                    int xo = (int)(dummy[0]+0.5);
                    int yo = (int)(dummy[1]+0.5);
                    int zo = (int)(dummy[2]+0.5);
                                
                    if ( (xo<boundCont[0] | xo>boundCont[1]) | (yo<boundCont[2] | yo>boundCont[3]) | (zo<boundCont[4] | zo>boundCont[5]) )
                    {
                        cInt[voxelK + voxelJ*nhood + voxelI*nhood*nhood] = nim->triInterpolate(i,j,k,component);
                    }
                    else
                    {
                        dummy[0] = (int)(dummy[0]+0.5);
                        dummy[1] = (int)(dummy[1]+0.5);
                        dummy[2] = (int)(dummy[2]+0.5);
                        cInt[voxelK + voxelJ*nhood + voxelI*nhood*nhood] = nim->getVoxel(dummy[0],dummy[1],dummy[2],component);
                    }
                voxelK++;
                }
            voxelJ++;
            }
        voxelI++;
        }
        
        //  regional statistical analysis output
        maxVal = arrayMax(cInt,kernelN);
        minVal = arrayMin(cInt,kernelN);
        medVal = arrayMedian(cInt,kernelN);
        stdVal = arrayStd(cInt,kernelN);
        
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        //  regional analysis printing
        
        //  flush left
        cout << left;
        
        //  print out table header
        if (component==0)
        {
            cout << endl;
            
            cout << "Regional analysis using a kernel size of " << nhood << " voxels^3 at a voxel located at ( " << x << ", " << y << ", " << z << " )" << endl << endl;
            
            cout << setw(16) << "Component #" << setw(13) << "Voxel value" << setw(12) << "Maximum" << setw(12) << "Minimum" << setw(21) << "Standard Deviation" << setw(8) << "Median" << endl << endl;
        }
        
        //  table contents
        cout << setw(16) << (component+1) << setw(13) << voxInt << setw(12) << maxVal << setw(12) << minVal << setw(21) << stdVal << setw(8) << medVal << endl;
        
    }
return 0;
}

