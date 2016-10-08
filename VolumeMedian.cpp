//
//  VolumeMedian.cpp
//  CW3
//
//  Created by Felix Bragman on 20/12/2012.
//  Copyright (c) 2012 Felix Bragman. All rights reserved.
//

#include <iostream>
#include <string>
#include <sstream>
#include "Image.h"
#include "Stats.h"
#include <cmath>

using namespace std;

int main(int argc, const char * argv[])
{
    if (argc < 4)
    {
		cout << "Error - " << "command linke syntax - ";
        cout << "VolumeQuery nifti_filename output_nifti_filename kernel_size" << endl;
		return 1;
    }
    
    //  command line inputs
    string filename = argv[1];
    
    string outputFilename = argv[2];
    
    string argument = argv[3];
    int kernel;
    istringstream(argument) >> kernel;
    
    if (kernel%2 == 0 | kernel<0)
    {
        cout << "kernel dimension must be a positive integer" << endl;
        return 1;
    }
    
    //  load image
    Image *nim = new Image(filename);
    //  load dummy image
    Image *dummy = new Image(filename);
    cout << "loading done" << endl;
    
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    //  image properties
    
    int *size = nim->getSize();             // grid size
    int nCo = nim->getNoOfComponents();     // number of components per voxel
    
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    //  averaging over i,j,k coordinates for all components
    
    //  for each voxel: loop over all neighbouring voxels needed for averaging
    //                  let (i,j,k) be centre voxel
    //                  loop from (i-dx,j-dy,k-dz) to (i+dx,j+dy,k+dz)
    
    //  method: read in image second time, use one image for analysis, one image for dumping new data
    //  note:   memory intensive, use copy constructor for private member variable data in the future
    
    int kernelN = kernel*kernel*kernel;
    
    cout << "Filtering in progress ..." << endl;
    
    for (int n = 0; n < nCo; ++n)
    {
        cout << 100*n/nCo << "%" << endl;
        for (int i = 0; i < size[0]; ++i)
        {
            for (int j = 0; j < size[1]; ++j)
            {
                for (int k =0; k < size[2]; ++k)
                {
                    double *tmp = new double[kernelN];
                    double tmpMed = 0.0;
                    
                    //  iterators (voxel position in neighbourhood)
                    int xIt = 0;
                    for (int kernelX = (i-(kernel-1)/2); (kernelX <= (i+(kernel-1)/2)); kernelX++)
                    {
                        int yIt = 0;
                        for (int kernelY = (j-(kernel-1)/2); (kernelY <= (j+(kernel-1)/2)); kernelY++)
                        {
                            int zIt = 0;
                            for (int kernelZ = (k-(kernel-1)/2); (kernelZ <= (k+(kernel-1)/2)); kernelZ++)
                            {
                                tmp[(zIt) + (yIt)*kernel + (xIt)*kernel*kernel] = nim->getVoxel(kernelX,kernelY,kernelZ,n);
                                zIt++;
                            }
                            yIt++;
                        }
                        xIt++;
                    }
                    //  For each set of voxels within the kernel, array -> arrayMedian : returns median (double)
                    tmpMed = arrayMedian(tmp,kernelN);
                    //  Input median into dummy image;
                    dummy->setVoxC(i,j,k,n,tmpMed);
                    delete [] tmp;
                }
            }
        }
    }
    
    cout << "Filtering complete" << endl;
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    //  output of averaged volume (in dummy)
    delete nim;
    cout << "Saving new volume" << endl;
    dummy->saveVol(filename,outputFilename);
    
    return 0;
}