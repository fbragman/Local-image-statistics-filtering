#include "Image.h"
#include "Stats.h"
#include <iostream>
#include <string>
#include <cmath>

using namespace std;
//Users/felixbragman/Documents/UCL/Programming/C++/CW3/trunk/VolumeMedian.cpp

    // constructor

Image::Image(string filename)
{
    cout << "constructor" << endl;
    nifti_image *nim = nifti_image_read(filename.c_str(), 1);
    
    //  update private member variables
    
    //  grid dimension
    int ndim = nim->dim[0];
    size[0] = nim->nx;  //  x grid or dim[1]
    size[1] = nim->ny;  //  y grid or dim[2]
    if (ndim > 2)       //  volumetric data
    {
        size[2] = nim->nz;  //  z grid or dim[3]
    }
        
    //  voxel spacing (spatial resolution)
    //  xyzt_units = NIFTI_UNITS_MM
    spacing[0] = nim->dx;
    spacing[1] = nim->dy;   // equivalent to pixdim
    spacing[2] = nim->dz;
	
    //  physical location of voxel centre i.e. volume origin
    origin[0] = nim->qoffset_x;
    origin[1] = nim->qoffset_y;
    origin[2] = nim->qoffset_z;
    
    //  number of components per voxel
    //  noOfComponents = 1 -> scalar valued
    //  noOfComponents = 2 -> vector
    if (nim->nu > nim->nt)
    {
        noOfComponents = nim->nu;
    }
    else
    {
        noOfComponents = nim->nt;
    }
    
    //  save data as a 1D array
    int numel = size[0]*size[1]*size[2]*noOfComponents;
    
    if(nim->datatype == 16)
    {
        float *dummy = static_cast<float*>(nim->data);
        data = new double[numel];
        for (int i = 0; i < numel; ++i)
        {
            data[i] = dummy[i];
        }
        delete [] dummy;
    }
    
    //  free image
    
    //nifti_image_free(nim);
    
}

    //  copy constructor
    //  NOT WORKING - Segmentation Fault: 11 at line 89

Image::Image(const Image& rhs)
{
    cout << "copy constructor" << endl;
    
    for (int n = 0; n < (rhs.noOfComponents); ++n)
    {
        for (int i = 0; i < rhs.size[0]; ++i)
        {
            for (int j = 0; j < rhs.size[1]; ++j)
            {
                for (int k = 0; k < rhs.size[2]; ++k)
                {
                    this->setVoxC(i,j,k,n,rhs.data[rhs.sub2ind(i,j,k,n)]);
                }
            }
        }
    }
}

    //  destructor

Image::~Image()
{
	cout << "destructor" << endl;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // auxillary private member function

    //  index converter :: [i,j,k,nCo] -> idx
    int Image::sub2ind(const int i,const int j,const int k,const int nCo) const
    {
        int vol = size[0]*size[1]*size[2];
        int idx = ((i + j*size[0] + k*size[0]*size[1]) + (vol*nCo));
    return idx;
    }

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    //  accessor functions

    //  access grid dimensions / image size
    int *Image::getSize()
    {
        int *result = new int[3];
        
        for (int i = 0; i<3; ++i)
        {
            result[i] = size[i];
        }
    return result;
    }
    
    //  access voxel spacing / determine spatial resolution
    double *Image::getRes()
    {
        double *result = new double[3];
        
        for (int i = 0; i<3; ++i)
        {
            result[i] = spacing[i];
        }
    return result;
    }

    //  access voxel origin

    double *Image::getOrigin()
    {
        double *result = new double[3];
    
        for (int i = 0; i<3; ++i)
            {
                result[i] = origin[i];
            }
    return result;
    }

    //  access number of components at each voxel

    int Image::getNoOfComponents()
    {
    return noOfComponents;
    }

    int Image::getIdx(const int i,const int j, const int k)
    {
        int idx = (i + j*size[0] + k*size[0]*size[1]);
    return idx;
    }

    double *Image::physLoc(const int i, const int j, const int k)
    {
        double *centre = new double[3];
        centre[0] = origin[0] + i*spacing[0];
        centre[1] = origin[1] + j*spacing[1];
        centre[2] = origin[2] + k*spacing[2];
    return centre;
    }

    double *Image::getContLoc(const double x, const double y, const double z)
    {
        double location[3] = {x,y,z};
        double *contLoc = new double[3];
        for (int i = 0; i<3; ++i)
        {
            contLoc[i] = (location[i] - origin[i])/spacing[i];
        }
    return contLoc;
    }

    //  normal image index
    double Image::getVoxel(const int i, const int j, const int k, const int nCo) const
    {
            if ( (i < 0 | i > size[0]) | (j < 0 | j > size[1]) | (k < 0 | k > size[2]))
            {
                return 0.0;
            }
            
            else
            {
                int voxelIdx = (this->sub2ind(i,j,k,nCo));
                return this->data[voxelIdx];
            }
    }

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    //  interpolate

    double Image::triInterpolate(const int x, const int y, const int z,const int nCo)
    {
        //  step 1: coordinates lattice points of grid surrounding C
        //          e.g.    c0000 = (xo,yo,zo)
        //                  c1111 = (x1,y1,z1)
        double xo = floor(x);
        double x1 = ceil(x);
        double yo = floor(y);
        double y1 = ceil(y);
        double zo = floor(z);
        double z1 = ceil(z);
        
        //  step 2: calculate xd,yd,zd
        double xd = x - xo;
        double yd = y - yo;
        double zd = z - zo;
        
        //  step 3: nearest neighbour continuous voxel indices
        double *contIdx = new double[3];
        contIdx = this->getContLoc(x, y, z);
        
        int idxo = floor(contIdx[0]);
        int idx1 = ceil(contIdx[0]);
        int idyo = floor(contIdx[1]);
        int idy1 = ceil(contIdx[1]);
        int idzo = floor(contIdx[2]);
        int idz1 = ceil(contIdx[2]);
        delete [] contIdx;
        
        double voxInterpolate = 0.0;
            
        //for (int i = 0; i < noOfComponents; ++i)
        //{
            //  step 4: linear interpolation to find coo,co1,c1o,c11
            double c00 = this->getVoxel(idxo, idyo, idzo, nCo)*(1-xd) + this->getVoxel(idx1, idyo, idzo, nCo)*xd;
            double c10 = this->getVoxel(idxo, idy1, idzo, nCo)*(1-xd) + this->getVoxel(idx1, idy1, idzo, nCo)*xd;
            double c01 = this->getVoxel(idxo, idyo, idz1, nCo)*(1-xd) + this->getVoxel(idx1, idyo, idz1, nCo)*xd;
            double c11 = this->getVoxel(idxo, idy1, idz1, nCo)*(1-xd) + this->getVoxel(idx1, idy1, idz1, nCo)*xd;
            
            //  step 5: linear interpolation to find co,c1
            double c0 = c00*(1-yd) + c10*yd;
            double c1 = c01*(1-yd) + c11*yd;
            
            //  step 6: linear interpolation to find c
            voxInterpolate = c0*(1-zd) + c1*zd;
        //}
    return voxInterpolate;
    }

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    //  modifiers

    void Image::setVoxC(const int i, const int j, const int k, const int nCo,double setMag)
    {
        this->data[this->sub2ind(i,j,k,nCo)] = setMag;
    }

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    //  operator << overloading

    ostream& operator<< (ostream& out,const Image& nim)
    {
        out << "Number of components per voxel: " << nim.noOfComponents << endl;
        out << "Image size (voxels): [ " << nim.size[0] << " " << nim.size[1] << " " << nim.size[2] << " ]" << endl;
        out << "Voxel dimensions (mm): [ " << nim.spacing[0] << " " << nim.spacing[1] << " " << nim.spacing[2] << " ]" << endl;
        out << "Origin voxel centre (mm): [ " << nim.origin[0] << " " << nim.origin[1] << " " << nim.origin[2] << " ]" << endl;
    return out;
    }

    //  operator += overloading
    Image& Image::operator+= (const Image& rhs)
    {
        cout << "+= overloading " << endl;
        for (int nCo = 0; nCo < (rhs.noOfComponents-1); ++nCo)
        {
            for (int i = 0; i <= rhs.size[0]; ++i)
            {
                for (int j = 0; j <= rhs.size[1]; ++j)
                {
                    for (int k = 0; k <= rhs.size[2]; ++k)
                    {
                        this->data[this->sub2ind(i,j,k,nCo)] += rhs.data[rhs.sub2ind(i,j,k,nCo)];
                    }
                }
            }
        }
    return *this;
    }

    //  operator + overloading
    Image operator+ (const Image& rhs, const Image& lhs)
    {
        //  operate over private member variable data in each object rhs/lhs
        //  dimension and number of scalar components of lhs equals that of rhs
        //  data : 1D array for all components
        cout << " + overloading " << endl;
        Image tmp(lhs);
        tmp += rhs;
        return tmp;
    }
     
    //  operator = overloading
    Image& Image::operator= (const Image& rhs)
    {
        cout << " = overloading " << endl;

        this->data = rhs.data;
        return *this;
    }

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    //  saving image

    //  method: read in input
    //          change data values

    void Image::saveVol(const string & filename, const string & outputFilename)
    {
        //  reading in original image
        nifti_image *nim = nifti_image_read(filename.c_str(), 1);
        
        //  total number of voxels in the n-D data
        int numel = size[0]*size[1]*size[2]*noOfComponents;
        
        double *dummy = static_cast<double*>(this->data);
        float* ndata = new float[numel];
        for (int i = 0; i < numel; ++i)
        {
            ndata[i] = dummy[i];
        }
        nim->data = ndata;
        delete [] dummy;
        
        nim->fname = new char[outputFilename.size() + 1];
        strcpy(nim->fname, outputFilename.c_str());
        
        //  image write
        nifti_image_write(nim);
        
        //  image unloading
        nifti_image_free(nim);
    }
