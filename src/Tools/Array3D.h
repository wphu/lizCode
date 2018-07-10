#ifndef ARRAY3D_H
#define ARRAY3D_H

#include <vector>
#include <iostream>

// Class Array3D: template class as 3D array
template<typename T>
class Array3D
{

public:

    // Constructor for Array3D: with no input argument
    Array3D() 
    {
        data = NULL;
    };

    //! Destructor for Array3D
    ~Array3D() 
    {
    } ;

    void allocate_dims(std::vector<unsigned int> dims_temp)
    {
        dims = dims_temp;
        if(dims.size() != 3) 
        {
            ERROR("Alloc error must be 4 : " << dims.size());
        }

        if(data)
        {
            delete [] data;
        }

        data = new T[dims[0]*dims[1]*dims[2]];
        data_3D= new T**[dims[0]];
        for (int i = 0; i < dims[0]; i++)
        {
            data_3D[i] = new T*[dims[1]];
            for (int j = 0; j < dims[1]; j++)
            {
                data_3D[i][j] = data + i * dims[1] * dims[2] + j * dims[2];
            }
        }

        //DEBUG(10,"Fields 3D created: " << dims_[0] << "x" << dims_[1] << "x" << dims_[2]);
        dim_global = dims[0] * dims[1] * dims[2];
    };


    // returns the dimension of the Array3D
	inline std::vector<unsigned int> get_dims () const
    {
        return dims;
    }

    // method used to put all entry of a Array3D at a given value val
    inline void put_to(T val)
    {
        if(data)
        {
            for(int i = 0; i < dim_global; i++)
            {
                data[i] = val;
            }
        }
            
    }


    //! 3D reference access to the linearized array (with check in DEBUG mode)
    inline T& operator () (unsigned int i, unsigned int j, unsigned int k)
    {
        return data_3D[i][j][k];
    };
    //! 3D access to the linearized array (with check in DEBUG mode)
    inline double operator () (unsigned int i, unsigned int j,unsigned int k) const
    {
        return data_3D[i][j][k];
    };


protected:

private:
    // pointer to the linearized array
    T* data;

    // data_3D[nx][ny][nz][nv]
    T ***data_3D;

    // vector containing the dimensions of the Array3D
    std::vector<unsigned int> dims;

    // All arrays may be viewed as a 1D array
    // Linearized diags
    int dim_global;

};

#endif
