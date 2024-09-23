// EMTG: Evolutionary Mission Trajectory Generator
// An open-source global optimization tool for preliminary mission design
// Provided by NASA Goddard Space Flight Center
//
// Copyright (c) 2013 - 2020 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// All Other Rights Reserved.

// Licensed under the NASA Open Source License (the "License"); 
// You may not use this file except in compliance with the License. 
// You may obtain a copy of the License at:
// https://opensource.org/licenses/NASA-1.3
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either 
// express or implied.   See the License for the specific language
// governing permissions and limitations under the License.

//EMTG Tensor class
//Jacob Englander 11/23/2018

#pragma once

#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <exception>
#include <vector>

#include "EMTG_Matrix.h"

#ifdef AD_INSTRUMENTATION
//#include "GSAD_A5.h"
#include "GSAD_2B.h"
#endif

namespace EMTG { namespace math {


    template <class T> class Tensor
    {
    public:
        //default constructor
        Tensor();

        //constructor for tensor when size is known but there is no input data
        Tensor(const size_t& n, const size_t& m, const size_t& p);

        //constructor for tensor when size is known but there is input data
        Tensor(const size_t& n, const size_t& m, const size_t& p, std::vector< Matrix<T> > input_data);

        //constructor for cubic tensors
        Tensor(const size_t& n, std::vector< Matrix<T> > input_data);

        //print methods
        virtual void print_to_screen();
        virtual void print_to_file();
        virtual void print_to_file(const std::string& filename);
        virtual void print_to_file(const std::string& filename, const std::vector< std::string >& rowheader);
        virtual void print_to_file(const std::string& filename, const std::vector< std::string >& rowheader, const std::vector< std::string >& columnheader);

        //assignment
        virtual void assign_all(std::vector< Matrix<T> >& input_data);
        virtual void assign_entry(const size_t& i, const size_t& j, const size_t& k, const T& value);
        virtual void assign_zeros();

        //get
        size_t get_n() const { return this->n; }
        size_t get_m() const { return this->m; }
        size_t get_p() const { return this->p; }
        Matrix<T>& getSlice(const size_t& k) { return this->values[k]; }
        virtual T operator() (const size_t& i, const size_t& j, const size_t& k) const { return this->getentry(i, j, k); };
        virtual T& operator() (const size_t& i, const size_t& j, const size_t& k) { return this->getreference(i, j, k); };

        inline T getentry(const size_t& i, const size_t& j, const size_t& k) const
        {
#ifndef FAST_EMTG_MATRIX
            if (i >= this->n || j >= this->m || k >= this->p || i < 0 || j < 0 || p < 0)
            {
                throw std::invalid_argument("EXCEPTION: Invalid index in[" + std::to_string(i) + ", " + std::to_string(j) + ", " + std::to_string(k) + "]. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }
#endif
            return this->values[k](i*this->m + j);
        }//end getentry()

        inline T& getreference(const size_t& i, const size_t& j, const size_t& k)
        {
#ifndef FAST_EMTG_MATRIX
            if (i >= this->n || j >= this->m || k >= this->p || i < 0 || j < 0 || p < 0)
            {
                throw std::invalid_argument("EXCEPTION: Invalid index in[" + std::to_string(i) + ", " + std::to_string(j) + ", " + std::to_string(k) + "]. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }
#endif
            return this->values[k](i*this->m + j);
        }//end getreference()

        //multiplication
        Matrix<T> bullet2_vector(const Matrix<T>& vector) const;

        //matrix addition
        virtual Tensor operator+ (const Tensor& OtherTensor) const;
        virtual Tensor& operator+= (const Tensor& OtherTensor);
        virtual Tensor operator+ (const T& scalar) const;
        virtual Tensor& operator+= (const T& scalar);

        //Tensor subtraction
        virtual Tensor operator- (const Tensor& OtherTensor) const;
        virtual Tensor& operator-= (const Tensor& OtherTensor);
        virtual Tensor operator- (const T& scalar) const;
        virtual Tensor& operator-= (const T& scalar);

        //multiplication(element-wise)
        virtual Tensor element_multiply(const Tensor& OtherTensor) const;
        virtual Tensor& element_multiply_equals(const Tensor& OtherTensor);

        //division (element-wise)
        virtual Tensor operator/ (const T& scalar) const;
        virtual Tensor element_divide(const Tensor& OtherTensor) const;
        virtual Tensor& operator/= (const T& scalar);
        virtual Tensor& element_divide_equals(const Tensor& OtherTensor);

    protected:
        // fields
        size_t n; //number of rows
        size_t m; //number of columns
        size_t p; //number of slices
        bool iscube;
        std::vector< Matrix<T> > values;
    }; //end class definition: Tensor
}}// close namespace