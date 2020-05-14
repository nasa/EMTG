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

//EMTG Matrix class
//Jacob Englander 1/13/2013

#pragma once

#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <exception>
#include <vector>

#ifdef AD_INSTRUMENTATION
//#include "GSAD_A5.h"
#include "GSAD_2B.h"
#endif

namespace EMTG { namespace math {

    enum MatrixType {standard, identity, Rxhat, Ryhat, Rzhat, Rahat, Rxhatdot, Ryhatdot, Rzhatdot, diagonal, skewsymmetric};

    template <class T> class Matrix
    {
    public:
        //default constructor
        Matrix();

        //constructor for standard matrices when size is known but there is no initialization data
        Matrix(const size_t& n_in, const size_t& m_in);

        //constructor for standard matrices when size is known and an array of data is supplied
        Matrix(const size_t& n_in, const size_t& m_in, const T* input_data);
        Matrix(const size_t& n_in, const size_t& m_in, const std::vector<T>& input_data);

        //constructor for standard matrices when size is known and a scalar is to be assigned to all values
        Matrix(const size_t& n_in, const size_t& m_in, const T& value);

        //constructor for special matrix types, all of which happen to be square
        Matrix(const size_t& n_in, const T& input_data, const MatrixType& type_in);
        Matrix(const size_t& n_in, const std::vector<T>& input_data, const MatrixType& type_in);
        Matrix(const size_t& n_in, const T& input_data, const T& input_data_derivative, const MatrixType& type_in);
        Matrix(const Matrix<T>& input_vector, const MatrixType& type_in);

        //constructor for identity
        Matrix(const size_t& n_in, const MatrixType& type_in);

        //constructor for rotation matrix about an arbitrary axis
        Matrix(const size_t& n_in, const T& input_data, const MatrixType& type_in, const Matrix<T> rotaxis);

        //destructor
        virtual ~Matrix();

        //********************************************************************
        //public methods. These methods are virtual so that other classes can overload them
        
        //value assignment functions
        virtual void assign_all(const T* input_data);
        virtual void assign_all(const std::vector<T>& input_data);
        virtual void assign_entry(const size_t& i, const size_t& j, const T& value);
        virtual void assign_zeros();
        virtual void assign_constant(const T& value);
        virtual void assign_row(const size_t& i, const T* rowvalues);
        virtual void assign_row(const size_t& i, const Matrix& TheRow);
        virtual void assign_column(const size_t& j, const T* columnvalues);
        virtual void assign_column(const size_t& j, const Matrix& TheColumn);
        virtual void assign_diagonal(const std::vector<T>& input_data);
        virtual void assign_diagonal(const T& input_value);
        virtual void clear();
        virtual void resize(const size_t& n_in, const size_t& m_in);
        virtual void resize(const size_t& n_in, const size_t& m_in, const T& value);
        virtual void resize(const size_t& n_in, const size_t& m_in, const MatrixType& type_in);
        virtual void construct_rotation_matrix(const T& angle);
        virtual void construct_rotation_matrix_derivative(const T& angle, const T& angledot);
        virtual void construct_identity_matrix();
        virtual Matrix& operator= (const T& scalar);
        virtual Matrix& operator= (const Matrix& OtherMatrix);
        virtual void shallow_copy(const Matrix& OtherMatrix);
        virtual void shallow_copy(const Matrix& OtherMatrix, const size_t& num_entries);
        

        //get methods
        size_t get_n() const { return this->n; }
        size_t get_m() const { return this->m; }
        //virtual T getentry(const size_t& i, const size_t& j) const;
        virtual Matrix getSubMatrix1D(const size_t& firstEntry, const size_t& lastEntry) const;
        virtual Matrix getrow(const size_t& i) const;
        virtual Matrix getcolumn(const size_t& j) const;
        virtual T operator() (const size_t& i, const size_t& j) const;
        virtual T operator() (const size_t& i) const;
        //virtual T& getreference(const size_t& i, const size_t& j);
        virtual T& operator() (const size_t& i, const size_t& j);
        virtual T& operator() (const size_t& i);
        virtual std::vector<T> getValues() const { return this->values; }


        inline T getentry(const size_t& i, const size_t& j) const
        {
#ifndef FAST_EMTG_MATRIX
            if (i > n - 1 || j > m - 1 || i < 0 || j < 0)
            {
                throw std::invalid_argument("EXCEPTION: Invalid index in [" << i << ", " << j << "]. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }
#endif
            return values[i*m + j];
        }

        inline T& getreference(const size_t& i, const size_t& j)
        {
#ifndef FAST_EMTG_MATRIX
            if (i > n - 1 || j > m - 1 || i < 0 || j < 0)
            {
                throw std::invalid_argument("EXCEPTION: Invalid index in [" << i << ", " << j << "]. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }
#endif
            return values[i*m + j];
        }

        //print methods
        virtual void print_to_screen();
        virtual void print_to_file();
        virtual void print_to_file(const std::string& filename, const bool newFile = true);
        virtual void print_to_file(const std::string& filename, const std::vector< std::string >& rowheader, const bool newFile = true);
        virtual void print_to_file(const std::string& filename, const std::vector< std::string >& rowheader, const std::vector< std::string >& columnheader, const bool newFile = true);

        //basic math operations

        //sign change
        virtual Matrix operator- ();

        //matrix addition
        virtual Matrix operator+ (const Matrix& OtherMatrix) const;
        virtual Matrix& operator+= (const Matrix& OtherMatrix);
        virtual Matrix operator+ (const T& scalar) const;
        virtual Matrix& operator+= (const T& scalar);

        //matrix subtraction
        virtual Matrix operator- (const Matrix& OtherMatrix) const;
        virtual Matrix& operator-= (const Matrix& OtherMatrix);
        virtual Matrix operator- (const T& scalar) const;
        virtual Matrix& operator-= (const T& scalar);

        //right multiplication (left multiplication is done by using the other object's multiplier)
        virtual Matrix operator* (const Matrix& OtherMatrix) const;
        virtual Matrix operator* (const T& scalar) const;
        virtual Matrix& operator*= (const Matrix& OtherMatrix);
        virtual Matrix& operator*= (const T& scalar);
        
        //comparisons
        virtual bool operator== (const Matrix& OtherMatrix) const;
        virtual bool operator!= (const Matrix& OtherMatrix) const;
#ifdef _EMTG_thruth_table
        virtual Matrix<bool> operator<= (const Matrix& OtherMatrix);
        virtual Matrix<bool> operator>= (const Matrix& OtherMatrix);
        virtual Matrix<bool> operator< (const Matrix& OtherMatrix);
        virtual Matrix<bool> operator> (const Matrix& OtherMatrix);
#endif
        
        //multiplication(element-wise)
        virtual Matrix element_multiply(const Matrix& OtherMatrix) const;
        virtual Matrix& element_multiply_equals (const Matrix& OtherMatrix);

        //division (element-wise)
        virtual Matrix operator/ (const T& scalar) const;
        virtual Matrix element_divide(const Matrix& OtherMatrix) const;
        virtual Matrix& operator/= (const T& scalar);
        virtual Matrix& element_divide_equals (const Matrix& OtherMatrix);
        
        //special matrix math
        virtual Matrix transpose() const;
        virtual T determinant() const;
        virtual T determinant(const size_t& n) const;
        virtual Matrix Cofactor(const size_t& p, const size_t& q, const size_t& n) const;
        virtual Matrix adjoint(const size_t& n) const;
        virtual Matrix inverse() const;
        virtual Matrix horz_cat(const Matrix& OtherMatrix) const;
        virtual Matrix vert_cat(const Matrix& OtherMatrix) const;
        virtual Matrix reshape(const size_t& new_n, const size_t& new_m) const;

        //vector math
        virtual T norm() const;
        virtual T norm(const size_t& FirstEntry, const size_t& LastEntry) const;
        virtual T hypnorm(const double& hypconstant) const;
        virtual Matrix unitize() const;
        virtual T dot(const Matrix& OtherMatrix) const;
        virtual Matrix cross(const Matrix& OtherMatrix) const;
        virtual Matrix unitcross(const Matrix& OtherMatrix) const;
        virtual Matrix unitDerivative(const Matrix& OtherMatrix) const;
        virtual Matrix crossDerivative(const Matrix& OtherMatrix) const;
        virtual void cross_in_place(const Matrix& OtherMatrix, Matrix& TargetMatrix);

#ifdef AD_INSTRUMENTATION
        friend Matrix<GSAD::adouble> operator* (const Matrix<double>& otherMatrix, const Matrix<GSAD::adouble> thisMatrix);
#endif        

    protected:
        // fields
        size_t n; //number of rows
        size_t m; //number of columns
        MatrixType type;
        bool issquare;
        bool isvector;
        std::vector<T> values;
    }; //end class definition: Matrix
}}// close namespace