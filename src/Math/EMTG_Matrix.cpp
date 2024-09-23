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

#include "EMTG_Matrix.h"
#include <exception>

namespace EMTG 
{ 
    namespace math
    {
        //***************************************************
        //method definitions
        template <class T> Matrix<T>::Matrix() :
            n(1),
            m(1),
            values(1, 0.0),
            type(standard),
            issquare(true),
            isvector(true)
        {
            //default constructor constructs a 1x1 matrix initialized to zero
        }

        //constructor when size is known but there is no initialization data
        template <class T> Matrix<T>::Matrix(const size_t& n_in, const size_t& m_in) : 
            n(n_in),
            m(m_in),
            values(n * m, 0.0),
            type(standard),
            isvector(n == 1 || m == 1 ? true : false),
            issquare(n == m ? true : false)
        {
        }

        //constructor when size is known and an array of data is supplied
        template <class T> Matrix<T>::Matrix(const size_t& n_in, const size_t& m_in, const T* input_data) : 
            n(n_in), 
            m(m_in),
            values(n * m, 0.0),
            type(standard),
            isvector(n == 1 || m == 1 ? true : false),
            issquare(n == m ? true : false)
        {
            assign_all(input_data);
        }

        template <class T> Matrix<T>::Matrix(const size_t& n_in, const size_t& m_in, const std::vector<T>& input_data) :
            n(n_in),
            m(m_in),
            values(n * m, 0.0),
            type(standard),
            isvector(n == 1 || m == 1 ? true : false),
            issquare(n == m ? true : false)
        {
            assign_all(input_data);
        }

        template <class T> Matrix<T>::Matrix(const size_t& n_in, const size_t& m_in, const T& value) :
            n(n_in),
            m(m_in),
            values(n * m, value),
            type(standard),
            isvector(n == 1 || m == 1 ? true : false),
            issquare(n == m ? true : false)
        {
        }

    
        template <class T> Matrix<T>::Matrix(const size_t& n_in, const T& input_data, const MatrixType& type_in) :
            n(n_in), 
            m(n_in),
            values(n * m, 0.0),
            type(type_in),
            isvector(n == 1 || m == 1 ? true : false),
            issquare(true)
        {
            switch (type)
            {
            case standard:
                {
                    throw std::invalid_argument("EXCEPTION: Do not use the special matrix type constructor to build a standard matrix! Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
                    break;
                }
            case identity:
                {
                    construct_identity_matrix();
                    break;
                }
            case diagonal:
            {
                this->assign_diagonal(input_data);
                break;
            }
            default: //currently this is all of the rotation matrices
                {
                    construct_rotation_matrix(input_data);
                    break;
                }
            }
        }
        
        template <class T> Matrix<T>::Matrix(const size_t& n_in, const std::vector<T>& input_data, const MatrixType& type_in) :
            n(n_in),
            m(n_in),
            values(n * m, 0.0),
            type(type_in),
            isvector(n == 1 || m == 1 ? true : false),
            issquare(true)
        {
            switch (type)
            {
                case diagonal:
                {
                    this->assign_diagonal(input_data);
                    break;
                }
                default: //currently this is all of the rotation matrices
                {
                    throw std::invalid_argument("I do not recognize this matrix type when passing in a vector of input data: " + std::to_string(type_in) + ". Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
                }
            }
        }

        template <class T> Matrix<T>::Matrix(const size_t& n_in, const T& input_data, const T& input_data_derivative, const MatrixType& type_in) :
            n(n_in),
            m(n_in),
            values(n * m, 0.0),
            type(type_in),
            isvector(n == 1 || m == 1 ? true : false),
            issquare(true)
        {
            this->construct_rotation_matrix_derivative(input_data, input_data_derivative);
        }

        template <class T> Matrix<T>::Matrix(const size_t& n_in, const MatrixType& type_in) :
            n(n_in),
            m(n_in),
            values(n * m, 0.0),
            type(type_in),
            isvector(n == 1 || m == 1 ? true : false),
            issquare(true)
        {
            if (type == identity)
                construct_identity_matrix();
            else
            {
                throw std::invalid_argument("EXCEPTION: identity constructor may only be used to construct the identity matrix. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }
        }

        template <class T> Matrix<T>::Matrix(const Matrix<T>& input_vector, const MatrixType& type_in) :
            n(3),
            m(3),
            values(n * m, 0.0),
            type(type_in),
            isvector(false),
            issquare(true)
        {
            switch (type)
            {
            case skewsymmetric:
            {
#ifndef FAST_EMTG_MATRIX
                if (input_vector.get_n() != 3 || input_vector.get_m() != 1)
                    throw std::invalid_arguement("skew-symmetric matrix constructor must be passed a 3-vector. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
#endif
                this->assign_entry(0, 0,  0.0);
                this->assign_entry(0, 1, -input_vector(2));
                this->assign_entry(0, 2,  input_vector(1));
                this->assign_entry(1, 0,  input_vector(2));
                this->assign_entry(1, 1,  0.0);
                this->assign_entry(1, 2, -input_vector(0));
                this->assign_entry(2, 0, -input_vector(1));
                this->assign_entry(2, 1,  input_vector(0));
                this->assign_entry(2, 2,  0.0);

                break;
            }
            default: //currently this is all of the rotation matrices
            {
                throw std::invalid_argument("I do not recognize this matrix type when passing in an input vector and a type: " + std::to_string(type_in) + ". Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }
            }
        }


        //constructor for a rotation matrix about an arbitary axis
        template <class T> Matrix<T>::Matrix(const size_t& n_in, const T& input_data, const MatrixType& type_in, const Matrix<T> rotaxis) : 
            n(n_in), 
            m(n_in),
            values(n * m, 0.0),
            type(type_in),
            isvector(n == 1 || m == 1 ? true : false),
            issquare(true)
        {
            //assign the matrix type
            type = type_in;
            if (this->type != Rahat)
            {
                throw std::invalid_argument("EXCEPTION: Do not use the rotation matrix about arbitrary axis constructor for any other type of matrix. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }

            //get the rotation angle
            T ctheta = cos(input_data);
            T stheta = sin(input_data);
            T ux = rotaxis(0);
            T uy = rotaxis(1);
            T uz = rotaxis(2);

            this->assign_entry(0, 0, ctheta + ux*ux*(1 - ctheta));
            this->assign_entry(0, 1, ux*uy*(1 - ctheta) - uz*stheta);
            this->assign_entry(0, 2, ux*uz*(1 - ctheta) + uy*stheta);
            this->assign_entry(1, 0, uy*ux*(1 - ctheta) + uz*stheta);
            this->assign_entry(1, 1, ctheta + uy*uy*(1 - ctheta));
            this->assign_entry(1, 2, uy*uz*(1 - ctheta) - ux*stheta);
            this->assign_entry(2, 0, uz*ux*(1 - ctheta) - uy*stheta);
            this->assign_entry(2, 1, uz*uy*(1 - ctheta) + ux*stheta);
            this->assign_entry(2, 2, ctheta + uz*uz*(1 - ctheta));
        }


        //***************************************************
        //destructor
        template <class T> Matrix<T>::~Matrix()
        {
        }

        //***************************************************
        //assignment functions
        template<class T> void Matrix<T>::assign_all(const T* input_data)
        {
            for (size_t i = 0; i < n; ++i)
            {
                for (size_t j = 0; j < m; ++j)
                    values[i*m + j] = input_data[i*m + j];
            }
        }

        template<class T> void Matrix<T>::assign_all(const std::vector<T>& input_data)
        {
            if (values.size() == input_data.size())
                values = input_data;
            else
            {
                throw std::invalid_argument("EXCEPTION: sizes of input data vector and internal values vector must match. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }
        }

        template<class T> void Matrix<T>::assign_entry(const size_t& i, const size_t& j, const T& value)
        {
            if (i >= n || j >= m)
            {
                throw std::invalid_argument("EXCEPTION: Invalid index in [" + std::to_string(i) + ", " + std::to_string(j) + "]. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }
            values[i*m + j] = value;
        }

        template<class T> void Matrix<T>::assign_zeros()
        {
            for (size_t i = 0; i < n; ++i)
            {
                for (size_t j = 0; j < m; ++j)
                    values[i*m + j] = 0;
            }
        }

        template <class T> void Matrix<T>::assign_constant(const T& value)
        {
            for (size_t i = 0; i < n; ++i)
            {
                for (size_t j = 0; j < m; ++j)
                    values[i*m + j] = value;
            }
        }

        template <class T> void Matrix<T>::assign_row(const size_t& i, const T* rowvalues)
        {
            for (size_t j = 0; j < m; ++j)
                values[i*m+j] = rowvalues[j];
        }

        template <class T> void Matrix<T>::assign_row(const size_t& i, const Matrix& TheRow)
        {
            if (!(TheRow.m == m && TheRow.n == 1))
            {
                throw std::invalid_argument("EXCEPTION: Cannot assign row because dimension is incorrect. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }

            for (size_t j = 0; j < m; ++j)
                values[i*m+j] = TheRow(j);
        }

        template <class T> void Matrix<T>::assign_column(const size_t& j, const T* columnvalues)
        {
            for (size_t i = 0; i < n; ++i)
                values[i*m+j] = columnvalues[i];
        }

        template <class T> void Matrix<T>::assign_column(const size_t& j, const Matrix& TheColumn)
        {
            if (!(TheColumn.n == n && TheColumn.m == 1))
            {
                throw std::invalid_argument("EXCEPTION: Cannot assign column because dimension is incorrect. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }

            for (size_t i = 0; i < n; ++i)
                values[i*m+j] = TheColumn(i);
        }

        template <class T> void Matrix<T>::assign_diagonal(const std::vector<T>& input_data)
        {
#ifndef FAST_EMTG_MATRIX
            if (!this->issquare)
            {
                throw std::invalid_argument("You cannot assign the diagonal of a non-square matrix. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }
            else if (this->n != input_data.size())
            {
                throw std::invalid_argument("You have attempted to assign the diagonal values of a matrix but your diagonal vector is not the same dimension as the matrix. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }
#endif

            for (size_t i = 0; i < this->n; ++i)
            {
                this->values[i * this->m] = input_data[i];
            }
        }

        template <class T> void Matrix<T>::assign_diagonal(const T& input_value)
        {
#ifndef FAST_EMTG_MATRIX
            if (!this->issquare)
            {
                throw std::invalid_argument("You cannot assign the diagonal of a non-square matrix. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }
            else if (this->n != input_data.size())
            {
                throw std::invalid_argument("You have attempted to assign the diagonal values of a matrix but your diagonal vector is not the same dimension as the matrix. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }
#endif

            for (size_t i = 0; i < this->n; ++i)
            {
                this->values[i * this->m] = input_value;
            }
        }

        template <class T> void Matrix<T>::construct_identity_matrix()
        {
            for (size_t i = 0; i < n; ++i)
            {
                for (size_t j = 0; j < m; ++j)
                {
                    values[i*m+j] = (i == j) ? 1 : 0;
                }
            }
        }

        template <class T> void Matrix<T>::construct_rotation_matrix(const T& angle)
        {
            //pre-compute trig values
            T cangle = cos(angle);
            T sangle = sin(angle);

            switch (type)
            {
                case Rxhat:
                {
                    assign_entry(0, 0,  1.0);
                    assign_entry(0, 1,  0.0);
                    assign_entry(0, 2,  0.0);
                    assign_entry(1, 0,  0.0);
                    assign_entry(1, 1,  cangle);
                    assign_entry(1, 2, -sangle);
                    assign_entry(2, 0,  0);
                    assign_entry(2, 1,  sangle);
                    assign_entry(2, 2,  cangle);
                    break;
                }
                case Ryhat:
                {
                    assign_entry(0, 0,  cangle);
                    assign_entry(0, 1,  0.0);
                    assign_entry(0, 2,  sangle);
                    assign_entry(1, 0,  0.0);
                    assign_entry(1, 1,  1.0);
                    assign_entry(1, 2,  0.0);
                    assign_entry(2, 0, -sangle);
                    assign_entry(2, 1,  0.0);
                    assign_entry(2, 2,  cangle);
                    break;
                }
                case Rzhat:
                {
                    assign_entry(0, 0,  cangle);
                    assign_entry(0, 1, -sangle);
                    assign_entry(0, 2,  0.0);
                    assign_entry(1, 0,  sangle);
                    assign_entry(1, 1,  cangle);
                    assign_entry(1, 2,  0.0);
                    assign_entry(2, 0,  0.0);
                    assign_entry(2, 1,  0.0);
                    assign_entry(2, 2,  1.0);
                    break;
                }
                default:
                    std::cout<<"type "<<type<<" not handled"<< std::endl;
                    break;            
            }
        }

        template<class T> void Matrix<T>::construct_rotation_matrix_derivative(const T& angle, const T& angledot)
        {
            //pre-compute trig values
            T cangle = cos(angle);
            T sangle = sin(angle);

            switch (type)
            {
                case Rxhatdot:
                {
                    assign_entry(0, 0, 0.0);
                    assign_entry(0, 1, 0.0);
                    assign_entry(0, 2, 0.0);
                    assign_entry(1, 0, 0.0);
                    assign_entry(1, 1, -sangle * angledot);
                    assign_entry(1, 2, -cangle * angledot);
                    assign_entry(2, 0, 0);
                    assign_entry(2, 1, cangle * angledot);
                    assign_entry(2, 2, -sangle * angledot);
                    break;
                }
                case Ryhatdot:
                {
                    assign_entry(0, 0, -sangle * angledot);
                    assign_entry(0, 1, 0.0);
                    assign_entry(0, 2, -cangle * angledot);
                    assign_entry(1, 0, 0.0);
                    assign_entry(1, 1, 0.0);
                    assign_entry(1, 2, 0.0);
                    assign_entry(2, 0, cangle * angledot);
                    assign_entry(2, 1, 0.0);
                    assign_entry(2, 2, -sangle * angledot);
                    break;
                }
                case Rzhatdot:
                {
                    assign_entry(0, 0, -sangle * angledot);
                    assign_entry(0, 1, -cangle * angledot);
                    assign_entry(0, 2, 0.0);
                    assign_entry(1, 0, cangle * angledot);
                    assign_entry(1, 1, -sangle * angledot);
                    assign_entry(1, 2, 0.0);
                    assign_entry(2, 0, 0.0);
                    assign_entry(2, 1, 0.0);
                    assign_entry(2, 2, 0.0);
                    break;
                }
                default:
                    std::cout << "type " << type << " not handled" << std::endl;
                    break;
            }
        }

        template<class T> void Matrix<T>::clear()
        {
            n = 1;
            m = 1;

            issquare = true;
            isvector = true;

            values.resize(1, 0.0);
        }
    
        template<class T> void Matrix<T>::resize(const size_t& n_in, const size_t& m_in)
        {
            n = n_in;
            m = m_in;

            issquare = (n == m ? true : false);
            isvector = (n == 1 || m == 1 ? true : false);

            values.resize(n * m, 0.0);
        }

        template<class T> void Matrix<T>::resize(const size_t& n_in, const size_t& m_in, const T& value)
        {
            n = n_in;
            m = m_in;

            issquare = (n == m ? true : false);
            isvector = (n == 1 || m == 1 ? true : false);

            values.resize(n * m, value);
        }

        template<class T> void Matrix<T>::resize(const size_t& n_in, const size_t& m_in, const MatrixType& type_in)
        {
            this->resize(n_in, m_in);

            this->type = type_in;

            switch (type_in)
            {
                case standard:
                {
                    throw std::invalid_argument("EXCEPTION: Do not use the special matrix type constructor to build a standard matrix. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
                    break;
                }
                case identity:
                {
                    construct_identity_matrix();
                    break;
                }
                default: //currently this is all of the rotation matrices
                {
                    throw std::invalid_argument("EXCEPTION: You can't resize to a rotation matrix. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
                }
            }
        }

        template<class T> Matrix<T>& Matrix<T>::operator= (const T& scalar)
        {
            for (size_t i = 0; i < n; ++i)
            {
                for (size_t j = 0; j < m; ++j)
                    values[i*m+j] = scalar;
            }

            return *this;
        }
    
        template<class T> Matrix<T>& Matrix<T>::operator= (const Matrix& OtherMatrix)
        {
            n = OtherMatrix.n;
            m = OtherMatrix.m;

            values = OtherMatrix.values;
            type = OtherMatrix.type;
            issquare = OtherMatrix.issquare;
            isvector = OtherMatrix.isvector;

            return *this;
        }

        template<class T> void Matrix<T>::shallow_copy(const Matrix& OtherMatrix)
        {
    #ifndef FAST_EMTG_MATRIX
            if (!(this->n == OtherMatrix.n) || !(this->m == OtherMatrix.m))
            {
                throw std::invalid_argument("ERROR: Matrix::shallow_copy() requires matrices of the same dimensions. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }
    #endif
            size_t size = this->n * this->m;

            for (size_t k = 0; k < size; ++k)
                this->values[k] = OtherMatrix.values[k];
        }


        template<class T> void Matrix<T>::shallow_copy(const Matrix& OtherMatrix, const size_t& num_entries)
        {
            for (size_t k = 0; k < num_entries; ++k)
                this->values[k] = OtherMatrix.values[k];
        }

        //***************************************************
        //get methods

        template<class T> Matrix<T> Matrix<T>::getSubMatrix1D(const size_t& firstEntry, const size_t& lastEntry) const
        {
            size_t num_entries = lastEntry - firstEntry + 1;

            Matrix<T> TheSubMatrix(num_entries, 1);

            for (size_t entryIndex = 0; entryIndex < num_entries; ++entryIndex)
                TheSubMatrix(entryIndex) = operator()(entryIndex + firstEntry);

            return TheSubMatrix;
        }

        template<class T> Matrix<T> Matrix<T>::getrow(const size_t& i) const
        {
            Matrix<T> TheRow(1, m);

            for (size_t j = 0; j < m; ++j)
                TheRow(0, j) = getentry(i, j);

            return TheRow;
        }

        template<class T> Matrix<T> Matrix<T>::getcolumn(const size_t& j) const
        {
            Matrix<T> TheColumn(n, 1);

            for (size_t i = 0; i < n; ++i)
                TheColumn(i, 0) = getentry(i, j);

            return TheColumn;
        }

        template<class T> T Matrix<T>::operator() (const size_t& i, const size_t& j) const
        {
            return getentry(i, j);
        }

        template<class T> T Matrix<T>::operator() (const size_t& i) const
        {
#ifndef FAST_EMTG_MATRIX
            if (!isvector)
            {
                throw std::invalid_argument("EXCEPTION: Cannot access a non-vector matrix using only one index. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }
#endif

            return m == 1 ? getentry(i, 0) : getentry(0, i);
        }
        
        template<class T> T& Matrix<T>::operator() (const size_t& i, const size_t& j)
        {
    #ifdef FAST_EMTG_MATRIX
            return *(this->values.data() + i * this->m + j);
    #else
            return getreference(i, j);
    #endif
        }

        template<class T> T& Matrix<T>::operator() (const size_t& i)
        {
#ifndef FAST_EMTG_MATRIX
            if (!isvector)
            {
                throw std::invalid_argument("EXCEPTION: Cannot access a non-vector matrix using only one index. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }
#endif

            return m == 1 ? getreference(i, 0) : getreference(0, i);
        }

        //***************************************************
        //print methods
        template<class T> void Matrix<T>::print_to_screen()
        {
            std::cout.width(15);
            std::cout.precision(20);

            for (size_t i = 0; i < n; ++i)
            {
                std::cout.precision(20);
                std::cout << values[i*m];
                for (size_t j = 1; j < m; ++j)
                {
                    std::cout.width(25);
                    std::cout.precision(20);
                    std::cout << values[i*m+j];
                }
                std::cout << std::endl;
            }
        }

        template<class T> void Matrix<T>::print_to_file()
        {
            this->print_to_file("Matrix.txt");
        }

        template<class T> void Matrix<T>::print_to_file(const std::string& filename, const bool newFile)
        {
            std::ofstream outputfile;
            if (newFile)
                outputfile = std::ofstream(filename.c_str(), std::ios::trunc);
            else
                outputfile = std::ofstream(filename.c_str(), std::ios::app);

            for (size_t i = 0; i < n; ++i)
            {
                outputfile.width(30);
                outputfile.precision(20);
                outputfile << values[i*m];
                for (size_t j = 1; j < m; ++j)
                {
                    outputfile.width(30);
                    outputfile.precision(20);
                    outputfile << values[i*m+j];
                }
                outputfile << std::endl;
            }

            outputfile.close();
        }

        template<class T> void Matrix<T>::print_to_file(const std::string& filename, const std::vector< std::string >& columnheader, const bool newFile)
        {
            std::ofstream outputfile;
            if (newFile)
                outputfile = std::ofstream(filename.c_str(), std::ios::trunc);
            else
                outputfile = std::ofstream(filename.c_str(), std::ios::app);
        
            for (size_t j = 0; j < m; ++j)
            {
                outputfile.width(30);
                outputfile.precision(20);
                if (j < columnheader.size())
                    outputfile << columnheader[j];
                else
                    outputfile << "N/A";
            }
            outputfile << std::endl;
            for (size_t i = 0; i < n; ++i)
            {
                outputfile.width(30);
                outputfile.precision(20);
                outputfile << values[i*m];
                for (size_t j = 1; j < m; ++j)
                {
                    outputfile.width(30);
                    outputfile.precision(20);
                    outputfile << values[i*m + j];
                }
                outputfile << std::endl;
            }

            outputfile.close();
        }
        template<class T> void Matrix<T>::print_to_file(const std::string& filename, const std::vector< std::string >& columnheader, const std::vector< std::string >& rowheader, const bool newFile)
        {
            std::ofstream outputfile;
            if (newFile)
                outputfile = std::ofstream(filename.c_str(), std::ios::trunc);
            else
                outputfile = std::ofstream(filename.c_str(), std::ios::app);

            outputfile.width(30);
            outputfile.precision(20);
            outputfile << std::left << "rows";
            for (size_t j = 0; j < m; ++j)
            {
                outputfile.width(30);
                outputfile.precision(20);
                if (j < columnheader.size())
                    outputfile << columnheader[j];
                else
                    outputfile << "N/A";
            }
            outputfile << std::endl;
            for (size_t i = 0; i < n; ++i)
            {
                outputfile.width(30);
                if (i < rowheader.size())
                    outputfile << std::left << rowheader[i];
                else
                    outputfile << std::left << "N/A";
                for (size_t j = 0; j < m; ++j)
                {
                    outputfile.width(30);
                    outputfile.precision(20);
                    outputfile << values[i*m + j];
                }
                outputfile << std::endl;
            }

            outputfile.close();
        }
        //***************************************************
        //basic math operations

        //sign change
        template<class T> Matrix<T> Matrix<T>::operator- ()
        {
            math::Matrix<T> newMatrix = *this;

            for (size_t i = 0; i < n; ++i)
            {
                for (size_t j = 0; j < m; ++j)
                    newMatrix.values[i * m + j] *= -1;
            }

            return newMatrix;
        }

        //***************************************************
        //matrix addition
        template<class T> Matrix<T> Matrix<T>::operator+ (const Matrix& OtherMatrix) const
        {
            //first check to see if the matrices are the same size
            if (!(n == OtherMatrix.n) || !(m == OtherMatrix.m))
            {
                throw std::invalid_argument("EXCEPTION: Matrix sizes do not match (operator+). Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));;
            }

            Matrix<T> NewMatrix(n, m);

            for (size_t i = 0; i < n; ++i)
            {
                for (size_t j = 0; j < m; ++j)
                {
    #ifdef FAST_EMTG_MATRIX
                    *(NewMatrix.values.data() + i * NewMatrix.m + j) = *(this->values.data() + i * this->m + j) + *(OtherMatrix.values.data() + i * OtherMatrix.m + j);
    #else
                    NewMatrix.assign_entry(i, j, getentry(i, j) + OtherMatrix.getentry(i, j));
    #endif
                }
            }

            return NewMatrix;
        }

        template<class T> Matrix<T> Matrix<T>::operator+ (const T& scalar) const
        {
            Matrix<T> NewMatrix(n, m, values);

            NewMatrix += scalar;

            return NewMatrix;
        }

        template<class T> Matrix<T>& Matrix<T>::operator+= (const Matrix& OtherMatrix)
        {
            //first check to see if the matrices are the same size
            if (!(n == OtherMatrix.n) || !(m == OtherMatrix.m))
            {
                throw std::invalid_argument("EXCEPTION: Matrix sizes do not match (operator+=). Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }

            for (size_t i = 0; i < n; ++i)
            {
                for (size_t j = 0; j < m; ++j)
                {
                    *(this->values.data() + i * this->m + j) += *(OtherMatrix.values.data() + i * m + j);
                }
            }

            return *this;
        }

        template<class T> Matrix<T>& Matrix<T>::operator+= (const T& scalar)
        {
            for (size_t i = 0; i < n; ++i)
            {
                for (size_t j = 0; j < m; ++j)
                {
                    *(this->values.data() + i * this->m + j) += scalar;
                }
            }

            return *this;
        }

        //***************************************************
        //matrix subtraction
        template<class T> Matrix<T> Matrix<T>::operator- (const Matrix& OtherMatrix) const
        {
            //first check to see if the matrices are the same size
            if (!(n == OtherMatrix.n) || !(m == OtherMatrix.m))
            {
                throw std::invalid_argument("EXCEPTION: Matrix sizes do not match (operator-). Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }

            Matrix<T> NewMatrix(n, m);

            for (size_t i = 0; i < n; ++i)
            {
                for (size_t j = 0; j < m; ++j)
                {
    #ifdef FAST_EMTG_MATRIX
                    *(NewMatrix.values.data() + i * NewMatrix.m + j) = *(this->values.data() + i * this->m + j) - *(OtherMatrix.values.data() + i * OtherMatrix.m + j);
    #else
                    NewMatrix.assign_entry(i, j, getentry(i,j) - OtherMatrix.getentry(i,j));
    #endif
                }
            }

            return NewMatrix;
        }

        template<class T> Matrix<T> Matrix<T>::operator- (const T& scalar) const
        {
            Matrix<T> NewMatrix(n, m, values);

            NewMatrix -= scalar;

            return NewMatrix;
        }

        template<class T> Matrix<T>& Matrix<T>::operator-= (const Matrix& OtherMatrix)
        {
            //first check to see if the matrices are the same size
            if (!(n == OtherMatrix.n) || !(m == OtherMatrix.m))
            {
                throw std::invalid_argument("EXCEPTION: Matrix sizes do not match (operator-=). Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }

            for (size_t i = 0; i < n; ++i)
            {
                for (size_t j = 0; j < m; ++j)
                {
                    *(this->values.data() + i * this->m + j) -= *(OtherMatrix.values.data() + i * OtherMatrix.m + j);
                }
            }

            return *this;
        }

        template<class T> Matrix<T>& Matrix<T>::operator-= (const T& scalar)
        {
            for (size_t i = 0; i < n; ++i)
            {
                for (size_t j = 0; j < n; ++j)
                {
                    *(this->values.data() + i * this->m + j) -= scalar;
                }
            }

            return *this;
        }

        //right multiplication (left multiplication is done by using the other object's left multiplier
        template<class T> Matrix<T> Matrix<T>::operator* (const Matrix& OtherMatrix) const
        {
            //first check to see if the matrices are compatible
            if (!(m == OtherMatrix.n))
            {
                throw std::invalid_argument("EXCEPTION: matrix sizes do not match (operator*). Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }

            Matrix<T> NewMatrix(n, OtherMatrix.m);

            for (size_t i = 0; i < NewMatrix.n; ++i)
            {
                for (size_t j = 0; j < NewMatrix.m; ++j)
                {
                    for (size_t k = 0; k < m; ++k)
    #ifdef FAST_EMTG_MATRIX
                        *(NewMatrix.values.data() + i * NewMatrix.m + j) += *(this->values.data() + i * this->m + k) * *(OtherMatrix.values.data() + k * OtherMatrix.m + j);
    #else
                        NewMatrix(i, j) += getentry(i, k) * OtherMatrix(k, j);
    #endif
                }
            }


            return NewMatrix;
        }

        template<class T> Matrix<T> Matrix<T>::operator* (const T& scalar) const
        {
            Matrix<T> NewMatrix(n, m, values);

            NewMatrix *= scalar;

            return NewMatrix;
        }

        template<class T> Matrix<T>& Matrix<T>::operator*= (const Matrix& OtherMatrix)
        {
            Matrix<T> NewMatrix = *this * OtherMatrix;

            assign_all(NewMatrix.values);

            return *this;
        }
    
        template<class T> Matrix<T>& Matrix<T>::operator*= (const T& scalar)
        {
            for (size_t i = 0; i < n; ++i)
            {
                for (size_t j = 0; j < m; ++j)
                    *(values.data() + i*m + j) *= scalar;
            }

            return *this;
        }

        //comparators
        template<class T> bool Matrix<T>::operator== (const Matrix& OtherMatrix) const
        {
            if (!(m == OtherMatrix.m))
                return false;

            for (size_t i = 0; i < n; ++i)
            {
                for (size_t j = 0; j < m; ++j)
                {
                    if (!(*(this->values.data() + i*m+j) == *(OtherMatrix.values.data() + i*OtherMatrix.m+j)))
                        return false;
                }
            }
        
            return true;
        }

        template<class T> bool Matrix<T>::operator!= (const Matrix& OtherMatrix) const
        {
            return !(*this == OtherMatrix);
        }

    #ifdef _EMTG_thruth_table
        template<class T> Matrix<bool> Matrix<T>::operator<= (const Matrix& OtherMatrix)
        {
            //first check to see if the matrices are compatible
            if (!(n == OtherMatrix.n && m == OtherMatrix.m))
            {
                throw std::invalid_argument("EXCEPTION: matrix sizes do not match (operator<=). Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }

            Matrix<bool> TruthTable(n, m);
        
            for (size_t i = 0; i < n; ++i)
            {
                for (size_t j = 0; j < m; ++j)
                    TruthTable(i,j) = values[i*m+j] <= OtherMatrix.values[i*m+j] ? true : false;
            }

            return TruthTable;
        }

        template<class T> Matrix<bool> Matrix<T>::operator>= (const Matrix& OtherMatrix)
        {
            //first check to see if the matrices are compatible
            if (!(n == OtherMatrix.n && m == OtherMatrix.m))
            {
                throw std::invalid_argument("EXCEPTION: matrix sizes do not match (operator>=). Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }

            Matrix<bool> TruthTable(n, m);

            for (size_t i = 0; i < n; ++i)
            {
                for (size_t j = 0; j < m; ++j)
                    TruthTable(i,j) = values[i*m+j] >= OtherMatrix.values[i*m+j] ? true : false;
            }

            return TruthTable;
        }

        template<class T> Matrix<bool> Matrix<T>::operator< (const Matrix& OtherMatrix)
            {
            //first check to see if the matrices are compatible
            if (!(n == OtherMatrix.n && m == OtherMatrix.m))
            {
                throw std::invalid_argument("EXCEPTION: matrix sizes do not match (operator<). Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }

            Matrix<bool> TruthTable(n, m);

            for (size_t i = 0; i < n; ++i)
            {
                for (size_t j = 0; j < m; ++j)
                    TruthTable(i,j) = values[i*m+j] < OtherMatrix.values[i*m+j] ? true : false;
            }

            return TruthTable;
        }

        template<class T> Matrix<bool> Matrix<T>::operator> (const Matrix& OtherMatrix)
        {
            //first check to see if the matrices are compatible
            if (!(n == OtherMatrix.n && m == OtherMatrix.m))
            {
                throw std::invalid_argument("EXCEPTION: matrix sizes do not match (operator>). Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }

            Matrix<bool> TruthTable(n, m);

            for (size_t i = 0; i < n; ++i)
            {
                for (size_t j = 0; j < m; ++j)
                    TruthTable(i,j) = values[i*m+j] > OtherMatrix.values[i*m+j] ? true : false;
            }

            return TruthTable;
        }
    #endif //_EMTG_thruth_table

        //multiplication(element-wise)
        template<class T> Matrix<T> Matrix<T>::element_multiply(const Matrix& OtherMatrix) const
        {
            //first check to make sure that the operation can be done
            //is the number of rows the same?
            if (n == OtherMatrix.n)
            {
                //is the number of columns the same?
                if (m == OtherMatrix.m)
                {
                    //the number of rows and columns is the same, so operate on the entire matrix
                    Matrix<T> NewMatrix(n,m);

                    for (size_t i = 0; i < n; ++i)
                    {
                        for (size_t j = 0; j < m; ++j)
    #ifdef FAST_EMTG_MATRIX
                            *(NewMatrix.values.data() + i * NewMatrix.m + j) = *(this->values.data() + i * this->m + j) * *(OtherMatrix.values.data() + i * OtherMatrix.m + j);
    #else
                            NewMatrix(i,j) = values[i*m+j] * OtherMatrix(i,j);
    #endif
                    }

                    return NewMatrix;
                }
                else if (OtherMatrix.m == 1)
                {
                    //we actually want to operate through by a column vector
                    Matrix<T> NewMatrix(n,m);
                
                    for (size_t i = 0; i < n; ++i)
                    {
                        for (size_t j = 0; j < m; ++j)
    #ifdef FAST_EMTG_MATRIX
                            *(NewMatrix.values.data() + i * NewMatrix.m + j) = *(this->values.data() + i * this->m + j) * *(OtherMatrix.values.data() + i);
    #else
                            NewMatrix(i,j) = values[i*m+j] * OtherMatrix.values[i];
    #endif
                    }

                    return NewMatrix;
                }
                else
                {
                    //can't do this
                    throw std::invalid_argument("EXCEPTION: matrix sizes do not match (element multiply). Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
                }
            }
            else if (m == OtherMatrix.m)
            {
                if (OtherMatrix.n == 1)
                {
                    //we actually want to operate through by a row vector
                    Matrix<T> NewMatrix(n,m);
                
                    for (size_t i = 0; i < n; ++i)
                    {
                        for (size_t j = 0; j < m; ++j)
                            NewMatrix(i,j) = values[i*m+j] * OtherMatrix.values[j];
                    }

                    return NewMatrix;
                }
                else
                {
                    //can't do this
                    throw std::invalid_argument("EXCEPTION: matrix sizes do not match (element multiply). Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
                }
            }
            else
            {
                //can't do this
                throw std::invalid_argument("EXCEPTION: matrix sizes do not match (element multiply). Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }
        }
    
        template<class T> Matrix<T>& Matrix<T>::element_multiply_equals (const Matrix& OtherMatrix)
        {
            //first check to make sure that the operation can be done
            //is the number of rows the same?
            if (n == OtherMatrix.n)
            {
                //is the number of columns the same?
                if (m == OtherMatrix.m)
                {
                    //the number of rows and columns is the same, so operate on the entire matrix

                    for (size_t i = 0; i < n; ++i)
                    {
                        for (size_t j = 0; j < m; ++j)
    #ifdef FAST_EMTG_MATRIX
                            *(this->values.data() + i * this->m + j) *= *(OtherMatrix.values.data() + i * OtherMatrix.m + j);
    #else
                            values[i*m+j] *= OtherMatrix(i,j);
    #endif
                    }

                    return *this;
                }
                else if (OtherMatrix.m == 1)
                {
                    //we actually want to operate through by a column vector

                    for (size_t i = 0; i < n; ++i)
                    {
                        for (size_t j = 0; j < m; ++j)
    #ifdef FAST_EMTG_MATRIX
                            *(this->values.data() + i * this->m + j) *= *(OtherMatrix.values.data() + i);
    #else
                            values[i*m+j] *= OtherMatrix.values[i];
    #endif
                    }
                
                    return *this;
                }
                else
                {
                    //can't do this
                    throw std::invalid_argument("EXCEPTION: matrix sizes do not match (element multiply-equals). Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
                }
            }
            else if (m == OtherMatrix.m)
            {
                if (OtherMatrix.n == 1)
                {
                    //we actually want to operate through by a row vector
                    for (size_t i = 0; i < n; ++i)
                    {
                        for (size_t j = 0; j < m; ++j)
                            *(this->values.data() + i*this->m + j) *= *(OtherMatrix.values.data() + j);
                    }

                    return *this;
                }
                else
                {
                    //can't do this
                    throw std::invalid_argument("EXCEPTION: matrix sizes do not match (element multiply-equals). Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
                }
            }
            else
            {
                //can't do this
                throw std::invalid_argument("EXCEPTION: matrix sizes do not match (element multiply-equals). Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }
        }

        //division (element-wise)
        template<class T> Matrix<T> Matrix<T>::operator/ (const T& scalar) const
        {
            Matrix<T> NewMatrix(n,m, values);

            NewMatrix /= scalar;

            return NewMatrix;
        }
    
        template<class T> Matrix<T> Matrix<T>::element_divide(const Matrix& OtherMatrix) const
        {
            //first check to make sure that the operation can be done
            //is the number of rows the same?
            if (n == OtherMatrix.n)
            {
                //is the number of columns the same?
                if (m == OtherMatrix.m)
                {
                    //the number of rows and columns is the same, so operate on the entire matrix
                    Matrix<T> NewMatrix(n,m);

                    for (size_t i = 0; i < n; ++i)
                    {
                        for (size_t j = 0; j < m; ++j)
                            NewMatrix(i,j) = values[i*m+j] / OtherMatrix(i,j);
                    }

                    return NewMatrix;
                }
                else if (OtherMatrix.m == 1)
                {
                    //we actually want to operate through by a column vector
                    Matrix<T> NewMatrix(n,m);
                
                    for (size_t i = 0; i < n; ++i)
                    {
                        for (size_t j = 0; j < m; ++j)
                            NewMatrix(i,j) = values[i*m+j] / OtherMatrix.values[i];
                    }

                    return NewMatrix;
                }
                else
                {
                    //can't do this
                    throw std::invalid_argument("EXCEPTION: matrix sizes do not match (element divide). Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
                }
            }
            else if (m == OtherMatrix.m)
            {
                if (OtherMatrix.n == 1)
                {
                    //we actually want to operate through by a row vector
                    Matrix<T> NewMatrix(n,m);
                
                    for (size_t i = 0; i < n; ++i)
                    {
                        for (size_t j = 0; j < m; ++j)
                            NewMatrix(i,j) = values[i*m+j] / OtherMatrix.values[j];
                    }

                    return NewMatrix;
                }
                else
                {
                    //can't do this
                    throw std::invalid_argument("EXCEPTION: matrix sizes do not match (element divide). Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));;
                }
            }
            else
            {
                //can't do this
                throw std::invalid_argument("EXCEPTION: matrix sizes do not match (element divide). Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }
        }

        template<class T> Matrix<T>& Matrix<T>::operator/= (const T& scalar)
        {
            T inverse_scalar = 1 / scalar;

            for (size_t i = 0; i < n; ++i)
            {
                for (size_t j = 0; j < m; ++j)
                    *(this->values.data() + i * this->m + j) *= inverse_scalar;
            }

            return *this;
        }

        template<class T> Matrix<T>& Matrix<T>::element_divide_equals (const Matrix& OtherMatrix)
        {
            //first check to make sure that the operation can be done
            //is the number of rows the same?
            if (n == OtherMatrix.n)
            {
                //is the number of columns the same?
                if (m == OtherMatrix.m)
                {
                    //the number of rows and columns is the same, so operate on the entire matrix
                    Matrix<T> NewMatrix(n,m);

                    for (size_t i = 0; i < n; ++i)
                    {
                        for (size_t j = 0; j < m; ++j)
                        {
                            values[i*m + j] /= OtherMatrix(i, j);
                            NewMatrix(i, j) = values[i*m + j];
                        }
                    }

                    return *this;
                }
                else if (OtherMatrix.m == 1)
                {
                    //we actually want to operate through by a column vector
                    Matrix<T> NewMatrix(n,m);
                
                    for (size_t i = 0; i < n; ++i)
                    {
                        for (size_t j = 0; j < m; ++j)
                        {
                            values[i*m + j] /= OtherMatrix.values[i];
                            NewMatrix(i, j) = values[i*m + j];
                        }
                    }

                    return *this;
                }
                else
                {
                    //can't do this
                    throw std::invalid_argument("EXCEPTION: matrix sizes do not match (element divide-equals). Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
                }
            }
            else if (m == OtherMatrix.m)
            {
                if (OtherMatrix.n == 1)
                {
                    //we actually want to operate through by a row vector
                    Matrix<T> NewMatrix(n,m);
                
                    for (size_t i = 0; i < n; ++i)
                    {
                        for (size_t j = 0; j < m; ++j)
                        {
                            values[i*m + j] /= OtherMatrix.values[j];
                            NewMatrix(i, j) = values[i*m + j];
                        }
                    }

                    return *this;
                }
                else
                {
                    //can't do this
                    throw std::invalid_argument("EXCEPTION: matrix sizes do not match (element divide-equals). Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
                }
            }
            else
            {
                //can't do this
                throw std::invalid_argument("EXCEPTION: matrix sizes do not match (element divide-equals). Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }
        }

        //*****************************************************
        //Special Matrix Math Functions
        template<class T> Matrix<T> Matrix<T>::transpose() const
        {
            Matrix<T> NewMatrix(m, n);

            for (size_t i = 0; i < n; ++i)
            {
                for (size_t j = 0; j < m; ++j)
                    NewMatrix(j,i) = getentry(i,j);
            }
        
            return NewMatrix;
        }


        template<class T> T Matrix<T>::determinant() const
        {
            return this->determinant(this->n);
        }

        template<class T> Matrix<T> Matrix<T>::Cofactor(const size_t& p, const size_t& q, const size_t& n) const
        {
            size_t i = 0, j = 0;
            Matrix<T> cofactor(n, n, 0.0);

            // Looping for each element of the matrix
            for (size_t row = 0; row < n; ++row)
            {
                for (size_t col = 0; col < n; ++col)
                {
                    //  Copying into temporary matrix only those element
                    //  which are not in given row and column
                    if (row != p && col != q)
                    {
                        cofactor(i, j++) = this->getentry(row, col);

                        // Row is filled, so increase row index and
                        // reset col index
                        if (j == n - 1)
                        {
                            j = 0;
                            ++i;
                        }
                    }
                }
            }

            return cofactor;
        }//end Cofactor()

        template<class T> T Matrix<T>::determinant(const size_t& n) const
        {
            if (!issquare)
            {
                throw std::invalid_argument("EXCEPTION: Cannot take the determinant of a non-square matrix. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }

            if (n == 1)
            {
                //the determinant of a 1x1 is trivial

                return values[0];
            }
            else if (n == 2)
            {
                //the determinant of a 2x2 is ad - bc

                return values[0] * values[3] - values[1] * values[2];
            }
            else if (n == 3)
            {
                //the determinant of a 3x3 is aei + bfg + cdh - ceg - bdi - afh

                return values[0]*values[4]*values[8] + values[1]*values[5]*values[6] + values[2]*values[3]*values[7] - values[2]*values[4]*values[6] - values[1]*values[3]*values[8] - values[0]*values[5]*values[7];
            }
            else
            {
                T D = 0.0;
                Matrix<T> cofactor;
                int sign = 1;

                for (size_t f = 0; f < n; ++f)
                {
                    cofactor = this->Cofactor(0, f, n);
                    D += sign * this->getentry(0, f) + cofactor.determinant(n - 1);
                    sign = -sign; //terms have alternate sign
                }
                
                return D;
            }
        }//end determinant()

        template<class T> Matrix<T> Matrix<T>::adjoint(const size_t& n) const
        {
            Matrix<T> adjoint(n, n, 0.0);
            Matrix<T> cofactor;

            if (n == 1)
            {
                adjoint(0, 0) = 1;
                return adjoint;
            }

            int sign = 1;

            for (size_t i = 0; i < n; i++)
            {
                for (size_t j = 0; j < n; j++)
                {
                    // Get cofactor of A[i][j]
                    cofactor = this->Cofactor(i, j, n);

                    // sign of adjoint[j][i] positive if sum of row
                    // and column indexes is even.
                    sign = ((i + j) % 2 == 0) ? 1 : -1;

                    // Interchanging rows and columns to get the
                    // transpose of the cofactor matrix
                    adjoint(j, i) = sign * cofactor.determinant(n - 1);
                }
            }

            return adjoint;
        }//end adjoint()

        template<class T> Matrix<T> Matrix<T>::inverse() const
        {
            //can only invert square matrices
            if (!issquare)
            {
                throw std::invalid_argument("EXCEPTTION: Non-square matrices cannot be inverted.");
            }

            switch (type)
            {
                case standard:
                {
                    //first check if the matrix is invertible, i.e is the determinant nonzero?
                    T TheDeterminant = determinant();
                    if (TheDeterminant == 0)
                    {
                        throw std::invalid_argument("EXCEPTION: Cannot invert matrix because determinant is zero. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
                    }
                    //create a matrix to hold the inverse
                    Matrix<T> TheInverse(n, m);

                    switch (n)
                    {
                        case 1:
                        {
                            TheInverse(0,0) = values[0];
                            break;
                        }
                        case 2:
                        {
                            TheInverse(0,0) = values[3];
                            TheInverse(0,1) = -values[1];
                            TheInverse(1,0) = -values[2];
                            TheInverse(1,1) = values[0];
                            TheInverse /= TheDeterminant;
                            break;
                        }
                        case 3:
                        {
                            TheInverse(0,0) = values[4]*values[8] - values[5]*values[7];
                            TheInverse(0,1) = values[5]*values[6] - values[3]*values[8];
                            TheInverse(0,2) = values[3]*values[7] - values[4]*values[6];
                            TheInverse(1,0) = values[2]*values[7] - values[1]*values[8];
                            TheInverse(1,1) = values[0]*values[8] - values[2]*values[6];
                            TheInverse(1,2) = values[6]*values[1] - values[0]*values[7];
                            TheInverse(2,0) = values[1]*values[5] - values[2]*values[4];
                            TheInverse(2,1) = values[2]*values[3] - values[0]*values[5];
                            TheInverse(2,2) = values[0]*values[4] - values[1]*values[3];
                            TheInverse /= TheDeterminant;
                            break;
                        }
                        default: //n x n
                        {
                            T det = this->determinant();
                            if (det < 1.0e-20)
                            {
                                throw std::underflow_error("Singular matrix, can't find its inverse. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
                                break;
                            }

                            // Find adjoint
                            Matrix<T> adjoint = this->adjoint(this->n);

                            // Find Inverse using formula "inverse(A) = adj(A)/det(A)"
                            for (int i = 0; i < this->n; ++i)
                                for (int j = 0; j < this->n; ++j)
                                    TheInverse(i, j) = adjoint(i, j) / det;
                        }
                    }
                
                    return TheInverse.transpose();
                    break;
                }
                case identity:
                {
                    return *this;
                    break;
                }
                case Rxhat:
                {
                    return transpose();
                    break;
                }
                case Ryhat:
                {
                    return transpose();
                    break;
                }
                case Rzhat:
                {
                    return transpose();
                    break;
                }
                case Rahat:
                {
                    return transpose();
                    break;
                }
                default:
                    throw std::invalid_argument("EXCEPTION: Can't invert a matrix of type " + std::to_string(this->type) + ". Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
                    break;
            }

            return *this;
        }

        template<class T> Matrix<T> Matrix<T>::horz_cat(const Matrix& OtherMatrix) const
        {
            if (!(n == OtherMatrix.n))
            {
                throw std::invalid_argument("EXCEPTION: Cannot perform horizontal concatenation because vertical dimension does not match. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }

            Matrix NewMatrix(n, m + OtherMatrix.m);
        
            for (size_t i = 0; i < n; ++i)
            {
                for (size_t j1 = 0; j1 < m; ++j1)
                    NewMatrix(i, j1) = values[i*m+j1];

                for (size_t j2 = 0; j2 < OtherMatrix.m; ++j2)
                    NewMatrix(i,j2 + m) = OtherMatrix(i,j2);
            }

            return NewMatrix;
        }

        template<class T> Matrix<T> Matrix<T>::vert_cat(const Matrix& OtherMatrix) const
        {
            if (!(m == OtherMatrix.m))
            {
                throw std::invalid_argument("EXCEPTION: Cannot perform vertical concatenation because horizontal dimension does not match. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }

            Matrix NewMatrix(n + OtherMatrix.n, m);
        

            for (size_t j = 0; j < m; ++j)
            {
                for (size_t i1 = 0; i1 < n; ++i1)
                    NewMatrix(i1, j) = values[i1*m+j];

                for (size_t i2 = 0; i2 < OtherMatrix.n; ++i2)
                    NewMatrix(i2 + n,j) = OtherMatrix(i2,j);
            }

            return NewMatrix;
        }

        //reshape function; intended to work just like reshape in MATLAB
        template<class T> Matrix<T> Matrix<T>::reshape(const size_t& new_n, const size_t& new_m) const
        {
            if (!(new_n*new_m == n*m))
            {
                throw std::invalid_argument("EXCEPTION: Matrix sizes much match in a reshape. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));;
            }

            Matrix<T> NewMatrix(new_n, new_m, values);

            return NewMatrix;
        }
    
        //*************************************************
        //vector math
        template<class T> T Matrix<T>::norm() const
        {
            return sqrt(dot(*this));        
        }


        template<class T> T Matrix<T>::norm(const size_t& FirstEntry, const size_t& LastEntry) const
        {
            T TheNormSquared = 0.0;

            for (size_t i = FirstEntry; i < LastEntry; ++i)
            {
                T Term = this->operator()(i);
                TheNormSquared += Term * Term;
            }

            return sqrt(TheNormSquared);
        }

        template<class T> T Matrix<T>::hypnorm(const double& hypconstant) const
        {
            return sqrt(dot(*this) + hypconstant * hypconstant);
        }

        template<class T> Matrix<T> Matrix<T>::unitize() const
        {

            return (*this / (this->norm() + 1.0e-20));
        }

        template<class T> T Matrix<T>::dot(const Matrix& OtherMatrix) const
        {
            T TheDot = 0;

            if (!isvector || !OtherMatrix.isvector)
            {
                throw std::invalid_argument("EXCEPTION: Cannot take the norm of a non-vector. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }

            size_t dim = n == 1 ? m : n;
            size_t Otherdim = OtherMatrix.n == 1 ? OtherMatrix.m : OtherMatrix.n;

            if (!(dim == Otherdim))
            {
                throw std::invalid_argument("EXCEPTION: Can only take the dot product of equal length vectors. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }

            for (size_t k = 0; k < dim; ++k)
                TheDot += values[k]*OtherMatrix.values[k];

            return TheDot;    
        }

        template<class T> Matrix<T> Matrix<T>::cross(const Matrix& OtherMatrix) const
        {
            //cross product only works on 3-vectors
            if (!isvector || !OtherMatrix.isvector)
            {
                throw std::invalid_argument("EXCEPTION: Cannot take the cross product of a non-3-vector. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }

            size_t dim = n == 1 ? m : n;
            size_t Otherdim = OtherMatrix.n == 1 ? OtherMatrix.m : OtherMatrix.n;

            if (!(dim == 3 && Otherdim == 3))
            {
                throw std::invalid_argument("EXCEPTION: Can only take cross product of two 3-vectors. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }

            Matrix TheCross(3,1);

            TheCross(0,0) = values[1]*OtherMatrix.values[2] - values[2]*OtherMatrix.values[1];
            TheCross(1,0) = values[2]*OtherMatrix.values[0] - values[0]*OtherMatrix.values[2];
            TheCross(2,0) = values[0]*OtherMatrix.values[1] - values[1]*OtherMatrix.values[0];

            return TheCross;
        }

        template<class T> Matrix<T> Matrix<T>::unitcross(const Matrix& OtherMatrix) const
        {
            Matrix TheUnitCross = cross(OtherMatrix);

            return TheUnitCross / TheUnitCross.norm();
        }

        template<class T> Matrix<T> Matrix<T>::unitDerivative(const Matrix& OtherMatrix) const
        {
            //we make the assumption the OtherMatrix is the non-unit form of the current vector

            T magnitude = OtherMatrix.norm();
            Matrix I = Matrix(this->n, MatrixType::identity);
            Matrix OuterProduct = OtherMatrix * OtherMatrix.transpose();

            Matrix dx_unit_dx = (I - OuterProduct / magnitude / magnitude) / magnitude;

            return dx_unit_dx;
        }

        template<class T> Matrix<T> Matrix<T>::crossDerivative(const Matrix& OtherMatrix) const
        {
            //cross product only works on 3-vectors
            if (!isvector || !OtherMatrix.isvector)
            {
                throw std::invalid_argument("EXCEPTION: Cannot take the cross product of a non-3-vector, let alone the derivatives thereof. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }

            size_t dim = n == 1 ? m : n;
            size_t Otherdim = OtherMatrix.n == 1 ? OtherMatrix.m : OtherMatrix.n;

            if (!(dim == 3 && Otherdim == 3))
            {
                throw std::invalid_argument("EXCEPTION: Can only take the derivative of the cross product of two 3-vectors. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }

            Matrix<T> crossDerivative(6, 3, 0.0);
            crossDerivative(0, 0) = 0.0;
            crossDerivative(0, 1) = -OtherMatrix(2);
            crossDerivative(0, 2) = OtherMatrix(1);
            crossDerivative(1, 0) = OtherMatrix(2);
            crossDerivative(1, 1) = 0.0;
            crossDerivative(1, 2) = -OtherMatrix(0);
            crossDerivative(2, 0) = -OtherMatrix(1);
            crossDerivative(2, 1) = OtherMatrix(0);
            crossDerivative(2, 2) = 0.0;
            crossDerivative(3, 0) = 0.0;
            crossDerivative(3, 1) = this->values[2];
            crossDerivative(3, 2) = -this->values[1];
            crossDerivative(4, 0) = -this->values[2];
            crossDerivative(4, 1) = 0.0;
            crossDerivative(4, 2) = this->values[0];
            crossDerivative(5, 0) = this->values[1];
            crossDerivative(5, 1) = -this->values[0];
            crossDerivative(5, 2) = 0.0;

            return crossDerivative;
        }

        template<class T> void Matrix<T>::cross_in_place(const Matrix& OtherMatrix, Matrix& TargetMatrix)
        {
            size_t TargetDim = TargetMatrix.n == 1 ? TargetMatrix.m : TargetMatrix.n;
            if (!(TargetDim == 3))
            {
                throw std::invalid_argument("EXCEPTION: Target matrix of a cross product must be a 3-vector. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }

            TargetMatrix.values[0] = values[1]*OtherMatrix.values[2] - values[2]*OtherMatrix.values[1];
            TargetMatrix.values[1] = values[2]*OtherMatrix.values[0] - values[0]*OtherMatrix.values[2];
            TargetMatrix.values[2] = values[0]*OtherMatrix.values[1] - values[1]*OtherMatrix.values[0];
        }

        //left-handed operators

#ifdef AD_INSTRUMENTATION
        Matrix<GSAD::adouble> operator* (const Matrix<double>& otherMatrix, const Matrix<GSAD::adouble> thisMatrix)
        {
            //first check to see if the matrices are compatible
            if (!(thisMatrix.n == otherMatrix.m))
            {
                throw std::invalid_argument("EXCEPTION: matrix sizes do not match (operator*). Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }

            Matrix<GSAD::adouble> NewMatrix(otherMatrix.n, thisMatrix.m);

            for (size_t i = 0; i < NewMatrix.n; ++i)
            {
                for (size_t j = 0; j < NewMatrix.m; ++j)
                {
                    for (size_t k = 0; k < otherMatrix.m; ++k)
                        NewMatrix(i, j) += otherMatrix.getentry(i, k) * thisMatrix(k, j);
                }
            }

            return NewMatrix;
        }
#endif



        //explicit instantiation
#ifdef AD_INSTRUMENTATION
        template class Matrix<GSAD::adouble>;
#endif
        template class Matrix<double>;
    }//close namespace Math
}// close namespace EMTG