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

#include "EMTG_Tensor.h"
#include <exception>

namespace EMTG 
{ 
    namespace math
    {
        //default constructor
        template<class T> Tensor<T>::Tensor() :
            n(1),
            m(1),
            p(1),
            iscube(true),
            values(std::vector< Matrix<T> >(1, Matrix<T>(1, 1)))
        {}

        //constructor for tensor when size is known but there is no input data
        template<class T> Tensor<T>::Tensor(const size_t& n, const size_t& m, const size_t& p) :
            n(n),
            m(m),
            p(p),
            values(std::vector< Matrix<T> >(p, Matrix<T>(n, m, 0.0)))
        {}

        //constructor for tensor when size is known but there is input data
        template<class T> Tensor<T>::Tensor(const size_t& n, const size_t& m, const size_t& p, std::vector< Matrix<T> > input_data) :
            n(n),
            m(m),
            p(p),
            values(input_data)
        {
#ifndef FAST_EMTG_MATRIX
            if (input_data.size() != this->p)
            {
                throw std::invalid_arguement("tensor input data number of slices does not match tensor p dimension. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }
            else
            {
                for (Matrix<T> slice : this->values)
                {
                    if (slice.get_n() != this->n)
                    {
                        throw std::invalid_arguement("tensor input data n does not match tensor n dimension. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
                    }
                    else if (slice.get_m() != this->m)
                    {
                        throw std::invalid_arguement("tensor input data m does not match tensor m dimension. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
                    }
                }
            }//end error checking
#endif
            if (n == m && n == p)
                this->iscube = true;
            else
                this->iscube = false;
        }


        //constructor for cubic tensors
        template<class T> Tensor<T>::Tensor(const size_t& n, std::vector< Matrix<T> > input_data) :
            n(n),
            m(n),
            p(n),
            iscube(true),
            values(input_data)
        {
#ifndef FAST_EMTG_MATRIX
            if (input_data.size() != this->p)
            {
                throw std::invalid_arguement("tensor input data number of slices does not match tensor p dimension. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }
            else
            {
                for (Matrix<T> slice : this->values)
                {
                    if (slice.get_n() != this->n)
                    {
                        throw std::invalid_arguement("tensor input data n does not match tensor n dimension. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
                    }
                    else if (slice.get_m() != this->m)
                    {
                        throw std::invalid_arguement("tensor input data m does not match tensor m dimension. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
                    }
                }
            }//end error checking
#endif
            if (n == m && n == p)
                this->iscube = true;
            else
                this->iscube = false;
        }


        template<class T> void Tensor<T>::print_to_screen()
        {
            for (Matrix<T> slice : this->values)
                slice.print_to_screen();
        }

        template<class T> void Tensor<T>::print_to_file()
        {
            this->print_to_file("Tensor.txt");
        }

        template<class T> void Tensor<T>::print_to_file(const std::string& filename)
        {
            bool newFile = true;

            for (Matrix<T> slice : this->values)
            {
                slice.print_to_file(filename, newFile);
                newFile = false;
            }
        }

        template<class T> void Tensor<T>::print_to_file(const std::string& filename, const std::vector< std::string >& rowheader)
        {
            bool newFile = true;

            for (Matrix<T> slice : this->values)
            {
                slice.print_to_file(filename, rowheader, newFile);
                newFile = false;
            }
        }

        template<class T> void Tensor<T>::print_to_file(const std::string& filename, const std::vector< std::string >& rowheader, const std::vector< std::string >& columnheader)
        {
            bool newFile = true;

            for (Matrix<T> slice : this->values)
            {
                slice.print_to_file(filename, rowheader, columnheader, newFile);
                newFile = false;
            }
        }


        //assignment
        template<class T> void Tensor<T>::assign_all(std::vector< Matrix<T> >& input_data)
        {

#ifndef FAST_EMTG_MATRIX
            if (input_data.size() != this->p)
            {
                throw std::invalid_arguement("tensor input data number of slices does not match tensor p dimension. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));;
            }
            else
            {
                for (Matrix<T> slice : this->values)
                {
                    if (slice.get_n() != this->n)
                    {
                        throw std::invalid_arguement("tensor input data n does not match tensor n dimension. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
                    }
                    else if (slice.get_m() != this->m)
                    {
                        throw std::invalid_arguement("tensor input data m does not match tensor m dimension. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
                    }
                }
            }//end error checking
#endif
            this->values = input_data;
        }//end assign_all()

        template<class T> void Tensor<T>::assign_entry(const size_t& i, const size_t& j, const size_t& k, const T& value)
        {
            this->getreference(i, j, k) = value;
        }//end assign_entries()

        template<class T> void Tensor<T>::assign_zeros()
        {
            for (Matrix<T> slice : this->values)
                slice.assign_zeros();
        }//end assign_zeros()


        template<class T> Matrix<T> Tensor<T>::bullet2_vector(const Matrix<T>& vector) const
        {
#ifndef FAST_EMTG_MATRIX
            if (this->iscube != true)
            {
                throw std::invalid_arguement("bullet2_vector can only be called on a cubic tensor. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }
            if (vector.get_n() != this->n)
            {
                throw std::invalid_arguement("bullet2_vector only works if the length of the column vector matches the dimension of the tensor. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }
            if (vector.get_m() != this->1)
            {
                throw std::invalid_arguement("bullet2_vector only works with column vectors. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }
#endif

            Matrix<T> bullet(this->n, this->n, 0.0);

            for (size_t k = 0; k < this->p; ++k)
            {
                bullet.assign_column(k, this->values[k] * vector);
            }

            return bullet;
        }//end bullet2_vector()

        //***************************************************
        //Tensor addition
        template<class T> Tensor<T> Tensor<T>::operator+ (const Tensor& OtherTensor) const
        {
            //first check to see if the tensors are the same size
#ifndef FAST_EMTG_MATRIX
            if (!(n == OtherTensor.n) || !(m == OtherTensor.m) || !(p == OtherTensor.p))
            {
                throw std::invalid_arguement("Tensor sizes do not match(operator+). Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }
#endif
            std::vector< Matrix<T> > stack_o_matrices(this->p);

            for (size_t k = 0; k < this->p; ++k)
                stack_o_matrices[k] = this->values[k] + OtherTensor.values[k];

            Tensor<T> NewTensor(n, m, p, stack_o_matrices);

            return NewTensor;
        }

        template<class T> Tensor<T> Tensor<T>::operator+ (const T& scalar) const
        {
            Tensor<T> NewTensor(this->n, this->m, this->p, this->values);

            NewTensor += scalar;

            return NewTensor;
        }

        template<class T> Tensor<T>& Tensor<T>::operator+= (const Tensor& OtherTensor)
        {
            //first check to see if the matrices are the same size
#ifndef FAST_EMTG_MATRIX
            if (!(n == OtherTensor.n) || !(m == OtherTensor.m) || !(p == OtherTensor.p))
            {
                throw std::invalid_arguement("Tensor sizes do not match(operator+=). Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }
#endif

            for (size_t k = 0; k < this->p; ++k)
                this->values[k] += OtherTensor.values[k];

            return *this;
        }

        template<class T> Tensor<T>& Tensor<T>::operator+= (const T& scalar)
        {
            //first check to see if the matrices are the same size
#ifndef FAST_EMTG_MATRIX
            if (!(n == OtherTensor.n) || !(m == OtherTensor.m) || !(p == OtherTensor.p))
            {
                throw std::invalid_arguement("Tensor sizes do not match(operator+=). Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));;
            }
#endif

            for (size_t k = 0; k < this->p; ++k)
                this->values[k] += scalar;

            return *this;
        }

        //***************************************************
        //Tensor subtraction
        template<class T> Tensor<T> Tensor<T>::operator- (const Tensor& OtherTensor) const
        {
            //first check to see if the tensors are the same size
#ifndef FAST_EMTG_MATRIX
            if (!(n == OtherTensor.n) || !(m == OtherTensor.m) || !(p == OtherTensor.p))
            {
                throw std::invalid_arguement("Tensor sizes do not match(operator-). Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }
#endif
            std::vector< Matrix<T> > stack_o_matrices(this->p);

            for (size_t k = 0; k < this->p; ++k)
                stack_o_matrices[k] = this->values[k] - OtherTensor.values[k];

            Tensor<T> NewTensor(n, m, p, stack_o_matrices);

            return NewTensor;
        }

        template<class T> Tensor<T> Tensor<T>::operator- (const T& scalar) const
        {
            Tensor<T> NewTensor(this->n, this->m, this->p, this->values);

            NewTensor -= scalar;

            return NewTensor;
        }

        template<class T> Tensor<T>& Tensor<T>::operator-= (const Tensor& OtherTensor)
        {
            //first check to see if the matrices are the same size
#ifndef FAST_EMTG_MATRIX
            if (!(n == OtherTensor.n) || !(m == OtherTensor.m) || !(p == OtherTensor.p))
            {
                throw std::invalid_arguement("Tensor sizes do not match(operator-=). Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }
#endif

            for (size_t k = 0; k < this->p; ++k)
                this->values[k] -= OtherTensor.values[k];

            return *this;
        }

        template<class T> Tensor<T>& Tensor<T>::operator-= (const T& scalar)
        {
            //first check to see if the matrices are the same size
#ifndef FAST_EMTG_MATRIX
            if (!(n == OtherTensor.n) || !(m == OtherTensor.m) || !(p == OtherTensor.p))
            {
                throw std::invalid_arguement("Tensor sizes do not match(operator-=). Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }
#endif

            for (size_t k = 0; k < this->p; ++k)
                this->values[k] -= scalar;

            return *this;
        }



        //multiplication(element-wise)
        template<class T> Tensor<T> Tensor<T>::element_multiply(const Tensor& OtherTensor) const
        {
            //first check to see if the matrices are the same size
#ifndef FAST_EMTG_MATRIX
            if (!(n == OtherTensor.n) || !(m == OtherTensor.m) || !(p == OtherTensor.p))
            {
                throw std::invalid_arguement("Tensor sizes do not match(element multiply). Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }
#endif

            for (size_t k = 0; k < this->p; ++k)
                this->values[k].element_multiply(OtherTensor.values[k]);

            return *this;
        }

        template<class T> Tensor<T>& Tensor<T>::element_multiply_equals(const Tensor& OtherTensor)
        {
            //first check to see if the matrices are the same size
#ifndef FAST_EMTG_MATRIX
            if (!(n == OtherTensor.n) || !(m == OtherTensor.m) || !(p == OtherTensor.p))
            {
                throw std::invalid_arguement("Tensor sizes do not match(element multiply-equals). Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }
#endif

            for (size_t k = 0; k < this->p; ++k)
                this->values[k].element_multiply_equals(OtherTensor.values[k]);

            return *this;
        }

        //division (element-wise)
        template<class T> Tensor<T> Tensor<T>::operator/ (const T& scalar) const
        {
            //first check to see if the matrices are the same size
#ifndef FAST_EMTG_MATRIX
            if (!(n == OtherTensor.n) || !(m == OtherTensor.m) || !(p == OtherTensor.p))
            {
                throw std::invalid_arguement("Tensor sizes do not match(operator/). Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }
#endif


            std::vector< Matrix<T> > stack_o_matrices(this->p);

            for (size_t k = 0; k < this->p; ++k)
                stack_o_matrices[k] = this->values[k] / scalar;

            Tensor<T> NewTensor(n, m, p, stack_o_matrices);

            return NewTensor;
        }

        template<class T> Tensor<T> Tensor<T>::element_divide(const Tensor& OtherTensor) const
        {
            //first check to see if the matrices are the same size
#ifndef FAST_EMTG_MATRIX
            if (!(n == OtherTensor.n) || !(m == OtherTensor.m) || !(p == OtherTensor.p))
            {
                throw std::invalid_arguement("Tensor sizes do not match(element divide). Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }
#endif

            for (size_t k = 0; k < this->p; ++k)
                this->values[k].element_divide(OtherTensor.values[k]);

            return *this;
        }

        template<class T> Tensor<T>& Tensor<T>::operator/= (const T& scalar)
        {

            for (size_t k = 0; k < this->p; ++k)
                this->values[k] /= scalar;

            return *this;
        }

        template<class T> Tensor<T>& Tensor<T>::element_divide_equals(const Tensor& OtherTensor)
        {
            //first check to see if the matrices are the same size
#ifndef FAST_EMTG_MATRIX
            if (!(n == OtherTensor.n) || !(m == OtherTensor.m) || !(p == OtherTensor.p))
            {
                throw std::invalid_arguement("Tensor sizes do not match(element divide-equals). Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }
#endif

            for (size_t k = 0; k < this->p; ++k)
                this->values[k].element_divide_equals(OtherTensor.values[k]);

            return *this;
        }

        //explicit instantiation
#ifdef AD_INSTRUMENTATION
        template class Tensor<GSAD::adouble>;
#endif
        template class Tensor<double>;
    }//close namespace Math
}// close namespace EMTG