// EMTG: Evolutionary Mission Trajectory Generator
// An open-source global optimization tool for preliminary mission design
// Provided by NASA Goddard Space Flight Center
//
// Copyright (c) 2013 - 2024 United States Government as represented by the
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

// harmonic gravity field model

#ifndef HARMONIC_GRAVITY_FIELD_H
#define HARMONIC_GRAVITY_FIELD_H

#include <vector>
#include "doubleType.h"
#include "EMTG_Matrix.h"

namespace EMTG
{
    namespace Astrodynamics
    {
        class HarmonicGravityField
        {
        public:
            // constructors
            HarmonicGravityField();
            HarmonicGravityField(const int & degree_in, const int & order_in);

            virtual ~HarmonicGravityField();

            // methods
            inline void setCentralBodyGM(const double & mu) { this->mu = mu; };
            inline void setReferenceRadius(const double & ref_radius) { this->harmonic_field_reference_radius = ref_radius; };
            inline void setDegreeAndOrder(const size_t & degree_in, const size_t & order_in) 
            {
                if (degree_in < order_in)
                {
                    throw std::invalid_argument("Error in HarmonicGravityField::setDegreeAndOrder. Degree must be >= Order. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
                }

                this->degree = degree_in;
                this->order = order_in;
                this->degree_order_set = true;
                this->initializeContainers();
            }
            inline void setMassCoefficients(const std::vector< std::vector<double>> & Cnm, const std::vector< std::vector<double>> & Snm, const bool & are_normalized) 
            { 
                if (!this->degree_order_set)
                {
                    throw std::runtime_error("Error in HarmonicGravityField::setMassCoefficients. HarmonicGravityField::setDegreeAndOrder must be called. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
                }
                this->C_coeffs = Cnm; 
                this->S_coeffs = Snm; 
                if (!are_normalized)
                {
                    for (size_t n = 0; n < this->degree + 3; ++n)
                    {
                        for (size_t m = 0; m < this->order + 3; ++m)
                        {
                            this->C_coeffs[n][m] /= this->Nnm[n][m];
                            this->S_coeffs[n][m] /= this->Nnm[n][m];
                        }
                    }
                }
            };

            math::Matrix<doubleType> getBCFaccel() const { return this->spherical_harmonic_term_acceleration_BCF; };
            math::Matrix<doubleType> getBCFpositionJacobian() const { return this->accel_position_Jacobian_BCF; };
            double getCentralBodyGM() const { return this->mu; };
            double getReferenceRadius() const { return this->harmonic_field_reference_radius; };

            void computeHarmonicGravityField(const math::Matrix<doubleType> & r_body2sc_BCF, 
                                             const bool & generate_derivatives);

            void parseSTKgrvFile(const std::string & gravity_file_string_in);

        protected:

            void readMassCoefficientLine(const std::vector<std::string> & linecell,
                                            const bool& is_normalized);

            void readLoveNumbers(const std::vector<std::string> & linecell);

            void initializeContainers();
            void initializeNormalization();
            void computeDerivedLegendrePolynomials();

            // fields
            size_t degree;
            size_t order;
            bool degree_order_set;
            size_t max_degree;
            size_t max_order;

            double mu;
            double harmonic_field_reference_radius;

            std::string gravity_file_string;
            std::vector<double> love_numbers;
            std::vector<double> J_coeffs;
            std::vector< std::vector<double> > C_coeffs;
            std::vector< std::vector<double> > S_coeffs;
            std::vector< std::vector<double> > Nnm;
            std::vector< std::vector<double> > Beta;
            std::vector< std::vector<double> > Alpha;
            std::vector< std::vector<double> > NA_n_mp1;
            std::vector< std::vector<double> > NA_np1_mp1;
            std::vector< std::vector<double> > NA_n_mp2;
            std::vector< std::vector<double> > NA_np1_mp2;
            std::vector< std::vector<double> > NA_np2_mp2;
            
            doubleType S;
            doubleType T;
            doubleType U;
            std::vector<doubleType> RM;
            std::vector<doubleType> IM;
            std::vector< std::vector<doubleType> > ANM;
            math::Matrix<doubleType> r_body2sc_BCF;
            doubleType r_body2sc_BCF_norm;
            math::Matrix<doubleType> spherical_harmonic_term_acceleration_BCF;
            math::Matrix<doubleType> accel_position_Jacobian_BCF;

        };
    } // end Astrodynamics namespace
} // end EMTG namespace
#endif // HARMONIC_GRAVITY_FIELD_H