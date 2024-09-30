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

#include <algorithm>
#include "boost/algorithm//string.hpp"

#include "HarmonicGravityField.h"
#include "file_utilities.h"

namespace EMTG
{
    namespace Astrodynamics
    {
        //constructors
        HarmonicGravityField::HarmonicGravityField()
        {
            this->spherical_harmonic_term_acceleration_BCF.resize(3, 1, 0.0);
            this->accel_position_Jacobian_BCF.resize(3, 3, 0.0);

            this->degree_order_set = false;
        }

        HarmonicGravityField::HarmonicGravityField(const int & degree_in, const int & order_in)
        {
            this->spherical_harmonic_term_acceleration_BCF.resize(3, 1, 0.0);
            this->accel_position_Jacobian_BCF.resize(3, 3, 0.0);

            this->setDegreeAndOrder(degree_in, order_in);
        }

        HarmonicGravityField::~HarmonicGravityField() {}

        //methods
        void HarmonicGravityField::initializeContainers()
        {
            if (!this->degree_order_set)
            {
                throw std::runtime_error("Error in HarmonicGravityField::initializeContainers. HarmonicGravityField::setDegreeAndOrder must be called. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }

            // initialize everything based on degree because degree n >= order m
            size_t degp1 = this->degree + 1;
            size_t degp4 = this->degree + 4;
            this->J_coeffs.resize(degp4, 0.0);
            this->C_coeffs.resize(degp4, std::vector<double>(degp4, 0.0));
            this->S_coeffs.resize(degp4, std::vector<double>(degp4, 0.0));
            this->ANM.resize(degp4, std::vector<doubleType>(degp4, 0.0));
            this->Nnm.resize(degp4, std::vector<double>(degp4, 0.0));
            this->Beta.resize(degp4, std::vector<double>(degp4, 0.0));
            this->Alpha.resize(degp4, std::vector<double>(degp4, 0.0));
            this->RM.resize(degp4, 0.0);
            this->IM.resize(degp4, 0.0);
            this->NA_n_mp1.resize(degp1, std::vector<double>(degp1, 0.0));
            this->NA_np1_mp1.resize(degp1, std::vector<double>(degp1, 0.0));
            this->NA_n_mp2.resize(degp1, std::vector<double>(degp1, 0.0));
            this->NA_np1_mp2.resize(degp1, std::vector<double>(degp1, 0.0));
            this->NA_np2_mp2.resize(degp1, std::vector<double>(degp1, 0.0));
            this->initializeNormalization();
        }

        void HarmonicGravityField::initializeNormalization()
        {
            // diagonal recursion
            // diagonal derived Legendre Polynomials do not depend on the spacecraft state
            // these can be initialized at construction time
            this->ANM[0][0] = 1.0;
            for (size_t n = 1; n < this->degree + 3; ++n)
            {
                this->ANM[n][n] = sqrt((2.0*n + 1.0) / (2.0*n)) * this->ANM[n - 1][n - 1];
            }

            // compute mass coefficient normalization constants
            //   N_n,0 = sqrt(2n+1)
            //   N_n,m = sqrt(2(2n+1) * (n-m)! / (n+m)!)
            //   
            //   N_n,m = N_n,m-1 / sqrt( (n+m) * (n-m+1) )
            for (size_t n = 0; n < this->degree + 3; ++n)
            {
                // Temporary, to make following loop work
                this->Nnm[n][0] = sqrt((2.0 * (2.0 * n + 1.0)));
                for (size_t m = 1; m <= n + 2 && m < this->order + 3; ++m)
                {
                    this->Nnm[n][m] = this->Nnm[n][m - 1] / sqrt(double(n + m) * (n - m + 1.0));
                }
                // Now set true value
                this->Nnm[n][0] = sqrt((2.0 * n + 1.0));
            }

            // Pre-initialize the derived Legendre Polynomial normalization coefficients
            // Here we are using the normalization scheme from: 
            // Lundberg, J. B., & Schutz, B. E. (1988). Recursion Formulas of Legendre Functions for Use with Nonsingular Geopotential Models. 
            // Journal of Guidance, Control, and Dynamics, 11(1), 31–38. https://doi.org/10.2514/3.20266
            for (size_t m = 0; m < this->order + 3; ++m)
            {
                for (size_t n = m + 2; n < this->degree + 3; ++n)
                {
                    this->Alpha[n][m] = sqrt(((2.0 * n + 1.0) * (2.0 * n - 1.0)) / ((n - m) * (n + m)));
                    this->Beta[n][m] = sqrt(((2.0 * n + 1.0) * (n - m - 1.0) * (n + m - 1.0)) / ((2.0 * n - 3.0) * (n + m) * (n - m)));
                }
            }

            // Normalization constants for further derivatives of the derived Legendre Polynomials
            double sqrt2 = sqrt(2.0);
            for (size_t n = 0; n <= this->degree; ++n)
            {
                for (size_t m = 0; m <= n && m <= this->order; ++m)
                {
                    double ndub = n;
                    NA_n_mp1[n][m]   = sqrt((n - m) * (n + m + 1.0));
                    NA_np1_mp1[n][m] = sqrt(((2.0 * n + 1.0) * (n + m + 2.0) * (n + m + 1.0)) / (2.0 * n + 3.0));
                    NA_n_mp2[n][m]   = sqrt((n - m) * (n - m - 1.0) * (n + m + 1.0) * (n + m + 2.0));
                    NA_np1_mp2[n][m] = sqrt((2.0 * n + 1.0) / (2.0 * n + 3.0) * ((n - m) * (n + m + 1.0) * (n + m + 2.0) * (n + m + 3.0)));
                    NA_np2_mp2[n][m] = sqrt((2.0 * n + 1.0) / (2.0 * n + 5.0) * ((n + m + 1.0) * (n + m + 2.0) * (n + m + 3.0) * (n + m + 4.0)));
                    if (m == 0)
                    {
                        NA_n_mp1[n][m]   /= sqrt2;
                        NA_np1_mp1[n][m] /= sqrt2;
                        NA_n_mp2[n][m]   /= sqrt2;
                        NA_np1_mp2[n][m] /= sqrt2;
                        NA_np2_mp2[n][m] /= sqrt2;
                    }
                }
            }
        }
        
        void HarmonicGravityField::computeDerivedLegendrePolynomials()
        {
            // off-diagonal recursion
            this->ANM[1][0] = this->U * sqrt(3.0);
            for (size_t n = 1; n < this->degree + 3; ++n)
            {
                this->ANM[n + 1][n] = this->U * sqrt(2.0 * n + 3.0) * this->ANM[n][n];
            }

            // column-fill recursion
            for (size_t m = 0; m < this->order + 3; ++m)
            {
                // Lundberg and Schutz Table 2, Row 1
                for (size_t n = m + 2; n < this->degree + 3; ++n)
                {
                    this->ANM[n][m] = this->U * this->Alpha[n][m] * this->ANM[n - 1][m] - this->Beta[n][m] * this->ANM[n - 2][m];
                }
                // Pines Eq.(24)
                this->RM[m] = m == 0 ? 1.0 : this->S * this->RM[m - 1] - this->T * this->IM[m - 1]; // real part of (this->s + i*this->T)^m
                this->IM[m] = m == 0 ? 0.0 : this->S * this->IM[m - 1] + this->T * this->RM[m - 1]; // imaginary part of (this->s + i*this->T)^m
            }
        }

        void HarmonicGravityField::computeHarmonicGravityField(const math::Matrix<doubleType> & r_body2sc_BCF,
                                                               const bool & generate_derivatives)
        {
            this->r_body2sc_BCF_norm = r_body2sc_BCF.norm();
            this->S = r_body2sc_BCF(0) / this->r_body2sc_BCF_norm;
            this->T = r_body2sc_BCF(1) / this->r_body2sc_BCF_norm;
            this->U = r_body2sc_BCF(2) / this->r_body2sc_BCF_norm;

            this->computeDerivedLegendrePolynomials();

            // Initialize recursion in Pines Eq. 26
            doubleType rho = this->harmonic_field_reference_radius / this->r_body2sc_BCF_norm;
            doubleType rho0 = this->mu / this->r_body2sc_BCF_norm;
            doubleType rho_np1 = rho * rho0;
            doubleType rho_np2 = rho_np1 * rho;

            // Initialize Pines Eq. 30
            doubleType a1 = 0.0;
            doubleType a2 = 0.0;
            doubleType a3 = 0.0;
            doubleType a4 = 0.0;
            doubleType a11 = 0.0;
            doubleType a12 = 0.0;
            doubleType a13 = 0.0;
            doubleType a14 = 0.0;
            doubleType a23 = 0.0;
            doubleType a24 = 0.0;
            doubleType a33 = 0.0;
            doubleType a34 = 0.0;
            doubleType a44 = 0.0;
            double sqrt2 = sqrt(2.0);
            for (size_t n = 1; n <= this->degree; ++n)
            {
                rho_np1 *= rho;
                doubleType a1_sum = 0.0;
                doubleType a2_sum = 0.0;
                doubleType a3_sum = 0.0;
                doubleType a4_sum = 0.0;
                doubleType a11_sum = 0.0;
                doubleType a12_sum = 0.0;
                doubleType a13_sum = 0.0;
                doubleType a14_sum = 0.0;
                doubleType a23_sum = 0.0;
                doubleType a24_sum = 0.0;
                doubleType a33_sum = 0.0;
                doubleType a34_sum = 0.0;
                doubleType a44_sum = 0.0;

                for (size_t m = 0; m <= n && m <= this->order; ++m)
                {
                    double Cnm = this->C_coeffs[n][m];
                    double Snm = this->S_coeffs[n][m];
                    // Pines Eq. 27
                    doubleType D = (Cnm * this->RM[m] + Snm * this->IM[m]) * sqrt2;
                    doubleType E = m == 0 ? 0 : (Cnm * this->RM[m - 1] + Snm * this->IM[m - 1]) * sqrt2;
                    doubleType F = m == 0 ? 0 : (Snm * this->RM[m - 1] - Cnm * this->IM[m - 1]) * sqrt2;
                    // Correct for normalization
                    doubleType Avv00 =                          this->ANM[n][m];
                    doubleType Avv01 = this->NA_n_mp1[n][m]   * this->ANM[n][m + 1];
                    doubleType Avv11 = this->NA_np1_mp1[n][m] * this->ANM[n + 1][m + 1];
                    // Pines Eq. 30
                    a1_sum += m * Avv00 * E;
                    a2_sum += m * Avv00 * F;
                    a3_sum +=     Avv01 * D;
                    a4_sum +=     Avv11 * D; // Pines Eq. 30b

                    if (generate_derivatives)
                    {
                        // Pines Eq. 27
                        doubleType G = m <= 1 ? 0 : (Cnm * this->RM[m - 2] + Snm * this->IM[m - 2]) * sqrt2;
                        doubleType H = m <= 1 ? 0 : (Snm * this->RM[m - 2] - Cnm * this->IM[m - 2]) * sqrt2;
                        
                        // Correct for normalization
                        doubleType Avv02 = this->NA_n_mp2[n][m]   * this->ANM[n][m + 2];
                        doubleType Avv12 = this->NA_np1_mp2[n][m] * this->ANM[n + 1][m + 2];
                        doubleType Avv22 = this->NA_np2_mp2[n][m] * this->ANM[n + 2][m + 2];

                        // Pines Eq. 36
                        a11_sum += m * (m - 1) * Avv00 * G;
                        a12_sum += m * (m - 1) * Avv00 * H;
                        a13_sum += m * Avv01 * E;
                        a14_sum += m * Avv11 * E;
                        a23_sum += m * Avv01 * F;
                        a24_sum += m * Avv11 * F;
                        a33_sum += Avv02 * D;
                        a34_sum += Avv12 * D;
                        a44_sum += Avv22 * D;
                    }

                } // end loop over gravity field orders

                // Pines Eq. 30 inner sum constant multiply
                doubleType rho_np1_over_ref_rad = rho_np1 / this->harmonic_field_reference_radius;
                a1 += rho_np1_over_ref_rad * a1_sum;
                a2 += rho_np1_over_ref_rad * a2_sum;
                a3 += rho_np1_over_ref_rad * a3_sum;
                a4 -= rho_np1_over_ref_rad * a4_sum; // Pines Eq. 30b

                // Pines Eq. 31 
                this->spherical_harmonic_term_acceleration_BCF(0) = a1 + a4 * this->S;
                this->spherical_harmonic_term_acceleration_BCF(1) = a2 + a4 * this->T;
                this->spherical_harmonic_term_acceleration_BCF(2) = a3 + a4 * this->U;

                if (generate_derivatives)
                {
                    rho_np2 *= rho;

                    // Pines Eq. 36
                    doubleType rho_np2_over_ref_rad2 = rho_np2 / (this->harmonic_field_reference_radius * this->harmonic_field_reference_radius);
                    a11 += rho_np2_over_ref_rad2 * a11_sum;
                    a12 += rho_np2_over_ref_rad2 * a12_sum;
                    a13 += rho_np2_over_ref_rad2 * a13_sum;
                    a14 -= rho_np2_over_ref_rad2 * a14_sum;
                    a23 += rho_np2_over_ref_rad2 * a23_sum;
                    a24 -= rho_np2_over_ref_rad2 * a24_sum;
                    a33 += rho_np2_over_ref_rad2 * a33_sum;
                    a34 -= rho_np2_over_ref_rad2 * a34_sum;
                    a44 += rho_np2_over_ref_rad2 * a44_sum;

                    // Pines Eq. 37
                    // symmetric matrix, so cross-diagonal entries are equivalent
                    this->accel_position_Jacobian_BCF(0, 0) =  a11 + this->S * this->S * a44 + a4 / this->r_body2sc_BCF_norm + 2.0 * this->S * a14;
                    this->accel_position_Jacobian_BCF(1, 1) = -a11 + this->T * this->T * a44 + a4 / this->r_body2sc_BCF_norm + 2.0 * this->T * a24;
                    this->accel_position_Jacobian_BCF(2, 2) =  a33 + this->U * this->U * a44 + a4 / this->r_body2sc_BCF_norm + 2.0 * this->U * a34;
                    this->accel_position_Jacobian_BCF(0, 1) =
                    this->accel_position_Jacobian_BCF(1, 0) =  a12 + this->S * this->T * a44 + this->S * a24 + this->T * a14;
                    this->accel_position_Jacobian_BCF(0, 2) =
                    this->accel_position_Jacobian_BCF(2, 0) =  a13 + this->S * this->U * a44 + this->S * a34 + this->U * a14;
                    this->accel_position_Jacobian_BCF(1, 2) =
                    this->accel_position_Jacobian_BCF(2, 1) =  a23 + this->T * this->U * a44 + this->U * a24 + this->T * a34;
                }

            } // end loop over gravity field degrees

        }

        void HarmonicGravityField::parseSTKgrvFile(const std::string & gravity_file_string_in)
        {
            std::ifstream grav_file;
            grav_file.open(gravity_file_string_in, std::ios::in);

            if (!grav_file)
            {
                throw std::runtime_error("Gravity file: " + gravity_file_string + " was not opened successfully. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }

            this->gravity_file_string = gravity_file_string_in;

            std::string line;

            std::string Body;

            bool is_normalized;
            bool reading_mass_coefficients = false;
            bool reading_love_numbers = false;
            while (EMTG::file_utilities::safeGetline(grav_file, line))
            {
                if (line.size() > 0) // skip blank lines
                {
                    std::vector<std::string> linecell;
                    std::replace(line.begin(), line.end(),  '\t', ' ');
                    boost::split(linecell, line, boost::is_any_of(" "), boost::token_compress_on);
                    linecell.erase(
                        std::remove_if(
                            linecell.begin(),
                            linecell.end(),
                            [](std::string const & s) { return s.empty(); }),
                        linecell.end());

                    if (linecell.front() == "#")
                    {
                        continue;
                    }
                    else if (linecell.front() == "CentralBody")
                    {
                        Body = linecell[1];
                    }
                    else if (linecell.front() == "Degree")
                    {
                        this->max_degree = std::stoi(linecell[1]);
                    }
                    else if (linecell.front() == "Order")
                    {
                        this->max_order = std::stoi(linecell[1]);
                    }
                    else if (linecell.front() == "Gm")
                    {
                        this->mu = std::stod(linecell[1]);
                        this->mu /= 1.0E+09; // m^3/s^2 to km^3/s^2 conversion
                    }
                    else if (linecell.front() == "RefDistance")
                    {
                        this->harmonic_field_reference_radius = std::stod(linecell[1]);
                        this->harmonic_field_reference_radius /= 1000.0; // m to km conversion
                    }
                    else if (linecell.front() == "Normalized")
                    {
                        bool a_decision_was_made = false;
                        for (size_t k = 0; k < linecell.size(); ++k)
                        {
                            if (linecell[k] == "Yes")
                            {
                                a_decision_was_made = true;
                                is_normalized = true;
                                break;
                            }
                            else if (linecell[k] == "No")
                            {
                                a_decision_was_made = true;
                                is_normalized = false;
                                break;
                            }
                        }
                        if (!a_decision_was_made)
                        {
                            throw std::invalid_argument("Gravity file " + gravity_file_string + " must specify 'Normalized Yes/No'");
                        }
                        continue;
                    }
                    else if (linecell.front() == "BEGIN" && linecell.back() == "Coefficients")
                    {
                        reading_mass_coefficients = true;
                        continue;
                    }
                    else if (linecell.front() == "END" && linecell.back() == "Coefficients")
                    {
                        reading_mass_coefficients = false;
                        continue;
                    }
                    else if (linecell.front() == "Begin" && linecell.back() == "LoveNumbers")
                    {
                        reading_love_numbers = true;
                        continue;
                    }
                    else if (linecell.front() == "End" && linecell.back() == "LoveNumbers")
                    {
                        reading_love_numbers = false;
                        continue;
                    }

                    // check to see if we're currently in the process of reading something
                    if (reading_mass_coefficients)
                    {
                        this->readMassCoefficientLine(linecell, is_normalized);
                    }
                    else if (reading_love_numbers)
                    {
                        this->readLoveNumbers(linecell);
                    }
                }
            } // end loop through gravity file lines

        } // end parse_STK_grv_file

        void HarmonicGravityField::readMassCoefficientLine(const std::vector<std::string> & linecell,
                                                           const bool & is_normalized)
        {
            int n = std::stoi(linecell[0]);
            int m = std::stoi(linecell[1]);
            double Ccoeff = std::stod(linecell[2]);
            double Scoeff = std::stod(linecell[3]);

            // if we have gone beyond what the user wants to model, then don't try loading these
            // coefficients into memory
            if (n > this->degree || m > this->order)
            {
                return;
            }

            // if the coefficients are not normalized...we must normalize them!
            if (!is_normalized)
            {
                Ccoeff /= this->Nnm[n][m];
                Scoeff /= this->Nnm[n][m];
            }

            if (m == 0)
            {
                // by convention, J[n] = -C[n][0]
                this->J_coeffs[n] = -Ccoeff;
                this->C_coeffs[n][0] = Ccoeff;
            }
            else
            {
                this->C_coeffs[n][m] = Ccoeff;
                this->S_coeffs[n][m] = Scoeff;
            }
        }

        void HarmonicGravityField::readLoveNumbers(const std::vector<std::string> & linecell)
        {
            for (size_t k = 0; k < linecell.size(); ++k)
            {
                this->love_numbers.push_back(std::stod(linecell[k]));
            }
        }
        
    }//close namespace Astrodynamics
}//close namespace EMTG