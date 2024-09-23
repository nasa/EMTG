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

//EMTG frame class
//Jacob Englander 1/3/2013
#include "frame.h"

#include "EMTG_Matrix.h"
#include "EMTG_Tensor.h"

namespace EMTG 
{ 
    namespace Astrodynamics
    {

        //default constructor initializes to default Earth and then constructs dummy rotation matrices for J2000
        frame::frame(void)
        {
            this->initialize();
        }

        //constructor with data
        frame::frame(const double& alpha0_in, const double& alphadot_in, const double& delta0_in, const double& deltadot_in, const double& W_in, const double& Wdot_in)
        {
            this->initialize(alpha0_in, alphadot_in, delta0_in, deltadot_in, W_in, Wdot_in);
        }

        frame::frame(const double& alpha0_in,
            const double& alphadot_in,
            const double& delta0_in,
            const double& deltadot_in,
            const double& W_in,
            const double& Wdot_in,
            const double& PA1_0,
            const double& PA1dot,
            const double& PA2_0,
            const double& PA2dot,
            const double& PA3_0,
            const double& PA3dot)
        {
            this->initialize(alpha0_in, 
                alphadot_in,
                delta0_in, 
                deltadot_in,
                W_in,
                Wdot_in,
                PA1_0,
                PA1dot,
                PA2_0,
                PA2dot,
                PA3_0,
                PA3dot);
        }

        //destructor
        frame::~frame(void)
        {
        }

        //*************************************************************
        //methods

        //initialization method
        void frame::initialize(const double& alpha0,
            const double& alphadot,
            const double& delta0,
            const double& deltadot,
            const double& W,
            const double& Wdot,
            const double& PA1_0,
            const double& PA1dot,
            const double& PA2_0,
            const double& PA2dot,
            const double& PA3_0,
            const double& PA3dot)
        {
            this->alpha0 = alpha0;
            this->alphadot = alphadot;
            this->delta0 = delta0;
            this->deltadot = deltadot;
            this->W0 = W;
            this->Wdot = Wdot;
            this->PA1_0 = PA1_0;
            this->PA1dot = PA1dot;
            this->PA2_0 = PA2_0;
            this->PA2dot = PA2dot;
            this->PA3_0 = PA3_0;
            this->PA3dot = PA3dot;

            xhat = math::Matrix<doubleType>(3, 1, 0.0);
            yhat = xhat;
            zhat = xhat;
            xhat(0) = 1.0;
            yhat(1) = 1.0;
            zhat(2) = 1.0;
            this->semi_axis_a = 1.0;
            this->semi_axis_b = 1.0;
            this->semi_axis_c = 1.0;

            this->construct_rotation_matrices_J2000();

            this->dR_dreferenceVector = math::Tensor<doubleType>(3, 3, 3);
        }

        //initialization method for Earth J2000 frame, i.e. for when rotation matrices should be the identity matrix
        void frame::initialize()
        {
            this->alpha0 = -math::PIover2;
            this->alphadot = 0.0;
            this->delta0 = math::PIover2;
            this->deltadot = 0.0;
            this->W0 = 3.3187;
            this->Wdot = 6.3004;

            xhat = math::Matrix<doubleType>(3, 1, 0.0);
            yhat = xhat;
            zhat = xhat;
            xhat(0) = 1.0;
            yhat(1) = 1.0;
            zhat(2) = 1.0;

            this->construct_rotation_matrices_J2000();

            this->dR_dreferenceVector = math::Tensor<doubleType>(3, 3, 3);
        }//end initialize()

        //method to initialize with ecliptic reference angles
        void frame::initialize_ecliptic()
        {
            this->initialize(-90.0 * math::deg2rad,
                0.0 * math::deg2rad,
                66.560709000000003 * math::deg2rad,
                0.0 * math::deg2rad,
                84.176 * math::deg2rad,
                14.1844000 * math::deg2rad);
        }//end initialize_ecliptic()

        //construct the rotation matrices between ICRF and the local frame
        //epoch is in ET MJD seconds (i.e., seconds past the J2000 epoch, i.e., the time state of the state vector)
        void frame::construct_rotation_matrices(const doubleType& ETepoch, const bool& GenerateDerivatives)
        {
            doubleType days_since_reference_epoch = ETepoch / 86400.0 - 51544.5;
            doubleType centuries_since_reference_epoch = days_since_reference_epoch / 36525.0;

            //compute the current values of the angles alpha and delta
            this->alpha = this->alpha0 + this->alphadot * centuries_since_reference_epoch;
            this->delta = this->delta0 + this->deltadot * centuries_since_reference_epoch;
            this->W = this->W0 + this->Wdot * days_since_reference_epoch;

            //next, construct the rotation matrix and the derivative of the spin rotation
            math::Matrix<doubleType> Rx(3, (math::PIover2 - delta), math::Rxhat);
            math::Matrix<doubleType> Rz(3, (math::PIover2 + alpha), math::Rzhat);
            this->R_from_BCI_to_ICRF = Rz*Rx;
            this->R_from_ICRF_to_BCI = this->R_from_BCI_to_ICRF.transpose();
            this->R_from_BCF_to_BCI = math::Matrix<doubleType>(3, this->W, math::Rzhat);
            this->R_from_BCI_to_BCF = this->R_from_BCF_to_BCI.transpose();
                       
            //J2000 spin matrix
            this->R_from_J2000BCF_to_J2000BCI = math::Matrix<doubleType>(3, this->W, math::MatrixType::Rzhat);
            this->R_from_J2000BCI_to_J2000BCF = this->R_from_J2000BCF_to_J2000BCI.transpose();

            //and the full transformation
            this->R_from_BCF_to_ICRF = this->R_from_BCI_to_ICRF * this->R_from_BCF_to_BCI;
            this->R_from_ICRF_to_BCF = this->R_from_BCF_to_ICRF.transpose();

            //J2000 version
            this->R_from_J2000BCF_to_ICRF = this->R_from_J2000BCI_to_ICRF * this->R_from_J2000BCF_to_J2000BCI;
            this->R_from_ICRF_to_J2000BCF = this->R_from_J2000BCF_to_ICRF.transpose();

            if (GenerateDerivatives)
            {
                //derivative of Rx and Rx transpose, remember that deltadot is in radians/century
                math::Matrix<doubleType> dRx_dt(3, (math::PIover2 - this->delta), -this->deltadot / 36525.0 / 86400.0, math::Rxhatdot);

                //derivative of Rz and Rz transpose, remember that alphadot is in radians/century
                math::Matrix<doubleType> dRz_dt(3, (math::PIover2 + this->alpha), this->alphadot / 36525.0 / 86400.0, math::Rzhatdot);

                //chained derivative from ICRF to BCI
                this->dR_from_BCI_to_ICRF_dt = Rz * dRx_dt + dRz_dt * Rx;

                //chained derivative from BCI to ICRF
                this->dR_from_ICRF_to_BCI_dt = this->dR_from_BCI_to_ICRF_dt.transpose();

                //derivative of spin matrix
                this->dR_from_BCF_to_BCI_dt = math::Matrix<doubleType>(3, this->W, this->Wdot / 86400.0, math::Rzhatdot);
                this->dR_from_BCI_to_BCF_dt = this->dR_from_BCF_to_BCI_dt.transpose();

                //total derivative
                this->dR_from_BCF_to_ICRF_dt = this->R_from_BCI_to_ICRF * this->dR_from_BCF_to_BCI_dt + this->dR_from_BCI_to_ICRF_dt * this->R_from_BCF_to_BCI;
                this->dR_from_ICRF_to_BCF_dt = this->dR_from_BCF_to_ICRF_dt.transpose();


                //J2000 version
                //derivative of spin matrix
                this->dR_from_J2000BCF_to_J2000BCI_dt = math::Matrix<doubleType>(3, this->W, this->Wdot / 86400.0, math::Rzhatdot);
                this->dR_from_J2000BCI_to_J2000BCF_dt = this->dR_from_J2000BCF_to_J2000BCI_dt.transpose();

                //total derivative
                this->dR_from_J2000BCF_to_ICRF_dt = this->R_from_J2000BCI_to_ICRF * this->dR_from_J2000BCF_to_J2000BCI_dt;//this next term is always zero + this->dR_from_J2000BCI_to_ICRF_dt * this->R_from_J2000BCF_to_J2000BCI;
                this->dR_from_ICRF_to_J2000BCF_dt = this->dR_from_J2000BCF_to_ICRF_dt.transpose();
            }
        }//end construct_rotation_matrices


        void frame::construct_PrincipleAxes_to_BCF_rotation(const doubleType& ETepoch, 
            const bool& GenerateDerivatives)
        {
            doubleType days_since_reference_epoch = ETepoch / 86400.0 - 51544.5;

            //PrincipleAxes
            this->PA1 = this->PA1_0 + this->PA1dot * days_since_reference_epoch;
            this->PA2 = this->PA2_0 + this->PA2dot * days_since_reference_epoch;
            this->PA3 = this->PA3_0 + this->PA3dot * days_since_reference_epoch;
            math::Matrix<doubleType> R_from_BCF_to_Fprime = math::Matrix<doubleType>(3, this->PA3, math::Rzhat);
            math::Matrix<doubleType> R_from_Fprime_to_Fdprime = math::Matrix<doubleType>(3, this->PA2, math::Rxhat);
            math::Matrix<doubleType> R_from_Fdprime_to_PrincipleAxes = math::Matrix<doubleType>(3, this->PA1, math::Rzhat);
            this->R_from_BCF_to_PrincipleAxes = R_from_Fdprime_to_PrincipleAxes * R_from_Fprime_to_Fdprime * R_from_BCF_to_Fprime;
            this->R_from_PrincipleAxes_to_BCF = this->R_from_BCF_to_PrincipleAxes.transpose();

            if (GenerateDerivatives)
            {
                math::Matrix<doubleType> dBCF_to_Fprime_dt(3, this->PA3, this->PA3dot / 86400.0, math::Rzhatdot);
                math::Matrix<doubleType> dFprime_to_Fdprime_dt(3, this->PA2, this->PA2dot / 86400.0, math::Rxhatdot);
                math::Matrix<doubleType> dFdprime_to_PrincipleAxes_dt(3, this->PA1, this->PA1dot / 86400.0, math::Rzhatdot);
                math::Matrix<doubleType> d_R_d_PA1 = R_from_Fdprime_to_PrincipleAxes * R_from_Fprime_to_Fdprime * dBCF_to_Fprime_dt;
                math::Matrix<doubleType> d_R_d_PA2 = R_from_Fdprime_to_PrincipleAxes * dFprime_to_Fdprime_dt * R_from_BCF_to_Fprime;
                math::Matrix<doubleType> d_R_d_PA3 = dFdprime_to_PrincipleAxes_dt * R_from_Fprime_to_Fdprime * R_from_BCF_to_Fprime;
                this->dR_from_BCF_to_PrincipleAxes_dt = d_R_d_PA1 + d_R_d_PA2 + d_R_d_PA3;
                this->dR_from_PrincipleAxes_to_BCF_dt = this->dR_from_BCF_to_PrincipleAxes_dt.transpose();
            }
        }//end construct_Topocentric_to_BCF_rotation()

        void frame::construct_Topocentric_to_BCF_rotation(const doubleType& ETepoch, 
            math::Matrix<doubleType> rBCF,
            const bool& GenerateDerivatives)
        {
            doubleType days_since_reference_epoch = ETepoch / 86400.0 - 51544.5;

            //Step 1: calculate unit vectors of Topocentric frame in BCF frame
            //Step 1.1: calculate unit vectors of outward normal vector in BCF
            //Step 1.1.1: calculate the outward normal vector in PrincipleAxes
            this->construct_PrincipleAxes_to_BCF_rotation(ETepoch, GenerateDerivatives);
            math::Matrix<doubleType> rPrincipleAxes = this->R_from_BCF_to_PrincipleAxes * rBCF;
            math::Matrix<doubleType> ellipseMatrix = math::Matrix<doubleType>(3, std::vector<doubleType>({ this->semi_axis_a * this->semi_axis_a,
                                                                                                           this->semi_axis_b * this->semi_axis_b,
                                                                                                           this->semi_axis_c * this->semi_axis_a }),
                                                                                                           math::diagonal);

            math::Matrix<doubleType> nPrincipleAxes = ellipseMatrix * rPrincipleAxes * 2.0;

            //Step 1.1.2: rotate the outward normal "up" vector to BCF, and create a unit vector
            math::Matrix<doubleType> uBCF = this->R_from_PrincipleAxes_to_BCF * nPrincipleAxes;
            math::Matrix<doubleType> unit_uBCF = uBCF.unitize();

            //Step 1.2: e (east) vector in BCF frame
            math::Matrix<doubleType> eBCF = this->zhat.cross(this->R_from_PrincipleAxes_to_BCF * ellipseMatrix * this->R_from_BCF_to_PrincipleAxes * rBCF * 2.0);
            math::Matrix<doubleType> unit_eBCF = eBCF.unitize();

            //Step 1.3: pseudo-s (south) vector in BCF frame
            math::Matrix<doubleType> sBCF = eBCF.cross(uBCF);
            math::Matrix<doubleType> unit_sBCF = sBCF.unitize();

            //Step 2: form the transformation matrix
            this->R_from_BCF_to_Topocentric = unit_sBCF.transpose().vert_cat(unit_eBCF.transpose().vert_cat(unit_uBCF.transpose()));
            this->R_from_Topocentric_to_BCF = this->R_from_BCF_to_Topocentric.transpose();

            if (GenerateDerivatives)
            {
                //Step 3.1: generate the derivatives of the unit vectors
                //Step 3.1.1: uBCF
                math::Matrix<doubleType> d_nPrincipleAxes_drPrincipleAxes = ellipseMatrix * 2.0;
                math::Matrix<doubleType> d_uBCF_drPrincipleAxes = this->R_from_PrincipleAxes_to_BCF * d_nPrincipleAxes_drPrincipleAxes;
                math::Matrix<doubleType> d_unit_uBCF_duBCF = unit_uBCF.unitDerivative(uBCF);

                //Step 3.1.1.1: uBCF derivative with respect to time
                math::Matrix<doubleType> d_uBCF_dt = this->dR_from_PrincipleAxes_to_BCF_dt * nPrincipleAxes +
                                                     this->R_from_PrincipleAxes_to_BCF * ellipseMatrix * this->dR_from_BCF_to_PrincipleAxes_dt * rBCF * 2.0;
                math::Matrix<doubleType> d_unit_uBCF_dt = d_unit_uBCF_duBCF * d_uBCF_dt;
                    
                //Step 3.1.1.2: uBCF derivative with respect to rBCF
                math::Matrix<doubleType> d_uBCF_d_rBCF = d_uBCF_drPrincipleAxes * this->R_from_BCF_to_PrincipleAxes;
                math::Matrix<doubleType> d_unit_uBCF_d_rBCF = d_unit_uBCF_duBCF * d_uBCF_d_rBCF;

                //Step 3.1.1.3: sBCF derivative with respect to vBCF
                math::Matrix<doubleType> d_uBCF_dvBCF = math::Matrix<doubleType>(3, 3, 0.0);

                //Step 3.1.2: eBCF
                math::Matrix<doubleType> d_unit_eBCF_deBCF = unit_eBCF.unitDerivative(eBCF);

                //Step 3.1.2.1: eBCF derivative with respect to time
                math::Matrix<doubleType> d_eBCF_dt = zhat.cross((this->dR_from_PrincipleAxes_to_BCF_dt * ellipseMatrix * this->R_from_BCF_to_PrincipleAxes + this->R_from_PrincipleAxes_to_BCF * ellipseMatrix * this->dR_from_BCF_to_PrincipleAxes_dt) * rBCF) * 2.0;
                math::Matrix<doubleType> d_unit_eBCF_dt = d_unit_eBCF_deBCF * d_eBCF_dt;

                //Step 3.1.2.2: eBCF derivative with respect to rBCF
                math::Matrix<doubleType> I3 = math::Matrix<doubleType>(3, math::MatrixType::identity);
                math::Matrix<doubleType> kcross = math::Matrix<doubleType>(zhat, math::MatrixType::skewsymmetric);
                math::Matrix<doubleType> eBCFcross = math::Matrix<doubleType>(eBCF, math::MatrixType::skewsymmetric);
                math::Matrix<doubleType> d_eBCF_drBCF = kcross * this->R_from_PrincipleAxes_to_BCF * ellipseMatrix * this->R_from_BCF_to_PrincipleAxes * I3 * 2.0;
                math::Matrix<doubleType> d_unit_eBCF_d_rBCF = d_unit_eBCF_deBCF * d_eBCF_drBCF;

                //Step 3.1.2.3: eBCF derivative with respect to vBCF
                math::Matrix<doubleType> d_eBCF_dvBCF = math::Matrix<doubleType>(3, 3, 0.0);

                //Step 3.1.3: sBCF
                math::Matrix<doubleType> d_unit_sBCF_dsBCF = unit_sBCF.unitDerivative(sBCF);

                //Step 3.1.3.1: sBCF derivative with respect to time
                math::Matrix<doubleType> d_sBCF_dt = d_eBCF_dt.cross(uBCF) + eBCF.cross(d_uBCF_dt);
                math::Matrix<doubleType> d_unit_sBCF_dt = d_unit_sBCF_dsBCF * d_sBCF_dt;

                //Step 3.1.3.2: sBCF derivative with respect to rBCF
                math::Matrix<doubleType> d_eBCF_cross_d_rBCF_x = math::Matrix<doubleType>(3, 3, std::vector<doubleType>({0., -d_eBCF_drBCF(2, 0), d_eBCF_drBCF(1, 0), d_eBCF_drBCF(2, 0), 0.0, -d_eBCF_drBCF(0, 0), -d_eBCF_drBCF(1, 0), d_eBCF_drBCF(0, 0), 0.0}));
                math::Matrix<doubleType> d_eBCF_cross_d_rBCF_y = math::Matrix<doubleType>(3, 3, std::vector<doubleType>({0., -d_eBCF_drBCF(2, 1), d_eBCF_drBCF(1, 1), d_eBCF_drBCF(2, 1), 0.0, -d_eBCF_drBCF(0, 1), -d_eBCF_drBCF(1, 1), d_eBCF_drBCF(0, 1), 0.0}));
                math::Matrix<doubleType> d_eBCF_cross_d_rBCF_z = math::Matrix<doubleType>(3, 3, std::vector<doubleType>({0., -d_eBCF_drBCF(2, 2), d_eBCF_drBCF(1, 2), d_eBCF_drBCF(2, 2), 0.0, -d_eBCF_drBCF(0, 2), -d_eBCF_drBCF(1, 2), d_eBCF_drBCF(0, 2), 0.0}));
                math::Tensor<doubleType> d_eBCF_cross_d_rBCF = math::Tensor<doubleType>(3, std::vector< math::Matrix<doubleType> >({ d_eBCF_cross_d_rBCF_x , d_eBCF_cross_d_rBCF_y , d_eBCF_cross_d_rBCF_z }));
                
                math::Matrix<doubleType> d_sBCF_d_rBCF = d_eBCF_cross_d_rBCF.bullet2_vector(uBCF) + eBCFcross * d_uBCF_d_rBCF;

                math::Matrix<doubleType> d_unit_sBCF_d_rBCF = d_unit_sBCF_dsBCF * d_sBCF_d_rBCF;
                
                //Step 3.1.3.3: sBCF derivative with respect to vBCF
                math::Matrix<doubleType> d_sBCF_dvBCF = eBCFcross * d_uBCF_dvBCF;
                math::Matrix<doubleType> d_unit_sBCF_dvBCF = d_unit_sBCF_dsBCF * d_sBCF_dvBCF;

                //Step 3.2: derivative of the transformation matrix
                //Step 3.2.1: with respect to time
                this->dR_from_BCF_to_Topocentric_dt = d_unit_sBCF_dt.transpose().vert_cat(d_unit_eBCF_dt.transpose().vert_cat(d_unit_uBCF_dt.transpose()));
                this->dR_from_Topocentric_to_BCF_dt = this->dR_from_BCF_to_Topocentric_dt.transpose();

                //Step 3.2.2: with respect to position
                this->dR_from_BCF_to_Topocentric_drBCF = math::Tensor<doubleType>(3, 3, 3, std::vector< math::Matrix<doubleType> >({d_unit_sBCF_d_rBCF, d_unit_eBCF_d_rBCF, d_unit_uBCF_d_rBCF}));
                //this->dR_from_Topocentric_to_BCF_drBCF = this->dR_from_BCF_to_Topocentric_drBCF.transpose();      I don't know how to do this!

                //Step 3.2.3: with respect to velocity
                this->dR_from_Topocentric_to_BCF_dvBCF = math::Tensor<doubleType>(3, 3, 3, std::vector< math::Matrix<doubleType> >(3, math::Matrix<doubleType>(3, 1, 0.0)));
                this->dR_from_BCF_to_Topocentric_dvBCF = math::Tensor<doubleType>(3, 3, 3, std::vector< math::Matrix<doubleType> >(3, math::Matrix<doubleType>(3, 1, 0.0)));
            }//end derivatives
        }//end construct_Topocentric_to_BCF_rotation()

        void frame::construct_Polar_to_BCF_rotation(const doubleType& ETepoch,
            const bool& GenerateDerivatives)
        {
            doubleType days_since_reference_epoch = ETepoch / 86400.0 - 51544.5;

            if (GenerateDerivatives)
            {

            }
        }//end construct_Polar_to_BCF_rotation()

        void frame::construct_SAM_to_J2000BCI_rotation(const doubleType& ETepoch,
            math::Matrix<doubleType> R_to_Sun,
            math::Matrix<doubleType> dR_to_Sun_dt,
            const bool& GenerateDerivatives)
        {
            doubleType r = R_to_Sun.norm();
            doubleType dr_dt = (R_to_Sun(0) * dR_to_Sun_dt(0) + R_to_Sun(1) * dR_to_Sun_dt(1) + R_to_Sun(2) * dR_to_Sun_dt(2)) / r;
            math::Matrix<doubleType> R = R_to_Sun / r;
            math::Matrix<doubleType> Rdot = dR_to_Sun_dt / r - R_to_Sun * dr_dt / r / r;
            doubleType rx = R_to_Sun(0);
            doubleType ry = R_to_Sun(1);
            doubleType rz = R_to_Sun(2);
            doubleType drx_dt = dR_to_Sun_dt(0);
            doubleType dry_dt = dR_to_Sun_dt(1);
            doubleType drz_dt = dR_to_Sun_dt(2);

            math::Matrix<doubleType> khat(3, 1, std::vector<doubleType>({ 0.0, 0.0, 1.0 }));
            math::Matrix<doubleType> X = R;
            math::Matrix<doubleType> Y = R.cross(khat);
            math::Matrix<doubleType> Z = X.cross(Y);

            math::Matrix<doubleType> Xhat = X.unitize();
            math::Matrix<doubleType> Yhat = Y.unitize();
            math::Matrix<doubleType> Zhat = Z.unitize();

            //rotation
            this->R_from_J2000BCI_to_SAM = math::Matrix<doubleType>(3, 3, 0.0);
            this->R_from_J2000BCI_to_SAM(0, 0) = Xhat(0);
            this->R_from_J2000BCI_to_SAM(0, 1) = Xhat(1);
            this->R_from_J2000BCI_to_SAM(0, 2) = Xhat(2);
            this->R_from_J2000BCI_to_SAM(1, 0) = Yhat(0);
            this->R_from_J2000BCI_to_SAM(1, 1) = Yhat(1);
            this->R_from_J2000BCI_to_SAM(1, 2) = Yhat(2);
            this->R_from_J2000BCI_to_SAM(2, 0) = Zhat(0);
            this->R_from_J2000BCI_to_SAM(2, 1) = Zhat(1);
            this->R_from_J2000BCI_to_SAM(2, 2) = Zhat(2);

            this->R_from_SAM_to_J2000BCI = this->R_from_J2000BCI_to_SAM.transpose();

            //time derivative
            doubleType rxy3 = sqrt(rx*rx + ry * ry) * sqrt(rx*rx + ry * ry) * sqrt(rx*rx + ry * ry);
            doubleType r3 = r * r * r;
            doubleType rx2 = rx * rx;
            doubleType rx3 = rx2 * rx;
            doubleType rx4 = rx2 * rx2;
            doubleType ry2 = ry * ry;
            doubleType ry3 = ry2 * ry;
            doubleType ry4 = ry2 * ry2;
            doubleType rz2 = rz * rz;
            doubleType A = rx2 * rz2 + ry2 * rz2 + (rx2 + ry2)*(rx2 + ry2);
            doubleType dA_drx = 2.0 * rx * rz2 + 4.0 * rx3 + 4.0 * rx * ry2;
            doubleType dA_dry = 2.0 * ry * rz2 + 4.0 * rx2 * ry + 4.0 * ry3;
            doubleType dA_drz = 2.0 * rx2 * rz + 2.0 * ry2 * rz;
            doubleType dA_dt = dA_drx * drx_dt + dA_dry * dry_dt + dA_drz * drz_dt;

            this->dR_from_J2000BCI_to_SAM_dt = math::Matrix<doubleType>(3, 3, 0.0);
            this->dR_from_J2000BCI_to_SAM_dt(0, 0) = (-rx * ry*dry_dt - rx * rz*drz_dt + ry2 * drx_dt + rz2 * drx_dt) / r3;
            this->dR_from_J2000BCI_to_SAM_dt(0, 1) = (rx2 * dry_dt - rx * ry*drx_dt - ry * rz*drz_dt + rz2 * dry_dt) / r3;
            this->dR_from_J2000BCI_to_SAM_dt(0, 2) = (rx2 * drz_dt - rx * rz*drx_dt + ry2 * drz_dt - ry * rz*dry_dt) / r3;
            this->dR_from_J2000BCI_to_SAM_dt(1, 0) = -rx * (ry * drx_dt - rx * dry_dt) / rxy3;
            this->dR_from_J2000BCI_to_SAM_dt(1, 1) = -ry * (ry * drx_dt - rx * dry_dt) / rxy3;
            this->dR_from_J2000BCI_to_SAM_dt(1, 2) = 0.0;
            this->dR_from_J2000BCI_to_SAM_dt(2, 0) = -(-(rx*drz_dt + rz * drx_dt)*((rx2 + ry2)*(rx2 + ry2) + rx2 * rz2 + ry2 * rz2) + (2 * (rx*drx_dt + ry * dry_dt)*(rx2 + ry2) + rx2 * rz*drz_dt + rx * rz2 * drx_dt + ry2 * rz*drz_dt + ry * rz2 * dry_dt)*rx*rz) / pow((rx2 + ry2)*(rx2 + ry2) + rx2 * rz2 + ry2 * rz2, 3.0/2.0);
            this->dR_from_J2000BCI_to_SAM_dt(2, 1) = -(-(ry*drz_dt + rz * dry_dt)*((rx2 + ry2)*(rx2 + ry2) + rx2 * rz2 + ry2 * rz2) + (2 * (rx*drx_dt + ry * dry_dt)*(rx2 + ry2) + rx2 * rz*drz_dt + rx * rz2 * drx_dt + ry2 * rz*drz_dt + ry * rz2 * dry_dt)*ry*rz) / pow((rx2 + ry2)*(rx2 + ry2) + rx2 * rz2 + ry2 * rz2, 3.0/2.0);
            this->dR_from_J2000BCI_to_SAM_dt(2, 2) = -(-rx2 * drz_dt + rx * rz*drx_dt - ry2 * drz_dt + ry * rz*dry_dt)*rz / ((rx2 + ry2 + rz2)*sqrt(rx4 + 2 * rx2 * ry2 + rx2 * rz2 + ry4 + ry2 * rz2));
            
            this->dR_from_SAM_to_J2000BCI_dt = this->dR_from_J2000BCI_to_SAM_dt.transpose();


#ifdef AD_INSTRUMENTATION
            std::ofstream out("tests/SAM_time_derivatives.csv", std::ios::trunc);
            out << "i, j, Analytical, Algorithmic, Absolute error, Relative error, an/al, al/an" << std::endl;
            for (size_t i : {0, 1, 2})
            {
                for (size_t j : {0, 1, 2})
                {
                    double algorithmic = R_from_J2000BCI_to_SAM(i, j).getDerivative(394) / 2780.3088828875989;
                    double analytical = dR_from_J2000BCI_to_SAM_dt(i, j)_GETVALUE;
                    double abserror = analytical - algorithmic;
                    double relerror = abserror / algorithmic;
                    double an_al = analytical / algorithmic;
                    double al_an = algorithmic / analytical;
                    out << i << "," << j << "," << analytical << "," << algorithmic << "," << abserror << "," << relerror << "," << an_al << "," << al_an << std::endl;
                }
            }
            out.close();
#endif

        }//end construct_SAM_to_J2000BCI_rotation()

        void frame::construct_ObjectReferenced_to_ICRF_rotation(const doubleType& ETepoch,
            math::Matrix<doubleType> referenceVector,
            math::Matrix<doubleType> dreferenceVector_dt,
            const bool& GenerateDerivatives)
        {
            math::Matrix<doubleType> R = referenceVector.getSubMatrix1D(0, 2);
            math::Matrix<doubleType> V = referenceVector.getSubMatrix1D(3, 5);
            math::Matrix<doubleType> dRdt = dreferenceVector_dt.getSubMatrix1D(0, 2);
            math::Matrix<doubleType> dVdt = dreferenceVector_dt.getSubMatrix1D(3, 5);
            doubleType r = R.norm();
            doubleType v = V.norm();
            doubleType dr_dt = (R(0) * dRdt(0) + R(1) * dRdt(1) + R(2) * dRdt(2)) / r;
            doubleType dv_dt = (V(0) * dVdt(0) + V(1) * dVdt(1) + V(2) * dVdt(2)) / v;
            math::Matrix<doubleType> Rhat = R / r;
            math::Matrix<doubleType> dRhat_dt = dRdt / r - R * dr_dt / r / r;
            math::Matrix<doubleType> Vhat = V / v;
            math::Matrix<doubleType> dVhat_dt = dVdt / v - V * dr_dt / v / v;
            math::Matrix<doubleType> H = R.cross(V);
            doubleType h = H.norm();
            math::Matrix<doubleType> Hhat = H / h;

            doubleType rx = referenceVector(0);
            doubleType ry = referenceVector(1);
            doubleType rz = referenceVector(2);
            doubleType vx = referenceVector(3);
            doubleType vy = referenceVector(4);
            doubleType vz = referenceVector(5);
            doubleType drx_dt = dreferenceVector_dt(0);
            doubleType dry_dt = dreferenceVector_dt(1);
            doubleType drz_dt = dreferenceVector_dt(2);
            doubleType dvx_dt = dreferenceVector_dt(3);
            doubleType dvy_dt = dreferenceVector_dt(4);
            doubleType dvz_dt = dreferenceVector_dt(5);

            math::Matrix<doubleType> Xhat = Rhat;
            math::Matrix<doubleType> Yhat = Vhat;
            math::Matrix<doubleType> Zhat = Hhat;

            //rotation
            this->R_from_ICRF_to_ObjectReferenced = math::Matrix<doubleType>(3, 3, 0.0);
            this->R_from_ICRF_to_ObjectReferenced(0, 0) = Xhat(0);
            this->R_from_ICRF_to_ObjectReferenced(0, 1) = Xhat(1);
            this->R_from_ICRF_to_ObjectReferenced(0, 2) = Xhat(2);
            this->R_from_ICRF_to_ObjectReferenced(1, 0) = Yhat(0);
            this->R_from_ICRF_to_ObjectReferenced(1, 1) = Yhat(1);
            this->R_from_ICRF_to_ObjectReferenced(1, 2) = Yhat(2);
            this->R_from_ICRF_to_ObjectReferenced(2, 0) = Zhat(0);
            this->R_from_ICRF_to_ObjectReferenced(2, 1) = Zhat(1);
            this->R_from_ICRF_to_ObjectReferenced(2, 2) = Zhat(2);

            this->R_from_ObjectReferenced_to_ICRF = this->R_from_ICRF_to_ObjectReferenced.transpose();

            //time derivative
            doubleType r3 = r * r * r;
            doubleType v3 = v * v * v;

            math::Matrix<doubleType> dZhat_dt = R.cross(dVdt) / h - Zhat / h * (R.cross(dVdt).dot(Zhat));

            this->dR_from_ICRF_to_ObjectReferenced_dt = math::Matrix<doubleType>(3, 3, 0.0);
            this->dR_from_ICRF_to_ObjectReferenced_dt(0, 0) = ((ry*ry + rz*rz) * drx_dt - rx * ry * dry_dt - rx * rz * drz_dt) / r3;
            this->dR_from_ICRF_to_ObjectReferenced_dt(0, 1) = (-rx * ry * drx_dt + (rx*rx + rz*rz) * dry_dt - ry * rz * drz_dt) / r3;
            this->dR_from_ICRF_to_ObjectReferenced_dt(0, 2) = (-rx * rz * drx_dt - ry * rz * dry_dt + (rx*rx + ry*ry) * drz_dt) / r3;
            this->dR_from_ICRF_to_ObjectReferenced_dt(1, 0) = ((vy*vy + vz * vz) * dvx_dt - vx * vy * dvy_dt - vx * vz * dvz_dt) / v3;
            this->dR_from_ICRF_to_ObjectReferenced_dt(1, 1) = (-vx * vy * dvx_dt + (vx*vx + vz * vz) * dvy_dt - vy * vz * dvz_dt) / v3;
            this->dR_from_ICRF_to_ObjectReferenced_dt(1, 2) = (-vx * vz * dvx_dt - vy * vz * dvy_dt + (vx*vx + vy * vy) * dvz_dt) / v3;
            this->dR_from_ICRF_to_ObjectReferenced_dt(2, 0) = dZhat_dt(0);
            this->dR_from_ICRF_to_ObjectReferenced_dt(2, 1) = dZhat_dt(1);
            this->dR_from_ICRF_to_ObjectReferenced_dt(2, 2) = dZhat_dt(2);

            this->dR_from_ObjectReferenced_to_ICRF_dt = this->dR_from_ICRF_to_ObjectReferenced_dt.transpose();


#ifdef AD_INSTRUMENTATION
            std::ofstream out("tests/ObjectReferenced_time_derivatives.csv", std::ios::trunc);
            out << "i, j, Analytical, Algorithmic, Absolute error, Relative error, an/al, al/an" << std::endl;
            for (size_t i : {0, 1, 2})
            {
                for (size_t j : {0, 1, 2})
                {
                    double algorithmic = R_from_ICRF_to_ObjectReferenced(i, j).getDerivative(7) / 13713.35828;
                    double analytical = dR_from_ICRF_to_ObjectReferenced_dt(i, j)_GETVALUE;
                    double abserror = analytical - algorithmic;
                    double relerror = abserror / algorithmic;
                    double an_al = analytical / algorithmic;
                    double al_an = algorithmic / analytical;
                    out << i << "," << j << "," << analytical << "," << algorithmic << "," << abserror << "," << relerror << "," << an_al << "," << al_an << std::endl;
                }
            }
            out.close();
#endif

        }//end construct_ObjectReferenced_to_ICRF_rotation()

        //construct rotation frames for J2000
        void frame::construct_rotation_matrices_J2000(const bool& GenerateDerivatives)
        {
            //construct the J2000 rotations, all of whom have zero derivatives
            doubleType delta0 = this->delta0;
            doubleType alpha0 = this->alpha0;
            math::Matrix<doubleType> Rx(3, (math::PIover2 - delta0), math::MatrixType::Rxhat);
            math::Matrix<doubleType> Rz(3, (math::PIover2 + alpha0), math::MatrixType::Rzhat);

            this->R_from_J2000BCI_to_ICRF = Rz*Rx;

            this->R_from_ICRF_to_J2000BCI = this->R_from_J2000BCI_to_ICRF.transpose();
        }//end construct_rotation_matrices_J2000()


        math::Matrix<doubleType> frame::get_R(const ReferenceFrame& OriginalReferenceFrame, const ReferenceFrame& RotatedReferenceFrame) const
        {
            math::Matrix<doubleType> R;

            if (OriginalReferenceFrame == RotatedReferenceFrame)
            {
                R = math::Matrix<doubleType>(3, math::identity);
            }

            else if (OriginalReferenceFrame == EMTG::ReferenceFrame::ICRF)
            {
                if (RotatedReferenceFrame == EMTG::ReferenceFrame::TrueOfDate_BCI)
                {
                    R = this->R_from_ICRF_to_BCI;
                }
                else if (RotatedReferenceFrame == EMTG::ReferenceFrame::TrueOfDate_BCF)
                {
                    R = this->R_from_ICRF_to_BCF;
                }
                else if (RotatedReferenceFrame == EMTG::ReferenceFrame::J2000_BCI)
                {
                    R = this->R_from_ICRF_to_J2000BCI;
                }
                else if (RotatedReferenceFrame == EMTG::ReferenceFrame::J2000_BCF)
                {
                    R = this->R_from_ICRF_to_J2000BCF;
                }
                else if (RotatedReferenceFrame == EMTG::ReferenceFrame::PrincipleAxes)
                {
                    R = this->R_from_BCF_to_PrincipleAxes * this->get_R(EMTG::ReferenceFrame::ICRF, EMTG::ReferenceFrame::TrueOfDate_BCF);
                }
                else if (RotatedReferenceFrame == EMTG::ReferenceFrame::Topocentric)
                {
                    R = this->R_from_BCF_to_Topocentric * this->get_R(EMTG::ReferenceFrame::ICRF, EMTG::ReferenceFrame::TrueOfDate_BCF);
                }
                else if (RotatedReferenceFrame == ReferenceFrame::SAM)
                {
                    R = this->R_from_J2000BCI_to_SAM * this->R_from_ICRF_to_J2000BCI;
                }
                else if (RotatedReferenceFrame == ReferenceFrame::ObjectReferenced)
                {
                    R = this->R_from_ICRF_to_ObjectReferenced;
                }
            }
            else if (OriginalReferenceFrame == EMTG::ReferenceFrame::TrueOfDate_BCI)
            {
                if (RotatedReferenceFrame == EMTG::ReferenceFrame::ICRF)
                {
                    R = this->R_from_BCI_to_ICRF;
                }
                else if (RotatedReferenceFrame == EMTG::ReferenceFrame::TrueOfDate_BCF)
                {
                    R = this->R_from_BCI_to_BCF;
                }
            }
            else if (OriginalReferenceFrame == EMTG::ReferenceFrame::TrueOfDate_BCF)
            {
                if (RotatedReferenceFrame == EMTG::ReferenceFrame::ICRF)
                {
                    R = this->R_from_BCF_to_ICRF;
                }
                else if (RotatedReferenceFrame == EMTG::ReferenceFrame::TrueOfDate_BCI)
                {
                    R = this->R_from_BCF_to_BCI;
                }
                else
                {
                    R = this->get_R(EMTG::ReferenceFrame::ICRF, RotatedReferenceFrame) * this->R_from_BCF_to_ICRF;
                }
            }
            else if (OriginalReferenceFrame == EMTG::ReferenceFrame::J2000_BCI)
            {
                if (RotatedReferenceFrame == EMTG::ReferenceFrame::ICRF)
                {
                    R = this->R_from_J2000BCI_to_ICRF;
                }
                else if (RotatedReferenceFrame == EMTG::ReferenceFrame::J2000_BCF)
                {
                    R = this->R_from_J2000BCI_to_J2000BCF;
                }
            }
            else if (OriginalReferenceFrame == EMTG::ReferenceFrame::J2000_BCF)
            {
                if (RotatedReferenceFrame == EMTG::ReferenceFrame::ICRF)
                {
                    R = this->R_from_J2000BCF_to_ICRF;
                }
                else if (RotatedReferenceFrame == EMTG::ReferenceFrame::J2000_BCI)
                {
                    R = this->R_from_J2000BCF_to_J2000BCI;
                }
            }
            else if (OriginalReferenceFrame == EMTG::ReferenceFrame::PrincipleAxes)
            {
                if (RotatedReferenceFrame == EMTG::ReferenceFrame::TrueOfDate_BCF)
                {
                    R = this->R_from_PrincipleAxes_to_BCF;
                }
                else
                {
                    R = this->get_R(EMTG::ReferenceFrame::TrueOfDate_BCF, RotatedReferenceFrame) * this->R_from_PrincipleAxes_to_BCF;
                }
            }
            else if (OriginalReferenceFrame == EMTG::ReferenceFrame::Topocentric)
            {
                if (RotatedReferenceFrame == EMTG::ReferenceFrame::TrueOfDate_BCF)
                {
                    R = this->R_from_Topocentric_to_BCF;
                }
                else
                {
                    R = this->get_R(EMTG::ReferenceFrame::TrueOfDate_BCF, RotatedReferenceFrame) * this->R_from_Topocentric_to_BCF;
                }
            }
            else if (OriginalReferenceFrame == ReferenceFrame::SAM)
            {
                R = this->get_R(ReferenceFrame::J2000_BCI, RotatedReferenceFrame) * this->R_from_SAM_to_J2000BCI;
            }
            else if (OriginalReferenceFrame == ReferenceFrame::ObjectReferenced)
            {
                R = this->get_R(ReferenceFrame::ICRF, RotatedReferenceFrame) * this->R_from_ObjectReferenced_to_ICRF;
            }

            return R;
        }//end get_R()

        math::Matrix<doubleType> frame::get_dRdt(const ReferenceFrame& OriginalReferenceFrame, const ReferenceFrame& RotatedReferenceFrame) const
        {
            math::Matrix<doubleType> Rdot;

            if (OriginalReferenceFrame == RotatedReferenceFrame)
            {
                Rdot = math::Matrix<doubleType>(3, 3, 0.0);
            }

            if (OriginalReferenceFrame == EMTG::ReferenceFrame::ICRF)
            {
                if (RotatedReferenceFrame == EMTG::ReferenceFrame::TrueOfDate_BCI)
                {
                    Rdot = this->dR_from_ICRF_to_BCI_dt;
                }
                else if (RotatedReferenceFrame == EMTG::ReferenceFrame::TrueOfDate_BCF)
                {
                    Rdot = this->dR_from_ICRF_to_BCF_dt;
                }
                else if (RotatedReferenceFrame == EMTG::ReferenceFrame::J2000_BCI)
                {
                    Rdot = math::Matrix<doubleType>(3, 3, 0.0);
                }
                else if (RotatedReferenceFrame == EMTG::ReferenceFrame::J2000_BCF)
                {
                    Rdot = this->dR_from_ICRF_to_J2000BCF_dt;
                }
                else if (RotatedReferenceFrame == EMTG::ReferenceFrame::PrincipleAxes)
                {
                    Rdot = this->dR_from_BCF_to_PrincipleAxes_dt * this->get_R(EMTG::ReferenceFrame::ICRF, EMTG::ReferenceFrame::TrueOfDate_BCF)
                         + this->R_from_BCF_to_PrincipleAxes * this->get_dRdt(EMTG::ReferenceFrame::ICRF, EMTG::ReferenceFrame::TrueOfDate_BCF);
                }
                else if (RotatedReferenceFrame == EMTG::ReferenceFrame::Topocentric)
                {
                    Rdot = this->dR_from_BCF_to_Topocentric_dt * this->get_R(EMTG::ReferenceFrame::ICRF, EMTG::ReferenceFrame::TrueOfDate_BCF)
                        + this->R_from_BCF_to_Topocentric * this->get_dRdt(EMTG::ReferenceFrame::ICRF, EMTG::ReferenceFrame::TrueOfDate_BCF);
                }
                else if (RotatedReferenceFrame == ReferenceFrame::SAM)
                {
                    Rdot = this->dR_from_J2000BCI_to_SAM_dt * this->R_from_ICRF_to_J2000BCI;
                }
                else if (RotatedReferenceFrame == ReferenceFrame::ObjectReferenced)
                {
                    Rdot = this->dR_from_ICRF_to_ObjectReferenced_dt;
                }
            }
            else if (OriginalReferenceFrame == EMTG::ReferenceFrame::TrueOfDate_BCI)
            {
                if (RotatedReferenceFrame == EMTG::ReferenceFrame::ICRF)
                {
                    Rdot = this->dR_from_BCI_to_ICRF_dt;
                }
                else if (RotatedReferenceFrame == EMTG::ReferenceFrame::TrueOfDate_BCF)
                {
                    Rdot = this->dR_from_BCI_to_BCF_dt;
                }
            }
            else if (OriginalReferenceFrame == EMTG::ReferenceFrame::TrueOfDate_BCF)
            {
                if (RotatedReferenceFrame == EMTG::ReferenceFrame::ICRF)
                {
                    Rdot = this->dR_from_BCF_to_ICRF_dt;
                }
                else if (RotatedReferenceFrame == EMTG::ReferenceFrame::TrueOfDate_BCI)
                {
                    Rdot = this->dR_from_BCF_to_BCI_dt;
                }
                else if (RotatedReferenceFrame == EMTG::ReferenceFrame::Topocentric)
                {
                    Rdot = this->dR_from_BCF_to_Topocentric_dt;
                }
            }
            else if (OriginalReferenceFrame == EMTG::ReferenceFrame::J2000_BCI)
            {
                if (RotatedReferenceFrame == EMTG::ReferenceFrame::ICRF)
                {
                    Rdot = math::Matrix<doubleType>(3, 3, 0.0);
                }
                else if (RotatedReferenceFrame == EMTG::ReferenceFrame::J2000_BCF)
                {
                    Rdot = this->dR_from_J2000BCI_to_J2000BCF_dt;
                }
            }
            else if (OriginalReferenceFrame == EMTG::ReferenceFrame::J2000_BCF)
            {
                if (RotatedReferenceFrame == EMTG::ReferenceFrame::ICRF)
                {
                    Rdot = this->dR_from_J2000BCF_to_ICRF_dt;
                }
                else if (RotatedReferenceFrame == EMTG::ReferenceFrame::J2000_BCI)
                {
                    Rdot = this->dR_from_J2000BCF_to_J2000BCI_dt;
                }
            }
            else if (OriginalReferenceFrame == EMTG::ReferenceFrame::PrincipleAxes)
            {
                if (RotatedReferenceFrame == EMTG::ReferenceFrame::TrueOfDate_BCF)
                {
                    Rdot = this->dR_from_PrincipleAxes_to_BCF_dt;
                }
                else
                {
                    Rdot = this->get_dRdt(EMTG::ReferenceFrame::TrueOfDate_BCF, RotatedReferenceFrame) * this->R_from_PrincipleAxes_to_BCF
                         + this->get_R(EMTG::ReferenceFrame::TrueOfDate_BCF, RotatedReferenceFrame) * this->dR_from_PrincipleAxes_to_BCF_dt;
                }
            }
            else if (OriginalReferenceFrame == EMTG::ReferenceFrame::Topocentric)
            {
                if (RotatedReferenceFrame == EMTG::ReferenceFrame::TrueOfDate_BCF)
                {
                    Rdot = this->dR_from_Topocentric_to_BCF_dt;
                }
                else
                {
                    Rdot = this->get_dRdt(EMTG::ReferenceFrame::TrueOfDate_BCF, RotatedReferenceFrame) * this->R_from_Topocentric_to_BCF
                        + this->get_R(EMTG::ReferenceFrame::TrueOfDate_BCF, RotatedReferenceFrame) * this->dR_from_Topocentric_to_BCF_dt;
                }
            }
            else if (OriginalReferenceFrame == ReferenceFrame::SAM)
            {
                Rdot = this->get_R(ReferenceFrame::J2000_BCI, RotatedReferenceFrame) * this->dR_from_SAM_to_J2000BCI_dt;
            }
            else if (OriginalReferenceFrame == ReferenceFrame::ObjectReferenced)
            {
                Rdot = this->get_R(ReferenceFrame::ICRF, RotatedReferenceFrame) * this->dR_from_ObjectReferenced_to_ICRF_dt;
            }

            return Rdot;
        }//end get_dRdt()

        //method to rotate a state, assumes that the frame class has already been initialized to the correct epoch
        void frame::rotate_frame_to_frame(const EMTG::ReferenceFrame& OriginalReferenceFrame,
            const math::Matrix<doubleType>& state,
            const math::Matrix<doubleType>& dstate_dt,
            const math::Matrix<doubleType>& referenceVector,
            const math::Matrix<doubleType>& dreferenceVector_dt,
            const EMTG::ReferenceFrame& RotatedReferenceFrame,
            math::Matrix<doubleType>& Rstate,
            math::Matrix<doubleType>& dRstate_dstate,
            math::Matrix<doubleType>& dRstate_dt,
            math::Matrix<doubleType>& dRstate_dreferenceVector,
            const bool& GenerateDerivatives)
        {
            if (OriginalReferenceFrame == RotatedReferenceFrame)
            {
                Rstate = state;
                if (GenerateDerivatives)
                {
                    dRstate_dt = dreferenceVector_dt;

                    if (state == referenceVector)
                    {
                        dRstate_dreferenceVector = math::Matrix<doubleType>(3, math::MatrixType::identity);
                    }
                    else
                    {
                        dRstate_dreferenceVector = math::Matrix<doubleType>(3, 3, 0.0);
                    }
                }
                return;
            }
            else
            {
                math::Matrix<doubleType> R = this->get_R(OriginalReferenceFrame, RotatedReferenceFrame);
                math::Matrix<doubleType> Rdot = this->get_dRdt(OriginalReferenceFrame, RotatedReferenceFrame);

                Rstate = R * state;
                if (GenerateDerivatives)
                {
                    if (OriginalReferenceFrame == ReferenceFrame::ObjectReferenced || RotatedReferenceFrame == ReferenceFrame::ObjectReferenced)
                    {
                        dRstate_dt = R * dstate_dt + Rdot * state;
                        dRstate_dstate = R;
                    }
                    else
                    {
                        dRstate_dt = R * dstate_dt + (Rdot + this->dR_dreferenceVector.bullet2_vector(dreferenceVector_dt)) * state;
                        dRstate_dreferenceVector = this->dR_dreferenceVector.bullet2_vector(state);
                        if (state == referenceVector)
                        {
                            dRstate_dstate = R + dRstate_dreferenceVector;
                        }
                        else
                        {
                            dRstate_dstate = R;
                        }
                    }
                }

                return;
            }
        }//end rotate_frame_to_frame()

        //method to rotate a state if the frame class is not yet initialized
        void frame::rotate_frame_to_frame(
            const EMTG::ReferenceFrame& OriginalReferenceFrame,
            const math::Matrix<doubleType>& state,
            const math::Matrix<doubleType>& dstate_dt,
            const math::Matrix<doubleType>& referenceVector,
            const math::Matrix<doubleType>& dreferenceVector_dt,
            const EMTG::ReferenceFrame& RotatedReferenceFrame,
            math::Matrix<doubleType>& Rstate,
            math::Matrix<doubleType>& dRstate_dstate,
            math::Matrix<doubleType>& dRstate_dt,
            math::Matrix<doubleType>& dRstate_dreferenceVector,
            const doubleType& ETsecSinceJ2000,
            const bool& GenerateDerivatives)
        {
            this->construct_rotation_matrices_J2000();
            this->construct_rotation_matrices(ETsecSinceJ2000, GenerateDerivatives);

            if (OriginalReferenceFrame == ReferenceFrame::PrincipleAxes || RotatedReferenceFrame == ReferenceFrame::PrincipleAxes)
                this->construct_PrincipleAxes_to_BCF_rotation(ETsecSinceJ2000, GenerateDerivatives);
            if (OriginalReferenceFrame == ReferenceFrame::Topocentric || RotatedReferenceFrame == ReferenceFrame::Topocentric)
                this->construct_Topocentric_to_BCF_rotation(ETsecSinceJ2000, referenceVector, GenerateDerivatives);
            if (OriginalReferenceFrame == ReferenceFrame::Polar || RotatedReferenceFrame == ReferenceFrame::Polar)
                this->construct_Polar_to_BCF_rotation(ETsecSinceJ2000, GenerateDerivatives);
            if (OriginalReferenceFrame == ReferenceFrame::SAM || RotatedReferenceFrame == ReferenceFrame::SAM)
                this->construct_SAM_to_J2000BCI_rotation(ETsecSinceJ2000, referenceVector, dreferenceVector_dt, GenerateDerivatives);
            if (OriginalReferenceFrame == ReferenceFrame::ObjectReferenced || RotatedReferenceFrame == ReferenceFrame::ObjectReferenced)
                this->construct_ObjectReferenced_to_ICRF_rotation(ETsecSinceJ2000, referenceVector, dreferenceVector_dt, GenerateDerivatives);

            this->rotate_frame_to_frame( OriginalReferenceFrame,
                                         state,
                                         dstate_dt,
                                         referenceVector,
                                         dreferenceVector_dt,
                                         RotatedReferenceFrame,
                                         Rstate,
                                         dRstate_dstate,
                                         dRstate_dt,
                                         dRstate_dreferenceVector,
                                         GenerateDerivatives);
        }//end rotate_frame_to_frame()

        //alias for if you want the rotation matrix but not the rest of the derivatives
        void frame::rotate_frame_to_frame(
            const EMTG::ReferenceFrame& OriginalReferenceFrame,
            const math::Matrix<doubleType>& state,
            const EMTG::ReferenceFrame& RotatedReferenceFrame,
            math::Matrix<doubleType>& Rstate,
            math::Matrix<doubleType>& dRstate_dstate,
            const doubleType& ETsecSinceJ2000)
        {
            //make up dummies for a lot of the entries that we won't use because we don't care about derivatives
            static math::Matrix<doubleType> referenceVector(3, 1, 0.0);
            static math::Matrix<doubleType> dstate_dt(3, 1, 0.0);
            static math::Matrix<doubleType> dreferenceVector_dt(3, 1, 0.0);
            static math::Matrix<doubleType> dRstate_dt(3, 1, 0.0);
            static math::Matrix<doubleType> dRstate_dreferenceVector(3, 3, 0.0);

            this->rotate_frame_to_frame(OriginalReferenceFrame,
                state,
                dstate_dt,
                referenceVector,
                dreferenceVector_dt,
                RotatedReferenceFrame,
                Rstate,
                dRstate_dstate,
                dRstate_dt,
                dRstate_dreferenceVector,
                ETsecSinceJ2000,
                false);

            dRstate_dstate = this->get_R(OriginalReferenceFrame, RotatedReferenceFrame);

        }//end rotate_frame_to_frame() if you only want the rotation matrix

        //alias for if you don't want the derivatives
        void frame::rotate_frame_to_frame(
            const EMTG::ReferenceFrame& OriginalReferenceFrame,
            const math::Matrix<doubleType>& state,
            const EMTG::ReferenceFrame& RotatedReferenceFrame,
            math::Matrix<doubleType>& Rstate,
            const doubleType& ETsecSinceJ2000)
        {
            //make up dummies for a lot of the entries that we won't use because we don't care about derivatives
            static math::Matrix<doubleType> referenceVector(3, 1, 0.0);
            static math::Matrix<doubleType> dstate_dt(3, 1, 0.0);
            static math::Matrix<doubleType> dreferenceVector_dt(3, 1, 0.0);
            static math::Matrix<doubleType> dRstate_dstate(3, 3, 0.0);
            static math::Matrix<doubleType> dRstate_dt(3, 1, 0.0);
            static math::Matrix<doubleType> dRstate_dreferenceVector(3, 3, 0.0);

            this->rotate_frame_to_frame(OriginalReferenceFrame,
                                        state,
                                        dstate_dt,
                                        referenceVector,
                                        dreferenceVector_dt,
                                        RotatedReferenceFrame,
                                        Rstate,
                                        dRstate_dstate,
                                        dRstate_dt,
                                        dRstate_dreferenceVector,
                                        ETsecSinceJ2000,
                                        false);
        }//end rotate_frame_to_frame() if you don't want to return derivatives


        void frame::rotate_frame_to_frame(
            const EMTG::ReferenceFrame& OriginalReferenceFrame,
            const math::Matrix<doubleType>& state,
            const math::Matrix<doubleType>& referenceVector,
            const EMTG::ReferenceFrame& RotatedReferenceFrame,
            math::Matrix<doubleType>& Rstate,
            const doubleType& ETsecSinceJ2000)
        {
            //make up dummies for a lot of the entries that we won't use because we don't care about derivatives
            static math::Matrix<doubleType> dstate_dt(3, 1, 0.0);
            static math::Matrix<doubleType> dreferenceVector_dt(3, 1, 0.0);
            static math::Matrix<doubleType> dRstate_dstate(3, 3, 0.0);
            static math::Matrix<doubleType> dRstate_dt(3, 1, 0.0);
            static math::Matrix<doubleType> dRstate_dreferenceVector(3, 3, 0.0);

            this->rotate_frame_to_frame(OriginalReferenceFrame,
                state,
                dstate_dt,
                referenceVector,
                dreferenceVector_dt,
                RotatedReferenceFrame,
                Rstate,
                dRstate_dstate,
                dRstate_dt,
                dRstate_dreferenceVector,
                ETsecSinceJ2000,
                false);
        }//end rotate_frame_to_frame() if you don't want to return derivatives, but you do have a reference vector


		// Getters for reference angles and rates
		doubleType frame::getAlpha(const doubleType& ReferenceEpochMJDSeconds)
		{
			doubleType days_since_reference_epoch = ReferenceEpochMJDSeconds / 86400.0 - 51544.5;
			doubleType centuries_since_reference_epoch = days_since_reference_epoch / 36525.0;
			return this->alpha0 + this->alphadot * centuries_since_reference_epoch;
		}

		doubleType frame::getDelta(const doubleType& ReferenceEpochMJDSeconds)
		{
			doubleType days_since_reference_epoch = ReferenceEpochMJDSeconds / 86400.0 - 51544.5;
			doubleType centuries_since_reference_epoch = days_since_reference_epoch / 36525.0;
			return this->delta0 + this->deltadot * centuries_since_reference_epoch;
		}

		doubleType frame::getW(const doubleType& ReferenceEpochMJDSeconds)
		{
			doubleType days_since_reference_epoch = ReferenceEpochMJDSeconds / 86400.0 - 51544.5;
			doubleType centuries_since_reference_epoch = days_since_reference_epoch / 36525.0;
			return this->W0 + this->Wdot * days_since_reference_epoch;
		}

		doubleType frame::getAlpha0()
		{
			return this->alpha0;
		}

		doubleType frame::getDelta0()
		{
			return this->delta0;
		}

		doubleType frame::getW0()
		{
			return this->W0;
		}

		doubleType frame::getAlphadot()
		{
			return this->alphadot;
		}

		doubleType frame::getDeltadot()
		{
			return this->deltadot;
		}

		doubleType frame::getWdot()
		{
			return this->Wdot;
		}

    }//close namespace Astrodynamics
} //close namespace EMTG