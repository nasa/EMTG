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
//templated 2/16/2017

#pragma once

#include <cmath>

#include "EMTG_Matrix.h"
#include "EMTG_Tensor.h"
#include "EMTG_math.h"
#include "EMTG_enums.h"

namespace EMTG 
{ 
    namespace Astrodynamics
    {
        class frame
        {
        public:
            //constructor
            frame(void);

            frame(const double& alpha0_in,
                const double& alphadot_in,
                const double& delta0_in,
                const double& deltadot_in,
                const double& W_in,
                const double& Wdot_in);

            frame(const double& alpha0_in,
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
                const double& PA3dot);

            //destructor
            virtual ~frame(void);

            //methods
            void initialize();
            void initialize(const double& alpha0,
                const double& alphadot,
                const double& delta0,
                const double& deltadot,
                const double& W,
                const double& Wdot,
                const double& PA1_0 = 0.0,
                const double& PA1dot = 0.0,
                const double& PA2_0 = 0.0,
                const double& PA2dot = 0.0,
                const double& PA3_0 = 0.0,
                const double& PA3dot = 0.0);
            void initialize_ecliptic();

            void rotate_frame_to_frame(
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
                const bool& GenerateDerivatives = false);

            void rotate_frame_to_frame(
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
                const doubleType& ReferenceEpochJD = 2451545.0,
                const bool& GenerateDerivatives = false);

            void rotate_frame_to_frame(
                const EMTG::ReferenceFrame& OriginalReferenceFrame,
                const math::Matrix<doubleType>& state,
                const EMTG::ReferenceFrame& RotatedReferenceFrame,
                math::Matrix<doubleType>& Rstate,
                math::Matrix<doubleType>& dRstate_dstate,
                const doubleType& ReferenceEpochJD);

            void rotate_frame_to_frame(
                const EMTG::ReferenceFrame& OriginalReferenceFrame,
                const math::Matrix<doubleType>& state,
                const EMTG::ReferenceFrame& RotatedReferenceFrame,
                math::Matrix<doubleType>& Rstate,
                const doubleType& ReferenceEpochJD = 2451545.0);
            
            void rotate_frame_to_frame(
                const EMTG::ReferenceFrame& OriginalReferenceFrame,
                const math::Matrix<doubleType>& state,
                const math::Matrix<doubleType>& referenceVector,
                const EMTG::ReferenceFrame& RotatedReferenceFrame,
                math::Matrix<doubleType>& Rstate,
                const doubleType& ReferenceEpochJD = 2451545.0);

            void construct_rotation_matrices(const doubleType& ETepoch, const bool& GenerateDerivatives = false);
            void construct_rotation_matrices_J2000(const bool& GenerateDerivatives = false);
            void construct_PrincipleAxes_to_BCF_rotation(const doubleType& ETepoch, const bool& GenerateDerivatives = false);
            void construct_Topocentric_to_BCF_rotation(const doubleType& ETepoch, math::Matrix<doubleType> rBCF, const bool& GenerateDerivatives = false);
            void construct_Polar_to_BCF_rotation(const doubleType& ETepoch, const bool& GenerateDerivatives = false);
            void construct_SAM_to_J2000BCI_rotation(const doubleType& ETepoch, math::Matrix<doubleType> R_to_Sun, math::Matrix<doubleType> dR_to_Sun_dt, const bool& GenerateDerivatives = false);
            void construct_ObjectReferenced_to_ICRF_rotation(const doubleType& ETepoch, math::Matrix<doubleType> referenceVector, math::Matrix<doubleType> dreferenceVector_dt, const bool& GenerateDerivatives = false);

            //get methods
            inline double get_alpha0() const { return this->alpha0; };
            inline double get_alphadot() const { return this->alphadot; };
            inline double get_delta0() const { return this->delta0; };
            inline double get_deltadot() const { return this->deltadot; };
            inline double getW0() const { return this->W0; };
            inline double getWdot() const { return this->Wdot; };
            inline math::Matrix<doubleType> get_xhat() const { return this->xhat; };
            inline math::Matrix<doubleType> get_yhat() const { return this->zhat; };
            inline math::Matrix<doubleType> get_zhat() const { return this->zhat; };
            inline math::Matrix<doubleType> get_R_from_ICRF_to_J2000BCI() const { return this->R_from_ICRF_to_J2000BCI; };
            inline math::Matrix<doubleType> get_R_from_J2000BCI_to_ICRF() const { return this->R_from_J2000BCI_to_ICRF; };
            inline math::Matrix<doubleType> get_R_from_J2000BCI_to_J2000BCF() const { return this->R_from_J2000BCI_to_J2000BCF; };
            inline math::Matrix<doubleType> get_R_from_J2000BCF_to_J2000BCI() const { return this->R_from_J2000BCF_to_J2000BCI; };
            inline math::Matrix<doubleType> get_R_from_J2000BCF_to_ICRF() const { return this->R_from_J2000BCF_to_ICRF; };
            inline math::Matrix<doubleType> get_R_from_ICRF_to_J2000BCF() const { return this->R_from_ICRF_to_J2000BCF; };
            inline math::Matrix<doubleType> get_R_from_BCI_to_ICRF     () const { return this->R_from_BCI_to_ICRF     ; };
            inline math::Matrix<doubleType> get_R_from_ICRF_to_BCI     () const { return this->R_from_ICRF_to_BCI     ; };
            inline math::Matrix<doubleType> get_R_from_BCI_to_BCF      () const { return this->R_from_BCI_to_BCF      ; };
            inline math::Matrix<doubleType> get_R_from_BCF_to_BCI      () const { return this->R_from_BCF_to_BCI      ; };
            inline math::Matrix<doubleType> get_R_from_BCF_to_ICRF     () const { return this->R_from_BCF_to_ICRF     ; };
            inline math::Matrix<doubleType> get_R_from_ICRF_to_BCF     () const { return this->R_from_ICRF_to_BCF     ; };
            inline math::Matrix<doubleType> get_dR_from_BCI_to_ICRF_dt () const { return this->dR_from_BCI_to_ICRF_dt ; };
            inline math::Matrix<doubleType> get_dR_from_ICRF_to_BCI_dt () const { return this->dR_from_ICRF_to_BCI_dt ; };
            inline math::Matrix<doubleType> get_dR_from_BCI_to_BCF_dt  () const { return this->dR_from_BCI_to_BCF_dt  ; };
            inline math::Matrix<doubleType> get_dR_from_BCF_to_BCI_dt  () const { return this->dR_from_BCF_to_BCI_dt  ; };
            inline math::Matrix<doubleType> get_dR_from_BCF_to_ICRF_dt () const { return this->dR_from_BCF_to_ICRF_dt ; };
            inline math::Matrix<doubleType> get_dR_from_ICRF_to_BCF_dt () const { return this->dR_from_ICRF_to_BCF_dt ; };
            inline math::Matrix<doubleType> get_dR_from_J2000BCI_to_J2000BCF_dt() const { return this->dR_from_J2000BCI_to_J2000BCF_dt; };
            inline math::Matrix<doubleType> get_dR_from_J2000BCF_to_J2000BCI_dt() const { return this->dR_from_J2000BCF_to_J2000BCI_dt; };
            inline math::Matrix<doubleType> get_dR_from_J2000BCF_to_ICRF_dt() const { return this->dR_from_J2000BCF_to_ICRF_dt; };
            inline math::Matrix<doubleType> get_dR_from_ICRF_to_J2000BCF_dt() const { return this->dR_from_ICRF_to_J2000BCF_dt; };

            inline math::Matrix<doubleType> get_R_from_PrincipleAxes_to_BCF() const { return this->R_from_PrincipleAxes_to_BCF; };
            inline math::Matrix<doubleType> get_R_from_BCF_to_PrincipleAxes() const { return this->R_from_BCF_to_PrincipleAxes; };
            inline math::Matrix<doubleType> get_dR_from_PrincipleAxes_to_BCF_dt() const { return this->dR_from_PrincipleAxes_to_BCF_dt; };
            inline math::Matrix<doubleType> get_dR_from_BCF_to_PrincipleAxes_dt() const { return this->dR_from_BCF_to_PrincipleAxes_dt; };
            inline math::Matrix<doubleType> get_R_from_Topocentric_to_BCF() const { return this->R_from_Topocentric_to_BCF; };
            inline math::Matrix<doubleType> get_R_from_BCF_to_Topocentric() const { return this->R_from_BCF_to_Topocentric; };
            inline math::Matrix<doubleType> get_dR_from_Topocentric_to_BCF_dt() const { return this->dR_from_Topocentric_to_BCF_dt; };
            inline math::Matrix<doubleType> get_dR_from_BCF_to_Topocentric_dt() const { return this->dR_from_BCF_to_Topocentric_dt; };

            inline math::Tensor<doubleType> get_dR_from_Topocentric_to_BCF_drBCF() const {return this->dR_from_Topocentric_to_BCF_drBCF; };
            inline math::Tensor<doubleType> get_dR_from_BCF_to_Topocentric_drBCF() const {return this->dR_from_BCF_to_Topocentric_drBCF; };

            inline math::Matrix<doubleType> get_R_from_Polar_to_BCF() const { return this->R_from_Polar_to_BCF; };
            inline math::Matrix<doubleType> get_R_from_BCF_to_Polar() const { return this->R_from_BCF_to_Polar; };
            inline math::Matrix<doubleType> get_dR_from_Polar_to_BCF_dt() const { return this->dR_from_Polar_to_BCF_dt; };
            inline math::Matrix<doubleType> get_dR_from_BCF_to_Polar_dt() const { return this->dR_from_BCF_to_Polar_dt; };


            inline math::Matrix<doubleType> get_R_from_SAM_to_J2000BCI() const { return this->R_from_SAM_to_J2000BCI; };
            inline math::Matrix<doubleType> get_R_from_J2000BCI_to_SAM() const {return this->R_from_J2000BCI_to_SAM; };
            inline math::Matrix<doubleType> get_dR_from_SAM_to_J2000BCI_dt() const {return this->dR_from_SAM_to_J2000BCI_dt; };
            inline math::Matrix<doubleType> get_dR_from_J2000BCI_to_SAM_dt() const { return this->dR_from_J2000BCI_to_SAM_dt; };

            inline math::Matrix<doubleType> get_R_from_ObjectReferenced_to_ICRF() const { return this->R_from_ObjectReferenced_to_ICRF; };
            inline math::Matrix<doubleType> get_R_from_ICRF_to_ObjectReferenced() const { return this->R_from_ICRF_to_ObjectReferenced; };
            inline math::Matrix<doubleType> get_dR_from_ObjectReferenced_to_ICRF_dt() const { return this->dR_from_ObjectReferenced_to_ICRF_dt; };
            inline math::Matrix<doubleType> get_dR_from_ICRF_to_ObjectReferenced_dt() const { return this->dR_from_ICRF_to_ObjectReferenced_dt; };

            math::Matrix<doubleType> get_R(const ReferenceFrame& OriginalReferenceFrame, const ReferenceFrame& RotatedReferenceFrame) const;
            math::Matrix<doubleType> get_dRdt(const ReferenceFrame& OriginalReferenceFrame, const ReferenceFrame& RotatedReferenceFrame) const;
            math::Tensor<doubleType> get_dR_dreferenceVector() const { return this->dR_dreferenceVector; }

            //set methods
            inline void set_semi_axis_a(const double& semi_axis_a) { this->semi_axis_a = semi_axis_a; }
            inline void set_semi_axis_b(const double& semi_axis_b) { this->semi_axis_b = semi_axis_b; }
            inline void set_semi_axis_c(const double& semi_axis_c) { this->semi_axis_c = semi_axis_c; }

			// get methods

			// get the angles
			doubleType getAlpha(const doubleType& ReferenceEpochJD);
			doubleType getDelta(const doubleType& ReferenceEpochJD);
			doubleType getW(const doubleType& ReferenceEpochJD);
			doubleType getAlpha0();
			doubleType getDelta0();
			doubleType getW0();
			doubleType getAlphadot();
			doubleType getDeltadot();
			doubleType getWdot();
                                        
        private:
            //fields
            double alpha0, alphadot, delta0, deltadot, W0, Wdot; //BCI and BCF rotation angles
            double PA1_0, PA2_0, PA3_0; //principle axes offset from BCF, defaults to zero
            double PA1dot, PA2dot, PA3dot; //rates for the principle axes offset
            double semi_axis_a, semi_axis_b, semi_axis_c; //ellipsoid axes
            doubleType alpha, delta, W, PA1, PA2, PA3;
            math::Matrix<doubleType> R_from_ICRF_to_J2000BCI;
            math::Matrix<doubleType> R_from_J2000BCI_to_ICRF;
            math::Matrix<doubleType> R_from_J2000BCI_to_J2000BCF;
            math::Matrix<doubleType> R_from_J2000BCF_to_J2000BCI;
            math::Matrix<doubleType> R_from_ICRF_to_J2000BCF;
            math::Matrix<doubleType> R_from_J2000BCF_to_ICRF;
            math::Matrix<doubleType> R_from_BCI_to_ICRF;
            math::Matrix<doubleType> R_from_ICRF_to_BCI;
            math::Matrix<doubleType> R_from_BCI_to_BCF;
            math::Matrix<doubleType> R_from_BCF_to_BCI;
            math::Matrix<doubleType> R_from_BCF_to_ICRF;
            math::Matrix<doubleType> R_from_ICRF_to_BCF;
            math::Matrix<doubleType> dR_from_BCI_to_ICRF_dt;
            math::Matrix<doubleType> dR_from_ICRF_to_BCI_dt;
            math::Matrix<doubleType> dR_from_BCI_to_BCF_dt;
            math::Matrix<doubleType> dR_from_BCF_to_BCI_dt;
            math::Matrix<doubleType> dR_from_BCF_to_ICRF_dt;
            math::Matrix<doubleType> dR_from_ICRF_to_BCF_dt;
            math::Matrix<doubleType> dR_from_J2000BCI_to_J2000BCF_dt;
            math::Matrix<doubleType> dR_from_J2000BCF_to_J2000BCI_dt;
            math::Matrix<doubleType> dR_from_J2000BCF_to_ICRF_dt;
            math::Matrix<doubleType> dR_from_ICRF_to_J2000BCF_dt;

            math::Matrix<doubleType> R_from_PrincipleAxes_to_BCF;
            math::Matrix<doubleType> R_from_BCF_to_PrincipleAxes;
            math::Matrix<doubleType> dR_from_PrincipleAxes_to_BCF_dt;
            math::Matrix<doubleType> dR_from_BCF_to_PrincipleAxes_dt;

            math::Matrix<doubleType> R_from_Topocentric_to_BCF;
            math::Matrix<doubleType> R_from_BCF_to_Topocentric;
            math::Matrix<doubleType> dR_from_Topocentric_to_BCF_dt;
            math::Matrix<doubleType> dR_from_BCF_to_Topocentric_dt;
            math::Tensor<doubleType> dR_from_Topocentric_to_BCF_drBCF;
            math::Tensor<doubleType> dR_from_BCF_to_Topocentric_drBCF;
            math::Tensor<doubleType> dR_from_Topocentric_to_BCF_dvBCF;
            math::Tensor<doubleType> dR_from_BCF_to_Topocentric_dvBCF;

            math::Matrix<doubleType> R_from_Polar_to_BCF;
            math::Matrix<doubleType> R_from_BCF_to_Polar;
            math::Matrix<doubleType> dR_from_Polar_to_BCF_dt;
            math::Matrix<doubleType> dR_from_BCF_to_Polar_dt;

            math::Matrix<doubleType> R_from_SAM_to_J2000BCI;
            math::Matrix<doubleType> R_from_J2000BCI_to_SAM;
            math::Matrix<doubleType> dR_from_SAM_to_J2000BCI_dt;
            math::Matrix<doubleType> dR_from_J2000BCI_to_SAM_dt;

            math::Matrix<doubleType> R_from_ObjectReferenced_to_ICRF;
            math::Matrix<doubleType> R_from_ICRF_to_ObjectReferenced;
            math::Matrix<doubleType> dR_from_ObjectReferenced_to_ICRF_dt;
            math::Matrix<doubleType> dR_from_ICRF_to_ObjectReferenced_dt;

            math::Tensor<doubleType> dR_dreferenceVector;

            math::Matrix<doubleType> xhat, yhat, zhat;
        };

    }//close namespace Astrodynamics
} //close namespace EMTG