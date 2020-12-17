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

// aerodynamic drag term

#ifndef AERODYNAMIC_DRAG_TERM_H
#define AERODYNAMIC_DRAG_TERM_H

#include "SpacecraftAccelerationModel.h"
#include "SpacecraftAccelerationModelTerm.h"
#include "doubleType.h"
#include "ExponentialAtmosphere.h"
#include "BodydeticConversions.h"

namespace EMTG
{
	namespace Astrodynamics
	{
		class AerodynamicDragTerm : public SpacecraftAccelerationModelTerm
		{
		public:
			//constructor
			AerodynamicDragTerm(SpacecraftAccelerationModel * acceleration_model_in);

			virtual ~AerodynamicDragTerm();

			//methods
			virtual void computeAccelerationTerm();
			virtual void computeAccelerationTerm(const bool & generate_derivatives);
			void populateInstrumentationFile(std::ofstream & acceleration_model_file) override;
			inline AerodynamicDragTerm* clone() const { return new AerodynamicDragTerm(*this); }

			//fields
			doubleType Cd; // coefficient of drag [0,inf)
			doubleType sc_area; // spacecraft area normal to wind direction, km^2
			doubleType density; // atmospheric density, kg/m^3
			doubleType DdensityDh; // derivative of density w.r.t. altitude, [kg/m^3] / km
			doubleType vRelativeMag; // norm of vRelative, km/s
			doubleType drag_force_norm; // kN
			math::Matrix<doubleType> omega; // angular velocity of central body about its spin axis in icrf, rad/s
			math::Matrix<doubleType> omega_BCI; // angular velocity of central body about its spin axis in true of date BCI, rad/s
			math::Matrix<doubleType> vAtmosphere;  // velocity vector of atmosphere in icrf, km/s
			math::Matrix<doubleType> vRelative; // v_spacecraft - vAtmosphere
			math::Matrix<doubleType> rBCF; // position vector in true of date body fixed frame (km)
			math::Matrix<doubleType> dLLAdBCF; // partials of lat/lon/alt w.r.t. true of date body fixed position vector
			math::Matrix <doubleType> lla_output; // latitude/longitude/altitude

			math::Matrix<double> DhDrBCF; // derivative of altitude w.r.t. true-of-date body-fixed
			math::Matrix<double> DrBCFDr; // derivative of position of true-of-date body-fixed position w.r.t. ICRF position; same as transformation matrix from ICRF to true-of-date, body-fixed frame
			math::Matrix<double> DhDr; // derivative of altitude w.r.t. ICRF r
			math::Matrix<double> Ddensity_kmDr; // derivative of density w.r.t. ICRF r
			math::Matrix<double> Ddensity_kmDt; // derivative of density w.r.t. ICRF t
			math::Matrix<double> DaccDstate; // derivative of acceleration vector w.r.t. [r, v, m, t]
			math::Matrix<double> DmInvDstate; // derivative of 1/mass w.r.t. [r, v, m, t]
			math::Matrix<double> Ddensity_kmDstate; // derivative of density w.r.t. [r, v, m, t]
			math::Matrix<double> DvRelativeMagDstate; // derivative of norm(vRelative) w.r.t. [r, v, m, t]
			math::Matrix<double> DvRelativeDstate; // derivative of vRelative (ICRF) w.r.t. [r, v, m, t]
			math::Matrix<double> DrDstate; // derivative of r (ICRF) w.r.t. [r, v, m, t]
			math::Matrix<double> DvDstate; // derivative of v (ICRF) w.r.t. [r, v, m, t]
			math::Matrix<double> DcrossTermDx; // derivative of omega \cross r w.r.t. [r, v, m, t]
			math::Matrix<double> DrBCFDt; // derivative of r TOD BCF w.r.t t
			math::Matrix<double> vRelativeDouble; // double version of vRelative
			math::Matrix<double> DaccDpropVars; // derivative of acceleration w.r.t. prop vars
			math::Matrix<doubleType> DomegaDstate; // derivative of omega w.r.t. [r, v, m, t]
			math::Matrix<doubleType> DrBCFDt_doubleType; // derivative of r TOD BCF w.r.t t
			math::Matrix<doubleType> dICRFdBCI; // rotation matrix from BCI to ICRF
			math::Matrix<doubleType> referenceVector; // don't care
			math::Matrix<doubleType> dreferenceVector_dt; // don't care
			math::Matrix<doubleType> DomegaDt; // derivative of omega in ICRF w.r.t. t
			math::Matrix<doubleType> dRstate_dreferenceVector; // don't care
			math::Matrix<doubleType> DrBCFDr_doubleType; // derivative of r in the BCF TOD frame w.r.t. r ICRF as a doubleType
			math::Matrix<doubleType> temp_matrix_out;
		};
	} // end Astrodynamics namespace
} // end EMTG namespace

#endif