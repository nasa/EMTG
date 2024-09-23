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

#include "SphericalRADECStateRepresentation.h"
#include "EMTG_math.h"

namespace EMTG
{
    namespace Astrodynamics
    {
        SphericalRADECStateRepresentation::SphericalRADECStateRepresentation(const double& mu) : StateRepresentationBase::StateRepresentationBase(mu)
        {
            this->name = "SphericalRADEC";
            this->stateNames = std::vector<std::string>({ "r", "RA", "DEC", "v", "vRA", "DEC" });
        };

        math::Matrix<doubleType> SphericalRADECStateRepresentation::convertFromRepresentationToCartesian(const math::Matrix<doubleType>& StateVectorThisRepresentationIn, const bool& needG)
        {
            //Step 1: set the state vector
            this->setStateVectorThisRepresentation(StateVectorThisRepresentationIn);

            //Step 2: convert the state           
            const doubleType& r =    this->StateVectorThisRepresentation(0);
            const doubleType& RA =   this->StateVectorThisRepresentation(1);
            const doubleType& DEC =  this->StateVectorThisRepresentation(2);
            const doubleType& v =    this->StateVectorThisRepresentation(3);
            const doubleType& vRA =  this->StateVectorThisRepresentation(4);
            const doubleType& vDEC = this->StateVectorThisRepresentation(5);

            doubleType cosRA = cos(RA);
            doubleType sinRA = sin(RA);
            doubleType cosDEC = cos(DEC);
            doubleType sinDEC = sin(DEC);
            doubleType cosvRA = cos(vRA);
            doubleType sinvRA = sin(vRA);
            doubleType cosvDEC = cos(vDEC);
            doubleType sinvDEC = sin(vDEC);

            this->StateVectorCartesian(0) = r * cosRA * cosDEC;
            this->StateVectorCartesian(1) = r * sinRA * cosDEC;
            this->StateVectorCartesian(2) = r * sinDEC;
            this->StateVectorCartesian(3) = v * cosvRA * cosvDEC;
            this->StateVectorCartesian(4) = v * sinvRA * cosvDEC;
            this->StateVectorCartesian(5) = v * sinvDEC;
            

            //Step 3: construct the partial derivatives
            if (needG)
            {
                this->RepresentationToCartesianTransitionMatrix.assign_zeros();
                
                doubleType dx_dr = (cosRA * cosDEC);
                doubleType dy_dr = (sinRA * cosDEC);
                doubleType dz_dr = sinDEC;

                doubleType dx_dRA = (-r * cosDEC * sinRA);
                doubleType dy_dRA = (r * cosDEC * cosRA);

                doubleType dx_dDEC = (-r * cosRA*sinDEC);
                doubleType dy_dDEC = (-r * sinDEC*sinRA);
                doubleType dz_dDEC = (r*cosDEC);

                doubleType dvx_dv = (cosvRA * cosvDEC);
                doubleType dvy_dv = (sinvRA * cosvDEC);
                doubleType dvz_dv = sinvDEC;

                doubleType dvx_dvRA = (-v * cosvDEC * sinvRA);
                doubleType dvy_dvRA = (v * cosvDEC * cosvRA);

                doubleType dvx_dvDEC = (-v * cosvRA*sinvDEC);
                doubleType dvy_dvDEC = (-v * sinvDEC*sinvRA);
                doubleType dvz_dvDEC = (v*cosvDEC);

                //transition matrix has NewThing for rows and OldThing for columns
                this->RepresentationToCartesianTransitionMatrix(0, 0) = dx_dr;
                this->RepresentationToCartesianTransitionMatrix(1, 0) = dy_dr;
                this->RepresentationToCartesianTransitionMatrix(2, 0) = dz_dr;

                this->RepresentationToCartesianTransitionMatrix(0, 1) = dx_dRA;
                this->RepresentationToCartesianTransitionMatrix(1, 1) = dy_dRA;

                this->RepresentationToCartesianTransitionMatrix(0, 2) = dx_dDEC;
                this->RepresentationToCartesianTransitionMatrix(1, 2) = dy_dDEC;
                this->RepresentationToCartesianTransitionMatrix(2, 2) = dz_dDEC;

                this->RepresentationToCartesianTransitionMatrix(3, 3) = dvx_dv;
                this->RepresentationToCartesianTransitionMatrix(4, 3) = dvy_dv;
                this->RepresentationToCartesianTransitionMatrix(5, 3) = dvz_dv;

                this->RepresentationToCartesianTransitionMatrix(3, 4) = dvx_dvRA;
                this->RepresentationToCartesianTransitionMatrix(4, 4) = dvy_dvRA;

                this->RepresentationToCartesianTransitionMatrix(3, 5) = dvx_dvDEC;
                this->RepresentationToCartesianTransitionMatrix(4, 5) = dvy_dvDEC;
                this->RepresentationToCartesianTransitionMatrix(5, 5) = dvz_dvDEC;
            }

            //Step 4: return the converted vector
            return this->StateVectorCartesian;
        }//end convertFromRepresentationToCartesian()

        math::Matrix<doubleType> SphericalRADECStateRepresentation::convertFromCartesianToRepresentation(const math::Matrix<doubleType>& StateVectorCartesianIn, const bool& needG)
        {
            //Step 1: set the state vector
            this->setStateVectorCartesian(StateVectorCartesianIn);

            //Step 2: convert the state           
            //this is trivial for Cartesian
            const doubleType& x = this->StateVectorCartesian(0);
            const doubleType& y = this->StateVectorCartesian(1);
            const doubleType& z = this->StateVectorCartesian(2);
            const doubleType& vx = this->StateVectorCartesian(3);
            const doubleType& vy = this->StateVectorCartesian(4);
            const doubleType& vz = this->StateVectorCartesian(5);

            doubleType r = sqrt(x * x + y*y + z*z);
            doubleType RA = atan2(y, x);
            doubleType DEC = math::safe_asin(z / r);
            doubleType v = sqrt(vx * vx + vy*vy + vz*vz);
            doubleType vRA = atan2(vy, vx);
            doubleType vDEC = math::safe_asin(vz / v);

            this->StateVectorThisRepresentation(0) = r;
            this->StateVectorThisRepresentation(1) = RA;
            this->StateVectorThisRepresentation(2) = DEC;
            this->StateVectorThisRepresentation(3) = v;
            this->StateVectorThisRepresentation(4) = vRA;
            this->StateVectorThisRepresentation(5) = vDEC;

            //Step 3: construct the partial derivatives
            if (needG)
            {
                this->CartesianToRepresentationTransitionMatrix.assign_zeros();

                doubleType r3 = r * r * r;
                doubleType v3 = v * v * v;

                doubleType dr_dx = x / r;
                doubleType dr_dy = y / r;
                doubleType dr_dz = z / r;

                doubleType dRA_dx = -y / (x * x + y * y);
                doubleType dRA_dy = x / (x * x + y * y);

                doubleType dDEC_dx = -1.0 / sqrt(1.0 - (z*z) / (r*r)) * x * z / r3;
                doubleType dDEC_dy = -1.0 / sqrt(1.0 - (z*z) / (r*r)) * y * z / r3;
                doubleType dDEC_dz = 1.0 / sqrt(1.0 - (z*z) / (r*r)) * (x * x + y * y) / r3;

                doubleType dv_dvx = vx / v;
                doubleType dv_dvy = vy / v;
                doubleType dv_dvz = vz / v;

                doubleType dvRA_dvx = -vy / (vx * vx + vy * vy);
                doubleType dvRA_dvy = vx / (vx * vx + vy * vy);

                doubleType dvDEC_dvx = -1.0 / sqrt(1.0 - (vz*vz) / (v*v)) * vx * vz / v3;
                doubleType dvDEC_dvy = -1.0 / sqrt(1.0 - (vz*vz) / (v*v)) * vy * vz / v3;
                doubleType dvDEC_dvz = 1.0 / sqrt(1.0 - (vz*vz) / (v*v)) * (vx * vx + vy * vy) / v3;

                //transition matrix has NewThing for rows and OldThing for columns
                this->CartesianToRepresentationTransitionMatrix(0, 0) = dr_dx;
                this->CartesianToRepresentationTransitionMatrix(0, 1) = dr_dy;
                this->CartesianToRepresentationTransitionMatrix(0, 2) = dr_dz;

                this->CartesianToRepresentationTransitionMatrix(1, 0) = dRA_dx;
                this->CartesianToRepresentationTransitionMatrix(1, 1) = dRA_dy;

                this->CartesianToRepresentationTransitionMatrix(2, 0) = dDEC_dx;
                this->CartesianToRepresentationTransitionMatrix(2, 1) = dDEC_dy;
                this->CartesianToRepresentationTransitionMatrix(2, 2) = dDEC_dz;

                this->CartesianToRepresentationTransitionMatrix(3, 3) = dv_dvx;
                this->CartesianToRepresentationTransitionMatrix(3, 4) = dv_dvy;
                this->CartesianToRepresentationTransitionMatrix(3, 5) = dv_dvz;

                this->CartesianToRepresentationTransitionMatrix(4, 3) = dvRA_dvx;
                this->CartesianToRepresentationTransitionMatrix(4, 4) = dvRA_dvy;

                this->CartesianToRepresentationTransitionMatrix(5, 3) = dvDEC_dvx;
                this->CartesianToRepresentationTransitionMatrix(5, 4) = dvDEC_dvy;
                this->CartesianToRepresentationTransitionMatrix(5, 5) = dvDEC_dvz;
            }

            //Step 4: return the converted vector
            return this->StateVectorThisRepresentation;
        }//end convertFromCartesianToRepresentation()
    }//end namespace StateRepresentation
}//end namespace EMTG