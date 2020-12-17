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

#include "SphericalAZFPAStateRepresentation.h"
#include "EMTG_math.h"

namespace EMTG
{
    namespace Astrodynamics
    {
        SphericalAZFPAStateRepresentation::SphericalAZFPAStateRepresentation(const double& mu) : StateRepresentationBase::StateRepresentationBase(mu)
        {
            this->name = "SphericalAZFPA";
            this->stateNames = std::vector<std::string>({ "r", "RA", "DEC", "v", "AZ", "FPA" });
        };

        math::Matrix<doubleType> SphericalAZFPAStateRepresentation::convertFromRepresentationToCartesian(const math::Matrix<doubleType>& StateVectorThisRepresentationIn, const bool& needG)
        {
            //Step 1: set the state vector
            this->setStateVectorThisRepresentation(StateVectorThisRepresentationIn);

            //Step 2: convert the state           
            const doubleType& r = this->StateVectorThisRepresentation(0);
            const doubleType& RA = this->StateVectorThisRepresentation(1);
            const doubleType& DEC = this->StateVectorThisRepresentation(2);
            const doubleType& v = this->StateVectorThisRepresentation(3);
            const doubleType& AZ = this->StateVectorThisRepresentation(4);
            const doubleType& FPA = this->StateVectorThisRepresentation(5);

            doubleType cosRA = cos(RA);
            doubleType sinRA = sin(RA);
            doubleType cosDEC = cos(DEC);
            doubleType sinDEC = sin(DEC);
            doubleType cosAZ = cos(AZ);
            doubleType sinAZ = sin(AZ);
            doubleType cosFPA = cos(FPA);
            doubleType sinFPA = sin(FPA);

            this->StateVectorCartesian(0) = r * cosRA * cosDEC;
            this->StateVectorCartesian(1) = r * sinRA * cosDEC;
            this->StateVectorCartesian(2) = r * sinDEC;
            this->StateVectorCartesian(3) = -v * (sinFPA*(sinAZ*sinRA + cosAZ * cosRA*sinDEC) - cosFPA * cosDEC*cosRA);
            this->StateVectorCartesian(4) =  v * (sinFPA*(cosRA*sinAZ - cosAZ * sinDEC*sinRA) + cosFPA * cosDEC*sinRA);
            this->StateVectorCartesian(5) =  v * (cosFPA*sinDEC + cosDEC * cosAZ*sinFPA);
                       
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

                doubleType dvx_dRA = (-v * (sinFPA*(cosRA*sinAZ - cosAZ * sinDEC*sinRA) + cosFPA * cosDEC*sinRA));
                doubleType dvy_dRA = (-v * (sinFPA*(sinAZ*sinRA + cosAZ * cosRA*sinDEC) - cosFPA * cosDEC*cosRA));

                doubleType dvx_dDEC = (-v * cosRA*(cosFPA*sinDEC + cosDEC * cosAZ*sinFPA));
                doubleType dvy_dDEC = (-v * sinRA*(cosFPA*sinDEC + cosDEC * cosAZ*sinFPA));
                doubleType dvz_dDEC = (v*(cosFPA*cosDEC - cosAZ * sinFPA*sinDEC));

                doubleType dvx_dv = (-(sinFPA*(sinAZ*sinRA + cosAZ * cosRA*sinDEC) - cosFPA * cosDEC*cosRA));
                doubleType dvy_dv = ((sinFPA*(cosRA*sinAZ - cosAZ * sinDEC*sinRA) + cosFPA * cosDEC*sinRA));
                doubleType dvz_dv = ((cosFPA*sinDEC + cosDEC * cosAZ*sinFPA));

                doubleType dvx_dFPA = (-v * (cosFPA*(sinAZ*sinRA + cosAZ * cosRA*sinDEC) + cosDEC * sinFPA*cosRA));
                doubleType dvy_dFPA = (v*(cosFPA*(cosRA*sinAZ - cosAZ * sinDEC*sinRA) - cosDEC * sinFPA*sinRA));
                doubleType dvz_dFPA = (v*cosFPA*cosDEC*cosAZ - v * sinFPA*sinDEC);

                doubleType dvx_dAZ = (-v * sinFPA*(cosAZ*sinRA - cosRA * sinDEC*sinAZ));
                doubleType dvy_dAZ = (v*sinFPA*(cosAZ*cosRA + sinDEC * sinAZ*sinRA));
                doubleType dvz_dAZ = (-v * cosDEC*sinFPA*sinAZ);

                //transition matrix has NewThing for rows and OldThing for columns
                this->RepresentationToCartesianTransitionMatrix(0, 0) = dx_dr;
                this->RepresentationToCartesianTransitionMatrix(1, 0) = dy_dr;
                this->RepresentationToCartesianTransitionMatrix(2, 0) = dz_dr;

                this->RepresentationToCartesianTransitionMatrix(0, 1) = dx_dRA;
                this->RepresentationToCartesianTransitionMatrix(1, 1) = dy_dRA;
                this->RepresentationToCartesianTransitionMatrix(3, 1) = dvx_dRA;
                this->RepresentationToCartesianTransitionMatrix(4, 1) = dvy_dRA;

                this->RepresentationToCartesianTransitionMatrix(0, 2) = dx_dDEC;
                this->RepresentationToCartesianTransitionMatrix(1, 2) = dy_dDEC;
                this->RepresentationToCartesianTransitionMatrix(2, 2) = dz_dDEC;
                this->RepresentationToCartesianTransitionMatrix(3, 2) = dvx_dDEC;
                this->RepresentationToCartesianTransitionMatrix(4, 2) = dvy_dDEC;
                this->RepresentationToCartesianTransitionMatrix(5, 2) = dvz_dDEC;

                this->RepresentationToCartesianTransitionMatrix(3, 3) = dvx_dv;
                this->RepresentationToCartesianTransitionMatrix(4, 3) = dvy_dv;
                this->RepresentationToCartesianTransitionMatrix(5, 3) = dvz_dv;

                this->RepresentationToCartesianTransitionMatrix(3, 4) = dvx_dAZ;
                this->RepresentationToCartesianTransitionMatrix(4, 4) = dvy_dAZ;
                this->RepresentationToCartesianTransitionMatrix(5, 4) = dvz_dAZ;

                this->RepresentationToCartesianTransitionMatrix(3, 5) = dvx_dFPA;
                this->RepresentationToCartesianTransitionMatrix(4, 5) = dvy_dFPA;
                this->RepresentationToCartesianTransitionMatrix(5, 5) = dvz_dFPA;
            }

            //Step 4: return the converted vector
            return this->StateVectorCartesian;
        }//end convertFromRepresentationToCartesian()

        math::Matrix<doubleType> SphericalAZFPAStateRepresentation::convertFromCartesianToRepresentation(const math::Matrix<doubleType>& StateVectorCartesianIn, const bool& needG)
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

            doubleType r = sqrt(x * x + y * y + z * z);
            doubleType RA = atan2(y, x);
            doubleType DEC = math::safe_asin(z / r);
            doubleType v = sqrt(vx * vx + vy * vy + vz * vz);


            doubleType FPA = acos((x * vx + y * vy + z * vz) / r / v);

            //azimuth is complicated
            math::Matrix<doubleType> xhat(3, 1, { cos(RA)*cos(DEC), sin(RA)*cos(DEC), sin(DEC) });
            math::Matrix<doubleType> yhat(3, 1, { cos(RA + math::PIover2), sin(RA + math::PIover2), 0.0 });
            math::Matrix<doubleType> zhat(3, 1, { -cos(RA)*sin(DEC), -sin(RA)*sin(DEC), cos(DEC) });
            math::Matrix<doubleType> R = (xhat.horz_cat(yhat).horz_cat(zhat)).transpose();
            math::Matrix<doubleType> V(3, 1, { vx, vy, vz });
            math::Matrix<doubleType> Vprime = R * V;
            doubleType AZ = atan2(Vprime(1), Vprime(2));

            this->StateVectorThisRepresentation(0) = r;
            this->StateVectorThisRepresentation(1) = RA;
            this->StateVectorThisRepresentation(2) = DEC;
            this->StateVectorThisRepresentation(3) = v;
            this->StateVectorThisRepresentation(4) = AZ;
            this->StateVectorThisRepresentation(5) = FPA;

            //Step 3: construct the partial derivatives
            if (needG)
            {
                this->CartesianToRepresentationTransitionMatrix.assign_zeros();

                doubleType r3 = r * r * r;
                doubleType v3 = v * v * v;

                doubleType dRA_dx = -y / (x * x + y * y);
                doubleType dRA_dy = x / (x * x + y * y);

                doubleType dDEC_dx = -1.0 / sqrt(1.0 - (z*z) / (r*r)) * x * z / r3;
                doubleType dDEC_dy = -1.0 / sqrt(1.0 - (z*z) / (r*r)) * y * z / r3;
                doubleType dDEC_dz = 1.0 / sqrt(1.0 - (z*z) / (r*r)) * (x * x + y * y) / r3;
                
                doubleType cosRA = cos(RA);
                doubleType sinRA = sin(RA);
                doubleType cosDEC = cos(DEC);
                doubleType sinDEC = sin(DEC);

                doubleType dvpx_drx = (-sinDEC * cosRA * dDEC_dx - cosDEC * sinRA * dRA_dx) * vx + (-sinDEC * sinRA * dDEC_dx + cosDEC * cosRA * dRA_dx) * vy + cosDEC * dDEC_dx * vz;
                doubleType dvpx_dry = (-sinDEC * cosRA * dDEC_dy - cosDEC * sinRA * dRA_dy) * vx + (-sinDEC * sinRA * dDEC_dy + cosDEC * cosRA * dRA_dy) * vy + cosDEC * dDEC_dy * vz;
                doubleType dvpx_drz = (-sinDEC * cosRA * dDEC_dz)                           * vx + (-sinDEC * sinRA * dDEC_dz)                           * vy + cosDEC * dDEC_dz * vz;
                doubleType dvpy_drx = -cosRA * dRA_dx * vx - sinRA * dRA_dx * vy;
                doubleType dvpy_dry = -cosRA * dRA_dy * vx - sinRA * dRA_dy * vy;
                doubleType dvpy_drz = 0.0;
                doubleType dvpz_drx = (-cosDEC * cosRA * dDEC_dx + sinDEC * sinRA * dRA_dx) * vx + (-cosDEC * sinRA * dDEC_dx - sinDEC * cosRA * dRA_dx) * vy - sinDEC * dDEC_dx * vz;
                doubleType dvpz_dry = (-cosDEC * cosRA * dDEC_dy + sinDEC * sinRA * dRA_dy) * vx + (-cosDEC * sinRA * dDEC_dy - sinDEC * cosRA * dRA_dy) * vy - sinDEC * dDEC_dy * vz;
                doubleType dvpz_drz = (-cosDEC * cosRA * dDEC_dz)                           * vx + (-cosDEC * sinRA * dDEC_dz)                           * vy - sinDEC * dDEC_dz * vz;

                doubleType dvpx_dvx = cosDEC * cosRA;
                doubleType dvpx_dvy = cosDEC * sinRA;
                doubleType dvpx_dvz = sinDEC;
                doubleType dvpy_dvx = -sinRA;
                doubleType dvpy_dvy = cosRA;
                doubleType dvpy_dvz = 0.0;
                doubleType dvpz_dvx = -sinDEC * cosRA;
                doubleType dvpz_dvy = -sinDEC * sinRA;
                doubleType dvpz_dvz = cosDEC;

                doubleType vpx = Vprime(0);
                doubleType vpy = Vprime(1);
                doubleType vpz = Vprime(2);

                doubleType dAZ_dx = 1.0 / (vpy * vpy + vpz * vpz) * (-vpy * dvpz_drx + vpz * dvpy_drx);
                doubleType dAZ_dy = 1.0 / (vpy * vpy + vpz * vpz) * (-vpy * dvpz_dry + vpz * dvpy_dry);
                doubleType dAZ_dz = 1.0 / (vpy * vpy + vpz * vpz) * (-vpy * dvpz_drz + vpz * dvpy_drz);
                doubleType dAZ_dvx = 1.0 / (vpy * vpy + vpz * vpz) * (-vpy * dvpz_dvx + vpz * dvpy_dvx);
                doubleType dAZ_dvy = 1.0 / (vpy * vpy + vpz * vpz) * (-vpy * dvpz_dvy + vpz * dvpy_dvy);
                doubleType dAZ_dvz = 1.0 / (vpy * vpy + vpz * vpz) * (-vpy * dvpz_dvz + vpz * dvpy_dvz);

                doubleType rdotv = x * vx + y * vy + z * vz;

                doubleType dr_dx = x / r;
                doubleType dr_dy = y / r;
                doubleType dr_dz = z / r;
                doubleType dv_dvx = vx / v;
                doubleType dv_dvy = vy / v;
                doubleType dv_dvz = vz / v;
                doubleType cosFPA = (x * vx + y * vy + z * vz) / r / v;

                doubleType dcosFPA_dx  = vx / (r * v) - rdotv * dr_dx / (r * r * v);
                doubleType dcosFPA_dy  = vy / (r * v) - rdotv * dr_dy / (r * r * v);
                doubleType dcosFPA_dz  = vz / (r * v) - rdotv * dr_dz / (r * r * v);
                doubleType dcosFPA_dvx =  x / (r * v) - rdotv * dv_dvx / (r * v * v);
                doubleType dcosFPA_dvy =  y / (r * v) - rdotv * dv_dvy / (r * v * v);
                doubleType dcosFPA_dvz =  z / (r * v) - rdotv * dv_dvz / (r * v * v);                

                doubleType dFPA_dx = -1.0 / sqrt( 1 - cosFPA * cosFPA) * dcosFPA_dx;
                doubleType dFPA_dy = -1.0 / sqrt( 1 - cosFPA * cosFPA) * dcosFPA_dy;
                doubleType dFPA_dz = -1.0 / sqrt( 1 - cosFPA * cosFPA) * dcosFPA_dz;
                doubleType dFPA_dvx = -1.0 / sqrt( 1 - cosFPA * cosFPA) * dcosFPA_dvx;
                doubleType dFPA_dvy = -1.0 / sqrt( 1 - cosFPA * cosFPA) * dcosFPA_dvy;
                doubleType dFPA_dvz = -1.0 / sqrt( 1 - cosFPA * cosFPA) * dcosFPA_dvz;

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

                this->CartesianToRepresentationTransitionMatrix(4, 0) = dAZ_dx;
                this->CartesianToRepresentationTransitionMatrix(4, 1) = dAZ_dy;
                this->CartesianToRepresentationTransitionMatrix(4, 2) = dAZ_dz;
                this->CartesianToRepresentationTransitionMatrix(4, 3) = dAZ_dvx;
                this->CartesianToRepresentationTransitionMatrix(4, 4) = dAZ_dvy;
                this->CartesianToRepresentationTransitionMatrix(4, 5) = dAZ_dvz;

                this->CartesianToRepresentationTransitionMatrix(5, 0) = dFPA_dx;
                this->CartesianToRepresentationTransitionMatrix(5, 1) = dFPA_dy;
                this->CartesianToRepresentationTransitionMatrix(5, 2) = dFPA_dz;
                this->CartesianToRepresentationTransitionMatrix(5, 3) = dFPA_dvx;
                this->CartesianToRepresentationTransitionMatrix(5, 4) = dFPA_dvy;
                this->CartesianToRepresentationTransitionMatrix(5, 5) = dFPA_dvz;
            }

            //Step 4: return the converted vector
            return this->StateVectorThisRepresentation;
        }//end convertFromCartesianToRepresentation()
    }//end namespace StateRepresentation
}//end namespace EMTG