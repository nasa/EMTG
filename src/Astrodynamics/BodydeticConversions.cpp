// EMTG: Evolutionary Mission Trajectory Generator
// An open-source global optimization tool for preliminary mission design
// Provided by NASA Goddard Space Flight Center
//
// Copyright (c) 2014 - 2016 United States Government as represented by the
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

// source file for bodydetic conversions
// Noble Hatten 11/19/2019

#include<vector>

#include "EMTG_math.h"
#include "EMTG_Matrix.h"
#include "BodydeticConversions.h"


namespace EMTG
{
    namespace Astrodynamics
    {

        /*
        Convert from latitude, longitude, altitude (LLA) to body-centered, body-fixed (BCF) Cartesian position vector.
        Uses equations that work if body is an oblate spheroid. Could work on triaxial ellipsoid with modification of
        input arguments and calculation of ee2.
        @param LLA 3x1; geodetic latitude (rad), longitude (rad), geodetic altitude (km)
        @param Re Equatorial radius of central body (km)
        @param f flattening coefficient of central body; f = (Re - Rp) / Re
        @param rbcf 3x1; Body-centered, body-fixed (BCF) Cartesian position vector (km)
        @param generateDerivatives Whether or not to generate derivatives; CURRENTLY DOES NOTHING, NOT IMPLENENTED UNLESS IT TURNS OUT WE NEED IT
        @param dBCFdLLA 3x3; If generating derivatives, the Jacobian is populated here
        */
        void LLA2BCF_oblate(const math::Matrix<doubleType>& LLA,
            const doubleType& Re,
            const doubleType& f,
            math::Matrix<doubleType>& rbcf,
            const bool& generateDerivatives,
            math::Matrix<doubleType>& dBCFdLLA)
        {
            doubleType lat = LLA(0);
            doubleType lon = LLA(1);
            doubleType alt = LLA(2);
            doubleType sLat = sin(lat);
            doubleType cLat = cos(lat);
            doubleType sLon = sin(lon);
            doubleType Rp = Re * (1. - f);
            doubleType ex2 = (pow(Re, 2) - pow(Rp, 2)) / pow(Re, 2);
            doubleType ee2 = 0.; // true for oblate spheroid, not true for triaxial ellipsoid
            doubleType v = Re / pow(1. - ex2 * pow(sLat, 2) - ee2 * pow(cLat, 2) * pow(sLon, 2), 0.5);
            rbcf(0) = (v + alt) * cLat * cos(lon);
            rbcf(1) = (v * (1. - ee2) + alt) * cLat * sLon;
            rbcf(2) = (v * (1. - ex2) + alt) * sLat;

        }
        /*
        Convert from body-centered, body-fixed (BCF) Cartesian position vector to latitude, longitude, altitude (LLA).
        Uses equations that work if body is an oblate spheroid.
        Based on Ryo Nakamura's code originally written in BoundaryDeticAltitudeConstraint.cpp
        @param rbcf 3x1; Body-centered, body-fixed (BCF) Cartesian position vector (km)
        @param Re Equatorial radius of central body (km)
        @param f flattening coefficient of central body; f = (Re - Rp) / Re
        @param LLA 3x1; geodetic latitude (rad), longitude (rad), geodetic altitude (km)
        @param generateDerivatives Whether or not to generate derivatives
        @param dLLAdBCF 3x3; If generating derivatives, the Jacobian is populated here
        */
        void BCF2LLA_oblate(const math::Matrix<doubleType>& rbcf,
            const doubleType& Re,
            const doubleType& f,
            math::Matrix<doubleType>& LLA,
            const bool& generateDerivatives,
            math::Matrix<doubleType>& dLLAdBCF)
        {
            //Step 3: compute deticAltitude magnitude
            //First, make sure that the probe is not at the pole because there is a singularity at the pole. 
            doubleType rx = rbcf(0);
            doubleType ry = rbcf(1);
            doubleType rz = rbcf(2);
            doubleType rxy = sqrt(rx * rx + ry * ry);
            if (rxy < math::SMALL)
            {
                throw std::runtime_error("The probe is too close to the pole to calculate detic altitude and its derivatives.");
            }

            doubleType f2 = (1.0 - f) * (1.0 - f);
            doubleType e2 = 1. - f2;

            doubleType r = rbcf.norm();

            doubleType lambda = atan2(rz, rxy);
            doubleType tanlambda = tan(lambda);
            doubleType tanlambda2 = tanlambda * tanlambda;
            doubleType xa = (1.0 - f) * Re / sqrt(tanlambda2 + f2);
            doubleType mua = atan2(tanlambda, f2);
            doubleType ra = xa / cos(lambda);
            doubleType l = r - ra;
            doubleType dellambda = mua - lambda;
            doubleType h = l * cos(dellambda);
            doubleType rhoa = Re * f2 / pow(1.0 - (2.0 * f - f * f) * pow(sin(mua), 2.0), 1.5);
            doubleType delmu1 = l * sin(dellambda) / (rhoa + h);
            doubleType delmu = atan(delmu1);
            doubleType deticLatitude = mua - delmu;

            doubleType mu = mua - delmu; // mu is deticLatitude

            doubleType clat = cos(mu);
            doubleType slat = sin(mu);

            doubleType N = Re / sqrt(1.0 - e2 * slat * slat);

            doubleType deticAltitude = rxy * clat + (rz + e2 * N * slat) * slat - N;

            doubleType longitude = atan2(ry, rx);

            LLA(0) = deticLatitude;
            LLA(1) = longitude;
            LLA(2) = deticAltitude;

            if (generateDerivatives)
            {
                doubleType dr_drx = rx / r;
                doubleType dr_dry = ry / r;
                doubleType dr_drz = rz / r;

                doubleType dlambda_drx = -(rx * pow(rxy, -3.0) * rz) / (rz * rz / rxy / rxy + 1.0);
                doubleType dlambda_dry = -(ry * pow(rxy, -3.0) * rz) / (rz * rz / rxy / rxy + 1.0);
                doubleType dlambda_drz = 1.0 / (rxy * (rz * rz / rxy / rxy + 1.0));

                doubleType dxa_dlambda = -(Re * (1.0 - f) * pow(1.0 / cos(lambda), 2.0) * tanlambda) / pow(tanlambda2 + f2, 1.5);
                doubleType dxa_drx = dxa_dlambda * dlambda_drx;
                doubleType dxa_dry = dxa_dlambda * dlambda_dry;
                doubleType dxa_drz = dxa_dlambda * dlambda_drz;

                doubleType dmua_dlambda = pow(1.0 / cos(lambda), 2.0) / (f2 * (tanlambda2 / f2 / f2 + 1.0));
                doubleType dmua_drx = dmua_dlambda * dlambda_drx;
                doubleType dmua_dry = dmua_dlambda * dlambda_dry;
                doubleType dmua_drz = dmua_dlambda * dlambda_drz;

                doubleType dra_drx = dxa_drx / cos(lambda) + xa * tanlambda / cos(lambda) * dlambda_drx;
                doubleType dra_dry = dxa_dry / cos(lambda) + xa * tanlambda / cos(lambda) * dlambda_dry;
                doubleType dra_drz = dxa_drz / cos(lambda) + xa * tanlambda / cos(lambda) * dlambda_drz;

                doubleType dl_drx = dr_drx - dra_drx;
                doubleType dl_dry = dr_dry - dra_dry;
                doubleType dl_drz = dr_drz - dra_drz;

                doubleType ddellambda_drx = dmua_drx - dlambda_drx;
                doubleType ddellambda_dry = dmua_dry - dlambda_dry;
                doubleType ddellambda_drz = dmua_drz - dlambda_drz;

                doubleType dh_drx = dl_drx * cos(dellambda) - l * sin(dellambda) * ddellambda_drx;
                doubleType dh_dry = dl_dry * cos(dellambda) - l * sin(dellambda) * ddellambda_dry;
                doubleType dh_drz = dl_drz * cos(dellambda) - l * sin(dellambda) * ddellambda_drz;

                doubleType drhoa_dmua = (3.0 * Re * f2 * (2.0 * f - f * f) * cos(mua) * sin(mua)) / pow(1.0 - (2.0 * f - f * f) * pow(sin(mua), 2.0), 2.5);
                doubleType drhoa_drx = drhoa_dmua * dmua_drx;
                doubleType drhoa_dry = drhoa_dmua * dmua_dry;
                doubleType drhoa_drz = drhoa_dmua * dmua_drz;

                doubleType ddelmu1_drx = l * cos(dellambda) / (rhoa + h) * ddellambda_drx + delmu1 * (dl_drx / l - (drhoa_drx + dh_drx) / (rhoa + h));
                doubleType ddelmu1_dry = l * cos(dellambda) / (rhoa + h) * ddellambda_dry + delmu1 * (dl_dry / l - (drhoa_dry + dh_dry) / (rhoa + h));
                doubleType ddelmu1_drz = l * cos(dellambda) / (rhoa + h) * ddellambda_drz + delmu1 * (dl_drz / l - (drhoa_drz + dh_drz) / (rhoa + h));
                doubleType ddelmu_drx = 1.0 / (delmu1 * delmu1 + 1) * ddelmu1_drx;
                doubleType ddelmu_dry = 1.0 / (delmu1 * delmu1 + 1) * ddelmu1_dry;
                doubleType ddelmu_drz = 1.0 / (delmu1 * delmu1 + 1) * ddelmu1_drz;

                doubleType dmu_drx = (dmua_drx - ddelmu_drx);
                doubleType dmu_dry = (dmua_dry - ddelmu_dry);
                doubleType dmu_drz = (dmua_drz - ddelmu_drz);

                doubleType dN_drx = Re * e2 * slat / pow(1.0 - e2 * slat * slat, 1.5) * clat * dmu_drx;
                doubleType dN_dry = Re * e2 * slat / pow(1.0 - e2 * slat * slat, 1.5) * clat * dmu_dry;
                doubleType dN_drz = Re * e2 * slat / pow(1.0 - e2 * slat * slat, 1.5) * clat * dmu_drz;

                doubleType ddetalt_drx = (rx / rxy * clat - rxy * slat * dmu_drx + (e2 * dN_drx * slat + e2 * N * clat * dmu_drx) * slat + (rz + e2 * N * slat) * clat * dmu_drx - dN_drx);
                doubleType ddetalt_dry = (ry / rxy * clat - rxy * slat * dmu_dry + (e2 * dN_dry * slat + e2 * N * clat * dmu_dry) * slat + (rz + e2 * N * slat) * clat * dmu_dry - dN_dry);
                doubleType ddetalt_drz = (-rxy * slat * dmu_drz + (1.0 + e2 * dN_drz * slat + e2 * N * clat * dmu_drz) * slat + (rz + e2 * N * slat) * clat * dmu_drz - dN_drz);

                doubleType dLonTerm = 1. / pow(rxy, 2);
                doubleType dLondX = -ry * dLonTerm;
                doubleType dLondY = rx * dLonTerm;

                // populate the matrix

                // latitude
                dLLAdBCF(0, 0) = dmu_drx;
                dLLAdBCF(0, 1) = dmu_dry;
                dLLAdBCF(0, 2) = dmu_drz;

                // longitude
                dLLAdBCF(1, 0) = dLondX;
                dLLAdBCF(1, 1) = dLondY;
                dLLAdBCF(1, 2) = 0.;

                // altitude
                dLLAdBCF(2, 0) = ddetalt_drx;
                dLLAdBCF(2, 1) = ddetalt_dry;
                dLLAdBCF(2, 2) = ddetalt_drz;
            }
        }
    }//end namespace Astrodynamics
}//end namespace EMTG