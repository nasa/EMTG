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

#ifndef OUTGOING_B_PLANE_STATE_REPRESENTATION_H
#define OUTGOING_B_PLANE_STATE_REPRESENTATION_H

#include "StateRepresentation.h"

namespace EMTG
{
    namespace Astrodynamics
    {
        class OutgoingBplaneStateRepresentation : public StateRepresentationBase
        {
        public:
            OutgoingBplaneStateRepresentation() : OutgoingBplaneStateRepresentation(1.0) {};
            OutgoingBplaneStateRepresentation(const double& mu);

            ~OutgoingBplaneStateRepresentation() {};

            //methods
            virtual math::Matrix<doubleType> convertFromRepresentationToCartesian(const math::Matrix<doubleType>& StateVectorThisRepresentationIn, const bool& needG = false);
            virtual math::Matrix<doubleType> convertFromCartesianToRepresentation(const math::Matrix<doubleType>& StateVectorCartesianIn, const bool& needG = false);

            doubleType getBdotR() const { return this->BdotR; }
            doubleType getBdotT() const { return this->BdotT; }
            math::Matrix<doubleType> get_dBdotR_dx_cartesian() const { return this->dBdotR_dx_cartesian; }
            math::Matrix<doubleType> get_dBdotT_dx_cartesian() const { return this->dBdotT_dx_cartesian; }

        private:
            doubleType BdotR, BdotT;
            math::Matrix<doubleType> dBdotR_dx_cartesian;
            math::Matrix<doubleType> dBdotT_dx_cartesian;
        };
    }//end namespace StateRepresentation
}//end namespace EMTG

#endif // OUTGOING_B_PLANE_STATE_REPRESENTATION_H