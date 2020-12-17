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

#ifndef CARTESIAN_STATE_REPRESENTATION_H
#define CARTESIAN_STATE_REPRESENTATION_H

#include "StateRepresentation.h"

namespace EMTG
{
    namespace Astrodynamics
    {
        class CartesianStateRepresentation : public StateRepresentationBase
        {
        public:
            CartesianStateRepresentation(const double& mu);

            ~CartesianStateRepresentation() {};

            //methods
            virtual math::Matrix<doubleType> convertFromRepresentationToCartesian(const math::Matrix<doubleType>& StateVectorThisRepresentationIn, const bool& needG = false);
            virtual math::Matrix<doubleType> convertFromCartesianToRepresentation(const math::Matrix<doubleType>& StateVectorCartesianIn, const bool& needG = false);
        };
    }//end namespace StateRepresentation
}//end namespace EMTG

#endif // CARTESIAN_STATE_REPRESENTATION_H