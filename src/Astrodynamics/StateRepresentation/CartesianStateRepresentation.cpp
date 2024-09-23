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

#include "CartesianStateRepresentation.h"

namespace EMTG
{
    namespace Astrodynamics
    {
        CartesianStateRepresentation::CartesianStateRepresentation(const double& mu) : StateRepresentationBase::StateRepresentationBase(mu)
        {
            this->name = "Cartesian";
            this->stateNames = std::vector<std::string>({ "x", "y", "z", "vx", "vy", "vz" });
        };

        math::Matrix<doubleType> CartesianStateRepresentation::convertFromRepresentationToCartesian(const math::Matrix<doubleType>& StateVectorThisRepresentationIn, const bool& needG)
        {
            //Step 1: set the state vector
            this->setStateVectorThisRepresentation(StateVectorThisRepresentationIn);

            //Step 2: convert the state           
            //this is trivial for Cartesian
            this->setStateVectorCartesian(this->StateVectorThisRepresentation);

            //Step 3: construct the partial derivatives
            //nothing to do here because the cartesian transition matrix is an identity

            //Step 4: return the converted vector
            return this->StateVectorCartesian;
        }

        math::Matrix<doubleType> CartesianStateRepresentation::convertFromCartesianToRepresentation(const math::Matrix<doubleType>& StateVectorCartesianIn, const bool& needG)
        {
            //Step 1: set the state vector
            this->setStateVectorCartesian(StateVectorCartesianIn);

            //Step 2: convert the state           
            //this is trivial for Cartesian
            this->setStateVectorThisRepresentation(this->StateVectorCartesian);

            //Step 3: construct the partial derivatives
            //nothing to do here because the cartesian transition matrix is an identity

            //Step 4: return the converted vector
            return this->StateVectorThisRepresentation;
        }
    }//end namespace StateRepresentation
}//end namespace EMTG*