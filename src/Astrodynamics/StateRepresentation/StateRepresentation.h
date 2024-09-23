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

#ifndef STATE_REPRESENTATION_H
#define STATE_REPRESENTATION_H

#include <vector>

#include "doubleType.h"
#include <EMTG_Matrix.h>

namespace EMTG
{
    namespace Astrodynamics
    {
        class StateRepresentationBase
        {
        public:
            StateRepresentationBase(const double& mu) : StateVectorThisRepresentation(6, 1, 0.0),
                StateVectorCartesian(6, 1, 0.0),
                RepresentationToCartesianTransitionMatrix(6, math::identity),
                CartesianToRepresentationTransitionMatrix(6, math::identity),
                mu(mu)
            {};

            virtual ~StateRepresentationBase() {};

            //methods
            virtual math::Matrix<doubleType> convertFromRepresentationToCartesian(const math::Matrix<doubleType>& StateVectorThisRepresentationIn, const bool& needG = false) = 0;
            virtual math::Matrix<doubleType> convertFromCartesianToRepresentation(const math::Matrix<doubleType>& StateVectorCartesianIn, const bool& needG = false) = 0;

            inline math::Matrix<doubleType> getStateVectorThisRepresentation() const { return this->StateVectorThisRepresentation; }
            inline math::Matrix<doubleType> getStateVectorCartesian() const { return this->StateVectorCartesian; }
            inline math::Matrix<doubleType> getRepresentationToCartesianTransitionMatrix() const { return this->RepresentationToCartesianTransitionMatrix; }
            inline math::Matrix<doubleType> getCartesianToRepresentationTransitionMatrix() const { return this->CartesianToRepresentationTransitionMatrix; }
            inline std::string getName() const { return this->name; }
            inline std::vector<std::string> getStateNames() { return this->stateNames; }

            inline void setmu(const double& mu) { this->mu = mu; }

            inline void setStateVectorThisRepresentation(const math::Matrix<doubleType>& StateVectorThisRepresentation)
            { 
                for (size_t stateIndex : {0, 1, 2, 3, 4, 5})
                    this->StateVectorThisRepresentation(stateIndex) = StateVectorThisRepresentation(stateIndex);
            }
            inline void setStateVectorCartesian(const math::Matrix<doubleType>& StateVectorCartesian)
            {
                for (size_t stateIndex : {0, 1, 2, 3, 4, 5})
                    this->StateVectorCartesian(stateIndex) = StateVectorCartesian(stateIndex);
            }

            //data
        protected:
            std::string name;
            std::vector<std::string> stateNames;
            math::Matrix<doubleType> StateVectorThisRepresentation;
            math::Matrix<doubleType> StateVectorCartesian;
            math::Matrix<doubleType> RepresentationToCartesianTransitionMatrix;
            math::Matrix<doubleType> CartesianToRepresentationTransitionMatrix;
            double mu;
        };
    }//end namespace StateRepresentation
}//end namespace EMTG

#endif // STATE_REPRESENTATION_H