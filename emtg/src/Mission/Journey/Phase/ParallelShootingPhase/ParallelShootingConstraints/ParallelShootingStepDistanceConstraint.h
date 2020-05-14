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

//parallel shooting distance constraint
//Jacob Englander 2/21/2018
#pragma once

#include <tuple>
#include <string>

#include "missionoptions.h"
#include "universe.h"

#include "sparsey_thing.h"

#include "ParallelShootingStep.h"

namespace EMTG
{
    namespace Phases
    {
        //forward declare step
        class ParallelShootingStep;

        class ParallelShootingStepDistanceConstraint : public sparsey_thing
        {
        public:
            ParallelShootingStepDistanceConstraint();
            ParallelShootingStepDistanceConstraint(const std::string& ConstraintDefinition,
                const size_t& journeyIndex,
                const size_t& phaseIndex,
                const size_t& stepIndex,
                ParallelShootingStep* myStep,
                missionoptions* myOptions,
                Astrodynamics::universe* myUniverse);

            virtual ~ParallelShootingStepDistanceConstraint() {};

            //clone
            virtual ParallelShootingStepDistanceConstraint* clone() const { return new ParallelShootingStepDistanceConstraint(*this); }

            void calcbounds();

            void process(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

        protected:
            void parse_constraint_definition(const std::string& ConstraintDefinition);


            //fields
            std::string name;
            ParallelShootingStep* myStep;
            Astrodynamics::universe* myUniverse;
            missionoptions* myOptions;
            JourneyOptions* myJourneyOptions;
            size_t stepIndex;
            size_t journeyIndex;
            size_t phaseIndex;

            Astrodynamics::body* myBody;

            bool isDefinedRelativeToCentralBody;

            std::tuple<int, double, double> distance_constraint_definition; //body, lower bound in km, upper bound in km
            math::Matrix<doubleType> distance_constraint_relative_position;
            math::Matrix<double> distance_constraint_body_position_time_derivatives;
            doubleType distance_from_body;
            double lowerBound, upperBound;

            std::vector< std::vector<size_t> > dIndex_distance_constraints_wrt_StepLeftPosition; //varIndex, dIndex
            std::vector<size_t> G_indices_distance_constraints_wrt_StepLeftPosition;
            std::vector<size_t> G_indices_distance_constraints_wrt_StepLeftTime;
        };//end class ParallelShootingStepDistanceConstraint


        inline ParallelShootingStepDistanceConstraint * new_clone(ParallelShootingStepDistanceConstraint const & other)
        {
            return other.clone();
        }
    }//end namespace Phases
}//end namespace EMTG