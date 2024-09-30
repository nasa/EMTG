// EMTG: Evolutionary Mission Trajectory Generator
// An open-source global optimization tool for preliminary mission design
// Provided by NASA Goddard Space Flight Center
//
// Copyright (c) 2013 - 2024 United States Government as represented by the
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

#pragma once

#include "SpecializedBoundaryConstraintBase.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        class BoundaryEventBase;

        namespace SpecializedConstraints
        {

            class BoundaryStateInTwoBodyRotatingFrameConstraint : virtual public SpecializedBoundaryConstraintBase
            {
            public:
                //constructors
				BoundaryStateInTwoBodyRotatingFrameConstraint() {};

				BoundaryStateInTwoBodyRotatingFrameConstraint(const std::string& name,
                    const size_t& journeyIndex,
                    const size_t& phaseIndex,
                    const size_t& stageIndex,
                    Astrodynamics::universe* Universe,
                    HardwareModels::Spacecraft* mySpacecraft,
                    missionoptions* myOptions,
                    BoundaryEventBase* myBoundaryEvent,
                    const std::string& constraintDefinition);

                //destructor
                virtual ~BoundaryStateInTwoBodyRotatingFrameConstraint() {};


                //public methods
                virtual void output(std::ofstream& outputfile);

                virtual void calcbounds();

                //process
                virtual void process_constraint(const std::vector<doubleType>& X,
                    size_t& Xindex,
                    std::vector<doubleType>& F,
                    size_t& Findex,
                    std::vector<double>& G,
                    const bool& needG);

            protected:
                //fields
				size_t stateConstrainedIndex; // index of the state being constrained: 0 to 5 for x to vz
				size_t originOfConstraintFrame; // 1: body1. 2: body2. only body2 supported currently
				std::string nameOfOriginOfConstraintFrame;

				doubleType constrainedStateValue; // the value of the constrained state

				bool body1IsCB; // is body1 the central body?
				bool body2IsCB; // is body2 the central body?
				std::string body1Name;
				std::string body2Name;
				Astrodynamics::body* myBody1;
				Astrodynamics::body* myBody2;

				std::vector<size_t> Gindex_constraint_wrt_time_variables;
				std::vector< std::vector<size_t> > Gindex_constraint_position_wrt_StateAroundEvent_variables;//stateIndex, Gindex
				std::vector< std::vector<size_t> > Gindex_constraint_position_wrt_StateAroundEvent_time_variables;//stateIndex, Gindex

				std::vector< std::vector<size_t> > dIndex_constraint_position_wrt_StateAroundEvent;//stateIndex, dIndex
				std::vector< std::vector<size_t> > dIndex_constraint_position_wrt_StateAroundEvent_wrt_Time;//stateIndex, dIndex

				std::vector< std::vector<size_t> > Gindex_constraint_velocity_wrt_StateAroundEvent_variables;//stateIndex, Gindex
				std::vector< std::vector<size_t> > Gindex_constraint_velocity_wrt_StateAroundEvent_time_variables;//stateIndex, Gindex

				std::vector< std::vector<size_t> > dIndex_constraint_velocity_wrt_StateAroundEvent;//stateIndex, dIndex
				std::vector< std::vector<size_t> > dIndex_constraint_velocity_wrt_StateAroundEvent_wrt_Time;//stateIndex, dIndex

                std::vector<size_t> Gindex_distance_constraint_wrt_time_variables;
                std::vector< std::vector<size_t> > Gindex_distance_constraint_wrt_StateAfterEvent_variables;//stateIndex, Gindex
                std::vector< std::vector<size_t> > Gindex_distance_constraint_wrt_StateAfterEvent_time_variables;//stateIndex, Gindex

                std::vector< std::vector<size_t> > dIndex_distance_with_respect_to_StateAfterEvent;//stateIndex, dIndex
                std::vector< std::vector<size_t> > dIndex_distance_with_respect_to_StateAfterEvent_wrt_Time;//stateIndex, dIndex
            };
        }//end namespace SpecializedConstraints
    }//end namespace BoundaryEvents
}//end namespace EMTG