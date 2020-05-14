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

#ifndef PERIAPSE_DISTANCE_CONSTRAINT_H
#define PERIAPSE_DISTANCE_CONSTRAINT_H

#include "OrbitElementConstraintBase.h"

namespace EMTG
{
    namespace BoundaryEvents
    {

        namespace SpecializedConstraints
        {

            class PeriapseDistanceConstraint : public OrbitElementConstraintBase
            {
            public:
                //constructors
				PeriapseDistanceConstraint() {};

				PeriapseDistanceConstraint(const std::string& name,
                                        const size_t& journeyIndex,
                                        const size_t& phaseIndex,
                                        const size_t& stageIndex,
                                        Astrodynamics::universe* Universe,
                                        HardwareModels::Spacecraft* mySpacecraft,
                                        missionoptions* myOptions,
                                        BoundaryEventBase* myBoundaryEvent,
                                        const std::string& constraintDefinition);

                //destructor
                virtual ~PeriapseDistanceConstraint() {};


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
                // fields
				doubleType semimajor_axis;
				doubleType eccentricity;
                doubleType periapse_distance;
                
                std::vector< std::vector<size_t> > Gindex_periapse_distance_position_wrt_StateAfterEvent_variables;//stateIndex, Gindex
                std::vector< std::vector<size_t> > Gindex_periapse_distance_position_wrt_StateAfterEvent_time_variables;//stateIndex, Gindex - I suspect this will be empty
                std::vector< std::vector<size_t> > Gindex_periapse_distance_velocity_wrt_StateAfterEvent_variables;//stateIndex, Gindex
                std::vector< std::vector<size_t> > Gindex_periapse_distance_velocity_wrt_StateAfterEvent_time_variables;//stateIndex, Gindex - I suspect this will be empty
				
				std::vector< std::vector<size_t> > dIndex_periapse_distance_position_wrt_StateAfterEvent;//stateIndex, dIndex
                std::vector< std::vector<size_t> > dIndex_periapse_distance_position_wrt_StateAfterEvent_wrt_Time;//stateIndex, dIndex
                std::vector< std::vector<size_t> > dIndex_periapse_distance_velocity_wrt_StateAfterEvent;//stateIndex, dIndex
                std::vector< std::vector<size_t> > dIndex_periapse_distance_velocity_wrt_StateAfterEvent_wrt_Time;//stateIndex, dIndex
            };
        } // end namespace SpecializedConstraints
    } // end namespace BoundaryEvents
} // end namespace EMTG

#endif // PERIAPSE_DISTANCE_CONSTRAINT_H