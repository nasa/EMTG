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

#ifndef ORBIT_ELEMENT_CONSTRAINT_BASE_H
#define ORBIT_ELEMENT_CONSTRAINT_BASE_H

#include "EMTG_enums.h"
#include "SpecializedBoundaryConstraintBase.h"


namespace EMTG
{
    namespace BoundaryEvents
    {

        namespace SpecializedConstraints
        {

            class OrbitElementConstraintBase : public SpecializedBoundaryConstraintBase
            {
            public:
                //constructors
				OrbitElementConstraintBase() {};

				OrbitElementConstraintBase(const std::string& name,
										   const size_t& journeyIndex,
										   const size_t& phaseIndex,
										   const size_t& stageIndex,
										   Astrodynamics::universe* Universe,
										   HardwareModels::Spacecraft* mySpacecraft,
										   missionoptions* myOptions,
										   BoundaryEventBase* myBoundaryEvent,
										   const std::string& constraintDefinition);

                //destructor
                virtual ~OrbitElementConstraintBase() {};

                //public methods
                virtual void output(std::ofstream& outputfile) = 0;

                //calcbounds goes in the specialized event
                virtual void calcbounds() = 0;

                //process
                virtual void process_constraint(const std::vector<doubleType>& X,
                                                size_t& Xindex,
                                                std::vector<doubleType>& F,
                                                size_t& Findex,
                                                std::vector<double>& G,
                                                const bool& needG) = 0;

            protected:
                //fields
				math::Matrix<doubleType> boundary_orbit_elements;
				math::Matrix<doubleType> boundary_orbit_elements_Jacobian;
                ReferenceFrame user_specified_frame;
            };
        }//end namespace SpecializedConstraints
    }//end namespace BoundaryEvents
}//end namespace EMTG

#endif // ORBIT_ELEMENT_CONSTRAINT_BASE_H