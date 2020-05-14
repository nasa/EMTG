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

#pragma once

#include "SpecializedBoundaryConstraintBase.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        class BoundaryEventBase;

        namespace SpecializedConstraints
        {

            class BoundaryDistanceConstraint : virtual public SpecializedBoundaryConstraintBase
            {
            public:
                //constructors
                BoundaryDistanceConstraint() {};

                BoundaryDistanceConstraint(const std::string& name,
                    const size_t& journeyIndex,
                    const size_t& phaseIndex,
                    const size_t& stageIndex,
                    Astrodynamics::universe* Universe,
                    HardwareModels::Spacecraft* mySpacecraft,
                    missionoptions* myOptions,
                    BoundaryEventBase* myBoundaryEvent,
                    const std::string& constraintDefinition);

                //destructor
                virtual ~BoundaryDistanceConstraint() {};


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
                doubleType Distance;
                math::Matrix<doubleType> r_body_CB;
                math::Matrix<doubleType> r_sc_CB;
                math::Matrix<doubleType> r_sc_body;
                Astrodynamics::body* myBody;

                bool isCentralBodyConstraint;
                std::string bodyName;
                math::Matrix<double> dr_sc_CB_dt;
                math::Matrix<double> dr_body_CB_dt;
                double dDistance_dt;

                std::vector<size_t> Gindex_distance_constraint_wrt_time_variables;
                std::vector< std::vector<size_t> > Gindex_distance_constraint_wrt_StateAfterEvent_variables;//stateIndex, Gindex
                std::vector< std::vector<size_t> > Gindex_distance_constraint_wrt_StateAfterEvent_time_variables;//stateIndex, Gindex

                std::vector< std::vector<size_t> > dIndex_distance_with_respect_to_StateAfterEvent;//stateIndex, dIndex
                std::vector< std::vector<size_t> > dIndex_distance_with_respect_to_StateAfterEvent_wrt_Time;//stateIndex, dIndex
            };
        }//end namespace SpecializedConstraints
    }//end namespace BoundaryEvents
}//end namespace EMTG