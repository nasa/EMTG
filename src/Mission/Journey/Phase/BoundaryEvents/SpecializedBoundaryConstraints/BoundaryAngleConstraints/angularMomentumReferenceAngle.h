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

            class angularMomentumReferenceAngle : virtual public SpecializedBoundaryConstraintBase
            {
            public:
                //constructors
                angularMomentumReferenceAngle() {};

                angularMomentumReferenceAngle(const std::string& name,
                    const size_t& journeyIndex,
                    const size_t& phaseIndex,
                    const size_t& stageIndex,
                    Astrodynamics::universe* Universe,
                    HardwareModels::Spacecraft* mySpacecraft,
                    missionoptions* myOptions,
                    BoundaryEventBase* myBoundaryEvent,
                    const std::string& constraintDefinition);

                //destructor
                virtual ~angularMomentumReferenceAngle() {};


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
                doubleType Angle;
                math::Matrix<doubleType> hvec;
                math::Matrix<doubleType> R_probe_ref;
                math::Matrix<doubleType> R_cb_ref;
                Astrodynamics::body* myReference;

                bool refIsCentralBody;
                std::string refName;


                math::Matrix<double> dR_cb_ref_dt;
                math::Matrix<double> dR_cb_probe_dt;
                math::Matrix<double> dV_cb_probe_dt;

                //this constraint has derivatives with respect to any non-time variable that affects spacecraft position and velocity
                //spacecraft position is used for "Body" vector
                //spacecraft velocity is used for "Probe" vector
                //this constraint has derivatives with respect to all time variables because they affect both spacecraft and reference
                std::vector<size_t> Gindex_constraint_wrt_time_variables;
                std::vector< std::vector<size_t> > Gindex_constraint_wrt_StateAfterEvent_variables;//stateIndex, Gindex
                std::vector< std::vector<size_t> > Gindex_constraint_wrt_StateAfterEvent_time_variables;//stateIndex, Gindex

                std::vector< std::vector<size_t> > dIndex_with_respect_to_StateAfterEvent;//stateIndex, dIndex
                std::vector< std::vector<size_t> > dIndex_with_respect_to_StateAfterEvent_wrt_Time;//stateIndex, dIndex
            };
        }//end namespace SpecializedConstraints
    }//end namespace BoundaryEvents
}//end namespace EMTG