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

#include "BoundaryEventBase.h"
#include "PropagatorBase.h"
#include "StateRepresentation.h"

#include "IntegrationScheme.h"
#include "TimeDomainSpacecraftEOM.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        class FreePointBoundary : public virtual BoundaryEventBase
        {
        public:
            //constructor
            FreePointBoundary() ;
            FreePointBoundary(const std::string& name,
                                const size_t& journeyIndex,
                                const size_t& phaseIndex,
                                size_t& stageIndex,
                                Astrodynamics::universe* Universe,
                                HardwareModels::Spacecraft* mySpacecraft,
                                missionoptions* myOptions);

            virtual void initialize(const std::string& name,
                                    const size_t& journeyIndex,
                                    const size_t& phaseIndex,
                                    size_t& stageIndex,
                                    Astrodynamics::universe* Universe,
                                    HardwareModels::Spacecraft* mySpacecraft,
                                    missionoptions* myOptions);
            
            //destructor
            ~FreePointBoundary();

            //output
            virtual void output(std::ofstream& outputfile,
                                const double& launchdate,
                                size_t& eventcount) = 0;
            //calcbounds

        protected:
            virtual void calcbounds_event_left_side(const std::vector< std::tuple<double, double> >& StateBounds, std::vector<size_t> timeVariables);
            virtual void calcbounds_event_right_side() = 0;

            virtual void process_event_left_side(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG) = 0;

            virtual void process_event_right_side(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG) = 0;

            //fields

            //propagation
            bool AllowStateToPropagate;
            Astrodynamics::PropagatorBase* myPropagator;
            math::Matrix<doubleType> StateBeforeEventBeforePropagation;
            math::Matrix<double> STM;
            math::Matrix<double> dPropagatedStatedIndependentVariable;
            double dPropagationTime_dIndependentVariable;
            double ReferenceEpoch;

            //integrator pointers
            Astrodynamics::SpacecraftAccelerationModel* mySpacecraftAccelerationModel;
            Astrodynamics::TimeDomainSpacecraftEOM myEOM;
            Integration::IntegrationScheme* myIntegrationScheme;
            size_t total_number_of_states_to_integrate;

            //states
            StateRepresentation myStateRepresentationEnum;
            Astrodynamics::StateRepresentationBase* myStateRepresentation;
            std::vector<bool> encodedStateIsAngle;
            math::Matrix<doubleType> StateBeforeEventEncoded;
            std::vector<size_t> Xindex_encoded_state;
            ReferenceFrame myEncodedReferenceFrame;
            std::vector< std::vector<size_t> > dIndex_StateBeforeEvent_wrt_EncodedState;//encoded state, state before event
            std::vector< std::vector<size_t> > dIndex_StateBeforeEvent_wrt_Time;//state, timeVariable

            //8x8 matrix describing the derivatives of the state before event with respect to decision variables
            //this bakes in EVERYTHING, including ephemeris lookup vs encoded state, free point rotation and propagation, etc.
            math::Matrix<doubleType> dStateBeforeEvent_dStateBeforeEventEncoded;

            //derivatives of state before event with respect to decision variables, when chained through rotation and propagation
            math::Matrix<doubleType> dStateBeforeEventEncoded_dEncodedElements;//if cartesian state then this is 8x8 identity
                                                                               //otherwise (COE) it has a 1 in (6,6) and (7,7),
                                                                               //the upper 6x6 is orbit element conversions, and everything else is zeros
            math::Matrix<doubleType> dStateBeforeEventICRF_dStateBeforeEventEncoded;//if encoded state then this is 8x8 identity
                                                                                   //otherwise (ephemeris) it has a 1 in (6,6) and (7,7),
                                                                                   //  the (i=0,5, j=7) column is dense, and otherwise zeros 
            math::Matrix<doubleType> dStateBeforeEvent_dStateBeforeEventICRF;//upper 6x6 is rotation matrix, (6,6) and (7,7) are 1,
                                                                                   //(i=0,5, j=7) column is derivatives of rotated state wrt epoch
                                                                                               
            //state rotation matrix 8 x 8
            math::Matrix<doubleType> R_from_local_to_ICRF;
            math::Matrix<doubleType> dR_from_local_to_ICRF_dt;

            //tracking of derivatives for central body exclusion constraint
            std::vector<size_t> Gindices_dRBeforeEvent_dStateBeforeEventEncoded;//stateindex
        };
    }//close namespace events
}//close namespace EMTG