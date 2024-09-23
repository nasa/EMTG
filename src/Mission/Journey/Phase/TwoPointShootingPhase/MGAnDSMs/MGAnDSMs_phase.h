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

//EMTGv9 MGAnDSMs phase
//Jacob Englander 8-25-2017

#pragma once

#include "TwoPointShootingPhase.h"

#include "SpacecraftAccelerationModel.h"
#include "TimeDomainSpacecraftEOM.h"
#include "IntegrationScheme.h"

#include "boost/ptr_container/ptr_vector.hpp"

#include "Forward_MGAnDSMs_subphase.h"
#include "Backward_MGAnDSMs_subphase.h"

namespace EMTG
{
    namespace Phases
    {
        class MGAnDSMs_phase : public TwoPointShootingPhase
        {
        public:
            //constructor
            MGAnDSMs_phase() : TwoPointShootingPhase::TwoPointShootingPhase() {};
            MGAnDSMs_phase(const std::string& name,
                const size_t& journeyIndex,
                const size_t& phaseIndex,
                size_t& stageIndex,
                phase* previousPhase,
                Astrodynamics::universe* Universe,
                HardwareModels::Spacecraft* mySpacecraft,
                HardwareModels::LaunchVehicle* myLaunchVehicle,
                missionoptions* myOptions);

            //destructor
            virtual ~MGAnDSMs_phase();

            //clone
            virtual MGAnDSMs_phase* clone() const { return new MGAnDSMs_phase(*this); }
            
            //output
            virtual void output(std::ofstream& outputfile,
                size_t& eventcount);

            virtual void output_STMs();

            virtual void output_maneuver_and_target_spec(std::ofstream& maneuver_spec_file, std::ofstream& target_spec_file, bool& haveManeuverNeedTarget);

            //calcbounds
            virtual void setup_calcbounds(std::vector<double>* Xupperbounds,
                std::vector<double>* Xlowerbounds,
                std::vector<double>* X_scale_factors,
                std::vector<double>* Fupperbounds,
                std::vector<double>* Flowerbounds,
				std::vector<double>* F_scale_factors,
                std::vector<std::string>* Xdescriptions,
                std::vector<std::string>* Fdescriptions,
                std::vector<size_t>* iGfun,
                std::vector<size_t>* jGvar,
                std::vector<std::string>* Gdescriptions,
                std::vector<size_t>* iAfun,
                std::vector<size_t>* jAvar,
                std::vector<std::string>* Adescriptions,
                std::vector<double>* A);

            virtual void calcbounds();

            virtual void output_ephemeris(std::ofstream& outputfile, std::ofstream& acceleration_model_file);

            //process goes in the specialized phase
            virtual void process_phase(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

        protected:
            //calcbounds
            virtual void calcbounds_phase_main();

            virtual void calcbounds_subphases();

            virtual void calcbounds_match_point_constraints();

            void calcbounds_burnindex_sum_constraint();

            void calcbounds_virtual_propellant_tanks();

            void calcbounds_deltav_contribution();

            //process
            virtual void process_phase_main(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            virtual void process_forward_half_phase(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            virtual void process_backward_half_phase(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            virtual void process_match_point_constraints(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            void process_burnindex_sum_constraint(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            void process_virtual_propellant_tanks(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            void process_deltav_contribution(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            void output_phase_main(std::ofstream& outputfile,
                size_t& eventcount);

            //fields
            boost::ptr_vector< Forward_MGAnDSMs_subphase > ForwardSubPhases;
            boost::ptr_vector< Backward_MGAnDSMs_subphase > BackwardSubPhases;

            std::vector< math::Matrix<doubleType> > SpacecraftStateForward;
            std::vector< math::Matrix<doubleType> > SpacecraftStateBackward;

            std::vector< std::vector< std::vector<size_t> > > Gindices_match_point_constraint_ForwardControl;//forwardSubPhaseIndex, constraint, control variable
            std::vector< std::vector< std::vector<size_t> > > Gindices_match_point_constraint_BackwardControl;//backwardSubPhaseIndex, constraint, control variable
            std::vector< std::vector<size_t> > Gindices_match_point_constraint_ForwardBurnIndex;//forwardSubPhaseIndex, constraint
            std::vector< std::vector<size_t> > Gindices_match_point_constraint_BackwardBurnIndex;//forwardSubPhaseIndex, constraint

            std::vector< math::Matrix<double> > ForwardSPTM;
            std::vector< math::Matrix<double> > BackwardSPTM;
            std::vector< math::Matrix<double> > CumulativeForwardSPTM;
            std::vector< math::Matrix<double> > CumulativeBackwardSPTM;
            math::Matrix<double> Forward_dPropagatedState_dIndependentVariable;
            math::Matrix<double> Backward_dPropagatedState_dIndependentVariable;

            size_t numberOfDSMs;
            size_t numberOfForwardSubphases;
            size_t numberOfBackwardSubphases;
            size_t subphaseMatchPointIndex;

            std::vector< std::vector<size_t> > Gindices_deltav_wrt_ForwardControl;
            std::vector< std::vector<size_t> > Gindices_deltav_wrt_BackwardControl;

            std::vector<size_t> Gindices_BurnIndexSumConstraint;
        };
    }//close namespace Phases
}//close namespace EMTG