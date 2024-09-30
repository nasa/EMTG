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

//EMTGv9 ProbeEntryPhase
//Jacob Englander 2-28-2020

#pragma once

#include "MGAnDSMs_phase.h"
#include "FreePointLTRendezvous.h"
#include "FreePointFreeDirectDeparture.h"

namespace EMTG 
{
    namespace Phases
    {
        class ProbeEntryPhase : public MGAnDSMs_phase
        {
        public:
            //constructor
            ProbeEntryPhase() : MGAnDSMs_phase::MGAnDSMs_phase() {};

            ProbeEntryPhase(const std::string& name,
                const size_t& journeyIndex,
                const size_t& phaseIndex,
                size_t& stageIndex,
                phase* previousPhase,
                Astrodynamics::universe* Universe,
                HardwareModels::Spacecraft* mySpacecraft,
                HardwareModels::LaunchVehicle* myLaunchVehicle,
                missionoptions* myOptions);

            //destructor
            virtual ~ProbeEntryPhase();

            //output
            virtual void output(std::ofstream& outputfile,
                size_t& eventcount);

            void output_probe();

            void output_STMs();

            void output_maneuver_and_target_spec(std::ofstream& maneuver_spec_file, std::ofstream& target_spec_file, bool& haveManeuverNeedTarget);

            void output_ephemeris(std::ofstream& outputfile, std::ofstream& acceleration_model_file);

            void output_probe_ephemeris();

            virtual void calcbounds();

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

            //process goes in the specialized phase
            virtual void process_phase(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

        protected:
            virtual void calcbounds_phase_main();

            virtual void calcbounds_probe_descent_time();

            virtual void calcbounds_probe_boundary_events();

            virtual void calcbounds_separation();

            virtual void calcbounds_spacecraft_match_point_derivatives_due_to_probe_release();

            virtual void calcbounds_probe_match_point_constraints_from_separation_to_AEI();

            virtual void calcbounds_probe_match_point_constraints_from_AEI_to_end();

            virtual void calcbounds_communication_distance_constraint();

            void process_phase_main(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            virtual void process_probe_descent_time(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            virtual void process_probe_boundary_events(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            virtual void process_separation(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            virtual void process_spacecraft_match_point_derivatives_due_to_probe_release(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            virtual void process_probe_separation_to_AEI_forward_half_phase(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            virtual void process_probe_separation_to_AEI_backward_half_phase(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            virtual void process_probe_AEI_to_end_forward_half_phase(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            virtual void process_probe_AEI_to_end_backward_half_phase(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            virtual void process_probe_match_point_constraints_from_separation_to_AEI(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            virtual void process_probe_match_point_constraints_from_AEI_to_end(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            virtual void process_communication_distance_constraint(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            //fields
            doubleType probeDescentTime;
            size_t Xindex_probeDescentTime;

            BoundaryEvents::FreePointLTRendezvous* probeEntryArrivalEvent;
            BoundaryEvents::FreePointFreeDirectDeparture* probeEntryDepartureEvent;
            BoundaryEvents::FreePointLTRendezvous* probeEndEvent;
            missionoptions probe_separation_to_AEI_MissionOptions;//doctored journeyOptions for the probe before AEI, fed to its arrival event
            missionoptions probe_AEI_to_end_MissionOptions;//doctored journeyOptions for the probe after AEI, fed to its arrival event

            math::Matrix<doubleType> probe_first_subphase_match_point_state_minus;
            math::Matrix<doubleType> probe_first_subphase_match_point_state_plus;
            math::Matrix<doubleType> probe_second_subphase_match_point_state_minus;
            math::Matrix<doubleType> probe_second_subphase_match_point_state_plus;
            math::Matrix<doubleType> spacecraft_state_after_separation;
            math::Matrix<doubleType> probe_state_after_separation;
            math::Matrix<doubleType> probe_state_at_AEI_arrival;
            math::Matrix<doubleType> probe_state_at_AEI_departure;
            math::Matrix<doubleType> probe_state_at_end_of_phase;

            math::Matrix<doubleType> spacecraft_separation_deltaV;
            math::Matrix<doubleType> probe_separation_deltaV;

            math::Matrix<double> probeReleaseTM;//transition matrix for the probe from the left-hand side of the phase (unseparated, pre-TCM) to the post-TCM, post separation state

            math::Matrix<double> probe_separation_to_AEI_ForwardHPTM;
            math::Matrix<double> probe_separation_to_AEI_BackwardHPTM;
            math::Matrix<double> probe_AEI_to_end_ForwardHPTM;
            math::Matrix<double> probe_AEI_to_end_BackwardHPTM;
                        
            Astrodynamics::SpacecraftAccelerationModel* probe_separation_to_AEI_SpacecraftAccelerationModel;
            Astrodynamics::SpacecraftAccelerationModel* probe_AEI_to_end_SpacecraftAccelerationModel;

            Integration::IntegrationScheme* probe_separation_to_AEI_IntegrationScheme;
            Integration::IntegrationScheme* probe_AEI_to_end_IntegrationScheme;

            Astrodynamics::PropagatorBase* probe_separation_to_AEI_ForwardHalfPhasePropagator;
            Astrodynamics::PropagatorBase* probe_separation_to_AEI_BackwardHalfPhasePropagator;
            Astrodynamics::PropagatorBase* probe_AEI_to_end_ForwardHalfPhasePropagator;
            Astrodynamics::PropagatorBase* probe_AEI_to_end_BackwardHalfPhasePropagator;
            double probe_separation_to_AEI_ForwardIntegrationStepLength, probe_separation_to_AEI_BackwardIntegrationStepLength;
            double probe_AEI_to_end_ForwardIntegrationStepLength, probe_AEI_to_end_BackwardIntegrationStepLength;
            double probe_separation_to_AEI_MatchPointFraction;
            double probe_AEI_to_end_MatchPointFraction;

            size_t numMatchConstraints_probe;

            size_t total_number_of_states_to_integrate;

            math::Matrix<double> probe_separation_to_AEI_ForwardSTM, probe_separation_to_AEI_BackwardSTM, probe_separation_to_AEI_ForwardSPTM, probe_separation_to_AEI_BackwardSPTM;
            math::Matrix<double> probe_AEI_to_end_ForwardSTM, probe_AEI_to_end_BackwardSTM, probe_AEI_to_end_ForwardSPTM, probe_AEI_to_end_BackwardSPTM;
            math::Matrix<double> probe_separation_to_AEI_Forward_dPropagatedStatedIndependentVariable, probe_separation_to_AEI_Backward_dPropagatedStatedIndependentVariable;
            math::Matrix<double> probe_AEI_to_end_Forward_dPropagatedStatedIndependentVariable, probe_AEI_to_end_Backward_dPropagatedStatedIndependentVariable;
            double probe_separation_to_AEI_dForwardStepSize_dPropagationVariable;
            double probe_AEI_to_end_dForwardStepSize_dPropagationVariable;
            double probe_separation_to_AEI_dBackwardStepSize_dPropagationVariable;
            double probe_AEI_to_end_dBackwardStepSize_dPropagationVariable;

            Astrodynamics::TimeDomainSpacecraftEOM probe_separation_to_AEI_EOM;
            Astrodynamics::TimeDomainSpacecraftEOM probe_AEI_to_end_EOM;

            double probeSpacecraftCommunicationDistance;

            bool isKeplerian;
        
            //constraint indices
            std::vector<size_t> Findices_match_point_constraints_probe_separation_to_AEI;
            std::vector<size_t> Findices_match_point_constraints_probe_AEI_to_end;

            //derivative indices
            //separation impulse direction constraint
            std::vector<size_t> Xindices_separation_direction_vector;
            std::vector<size_t> Gindices_separation_direction_vector_wrt_self;
            std::vector< std::vector<size_t> > Gindex_separation_vector_velocity_match_constraint;//vIndex, Xindex
            std::vector< std::vector<size_t> > Gindex_separation_vector_velocity_match_constraint_wrt_Time;//vIndex, Xindex
            std::vector< std::vector<size_t> > dIndex_separation_vector_velocity_match_constraint;//vIndex, dIndex
            std::vector< std::vector<size_t> > dIndex_separation_vector_velocity_match_constraint_wrt_Time;//vIndex, dIndex

            //derivatives of probe separation to AEI match point due to boundary variables
            std::vector<size_t> ListOfVariablesAffectingRightBoundary_probe_separation_to_AEI;
            std::vector< std::vector<size_t> > DerivativesOfLeftBoundaryByVariable_probe_separation_to_AEI; //localXindex, dIndex
            std::vector< std::vector<size_t> > DerivativesOfRightBoundaryByVariable_probe_separation_to_AEI; //localXindex, dIndex
            std::vector< std::vector<size_t> > G_indices_match_point_constraints_wrt_LeftBoundaryVariables_probe_separation_to_AEI;
            std::vector< std::vector<size_t> > G_indices_match_point_constraints_wrt_RightBoundaryVariables_probe_separation_to_AEI;
            std::vector<size_t> ListOfTimeVariablesAffectingRightBoundary_probe_separation_to_AEI;
            std::vector< std::vector<size_t> > DerivativesOfLeftBoundaryByTimeVariable_probe_separation_to_AEI; //localXindex, dIndex
            std::vector< std::vector<size_t> > DerivativesOfRightBoundaryByTimeVariable_probe_separation_to_AEI; //localXindex, dIndex
            std::vector< std::vector<size_t> > G_indices_match_point_constraints_wrt_LeftBoundaryTimeVariables_probe_separation_to_AEI;
            std::vector< std::vector<size_t> > G_indices_match_point_constraints_wrt_RightBoundaryTimeVariables_probe_separation_to_AEI;
            std::vector<size_t> G_indices_match_point_constraints_wrt_PhaseFlightTime_probe_separation_to_AEI;

            std::vector< std::vector<bool> > LeftBoundaryVariableAffectsMatchConstraint_probe_separation_to_AEI;//constraintIndex, entry in ListOfVariablesAffectingLeftBoundary_probe
            std::vector< std::vector<bool> > LeftBoundaryTimeVariableAffectsMatchConstraint_probe_separation_to_AEI;//constraintIndex, entry in ListOfTimeVariablesAffectingLeftBoundary_probe
            std::vector< std::vector<bool> > RightBoundaryVariableAffectsMatchConstraint_probe_separation_to_AEI;//constraintIndex, entry in ListOfVariablesAffectingRightBoundary_probe
            std::vector< std::vector<bool> > RightBoundaryTimeVariableAffectsMatchConstraint_probe_separation_to_AEI;//constraintIndex, entry in ListOfTimeVariablesAffectingRightBoundary_probe


            //derivatives of probe separation to AEI match point due to boundary variables
            std::vector<size_t> ListOfVariablesAffectingLeftBoundary_probe_AEI_to_end;
            std::vector<size_t> ListOfVariablesAffectingRightBoundary_probe_AEI_to_end;
            std::vector< std::vector<size_t> > DerivativesOfLeftBoundaryByVariable_probe_AEI_to_end; //localXindex, dIndex
            std::vector< std::vector<size_t> > DerivativesOfRightBoundaryByVariable_probe_AEI_to_end; //localXindex, dIndex
            std::vector< std::vector<size_t> > G_indices_match_point_constraints_wrt_LeftBoundaryVariables_probe_AEI_to_end;
            std::vector< std::vector<size_t> > G_indices_match_point_constraints_wrt_RightBoundaryVariables_probe_AEI_to_end;
            std::vector<size_t> ListOfTimeVariablesAffectingLeftBoundary_probe_AEI_to_end;
            std::vector<size_t> ListOfTimeVariablesAffectingRightBoundary_probe_AEI_to_end;
            std::vector< std::vector<size_t> > DerivativesOfLeftBoundaryByTimeVariable_probe_AEI_to_end; //localXindex, dIndex
            std::vector< std::vector<size_t> > DerivativesOfRightBoundaryByTimeVariable_probe_AEI_to_end; //localXindex, dIndex
            std::vector< std::vector<size_t> > G_indices_match_point_constraints_wrt_LeftBoundaryTimeVariables_probe_AEI_to_end;
            std::vector< std::vector<size_t> > G_indices_match_point_constraints_wrt_RightBoundaryTimeVariables_probe_AEI_to_end;
            std::vector<size_t> G_indices_match_point_constraints_wrt_PhaseFlightTime_probe_AEI_to_end;
            std::vector< std::vector<size_t> > TruthTable_Forward_MatchConstraints_Derivative_wrt_EncodedStates_first_subphase;

            std::vector< std::vector<bool> > LeftBoundaryVariableAffectsMatchConstraint_probe_AEI_to_end;//constraintIndex, entry in ListOfVariablesAffectingLeftBoundary_probe
            std::vector< std::vector<bool> > LeftBoundaryTimeVariableAffectsMatchConstraint_probe_AEI_to_end;//constraintIndex, entry in ListOfTimeVariablesAffectingLeftBoundary_probe
            std::vector< std::vector<bool> > RightBoundaryVariableAffectsMatchConstraint_probe_AEI_to_end;//constraintIndex, entry in ListOfVariablesAffectingRightBoundary_probe
            std::vector< std::vector<bool> > RightBoundaryTimeVariableAffectsMatchConstraint_probe_AEI_to_end;//constraintIndex, entry in ListOfTimeVariablesAffectingRightBoundary_probe

            //derivatives due to separation impulse
            std::vector< std::vector<size_t> > Gindex_spacecraft_match_point_constraints_wrt_separation_impulse;//stateIndex, velocityIndex
            std::vector< size_t > Gindex_spacecraft_match_point_constraints_wrt_left_mass;//stateIndex, velocityIndex
            std::vector< size_t > dIndex_spacecraft_left_mass;//stateIndex, velocityIndex
            std::vector< std::vector<size_t> > Gindex_probe_match_point_constraints_wrt_separation_impulse;//stateIndex, velocityIndex

            //derivatives of communications constraint
            std::vector< std::vector<size_t> > Gindex_probe_communication_distanceconstraint_wrt_spacecraftStateAfterPhase_variables;//stateIndex, Gindex
            std::vector< std::vector<size_t> > Gindex_probe_communication_distanceconstraint_wrt_time_variables;//stateIndex, Gindex
            std::vector< std::vector<size_t> > Gindex_probe_communication_distanceconstraint_wrt_probeStateAtAEntry_variables;//stateIndex, Gindex

            std::vector< std::vector<size_t> > dIndex_probe_communication_distanceconstraint_with_respect_to_spacecraftStateAfterPhase;//stateIndex, dIndex
            std::vector< std::vector<size_t> > dIndex_probe_communication_distanceconstraint_with_respect_to_time_variables;//stateIndex, dIndex
            std::vector< std::vector<size_t> > dIndex_probe_communication_distanceconstraint_with_respect_to_probeStateAtEntry;//stateIndex, dIndex
            
        };
    }//close namespace Phases
}//close namespace EMTG