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

//EMTGv9 ProbeEntryPhase
//Jacob Englander 2-28-2020

#include "ProbeEntryPhase.h"
#include "PropagatorFactory.h"
#include "IntegrationSchemeFactory.h"
#include "EMTG_string_utilities.h"

#include "boost/algorithm/string.hpp"

namespace EMTG
{
    namespace Phases
    {
        ProbeEntryPhase::ProbeEntryPhase(const std::string& name,
            const size_t& journeyIndex,
            const size_t& phaseIndex,
            size_t& stageIndex,
            phase* previousPhase,
            Astrodynamics::universe* Universe,
            HardwareModels::Spacecraft* mySpacecraft,
            HardwareModels::LaunchVehicle* myLaunchVehicle,
            missionoptions* myOptions) :
            MGAnDSMs_phase::MGAnDSMs_phase(name,
                journeyIndex,
                phaseIndex,
                stageIndex,
                previousPhase,
                Universe,
                mySpacecraft,
                myLaunchVehicle,
                myOptions)
        {
            //is this phase Keplerian?
            if (this->myOptions->propagatorType == PropagatorType::KeplerPropagator
                || (this->myJourneyOptions->override_PropagatorType && this->myJourneyOptions->propagatorType == PropagatorType::KeplerPropagator))
            {
                this->isKeplerian = true;
            }
            else
            {
                this->isKeplerian = false;
            }

            //derivative truth table override
            {
                //Probe separation impulse is dependent on spacecraft mass. Therefore match point constraints are dependent on spacecraft mass.
                for (size_t stateIndex = 0; stateIndex < 6; ++stateIndex)
                {
                    this->TruthTable_Forward_MatchConstraints_Derivative_wrt_EncodedStates[6][stateIndex] = true;
                    this->TruthTable_Backward_MatchConstraints_Derivative_wrt_EncodedStates[6][stateIndex] = true;
                }
                this->TruthTable_Forward_MatchConstraints_Derivative_wrt_EncodedStates[6][7] = true; //epoch
                this->TruthTable_Backward_MatchConstraints_Derivative_wrt_EncodedStates[6][7] = true; //epoch

                //mass variables do not affect other constraints
                for (size_t constraintIndex = 0; constraintIndex < 6; ++constraintIndex)
                {
                    this->TruthTable_Forward_MatchConstraints_Derivative_wrt_EncodedStates[constraintIndex][6] = true; //mass variables
                    this->TruthTable_Backward_MatchConstraints_Derivative_wrt_EncodedStates[constraintIndex][6] = true; //mass variables
                }
            }//end truth table



            //we need to assign the appropriate constraints to the spacecraft's and probe's arrival events
            std::vector<std::string> spacecraft_arrival_constraints;
            std::vector<std::string> probe_AEI_arrival_constraints;
            std::vector<std::string> probe_end_arrival_constraints;

            for (std::string& constraint : this->myJourneyOptions->BoundaryConstraintDefinitions)
            {
                //does the string contain the substring, "_probe"? If so, goes to the probe
                if (string_utilities::string_contains_substring(constraint, "_probe_AEI"))
                {
                    probe_AEI_arrival_constraints.push_back(boost::erase_all_copy(constraint, "_probe_AEI"));
                }
                else if (string_utilities::string_contains_substring(constraint, "_probe_end"))
                {
                    probe_end_arrival_constraints.push_back(boost::erase_all_copy(constraint, "_probe_end"));
                }
                else //otherwise it goes to the spacecraft
                {
                    spacecraft_arrival_constraints.push_back(constraint);
                }
            }

            //create the probe AEI right boundary event
            //the probe entry needs its own missionoptions object. This is necessary so that the FreePointBoundary constructor can work as intended.
            this->probe_separation_to_AEI_MissionOptions = *this->myOptions;
            this->probe_separation_to_AEI_MissionOptions.Journeys[this->journeyIndex].arrival_elements_vary_flag = this->probe_separation_to_AEI_MissionOptions.Journeys[this->journeyIndex].Probe_AEI_elements_vary_flag;
            this->probe_separation_to_AEI_MissionOptions.Journeys[this->journeyIndex].arrival_elements = this->probe_separation_to_AEI_MissionOptions.Journeys[this->journeyIndex].Probe_AEI_elements;
            this->probe_separation_to_AEI_MissionOptions.Journeys[this->journeyIndex].arrival_elements_bounds = this->probe_separation_to_AEI_MissionOptions.Journeys[this->journeyIndex].Probe_AEI_elements_bounds;
            this->probe_separation_to_AEI_MissionOptions.Journeys[this->journeyIndex].arrival_elements_reference_epoch = this->probe_separation_to_AEI_MissionOptions.Journeys[this->journeyIndex].Probe_AEI_elements_reference_epoch;
            this->probe_separation_to_AEI_MissionOptions.Journeys[this->journeyIndex].arrival_elements_state_representation = this->probe_separation_to_AEI_MissionOptions.Journeys[this->journeyIndex].Probe_AEI_elements_state_representation;
            this->probe_separation_to_AEI_MissionOptions.Journeys[this->journeyIndex].arrival_elements_frame = this->probe_separation_to_AEI_MissionOptions.Journeys[this->journeyIndex].Probe_AEI_elements_frame;
            this->probe_separation_to_AEI_MissionOptions.Journeys[this->journeyIndex].AllowJourneyFreePointArrivalToPropagate = this->probe_separation_to_AEI_MissionOptions.Journeys[this->journeyIndex].AllowJourneyProbeAEIToPropagate;
            this->probe_separation_to_AEI_MissionOptions.Journeys[this->journeyIndex].perturb_drag = this->probe_separation_to_AEI_MissionOptions.Journeys[this->journeyIndex].perturb_drag_probe_separation_to_AEI;
            this->probe_separation_to_AEI_MissionOptions.Journeys[this->journeyIndex].spacecraft_drag_area = this->probe_separation_to_AEI_MissionOptions.Journeys[this->journeyIndex].probe_drag_area_probe_separation_to_AEI;
            this->probe_separation_to_AEI_MissionOptions.Journeys[this->journeyIndex].coefficient_of_drag = this->probe_separation_to_AEI_MissionOptions.Journeys[this->journeyIndex].probe_coefficient_of_drag_probe_separation_to_AEI;
            this->probe_separation_to_AEI_MissionOptions.Journeys[this->journeyIndex].final_velocity = this->probe_separation_to_AEI_MissionOptions.Journeys[this->journeyIndex].probe_AEI_velocity;

            //initialize the boundary itself
            this->probeEntryArrivalEvent = new BoundaryEvents::FreePointIntercept(name + "ProbeAEIFreePointIntercept",
                this->journeyIndex,
                this->phaseIndex,
                stageIndex,
                this->myUniverse,
                this->mySpacecraft,
                &this->probe_separation_to_AEI_MissionOptions);

            this->myArrivalEvent->construct_boundary_constraints(spacecraft_arrival_constraints);
            this->probeEntryArrivalEvent->construct_boundary_constraints(probe_AEI_arrival_constraints);

            //create the second subphase if the user has requested it
            if (this->myJourneyOptions->ModelProbeSecondPhase)
            {
                //create the probe post-AEI "departure" event
                this->probeEntryDepartureEvent = new BoundaryEvents::FreePointFreeDirectDeparture(this->name + "ProbeAEIFreePointFreeDirectDeparture",
                    this->journeyIndex,
                    this->phaseIndex,
                    stageIndex,
                    this->myUniverse,
                    this->mySpacecraft,
                    &this->probe_separation_to_AEI_MissionOptions,
                    this->probeEntryArrivalEvent);


                //create the probe end of mission right boundary event
                //the probe entry needs its own missionoptions object. This is necessary so that the FreePointBoundary constructor can work as intended.
                this->probe_AEI_to_end_MissionOptions = *this->myOptions;
                this->probe_AEI_to_end_MissionOptions.Journeys[this->journeyIndex].arrival_elements_vary_flag = this->probe_AEI_to_end_MissionOptions.Journeys[this->journeyIndex].Probe_End_elements_vary_flag;
                this->probe_AEI_to_end_MissionOptions.Journeys[this->journeyIndex].arrival_elements = this->probe_AEI_to_end_MissionOptions.Journeys[this->journeyIndex].Probe_End_elements;
                this->probe_AEI_to_end_MissionOptions.Journeys[this->journeyIndex].arrival_elements_bounds = this->probe_AEI_to_end_MissionOptions.Journeys[this->journeyIndex].Probe_End_elements_bounds;
                this->probe_AEI_to_end_MissionOptions.Journeys[this->journeyIndex].arrival_elements_reference_epoch = this->probe_AEI_to_end_MissionOptions.Journeys[this->journeyIndex].Probe_End_elements_reference_epoch;
                this->probe_AEI_to_end_MissionOptions.Journeys[this->journeyIndex].arrival_elements_state_representation = this->probe_AEI_to_end_MissionOptions.Journeys[this->journeyIndex].Probe_End_elements_state_representation;
                this->probe_AEI_to_end_MissionOptions.Journeys[this->journeyIndex].arrival_elements_frame = this->probe_AEI_to_end_MissionOptions.Journeys[this->journeyIndex].Probe_End_elements_frame;
                this->probe_AEI_to_end_MissionOptions.Journeys[this->journeyIndex].AllowJourneyFreePointArrivalToPropagate = this->probe_AEI_to_end_MissionOptions.Journeys[this->journeyIndex].AllowJourneyProbeEndToPropagate;
                this->probe_AEI_to_end_MissionOptions.Journeys[this->journeyIndex].perturb_drag = this->probe_AEI_to_end_MissionOptions.Journeys[this->journeyIndex].perturb_drag_probe_AEI_to_end;
                this->probe_AEI_to_end_MissionOptions.Journeys[this->journeyIndex].spacecraft_drag_area = this->probe_AEI_to_end_MissionOptions.Journeys[this->journeyIndex].probe_drag_area_probe_AEI_to_end;
                this->probe_AEI_to_end_MissionOptions.Journeys[this->journeyIndex].coefficient_of_drag = this->probe_AEI_to_end_MissionOptions.Journeys[this->journeyIndex].probe_coefficient_of_drag_probe_AEI_to_end;
                this->probe_AEI_to_end_MissionOptions.Journeys[this->journeyIndex].final_velocity = this->probe_AEI_to_end_MissionOptions.Journeys[this->journeyIndex].probe_end_velocity;

                //initialize the boundary itself
                this->probeEndEvent = new BoundaryEvents::FreePointIntercept(name + "ProbeEndOfMissionFreePointIntercept",
                    this->journeyIndex,
                    this->phaseIndex,
                    stageIndex,
                    this->myUniverse,
                    this->mySpacecraft,
                    &this->probe_AEI_to_end_MissionOptions);

                this->probeEndEvent->construct_boundary_constraints(probe_end_arrival_constraints);
            }

            //the probe only has 8 states - position, velocity, mass, epoch. No propellant tanks. No match point on epoch.
            this->numMatchConstraints_probe = 7;

            //size state vectors
            this->probe_first_subphase_match_point_state_minus.resize(this->numStatesToPropagate, 1, 0.0);
            this->probe_first_subphase_match_point_state_plus.resize(this->numStatesToPropagate, 1, 0.0);
            this->probe_second_subphase_match_point_state_minus.resize(this->numStatesToPropagate, 1, 0.0);
            this->probe_second_subphase_match_point_state_plus.resize(this->numStatesToPropagate, 1, 0.0);
            this->spacecraft_state_after_separation.resize(this->numStatesToPropagate, 1, 0.0);
            this->probe_state_after_separation.resize(this->numStatesToPropagate, 1, 0.0);
            this->probe_state_at_AEI_arrival.resize(this->numStatesToPropagate, 1, 0.0);
            this->probe_state_at_AEI_departure.resize(this->numStatesToPropagate, 1, 0.0);
            this->probe_state_at_end_of_phase.resize(this->numStatesToPropagate, 1, 0.0);

            //separation maneuver
            this->spacecraft_separation_deltaV.resize(3, 1, 0.0);
            this->probe_separation_deltaV.resize(3, 1, 0.0);

            //size STMs
            this->TCMTM = math::Matrix<double>(this->numStatesToPropagate + 2, math::identity);
            this->probe_separation_to_AEI_ForwardSPTM = math::Matrix<double>(this->numStatesToPropagate + 2, math::identity);
            this->probe_separation_to_AEI_BackwardSPTM = math::Matrix<double>(this->numStatesToPropagate + 2, math::identity);
            this->probe_AEI_to_end_ForwardSPTM = math::Matrix<double>(this->numStatesToPropagate + 2, math::identity);
            this->probe_AEI_to_end_BackwardSPTM = math::Matrix<double>(this->numStatesToPropagate + 2, math::identity);
            this->probeReleaseTM = this->TCMTM;
            
            //size Gindex vectors for the probe
            this->G_indices_match_point_constraints_wrt_LeftBoundaryVariables_probe_separation_to_AEI.resize(this->numMatchConstraints_probe);
            this->G_indices_match_point_constraints_wrt_RightBoundaryVariables_probe_separation_to_AEI.resize(this->numMatchConstraints_probe);
            this->G_indices_match_point_constraints_wrt_LeftBoundaryTimeVariables_probe_separation_to_AEI.resize(this->numMatchConstraints_probe);
            this->G_indices_match_point_constraints_wrt_RightBoundaryTimeVariables_probe_separation_to_AEI.resize(this->numMatchConstraints_probe);
            this->G_indices_match_point_constraints_wrt_PhaseFlightTime_probe_separation_to_AEI.resize(this->numMatchConstraints_probe);
            this->G_indices_match_point_constraints_wrt_LeftBoundaryVariables_probe_AEI_to_end.resize(this->numMatchConstraints_probe);
            this->G_indices_match_point_constraints_wrt_RightBoundaryVariables_probe_AEI_to_end.resize(this->numMatchConstraints_probe);
            this->G_indices_match_point_constraints_wrt_LeftBoundaryTimeVariables_probe_AEI_to_end.resize(this->numMatchConstraints_probe);
            this->G_indices_match_point_constraints_wrt_RightBoundaryTimeVariables_probe_AEI_to_end.resize(this->numMatchConstraints_probe);
            this->G_indices_match_point_constraints_wrt_PhaseFlightTime_probe_AEI_to_end.resize(this->numMatchConstraints_probe);
            this->Gindex_separation_vector_velocity_match_constraint.resize(3);
            this->Gindex_separation_vector_velocity_match_constraint_wrt_Time.resize(3);
            this->dIndex_separation_vector_velocity_match_constraint.resize(3);
            this->dIndex_separation_vector_velocity_match_constraint_wrt_Time.resize(3);

            //redirect the first propagator for the spacecraft to point to the state after separation rather than after the initial TCM
            this->ForwardSubPhases[0].set_spacecraft_state_minus(this->spacecraft_state_after_separation);

            //configure propagators for the probe
            this->probe_separation_to_AEI_MatchPointFraction = this->myJourneyOptions->ProbeSeparationToAEI_MatchPointFraction;
            this->probe_separation_to_AEI_dForwardStepSize_dPropagationVariable = this->probe_separation_to_AEI_MatchPointFraction;
            this->probe_separation_to_AEI_dBackwardStepSize_dPropagationVariable = 1.0 - this->probe_separation_to_AEI_MatchPointFraction;
            this->probe_separation_to_AEI_ForwardIntegrationStepLength = this->myJourneyOptions->ProbeSeparationToAEI_ForwardIntegrationStepLength;
            this->probe_separation_to_AEI_BackwardIntegrationStepLength = this->myJourneyOptions->ProbeSeparationToAEI_BackwardIntegrationStepLength;

            this->probe_AEI_to_end_MatchPointFraction = this->myJourneyOptions->ProbeAEI_to_end_MatchPointFraction;
            this->probe_AEI_to_end_dForwardStepSize_dPropagationVariable = this->probe_AEI_to_end_MatchPointFraction;
            this->probe_AEI_to_end_dBackwardStepSize_dPropagationVariable = 1.0 - this->probe_AEI_to_end_MatchPointFraction;
            this->probe_AEI_to_end_ForwardIntegrationStepLength = this->myJourneyOptions->ProbeAEI_to_end_ForwardIntegrationStepLength;
            this->probe_AEI_to_end_BackwardIntegrationStepLength = this->myJourneyOptions->ProbeAEI_to_end_BackwardIntegrationStepLength;

            if (this->isKeplerian)
            {
                //create the propagators for the first subphase
                this->probe_separation_to_AEI_ForwardSTM = math::Matrix<double>(6, math::identity);
                this->probe_separation_to_AEI_BackwardSTM = math::Matrix<double>(6, math::identity);
                this->probe_separation_to_AEI_Forward_dPropagatedStatedIndependentVariable = math::Matrix<double>(6, 1, 0.0);
                this->probe_separation_to_AEI_Backward_dPropagatedStatedIndependentVariable = math::Matrix<double>(6, 1, 0.0);

                this->probe_separation_to_AEI_ForwardHalfPhasePropagator = Astrodynamics::CreatePropagator(&this->probe_separation_to_AEI_MissionOptions,
                    this->myUniverse,
                    6,
                    this->probe_state_after_separation,
                    this->probe_first_subphase_match_point_state_minus,
                    this->probe_separation_to_AEI_ForwardSTM,
                    this->probe_separation_to_AEI_Forward_dPropagatedStatedIndependentVariable,
                    &this->probe_separation_to_AEI_dForwardStepSize_dPropagationVariable);

                this->probe_separation_to_AEI_BackwardHalfPhasePropagator = Astrodynamics::CreatePropagator(&this->probe_separation_to_AEI_MissionOptions,
                    this->myUniverse,
                    6,
                    this->probe_state_at_AEI_arrival,
                    this->probe_first_subphase_match_point_state_plus,
                    this->probe_separation_to_AEI_BackwardSTM,
                    this->probe_separation_to_AEI_Backward_dPropagatedStatedIndependentVariable,
                    &this->probe_separation_to_AEI_dBackwardStepSize_dPropagationVariable);

                //create the propagators for the second subphase
                if (this->myJourneyOptions->ModelProbeSecondPhase)
                {
                    this->probe_AEI_to_end_ForwardSTM = math::Matrix<double>(6, math::identity);
                    this->probe_AEI_to_end_BackwardSTM = math::Matrix<double>(6, math::identity);
                    this->probe_AEI_to_end_Forward_dPropagatedStatedIndependentVariable = math::Matrix<double>(6, 1, 0.0);
                    this->probe_AEI_to_end_Backward_dPropagatedStatedIndependentVariable = math::Matrix<double>(6, 1, 0.0);

                    this->probe_AEI_to_end_ForwardHalfPhasePropagator = Astrodynamics::CreatePropagator(&this->probe_AEI_to_end_MissionOptions,
                        this->myUniverse,
                        6,
                        this->probe_state_at_AEI_departure,
                        this->probe_second_subphase_match_point_state_minus,
                        this->probe_AEI_to_end_ForwardSTM,
                        this->probe_AEI_to_end_Forward_dPropagatedStatedIndependentVariable,
                        &this->probe_AEI_to_end_dForwardStepSize_dPropagationVariable);

                    this->probe_AEI_to_end_BackwardHalfPhasePropagator = Astrodynamics::CreatePropagator(&this->probe_AEI_to_end_MissionOptions,
                        this->myUniverse,
                        6,
                        this->probe_state_at_end_of_phase,
                        this->probe_second_subphase_match_point_state_plus,
                        this->probe_AEI_to_end_BackwardSTM,
                        this->probe_AEI_to_end_Backward_dPropagatedStatedIndependentVariable,
                        &this->probe_AEI_to_end_dBackwardStepSize_dPropagationVariable);
                }
            }
            else //integrated propagator
            {
                //separation to AEI
                this->probe_separation_to_AEI_ForwardSTM = math::Matrix<double>(14, math::identity);
                this->probe_separation_to_AEI_BackwardSTM = math::Matrix<double>(14, math::identity);
                this->probe_separation_to_AEI_Forward_dPropagatedStatedIndependentVariable.resize(this->numStatesToPropagate, 2, 0.0);
                this->probe_separation_to_AEI_Backward_dPropagatedStatedIndependentVariable.resize(this->numStatesToPropagate, 2, 0.0);

                //acceleration model object
                this->probe_separation_to_AEI_SpacecraftAccelerationModel = new Astrodynamics::SpacecraftAccelerationModel(&this->probe_separation_to_AEI_MissionOptions,
                    &this->probe_separation_to_AEI_MissionOptions.Journeys[this->journeyIndex],
                    this->myUniverse,
                    this->Xdescriptions,
                    this->mySpacecraft,
                    11);

                this->probe_separation_to_AEI_SpacecraftAccelerationModel->setDutyCycle(this->PhaseDutyCycle);

                //EOM
                this->probe_separation_to_AEI_EOM.setSpacecraftAccelerationModel(this->probe_separation_to_AEI_SpacecraftAccelerationModel);

                //integration scheme
                this->probe_separation_to_AEI_IntegrationScheme = CreateIntegrationScheme(&this->probe_separation_to_AEI_EOM, 10, 11);

                this->probe_separation_to_AEI_ForwardHalfPhasePropagator = CreatePropagator(&this->probe_separation_to_AEI_MissionOptions,
                    this->myUniverse,
                    10,
                    11,
                    this->probe_state_after_separation,
                    this->probe_first_subphase_match_point_state_minus,
                    this->probe_separation_to_AEI_ForwardSTM,
                    this->probe_separation_to_AEI_Forward_dPropagatedStatedIndependentVariable,
                    (Integration::Integrand*) &this->probe_separation_to_AEI_EOM,
                    this->probe_separation_to_AEI_IntegrationScheme,
                    &this->probe_separation_to_AEI_dForwardStepSize_dPropagationVariable,
                    this->probe_separation_to_AEI_ForwardIntegrationStepLength);

                this->probe_separation_to_AEI_BackwardHalfPhasePropagator = CreatePropagator(&this->probe_separation_to_AEI_MissionOptions,
                    this->myUniverse,
                    10,
                    11,
                    this->probe_state_at_AEI_arrival,
                    this->probe_first_subphase_match_point_state_plus,
                    this->probe_separation_to_AEI_BackwardSTM,
                    this->probe_separation_to_AEI_Backward_dPropagatedStatedIndependentVariable,
                    (Integration::Integrand*) &this->probe_separation_to_AEI_EOM,
                    this->probe_separation_to_AEI_IntegrationScheme,
                    &this->probe_separation_to_AEI_dBackwardStepSize_dPropagationVariable,
                    this->probe_separation_to_AEI_BackwardIntegrationStepLength);

                //AEI to mission end
                if (this->myJourneyOptions->ModelProbeSecondPhase)
                {
                    this->probe_AEI_to_end_ForwardSTM = math::Matrix<double>(14, math::identity);
                    this->probe_AEI_to_end_BackwardSTM = math::Matrix<double>(14, math::identity);
                    this->probe_AEI_to_end_Forward_dPropagatedStatedIndependentVariable.resize(this->numStatesToPropagate, 2, 0.0);
                    this->probe_AEI_to_end_Backward_dPropagatedStatedIndependentVariable.resize(this->numStatesToPropagate, 2, 0.0);

                    //acceleration model object
                    this->probe_AEI_to_end_SpacecraftAccelerationModel = new Astrodynamics::SpacecraftAccelerationModel(&this->probe_separation_to_AEI_MissionOptions,
                        &this->probe_AEI_to_end_MissionOptions.Journeys[this->journeyIndex],
                        this->myUniverse,
                        this->Xdescriptions,
                        this->mySpacecraft,
                        10);
                    this->probe_AEI_to_end_SpacecraftAccelerationModel->setDutyCycle(this->PhaseDutyCycle);

                    //EOM
                    this->probe_AEI_to_end_EOM.setSpacecraftAccelerationModel(this->probe_AEI_to_end_SpacecraftAccelerationModel);

                    //integration scheme
                    this->probe_AEI_to_end_IntegrationScheme = CreateIntegrationScheme(&this->probe_AEI_to_end_EOM, 9, 10);

                    this->probe_AEI_to_end_ForwardHalfPhasePropagator = CreatePropagator(&this->probe_AEI_to_end_MissionOptions,
                        this->myUniverse,
                        9,
                        10,
                        this->probe_state_at_AEI_departure,
                        this->probe_second_subphase_match_point_state_minus,
                        this->probe_AEI_to_end_ForwardSTM,
                        this->probe_AEI_to_end_Forward_dPropagatedStatedIndependentVariable,
                        (Integration::Integrand*) &this->probe_AEI_to_end_EOM,
                        this->probe_AEI_to_end_IntegrationScheme,
                        &this->probe_AEI_to_end_dForwardStepSize_dPropagationVariable,
                        this->probe_AEI_to_end_ForwardIntegrationStepLength);

                    this->probe_AEI_to_end_BackwardHalfPhasePropagator = CreatePropagator(&this->probe_AEI_to_end_MissionOptions,
                        this->myUniverse,
                        9,
                        10,
                        this->probe_state_at_end_of_phase,
                        this->probe_second_subphase_match_point_state_plus,
                        this->probe_AEI_to_end_BackwardSTM,
                        this->probe_AEI_to_end_Backward_dPropagatedStatedIndependentVariable,
                        (Integration::Integrand*) &this->probe_AEI_to_end_EOM,
                        this->probe_AEI_to_end_IntegrationScheme,
                        &this->probe_AEI_to_end_dBackwardStepSize_dPropagationVariable,
                        this->probe_AEI_to_end_BackwardIntegrationStepLength);
                }
            }//done creating propagators

            //output stuff
            this->output_state.resize(this->numStatesToPropagate + 13 * 13, 1, 0.0);
        }//end constructor

        //destructor
        ProbeEntryPhase::~ProbeEntryPhase()
        {

            delete this->probeEntryArrivalEvent;

            delete this->probe_separation_to_AEI_ForwardHalfPhasePropagator;
            delete this->probe_separation_to_AEI_BackwardHalfPhasePropagator;

            if (this->myPropagatorType == PropagatorType::IntegratedPropagator)
            {
                delete this->probe_separation_to_AEI_IntegrationScheme;
                delete this->probe_separation_to_AEI_SpacecraftAccelerationModel;
            }


            if (this->myJourneyOptions->ModelProbeSecondPhase)
            {
                delete this->probeEntryDepartureEvent;
                delete this->probeEndEvent;

                delete this->probe_AEI_to_end_ForwardHalfPhasePropagator;
                delete this->probe_AEI_to_end_BackwardHalfPhasePropagator;

                if (this->myPropagatorType == PropagatorType::IntegratedPropagator)
                {
                    delete this->probe_AEI_to_end_IntegrationScheme;
                    delete this->probe_AEI_to_end_SpacecraftAccelerationModel;
                }
            }
        }//end destructor

        //******************************************calcbounds methods
        void ProbeEntryPhase::setup_calcbounds(std::vector<double>* Xupperbounds,
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
            std::vector<double>* A)
        {
            //base class
            this->MGAnDSMs_phase::setup_calcbounds(Xupperbounds,
                Xlowerbounds,
                X_scale_factors,
                Fupperbounds,
                Flowerbounds,
                F_scale_factors,
                Xdescriptions,
                Fdescriptions,
                iGfun,
                jGvar,
                Gdescriptions,
                iAfun,
                jAvar,
                Adescriptions,
                A);

            //probe entry arrival event
            this->probeEntryArrivalEvent->setup_calcbounds(Xupperbounds,
                Xlowerbounds,
                X_scale_factors,
                Fupperbounds,
                Flowerbounds,
                F_scale_factors,
                Xdescriptions,
                Fdescriptions,
                iGfun,
                jGvar,
                Gdescriptions,
                iAfun,
                jAvar,
                Adescriptions,
                A);

            //if modeling the second subphase
            if (this->myJourneyOptions->ModelProbeSecondPhase)
            {
                //probe entry departure event
                this->probeEntryDepartureEvent->setup_calcbounds(Xupperbounds,
                    Xlowerbounds,
                    X_scale_factors,
                    Fupperbounds,
                    Flowerbounds,
                    F_scale_factors,
                    Xdescriptions,
                    Fdescriptions,
                    iGfun,
                    jGvar,
                    Gdescriptions,
                    iAfun,
                    jAvar,
                    Adescriptions,
                    A);

                //probe end event
                this->probeEndEvent->setup_calcbounds(Xupperbounds,
                    Xlowerbounds,
                    X_scale_factors,
                    Fupperbounds,
                    Flowerbounds,
                    F_scale_factors,
                    Xdescriptions,
                    Fdescriptions,
                    iGfun,
                    jGvar,
                    Gdescriptions,
                    iAfun,
                    jAvar,
                    Adescriptions,
                    A);
            }
        }//end setup_calcbounds()

        void ProbeEntryPhase::calcbounds()
        {
            this->calcbounds_phase_left_boundary();

            this->calcbounds_phase_flight_time();

            this->calcbounds_phase_right_boundary();

            //if modeling the second subphase
            if (this->myJourneyOptions->ModelProbeSecondPhase)
            {
                this->calcbounds_probe_descent_time();
            }

            this->calcbounds_probe_boundary_events();

            this->calcbounds_phase_main();
        }//end calcbounds

        void ProbeEntryPhase::calcbounds_probe_descent_time()
        {
            //for now we'll give it the same bounds as the phase flight time
            //this->Xupperbounds->push_back(this->Xupperbounds->operator[](this->Xindex_PhaseFlightTime));
            //this->Xlowerbounds->push_back(this->Xlowerbounds->operator[](this->Xindex_PhaseFlightTime));
            this->Xupperbounds->push_back(86400.0);
            this->Xlowerbounds->push_back(0.0);
            this->Xdescriptions->push_back(this->prefix + "probe descent time");
            this->X_scale_factors->push_back(this->myUniverse->TU);
            this->Xindex_probeDescentTime = this->Xdescriptions->size() - 1;
        }//end calcbounds_probe_descent_time

        void ProbeEntryPhase::calcbounds_probe_boundary_events()
        {
            //The purpose of this method is to call calcbounds on that probe's boundary events

            //AEI event has the same dependencies as the carrier's arrival event
            std::vector<size_t> timeVariables = this->myArrivalEvent->get_Xindices_EventRightEpoch();
            this->probeEntryArrivalEvent->calcbounds(timeVariables);


            //if modeling the second subphase
            if (this->myJourneyOptions->ModelProbeSecondPhase)
            {
                //probe entry departure depends on probe entry arrival
                this->probeEntryDepartureEvent->calcbounds(this->probeEntryArrivalEvent->get_Xindices_EventRightEpoch());

                //probe end event is also dependent on probe descent time
                timeVariables.push_back(this->Xindex_probeDescentTime);

                this->probeEndEvent->calcbounds(timeVariables);
            }
        }//end calcbounds_probe_boundary_events()

        void ProbeEntryPhase::calcbounds_phase_main()
        {
            //Step 1: separation impulse
            this->calcbounds_separation();

            //Step 2: call calcbounds() on the base MGAnDSMs phase
            this->MGAnDSMs_phase::calcbounds_phase_main();

            //Step 3: additional match point derivative entries for the spacecraft, due to the separation impulse
            this->calcbounds_spacecraft_match_point_derivatives_due_to_probe_release();

            //Step 4: probe match points from separation to AEI
            this->calcbounds_probe_match_point_constraints_from_separation_to_AEI();

            //Step 5: communication constraint
            this->calcbounds_communication_distance_constraint();

            //Step 6: probe match points from AEI to end
            if (this->myJourneyOptions->ModelProbeSecondPhase)
            {
                this->calcbounds_probe_match_point_constraints_from_AEI_to_end();
            }
        }//end calcbounds_phase_main()

        void ProbeEntryPhase::calcbounds_spacecraft_match_point_derivatives_due_to_probe_release()
        {
            //here we account for spacecraft match point entries due to the probe separation
            //there will be derivatives with respect to the separation impulse unit vector
            //and also due to spacecraft initial mass (although this should already be present?)

            //Step 1: derivatives due to probe separation impulse
            for (size_t constraintIndex = 0; constraintIndex < this->numMatchConstraints; ++constraintIndex)
            {
                std::vector<size_t> stateGindex_spacecraft_match_point_constraints_wrt_separation_impulse;

                for (size_t velocityIndex = 0; velocityIndex < 3; ++velocityIndex)
                {
                    size_t Xindex = this->Xindices_separation_direction_vector[velocityIndex];

                    this->create_sparsity_entry(this->Findices_match_point_constraints[constraintIndex],
                        Xindex,
                        stateGindex_spacecraft_match_point_constraints_wrt_separation_impulse);
                }
                this->Gindex_spacecraft_match_point_constraints_wrt_separation_impulse.push_back(stateGindex_spacecraft_match_point_constraints_wrt_separation_impulse);
            }//end loop over match point constraints

            //Step 2: track derivatives wrt left mass
            std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateBeforePhase = this->myDepartureEvent->get_Derivatives_of_StateAfterEvent();//Xindex, stateIndex, derivative value
            
            for (size_t dIndex = 0; dIndex < Derivatives_of_StateBeforePhase.size(); ++dIndex)
            {
                size_t stateIndex = std::get<1>(Derivatives_of_StateBeforePhase[dIndex]);

                if (stateIndex == 6)//this is a mass variable
                {
                    this->dIndex_spacecraft_left_mass.push_back(dIndex);
                    size_t Xindex = std::get<0>(Derivatives_of_StateBeforePhase[dIndex]);

                    for (size_t constraintIndex = 0; constraintIndex < this->numMatchConstraints; ++constraintIndex)
                    {
                        //these entries should already exist, so all we're doing here is storing their Gindex
                        this->create_sparsity_entry(this->Findices_match_point_constraints[constraintIndex],
                            Xindex,
                            this->Gindex_spacecraft_match_point_constraints_wrt_left_mass);
                    }//end loop over match point constraints
                }
            }
        }//end calcbounds_spacecraft_match_point_derivatives_due_to_probe_release()

        void ProbeEntryPhase::calcbounds_separation()
        {
            std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_probe_StateAtEnd = this->probeEndEvent->get_Derivatives_of_StateBeforeEvent();//Xindex, stateIndex, derivative value
            std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_probe_StateAtEnd_wrt_Time = this->probeEndEvent->get_Derivatives_of_StateBeforeEvent_wrt_Time();//Xindex, stateIndex, derivative value
            for (size_t velocityIndex : {0, 1, 2})
            {
                //Step 1: variables defining the separation impulse direction
                this->Xlowerbounds->push_back(-1.0);
                this->Xupperbounds->push_back(1.0);
                this->Xdescriptions->push_back(this->prefix + "separation unit vector " + stateVectorNames[velocityIndex]);
                this->X_scale_factors->push_back(1.0);
                this->Xindices_separation_direction_vector.push_back(this->Xdescriptions->size() - 1);

                //Step 2: constraint to match the inertial velocity of the probe at entry
                this->Flowerbounds->push_back(0.0);
                this->Fupperbounds->push_back(0.0);
                this->Fdescriptions->push_back(this->prefix + "separation direction vector match probe velocity unit vector " + stateVectorNames[velocityIndex]);

                //Step 3: derivative w.r.t. this velocity component
                this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                    this->Xindices_separation_direction_vector[velocityIndex],
                    this->Gindices_separation_direction_vector_wrt_self);

                //Step 4: derivatives of velocity match w.r.t. non-time variables
                for (size_t dIndex = 0; dIndex < Derivatives_of_probe_StateAtEnd.size(); ++dIndex)
                {
                    size_t Xindex = std::get<0>(Derivatives_of_probe_StateAtEnd[dIndex]);
                    size_t stateIndex = std::get<1>(Derivatives_of_probe_StateAtEnd[dIndex]);

                    if (stateIndex > 2 && stateIndex < 6)//i.e. does this variable influence velocity
                    {
                        this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                            Xindex,
                            this->Gindex_separation_vector_velocity_match_constraint[velocityIndex]);

                        this->dIndex_separation_vector_velocity_match_constraint[velocityIndex].push_back(dIndex);
                    }
                }//end non-time Jacobian sparsity pattern

                //Step 5: derivatives of velocity match w.r.t. time variables
                for (size_t dIndex = 0; dIndex < Derivatives_of_probe_StateAtEnd_wrt_Time.size(); ++dIndex)
                {
                    size_t Xindex = std::get<0>(Derivatives_of_probe_StateAtEnd_wrt_Time[dIndex]);
                    size_t stateIndex = std::get<1>(Derivatives_of_probe_StateAtEnd_wrt_Time[dIndex]);

                    if (stateIndex > 2 && stateIndex < 6)//i.e. does this variable influence velocity
                    {

                        this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                            Xindex,
                            this->Gindex_separation_vector_velocity_match_constraint_wrt_Time[velocityIndex]);

                        this->dIndex_separation_vector_velocity_match_constraint_wrt_Time[velocityIndex].push_back(dIndex);
                    }
                }//end time Jacobian sparsity pattern
            }//end loop over separation unit vector components
        }//end calcbounds_separation()

        void ProbeEntryPhase::calcbounds_probe_match_point_constraints_from_separation_to_AEI()
        {
            std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateBeforePhase = this->myDepartureEvent->get_Derivatives_of_StateAfterEvent();//Xindex, stateIndex, derivative value
            std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_probe_StateAtAEI = this->probeEntryArrivalEvent->get_Derivatives_of_StateBeforeEvent();//Xindex, stateIndex, derivative value
            std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateBeforePhase_wrt_Time = this->myDepartureEvent->get_Derivatives_of_StateAfterEvent_wrt_Time();//Xindex, stateIndex, derivative value
            std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_probe_StateAtAEI_wrt_Time = this->probeEntryArrivalEvent->get_Derivatives_of_StateBeforeEvent_wrt_Time();//Xindex, stateIndex, derivative value

            //Step 1: make a list of variables affecting the right boundary of the subphase from separation to AEI
            {
                //Step 1.1: The left boundary is the same as for the spacecraft so we don't need to do anything new
                this->DerivativesOfLeftBoundaryByVariable_probe_separation_to_AEI = this->DerivativesOfLeftBoundaryByVariable;

                //Step 1.2: right boundary. This is special because the probe has its own arrival event
                for (size_t dIndex = 0; dIndex < Derivatives_of_probe_StateAtAEI.size(); ++dIndex)
                {
                    bool alreadyListed = false;
                    size_t Xindex = std::get<0>(Derivatives_of_probe_StateAtAEI[dIndex]);
                    for (size_t listIndex = 0; listIndex < this->ListOfVariablesAffectingRightBoundary_probe_separation_to_AEI.size(); ++listIndex)
                    {
                        if (this->ListOfVariablesAffectingRightBoundary_probe_separation_to_AEI[listIndex] == Xindex)
                        {
                            alreadyListed = true;
                            break;
                        }
                    }
                    if (!alreadyListed)
                    {
                        this->ListOfVariablesAffectingRightBoundary_probe_separation_to_AEI.push_back(Xindex);
                    }
                }

                this->DerivativesOfRightBoundaryByVariable_probe_separation_to_AEI.resize(ListOfVariablesAffectingRightBoundary_probe_separation_to_AEI.size());

                for (size_t listIndex = 0; listIndex < this->ListOfVariablesAffectingRightBoundary_probe_separation_to_AEI.size(); ++listIndex)
                {
                    for (size_t dIndex = 0; dIndex < Derivatives_of_probe_StateAtAEI.size(); ++dIndex)
                    {
                        size_t Xindex = std::get<0>(Derivatives_of_probe_StateAtAEI[dIndex]);

                        if (Xindex == this->ListOfVariablesAffectingRightBoundary_probe_separation_to_AEI[listIndex])
                        {
                            this->DerivativesOfRightBoundaryByVariable_probe_separation_to_AEI[listIndex].push_back(dIndex);
                        }
                    }
                }
            }//end non-time variables

            //Step 2: make a list of time variables affecting both boundaries of the subphase from separation to AEI
            {
                //Step 2.1: Left Boundary is the same as the spacecraft
                this->DerivativesOfLeftBoundaryByTimeVariable_probe_separation_to_AEI = this->DerivativesOfLeftBoundaryByTimeVariable;
                
                //Step 2.2: Right Boundary
                for (size_t dIndex = 0; dIndex < Derivatives_of_probe_StateAtAEI_wrt_Time.size(); ++dIndex)
                {
                    bool alreadyListed = false;
                    size_t Xindex = std::get<0>(Derivatives_of_probe_StateAtAEI_wrt_Time[dIndex]);
                    for (size_t listIndex = 0; listIndex < this->ListOfTimeVariablesAffectingRightBoundary_probe_separation_to_AEI.size(); ++listIndex)
                    {
                        if (this->ListOfTimeVariablesAffectingRightBoundary_probe_separation_to_AEI[listIndex] == Xindex)
                        {
                            alreadyListed = true;
                            break;
                        }
                    }
                    if (!alreadyListed)
                    {
                        this->ListOfTimeVariablesAffectingRightBoundary_probe_separation_to_AEI.push_back(Xindex);
                    }
                }

                this->DerivativesOfRightBoundaryByTimeVariable_probe_separation_to_AEI.resize(ListOfTimeVariablesAffectingRightBoundary_probe_separation_to_AEI.size());

                for (size_t listIndex = 0; listIndex < this->ListOfTimeVariablesAffectingRightBoundary_probe_separation_to_AEI.size(); ++listIndex)
                {
                    for (size_t dIndex = 0; dIndex < Derivatives_of_probe_StateAtAEI_wrt_Time.size(); ++dIndex)
                    {
                        size_t Xindex = std::get<0>(Derivatives_of_probe_StateAtAEI_wrt_Time[dIndex]);

                        if (Xindex == this->ListOfTimeVariablesAffectingRightBoundary_probe_separation_to_AEI[listIndex])
                        {
                            this->DerivativesOfRightBoundaryByTimeVariable_probe_separation_to_AEI[listIndex].push_back(dIndex);
                        }
                    }
                }
            }

            //Step 3: construct the match point constraints for the probe from separation to AEI
            for (size_t constraintIndex = 0; constraintIndex < this->numMatchConstraints_probe; ++constraintIndex)
            {
                this->Flowerbounds->push_back(-math::SMALL);
                this->Fupperbounds->push_back(math::SMALL);
                this->Fdescriptions->push_back(prefix + "probe separation to AEI match point " + this->matchPointConstraintNames[constraintIndex]);
                this->Findices_match_point_constraints_probe_separation_to_AEI.push_back(Fdescriptions->size() - 1);

                //derivatives with respect to phase flight time
                if (this->TruthTable_MatchConstraints_Derivative_wrt_PropagationTime[constraintIndex])
                {
                    this->create_sparsity_entry(this->Findices_match_point_constraints_probe_separation_to_AEI[constraintIndex],
                        this->Xindex_PhaseFlightTime,
                        this->G_indices_match_point_constraints_wrt_PhaseFlightTime_probe_separation_to_AEI[constraintIndex]);
                }


                //derivatives with respect to decision variables that affect boundary states
                //with respect to left boundary
                std::vector<bool> constraint_LeftBoundaryVariableAffectsMatchConstraint(this->ListOfVariablesAffectingLeftBoundary.size(), false);
                for (size_t varIndex = 0; varIndex < this->ListOfVariablesAffectingLeftBoundary.size(); ++varIndex)
                {
                    size_t Xindex = this->ListOfVariablesAffectingLeftBoundary[varIndex];

                    for (size_t entryIndex = 0; entryIndex < this->DerivativesOfLeftBoundaryByVariable_probe_separation_to_AEI[varIndex].size(); ++entryIndex)
                    {
                        size_t dIndex = this->DerivativesOfLeftBoundaryByVariable_probe_separation_to_AEI[varIndex][entryIndex];
                        size_t stateIndex = std::get<1>(Derivatives_of_StateBeforePhase[dIndex]);

                        //hooray, this variable affects the constraint. Make the Jacobian entry!
                        if (this->TruthTable_Forward_MatchConstraints_Derivative_wrt_EncodedStates[constraintIndex][stateIndex])
                        {
                            constraint_LeftBoundaryVariableAffectsMatchConstraint[varIndex] = true;
                            break;
                        }
                    }

                    if (constraint_LeftBoundaryVariableAffectsMatchConstraint[varIndex])
                    {
                        this->create_sparsity_entry(this->Findices_match_point_constraints_probe_separation_to_AEI[constraintIndex],
                            Xindex,
                            this->G_indices_match_point_constraints_wrt_LeftBoundaryVariables_probe_separation_to_AEI[constraintIndex]);
                    }
                    else
                        this->G_indices_match_point_constraints_wrt_LeftBoundaryVariables_probe_separation_to_AEI[constraintIndex].push_back(0);
                }
                this->LeftBoundaryVariableAffectsMatchConstraint_probe_separation_to_AEI.push_back(constraint_LeftBoundaryVariableAffectsMatchConstraint);

                //with respect to right boundary
                std::vector<bool> constraint_RightBoundaryVariableAffectsMatchConstraint(this->ListOfVariablesAffectingRightBoundary_probe_separation_to_AEI.size(), false);
                for (size_t varIndex = 0; varIndex < this->ListOfVariablesAffectingRightBoundary_probe_separation_to_AEI.size(); ++varIndex)
                {
                    size_t Xindex = this->ListOfVariablesAffectingRightBoundary_probe_separation_to_AEI[varIndex];

                    for (size_t entryIndex = 0; entryIndex < this->DerivativesOfRightBoundaryByVariable_probe_separation_to_AEI[varIndex].size(); ++entryIndex)
                    {
                        size_t dIndex = this->DerivativesOfRightBoundaryByVariable_probe_separation_to_AEI[varIndex][entryIndex];
                        size_t stateIndex = std::get<1>(Derivatives_of_probe_StateAtAEI[dIndex]);

                        //hooray, this variable affects the constraint. Make the Jacobian entry!
                        if (this->TruthTable_Backward_MatchConstraints_Derivative_wrt_EncodedStates[constraintIndex][stateIndex])
                        {
                            constraint_RightBoundaryVariableAffectsMatchConstraint[varIndex] = true;
                            break;
                        }
                    }

                    if (constraint_RightBoundaryVariableAffectsMatchConstraint[varIndex])
                    {
                        this->create_sparsity_entry(this->Findices_match_point_constraints_probe_separation_to_AEI[constraintIndex],
                            Xindex,
                            this->G_indices_match_point_constraints_wrt_RightBoundaryVariables_probe_separation_to_AEI[constraintIndex]);
                    }
                    else
                        this->G_indices_match_point_constraints_wrt_RightBoundaryVariables_probe_separation_to_AEI[constraintIndex].push_back(0);
                }
                this->RightBoundaryVariableAffectsMatchConstraint_probe_separation_to_AEI.push_back(constraint_RightBoundaryVariableAffectsMatchConstraint);

                //derivatives with respect to time variables that affect boundary states
                //with respect to left boundary
                std::vector<bool> constraint_LeftBoundaryTimeVariableAffectsMatchConstraint(this->ListOfTimeVariablesAffectingLeftBoundary.size(), false);
                for (size_t varIndex = 0; varIndex < this->ListOfTimeVariablesAffectingLeftBoundary.size(); ++varIndex)
                {
                    size_t Xindex = this->ListOfTimeVariablesAffectingLeftBoundary[varIndex];

                    for (size_t entryIndex = 0; entryIndex < this->DerivativesOfLeftBoundaryByTimeVariable_probe_separation_to_AEI[varIndex].size(); ++entryIndex)
                    {
                        size_t dIndex = this->DerivativesOfLeftBoundaryByTimeVariable_probe_separation_to_AEI[varIndex][entryIndex];
                        size_t stateIndex = std::get<1>(Derivatives_of_StateBeforePhase_wrt_Time[dIndex]);

                        //hooray, this variable affects the constraint. Make the Jacobian entry!
                        if (this->TruthTable_Forward_MatchConstraints_Derivative_wrt_EncodedStates[constraintIndex][stateIndex])
                        {
                            constraint_LeftBoundaryTimeVariableAffectsMatchConstraint[varIndex] = true;
                            break;
                        }
                    }

                    if (constraint_LeftBoundaryTimeVariableAffectsMatchConstraint[varIndex])
                    {
                        this->create_sparsity_entry(this->Findices_match_point_constraints_probe_separation_to_AEI[constraintIndex],
                            Xindex,
                            this->G_indices_match_point_constraints_wrt_LeftBoundaryTimeVariables_probe_separation_to_AEI[constraintIndex]);
                    }
                    else
                        this->G_indices_match_point_constraints_wrt_LeftBoundaryVariables_probe_separation_to_AEI[constraintIndex].push_back(0);
                }
                this->LeftBoundaryTimeVariableAffectsMatchConstraint_probe_separation_to_AEI.push_back(constraint_LeftBoundaryTimeVariableAffectsMatchConstraint);

                //with respect to right boundary
                std::vector<bool> constraint_RightBoundaryTimeVariableAffectsMatchConstraint(this->ListOfTimeVariablesAffectingRightBoundary_probe_separation_to_AEI.size(), false);
                for (size_t varIndex = 0; varIndex < this->ListOfTimeVariablesAffectingRightBoundary_probe_separation_to_AEI.size(); ++varIndex)
                {
                    size_t Xindex = this->ListOfTimeVariablesAffectingRightBoundary_probe_separation_to_AEI[varIndex];

                    for (size_t entryIndex = 0; entryIndex < this->DerivativesOfRightBoundaryByTimeVariable_probe_separation_to_AEI[varIndex].size(); ++entryIndex)
                    {
                        size_t dIndex = this->DerivativesOfRightBoundaryByTimeVariable_probe_separation_to_AEI[varIndex][entryIndex];
                        size_t stateIndex = std::get<1>(Derivatives_of_probe_StateAtAEI_wrt_Time[dIndex]);

                        //hooray, this variable affects the constraint. Make the Jacobian entry!
                        if (this->TruthTable_Backward_MatchConstraints_Derivative_wrt_EncodedStates[constraintIndex][stateIndex])
                        {
                            constraint_RightBoundaryTimeVariableAffectsMatchConstraint[varIndex] = true;
                            break;
                        }
                    }

                    if (constraint_RightBoundaryTimeVariableAffectsMatchConstraint[varIndex])
                    {
                        this->create_sparsity_entry(this->Findices_match_point_constraints_probe_separation_to_AEI[constraintIndex],
                            Xindex,
                            this->G_indices_match_point_constraints_wrt_RightBoundaryTimeVariables_probe_separation_to_AEI[constraintIndex]);
                    }
                    else
                        this->G_indices_match_point_constraints_wrt_RightBoundaryVariables_probe_separation_to_AEI[constraintIndex].push_back(0);
                }
                this->RightBoundaryTimeVariableAffectsMatchConstraint_probe_separation_to_AEI.push_back(constraint_RightBoundaryTimeVariableAffectsMatchConstraint);
            }

            //Step 4: derivatives of the match point from separation to AEI due to probe separation impulse
            for (size_t constraintIndex = 0; constraintIndex < this->numMatchConstraints_probe; ++constraintIndex)
            {
                std::vector<size_t> stateGindex_probe_match_point_constraints_wrt_separation_impulse;

                for (size_t velocityIndex = 0; velocityIndex < 3; ++velocityIndex)
                {
                    size_t Xindex = this->Xindices_separation_direction_vector[velocityIndex];

                    this->create_sparsity_entry(this->Findices_match_point_constraints_probe_separation_to_AEI[constraintIndex],
                        Xindex,
                        stateGindex_probe_match_point_constraints_wrt_separation_impulse);
                }
                this->Gindex_probe_match_point_constraints_wrt_separation_impulse.push_back(stateGindex_probe_match_point_constraints_wrt_separation_impulse);
            }
        }//end calcbounds_probe_match_point_constraints_from_separation_to_AEI()

        void ProbeEntryPhase::calcbounds_probe_match_point_constraints_from_AEI_to_end()
        {
            std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_probe_StateAfterAEI = this->probeEntryDepartureEvent->get_Derivatives_of_StateAfterEvent();//Xindex, stateIndex, derivative value
            std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_probe_StateAtEnd = this->probeEndEvent->get_Derivatives_of_StateBeforeEvent();//Xindex, stateIndex, derivative value
            std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateAfterAEI_wrt_Time = this->probeEntryDepartureEvent->get_Derivatives_of_StateAfterEvent_wrt_Time();//Xindex, stateIndex, derivative value
            std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_probe_StateAtEnd_wrt_Time = this->probeEndEvent->get_Derivatives_of_StateBeforeEvent_wrt_Time();//Xindex, stateIndex, derivative value

            //Step 1: make a list of variables affecting the right boundary of the subphase from separation to AEI
            {
                //Step 1.1: The left boundary is drawn from the second subphase's departure event
                for (size_t dIndex = 0; dIndex < Derivatives_of_probe_StateAfterAEI.size(); ++dIndex)
                {
                    bool alreadyListed = false;
                    size_t Xindex = std::get<0>(Derivatives_of_probe_StateAtEnd[dIndex]);
                    for (size_t listIndex = 0; listIndex < this->ListOfVariablesAffectingLeftBoundary_probe_AEI_to_end.size(); ++listIndex)
                    {
                        if (this->ListOfVariablesAffectingLeftBoundary_probe_AEI_to_end[listIndex] == Xindex)
                        {
                            alreadyListed = true;
                            break;
                        }
                    }
                    if (!alreadyListed)
                    {
                        this->ListOfVariablesAffectingLeftBoundary_probe_AEI_to_end.push_back(Xindex);
                    }
                }

                this->DerivativesOfLeftBoundaryByVariable_probe_AEI_to_end.resize(ListOfVariablesAffectingLeftBoundary_probe_AEI_to_end.size());

                for (size_t listIndex = 0; listIndex < this->ListOfVariablesAffectingLeftBoundary_probe_AEI_to_end.size(); ++listIndex)
                {
                    for (size_t dIndex = 0; dIndex < Derivatives_of_probe_StateAfterAEI.size(); ++dIndex)
                    {
                        size_t Xindex = std::get<0>(Derivatives_of_probe_StateAfterAEI[dIndex]);

                        if (Xindex == this->ListOfVariablesAffectingLeftBoundary_probe_AEI_to_end[listIndex])
                        {
                            this->DerivativesOfLeftBoundaryByVariable_probe_AEI_to_end[listIndex].push_back(dIndex);
                        }
                    }
                }

                //Step 1.2: right boundary is drawn from the second subphase's arrival event
                for (size_t dIndex = 0; dIndex < Derivatives_of_probe_StateAtEnd.size(); ++dIndex)
                {
                    bool alreadyListed = false;
                    size_t Xindex = std::get<0>(Derivatives_of_probe_StateAtEnd[dIndex]);
                    for (size_t listIndex = 0; listIndex < this->ListOfVariablesAffectingRightBoundary_probe_AEI_to_end.size(); ++listIndex)
                    {
                        if (this->ListOfVariablesAffectingRightBoundary_probe_AEI_to_end[listIndex] == Xindex)
                        {
                            alreadyListed = true;
                            break;
                        }
                    }
                    if (!alreadyListed)
                    {
                        this->ListOfVariablesAffectingRightBoundary_probe_AEI_to_end.push_back(Xindex);
                    }
                }

                this->DerivativesOfRightBoundaryByVariable_probe_AEI_to_end.resize(ListOfVariablesAffectingRightBoundary_probe_AEI_to_end.size());

                for (size_t listIndex = 0; listIndex < this->ListOfVariablesAffectingRightBoundary_probe_AEI_to_end.size(); ++listIndex)
                {
                    for (size_t dIndex = 0; dIndex < Derivatives_of_probe_StateAtEnd.size(); ++dIndex)
                    {
                        size_t Xindex = std::get<0>(Derivatives_of_probe_StateAtEnd[dIndex]);

                        if (Xindex == this->ListOfVariablesAffectingRightBoundary_probe_AEI_to_end[listIndex])
                        {
                            this->DerivativesOfRightBoundaryByVariable_probe_AEI_to_end[listIndex].push_back(dIndex);
                        }
                    }
                }
            }//end non-time variables

            //Step 2: make a list of time variables affecting both boundaries of the subphase from separation to AEI
            {
                //Step 2.1: The left boundary is drawn from the second subphase's departure event
                for (size_t dIndex = 0; dIndex < Derivatives_of_StateAfterAEI_wrt_Time.size(); ++dIndex)
                {
                    bool alreadyListed = false;
                    size_t Xindex = std::get<0>(Derivatives_of_probe_StateAtEnd_wrt_Time[dIndex]);
                    for (size_t listIndex = 0; listIndex < this->ListOfTimeVariablesAffectingLeftBoundary_probe_AEI_to_end.size(); ++listIndex)
                    {
                        if (this->ListOfTimeVariablesAffectingLeftBoundary_probe_AEI_to_end[listIndex] == Xindex)
                        {
                            alreadyListed = true;
                            break;
                        }
                    }
                    if (!alreadyListed)
                    {
                        this->ListOfTimeVariablesAffectingLeftBoundary_probe_AEI_to_end.push_back(Xindex);
                    }
                }

                this->DerivativesOfLeftBoundaryByTimeVariable_probe_AEI_to_end.resize(ListOfTimeVariablesAffectingLeftBoundary_probe_AEI_to_end.size());

                for (size_t listIndex = 0; listIndex < this->ListOfTimeVariablesAffectingLeftBoundary_probe_AEI_to_end.size(); ++listIndex)
                {
                    for (size_t dIndex = 0; dIndex < Derivatives_of_StateAfterAEI_wrt_Time.size(); ++dIndex)
                    {
                        size_t Xindex = std::get<0>(Derivatives_of_StateAfterAEI_wrt_Time[dIndex]);

                        if (Xindex == this->ListOfTimeVariablesAffectingLeftBoundary_probe_AEI_to_end[listIndex])
                        {
                            this->DerivativesOfLeftBoundaryByTimeVariable_probe_AEI_to_end[listIndex].push_back(dIndex);
                        }
                    }
                }

                //Step 2.2: right boundary is drawn from the second subphase's arrival event
                for (size_t dIndex = 0; dIndex < Derivatives_of_probe_StateAtEnd_wrt_Time.size(); ++dIndex)
                {
                    bool alreadyListed = false;
                    size_t Xindex = std::get<0>(Derivatives_of_probe_StateAtEnd_wrt_Time[dIndex]);
                    for (size_t listIndex = 0; listIndex < this->ListOfTimeVariablesAffectingRightBoundary_probe_AEI_to_end.size(); ++listIndex)
                    {
                        if (this->ListOfTimeVariablesAffectingRightBoundary_probe_AEI_to_end[listIndex] == Xindex)
                        {
                            alreadyListed = true;
                            break;
                        }
                    }
                    if (!alreadyListed)
                    {
                        this->ListOfTimeVariablesAffectingRightBoundary_probe_AEI_to_end.push_back(Xindex);
                    }
                }

                this->DerivativesOfRightBoundaryByTimeVariable_probe_AEI_to_end.resize(ListOfTimeVariablesAffectingRightBoundary_probe_AEI_to_end.size());

                for (size_t listIndex = 0; listIndex < this->ListOfTimeVariablesAffectingRightBoundary_probe_AEI_to_end.size(); ++listIndex)
                {
                    for (size_t dIndex = 0; dIndex < Derivatives_of_probe_StateAtEnd_wrt_Time.size(); ++dIndex)
                    {
                        size_t Xindex = std::get<0>(Derivatives_of_probe_StateAtEnd_wrt_Time[dIndex]);

                        if (Xindex == this->ListOfTimeVariablesAffectingRightBoundary_probe_AEI_to_end[listIndex])
                        {
                            this->DerivativesOfRightBoundaryByTimeVariable_probe_AEI_to_end[listIndex].push_back(dIndex);
                        }
                    }
                }
            }

            //Step 3: construct the match point constraints for the probe from AEI to probe mission end
            for (size_t constraintIndex = 0; constraintIndex < this->numMatchConstraints_probe; ++constraintIndex)
            {
                this->Flowerbounds->push_back(-math::SMALL);
                this->Fupperbounds->push_back(math::SMALL);
                this->Fdescriptions->push_back(prefix + "probe AEI to end match point " + this->matchPointConstraintNames[constraintIndex]);
                this->Findices_match_point_constraints_probe_AEI_to_end.push_back(Fdescriptions->size() - 1);

                //derivatives with respect to phase flight time
                if (this->TruthTable_MatchConstraints_Derivative_wrt_PropagationTime[constraintIndex])
                {
                    this->create_sparsity_entry(this->Findices_match_point_constraints_probe_AEI_to_end[constraintIndex],
                        this->Xindex_probeDescentTime,
                        this->G_indices_match_point_constraints_wrt_PhaseFlightTime_probe_AEI_to_end[constraintIndex]);
                }


                //derivatives with respect to decision variables that affect boundary states
                //with respect to left boundary
                std::vector<bool> constraint_LeftBoundaryVariableAffectsMatchConstraint(this->ListOfVariablesAffectingLeftBoundary_probe_AEI_to_end.size(), false);
                for (size_t varIndex = 0; varIndex < this->ListOfVariablesAffectingLeftBoundary_probe_AEI_to_end.size(); ++varIndex)
                {
                    size_t Xindex = this->ListOfVariablesAffectingLeftBoundary_probe_AEI_to_end[varIndex];

                    for (size_t entryIndex = 0; entryIndex < this->DerivativesOfLeftBoundaryByVariable_probe_AEI_to_end[varIndex].size(); ++entryIndex)
                    {
                        size_t dIndex = this->DerivativesOfLeftBoundaryByVariable_probe_AEI_to_end[varIndex][entryIndex];
                        size_t stateIndex = std::get<1>(Derivatives_of_probe_StateAfterAEI[dIndex]);

                        //hooray, this variable affects the constraint. Make the Jacobian entry!
                        if (this->TruthTable_Forward_MatchConstraints_Derivative_wrt_EncodedStates[constraintIndex][stateIndex])
                        {
                            constraint_LeftBoundaryVariableAffectsMatchConstraint[varIndex] = true;
                            break;
                        }
                    }

                    if (constraint_LeftBoundaryVariableAffectsMatchConstraint[varIndex])
                    {
                        this->create_sparsity_entry(this->Findices_match_point_constraints_probe_AEI_to_end[constraintIndex],
                            Xindex,
                            this->G_indices_match_point_constraints_wrt_LeftBoundaryVariables_probe_AEI_to_end[constraintIndex]);
                    }
                    else
                        this->G_indices_match_point_constraints_wrt_LeftBoundaryVariables_probe_AEI_to_end[constraintIndex].push_back(0);
                }
                this->LeftBoundaryVariableAffectsMatchConstraint_probe_AEI_to_end.push_back(constraint_LeftBoundaryVariableAffectsMatchConstraint);

                //with respect to right boundary
                std::vector<bool> constraint_RightBoundaryVariableAffectsMatchConstraint(this->ListOfVariablesAffectingRightBoundary_probe_AEI_to_end.size(), false);
                for (size_t varIndex = 0; varIndex < this->ListOfVariablesAffectingRightBoundary_probe_AEI_to_end.size(); ++varIndex)
                {
                    size_t Xindex = this->ListOfVariablesAffectingRightBoundary_probe_AEI_to_end[varIndex];

                    for (size_t entryIndex = 0; entryIndex < this->DerivativesOfRightBoundaryByVariable_probe_AEI_to_end[varIndex].size(); ++entryIndex)
                    {
                        size_t dIndex = this->DerivativesOfRightBoundaryByVariable_probe_AEI_to_end[varIndex][entryIndex];
                        size_t stateIndex = std::get<1>(Derivatives_of_probe_StateAtEnd[dIndex]);

                        //hooray, this variable affects the constraint. Make the Jacobian entry!
                        if (this->TruthTable_Backward_MatchConstraints_Derivative_wrt_EncodedStates[constraintIndex][stateIndex])
                        {
                            constraint_RightBoundaryVariableAffectsMatchConstraint[varIndex] = true;
                            break;
                        }
                    }

                    if (constraint_RightBoundaryVariableAffectsMatchConstraint[varIndex])
                    {
                        this->create_sparsity_entry(this->Findices_match_point_constraints_probe_AEI_to_end[constraintIndex],
                            Xindex,
                            this->G_indices_match_point_constraints_wrt_RightBoundaryVariables_probe_AEI_to_end[constraintIndex]);
                    }
                    else
                        this->G_indices_match_point_constraints_wrt_RightBoundaryVariables_probe_AEI_to_end[constraintIndex].push_back(0);
                }
                this->RightBoundaryVariableAffectsMatchConstraint_probe_AEI_to_end.push_back(constraint_RightBoundaryVariableAffectsMatchConstraint);

                //derivatives with respect to time variables that affect boundary states
                //with respect to left boundary
                std::vector<bool> constraint_LeftBoundaryTimeVariableAffectsMatchConstraint(this->ListOfTimeVariablesAffectingLeftBoundary_probe_AEI_to_end.size(), false);
                for (size_t varIndex = 0; varIndex < this->ListOfTimeVariablesAffectingLeftBoundary_probe_AEI_to_end.size(); ++varIndex)
                {
                    size_t Xindex = this->ListOfTimeVariablesAffectingLeftBoundary_probe_AEI_to_end[varIndex];

                    for (size_t entryIndex = 0; entryIndex < this->DerivativesOfLeftBoundaryByTimeVariable_probe_AEI_to_end[varIndex].size(); ++entryIndex)
                    {
                        size_t dIndex = this->DerivativesOfLeftBoundaryByTimeVariable_probe_AEI_to_end[varIndex][entryIndex];
                        size_t stateIndex = std::get<1>(Derivatives_of_StateAfterAEI_wrt_Time[dIndex]);

                        //hooray, this variable affects the constraint. Make the Jacobian entry!
                        if (this->TruthTable_Forward_MatchConstraints_Derivative_wrt_EncodedStates[constraintIndex][stateIndex])
                        {
                            constraint_LeftBoundaryTimeVariableAffectsMatchConstraint[varIndex] = true;
                            break;
                        }
                    }

                    if (constraint_LeftBoundaryTimeVariableAffectsMatchConstraint[varIndex])
                    {
                        this->create_sparsity_entry(this->Findices_match_point_constraints_probe_AEI_to_end[constraintIndex],
                            Xindex,
                            this->G_indices_match_point_constraints_wrt_LeftBoundaryTimeVariables_probe_AEI_to_end[constraintIndex]);
                    }
                    else
                        this->G_indices_match_point_constraints_wrt_LeftBoundaryVariables_probe_AEI_to_end[constraintIndex].push_back(0);
                }
                this->LeftBoundaryTimeVariableAffectsMatchConstraint_probe_AEI_to_end.push_back(constraint_LeftBoundaryTimeVariableAffectsMatchConstraint);

                //with respect to right boundary
                std::vector<bool> constraint_RightBoundaryTimeVariableAffectsMatchConstraint(this->ListOfTimeVariablesAffectingRightBoundary_probe_AEI_to_end.size(), false);
                for (size_t varIndex = 0; varIndex < this->ListOfTimeVariablesAffectingRightBoundary_probe_AEI_to_end.size(); ++varIndex)
                {
                    size_t Xindex = this->ListOfTimeVariablesAffectingRightBoundary_probe_AEI_to_end[varIndex];

                    for (size_t entryIndex = 0; entryIndex < this->DerivativesOfRightBoundaryByTimeVariable_probe_AEI_to_end[varIndex].size(); ++entryIndex)
                    {
                        size_t dIndex = this->DerivativesOfRightBoundaryByTimeVariable_probe_AEI_to_end[varIndex][entryIndex];
                        size_t stateIndex = std::get<1>(Derivatives_of_probe_StateAtEnd_wrt_Time[dIndex]);

                        //hooray, this variable affects the constraint. Make the Jacobian entry!
                        if (this->TruthTable_Backward_MatchConstraints_Derivative_wrt_EncodedStates[constraintIndex][stateIndex])
                        {
                            constraint_RightBoundaryTimeVariableAffectsMatchConstraint[varIndex] = true;
                            break;
                        }
                    }

                    if (constraint_RightBoundaryTimeVariableAffectsMatchConstraint[varIndex])
                    {
                        this->create_sparsity_entry(this->Findices_match_point_constraints_probe_AEI_to_end[constraintIndex],
                            Xindex,
                            this->G_indices_match_point_constraints_wrt_RightBoundaryTimeVariables_probe_AEI_to_end[constraintIndex]);
                    }
                    else
                        this->G_indices_match_point_constraints_wrt_RightBoundaryVariables_probe_AEI_to_end[constraintIndex].push_back(0);
                }
                this->RightBoundaryTimeVariableAffectsMatchConstraint_probe_AEI_to_end.push_back(constraint_RightBoundaryTimeVariableAffectsMatchConstraint);
            }

            //Step 4: derivatives of the match point from separation to AEI due to probe separation impulse
            for (size_t constraintIndex = 0; constraintIndex < this->numMatchConstraints_probe; ++constraintIndex)
            {
                std::vector<size_t> stateGindex_probe_match_point_constraints_wrt_separation_impulse;

                for (size_t velocityIndex = 0; velocityIndex < 3; ++velocityIndex)
                {
                    size_t Xindex = this->Xindices_separation_direction_vector[velocityIndex];

                    this->create_sparsity_entry(this->Findices_match_point_constraints_probe_AEI_to_end[constraintIndex],
                        Xindex,
                        stateGindex_probe_match_point_constraints_wrt_separation_impulse);
                }
                this->Gindex_probe_match_point_constraints_wrt_separation_impulse.push_back(stateGindex_probe_match_point_constraints_wrt_separation_impulse);
            }
        }//end calcbounds_probe_match_point_constraints_from_AEI_to_end()

        void ProbeEntryPhase::calcbounds_communication_distance_constraint()
        {
            //Step 1: we will need the derivatives with respect to time and non-time variables of both the probe's entry state and the spacecraft's state at end of phase
            std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_spacecraft_StateAfterPhase = this->myArrivalEvent->get_Derivatives_of_StateAfterEvent();//Xindex, stateIndex, derivative value
            std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_probe_StateAtAEI = this->probeEntryArrivalEvent->get_Derivatives_of_StateBeforeEvent();//Xindex, stateIndex, derivative value
            std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_spacecraft_StateAfterPhase_wrt_Time = this->myArrivalEvent->get_Derivatives_of_StateAfterEvent_wrt_Time();//Xindex, stateIndex, derivative value
            std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_probe_StateAtAEI_wrt_Time = this->probeEntryArrivalEvent->get_Derivatives_of_StateBeforeEvent_wrt_Time();//Xindex, stateIndex, derivative value

            //Step 2: create the constraint
            this->Flowerbounds->push_back(this->myJourneyOptions->probe_communication_distance_bounds[0] / this->myJourneyOptions->probe_communication_distance_bounds[1] - 1.0);
            this->Fupperbounds->push_back(0.0);
            this->Fdescriptions->push_back(prefix + "probe communication distance constraint");

            //Step 2: loop over variables affecting position
            for (size_t stateIndex : {0, 1, 2})
            {
                //Step 2.1: derivatives due to the spacecraft
                {
                    std::vector<size_t> state_dIndex_with_respect_to_spacecraftStateAfterPhase;
                    std::vector<size_t> state_Gindex_constraint_wrt_spacecraftStateAfterPhase_variables;

                    //non-time variables
                    for (size_t dIndex = 0; dIndex < Derivatives_of_spacecraft_StateAfterPhase.size(); ++dIndex)
                    {
                        if (std::get<1>(Derivatives_of_spacecraft_StateAfterPhase[dIndex]) == stateIndex)
                        {
                            state_dIndex_with_respect_to_spacecraftStateAfterPhase.push_back(dIndex);

                            this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                                std::get<0>(Derivatives_of_spacecraft_StateAfterPhase[dIndex]),
                                state_Gindex_constraint_wrt_spacecraftStateAfterPhase_variables);
                        }
                    }
                    this->dIndex_probe_communication_distanceconstraint_with_respect_to_spacecraftStateAfterPhase.push_back(state_dIndex_with_respect_to_spacecraftStateAfterPhase);
                    this->Gindex_probe_communication_distanceconstraint_wrt_spacecraftStateAfterPhase_variables.push_back(state_Gindex_constraint_wrt_spacecraftStateAfterPhase_variables);

                    //time variables
                    std::vector<size_t> state_dIndex_with_respect_to_time_variables;
                    std::vector<size_t> state_Gindex_constraint_wrt_time_variables;
                    for (size_t dIndex = 0; dIndex < Derivatives_of_spacecraft_StateAfterPhase_wrt_Time.size(); ++dIndex)
                    {
                        if (std::get<1>(Derivatives_of_spacecraft_StateAfterPhase_wrt_Time[dIndex]) == stateIndex)
                        {
                            state_dIndex_with_respect_to_time_variables.push_back(dIndex);

                            this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                                std::get<0>(Derivatives_of_spacecraft_StateAfterPhase_wrt_Time[dIndex]),
                                state_Gindex_constraint_wrt_time_variables);
                        }
                    }
                    this->dIndex_probe_communication_distanceconstraint_with_respect_to_time_variables.push_back(state_dIndex_with_respect_to_time_variables);
                    this->Gindex_probe_communication_distanceconstraint_wrt_time_variables.push_back(state_Gindex_constraint_wrt_time_variables);
                }//end derivatives due to spacecraft

                //Step 2.2: derivatives due to the probe.
                //Note that we do not need new G entries due to time variables for the probe because it has the same time variables as the spacecraft.
                {
                    std::vector<size_t> state_dIndex_with_respect_to_probeStateAtAEI;
                    std::vector<size_t> state_Gindex_constraint_wrt_probeStateAtAEI_variables;

                    //non-time variables
                    for (size_t dIndex = 0; dIndex < Derivatives_of_probe_StateAtAEI.size(); ++dIndex)
                    {
                        if (std::get<1>(Derivatives_of_probe_StateAtAEI[dIndex]) == stateIndex)
                        {
                            state_dIndex_with_respect_to_probeStateAtAEI.push_back(dIndex);

                            this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                                std::get<0>(Derivatives_of_probe_StateAtAEI[dIndex]),
                                state_Gindex_constraint_wrt_probeStateAtAEI_variables);
                        }
                    }
                    this->dIndex_probe_communication_distanceconstraint_with_respect_to_probeStateAtEntry.push_back(state_dIndex_with_respect_to_probeStateAtAEI);
                    this->Gindex_probe_communication_distanceconstraint_wrt_probeStateAtAEntry_variables.push_back(state_Gindex_constraint_wrt_probeStateAtAEI_variables);
                }//end derivatives due to probe
            }//end loop over position derivatives
        }//end calcbounds_communication_distance_constraint
        
        //******************************************process methods
        void ProbeEntryPhase::process_phase(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //do we want to propagate STMs or not?
            this->total_number_of_states_to_integrate = needG
                ? this->numStatesToPropagate - 1 + 12 * 12
                : this->numStatesToPropagate - 1;

            //phase left boundary - on the combined stack
            this->process_phase_left_boundary(X, Xindex, F, Findex, G, needG);

            this->process_phase_flight_time(X, Xindex, F, Findex, G, needG);

            //spacecraft right boundary
            this->process_phase_right_boundary(X, Xindex, F, Findex, G, needG);

            //probe descent time
            if (this->myJourneyOptions->ModelProbeSecondPhase)
            {
                this->process_probe_descent_time(X, Xindex, F, Findex, G, needG);
            }

            //probe right boundary
            this->process_probe_boundary_events(X, Xindex, F, Findex, G, needG);

            //stuff that happens in the middle
            this->process_phase_main(X, Xindex, F, Findex, G, needG);
        }//end process_phase()

        void ProbeEntryPhase::process_probe_descent_time(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            this->probeDescentTime = X[Xindex++];
        }//end process_probe_descent_time()

        void ProbeEntryPhase::process_probe_boundary_events(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //Step 1: probe AEI arrival event
            //Step 1.1: call the event
            this->probeEntryArrivalEvent->reset_ETM();
            this->probeEntryArrivalEvent->process_event(X, Xindex, F, Findex, G, needG);
                       
            //Step 1.2: extract the states
            math::Matrix<doubleType>& AEI_arrival_state = this->probeEntryArrivalEvent->get_state_before_event();
            for (size_t stateIndex = 0; stateIndex < 8; ++stateIndex)
                this->probe_state_at_AEI_arrival(stateIndex) = AEI_arrival_state(stateIndex);

            //Step 2: probe AEI departure event
            //Step 2.1: call the event
            this->probeEntryDepartureEvent->reset_ETM();
            this->probeEntryDepartureEvent->process_event(X, Xindex, F, Findex, G, needG);

            //Step 2.2: extract the states
            math::Matrix<doubleType>& AEI_departure_state = this->probeEntryDepartureEvent->get_state_before_event();
            for (size_t stateIndex = 0; stateIndex < 8; ++stateIndex)
                this->probe_state_at_AEI_departure(stateIndex) = AEI_departure_state(stateIndex);

            //if modeling the second subphase
            if (this->myJourneyOptions->ModelProbeSecondPhase)
            {
                //Step 3: probe AEI departure event
                //Step 3.1: call the event
                this->probeEndEvent->reset_ETM();
                this->probeEndEvent->process_event(X, Xindex, F, Findex, G, needG);

                //Step 3.2: extract the states
                math::Matrix<doubleType>& probe_end_state = this->probeEndEvent->get_state_before_event();
                for (size_t stateIndex = 0; stateIndex < 8; ++stateIndex)
                    this->probe_state_at_end_of_phase(stateIndex) = probe_end_state(stateIndex);
            }
        }//end process_probe_entry()

        void ProbeEntryPhase::process_phase_main(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //Step 1: separate the spacecraft from the probe and apply the separation impulse
            this->process_separation(X, Xindex, F, Findex, G, needG);

            //Step 2: do spacecraft things
            this->MGAnDSMs_phase::process_phase_main(X, Xindex, F, Findex, G, needG);

            this->process_spacecraft_match_point_derivatives_due_to_probe_release(X, Xindex, F, Findex, G, needG);

            //Step 3: do probe things
            this->process_probe_separation_to_AEI_forward_half_phase(X, Xindex, F, Findex, G, needG);

            this->process_probe_separation_to_AEI_backward_half_phase(X, Xindex, F, Findex, G, needG);

            this->process_probe_match_point_constraints_from_separation_to_AEI(X, Xindex, F, Findex, G, needG);

            this->process_communication_distance_constraint(X, Xindex, F, Findex, G, needG);

            //if modeling the second subphase
            if (this->myJourneyOptions->ModelProbeSecondPhase)
            {
                this->process_probe_AEI_to_end_forward_half_phase(X, Xindex, F, Findex, G, needG);

                this->process_probe_AEI_to_end_backward_half_phase(X, Xindex, F, Findex, G, needG);

                this->process_probe_match_point_constraints_from_AEI_to_end(X, Xindex, F, Findex, G, needG);
            }
            
        }//end process_phase_main()

        void ProbeEntryPhase::process_separation(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //Step 1: initialize the state of the probe and spacecraft to be equal to the stack
            this->spacecraft_state_after_separation = this->state_at_beginning_of_phase;
            this->probe_state_after_separation = this->state_at_beginning_of_phase;

            //Step 2: figure out the mass of the spacecraft and probe
            this->probe_state_after_separation(6) = this->myJourneyOptions->probe_mass;
            this->spacecraft_state_after_separation(6) -= this->myJourneyOptions->probe_mass;

            //Step 3: compute and apply the separation impulse
            //Step 3.1: extract sparse representation of separation direction
            //and apply the velocity match constraints to force the separation impulse to match the probe entry velocity unit vector
            math::Matrix<doubleType> separation_unit_vector(3, 1, 0.0);
            math::Matrix<doubleType> probe_entry_velocity_unit_vector = this->probe_state_at_end_of_phase.getSubMatrix1D(3, 5).unitize();
            for (size_t velocityIndex : {0, 1, 2})
            {
                //Step 3.1.1: extract the separation impulse component
                separation_unit_vector(velocityIndex) = X[Xindex++];

                //Step 3.1.2: apply the velocity match constraint
                F[Findex++] = separation_unit_vector(velocityIndex) - probe_entry_velocity_unit_vector(velocityIndex);
            }//end loop over separation velocity unit vector components


            //Step 3.2: compute delta-v on spacecraft and probe
            this->spacecraft_separation_deltaV = -separation_unit_vector * this->myJourneyOptions->probe_separation_impulse / this->spacecraft_state_after_separation(6);
            this->probe_separation_deltaV = separation_unit_vector * this->myJourneyOptions->probe_separation_impulse / this->probe_state_after_separation(6);
            for (size_t vIndex : {0, 1, 2})
            {
                //spacecraft impulse is backwards
                this->spacecraft_state_after_separation(3 + vIndex) += this->spacecraft_separation_deltaV(vIndex);

                //probe impulse is forwards
                this->probe_state_after_separation(3 + vIndex) += this->probe_separation_deltaV(vIndex);
            }

            //Step 3.3: derivatives of separation impulse velocity match constraint
            if (needG)
            {
                //This constraint enforces that the separation impulse components match the *unit vector* of the probe's inertial entry velocity
                //because that's where the nose points, so the spacecraft slews to that inertial attitude before releasing the probe
                math::Matrix<doubleType> dProbeEntryVelocityUnit_dVelocityComponents = separation_unit_vector.unitDerivative(this->probe_state_at_end_of_phase.getSubMatrix1D(3, 5));
                std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_probe_StateAtEnd = this->probeEndEvent->get_Derivatives_of_StateBeforeEvent();//Xindex, stateIndex, derivative value
                std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_probe_StateAtEnd_wrt_Time = this->probeEndEvent->get_Derivatives_of_StateBeforeEvent_wrt_Time();//Xindex, stateIndex, derivative value
                
                for (size_t velocityIndex : {0, 1, 2})
                {
                    //Step 3.3.1: derivatives with respect to the separation vector itself
                    G[this->Gindices_separation_direction_vector_wrt_self[velocityIndex]] = 1.0;

                    //Step 3.3.2: derivatives with respect to variables affecting the probe entry velocity
                    //first clear everything
                    for (size_t Gindex : this->Gindex_separation_vector_velocity_match_constraint[velocityIndex])
                    {
                        G[Gindex] = 0.0;
                    }

                    for (size_t dIndexIndex = 0; dIndexIndex < this->dIndex_separation_vector_velocity_match_constraint[velocityIndex].size(); ++dIndexIndex)
                    {
                        size_t dIndex = this->dIndex_separation_vector_velocity_match_constraint[velocityIndex][dIndexIndex];
                        size_t Xindex = std::get<0>(Derivatives_of_probe_StateAtEnd[dIndex]);
                        size_t probeVelocityIndex = std::get<1>(Derivatives_of_probe_StateAtEnd[dIndex]) - 3;
                        size_t Gindex = this->Gindex_separation_vector_velocity_match_constraint[velocityIndex][dIndexIndex];

                        double TheDerivative = dProbeEntryVelocityUnit_dVelocityComponents(velocityIndex, probeVelocityIndex) _GETVALUE
                            * std::get<2>(Derivatives_of_probe_StateAtEnd[dIndex]);

                        G[Gindex] += -this->X_scale_factors->operator[](Xindex)
                            * TheDerivative;
                    }//end loop over derivative entries

                    //Step 3.3.3: derivatives with respect to time variables affecting the probe entry velocity
                    //first clear everything
                    for (size_t Gindex : this->Gindex_separation_vector_velocity_match_constraint_wrt_Time[velocityIndex])
                    {
                        G[Gindex] = 0.0;
                    }

                    for (size_t dIndexIndex = 0; dIndexIndex < this->dIndex_separation_vector_velocity_match_constraint_wrt_Time[velocityIndex].size(); ++dIndexIndex)
                    {
                        size_t dIndex = this->dIndex_separation_vector_velocity_match_constraint_wrt_Time[velocityIndex][dIndexIndex];
                        size_t Xindex = std::get<0>(Derivatives_of_probe_StateAtEnd_wrt_Time[dIndex]);
                        size_t probeVelocityIndex = std::get<1>(Derivatives_of_probe_StateAtEnd_wrt_Time[dIndex]) - 3;
                        size_t Gindex = this->Gindex_separation_vector_velocity_match_constraint_wrt_Time[velocityIndex][dIndexIndex];

                        double TheDerivative = dProbeEntryVelocityUnit_dVelocityComponents(velocityIndex, probeVelocityIndex) _GETVALUE
                            * std::get<2>(Derivatives_of_probe_StateAtEnd_wrt_Time[dIndex]);

                        G[Gindex] += -this->X_scale_factors->operator[](Xindex)
                            * TheDerivative;
                    }//end loop over derivative entries
                }//end loop over velocity components
            }//end derivatives of separation impulse

            //Step 4: copy the spacecraft state back into state_at_beginning_of_phase so that the parent MGAnDSMs_phase process() methods work properly
            this->state_at_beginning_of_phase = this->spacecraft_state_after_separation;
        }//end process_separation()

        void ProbeEntryPhase::process_spacecraft_match_point_derivatives_due_to_probe_release(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            if (needG)
            {
                math::Matrix<double> dStateNow_dDecisionVariable(this->numStatesToPropagate + 2, 1, 0.0);
                math::Matrix<double> dMatchPointState_dDecisionVariable(this->numStatesToPropagate + 2, 1, 0.0);

                //Step 1: derivatives due to impulse direction
                for (size_t Vindex = 0; Vindex < 3; ++Vindex)
                {
                    dStateNow_dDecisionVariable.assign_zeros();
                    dStateNow_dDecisionVariable(3 + Vindex) = (-this->myJourneyOptions->probe_separation_impulse / this->spacecraft_state_after_separation(6))_GETVALUE;

                    dMatchPointState_dDecisionVariable = this->ForwardHPTM * dStateNow_dDecisionVariable;

                    for (size_t constraintIndex = 0; constraintIndex < this->numMatchConstraints; ++constraintIndex)
                    {
                        size_t stateIndex = this->matchPointConstraintStateIndex[constraintIndex];
                        size_t Gindex = this->Gindex_spacecraft_match_point_constraints_wrt_separation_impulse[constraintIndex][Vindex];
                        size_t Xindex = this->jGvar->operator[](Gindex);

                        G[Gindex] = -this->X_scale_factors->operator[](Xindex)
                            * dMatchPointState_dDecisionVariable(stateIndex)
                            * this->continuity_constraint_scale_factors(constraintIndex);
                    }//end loop over states
                }//end loop over separation direction components

                //Step 2: derivatives due to spacecraft mass
                std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateBeforePhase = this->myDepartureEvent->get_Derivatives_of_StateAfterEvent();//Xindex, stateIndex, derivative value
                dStateNow_dDecisionVariable.assign_zeros();
                for (size_t vIndex : {0, 1, 2})
                {
                    dStateNow_dDecisionVariable(3 + vIndex) = -(X[this->Xindices_separation_direction_vector[vIndex]] * this->myJourneyOptions->probe_separation_impulse
                        / this->spacecraft_state_after_separation(6) / this->spacecraft_state_after_separation(6))_GETVALUE;
                }


                for (size_t dIndex : this->dIndex_spacecraft_left_mass)
                {
                    dMatchPointState_dDecisionVariable = this->ForwardHPTM * dStateNow_dDecisionVariable * std::get<2>(Derivatives_of_StateBeforePhase[dIndex]);

                    for (size_t constraintIndex = 0; constraintIndex < this->numMatchConstraints; ++constraintIndex)
                    {
                        size_t Gindex = this->Gindex_spacecraft_match_point_constraints_wrt_left_mass[constraintIndex];
                        size_t Xindex = this->jGvar->operator[](Gindex);
                        size_t stateIndex = this->matchPointConstraintStateIndex[constraintIndex];

                        G[Gindex] += this->X_scale_factors->operator[](Xindex)
                            * dMatchPointState_dDecisionVariable(stateIndex)
                            * this->continuity_constraint_scale_factors(constraintIndex);
                    }//end loop over states
                }

            }//end derivatives
        }//end process_spacecraft_match_point_derivatives_due_to_probe_release()

        void ProbeEntryPhase::process_probe_separation_to_AEI_forward_half_phase(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //Step 1: call the forward propagator
            this->probe_separation_to_AEI_Forward_dPropagatedStatedIndependentVariable.assign_zeros();
            this->probe_separation_to_AEI_ForwardHalfPhasePropagator->setCurrentEpoch(this->probe_state_after_separation(7));
            this->probe_separation_to_AEI_ForwardHalfPhasePropagator->setIndexOfEpochInStateVec(7);
            this->probe_separation_to_AEI_ForwardHalfPhasePropagator->setCurrentIndependentVariable(this->probe_state_after_separation(7));
            this->probe_separation_to_AEI_ForwardHalfPhasePropagator->propagate(this->PhaseFlightTime * this->probe_separation_to_AEI_MatchPointFraction, needG);

            if (this->isKeplerian) //otherwise this happens in the propagator
            {
                this->probe_first_subphase_match_point_state_minus(7) = this->probe_state_after_separation(7) + this->PhaseFlightTime * this->probe_separation_to_AEI_MatchPointFraction;
                this->probe_first_subphase_match_point_state_minus(6) = this->probe_state_after_separation(6);
                this->probe_first_subphase_match_point_state_minus(8) = this->probe_state_after_separation(8); //virtual fuel (this is zero anyway)
            }

            //Step 2: derivatives
            if (needG)
            {
                //Step 2.1: upper left 6x6 is the regular STM
                size_t stateMax = this->isKeplerian ? 6 : 10;
                for (size_t i = 0; i < stateMax; ++i)
                    for (size_t j = 0; j < stateMax; ++j)
                        this->probe_separation_to_AEI_ForwardSPTM(i, j) = this->probe_separation_to_AEI_ForwardSTM(i, j);

                //Step 2.2: turn the upper right 10 x 2 into the Phi_t terms developed by Lantoine
                for (size_t i = 0; i < stateMax; ++i)
                {
                    //this really ought to be cleaned up...
                    if (this->isKeplerian)
                    {
                        this->probe_separation_to_AEI_ForwardSPTM(i, this->stateIndex_phase_propagation_variable) = this->probe_separation_to_AEI_Forward_dPropagatedStatedIndependentVariable(i);
                    }
                    else //integrated propagator
                    {
                        this->probe_separation_to_AEI_ForwardSPTM(i, this->stateIndex_phase_propagation_variable) = this->probe_separation_to_AEI_ForwardSTM(i, 13);
                    }
                }

                //Step 2.3: form the HPTM
                this->probe_separation_to_AEI_ForwardHPTM = this->probe_separation_to_AEI_ForwardSPTM * this->TCMTM;
            }//end derivatives
        }//end process_probe_separation_to_AEI_forward_half_phase()

        void ProbeEntryPhase::process_probe_separation_to_AEI_backward_half_phase(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //Step 1: call the backward propagator
            this->probe_separation_to_AEI_Backward_dPropagatedStatedIndependentVariable.assign_zeros();
            this->probe_separation_to_AEI_BackwardHalfPhasePropagator->setCurrentEpoch(this->probe_state_at_AEI_arrival(7));
            this->probe_separation_to_AEI_BackwardHalfPhasePropagator->setIndexOfEpochInStateVec(7);
            this->probe_separation_to_AEI_BackwardHalfPhasePropagator->setCurrentIndependentVariable(this->probe_state_at_AEI_arrival(7));
            this->probe_separation_to_AEI_BackwardHalfPhasePropagator->propagate(-this->PhaseFlightTime * (1.0 - this->probe_separation_to_AEI_MatchPointFraction), needG);
            
            if (this->isKeplerian) //otherwise this happens in the propagator
            {
                this->probe_first_subphase_match_point_state_plus(6) = this->probe_state_at_AEI_arrival(6);
                this->probe_first_subphase_match_point_state_plus(7) = this->probe_state_at_AEI_arrival(7) - this->PhaseFlightTime * (1.0 - this->probe_separation_to_AEI_MatchPointFraction);
                this->probe_first_subphase_match_point_state_plus(8) = this->probe_state_at_AEI_arrival(8); //virtual fuel (this is zero anyway)
            }

            //Step 2: derivatives
            if (needG)
            {
                //Step 2.1: upper left 6x6 is the regular STM
                size_t stateMax = this->isKeplerian ? 6 : 10;
                for (size_t i = 0; i < stateMax; ++i)
                    for (size_t j = 0; j < stateMax; ++j)
                        this->probe_separation_to_AEI_BackwardSPTM(i, j) = this->probe_separation_to_AEI_BackwardSTM(i, j);
                //Step 2.2: turn the upper right 6 x 2 into the Phi_t terms developed by Lantoine
                for (size_t i = 0; i < 6; ++i)
                {
                    //this really ought to be cleaned up...
                    if (this->isKeplerian)
                    {
                        this->probe_separation_to_AEI_BackwardSPTM(i, this->stateIndex_phase_propagation_variable) = this->probe_separation_to_AEI_Backward_dPropagatedStatedIndependentVariable(i);
                    }
                    else //integrated propagator
                    {
                        this->probe_separation_to_AEI_BackwardSPTM(i, this->stateIndex_phase_propagation_variable) = this->probe_separation_to_AEI_BackwardSTM(i, 13);
                    }
                }

                //Step 2.3: form the HPTM
                this->probe_separation_to_AEI_BackwardHPTM = this->probe_separation_to_AEI_BackwardSPTM;
            }//end derivatives
        }//end process_probe_separation_to_AEI_backward_half_phase()

        void ProbeEntryPhase::process_probe_match_point_constraints_from_separation_to_AEI(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //Step 1: compute the probe match point constraints from separation to AEI
            for (size_t constraintIndex = 0; constraintIndex < this->numMatchConstraints_probe; ++constraintIndex)
            {
                size_t stateIndex = this->matchPointConstraintStateIndex[constraintIndex];
                F[Findex++] = (this->probe_first_subphase_match_point_state_plus(stateIndex) - this->probe_first_subphase_match_point_state_minus(stateIndex))
                    * this->continuity_constraint_scale_factors(constraintIndex);
            }


            //Step 2: derivatives
            if (needG)
            {
                std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateBeforePhase = this->myDepartureEvent->get_Derivatives_of_StateAfterEvent();//Xindex, stateIndex, derivative value
                std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_probe_StateAtAEI = this->probeEntryArrivalEvent->get_Derivatives_of_StateBeforeEvent();//Xindex, stateIndex, derivative value
                std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateBeforePhase_wrt_Time = this->myDepartureEvent->get_Derivatives_of_StateAfterEvent_wrt_Time();//Xindex, stateIndex, derivative value
                std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_probe_StateAtAEI_wrt_Time = this->probeEntryArrivalEvent->get_Derivatives_of_StateBeforeEvent_wrt_Time();//Xindex, stateIndex, derivative value

                math::Matrix<double> dBoundaryState_dDecisionVariable(this->probe_separation_to_AEI_ForwardHPTM.get_n(), 1, 0.0);
                math::Matrix<double> dMatchPointState_dDecisionVariable(this->probe_separation_to_AEI_ForwardHPTM.get_n(), 1, 0.0);
                
                for (size_t constraintIndex = 0; constraintIndex < this->numMatchConstraints_probe; ++constraintIndex)
                {
                    size_t stateIndex = this->matchPointConstraintStateIndex[constraintIndex];
                    //Step 2.1: left boundary
                    for (size_t varIndex = 0; varIndex < this->DerivativesOfLeftBoundaryByVariable_probe_separation_to_AEI.size(); ++varIndex)
                    {
                        if (this->LeftBoundaryVariableAffectsMatchConstraint_probe_separation_to_AEI[constraintIndex][varIndex])
                        {
                            dBoundaryState_dDecisionVariable.assign_zeros();

                            for (size_t entryIndex = 0; entryIndex < this->DerivativesOfLeftBoundaryByVariable_probe_separation_to_AEI[varIndex].size(); ++entryIndex)
                            {
                                size_t dIndex = this->DerivativesOfLeftBoundaryByVariable_probe_separation_to_AEI[varIndex][entryIndex];
                                size_t stateIndex = std::get<1>(Derivatives_of_StateBeforePhase[dIndex]);

                                dBoundaryState_dDecisionVariable(stateIndex) = std::get<2>(Derivatives_of_StateBeforePhase[dIndex]);
                            }

                            dMatchPointState_dDecisionVariable = this->probe_separation_to_AEI_ForwardHPTM * dBoundaryState_dDecisionVariable;


                            size_t Gindex = this->G_indices_match_point_constraints_wrt_LeftBoundaryVariables_probe_separation_to_AEI[constraintIndex][varIndex];
                            size_t Xindex = this->jGvar->operator[](Gindex);

                            G[Gindex] = this->X_scale_factors->operator[](Xindex)
                                * -dMatchPointState_dDecisionVariable(stateIndex)
                                * this->continuity_constraint_scale_factors(constraintIndex);
                        }
                    }

                    //Step 2.2: right boundary
                    for (size_t varIndex = 0; varIndex < this->DerivativesOfRightBoundaryByVariable_probe_separation_to_AEI.size(); ++varIndex)
                    {
                        if (this->RightBoundaryVariableAffectsMatchConstraint_probe_separation_to_AEI[constraintIndex][varIndex])
                        {
                            dBoundaryState_dDecisionVariable.assign_zeros();

                            for (size_t entryIndex = 0; entryIndex < this->DerivativesOfRightBoundaryByVariable_probe_separation_to_AEI[varIndex].size(); ++entryIndex)
                            {
                                size_t dIndex = this->DerivativesOfRightBoundaryByVariable_probe_separation_to_AEI[varIndex][entryIndex];
                                size_t stateIndex = std::get<1>(Derivatives_of_probe_StateAtAEI[dIndex]);

                                dBoundaryState_dDecisionVariable(stateIndex) = std::get<2>(Derivatives_of_probe_StateAtAEI[dIndex]);
                            }

                            dMatchPointState_dDecisionVariable = this->probe_separation_to_AEI_BackwardHPTM * dBoundaryState_dDecisionVariable;


                            size_t Gindex = this->G_indices_match_point_constraints_wrt_RightBoundaryVariables_probe_separation_to_AEI[constraintIndex][varIndex];
                            size_t Xindex = this->jGvar->operator[](Gindex);

                            G[Gindex] = this->X_scale_factors->operator[](Xindex)
                                * dMatchPointState_dDecisionVariable(stateIndex)
                                * this->continuity_constraint_scale_factors(constraintIndex);
                        }
                    }

                    //Step 2.3: left boundary epoch
                    for (size_t varIndex = 0; varIndex < this->DerivativesOfLeftBoundaryByTimeVariable_probe_separation_to_AEI.size(); ++varIndex)
                    {
                        if (this->LeftBoundaryTimeVariableAffectsMatchConstraint_probe_separation_to_AEI[constraintIndex][varIndex])
                        {
                            dBoundaryState_dDecisionVariable.assign_zeros();

                            for (size_t entryIndex = 0; entryIndex < this->DerivativesOfLeftBoundaryByTimeVariable_probe_separation_to_AEI[varIndex].size(); ++entryIndex)
                            {
                                size_t dIndex = this->DerivativesOfLeftBoundaryByTimeVariable_probe_separation_to_AEI[varIndex][entryIndex];
                                size_t stateIndex = std::get<1>(Derivatives_of_StateBeforePhase_wrt_Time[dIndex]);

                                dBoundaryState_dDecisionVariable(stateIndex) = std::get<2>(Derivatives_of_StateBeforePhase_wrt_Time[dIndex]);
                            }

                            dMatchPointState_dDecisionVariable = this->probe_separation_to_AEI_ForwardHPTM * dBoundaryState_dDecisionVariable;


                            size_t Gindex = this->G_indices_match_point_constraints_wrt_LeftBoundaryTimeVariables_probe_separation_to_AEI[constraintIndex][varIndex];
                            size_t Xindex = this->jGvar->operator[](Gindex);

                            G[Gindex] = this->X_scale_factors->operator[](Xindex)
                                * -dMatchPointState_dDecisionVariable(stateIndex)
                                * this->continuity_constraint_scale_factors(constraintIndex);
                        }
                    }

                    //Step 2.4: right boundary epoch
                    for (size_t varIndex = 0; varIndex < this->DerivativesOfRightBoundaryByTimeVariable_probe_separation_to_AEI.size() - 1; ++varIndex)
                    {
                        if (this->RightBoundaryTimeVariableAffectsMatchConstraint_probe_separation_to_AEI[constraintIndex][varIndex])
                        {
                            dBoundaryState_dDecisionVariable.assign_zeros();

                            for (size_t entryIndex = 0; entryIndex < this->DerivativesOfRightBoundaryByTimeVariable_probe_separation_to_AEI[varIndex].size(); ++entryIndex)
                            {
                                size_t dIndex = this->DerivativesOfRightBoundaryByTimeVariable_probe_separation_to_AEI[varIndex][entryIndex];
                                size_t stateIndex = std::get<1>(Derivatives_of_probe_StateAtAEI_wrt_Time[dIndex]);

                                dBoundaryState_dDecisionVariable(stateIndex) = std::get<2>(Derivatives_of_probe_StateAtAEI_wrt_Time[dIndex]);
                            }

                            dMatchPointState_dDecisionVariable = this->probe_separation_to_AEI_BackwardHPTM * dBoundaryState_dDecisionVariable;


                            size_t Gindex = this->G_indices_match_point_constraints_wrt_RightBoundaryTimeVariables_probe_separation_to_AEI[constraintIndex][varIndex];
                            size_t Xindex = this->jGvar->operator[](Gindex);

                            G[Gindex] += this->X_scale_factors->operator[](Xindex)
                                * dMatchPointState_dDecisionVariable(stateIndex)
                                * this->continuity_constraint_scale_factors(constraintIndex);
                        }
                    }

                    //Step 2.5: right boundary epoch due to current phase flight time
                    if (this->RightBoundaryTimeVariableAffectsMatchConstraint_probe_separation_to_AEI[constraintIndex].back())
                    {
                        dBoundaryState_dDecisionVariable.assign_zeros();

                        for (size_t entryIndex = 0; entryIndex < this->DerivativesOfRightBoundaryByTimeVariable_probe_separation_to_AEI.back().size(); ++entryIndex)
                        {
                            size_t dIndex = this->DerivativesOfRightBoundaryByTimeVariable_probe_separation_to_AEI.back()[entryIndex];
                            size_t stateIndex = std::get<1>(Derivatives_of_probe_StateAtAEI_wrt_Time[dIndex]);

                            dBoundaryState_dDecisionVariable(stateIndex) = std::get<2>(Derivatives_of_probe_StateAtAEI_wrt_Time[dIndex]);
                        }

                        dMatchPointState_dDecisionVariable = this->probe_separation_to_AEI_BackwardHPTM * dBoundaryState_dDecisionVariable;


                        size_t Gindex = this->G_indices_match_point_constraints_wrt_RightBoundaryTimeVariables_probe_separation_to_AEI[constraintIndex].back();
                        size_t Xindex = this->jGvar->operator[](Gindex);

                        G[Gindex] = this->X_scale_factors->operator[](Xindex)
                            * dMatchPointState_dDecisionVariable(stateIndex)
                            * this->continuity_constraint_scale_factors(constraintIndex);
                    }
                }//end loop over match point constraints

                //Step 2.6: derivatives wrt propagation time
                {
                    dBoundaryState_dDecisionVariable.assign_zeros();
                    dBoundaryState_dDecisionVariable(this->stateIndex_phase_propagation_variable) = 1.0;

                    //Forward
                    dMatchPointState_dDecisionVariable = this->probe_separation_to_AEI_ForwardHPTM * dBoundaryState_dDecisionVariable;

                    for (size_t constraintIndex = 0; constraintIndex < this->numMatchConstraints_probe; ++constraintIndex)
                    {
                        if (this->TruthTable_MatchConstraints_Derivative_wrt_PropagationTime[constraintIndex])
                        {
                            size_t stateIndex = this->matchPointConstraintStateIndex[constraintIndex];
                            size_t Gindex = this->G_indices_match_point_constraints_wrt_PhaseFlightTime_probe_separation_to_AEI[constraintIndex];
                            size_t Xindex = this->jGvar->operator[](Gindex);

                            G[Gindex] += -this->X_scale_factors->operator[](Xindex)
                                * dMatchPointState_dDecisionVariable(stateIndex)
                                * this->continuity_constraint_scale_factors(constraintIndex);
                        }
                    }

                    //Backward
                    dMatchPointState_dDecisionVariable = this->probe_separation_to_AEI_BackwardHPTM * dBoundaryState_dDecisionVariable;

                    for (size_t constraintIndex = 0; constraintIndex < this->numMatchConstraints_probe; ++constraintIndex)
                    {
                        if (this->TruthTable_MatchConstraints_Derivative_wrt_PropagationTime[constraintIndex])
                        {
                            size_t stateIndex = this->matchPointConstraintStateIndex[constraintIndex];
                            size_t Gindex = this->G_indices_match_point_constraints_wrt_PhaseFlightTime_probe_separation_to_AEI[constraintIndex];
                            size_t Xindex = this->jGvar->operator[](Gindex);

                            G[Gindex] += this->X_scale_factors->operator[](Xindex)
                                * dMatchPointState_dDecisionVariable(stateIndex)
                                * this->continuity_constraint_scale_factors(constraintIndex);
                        }
                    }
                }//end propagation time derivative

                //Step 2.7: derivatives due to separation impulse
                for (size_t Vindex = 0; Vindex < 3; ++Vindex)
                {
                    dBoundaryState_dDecisionVariable.assign_zeros();
                    dBoundaryState_dDecisionVariable(3 + Vindex) = (-this->myJourneyOptions->probe_separation_impulse / this->probe_state_after_separation(6))_GETVALUE;

                    dMatchPointState_dDecisionVariable = this->probe_separation_to_AEI_ForwardHPTM * dBoundaryState_dDecisionVariable;

                    for (size_t constraintIndex = 0; constraintIndex < this->numMatchConstraints_probe; ++constraintIndex)
                    {
                        size_t stateIndex = this->matchPointConstraintStateIndex[constraintIndex];
                        size_t Gindex = this->Gindex_probe_match_point_constraints_wrt_separation_impulse[constraintIndex][Vindex];
                        size_t Xindex = this->jGvar->operator[](Gindex);

                        G[Gindex] = this->X_scale_factors->operator[](Xindex)
                            * dMatchPointState_dDecisionVariable(stateIndex)
                            * this->continuity_constraint_scale_factors(constraintIndex);
                    }//end loop over states

                }//end loop over separation direction components
            }
        }//end process_probe_match_point_constraints_from_separation_to_AEI()

        void ProbeEntryPhase::process_probe_AEI_to_end_forward_half_phase(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //Step 1: call the forward propagator
            this->probe_AEI_to_end_Forward_dPropagatedStatedIndependentVariable.assign_zeros();
            this->probe_AEI_to_end_ForwardHalfPhasePropagator->setCurrentEpoch(this->probe_state_at_AEI_departure(7));
            this->probe_AEI_to_end_ForwardHalfPhasePropagator->setIndexOfEpochInStateVec(7);
            this->probe_AEI_to_end_ForwardHalfPhasePropagator->setCurrentIndependentVariable(this->probe_state_at_AEI_departure(7));
            this->probe_AEI_to_end_ForwardHalfPhasePropagator->propagate(this->probeDescentTime * this->probe_AEI_to_end_MatchPointFraction, needG);

            if (this->isKeplerian) //otherwise this happens in the propagator
            {
                this->probe_second_subphase_match_point_state_minus(7) = this->probe_state_at_AEI_departure(7) + this->probeDescentTime * this->probe_AEI_to_end_MatchPointFraction;
                this->probe_second_subphase_match_point_state_minus(6) = this->probe_state_at_AEI_departure(6);
                this->probe_second_subphase_match_point_state_minus(8) = this->probe_state_at_AEI_departure(8); //virtual fuel (this is zero anyway)
            }

            //Step 2: derivatives
            if (needG)
            {
                //Step 2.1: upper left 6x6 is the regular STM
                size_t stateMax = this->isKeplerian ? 6 : 10;
                for (size_t i = 0; i < stateMax; ++i)
                    for (size_t j = 0; j < stateMax; ++j)
                        this->probe_AEI_to_end_ForwardSPTM(i, j) = this->probe_AEI_to_end_ForwardSTM(i, j);

                //Step 2.2: turn the upper right 10 x 2 into the Phi_t terms developed by Lantoine
                for (size_t i = 0; i < stateMax; ++i)
                {
                    //this really ought to be cleaned up...
                    if (this->isKeplerian)
                    {
                        this->probe_AEI_to_end_ForwardSPTM(i, this->stateIndex_phase_propagation_variable) = this->probe_AEI_to_end_Forward_dPropagatedStatedIndependentVariable(i);
                    }
                    else //integrated propagator
                    {
                        this->probe_AEI_to_end_ForwardSPTM(i, this->stateIndex_phase_propagation_variable) = this->probe_AEI_to_end_ForwardSTM(i, 13);
                    }
                }

                //Step 2.3: form the HPTM
                this->probe_AEI_to_end_ForwardHPTM = this->probe_AEI_to_end_ForwardSPTM * this->TCMTM;
            }//end derivatives
        }//end process_probe_AEI_to_end_forward_half_phase()

        void ProbeEntryPhase::process_probe_AEI_to_end_backward_half_phase(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //Step 1: call the backward propagator
            this->probe_AEI_to_end_Backward_dPropagatedStatedIndependentVariable.assign_zeros();
            this->probe_AEI_to_end_BackwardHalfPhasePropagator->setCurrentEpoch(this->probe_state_at_end_of_phase(7));
            this->probe_AEI_to_end_BackwardHalfPhasePropagator->setIndexOfEpochInStateVec(7);
            this->probe_AEI_to_end_BackwardHalfPhasePropagator->setCurrentIndependentVariable(this->probe_state_at_end_of_phase(7));
            this->probe_AEI_to_end_BackwardHalfPhasePropagator->propagate(-this->probeDescentTime * (1.0 - this->probe_AEI_to_end_MatchPointFraction), needG);

            if (this->isKeplerian) //otherwise this happens in the propagator
            {
                this->probe_second_subphase_match_point_state_plus(6) = this->probe_state_at_end_of_phase(6);
                this->probe_second_subphase_match_point_state_plus(7) = this->probe_state_at_end_of_phase(7) - this->probeDescentTime * (1.0 - this->probe_AEI_to_end_MatchPointFraction);
                this->probe_second_subphase_match_point_state_plus(8) = this->probe_state_at_end_of_phase(8); //virtual fuel (this is zero anyway)
            }

            //Step 2: derivatives
            if (needG)
            {
                //Step 2.1: upper left 6x6 is the regular STM
                size_t stateMax = this->isKeplerian ? 6 : 10;
                for (size_t i = 0; i < stateMax; ++i)
                    for (size_t j = 0; j < stateMax; ++j)
                        this->probe_AEI_to_end_BackwardSPTM(i, j) = this->probe_AEI_to_end_BackwardSTM(i, j);
                //Step 2.2: turn the upper right 6 x 2 into the Phi_t terms developed by Lantoine
                for (size_t i = 0; i < 6; ++i)
                {
                    //this really ought to be cleaned up...
                    if (this->isKeplerian)
                    {
                        this->probe_AEI_to_end_BackwardSPTM(i, this->stateIndex_phase_propagation_variable) = this->probe_AEI_to_end_Backward_dPropagatedStatedIndependentVariable(i);
                    }
                    else //integrated propagator
                    {
                        this->probe_AEI_to_end_BackwardSPTM(i, this->stateIndex_phase_propagation_variable) = this->probe_AEI_to_end_BackwardSTM(i, 13);
                    }
                }

                //Step 2.3: form the HPTM
                this->probe_AEI_to_end_BackwardHPTM = this->probe_AEI_to_end_BackwardSPTM;
            }//end derivatives
        }//end process_probe_AEI_to_end_backward_half_phase()
        
        void ProbeEntryPhase::process_probe_match_point_constraints_from_AEI_to_end(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //Step 1: compute the probe match point constraints from AEI to end
            for (size_t constraintIndex = 0; constraintIndex < this->numMatchConstraints_probe; ++constraintIndex)
            {
                size_t stateIndex = this->matchPointConstraintStateIndex[constraintIndex];
                F[Findex++] = (this->probe_second_subphase_match_point_state_plus(stateIndex) - this->probe_second_subphase_match_point_state_minus(stateIndex))
                    * this->continuity_constraint_scale_factors(constraintIndex);
            }


            //Step 2: derivatives
            if (needG)
            {
                std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_probe_StateAfterAEI = this->probeEntryDepartureEvent->get_Derivatives_of_StateAfterEvent();//Xindex, stateIndex, derivative value
                std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_probe_StateAtEnd = this->probeEndEvent->get_Derivatives_of_StateBeforeEvent();//Xindex, stateIndex, derivative value
                std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateAfterAEI_wrt_Time = this->probeEntryDepartureEvent->get_Derivatives_of_StateAfterEvent_wrt_Time();//Xindex, stateIndex, derivative value
                std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_probe_StateAtEnd_wrt_Time = this->probeEndEvent->get_Derivatives_of_StateBeforeEvent_wrt_Time();//Xindex, stateIndex, derivative value

                math::Matrix<double> dBoundaryState_dDecisionVariable(this->probe_AEI_to_end_ForwardHPTM.get_n(), 1, 0.0);
                math::Matrix<double> dMatchPointState_dDecisionVariable(this->probe_AEI_to_end_ForwardHPTM.get_n(), 1, 0.0);

                for (size_t constraintIndex = 0; constraintIndex < this->numMatchConstraints_probe; ++constraintIndex)
                {
                    size_t stateIndex = this->matchPointConstraintStateIndex[constraintIndex];
                    //Step 2.1: left boundary
                    for (size_t varIndex = 0; varIndex < this->DerivativesOfLeftBoundaryByVariable_probe_AEI_to_end.size(); ++varIndex)
                    {
                        if (this->LeftBoundaryVariableAffectsMatchConstraint_probe_AEI_to_end[constraintIndex][varIndex])
                        {
                            dBoundaryState_dDecisionVariable.assign_zeros();

                            for (size_t entryIndex = 0; entryIndex < this->DerivativesOfLeftBoundaryByVariable_probe_AEI_to_end[varIndex].size(); ++entryIndex)
                            {
                                size_t dIndex = this->DerivativesOfLeftBoundaryByVariable_probe_AEI_to_end[varIndex][entryIndex];
                                size_t stateIndex = std::get<1>(Derivatives_of_probe_StateAfterAEI[dIndex]);

                                dBoundaryState_dDecisionVariable(stateIndex) = std::get<2>(Derivatives_of_probe_StateAfterAEI[dIndex]);
                            }

                            dMatchPointState_dDecisionVariable = this->probe_AEI_to_end_ForwardHPTM * dBoundaryState_dDecisionVariable;


                            size_t Gindex = this->G_indices_match_point_constraints_wrt_LeftBoundaryVariables_probe_AEI_to_end[constraintIndex][varIndex];
                            size_t Xindex = this->jGvar->operator[](Gindex);

                            G[Gindex] = this->X_scale_factors->operator[](Xindex)
                                * -dMatchPointState_dDecisionVariable(stateIndex)
                                * this->continuity_constraint_scale_factors(constraintIndex);
                        }
                    }

                    //Step 2.2: right boundary
                    for (size_t varIndex = 0; varIndex < this->DerivativesOfRightBoundaryByVariable_probe_AEI_to_end.size(); ++varIndex)
                    {
                        if (this->RightBoundaryVariableAffectsMatchConstraint_probe_AEI_to_end[constraintIndex][varIndex])
                        {
                            dBoundaryState_dDecisionVariable.assign_zeros();

                            for (size_t entryIndex = 0; entryIndex < this->DerivativesOfRightBoundaryByVariable_probe_AEI_to_end[varIndex].size(); ++entryIndex)
                            {
                                size_t dIndex = this->DerivativesOfRightBoundaryByVariable_probe_AEI_to_end[varIndex][entryIndex];
                                size_t stateIndex = std::get<1>(Derivatives_of_probe_StateAtEnd[dIndex]);

                                dBoundaryState_dDecisionVariable(stateIndex) = std::get<2>(Derivatives_of_probe_StateAtEnd[dIndex]);
                            }

                            dMatchPointState_dDecisionVariable = this->probe_AEI_to_end_BackwardHPTM * dBoundaryState_dDecisionVariable;


                            size_t Gindex = this->G_indices_match_point_constraints_wrt_RightBoundaryVariables_probe_AEI_to_end[constraintIndex][varIndex];
                            size_t Xindex = this->jGvar->operator[](Gindex);

                            G[Gindex] = this->X_scale_factors->operator[](Xindex)
                                * dMatchPointState_dDecisionVariable(stateIndex)
                                * this->continuity_constraint_scale_factors(constraintIndex);
                        }
                    }

                    //Step 2.3: left boundary epoch
                    for (size_t varIndex = 0; varIndex < this->DerivativesOfLeftBoundaryByTimeVariable_probe_AEI_to_end.size(); ++varIndex)
                    {
                        if (this->LeftBoundaryTimeVariableAffectsMatchConstraint_probe_AEI_to_end[constraintIndex][varIndex])
                        {
                            dBoundaryState_dDecisionVariable.assign_zeros();

                            for (size_t entryIndex = 0; entryIndex < this->DerivativesOfLeftBoundaryByTimeVariable_probe_AEI_to_end[varIndex].size(); ++entryIndex)
                            {
                                size_t dIndex = this->DerivativesOfLeftBoundaryByTimeVariable_probe_AEI_to_end[varIndex][entryIndex];
                                size_t stateIndex = std::get<1>(Derivatives_of_StateAfterAEI_wrt_Time[dIndex]);

                                dBoundaryState_dDecisionVariable(stateIndex) = std::get<2>(Derivatives_of_StateAfterAEI_wrt_Time[dIndex]);
                            }

                            dMatchPointState_dDecisionVariable = this->probe_AEI_to_end_ForwardHPTM * dBoundaryState_dDecisionVariable;


                            size_t Gindex = this->G_indices_match_point_constraints_wrt_LeftBoundaryTimeVariables_probe_AEI_to_end[constraintIndex][varIndex];
                            size_t Xindex = this->jGvar->operator[](Gindex);

                            G[Gindex] = this->X_scale_factors->operator[](Xindex)
                                * -dMatchPointState_dDecisionVariable(stateIndex)
                                * this->continuity_constraint_scale_factors(constraintIndex);
                        }
                    }

                    //Step 2.4: right boundary epoch
                    for (size_t varIndex = 0; varIndex < this->DerivativesOfRightBoundaryByTimeVariable_probe_AEI_to_end.size() - 1; ++varIndex)
                    {
                        if (this->RightBoundaryTimeVariableAffectsMatchConstraint_probe_AEI_to_end[constraintIndex][varIndex])
                        {
                            dBoundaryState_dDecisionVariable.assign_zeros();

                            for (size_t entryIndex = 0; entryIndex < this->DerivativesOfRightBoundaryByTimeVariable_probe_AEI_to_end[varIndex].size(); ++entryIndex)
                            {
                                size_t dIndex = this->DerivativesOfRightBoundaryByTimeVariable_probe_AEI_to_end[varIndex][entryIndex];
                                size_t stateIndex = std::get<1>(Derivatives_of_probe_StateAtEnd_wrt_Time[dIndex]);

                                dBoundaryState_dDecisionVariable(stateIndex) = std::get<2>(Derivatives_of_probe_StateAtEnd_wrt_Time[dIndex]);
                            }

                            dMatchPointState_dDecisionVariable = this->probe_AEI_to_end_BackwardHPTM * dBoundaryState_dDecisionVariable;


                            size_t Gindex = this->G_indices_match_point_constraints_wrt_RightBoundaryTimeVariables_probe_AEI_to_end[constraintIndex][varIndex];
                            size_t Xindex = this->jGvar->operator[](Gindex);

                            G[Gindex] += this->X_scale_factors->operator[](Xindex)
                                * dMatchPointState_dDecisionVariable(stateIndex)
                                * this->continuity_constraint_scale_factors(constraintIndex);
                        }
                    }

                    //Step 2.5: right boundary epoch due to probe descent time
                    if (this->RightBoundaryTimeVariableAffectsMatchConstraint_probe_AEI_to_end[constraintIndex].back())
                    {
                        dBoundaryState_dDecisionVariable.assign_zeros();

                        for (size_t entryIndex = 0; entryIndex < this->DerivativesOfRightBoundaryByTimeVariable_probe_AEI_to_end.back().size(); ++entryIndex)
                        {
                            size_t dIndex = this->DerivativesOfRightBoundaryByTimeVariable_probe_AEI_to_end.back()[entryIndex];
                            size_t stateIndex = std::get<1>(Derivatives_of_probe_StateAtEnd_wrt_Time[dIndex]);

                            dBoundaryState_dDecisionVariable(stateIndex) = std::get<2>(Derivatives_of_probe_StateAtEnd_wrt_Time[dIndex]);
                        }

                        dMatchPointState_dDecisionVariable = this->probe_AEI_to_end_BackwardHPTM * dBoundaryState_dDecisionVariable;


                        size_t Gindex = this->G_indices_match_point_constraints_wrt_RightBoundaryTimeVariables_probe_AEI_to_end[constraintIndex].back();
                        size_t Xindex = this->jGvar->operator[](Gindex);

                        G[Gindex] = this->X_scale_factors->operator[](Xindex)
                            * dMatchPointState_dDecisionVariable(stateIndex)
                            * this->continuity_constraint_scale_factors(constraintIndex);
                    }
                }//end loop over match point constraints

                //Step 2.6: derivatives wrt propagation time
                {
                    dBoundaryState_dDecisionVariable.assign_zeros();
                    dBoundaryState_dDecisionVariable(this->stateIndex_phase_propagation_variable) = 1.0;

                    //Forward
                    dMatchPointState_dDecisionVariable = this->probe_AEI_to_end_ForwardHPTM * dBoundaryState_dDecisionVariable;

                    for (size_t constraintIndex = 0; constraintIndex < this->numMatchConstraints_probe; ++constraintIndex)
                    {
                        if (this->TruthTable_MatchConstraints_Derivative_wrt_PropagationTime[constraintIndex])
                        {
                            size_t stateIndex = this->matchPointConstraintStateIndex[constraintIndex];
                            size_t Gindex = this->G_indices_match_point_constraints_wrt_PhaseFlightTime_probe_AEI_to_end[constraintIndex];
                            size_t Xindex = this->jGvar->operator[](Gindex);

                            G[Gindex] += -this->X_scale_factors->operator[](Xindex)
                                * dMatchPointState_dDecisionVariable(stateIndex)
                                * this->continuity_constraint_scale_factors(constraintIndex);
                        }
                    }

                    //Backward
                    dMatchPointState_dDecisionVariable = this->probe_AEI_to_end_BackwardHPTM * dBoundaryState_dDecisionVariable;

                    for (size_t constraintIndex = 0; constraintIndex < this->numMatchConstraints_probe; ++constraintIndex)
                    {
                        if (this->TruthTable_MatchConstraints_Derivative_wrt_PropagationTime[constraintIndex])
                        {
                            size_t stateIndex = this->matchPointConstraintStateIndex[constraintIndex];
                            size_t Gindex = this->G_indices_match_point_constraints_wrt_PhaseFlightTime_probe_AEI_to_end[constraintIndex];
                            size_t Xindex = this->jGvar->operator[](Gindex);

                            G[Gindex] += this->X_scale_factors->operator[](Xindex)
                                * dMatchPointState_dDecisionVariable(stateIndex)
                                * this->continuity_constraint_scale_factors(constraintIndex);
                        }
                    }
                }//end propagation time derivative
            }
        }//end process_probe_AEI_to_end_match_point_constraints()
        
        void ProbeEntryPhase::process_communication_distance_constraint(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //Step 1: compute the distance between the probe and spacecraft
            math::Matrix<doubleType> R_from_probe_to_spacecraft = this->state_at_end_of_phase.getSubMatrix1D(0, 2) - this->probe_state_at_AEI_arrival.getSubMatrix1D(0, 2);
            doubleType r_from_probe_to_spacecraft = R_from_probe_to_spacecraft.norm();

            //Step 2: apply the constraint
            F[Findex++] = r_from_probe_to_spacecraft / this->myJourneyOptions->probe_communication_distance_bounds[1] - 1.0;

            //Step 3: derivatives
            if (needG)
            {
                //Step 3.1: derivative components
                double r = r_from_probe_to_spacecraft _GETVALUE;
                double xsc = this->state_at_end_of_phase(0)_GETVALUE;
                double ysc = this->state_at_end_of_phase(1)_GETVALUE;
                double zsc = this->state_at_end_of_phase(2)_GETVALUE;
                double xp = this->probe_state_at_end_of_phase(0)_GETVALUE;
                double yp = this->probe_state_at_end_of_phase(1)_GETVALUE;
                double zp = this->probe_state_at_end_of_phase(2)_GETVALUE;
                double dr_dxsc = (-xp + xsc) / r;
                double dr_dysc = (-yp + ysc) / r;
                double dr_dzsc = (-zp + zsc) / r;
                double dr_dxp = (xp - xsc) / r;
                double dr_dyp = (yp - ysc) / r;
                double dr_dzp = (zp - zsc) / r;
                math::Matrix<double> dr_dRsc(3, 1, { dr_dxsc, dr_dysc, dr_dzsc });
                math::Matrix<double> dr_dRp(3, 1, { dr_dxp, dr_dyp, dr_dzp });

                std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_spacecraft_StateAfterPhase = this->myArrivalEvent->get_Derivatives_of_StateAfterEvent();//Xindex, stateIndex, derivative value
                std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_probe_StateAtAEI = this->probeEntryArrivalEvent->get_Derivatives_of_StateBeforeEvent();//Xindex, stateIndex, derivative value
                std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_spacecraft_StateAfterPhase_wrt_Time = this->myArrivalEvent->get_Derivatives_of_StateAfterEvent_wrt_Time();//Xindex, stateIndex, derivative value
                std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_probe_StateAtAEI_wrt_Time = this->probeEntryArrivalEvent->get_Derivatives_of_StateBeforeEvent_wrt_Time();//Xindex, stateIndex, derivative value

                //Step 3.2: zero all derivatives
                for (size_t stateIndex : {0, 1, 2})
                {
                    for (size_t Gindex : this->Gindex_probe_communication_distanceconstraint_wrt_spacecraftStateAfterPhase_variables[stateIndex])
                    {
                        G[Gindex] = 0.0;
                    }

                    for (size_t Gindex : this->Gindex_probe_communication_distanceconstraint_wrt_probeStateAtAEntry_variables[stateIndex])
                    {
                        G[Gindex] = 0.0;
                    }

                    for (size_t Gindex : this->Gindex_probe_communication_distanceconstraint_wrt_time_variables[stateIndex])
                    {
                        G[Gindex] = 0.0;
                    }
                }//end zeroing

                for (size_t stateIndex : {0, 1, 2})
                {
                    //Step 3.3: derivatives with respect to non-time variables that affect the spacecraft
                    for (size_t entryIndex = 0; entryIndex < this->dIndex_probe_communication_distanceconstraint_with_respect_to_spacecraftStateAfterPhase[stateIndex].size(); ++entryIndex)
                    {
                        size_t dIndex = this->dIndex_probe_communication_distanceconstraint_with_respect_to_spacecraftStateAfterPhase[stateIndex][entryIndex];
                        size_t Gindex = this->Gindex_probe_communication_distanceconstraint_wrt_spacecraftStateAfterPhase_variables[stateIndex][entryIndex];
                        size_t Xindex = std::get<0>(Derivatives_of_spacecraft_StateAfterPhase[dIndex]);

                        double TheDerivative = dr_dRsc(stateIndex) * std::get<2>(Derivatives_of_spacecraft_StateAfterPhase[dIndex]);

                        G[Gindex] += this->X_scale_factors->operator[](Xindex)
                            * TheDerivative
                            / this->myJourneyOptions->probe_communication_distance_bounds[1];
                    }//end loop over derivative entries

                    //Step 3.4: derivatives with respect to non-time varibales that affect the probe
                    for (size_t entryIndex = 0; entryIndex < this->dIndex_probe_communication_distanceconstraint_with_respect_to_probeStateAtEntry[stateIndex].size(); ++entryIndex)
                    {
                        size_t dIndex = this->dIndex_probe_communication_distanceconstraint_with_respect_to_probeStateAtEntry[stateIndex][entryIndex];
                        size_t Gindex = this->Gindex_probe_communication_distanceconstraint_wrt_probeStateAtAEntry_variables[stateIndex][entryIndex];
                        size_t Xindex = std::get<0>(Derivatives_of_probe_StateAtAEI[dIndex]);

                        double TheDerivative = dr_dRp(stateIndex) * std::get<2>(Derivatives_of_probe_StateAtAEI[dIndex]);

                        G[Gindex] += this->X_scale_factors->operator[](Xindex)
                            * TheDerivative
                            / this->myJourneyOptions->probe_communication_distance_bounds[1];
                    }//end loop over derivative entries

                    //Step 3.5: derivatives with respect to time variables that affect both the spacecraft and the probe
                    for (size_t entryIndex = 0; entryIndex < this->dIndex_probe_communication_distanceconstraint_with_respect_to_time_variables[stateIndex].size(); ++entryIndex)
                    {
                        size_t dIndex = this->dIndex_probe_communication_distanceconstraint_with_respect_to_time_variables[stateIndex][entryIndex];
                        size_t Gindex = this->Gindex_probe_communication_distanceconstraint_wrt_time_variables[stateIndex][entryIndex];
                        size_t Xindex = std::get<0>(Derivatives_of_probe_StateAtAEI_wrt_Time[dIndex]);

                        double TheDerivative = dr_dRsc(stateIndex) * std::get<2>(Derivatives_of_spacecraft_StateAfterPhase_wrt_Time[dIndex])
                                             - dr_dRp(stateIndex) * std::get<2>(Derivatives_of_probe_StateAtAEI_wrt_Time[dIndex]);

                        G[Gindex] += this->X_scale_factors->operator[](Xindex)
                            * TheDerivative
                            / this->myJourneyOptions->probe_communication_distance_bounds[1];
                    }//end loop over derivative entries
                }//end loop over states
            }//end derivatives
        }//end process_communication_distance_constraint()

        //*****************************************output methods
        void ProbeEntryPhase::output(std::ofstream& outputfile,
            size_t& eventcount)
        {
            //this is identical to MGAnDSMs except that (1) it includes the separation delta-v and (2) it calls output_probe() to write the probe output file

            //Step 1: output the departure event
            this->myDepartureEvent->output(outputfile, this->LaunchDate _GETVALUE, eventcount);
            this->mySpacecraft->setActiveStage(this->stageIndex);

            //Step 2: output the initial TCM if applicable
            if (this->hasInitialTCM)
            {
                this->output_initial_TCM(outputfile, eventcount);
            }

            //Step 3: output the separation impulse
            //we are relying on output_power and output_active_power already being set by output_initial_TCM()
            {
                phase::write_output_line(outputfile,//outputfile
                    eventcount,//eventcount
                    "separation",//event_type
                    "deep-space",//event_location
                    0.0,// timestep_size,
                    -1,//flyby_altitude,
                    0,//BdotR
                    0,//BdotT
                    0,//angle1
                    0,//angle2
                    0,//C3
                    this->spacecraft_state_after_separation,//state
                    this->spacecraft_separation_deltaV,//dV
                    math::Matrix<doubleType>(3, 1, 0.0),//ThrustVector
                    this->spacecraft_separation_deltaV.norm(),//dVmag
                    0.0,//Thrust
                    -1,//Isp
                    this->output_power,//AvailPower
                    0.0,//mdot
                    0,//number_of_active_engines
                    this->output_active_power,
                    "none");//active_power)
            }
            

            //Step 4: output the main part of the phase
            this->output_phase_main(outputfile, eventcount);

            //Step 5: print the arrival event
            this->myArrivalEvent->output(outputfile, this->LaunchDate _GETVALUE, eventcount);
            
            //Step 6: the probe gets *its own output file*
            this->output_probe();
        }//end output()

        void ProbeEntryPhase::output_probe()
        {
            //The probe gets its own output file
            size_t eventcount = 0;
            doubleType ForwardOutputTimestep_separation_to_AEI = this->PhaseFlightTime * this->probe_separation_to_AEI_MatchPointFraction / (this->myOptions->num_timesteps / 2);
            doubleType BackwardOutputTimestep_separation_to_AEI = this->PhaseFlightTime * (1.0 - this->probe_separation_to_AEI_MatchPointFraction) / (this->myOptions->num_timesteps / 2);
            doubleType ForwardOutputTimestep_AEI_to_end = this->probeDescentTime * this->probe_AEI_to_end_MatchPointFraction / (this->myOptions->num_timesteps / 2);
            doubleType BackwardOutputTimestep_AEI_to_end = this->probeDescentTime * (1.0 - this->probe_AEI_to_end_MatchPointFraction) / (this->myOptions->num_timesteps / 2);
            math::Matrix<doubleType> empty3(3, 1, 0.0);

            //Step 1: create the file
            std::ofstream outputfile(this->myOptions->working_directory + "/" + this->myOptions->mission_name + "_probe.emtg");

            //Step 2: header
            //copied from journey.cpp. I tried to do a cute pointer-to-parent thing to avoid the copy but it blew up
            {
                outputfile.precision(20);

                outputfile << std::endl;
                outputfile << "Journey: " << 0 << std::endl;
                outputfile << "Journey name: " << this->myJourneyOptions->journey_name << "_probe_trajectory" << std::endl;
                outputfile << "Central Body: " << this->myUniverse->central_body_name << std::endl;
                outputfile << "Radius (km): " << this->myUniverse->central_body_radius << std::endl;
                outputfile << "mu (km^3/s^2): " << this->myUniverse->mu << std::endl;
                outputfile << "Characteristic length unit (km): " << this->myUniverse->LU << std::endl;
                std::vector<std::string> FrameDefinitions({ "ICRF","J2000_BCI","J2000_BCF","TrueOfDate_BCI","TrueOfDate_BCF" });
                outputfile << "Frame: " << FrameDefinitions[this->myOptions->output_file_frame] << std::endl;

                if (this->myOptions->output_file_frame > ReferenceFrame::ICRF)
                {
                    outputfile << "alpha0: " << this->myUniverse->LocalFrame.get_alpha0() << std::endl;

                    if (this->myOptions->output_file_frame == ReferenceFrame::TrueOfDate_BCI || this->myOptions->output_file_frame == ReferenceFrame::TrueOfDate_BCF)
                        outputfile << "alphadot: " << this->myUniverse->LocalFrame.get_alphadot() << std::endl;

                    outputfile << "delta0: " << this->myUniverse->LocalFrame.get_delta0() << std::endl;

                    if (this->myOptions->output_file_frame == ReferenceFrame::TrueOfDate_BCI || this->myOptions->output_file_frame == ReferenceFrame::TrueOfDate_BCF)
                        outputfile << "deltadot: " << this->myUniverse->LocalFrame.get_deltadot() << std::endl;

                    if (this->myOptions->output_file_frame == ReferenceFrame::J2000_BCF || this->myOptions->output_file_frame == ReferenceFrame::TrueOfDate_BCF)
                    {
                        outputfile << "W0: " << this->myUniverse->LocalFrame.getW0() << std::endl;

                        if (this->myOptions->output_file_frame == ReferenceFrame::TrueOfDate_BCF)
                            outputfile << "Wdot: " << this->myUniverse->LocalFrame.getWdot() << std::endl;
                    }
                }
                outputfile << std::endl;

                //next, column headers

                //column headers line 1
                outputfile.width(5); outputfile << "#";
                outputfile.width(3); outputfile << " | ";
                outputfile.width(16); outputfile << "JulianDate";
                outputfile.width(3); outputfile << " | ";
                outputfile.width(11); outputfile << "MM/DD/YYYY";
                outputfile.width(3); outputfile << " | ";
                outputfile.width(12); outputfile << "event type";
                outputfile.width(3); outputfile << " | ";
                outputfile.width(25); outputfile << "location";
                outputfile.width(3); outputfile << " | ";
                outputfile.width(15); outputfile << "step size";
                outputfile.width(3); outputfile << " | ";
                outputfile.width(19); outputfile << "altitude";
                outputfile.width(3); outputfile << " | ";
                outputfile.width(19); outputfile << "BdotR";
                outputfile.width(3); outputfile << " | ";
                outputfile.width(19); outputfile << "BdotT";
                outputfile.width(3); outputfile << " | ";
                outputfile.width(8); outputfile << "RA";
                outputfile.width(3); outputfile << " | ";
                outputfile.width(8); outputfile << "DEC";
                outputfile.width(3); outputfile << " | ";
                outputfile.width(14); outputfile << "C3";
                outputfile.width(3); outputfile << " | ";
                outputfile.width(19); outputfile << " x";
                outputfile.width(3); outputfile << " | ";
                outputfile.width(19); outputfile << " y";
                outputfile.width(3); outputfile << " | ";
                outputfile.width(19); outputfile << " z";
                outputfile.width(3); outputfile << " | ";
                outputfile.width(19); outputfile << " xdot";
                outputfile.width(3); outputfile << " | ";
                outputfile.width(19); outputfile << " ydot";
                outputfile.width(3); outputfile << " | ";
                outputfile.width(19); outputfile << " zdot";
                outputfile.width(3); outputfile << " | ";
                outputfile.width(19); outputfile << " dV_x";
                outputfile.width(3); outputfile << " | ";
                outputfile.width(19); outputfile << " dV_y";
                outputfile.width(3); outputfile << " | ";
                outputfile.width(19); outputfile << " dV_z";
                outputfile.width(3); outputfile << " | ";
                outputfile.width(19); outputfile << " T_x";
                outputfile.width(3); outputfile << " | ";
                outputfile.width(19); outputfile << " T_y";
                outputfile.width(3); outputfile << " | ";
                outputfile.width(19); outputfile << " T_z";
                outputfile.width(3); outputfile << " | ";
                outputfile.width(17); outputfile << "|dV| (km/s)";
                outputfile.width(3); outputfile << " | ";
                outputfile.width(14); outputfile << "Avail. Thrust";
                outputfile.width(3); outputfile << " | ";
                outputfile.width(14); outputfile << "Isp";
                outputfile.width(3); outputfile << " | ";
                outputfile.width(14); outputfile << "Avail. Power";
                outputfile.width(3); outputfile << " | ";
                outputfile.width(19); outputfile << "Mass Flow";
                outputfile.width(3); outputfile << " | ";
                outputfile.width(14); outputfile << "mass";
                outputfile.width(3); outputfile << " | ";
                outputfile.width(14); outputfile << "number of";
                outputfile.width(3); outputfile << " | ";
                outputfile.width(14); outputfile << "active power";
                outputfile.width(3); outputfile << " | ";
                outputfile.width(14); outputfile << "throttle level";
                outputfile.width(3); outputfile << " | ";
                outputfile << std::endl;

                //column headers line 2
                outputfile.width(5); outputfile << "";
                outputfile.width(3); outputfile << " | ";
                outputfile.width(16); outputfile << " (ET/TDB)";
                outputfile.width(3); outputfile << " | ";
                outputfile.width(11); outputfile << "";
                outputfile.width(3); outputfile << " | ";
                outputfile.width(12); outputfile << "";
                outputfile.width(3); outputfile << " | ";
                outputfile.width(25); outputfile << "";
                outputfile.width(3); outputfile << " | ";
                outputfile.width(15); outputfile << "(days)";
                outputfile.width(3); outputfile << " | ";
                outputfile.width(19); outputfile << "(km)";
                outputfile.width(3); outputfile << " | ";
                outputfile.width(19); outputfile << "(km)";
                outputfile.width(3); outputfile << " | ";
                outputfile.width(19); outputfile << "(km)";
                outputfile.width(3); outputfile << " | ";
                outputfile.width(8); outputfile << "degrees";
                outputfile.width(3); outputfile << " | ";
                outputfile.width(8); outputfile << "degrees";
                outputfile.width(3); outputfile << " | ";
                outputfile.width(14); outputfile << "(km^2/s^2)";
                outputfile.width(3); outputfile << " | ";
                outputfile.width(19); outputfile << "(km)";
                outputfile.width(3); outputfile << " | ";
                outputfile.width(19); outputfile << "(km)";
                outputfile.width(3); outputfile << " | ";
                outputfile.width(19); outputfile << "(km)";
                outputfile.width(3); outputfile << " | ";
                outputfile.width(19); outputfile << "(km/s)";
                outputfile.width(3); outputfile << " | ";
                outputfile.width(19); outputfile << "(km/s)";
                outputfile.width(3); outputfile << " | ";
                outputfile.width(19); outputfile << "(km/s)";
                outputfile.width(3); outputfile << " | ";
                outputfile.width(19); outputfile << "(km/s)";
                outputfile.width(3); outputfile << " | ";
                outputfile.width(19); outputfile << "(km/s)";
                outputfile.width(3); outputfile << " | ";
                outputfile.width(19); outputfile << "(km/s)";
                outputfile.width(3); outputfile << " | ";
                outputfile.width(19); outputfile << "(N)";
                outputfile.width(3); outputfile << " | ";
                outputfile.width(19); outputfile << "(N)";
                outputfile.width(3); outputfile << " | ";
                outputfile.width(19); outputfile << "(N)";
                outputfile.width(3); outputfile << " | ";
                outputfile.width(17); outputfile << "throttle (0-1)";
                outputfile.width(3); outputfile << " | ";
                outputfile.width(14); outputfile << "(N)";
                outputfile.width(3); outputfile << " | ";
                outputfile.width(14); outputfile << "(s)";
                outputfile.width(3); outputfile << " | ";
                outputfile.width(14); outputfile << "(kW)";
                outputfile.width(3); outputfile << " | ";
                outputfile.width(19); outputfile << "rate (kg/s)";
                outputfile.width(3); outputfile << " | ";
                outputfile.width(14); outputfile << "(kg)";
                outputfile.width(3); outputfile << " | ";
                outputfile.width(14); outputfile << "active engines";
                outputfile.width(3); outputfile << " | ";
                outputfile.width(14); outputfile << "(kW)";
                outputfile.width(3); outputfile << " | ";
                outputfile.width(14); outputfile << "";
                outputfile.width(3); outputfile << " | ";
                outputfile << std::endl;


                for (size_t k = 0; k < 627; ++k)
                    outputfile << "-";
                outputfile << std::endl;
            }

            //Step 3: separation event
            {
                phase::write_output_line(outputfile,//outputfile
                    eventcount,//eventcount
                    "separation",//event_type
                    "deep-space",//event_location
                    0.0,// timestep_size,
                    -1,//flyby_altitude,
                    0,//BdotR
                    0,//BdotT
                    0,//angle1
                    0,//angle2
                    0,//C3
                    this->probe_state_after_separation,//state
                    this->probe_separation_deltaV,//dV
                    math::Matrix<doubleType>(3, 1, 0.0),//ThrustVector
                    this->probe_separation_deltaV.norm(),//dVmag
                    0.0,//Thrust
                    -1,//Isp
                    0.0,//AvailPower
                    0.0,//mdot
                    0,//number_of_active_engines
                    0.0,
                    "none");//active_power)
            }

            //Step 4: forward half-phase to AEI
            {
                this->probe_separation_to_AEI_ForwardHalfPhasePropagator->setStateRight(output_state);

                for (size_t step = 0; step < this->myOptions->num_timesteps / 2; ++step)
                {
                    //propagate to the halfway point of this step
                    this->probe_separation_to_AEI_Forward_dPropagatedStatedIndependentVariable.assign_zeros();
                    this->probe_separation_to_AEI_ForwardHalfPhasePropagator->setCurrentEpoch(this->probe_state_after_separation(7));
                    this->probe_separation_to_AEI_ForwardHalfPhasePropagator->setIndexOfEpochInStateVec(7);
                    this->probe_separation_to_AEI_ForwardHalfPhasePropagator->setCurrentIndependentVariable(this->probe_state_after_separation(7));
                    this->probe_separation_to_AEI_ForwardHalfPhasePropagator->propagate(ForwardOutputTimestep_separation_to_AEI * (step + 0.5), false);

                    if (this->isKeplerian)
                    {
                        output_state(6) = this->probe_state_after_separation(6);
                        output_state(7) = this->probe_state_after_separation(7) + ForwardOutputTimestep_separation_to_AEI * (step + 0.5);
                    }

                    //print
                    this->write_output_line(outputfile,//outputfile
                        eventcount,//eventcount
                        "coast",//event_type
                        "deep-space",//event_location
                        ForwardOutputTimestep_separation_to_AEI / 86400.0,// timestep_size,
                        -1,//flyby_altitude,
                        0,//BdotR
                        0,//BdotT
                        0,//angle1
                        0,//angle2
                        0,//C3
                        output_state,//state
                        empty3,//delta-v
                        empty3,//ThrustVector
                        0.0,//dVmag
                        0.0,//Thrust
                        0.0,//Isp
                        0.0,//AvailPower
                        0.0,//mdot
                        0,//number_of_active_engines
                        0.0, //active_power
                        "none");
                }

                //reset the propagator to go where it is supposed to
                this->probe_separation_to_AEI_ForwardHalfPhasePropagator->setStateRight(this->probe_first_subphase_match_point_state_minus);
            }

            //Step 5: match point to AEI
            {
                phase::write_output_line(outputfile,//outputfile
                    eventcount,//eventcount
                    "match_point",//event_type
                    "deep-space",//event_location
                    0.0,// timestep_size,
                    -1,//flyby_altitude,
                    0,//BdotR
                    0,//BdotT
                    0,//angle1
                    0,//angle2
                    0,//C3
                    this->probe_first_subphase_match_point_state_minus,//state
                    math::Matrix<doubleType>(3, 1, 0.0),//dV
                    math::Matrix<doubleType>(3, 1, 0.0),//ThrustVector
                    0.0,//dVmag
                    0.0,//Thrust
                    0.0,//Isp
                    0.0,//AvailPower
                    0.0,//mdot
                    0,//number_of_active_engines
                    0.0,
                    "none");//active_power
            }

            //Step 6: backward half-phase to AEI
            {
                //set the propagator temporarily to dump into the output state
                this->probe_separation_to_AEI_BackwardHalfPhasePropagator->setStateRight(output_state);

                for (size_t step = 0; step < this->myOptions->num_timesteps / 2; ++step)
                {
                    size_t backstep = this->myOptions->num_timesteps / 2 - step - 1;

                    //propagate to the halfway point of this step
                    this->probe_separation_to_AEI_Backward_dPropagatedStatedIndependentVariable.assign_zeros();
                    this->probe_separation_to_AEI_BackwardHalfPhasePropagator->setCurrentEpoch(this->probe_state_at_AEI_arrival(7));
                    this->probe_separation_to_AEI_BackwardHalfPhasePropagator->setIndexOfEpochInStateVec(7);
                    this->probe_separation_to_AEI_BackwardHalfPhasePropagator->setCurrentIndependentVariable(this->probe_state_at_AEI_arrival(7));
                    this->probe_separation_to_AEI_BackwardHalfPhasePropagator->propagate(-BackwardOutputTimestep_separation_to_AEI * (backstep + 0.5), false);

                    if (this->isKeplerian)
                    {
                        output_state(6) = this->probe_state_at_AEI_arrival(6);
                        output_state(7) = this->probe_state_at_AEI_arrival(7) - BackwardOutputTimestep_separation_to_AEI * (backstep + 0.5);
                    }

                    //print
                    this->write_output_line(outputfile,//outputfile
                        eventcount,//eventcount
                        "coast",//event_type
                        "deep-space",//event_location
                        BackwardOutputTimestep_separation_to_AEI / 86400.0,// timestep_size,
                        -1,//flyby_altitude,
                        0,//BdotR
                        0,//BdotT
                        0,//angle1
                        0,//angle2
                        0,//C3
                        output_state,//state
                        empty3,//delta-v
                        empty3,//ThrustVector
                        0.0,//dVmag
                        0.0,//Thrust
                        0.0,//Isp
                        0.0,//AvailPower
                        0.0,//mdot
                        0,//number_of_active_engines
                        0.0,
                        "none");//active_power
                }
                
                //reset the propagator to go where it is supposed to
                this->probe_separation_to_AEI_BackwardHalfPhasePropagator->setStateRight(this->probe_first_subphase_match_point_state_plus);
            }

            //Step 7: probe AEI arrival event
            this->probeEntryArrivalEvent->output(outputfile, this->LaunchDate _GETVALUE, eventcount);

            //Step 8: forward half phase for descent
            {
                this->probe_AEI_to_end_ForwardHalfPhasePropagator->setStateRight(output_state);

                for (size_t step = 0; step < this->myOptions->num_timesteps / 2; ++step)
                {
                    //propagate to the halfway point of this step
                    this->probe_AEI_to_end_Forward_dPropagatedStatedIndependentVariable.assign_zeros();
                    this->probe_AEI_to_end_ForwardHalfPhasePropagator->setCurrentEpoch(this->probe_state_at_AEI_departure(7));
                    this->probe_AEI_to_end_ForwardHalfPhasePropagator->setIndexOfEpochInStateVec(7);
                    this->probe_AEI_to_end_ForwardHalfPhasePropagator->setCurrentIndependentVariable(this->probe_state_at_AEI_departure(7));
                    this->probe_AEI_to_end_ForwardHalfPhasePropagator->propagate(ForwardOutputTimestep_AEI_to_end * (step + 0.5), false);

                    if (this->isKeplerian)
                    {
                        output_state(6) = this->probe_state_at_AEI_departure(6);
                        output_state(7) = this->probe_state_at_AEI_departure(7) + ForwardOutputTimestep_AEI_to_end * (step + 0.5);
                    }

                    //print
                    this->write_output_line(outputfile,//outputfile
                        eventcount,//eventcount
                        "coast",//event_type
                        "deep-space",//event_location
                        ForwardOutputTimestep_AEI_to_end / 86400.0,// timestep_size,
                        -1,//flyby_altitude,
                        0,//BdotR
                        0,//BdotT
                        0,//angle1
                        0,//angle2
                        0,//C3
                        output_state,//state
                        empty3,//delta-v
                        empty3,//ThrustVector
                        0.0,//dVmag
                        0.0,//Thrust
                        0.0,//Isp
                        0.0,//AvailPower
                        0.0,//mdot
                        0,//number_of_active_engines
                        0.0, //active_power
                        "none");
                }

                //reset the propagator to go where it is supposed to
                this->probe_AEI_to_end_ForwardHalfPhasePropagator->setStateRight(this->probe_second_subphase_match_point_state_minus);
            }

            //Step 9: descent match point
            {
                phase::write_output_line(outputfile,//outputfile
                    eventcount,//eventcount
                    "match_point",//event_type
                    "deep-space",//event_location
                    0.0,// timestep_size,
                    -1,//flyby_altitude,
                    0,//BdotR
                    0,//BdotT
                    0,//angle1
                    0,//angle2
                    0,//C3
                    this->probe_second_subphase_match_point_state_minus,//state
                    math::Matrix<doubleType>(3, 1, 0.0),//dV
                    math::Matrix<doubleType>(3, 1, 0.0),//ThrustVector
                    0.0,//dVmag
                    0.0,//Thrust
                    0.0,//Isp
                    0.0,//AvailPower
                    0.0,//mdot
                    0,//number_of_active_engines
                    0.0,
                    "none");//active_power
            }

            //Step 10: backward half phase for descent
            {
                //set the propagator temporarily to dump into the output state
                this->probe_AEI_to_end_BackwardHalfPhasePropagator->setStateRight(output_state);

                for (size_t step = 0; step < this->myOptions->num_timesteps / 2; ++step)
                {
                    size_t backstep = this->myOptions->num_timesteps / 2 - step - 1;

                    //propagate to the halfway point of this step
                    this->probe_AEI_to_end_Backward_dPropagatedStatedIndependentVariable.assign_zeros();
                    this->probe_AEI_to_end_BackwardHalfPhasePropagator->setCurrentEpoch(this->probe_state_at_end_of_phase(7));
                    this->probe_AEI_to_end_BackwardHalfPhasePropagator->setIndexOfEpochInStateVec(7);
                    this->probe_AEI_to_end_BackwardHalfPhasePropagator->setCurrentIndependentVariable(this->probe_state_at_end_of_phase(7));
                    this->probe_AEI_to_end_BackwardHalfPhasePropagator->propagate(-BackwardOutputTimestep_AEI_to_end * (backstep + 0.5), false);

                    if (this->isKeplerian)
                    {
                        output_state(6) = this->probe_state_at_end_of_phase(6);
                        output_state(7) = this->probe_state_at_end_of_phase(7) - BackwardOutputTimestep_AEI_to_end * (backstep + 0.5);
                    }

                    //print
                    this->write_output_line(outputfile,//outputfile
                        eventcount,//eventcount
                        "coast",//event_type
                        "deep-space",//event_location
                        BackwardOutputTimestep_AEI_to_end / 86400.0,// timestep_size,
                        -1,//flyby_altitude,
                        0,//BdotR
                        0,//BdotT
                        0,//angle1
                        0,//angle2
                        0,//C3
                        output_state,//state
                        empty3,//delta-v
                        empty3,//ThrustVector
                        0.0,//dVmag
                        0.0,//Thrust
                        0.0,//Isp
                        0.0,//AvailPower
                        0.0,//mdot
                        0,//number_of_active_engines
                        0.0,
                        "none");//active_power
                }

                //reset the propagator to go where it is supposed to
                this->probe_AEI_to_end_BackwardHalfPhasePropagator->setStateRight(this->probe_second_subphase_match_point_state_plus);
            }

            //Step 11: probe end event
            this->probeEndEvent->output(outputfile, this->LaunchDate _GETVALUE, eventcount);

            //Step 12: probe boundary states
            {
                outputfile << std::endl;

                //create a 3-element storage vector that will be used every time something needs to be rotated to the local frame
                math::Matrix<doubleType> r(3, 1);
                math::Matrix<doubleType> v(3, 1);
                math::Matrix<doubleType> disp_r(3, 1);
                math::Matrix<doubleType> disp_v(3, 1);

                //construct the rotation matrix from ICRF to the Universe J2000 BCI frame
                outputfile << "Boundary states:" << std::endl;
                outputfile.precision(8);

                //left boundary
                outputfile << "Boundary 1: ";
                r = this->state_after_initial_TCM.getSubMatrix1D(0, 2);
                v = this->state_after_initial_TCM.getSubMatrix1D(3, 5);

                this->myUniverse->LocalFrame.rotate_frame_to_frame(ReferenceFrame::ICRF, r, this->myOptions->output_file_frame, disp_r);
                this->myUniverse->LocalFrame.rotate_frame_to_frame(ReferenceFrame::ICRF, v, this->myOptions->output_file_frame, disp_v);

                for (int j = 0; j < 3; ++j)
                    outputfile << " " << disp_r(j) _GETVALUE;
                for (int j = 0; j < 3; ++j)
                    outputfile << " " << disp_v(j) _GETVALUE;
                outputfile << std::endl;

                //AEI
                outputfile << "Boundary 2: ";
                r = this->probe_state_at_AEI_arrival.getSubMatrix1D(0, 2);
                v = this->probe_state_at_AEI_arrival.getSubMatrix1D(3, 5);

                this->myUniverse->LocalFrame.rotate_frame_to_frame(ReferenceFrame::ICRF, r, this->myOptions->output_file_frame, disp_r);
                this->myUniverse->LocalFrame.rotate_frame_to_frame(ReferenceFrame::ICRF, v, this->myOptions->output_file_frame, disp_v);

                for (int j = 0; j < 3; ++j)
                    outputfile << " " << disp_r(j) _GETVALUE;
                for (int j = 0; j < 3; ++j)
                    outputfile << " " << disp_v(j) _GETVALUE;
                outputfile << std::endl;

                //End
                outputfile << "Boundary 3: ";
                r = this->probe_state_at_end_of_phase.getSubMatrix1D(0, 2);
                v = this->probe_state_at_end_of_phase.getSubMatrix1D(3, 5);

                this->myUniverse->LocalFrame.rotate_frame_to_frame(ReferenceFrame::ICRF, r, this->myOptions->output_file_frame, disp_r);
                this->myUniverse->LocalFrame.rotate_frame_to_frame(ReferenceFrame::ICRF, v, this->myOptions->output_file_frame, disp_v);

                for (int j = 0; j < 3; ++j)
                    outputfile << " " << disp_r(j) _GETVALUE;
                for (int j = 0; j < 3; ++j)
                    outputfile << " " << disp_v(j) _GETVALUE;
                outputfile << std::endl;
            }

            //Step 13: probe constraints
            outputfile << std::endl;
            this->probeEntryArrivalEvent->output_specialized_constraints(outputfile);
            this->probeEntryDepartureEvent->output_specialized_constraints(outputfile);
            this->probeEndEvent->output_specialized_constraints(outputfile);

            //Step 14: close the file
            outputfile.close();
        }//end output_probe

        void ProbeEntryPhase::output_ephemeris(std::ofstream& outputfile, std::ofstream& acceleration_model_file)
        {
            //Step 1: output an ephemeris for the spacecraft
            this->MGAnDSMs_phase::output_ephemeris(outputfile, acceleration_model_file);

            //Step 2: output an ephemeris file for the probe
            this->output_probe_ephemeris();
        }//end output_ephemeris()

        void ProbeEntryPhase::output_probe_ephemeris()
        {            
            std::ofstream probe_acceleration_model_file;
            const size_t bufsize = 128 * 1024;
            char buf1[bufsize];
            char buf2[bufsize];                                    

            // open a new ephemeris file for the probe
            if (this->myOptions->generate_acceleration_model_instrumentation_file)
            {
                std::string acceleration_model_string = this->myOptions->working_directory + "//" + this->myOptions->mission_name + "_probe.accelerations";
                probe_acceleration_model_file.open(acceleration_model_string, std::ios::out | std::ios::trunc);

                probe_acceleration_model_file.rdbuf()->pubsetbuf(buf1, bufsize);

                probe_acceleration_model_file.setf(std::ios::scientific, std::ios::floatfield);
                probe_acceleration_model_file.precision(16);
            }
            
            // headers
            // ephemeris file header
            std::ofstream ephemeris_file(this->myOptions->working_directory + "//" + this->myOptions->mission_name + "_probe.ephemeris", std::ios::out | std::ios::trunc);
            ephemeris_file.rdbuf()->pubsetbuf(buf2, bufsize);
            ephemeris_file << "#epoch, x(km), y(km), z(km), vx(km/s), vy(km/s), vz(km/s)" << std::endl;
            
            // acceleration file header
            if (this->myOptions->generate_acceleration_model_instrumentation_file)
            {
                if (journeyIndex > 0)
                {
                    probe_acceleration_model_file << "\n";
                }
                probe_acceleration_model_file << "#epoch (JD)";
                probe_acceleration_model_file << ", sc_distance_from_cb" << ", x" << ", y" << ", z";
                probe_acceleration_model_file << ", velocity_norm" << ", vx" << ", vy" << ", vz";
                probe_acceleration_model_file << ", mass";
                probe_acceleration_model_file << ", acceleration_norm" << ", ax" << ", ay" << ", az";;
                if (this->myOptions->perturb_SRP)
                {
                    probe_acceleration_model_file << ", SRP_norm, SRP_x, SRP_y, SRP_z";
                }
                for (size_t k = 0; k < this->myOptions->Journeys.size(); ++k)
                {
                    if (this->myOptions->Journeys[k].perturb_drag)
                    {
                        probe_acceleration_model_file << ", Drag_norm, Drag_x, Drag_y, Drag_z";
                        break;
                    }
                }
                probe_acceleration_model_file << ", cb_point_mass_gravity_norm";
                probe_acceleration_model_file << ", cb_point_mass_gravity_x";
                probe_acceleration_model_file << ", cb_point_mass_gravity_y";
                probe_acceleration_model_file << ", cb_point_mass_gravity_z";
                if (this->myOptions->perturb_J2)
                {
                    probe_acceleration_model_file << ", cb_J2_gravity_norm";
                    probe_acceleration_model_file << ", cb_J2_gravity_x";
                    probe_acceleration_model_file << ", cb_J2_gravity_y";
                    probe_acceleration_model_file << ", cb_J2_gravity_z";
                }
                if (this->myOptions->perturb_thirdbody)
                {
                    for (size_t k = 0; k < this->myOptions->Journeys[journeyIndex].perturbation_bodies.size(); ++k)
                    {
                        probe_acceleration_model_file << ", " + this->myUniverse->bodies[k].name + "_distance_from_cb";
                        probe_acceleration_model_file << ", " + this->myUniverse->bodies[k].name + "_x";
                        probe_acceleration_model_file << ", " + this->myUniverse->bodies[k].name + "_y";
                        probe_acceleration_model_file << ", " + this->myUniverse->bodies[k].name + "_z";
                        probe_acceleration_model_file << ", " + this->myUniverse->bodies[k].name + "_speed_wrt_cb";
                        probe_acceleration_model_file << ", " + this->myUniverse->bodies[k].name + "_vx";
                        probe_acceleration_model_file << ", " + this->myUniverse->bodies[k].name + "_vy";
                        probe_acceleration_model_file << ", " + this->myUniverse->bodies[k].name + "_vz";
                        probe_acceleration_model_file << ", sc_distance_from_" + this->myUniverse->bodies[k].name;
                        probe_acceleration_model_file << ", " + this->myUniverse->bodies[k].name + "_2_sc_x";
                        probe_acceleration_model_file << ", " + this->myUniverse->bodies[k].name + "_2_sc_y";
                        probe_acceleration_model_file << ", " + this->myUniverse->bodies[k].name + "_2_sc_z";
                        probe_acceleration_model_file << ", sc_speed_wrt_" + this->myUniverse->bodies[k].name;
                        probe_acceleration_model_file << ", " + this->myUniverse->bodies[k].name + "_2_sc_vx";
                        probe_acceleration_model_file << ", " + this->myUniverse->bodies[k].name + "_2_sc_vy";
                        probe_acceleration_model_file << ", " + this->myUniverse->bodies[k].name + "_2_sc_vz";
                        probe_acceleration_model_file << ", " + this->myUniverse->bodies[k].name + "_point_mass_gravity_norm";
                        probe_acceleration_model_file << ", " + this->myUniverse->bodies[k].name + "_point_mass_gravity_x";
                        probe_acceleration_model_file << ", " + this->myUniverse->bodies[k].name + "_point_mass_gravity_y";
                        probe_acceleration_model_file << ", " + this->myUniverse->bodies[k].name + "_point_mass_gravity_z";
                    }
                }

                probe_acceleration_model_file.flush();
            }

            // write the cmd file for the ephemeris
            {
                std::string cmdstring = this->myOptions->working_directory + "//" + this->myOptions->mission_name + "_probe.cmd";
                std::ofstream cmd_file(cmdstring, std::ios::out | std::ios::trunc);

                cmd_file << "\\begindata" << std::endl;
                cmd_file << "INPUT_DATA_TYPE = 'STATES'" << std::endl;
                cmd_file << "OUTPUT_SPK_TYPE = 9" << std::endl;
                cmd_file << "OBJECT_ID = " << this->myOptions->spacecraft_SPICE_ID + 1 << std::endl;
                cmd_file << "OBJECT_NAME = '" << this->myOptions->mission_name << "_probe'" << std::endl;
                cmd_file << "CENTER_ID = " << this->myOptions->forward_integrated_ephemeris_central_body_SPICE_ID << std::endl;
                cmd_file << "CENTER_NAME = 'CENTER'" << std::endl;
                cmd_file << "REF_FRAME_NAME = 'J2000'" << std::endl;
                cmd_file << "PRODUCER_ID = 'EMTGv9'" << std::endl;
                cmd_file << "DATA_ORDER = 'EPOCH X Y Z VX VY VZ'" << std::endl;
                cmd_file << "INPUT_DATA_UNITS = ('ANGLES=DEGREES' 'DISTANCES=km')" << std::endl;
                cmd_file << "TIME_WRAPPER = '# TDB'" << std::endl;
                cmd_file << "DATA_DELIMITER = ','" << std::endl;
                cmd_file << "LINES_PER_RECORD = 1" << std::endl;
                cmd_file << "IGNORE_FIRST_LINE = 0" << std::endl;
                cmd_file << "LEAPSECONDS_FILE = '" << boost::replace_all_copy(this->myOptions->universe_folder, "\\", "/") + "/ephemeris_files/" + this->myOptions->SPICE_leap_seconds_kernel << "'" << std::endl;
                cmd_file << "POLYNOM_DEGREE = 3" << std::endl;
                cmd_file << "SEGMENT_ID = 'SPK_STATES_09'" << std::endl;
                cmd_file << "\\begintext" << std::endl;

                cmd_file.close();

                //Step 3: write a Python file that makes and checks the .bsp
                std::ofstream pyfile(this->myOptions->working_directory + "//bspwriter_probe.py", std::ios::out | std::ios::trunc);

                pyfile << "import os" << std::endl;
                if (strcmp(this->myOptions->spice_utility_extension.c_str(), "\"\""))
                {
                    pyfile << "mkspk_path = '" << boost::replace_all_copy(this->myOptions->spice_utilities_path, "\\", "/") << "/mkspk'" << std::endl;
                    pyfile << "brief_path = '" << boost::replace_all_copy(this->myOptions->spice_utilities_path, "\\", "/") << "/brief'" << std::endl;
                }
                else
                {
                    pyfile << "mkspk_path = '" << boost::replace_all_copy(this->myOptions->spice_utilities_path, "\\", "/") << "/mkspk" << this->myOptions->spice_utility_extension << "'" << std::endl;
                    pyfile << "brief_path = '" << boost::replace_all_copy(this->myOptions->spice_utilities_path, "\\", "/") << "/brief" << this->myOptions->spice_utility_extension << "'" << std::endl;
                }

                pyfile << "import sys " << std::endl;
                pyfile << "sys.path.append('" << boost::replace_all_copy(this->myOptions->pyemtg_path, "\\", "/") << "')" << std::endl;
                pyfile << "sys.path.append('" << boost::replace_all_copy(this->myOptions->pyemtg_path, "\\", "/") << "/SimpleMonteCarlo')" << std::endl;
                pyfile << "sys.path.append('" << boost::replace_all_copy(this->myOptions->pyemtg_path, "\\", "/") << "/SpiceyPy_Utilities') " << std::endl;

                pyfile << "import clean_spiceicles " << std::endl;

                pyfile << "clean_spiceicles.do_the_stuff(['" << boost::replace_all_copy(this->myOptions->working_directory, "\\", "/") << "','" << this->myOptions->mission_name << "_probe.ephemeris','" << boost::replace_all_copy(this->myOptions->universe_folder, "\\", "/") << "/ephemeris_files/']) " << std::endl;


                pyfile << "if os.path.exists('" << this->myOptions->mission_name << "_probe.bsp') :" << std::endl;
                pyfile << "    os.remove('" << this->myOptions->mission_name << "_probe.bsp')" << std::endl;
                pyfile << "os.system(mkspk_path + ' -setup " << this->myOptions->mission_name << "_probe.cmd -input " << this->myOptions->mission_name << "_probe_clean.ephemeris -output " << this->myOptions->mission_name << "_probe.bsp')" << std::endl;
                pyfile << "os.system(brief_path + ' " << this->myOptions->mission_name << "_probe.bsp')" << std::endl;

                pyfile.close();
            }
            
            // output the departure event
            this->myDepartureEvent->output_ephemeris(ephemeris_file);

            if (this->myOptions->generate_acceleration_model_instrumentation_file && !this->isKeplerian)
            {
                this->temp_state = this->myDepartureEvent->get_state_before_event();
                this->probe_separation_to_AEI_SpacecraftAccelerationModel->populateInstrumentationFile(probe_acceleration_model_file, this->temp_state, this->temp_state(7));
            }

            // output the probe's trajectory
            {
                // we'll need an output vector
                this->temp_state = this->output_state;
                // temporarily assign the propagator to the output state
                this->probe_separation_to_AEI_ForwardHalfPhasePropagator->setStateLeft(this->temp_state);
                this->probe_separation_to_AEI_ForwardHalfPhasePropagator->setStateRight(output_state);
                this->probe_separation_to_AEI_ForwardHalfPhasePropagator->setIndexOfEpochInStateVec(7);
                this->temp_state = this->probe_state_after_separation;

                // propagate and print, skipping the first entry
                double timeToPropagate = this->EphemerisOutputResolution;
                double propagated_time = 0;
                while (propagated_time < this->PhaseFlightTime)
                {
                    // propagate
                    this->probe_separation_to_AEI_Forward_dPropagatedStatedIndependentVariable.assign_zeros();
                    this->probe_separation_to_AEI_ForwardHalfPhasePropagator->setCurrentEpoch(this->temp_state(7));
                    this->probe_separation_to_AEI_ForwardHalfPhasePropagator->setCurrentIndependentVariable(this->temp_state(7));
                    this->probe_separation_to_AEI_ForwardHalfPhasePropagator->propagate(timeToPropagate, false);
                    this->temp_state = output_state;
                    propagated_time += timeToPropagate;

                    if (this->isKeplerian)
                    {
                        output_state(6) = this->probe_state_after_separation(6);
                        output_state(7) = this->probe_state_after_separation(7) + propagated_time;
                    }

                    // convert to Sun-centered if necessary
                    if (this->myUniverse->central_body.spice_ID != this->myOptions->forward_integrated_ephemeris_central_body_SPICE_ID)
                    {
                        double LT_dump;
                        double bodyStateDouble[6];
                        spkez_c(this->myUniverse->central_body.spice_ID, output_state(7)_GETVALUE - (51544.5 * 86400.0), "J2000", "NONE", this->myOptions->forward_integrated_ephemeris_central_body_SPICE_ID, bodyStateDouble, &LT_dump);

                        for (size_t stateIndex = 0; stateIndex < 6; ++stateIndex)
                            output_state(stateIndex) += bodyStateDouble[stateIndex];
                    }

                    // populate probe acceleration model file if applicable
                    if (this->myOptions->generate_acceleration_model_instrumentation_file && !this->isKeplerian)
                    {
                        this->probe_separation_to_AEI_SpacecraftAccelerationModel->populateInstrumentationFile(probe_acceleration_model_file, this->temp_state, this->temp_state(7));
                    }

                    // write line to ephemeris file
                    this->write_ephemeris_line(ephemeris_file,
                        output_state,
                        math::Matrix<doubleType>(3, 1, 0.0),//control vector
                        0.0,
                        0.0,
                        0.0,
                        0,
                        0.0,
                        "none");

                    // increment propagatedEpoch
                    timeToPropagate = (this->PhaseFlightTime _GETVALUE - propagated_time) > this->EphemerisOutputResolution
                        ? this->EphemerisOutputResolution
                        : (this->PhaseFlightTime _GETVALUE - propagated_time);
                }

                // reset the propagator to its original state vectors
                this->probe_separation_to_AEI_ForwardHalfPhasePropagator->setStateLeft(this->probe_state_after_separation);
                this->probe_separation_to_AEI_ForwardHalfPhasePropagator->setStateRight(this->probe_first_subphase_match_point_state_minus);
            }

            // output the arrival event
            this->myArrivalEvent->output_ephemeris(ephemeris_file);
            this->temp_state = this->myArrivalEvent->get_state_before_event();
            if (this->myOptions->generate_acceleration_model_instrumentation_file && !this->isKeplerian)
            {
                this->probe_separation_to_AEI_SpacecraftAccelerationModel->populateInstrumentationFile(probe_acceleration_model_file, this->temp_state, this->temp_state(7));
            }

            //Step 6: close the files
            ephemeris_file.close();
            probe_acceleration_model_file.close();
        }//end output_probe_ephemeris()
        
        void ProbeEntryPhase::output_STMs()
        {
            //Step 1: spacecraft
            this->MGAnDSMs_phase::output_STMs();

            //Step 2: probe
            this->probe_separation_to_AEI_ForwardSTM.print_to_file(this->myOptions->working_directory + "//" + boost::replace_all_copy(this->prefix, ": ", "") + "_probe_separation_to_AEI_STM_0_forward.stm");
            this->probe_separation_to_AEI_BackwardSTM.print_to_file(this->myOptions->working_directory + "//" + boost::replace_all_copy(this->prefix, ": ", "") + "_probe_separation_to_AEI_STM_1_backward.stm");
            this->probe_AEI_to_end_ForwardSTM.print_to_file(this->myOptions->working_directory + "//" + boost::replace_all_copy(this->prefix, ": ", "") + "_probe_AEI_to_end_STM_0_forward.stm");
            this->probe_AEI_to_end_BackwardSTM.print_to_file(this->myOptions->working_directory + "//" + boost::replace_all_copy(this->prefix, ": ", "") + "_probe_AEI_to_end_STM_1_backward.stm");
        }//end output_STMs()

        void ProbeEntryPhase::output_maneuver_and_target_spec(std::ofstream& maneuver_spec_file, std::ofstream& target_spec_file, bool& haveManeuverNeedTarget)
        {
            //Step 1: Reset the spacecraft
            this->mySpacecraft->setActiveStage(this->stageIndex);

            //Step 2: should write out a maneuver/target spec for the departure maneuver if there is one
            this->myDepartureEvent->output_maneuver_and_target_spec(maneuver_spec_file, target_spec_file, haveManeuverNeedTarget);

            //Step 3: the target of the departure maneuver (or previous maneuver if there wasn't a departure maneuver) is the *probe entry state*
            if (haveManeuverNeedTarget)
            {
                //reset the flag that tells us we need a target spec line
                haveManeuverNeedTarget = false;

                //Step 1.1: configure target spec
                target_spec_line myTargetSpecLine(this->name,
                    "EME2000",
                    this->myUniverse->central_body.name,
                    this->probe_state_at_end_of_phase(7),
                    this->probe_state_at_end_of_phase);

                //Step 1.2: write target spec object
                myTargetSpecLine.write(target_spec_file);
            }

            //Step 4: forward maneuver and target space for the divert(s)
            for (size_t subPhaseIndex = 0; subPhaseIndex < this->numberOfForwardSubphases; ++subPhaseIndex)
            {
                this->ForwardSubPhases[subPhaseIndex].output_maneuver_and_target_spec(maneuver_spec_file, target_spec_file, haveManeuverNeedTarget);
            }

            //Step 5: backward maneuver and target spec for the divert(s)
            for (size_t subPhaseIndex = 0; subPhaseIndex < this->numberOfBackwardSubphases; ++subPhaseIndex)
            {
                size_t backSubPhaseIndex = this->numberOfBackwardSubphases - subPhaseIndex - 1;
                this->BackwardSubPhases[backSubPhaseIndex].output_maneuver_and_target_spec(maneuver_spec_file, target_spec_file, haveManeuverNeedTarget);
            }

            //should write out a maneuver/target spec for the arrival maneuver if there is one
            this->myArrivalEvent->output_maneuver_and_target_spec(maneuver_spec_file, target_spec_file, haveManeuverNeedTarget);
        }//end output_maneuver_and_target_spec()
    }//close namespace Phases
}//close namespace EMTG