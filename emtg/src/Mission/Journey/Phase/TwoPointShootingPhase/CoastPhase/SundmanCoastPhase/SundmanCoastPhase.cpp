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

//EMTGv9 SundmanCoastPhase
//Jacob Englander 4-5-2018

#include "SundmanCoastPhase.h"

#include "PropagatorFactory.h"
#include "IntegrationSchemeFactory.h"


namespace EMTG
{
    namespace Phases
    {
        SundmanCoastPhase::SundmanCoastPhase(const std::string& name,
            const size_t& journeyIndex,
            const size_t& phaseIndex,
            size_t& stageIndex,
            phase* previousPhase,
            Astrodynamics::universe* Universe,
            HardwareModels::Spacecraft* mySpacecraft,
            HardwareModels::LaunchVehicle* myLaunchVehicle,
            missionoptions* myOptions) :
            CoastPhase::CoastPhase(name,
                journeyIndex,
                phaseIndex,
                stageIndex,
                previousPhase,
                Universe,
                mySpacecraft,
                myLaunchVehicle,
                myOptions,
                9, //numStatesToPropagate
                9) //numMatchConstraints
        {
            this->G_indices_match_point_constraints_wrt_SundmanIndependentVariable.resize(this->numMatchConstraints);

            this->initialize();
        }//end constructor

        void SundmanCoastPhase::initialize()
        {
            this->num_timesteps = this->myJourneyOptions->override_num_steps
                ? this->myJourneyOptions->number_of_steps
                : this->myOptions->num_timesteps;

            this->MatchPointFraction = this->myJourneyOptions->CoastPhaseMatchPointFraction;
            this->dForwardStepSize_dPropagationVariable = this->MatchPointFraction;
            this->dBackwardStepSize_dPropagationVariable = 1.0 - this->MatchPointFraction;

            this->stateVectorNames.push_back("virtual chemical fuel");

            //name our match point constraints
            this->matchPointConstraintNames.push_back("virtual chemical fuel");
            this->matchPointConstraintStateIndex.push_back(8);
            this->continuity_constraint_scale_factors(7) = 1.0 / (1.0 * this->myJourneyOptions->maximum_mass);
            
            this->stateVectorNames.push_back("phase flight time");

            //match point constraint for time
            this->matchPointConstraintNames.push_back("phase flight time");
            this->matchPointConstraintStateIndex.push_back(7);
            this->continuity_constraint_scale_factors(8) = 1.0 / this->myUniverse->TU;

            this->stateVectorNames.push_back("phase Sundman independent variable");
            this->stateIndex_phase_propagation_variable = this->numStatesToPropagate;

            //no maneuvers
            this->hasBipropManeuver = false;
            this->hasMonopropManeuver = false;
            this->hasElectricManeuver = false;

            //size the TCM transition matrix
            this->TCMTM = math::Matrix<double>(this->numStatesToPropagate + 1, math::identity);

            //configure propagators
            this->dStepSize_dPropagationVariable = 0.5;

            //which propagator do we want?
            if (this->myPropagatorType == PropagatorType::KeplerPropagator)
            {
                throw std::invalid_argument("Keplerian propagation is not possible with SundmanCoastPhase, use CoastPhase instead...aborting");
            }
            else //integrated propagator
            {
                this->isKeplerian = false;
                this->ForwardSTM = math::Matrix<double>(10, math::identity);
                this->BackwardSTM = math::Matrix<double>(10, math::identity);
                this->Forward_dPropagatedStatedIndependentVariable = math::Matrix<double>(this->numStatesToPropagate, 2, 0.0);
                this->Backward_dPropagatedStatedIndependentVariable = math::Matrix<double>(this->numStatesToPropagate, 2, 0.0);

                //acceleration model object
                this->mySpacecraftAccelerationModel = new Astrodynamics::SpacecraftAccelerationModel(this->myOptions,
                    this->myJourneyOptions,
                    this->myUniverse,
                    this->Xdescriptions,
                    this->mySpacecraft,
                    10); // STM size
                this->mySpacecraftAccelerationModel->setDutyCycle(this->PhaseDutyCycle);

                //EOM
                this->myEOM.setSpacecraftAccelerationModel(this->mySpacecraftAccelerationModel);

                //integration scheme
                this->myIntegrationScheme = CreateIntegrationScheme(&this->myEOM, this->numStatesToPropagate, 10);

                //integration step length

                if (this->myJourneyOptions->override_integration_step_size)
                {
                    this->integration_step_length = this->myJourneyOptions->integration_step_size / 180.0 * math::PI * this->myUniverse->LU;
                }
                else
                {
                    this->integration_step_length = 0.25 / 180.0 * math::PI * this->myUniverse->LU;
                }

                this->ForwardHalfPhasePropagator = CreatePropagator(this->myOptions,
                    this->myUniverse,                 
                    this->numStatesToPropagate,
                    10,
                    this->state_after_initial_TCM,
                    this->match_point_state_minus,
                    this->ForwardSTM,
                    this->Forward_dPropagatedStatedIndependentVariable,
                    (Integration::Integrand*) &this->myEOM,
                    this->myIntegrationScheme,
                    &this->dForwardStepSize_dPropagationVariable,
                    this->integration_step_length);

                this->BackwardHalfPhasePropagator = CreatePropagator(this->myOptions,
                    this->myUniverse,
                    this->numStatesToPropagate,
                    10,
                    this->state_at_end_of_phase,
                    this->match_point_state_plus,
                    this->BackwardSTM,
                    this->Backward_dPropagatedStatedIndependentVariable,
                    (Integration::Integrand*) &this->myEOM,
                    this->myIntegrationScheme,
                    &this->dBackwardStepSize_dPropagationVariable,
                    this->integration_step_length);
            }

            //SPTMs
            this->ForwardSPTM = math::Matrix<double>(this->numStatesToPropagate + 1, math::identity);
            this->BackwardSPTM = math::Matrix<double>(this->numStatesToPropagate + 1, math::identity);

            //derivative truth table

            //unless SRP is enabled, only mass affects mass
            if (!this->myOptions->perturb_SRP)
            {
                //only mass variables affect mass constraints
                for (size_t stateIndex = 0; stateIndex < 6; ++stateIndex)
                {
                    this->TruthTable_Forward_MatchConstraints_Derivative_wrt_EncodedStates[6][stateIndex] = false;
                    this->TruthTable_Backward_MatchConstraints_Derivative_wrt_EncodedStates[6][stateIndex] = false;
                }
                this->TruthTable_Forward_MatchConstraints_Derivative_wrt_EncodedStates[6][7] = false; //epoch
                this->TruthTable_Backward_MatchConstraints_Derivative_wrt_EncodedStates[6][7] = false; //epoch
            }

            //always true - propellant tanks are only affected by themselves and mass
            for (size_t stateIndex = 0; stateIndex < 6; ++stateIndex)
            {
                this->TruthTable_Forward_MatchConstraints_Derivative_wrt_EncodedStates[7][stateIndex] = false; //fuel mass wrt state
                this->TruthTable_Backward_MatchConstraints_Derivative_wrt_EncodedStates[7][stateIndex] = false; //fuel mass wrt state
            }
            //backward propellant tanks are not affected by mass
            this->TruthTable_Backward_MatchConstraints_Derivative_wrt_EncodedStates[7][6] = false; //fuel mass wrt mass

            this->TruthTable_Forward_MatchConstraints_Derivative_wrt_EncodedStates[7][7] = false; //fuel mass wrt epoch
            this->TruthTable_Backward_MatchConstraints_Derivative_wrt_EncodedStates[7][7] = false; //fuel mass wrt epoch
            //always true - propellant tank states do not affect anything but themselves
            for (size_t constraintIndex = 0; constraintIndex < 7; ++constraintIndex)
            {
                this->TruthTable_Forward_MatchConstraints_Derivative_wrt_EncodedStates[constraintIndex][8] = false; //other states wrt fuel mass
                this->TruthTable_Backward_MatchConstraints_Derivative_wrt_EncodedStates[constraintIndex][8] = false; //other states wrt fuel mass
            }

            //derivatives with respect to phase independent variable
            //turn off derivatives of mass and fuel wrt time unless ACS tracking is enabled, and turn off derivatives of oxidizer wrt time in all cases
            this->TruthTable_MatchConstraints_Derivative_wrt_SundmanIndependentVariable.resize(this->numMatchConstraints, true);
            if (!this->myOptions->trackACS)
            {
                this->TruthTable_MatchConstraints_Derivative_wrt_SundmanIndependentVariable[6] = false;

                this->TruthTable_MatchConstraints_Derivative_wrt_SundmanIndependentVariable[7] = false;
            }

            //output stuff
            this->output_state.resize(this->numStatesToPropagate + 13 * 13, 1, 0.0);
        }//end initialize

        //************************************calcbounds methods
        void SundmanCoastPhase::calcbounds_phase_flight_time()
        {
            //Step 1: create the independent variable
            this->Xlowerbounds->push_back(math::SMALL);
            this->Xupperbounds->push_back((math::PI - math::SMALL) * this->myUniverse->LU);
            this->X_scale_factors->push_back(1.0);
            this->Xdescriptions->push_back(this->prefix + "phase Sundman independent variable");
            this->Xindex_SundmanIndependentVariable = this->Xdescriptions->size() - 1;

            //Step 2: base class
            phase::calcbounds_phase_flight_time();
        }

        void SundmanCoastPhase::calcbounds_match_point_constraints()
        {
            //Step 1: base class
            this->TwoPointShootingPhase::calcbounds_match_point_constraints();
            
            //Step 2: virtual chemical fuel
            this->create_sparsity_entry(this->Findices_match_point_constraints[7],
                this->X_index_of_first_decision_variable_in_this_phase,
                true,
                this->prefix + "virtual chemical fuel",
                this->Gindex_dVirtualChemicalFuelConstraint_dVirtualChemicalFuel);

            //Step 3: sparsity pattern with respect to independent variable
            for (size_t constraintIndex = 0; constraintIndex < this->numMatchConstraints; ++constraintIndex)
            {
                if (this->TruthTable_MatchConstraints_Derivative_wrt_SundmanIndependentVariable[constraintIndex])
                {
                    this->create_sparsity_entry(this->Findices_match_point_constraints[constraintIndex],
                        this->Xindex_SundmanIndependentVariable,
                        this->G_indices_match_point_constraints_wrt_SundmanIndependentVariable[constraintIndex]);
                }
            }
        }//end calcbounds_match_point_constraints()

        //************************************process methods
        void SundmanCoastPhase::process_phase_flight_time(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //Step 1: propagation variable
            this->SundmanIndependentVariable = X[Xindex++];

            //Step 2: base class
            this->phase::process_phase_flight_time(X, Xindex, F, Findex, G, needG);
        }//end process_phase_flight_time()

        void SundmanCoastPhase::process_forward_half_phase(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //Step 1: set the tank states
            this->state_after_initial_TCM(8) = this->chemical_fuel_used; //TCM propellant, usually zero
            
            //Step 3: propagate
            this->Forward_dPropagatedStatedIndependentVariable.assign_zeros();
            this->ForwardHalfPhasePropagator->setCurrentEpoch(this->state_after_initial_TCM(7));
            this->ForwardHalfPhasePropagator->setIndexOfEpochInStateVec(7);
            this->ForwardHalfPhasePropagator->setCurrentIndependentVariable(0.0);
            this->ForwardHalfPhasePropagator->propagate(this->MatchPointFraction * this->SundmanIndependentVariable, needG);

            //MTM entries for ACS
            if (needG)
            {
                //Step 3: form the SPTM
                //Step 3.1: upper left 6x6 is the regular STM
                size_t stateMax = 10;
                for (size_t i = 0; i < stateMax; ++i)
                    for (size_t j = 0; j < stateMax; ++j)
                        this->ForwardSPTM(i, j) = this->ForwardSTM(i, j);
                //Step 3.2: turn the upper right 6 x 2 into the Phi_t terms developed by Lantoine
                for (size_t i = 0; i < stateMax; ++i)
                {

                    //this really ought to be cleaned up...
                    if (this->myPropagatorType == PropagatorType::KeplerPropagator)
                    {
                        this->ForwardSPTM(i, this->stateIndex_phase_propagation_variable) = this->Forward_dPropagatedStatedIndependentVariable(i);
                    }
                    else //integrated propagator
                    {
                        this->ForwardSPTM(i, this->stateIndex_phase_propagation_variable) = this->ForwardSTM(i, 9);
                    }
                }
                //this->ForwardSPTM(7, this->stateIndex_phase_propagation_variable) = this->Forward_dPropagatedStatedIndependentVariable(7, 1);
            }
        }//end process_forward_half_phase()

        void SundmanCoastPhase::process_backward_half_phase(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //Step 1: set the tank states
            this->state_at_end_of_phase(8) = this->virtual_chemical_fuel_used;

            //Step 3: call the backward propagator
            this->Backward_dPropagatedStatedIndependentVariable.assign_zeros();
            this->BackwardHalfPhasePropagator->setCurrentEpoch(this->state_at_end_of_phase(7));
            this->BackwardHalfPhasePropagator->setIndexOfEpochInStateVec(7);
            this->BackwardHalfPhasePropagator->setCurrentIndependentVariable(0.0);
            this->BackwardHalfPhasePropagator->propagate(-(1.0 - this->MatchPointFraction) * this->SundmanIndependentVariable, needG);
                                                                                     
            if (needG)
            {
                //Step 3: form the SPTM
                //Step 3.1: upper left 6x6 is the regular STM
                size_t stateMax = 10;
                for (size_t i = 0; i < stateMax; ++i)
                    for (size_t j = 0; j < stateMax; ++j)
                        this->BackwardSPTM(i, j) = this->BackwardSTM(i, j);
                //Step 3.2: turn the upper right 6 x 2 into the Phi_t terms developed by Lantoine
                for (size_t i = 0; i < stateMax; ++i)
                {
                    //this really ought to be cleaned up...
                    if (this->myPropagatorType == PropagatorType::KeplerPropagator)
                    {
                        this->BackwardSPTM(i, this->stateIndex_phase_propagation_variable) = this->Backward_dPropagatedStatedIndependentVariable(i);
                    }
                    else //integrated propagator
                    {
                        this->BackwardSPTM(i, this->stateIndex_phase_propagation_variable) = this->BackwardSTM(i, 9);
                    }
                    
                }
                //this->BackwardSPTM(7, this->stateIndex_phase_propagation_variable) = this->Backward_dPropagatedStatedIndependentVariable(7, 1);
            }
        }//end process_backward_half_phase()

        void SundmanCoastPhase::process_match_point_constraints(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //Step 1: contruct the SPTM chains and the forward and backward HPTM
            if (needG)
            {
                this->ForwardHPTM = this->ForwardSPTM * this->TCMTM;
                this->BackwardHPTM = this->BackwardSPTM;
            }

            //Step 2: call the base class
            TwoPointShootingPhase::process_match_point_constraints(X, Xindex, F, Findex, G, needG);


            if (needG)
            {
                math::Matrix<double> dStateNow_dDecisionVariable(this->numStatesToPropagate + 1, 1, 0.0);
                math::Matrix<double> dMatchPointState_dDecisionVariable(this->numStatesToPropagate + 1, 1, 0.0);

                //Step 3: Derivatives with respect to propagation time

                // NOTE: TwoPointShootingPhase is already computing the match point partials for PhaseFlightTime
                // and those computations are applicable for SundmanCoastPhase as well, so we don't need to perform the calculation
                // here as well

                /*
                dStateNow_dDecisionVariable.assign_zeros();
                
                // PhaseFlightTime does NOT impact the left hand side of the match point
                dStateNow_dDecisionVariable(7) = 0.0;

                //Forward
                dMatchPointState_dDecisionVariable = this->ForwardHPTM * dStateNow_dDecisionVariable;

                for (size_t constraintIndex = 0; constraintIndex < this->numMatchConstraints; ++constraintIndex)
                {
                    if (this->TruthTable_MatchConstraints_Derivative_wrt_PropagationTime[constraintIndex])
                    {
                        size_t stateIndex = this->matchPointConstraintStateIndex[constraintIndex];
                        size_t Gindex = this->G_indices_match_point_constraints_wrt_PhaseFlightTime[constraintIndex];
                        size_t Xindex = this->jGvar->operator[](Gindex);

                        G[Gindex] += -this->X_scale_factors->operator[](Xindex)
                            * dMatchPointState_dDecisionVariable(stateIndex)
                            * this->continuity_constraint_scale_factors(constraintIndex);
                    }
                }

                //Backward
                dStateNow_dDecisionVariable.assign_zeros();
                dStateNow_dDecisionVariable(7) = 1.0;
                dMatchPointState_dDecisionVariable = this->BackwardHPTM * dStateNow_dDecisionVariable;

                for (size_t constraintIndex = 0; constraintIndex < this->numMatchConstraints; ++constraintIndex)
                {
                    if (this->TruthTable_MatchConstraints_Derivative_wrt_PropagationTime[constraintIndex])
                    {
                        size_t stateIndex = this->matchPointConstraintStateIndex[constraintIndex];
                        size_t Gindex = this->G_indices_match_point_constraints_wrt_PhaseFlightTime[constraintIndex];
                        size_t Xindex = this->jGvar->operator[](Gindex);

                        G[Gindex] = this->X_scale_factors->operator[](Xindex)
                            * dMatchPointState_dDecisionVariable(stateIndex)
                            * this->continuity_constraint_scale_factors(constraintIndex);
                    }
                }
                */

                //Step 3: Derivatives with respect to Sundman variable
                dStateNow_dDecisionVariable.assign_zeros();
                dStateNow_dDecisionVariable(this->stateIndex_phase_propagation_variable) = 1.0;


                //Forward
                dMatchPointState_dDecisionVariable = this->ForwardHPTM * dStateNow_dDecisionVariable;

                for (size_t constraintIndex = 0; constraintIndex < this->numMatchConstraints; ++constraintIndex)
                {
                    if (this->TruthTable_MatchConstraints_Derivative_wrt_SundmanIndependentVariable[constraintIndex])
                    {
                        size_t stateIndex = this->matchPointConstraintStateIndex[constraintIndex];
                        size_t Gindex = this->G_indices_match_point_constraints_wrt_SundmanIndependentVariable[constraintIndex];
                        size_t Xindex = this->jGvar->operator[](Gindex);

                        G[Gindex] += -this->X_scale_factors->operator[](Xindex)
                            * dMatchPointState_dDecisionVariable(stateIndex)
                            * this->continuity_constraint_scale_factors(constraintIndex);
                    }
                }

                //Backward
                dMatchPointState_dDecisionVariable = this->BackwardHPTM * dStateNow_dDecisionVariable;

                for (size_t constraintIndex = 0; constraintIndex < this->numMatchConstraints; ++constraintIndex)
                {
                    if (this->TruthTable_MatchConstraints_Derivative_wrt_SundmanIndependentVariable[constraintIndex])
                    {
                        size_t stateIndex = this->matchPointConstraintStateIndex[constraintIndex];
                        size_t Gindex = this->G_indices_match_point_constraints_wrt_SundmanIndependentVariable[constraintIndex];
                        size_t Xindex = this->jGvar->operator[](Gindex);

                        G[Gindex] += this->X_scale_factors->operator[](Xindex)
                            * dMatchPointState_dDecisionVariable(stateIndex)
                            * this->continuity_constraint_scale_factors(constraintIndex);
                    }
                }

                //virtual chemical fuel
                {
                    size_t Gindex = this->Gindex_dVirtualChemicalFuelConstraint_dVirtualChemicalFuel;
                    size_t Xindex = this->jGvar->operator[](Gindex);

                    G[Gindex] = this->X_scale_factors->operator[](Xindex)
                        * this->continuity_constraint_scale_factors(7);
                }
            }//end derivatives
        }//end process_match_point_constraints()

        void SundmanCoastPhase::output(std::ofstream& outputfile,
            size_t& eventcount)
        {
            doubleType output_step = this->SundmanIndependentVariable / this->num_timesteps;
            math::Matrix<doubleType> empty3(3, 1, 0.0);
            doubleType previous_epoch;

            //Step 1: output the departure event
            this->myDepartureEvent->output(outputfile, this->LaunchDate _GETVALUE, eventcount);
            this->mySpacecraft->setActiveStage(this->stageIndex);


            //Step 2: output the initial TCM if applicable
            this->output_initial_TCM(outputfile, eventcount);

            //Step 3: print the forward half phase
            //set the propagator temporarily to dump into the output state
            this->ForwardHalfPhasePropagator->setStateRight(output_state);

            for (size_t step = 0; step < this->num_timesteps / 2; ++step)
            {
                //propagate to the halfway point of this step
                previous_epoch = step == 0 ? this->state_after_initial_TCM(7) : output_state(7);
                this->Forward_dPropagatedStatedIndependentVariable.assign_zeros();
                this->ForwardHalfPhasePropagator->setCurrentEpoch(this->state_after_initial_TCM(7));
                this->ForwardHalfPhasePropagator->setIndexOfEpochInStateVec(7);
                this->ForwardHalfPhasePropagator->setCurrentIndependentVariable(0.0);

                if (step == this->num_timesteps / 2 - 1)
                {
                    this->ForwardHalfPhasePropagator->setStorePropagationHistory(true);
                    this->ForwardHalfPhasePropagator->propagate(output_step * (step + 0.5), false);
                    this->ForwardHalfPhasePropagator->setStorePropagationHistory(false);
                }
                else
                {
                    this->ForwardHalfPhasePropagator->propagate(output_step * (step + 0.5), false);
                }

                math::Matrix<doubleType> R_sc_Sun(3, 1, 0.0);
                if (this->myUniverse->central_body_SPICE_ID == 10)
                {
                    R_sc_Sun = this->output_state.getSubMatrix1D(0, 2);
                }
                else
                {
                    //where is the central body relative to the sun?
                    doubleType central_body_state_and_derivatives[12];
                    this->myUniverse->locate_central_body(this->output_state(7),
                        central_body_state_and_derivatives,
                        *this->myOptions,
                        false);

                    math::Matrix<doubleType> R_CB_Sun(3, 1, 0.0);
                    for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
                    {
                        R_CB_Sun(stateIndex) = central_body_state_and_derivatives[stateIndex];
                    }

                    R_sc_Sun = this->output_state.getSubMatrix1D(0, 2) + R_CB_Sun;
                }
                doubleType r_sc_sun_AU = R_sc_Sun.norm() / this->myOptions->AU;
                this->mySpacecraft->computePowerState(r_sc_sun_AU, this->output_state(7));

                doubleType power = this->mySpacecraft->getAvailablePower();
                doubleType active_power = this->mySpacecraft->getEPActivePower();

                //print
                this->write_output_line(outputfile,//outputfile
                    eventcount,//eventcount
                    "coast",//event_type
                    "deep-space",//event_location
                    (output_state(7) - previous_epoch) / 86400.0,// timestep_size,
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
                    power,//AvailPower
                    0.0,//mdot
                    0,//number_of_active_engines
                    active_power,
                    "none");//active_power
            } // end forward output

            //reset the propagator to go where it is supposed to
            this->ForwardHalfPhasePropagator->setStateRight(this->match_point_state_minus);
            std::vector<double> forward_propagation_history = this->ForwardHalfPhasePropagator->getPropagationHistory();
            this->ForwardHalfPhasePropagator->setStorePropagationHistory(false);
            this->ForwardHalfPhasePropagator->clearPropagationHistory();

            math::Matrix<doubleType> R_sc_Sun(3, 1, 0.0);
            if (this->myUniverse->central_body_SPICE_ID == 10)
            {
                R_sc_Sun = this->match_point_state_minus.getSubMatrix1D(0, 2);
            }
            else
            {
                //where is the central body relative to the sun?
                doubleType central_body_state_and_derivatives[12];
                this->myUniverse->locate_central_body(this->match_point_state_minus(7),
                    central_body_state_and_derivatives,
                    *this->myOptions,
                    false);

                math::Matrix<doubleType> R_CB_Sun(3, 1, 0.0);
                for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
                {
                    R_CB_Sun(stateIndex) = central_body_state_and_derivatives[stateIndex];
                }

                R_sc_Sun = this->match_point_state_minus.getSubMatrix1D(0, 2) + R_CB_Sun;
            }
            doubleType r_sc_sun_AU = R_sc_Sun.norm() / this->myOptions->AU;
            this->mySpacecraft->computePowerState(r_sc_sun_AU, this->match_point_state_minus(7));

            doubleType power = this->mySpacecraft->getAvailablePower();
            doubleType active_power = this->mySpacecraft->getEPActivePower();

            //Step 4: print the match point
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
                this->match_point_state_minus,//state
                math::Matrix<doubleType>(3, 1, 0.0),//dV
                math::Matrix<doubleType>(3, 1, 0.0),//ThrustVector
                0.0,//dVmag
                0.0,//Thrust
                0.0,//Isp
                power,//AvailPower
                0.0,//mdot
                0,//number_of_active_engines
                active_power,
                "none");//active_power

            //Step 5: print the backward half phase
            //set the propagator temporarily to dump into the output state
            this->BackwardHalfPhasePropagator->setStateRight(output_state);

            for (size_t step = 0; step < this->num_timesteps / 2; ++step)
            {
                //propagate to the halfway point of this step
                previous_epoch = step == 0 ? this->match_point_state_plus(7) : output_state(7);
                size_t backstep = this->num_timesteps / 2 - step - 1;

                this->Backward_dPropagatedStatedIndependentVariable.assign_zeros();
                this->BackwardHalfPhasePropagator->setCurrentEpoch(this->state_at_end_of_phase(7));
                this->BackwardHalfPhasePropagator->setIndexOfEpochInStateVec(7);
                this->BackwardHalfPhasePropagator->setCurrentIndependentVariable(0.0);

                if (step == 0)
                {
                    this->BackwardHalfPhasePropagator->setStorePropagationHistory(true);
                    this->BackwardHalfPhasePropagator->propagate(-output_step * (backstep + 0.5), false);
                    this->BackwardHalfPhasePropagator->setStorePropagationHistory(false);
                }
                else
                {
                    this->BackwardHalfPhasePropagator->propagate(-output_step * (backstep + 0.5), false);
                }

                math::Matrix<doubleType> R_sc_Sun(3, 1, 0.0);
                if (this->myUniverse->central_body_SPICE_ID == 10)
                {
                    R_sc_Sun = this->output_state.getSubMatrix1D(0, 2);
                }
                else
                {
                    //where is the central body relative to the sun?
                    doubleType central_body_state_and_derivatives[12];
                    this->myUniverse->locate_central_body(this->output_state(7),
                        central_body_state_and_derivatives,
                        *this->myOptions,
                        false);

                    math::Matrix<doubleType> R_CB_Sun(3, 1, 0.0);
                    for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
                    {
                        R_CB_Sun(stateIndex) = central_body_state_and_derivatives[stateIndex];
                    }

                    R_sc_Sun = this->output_state.getSubMatrix1D(0, 2) + R_CB_Sun;
                }
                doubleType r_sc_sun_AU = R_sc_Sun.norm() / this->myOptions->AU;
                this->mySpacecraft->computePowerState(r_sc_sun_AU, this->output_state(7));

                doubleType power = this->mySpacecraft->getAvailablePower();
                doubleType active_power = this->mySpacecraft->getEPActivePower();

                //print
                this->write_output_line(outputfile,//outputfile
                    eventcount,//eventcount
                    "coast",//event_type
                    "deep-space",//event_location
                    (output_state(7) - previous_epoch) / 86400.0,// timestep_size,
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
                    power,//AvailPower
                    0.0,//mdot
                    0,//number_of_active_engines
                    active_power,
                    "none");//active_power
            } // end backward output

            //reset the propagator to go where it is supposed to
            this->BackwardHalfPhasePropagator->setStateRight(this->match_point_state_plus);
            std::vector<double> backward_propagation_history = this->BackwardHalfPhasePropagator->getPropagationHistory();
            this->BackwardHalfPhasePropagator->setStorePropagationHistory(false);
            this->BackwardHalfPhasePropagator->clearPropagationHistory();

            //Step 6: print the arrival event
            this->myArrivalEvent->output(outputfile, this->LaunchDate _GETVALUE, eventcount);

            // print out Sundman step size information
            outputfile << std::endl;
            outputfile << "s-domain integration step length: " << this->integration_step_length << std::endl;
            outputfile << "Left boundary integration step length (s): " << forward_propagation_history[0] << std::endl;
            outputfile << "Right boundary integration step length (s): " << std::fabs(backward_propagation_history[0]) << std::endl;
        }

        void SundmanCoastPhase::output_ephemeris(std::ofstream& outputfile, std::ofstream& acceleration_model_file)
        {
            // Reset the spacecraft
            this->mySpacecraft->setActiveStage(this->stageIndex);

            //Step 1: output the departure event
            this->myDepartureEvent->output_ephemeris(outputfile);

            if (this->myOptions->generate_acceleration_model_instrumentation_file && !this->isKeplerian)
            {
                this->temp_state = this->myDepartureEvent->get_state_before_event();
                this->mySpacecraftAccelerationModel->populateInstrumentationFile(acceleration_model_file, this->temp_state, this->temp_state(7));
            }

            //Step 2: print the phase
            {
                //Step 2.1: we'll need an output vector

                //Step 2.2: temporarily assign the propagator to the output state
                this->ForwardHalfPhasePropagator->setStateRight(output_state);

                //Step 2.3: propagate and print, skipping the first entry and propagating one degree of true anomaly at a time
                double anomalyToPropagate = math::deg2rad;
                while (anomalyToPropagate < this->SundmanIndependentVariable)
                {
                    //Step 2.3.1: propagate
                    this->Forward_dPropagatedStatedIndependentVariable.assign_zeros();
                    this->ForwardHalfPhasePropagator->setCurrentEpoch(this->state_after_initial_TCM(7));
                    this->ForwardHalfPhasePropagator->setIndexOfEpochInStateVec(7);
                    this->ForwardHalfPhasePropagator->setCurrentIndependentVariable(0.0);
                    this->ForwardHalfPhasePropagator->propagate(anomalyToPropagate, false);
                    output_state(6) = this->state_after_initial_TCM(6);
                    this->temp_state = output_state;

                    //Step 2.3.2: convert to Sun-centered if necessary
                    if (this->myUniverse->central_body.spice_ID != this->myOptions->forward_integrated_ephemeris_central_body_SPICE_ID)
                    {
                        double LT_dump;
                        double bodyStateDouble[6];
                        spkez_c(this->myUniverse->central_body.spice_ID, output_state(7)_GETVALUE - (51544.5 * 86400.0), "J2000", "NONE", this->myOptions->forward_integrated_ephemeris_central_body_SPICE_ID, bodyStateDouble, &LT_dump);

                        for (size_t stateIndex = 0; stateIndex < 6; ++stateIndex)
                            output_state(stateIndex) += bodyStateDouble[stateIndex];
                    }

                    //Step 2.3.3: print
                    if (this->myOptions->generate_acceleration_model_instrumentation_file && !this->isKeplerian)
                    {
                        this->mySpacecraftAccelerationModel->populateInstrumentationFile(acceleration_model_file, this->temp_state, this->temp_state(7));
                    }

                    this->write_ephemeris_line(outputfile,
                        output_state,
                        math::Matrix<doubleType>(3, 1, 0.0),//control vector
                        0.0,
                        0.0,
                        0.0,
                        0,
                        0.0,
                        "none");

                    //Step 2.3.4: increment propagatedEpoch
                    anomalyToPropagate += (this->SundmanIndependentVariable _GETVALUE - anomalyToPropagate) > math::deg2rad
                        ? math::deg2rad
                        : (this->SundmanIndependentVariable _GETVALUE - anomalyToPropagate);
                }

                //Step 2.4: reset the propagator to its original state vectors
                this->ForwardHalfPhasePropagator->setStateRight(this->match_point_state_minus);
            }

            //Step 3: print the arrival event
            this->myArrivalEvent->output_ephemeris(outputfile);
            this->temp_state = this->myArrivalEvent->get_state_before_event();
            if (this->myOptions->generate_acceleration_model_instrumentation_file && !this->isKeplerian)
            {
                this->mySpacecraftAccelerationModel->populateInstrumentationFile(acceleration_model_file, this->temp_state, this->temp_state(7));
            }

        }//end output_ephemeris()

        void SundmanCoastPhase::output_maneuver_and_target_spec(std::ofstream& maneuver_spec_file, std::ofstream& target_spec_file, bool& haveManeuverNeedTarget)
        {
            // Reset the spacecraft
            this->mySpacecraft->setActiveStage(this->stageIndex);

            //should write out a maneuver/target spec for the departure maneuver if there is one
            this->myDepartureEvent->output_maneuver_and_target_spec(maneuver_spec_file, target_spec_file, haveManeuverNeedTarget);

            //coast phases don't write maneuver and target spec files anyway, so nothing for the phase body

            //should write out a maneuver/target spec for the arrival maneuver if there is one
            this->myArrivalEvent->output_maneuver_and_target_spec(maneuver_spec_file, target_spec_file, haveManeuverNeedTarget);
        }//end output_maneuver_and_target_spec()

    }//end namespace Phases
}//end namespace EMTG