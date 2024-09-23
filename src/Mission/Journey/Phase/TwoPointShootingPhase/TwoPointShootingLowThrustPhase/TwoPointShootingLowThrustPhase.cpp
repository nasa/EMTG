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

//EMTGv9 TwoPointShootingPhase
//Jacob Englander 6-22-2017

#include "TwoPointShootingLowThrustPhase.h"
#include "PhaseDistanceConstraintFactory.h"

namespace EMTG
{
    namespace Phases
    {
        TwoPointShootingLowThrustPhase::TwoPointShootingLowThrustPhase(const std::string& name,
            const size_t& journeyIndex,
            const size_t& phaseIndex,
            size_t& stageIndex,
            phase* previousPhase,
            Astrodynamics::universe* Universe,
            HardwareModels::Spacecraft* mySpacecraft,
            HardwareModels::LaunchVehicle* myLaunchVehicle,
            missionoptions* myOptions,
            const size_t& numStatesToPropagate,
            const size_t& numMatchConstraints) :
            TwoPointShootingPhase::TwoPointShootingPhase(name,
                journeyIndex,
                phaseIndex,
                stageIndex,
                previousPhase,
                Universe,
                mySpacecraft,
                myLaunchVehicle,
                myOptions,
                numStatesToPropagate,
                numMatchConstraints)
        {
            //low-thrust phases have an electric maneuver
            this->hasElectricManeuver = true;

            //time things
            this->num_timesteps = this->myJourneyOptions->override_num_steps
                ? this->myJourneyOptions->number_of_steps
                : this->myOptions->num_timesteps;

            //if we have an odd number of timesteps, TPSLTPs will crash, so make it an even number
            this->num_timesteps = this->num_timesteps % 2 == 0 ? this->num_timesteps : this->num_timesteps + 1;

            this->ThrustStepLengths.resize(this->num_timesteps);
            this->dThrustStepLengths_dPropagationVariable.resize(this->num_timesteps);

            //state
            this->spacecraft_state_event_minus.resize(this->num_timesteps, math::Matrix<doubleType>(this->numStatesToPropagate, 1, 0.0));
            this->spacecraft_state_event_plus.resize(this->num_timesteps, math::Matrix<doubleType>(this->numStatesToPropagate, 1, 0.0));

            //control
            SpacecraftThrusterMode myThrusterMode = this->mySpacecraft->getCurrentStageOptions().getElectricPropulsionSystemOptions().getThrusterMode();
            if (myThrusterMode == FixedEfficiencyVSI
                || myThrusterMode == Stepped2D
                || myThrusterMode == Poly2D)
            {
                this->isVSI = true;
                this->num_controls = 4;
            }
            else
            {
                this->isVSI = false;
                this->num_controls = 3;
            }


            this->ControlVector.resize(this->num_timesteps, math::Matrix<doubleType>(this->num_controls, 1, 0.0));
            this->throttle.resize(this->num_timesteps);
            TruthTable_MatchConstraints_Derivative_wrt_Control.resize(this->numMatchConstraints, true);

            //hardware things
            this->max_thrust.resize(this->num_timesteps);
            this->max_mass_flow_rate.resize(this->num_timesteps);
            this->Isp.resize(this->num_timesteps);
            this->power.resize(this->num_timesteps);
            this->active_power.resize(this->num_timesteps);
            this->number_of_active_engines.resize(this->num_timesteps, 0);
            this->ThrottleLevel.resize(this->num_timesteps, -1);
            this->ThrottleLevelString.resize(this->num_timesteps, "none");

            //delta-v
            this->stepDeltavMagnitude.resize(this->num_timesteps);
            this->stepDeltaV.resize(this->num_timesteps, math::Matrix<doubleType>(3, 1, 0.0));
            this->G_indices_control_magnitude_constraints_wrt_Control.resize(this->num_timesteps);

            //distance contraint
            std::string shortprefix = "p" + std::to_string(this->phaseIndex);

            for (std::string& constraint : this->myJourneyOptions->PhaseDistanceConstraintDefinitions)
            {
                if (constraint.find("#") != 0
                    && (constraint.find(shortprefix) < 1024 || (constraint.find("pEnd") < 1024 && this->isLastPhaseInJourney)))
                {
                    this->distance_constraint_definitions.push_back(CreatePhaseDistanceConstraint(constraint, this->myOptions, this->myUniverse));
                }
            }
            this->number_of_distance_constraints = this->distance_constraint_definitions.size();

            if (this->number_of_distance_constraints > 0)
            {

                this->distance_constraint_relative_position.resize(this->num_timesteps, std::vector< math::Matrix<doubleType> >(this->number_of_distance_constraints, math::Matrix<doubleType>(3, 1, 0.0))); //step, constraint, state variable
                this->distance_constraint_body_position_time_derivatives.resize(this->num_timesteps, std::vector< math::Matrix<double> >(this->number_of_distance_constraints, math::Matrix<double>(3, 1, 0.0))); //step, constraint, state 
                this->distance_from_body.resize(this->num_timesteps, std::vector<doubleType>(this->number_of_distance_constraints, 0.0)); //step, constraint

                //derivative indices
                this->G_indices_distance_constraints_wrt_ForwardControl.resize(this->num_timesteps / 2, std::vector< std::vector< std::vector<size_t> > >(this->number_of_distance_constraints)); //step, constraint, controlstep, control variable
                this->G_indices_distance_constraints_wrt_BackwardControl.resize(this->num_timesteps / 2, std::vector< std::vector< std::vector<size_t> > >(this->number_of_distance_constraints)); //step, constraint, controlstep, control variable
                this->G_indices_distance_constraints_wrt_LeftBoundaryState.resize(this->num_timesteps, std::vector< std::vector<size_t> >(this->number_of_distance_constraints)); //step, constraint, state variable
                this->G_indices_distance_constraints_wrt_RightBoundaryState.resize(this->num_timesteps, std::vector< std::vector<size_t> >(this->number_of_distance_constraints)); //step, constraint, state variable
                this->G_indices_distance_constraints_wrt_LeftBoundaryTime.resize(this->num_timesteps, std::vector< std::vector<size_t> >(this->number_of_distance_constraints)); //step, constraint, state variable
                this->G_indices_distance_constraints_wrt_RightBoundaryTime.resize(this->num_timesteps, std::vector< std::vector<size_t> >(this->number_of_distance_constraints)); //step, constraint, state variable
                this->G_indices_distance_constraints_wrt_PhaseFlightTime.resize(this->num_timesteps, std::vector<size_t>(this->number_of_distance_constraints));
            }
        }//end constructor

        //******************************************calcbounds methods
        void TwoPointShootingLowThrustPhase::calcbounds_phase_main()
        {
            this->calcbounds_virtual_propellant_tanks();

            this->calcbounds_control();

            this->calcbounds_match_point_constraints();

            if (this->number_of_distance_constraints > 0)
                this->calcbounds_distance_constraints();


            if (this->myOptions->objective_type == 0)
                this->calcbounds_deltav_contribution();
        }//end calcbounds_phase_main()

        void TwoPointShootingLowThrustPhase::calcbounds_control()
        {
            std::vector<std::string> controlNames({ "x", "y", "z", "command" });
            for (size_t step = 0; step < this->num_timesteps; ++step)
            {
                std::vector<size_t> step_control_Xindex;
                //control unit vector
                for (size_t controlIndex = 0; controlIndex < 3; ++controlIndex)
                {
                    this->Xlowerbounds->push_back(-1.0);
                    this->Xupperbounds->push_back(1.0);
                    this->X_scale_factors->push_back(1.0);
                    this->Xdescriptions->push_back(prefix + "step " + std::to_string(step) + " u_" + controlNames[controlIndex]);
                    step_control_Xindex.push_back(this->Xdescriptions->size() - 1);
                }

                //u_command
                if (this->isVSI)
                {
                    this->Xlowerbounds->push_back(0.0);
                    this->Xupperbounds->push_back(1.0);
                    this->X_scale_factors->push_back(1.0);
                    this->Xdescriptions->push_back(prefix + "step " + std::to_string(step) + " u_" + controlNames[3]);
                    step_control_Xindex.push_back(this->Xdescriptions->size() - 1);
                }
                this->Xindex_control.push_back(step_control_Xindex);

                //control magnitude constraint
                if (this->myJourneyOptions->force_unit_magnitude_control == ControlMagnitudeType::UnitMagnitude)
                {
                    this->Flowerbounds->push_back(1.0);
                }
                else
                {
                    this->Flowerbounds->push_back(0.0);
                }
                if (this->myJourneyOptions->force_unit_magnitude_control == ControlMagnitudeType::ZeroMagnitude)
                {
                    this->Fupperbounds->push_back(0.0);
                }
                else
                {
                    this->Fupperbounds->push_back(1.0);
                }
                this->Fdescriptions->push_back(prefix + "step " + std::to_string(step) + " control magnitude");

                //sparsity pattern
                this->create_sparsity_entry_vector(Fdescriptions->size() - 1,
                    this->Xdescriptions->size() - 1 - this->num_controls,
                    true,
                    3,
                    "step " + std::to_string(step) + " u",
                    this->G_indices_control_magnitude_constraints_wrt_Control[step]);

                //if requested by the user, force the control in all steps to be the same as the first step
                if (this->myJourneyOptions->force_fixed_inertial_control && step > 0)
                {
                    std::vector<size_t> step_Gindices_fixed_inertial_control_wrt_this_step_control;
                    std::vector<size_t> step_Gindices_fixed_inertial_control_wrt_first_step_control;

                    for (size_t controlIndex = 0; controlIndex < this->Xindex_control[step].size(); ++controlIndex)
                    {
                        //u_i - u_0 = 0.0
                        //implement as nonlinear constraint so that SNOPT can temporarily violate it during the solve
                        this->Flowerbounds->push_back(-math::SMALL);
                        this->Fupperbounds->push_back(math::SMALL);
                        this->Fdescriptions->push_back(this->prefix + "Step " + std::to_string(step) + " u_" + controlNames[controlIndex] + " must match step 0");
                        this->create_sparsity_entry(this->Fdescriptions->size() - 1, this->Xindex_control[step][controlIndex], step_Gindices_fixed_inertial_control_wrt_this_step_control);
                        this->create_sparsity_entry(this->Fdescriptions->size() - 1, this->Xindex_control[0][controlIndex], step_Gindices_fixed_inertial_control_wrt_first_step_control);
                    }

                    this->Gindices_fixed_inertial_control_wrt_this_step_control.push_back(step_Gindices_fixed_inertial_control_wrt_this_step_control);
                    this->Gindices_fixed_inertial_control_wrt_first_step_control.push_back(step_Gindices_fixed_inertial_control_wrt_first_step_control);
                }//end fixed control vector constraint
            }//end loop over steps
        }//end calcbounds_control_variables()
        
        void TwoPointShootingLowThrustPhase::calcbounds_virtual_propellant_tanks()
        {

            //Step 1: virtual chemical fuel
            {
                this->Xlowerbounds->push_back(0.0);
                this->Xupperbounds->push_back(this->myJourneyOptions->maximum_mass);
                this->X_scale_factors->push_back(this->myJourneyOptions->maximum_mass);
                this->Xdescriptions->push_back(prefix + "virtual chemical fuel");
                size_t Xindex = this->Xdescriptions->size() - 1;

                //tell the spacecraft about this variable
                //global constraint
                this->mySpacecraft->appendGlobalChemicalFuelTank_Xindices(Xindex);
                this->mySpacecraft->appendGlobalChemicalFuelTank_Xscaleranges(this->Xupperbounds->back(), this->Xlowerbounds->back());
                //stage constraint
                this->mySpacecraft->appendChemicalFuelTank_Xindices(this->stageIndex, Xindex);
                this->mySpacecraft->appendChemicalFuelTank_Xscaleranges(this->stageIndex, this->Xupperbounds->back(), this->Xlowerbounds->back());
            }//end virtual chemical fuel tank

            //Step 2: virtual electric propellant
            {   
                this->Xlowerbounds->push_back(0.0);
                this->Xupperbounds->push_back(this->myJourneyOptions->maximum_mass);
                this->X_scale_factors->push_back(this->myJourneyOptions->maximum_mass);
                this->Xdescriptions->push_back(prefix + "virtual electric propellant");
                size_t Xindex = this->Xdescriptions->size() - 1;

                //tell the spacecraft about this variable
                //global constraint
                this->mySpacecraft->appendGlobalElectricPropellantTank_Xindices(Xindex);
                this->mySpacecraft->appendGlobalElectricPropellantTank_Xscaleranges(this->Xupperbounds->back(), this->Xlowerbounds->back());
                //stage constraint
                this->mySpacecraft->appendElectricPropellantTank_Xindices(this->stageIndex, Xindex);
                this->mySpacecraft->appendElectricPropellantTank_Xscaleranges(this->stageIndex, this->Xupperbounds->back(), this->Xlowerbounds->back());
            }//end virtual electric propellant tank
        }//end calcbounds_virtual_propellant_tanks

        void TwoPointShootingLowThrustPhase::calcbounds_deltav_contribution()
        {
            //the phase delta-v has derivatives with respect to boundary, phase flight time, and control

            //with respect to left boundary
            for (size_t varIndex = 0; varIndex < this->ListOfVariablesAffectingLeftBoundary.size(); ++varIndex)
            {
                size_t Xindex = this->ListOfVariablesAffectingLeftBoundary[varIndex];

                this->create_sparsity_entry(0,
                    Xindex,
                    this->G_indices_deltav_wrt_LeftBoundaryState);
            }

            //with respect to right boundary
            for (size_t varIndex = 0; varIndex < this->ListOfVariablesAffectingRightBoundary.size(); ++varIndex)
            {
                size_t Xindex = this->ListOfVariablesAffectingRightBoundary[varIndex];

                this->create_sparsity_entry(0,
                    Xindex,
                    this->G_indices_deltav_wrt_RightBoundaryState);
            }

            //derivatives with respect to time variables that affect boundary states
            //with respect to left boundary
            for (size_t varIndex = 0; varIndex < this->ListOfTimeVariablesAffectingLeftBoundary.size(); ++varIndex)
            {
                size_t Xindex = this->ListOfTimeVariablesAffectingLeftBoundary[varIndex];

                this->create_sparsity_entry(0,
                    Xindex,
                    this->G_indices_deltav_wrt_LeftBoundaryState_wrt_Time);
            }

            //with respect to right boundary
            for (size_t varIndex = 0; varIndex < this->ListOfTimeVariablesAffectingRightBoundary.size(); ++varIndex)
            {
                size_t Xindex = this->ListOfTimeVariablesAffectingRightBoundary[varIndex];

                this->create_sparsity_entry(0,
                    Xindex,
                    this->G_indices_deltav_wrt_RightBoundaryState_wrt_Time);
            }

            //derivative with respect to phase flight time
            this->create_sparsity_entry(0,
                this->Xindex_PhaseFlightTime,
                this->G_indices_deltav_wrt_PhaseFlightTime);

            //derivatives with respect to control
            for (size_t step = 0; step < this->num_timesteps; ++step)
            {
                std::vector<size_t> temp_indices;

                this->create_sparsity_entry_vector(0,
                    this->X_index_of_first_decision_variable_in_this_phase,
                    true,
                    this->num_controls,
                    "step " + std::to_string(step) + " u",
                    temp_indices);

                this->G_indices_deltav_wrt_Control.push_back(temp_indices);
            }
        }//end calcbounds_deltav_contribution

        void TwoPointShootingLowThrustPhase::calcbounds_match_point_constraints()
        {
            //base class
            TwoPointShootingPhase::calcbounds_match_point_constraints();

            //sparsity pattern with respect to control
            for (size_t constraintIndex = 0; constraintIndex < this->numMatchConstraints; ++constraintIndex)
            {
                std::vector< std::vector<size_t> > state_G_indices_match_point_constraints_wrt_Control;

                for (size_t step = 0; step < this->num_timesteps; ++step)
                {
                    std::vector<size_t> step_G_indices_match_point_constraints_wrt_Control;

                    if (this->TruthTable_MatchConstraints_Derivative_wrt_Control[constraintIndex])
                    {
                        this->create_sparsity_entry_vector(this->Findices_match_point_constraints[constraintIndex],
                            this->X_index_of_first_decision_variable_in_this_phase,
                            true,
                            this->num_controls,
                            "step " + std::to_string(step) + " u",
                            step_G_indices_match_point_constraints_wrt_Control);
                    }
                    else
                        step_G_indices_match_point_constraints_wrt_Control = std::vector<size_t>(this->num_controls, 0);

                    state_G_indices_match_point_constraints_wrt_Control.push_back(step_G_indices_match_point_constraints_wrt_Control);
                }

                this->G_indices_match_point_constraints_wrt_Control.push_back(state_G_indices_match_point_constraints_wrt_Control);
            }


            //virtual chemical fuel
            this->create_sparsity_entry(this->Findices_match_point_constraints[7],
                this->X_index_of_first_decision_variable_in_this_phase,
                true,
                this->prefix + "virtual chemical fuel",
                this->Gindex_dVirtualChemicalFuelConstraint_dVirtualChemicalFuel);

            //virtual electric propellant
            this->create_sparsity_entry(this->Findices_match_point_constraints[8],
                this->X_index_of_first_decision_variable_in_this_phase,
                true,
                this->prefix + "virtual electric propellant",
                this->Gindex_dVirtualElectricPropellantConstraint_dVirtualElectricPropellant);
        }//end calcbounds_match_point_constraints()

        void TwoPointShootingLowThrustPhase::calcbounds_distance_constraints()
        {
            //create constraints
            for (size_t step = 0; step < this->num_timesteps; ++step)
            {
                for (size_t constraintIndex = 0; constraintIndex < this->number_of_distance_constraints; ++constraintIndex)
                {
                    //create the constraint
                    this->Flowerbounds->push_back((std::get<1>(this->distance_constraint_definitions[constraintIndex]) - std::get<2>(this->distance_constraint_definitions[constraintIndex])) / this->myUniverse->LU);
                    this->Fupperbounds->push_back(0.0);
                    int bodyIndex = std::get<0>(this->distance_constraint_definitions[constraintIndex]);
                    if (bodyIndex == -2)
                    {
                        this->Fdescriptions->push_back(prefix + " step " + std::to_string(step) + " distance from " + this->myUniverse->central_body.name);
                    }
                    else
                    {
                        this->Fdescriptions->push_back(prefix + " step " + std::to_string(step) + " distance from " + this->myUniverse->bodies[bodyIndex - 1].name);
                    }

                    //derivatives with respect to boundary states
                    for (size_t encodedStateIndex = 0; encodedStateIndex < 8; ++encodedStateIndex)
                    {
                        if (step < this->num_timesteps / 2) //left half-phase
                        {
                            for (size_t varIndex = 0; varIndex < this->ListOfVariablesAffectingLeftBoundary.size(); ++varIndex)
                            {
                                size_t Xindex = this->ListOfVariablesAffectingLeftBoundary[varIndex];

                                this->create_sparsity_entry(Fdescriptions->size() - 1,
                                    Xindex,
                                    this->G_indices_distance_constraints_wrt_LeftBoundaryState[step][constraintIndex]);
                            }

                            for (size_t varIndex = 0; varIndex < this->ListOfTimeVariablesAffectingLeftBoundary.size(); ++varIndex)
                            {
                                size_t Xindex = this->ListOfTimeVariablesAffectingLeftBoundary[varIndex];

                                this->create_sparsity_entry(Fdescriptions->size() - 1,
                                    Xindex,
                                    this->G_indices_distance_constraints_wrt_LeftBoundaryTime[step][constraintIndex]);
                            }
                        }
                        else //right half-phase
                        {
                            for (size_t varIndex = 0; varIndex < this->ListOfVariablesAffectingRightBoundary.size(); ++varIndex)
                            {
                                size_t Xindex = this->ListOfVariablesAffectingRightBoundary[varIndex];

                                this->create_sparsity_entry(Fdescriptions->size() - 1,
                                    Xindex,
                                    this->G_indices_distance_constraints_wrt_RightBoundaryState[step][constraintIndex]);
                            }

                            for (size_t varIndex = 0; varIndex < this->ListOfTimeVariablesAffectingRightBoundary.size(); ++varIndex)
                            {
                                size_t Xindex = this->ListOfTimeVariablesAffectingRightBoundary[varIndex];

                                this->create_sparsity_entry(Fdescriptions->size() - 1,
                                    Xindex,
                                    this->G_indices_distance_constraints_wrt_RightBoundaryTime[step][constraintIndex]);
                            }
                        }
                    }

                    //derivatives with respect to phase flight time
                    this->create_sparsity_entry(Fdescriptions->size() - 1,
                        this->Xindex_PhaseFlightTime,
                        this->G_indices_distance_constraints_wrt_PhaseFlightTime[step][constraintIndex]);

                    //derivatives with respect to control - this is a triangle!
                    if (step < this->num_timesteps / 2)
                    {
                        //derivative with respect to all previous control variables
                        for (size_t controlstep = 0; controlstep < step; ++controlstep)
                        {
                            std::vector<size_t> controlstep_G_indices_distance_constraints_wrt_Control;

                            this->create_sparsity_entry_vector(Fdescriptions->size() - 1,
                                this->X_index_of_first_decision_variable_in_this_phase,
                                true,
                                this->num_controls,
                                "step " + std::to_string(controlstep) + " u",
                                controlstep_G_indices_distance_constraints_wrt_Control);

                            this->G_indices_distance_constraints_wrt_ForwardControl[step][constraintIndex].push_back(controlstep_G_indices_distance_constraints_wrt_Control);
                        }
                    }
                    else
                    {
                        //derivative with respect to all previous control variables
                        for (size_t controlstep = this->num_timesteps - 1; controlstep > step; --controlstep)
                        {
                            std::vector<size_t> controlstep_G_indices_distance_constraints_wrt_Control;

                            this->create_sparsity_entry_vector(Fdescriptions->size() - 1,
                                this->X_index_of_first_decision_variable_in_this_phase,
                                true,
                                this->num_controls,
                                "step " + std::to_string(controlstep) + " u",
                                controlstep_G_indices_distance_constraints_wrt_Control);

                            this->G_indices_distance_constraints_wrt_BackwardControl[this->num_timesteps - step - 1][constraintIndex].push_back(controlstep_G_indices_distance_constraints_wrt_Control);
                        }
                    }
                }
            }
        }//end calcbounds_distance_constraints()

        //******************************************process methods
        void TwoPointShootingLowThrustPhase::
            process_control(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //extract all of the control variables - yay!
            for (size_t step = 0; step < this->num_timesteps; ++step)
            {
                for (size_t controlIndex = 0; controlIndex < this->num_controls; ++controlIndex)
                {
                    this->ControlVector[step](controlIndex) = X[Xindex++];

                }//end loop over controls

                //control magnitude constraint
                this->throttle[step] = sqrt(  this->ControlVector[step](0) * this->ControlVector[step](0)
                                            + this->ControlVector[step](1) * this->ControlVector[step](1)
                                            + this->ControlVector[step](2) * this->ControlVector[step](2)
                                            + 1.0e-25);
                F[Findex++] = this->throttle[step];


                //if requested by the user, force the control in all steps to be the same as the first step
                if (this->myJourneyOptions->force_fixed_inertial_control && step > 0)
                {
                    for (size_t controlIndex = 0; controlIndex < this->num_controls; ++controlIndex)
                    {
                        //u_i - u_0 = 0.0
                        //implement as nonlinear constraint so that SNOPT can temporarily violate it during the solve
                        F[Findex++] = X[this->Xindex_control[step][controlIndex]] - X[this->Xindex_control[0][controlIndex]];

                        G[this->Gindices_fixed_inertial_control_wrt_this_step_control[step - 1][controlIndex]] = 1.0;
                        G[this->Gindices_fixed_inertial_control_wrt_first_step_control[step - 1][controlIndex]] = -1.0;
                    }
                }//end fixed control vector constraint
            }//end loop over steps

            //derivatives
            if (needG)
            {
                for (size_t step = 0; step < this->num_timesteps; ++step)
                {
                    for (size_t controlIndex = 0; controlIndex < this->num_controls; ++controlIndex)
                    {
                        size_t Gindex = this->G_indices_control_magnitude_constraints_wrt_Control[step][controlIndex];
                        size_t Xindex = this->jGvar->operator[](Gindex);

                        G[Gindex] = this->X_scale_factors->operator[](Xindex)
                            * (this->ControlVector[step](controlIndex) / this->throttle[step]) _GETVALUE;
                    }//end loop over controls
                }//end loop over steps
            }//end derivatives
        }//end process_control()

        void TwoPointShootingLowThrustPhase::
            process_virtual_propellant_tanks(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            this->virtual_chemical_fuel_used = X[Xindex++];
            this->virtual_electric_propellant_used = X[Xindex++];
        }//end process_virtual_propellant_tanks()



        void TwoPointShootingLowThrustPhase::output_STMs()
        {
            // Reset the spacecraft
            this->mySpacecraft->setActiveStage(this->stageIndex);

            std::ofstream summaryfile(this->myOptions->working_directory + "//" + boost::replace_all_copy(this->prefix, ": ", "") + ".stm_description");
            for (size_t stateIndex = 0; stateIndex < this->stateVectorNames.size(); ++stateIndex)
                summaryfile << this->stateVectorNames[stateIndex] << std::endl;
            summaryfile.close();

            size_t STMindex = 0;

            this->STM_Augmented_initial_coast.print_to_file(this->myOptions->working_directory + "//" + boost::replace_all_copy(this->prefix, ": ", "") + "_STM_" + std::to_string(STMindex++) + "_forward.stm", true);

            for (size_t stepIndex = 0; stepIndex < this->num_timesteps / 2; ++stepIndex)
                this->ForwardAugmentedSTM[stepIndex].print_to_file(this->myOptions->working_directory + "//" + boost::replace_all_copy(this->prefix, ": ", "") + "_STM_" + std::to_string(STMindex++) + "_forward.stm", true);

            for (int stepIndex = this->num_timesteps / 2 - 1; stepIndex >= 0; --stepIndex)
                this->BackwardAugmentedSTM[stepIndex].print_to_file(this->myOptions->working_directory + "//" + boost::replace_all_copy(this->prefix, ": ", "") + "_STM_" + std::to_string(STMindex++) + "_backward.stm", true);
            
            this->STM_Augmented_terminal_coast.print_to_file(this->myOptions->working_directory + "//" + boost::replace_all_copy(this->prefix, ": ", "") + "_STM_" + std::to_string(STMindex++) + "_backward.stm", true);
        }//end output_STMs()
    }//close namespace Phases
}//close namespace EMTG