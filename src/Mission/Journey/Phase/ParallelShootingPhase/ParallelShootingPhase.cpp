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

//EMTGv9 ParallelShootingPhase base class
//Jacob Englander 2-21-2018

#include "ParallelShootingPhase.h"

namespace EMTG
{
    namespace Phases
    {
        ParallelShootingPhase::ParallelShootingPhase()
        {}

        ParallelShootingPhase::ParallelShootingPhase(const std::string& name,
            const size_t& journeyIndex,
            const size_t& phaseIndex,
            size_t& stageIndex,
            phase* previousPhase,
            Astrodynamics::universe* Universe,
            HardwareModels::Spacecraft* mySpacecraft,
            HardwareModels::LaunchVehicle* myLaunchVehicle,
            missionoptions* myOptions) :
            phase::phase(name,
                journeyIndex,
                phaseIndex,
                stageIndex,
                previousPhase,
                Universe,
                mySpacecraft,
                myLaunchVehicle,
                myOptions,
                10, //numStatesToPropagate
                9) //numMatchConstraints
        {
            this->initialize(name,
                journeyIndex,
                phaseIndex,
                stageIndex,
                previousPhase,
                Universe,
                mySpacecraft,
                myLaunchVehicle,
                myOptions);
        }//end constructor


        void ParallelShootingPhase::initialize(const std::string& name,
            const size_t& journeyIndex,
            const size_t& phaseIndex,
            size_t& stageIndex,
            phase* previousPhase,
            Astrodynamics::universe* Universe,
            HardwareModels::Spacecraft* mySpacecraft,
            HardwareModels::LaunchVehicle* myLaunchVehicle,
            missionoptions* myOptions)
        {
            //name our match point constraints
            for (size_t stateIndex = 0; stateIndex < 7; ++stateIndex)
            {
                this->matchPointConstraintNames.push_back(this->stateVectorNames[stateIndex]);
                this->matchPointConstraintStateIndex.push_back(stateIndex);
            }
            this->stateVectorNames.push_back("virtual chemical fuel");
            this->stateVectorNames.push_back("virtual electric propellant");
            this->matchPointConstraintNames.push_back("virtual chemical fuel");
            this->matchPointConstraintNames.push_back("virtual electric propellant");
            this->matchPointConstraintStateIndex.push_back(8); //note we skipped the encode of boundary epoch, it does not get a match point constraint
            this->matchPointConstraintStateIndex.push_back(9);
            this->continuity_constraint_scale_factors(7) = 1.0 / this->myJourneyOptions->maximum_mass;
            this->continuity_constraint_scale_factors(8) = 1.0 / this->myJourneyOptions->maximum_mass;
            this->stateIndex_phase_propagation_variable = 13;
            this->stateVectorNames.push_back("u_x");
            this->stateVectorNames.push_back("u_y");
            this->stateVectorNames.push_back("u_z");
            this->stateVectorNames.push_back("phase flight time");


            //own members
            this->num_steps = this->myJourneyOptions->override_num_steps
                ? this->myJourneyOptions->number_of_steps
                : this->myOptions->num_timesteps;

            this->state_after_initial_coast.resize(10, 1, 0.0);
            this->state_before_terminal_coast.resize(10, 1, 0.0);

            this->ForcedCoast_dStepSize_dPropagationVariable = 0.0;
            this->dStepTime_dPhaseFlightTime = 1.0 / this->num_steps;

            if (this->hasInitialCoast)
            {
                this->STM_initial_coast.resize(14, 14, math::identity);
                this->STM_Augmented_initial_coast.resize(14, 14, math::identity);
                this->InitialCoast_dStatedIndependentVariable.resize(10, 2, 0.0);
            }

            if (this->hasTerminalCoast)
            {
                this->STM_terminal_coast.resize(14, 14, math::identity);
                this->STM_Augmented_terminal_coast.resize(14, 14, math::identity);
                this->TerminalCoast_dStatedIndependentVariable.resize(10, 2, 0.0);
            }

            this->hasElectricManeuver = true;

            //output stuff
            this->output_state.resize(this->numStatesToPropagate + 13 * 13, 1, 0.0);
        }//end initialize

        void ParallelShootingPhase::setup_calcbounds(std::vector<double>* Xupperbounds,
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
            //self
            this->phase::setup_calcbounds(Xupperbounds,
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

            //steps
            for (size_t stepIndex = 0; stepIndex < this->num_steps; ++stepIndex)
                this->mySteps[stepIndex].setup_calcbounds(Xupperbounds,
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
        }//end setup_calcbounds

        ParallelShootingPhase::~ParallelShootingPhase()
        {
        }//end destructor

        void ParallelShootingPhase::calcbounds()
        {
            this->calcbounds_phase_left_boundary();

            this->calcbounds_phase_flight_time();

            this->calcbounds_phase_right_boundary();

            this->calcbounds_virtual_propellant_tanks();

            this->calcbounds_phase_main();

            //do we want fixed inertial control?
            if (this->myJourneyOptions->force_fixed_inertial_control)
                this->calcbounds_fixed_inertial_control();
        }//end calcbounds
        
        void ParallelShootingPhase::calcbounds_virtual_propellant_tanks()
        {

            //Step 1: virtual chemical fuel
            {
                this->Xlowerbounds->push_back(0.0);
                this->Xupperbounds->push_back(this->myJourneyOptions->maximum_mass);
                this->X_scale_factors->push_back(this->myJourneyOptions->maximum_mass);
                this->Xdescriptions->push_back(prefix + "virtual chemical fuel");
                this->Xindex_virtual_chemical_fuel_tank = this->Xdescriptions->size() - 1;

                //tell the spacecraft about this variable
                //global constraint
                this->mySpacecraft->appendGlobalChemicalFuelTank_Xindices(this->Xindex_virtual_chemical_fuel_tank);
                this->mySpacecraft->appendGlobalChemicalFuelTank_Xscaleranges(this->Xupperbounds->back(), this->Xlowerbounds->back());
                //stage constraint
                this->mySpacecraft->appendChemicalFuelTank_Xindices(this->stageIndex, this->Xindex_virtual_chemical_fuel_tank);
                this->mySpacecraft->appendChemicalFuelTank_Xscaleranges(this->stageIndex, this->Xupperbounds->back(), this->Xlowerbounds->back());
            }//end virtual chemical fuel tank

             //Step 2: virtual electric propellant
            {
                this->Xlowerbounds->push_back(0.0);
                this->Xupperbounds->push_back(this->myJourneyOptions->maximum_mass);
                this->X_scale_factors->push_back(this->myJourneyOptions->maximum_mass);
                this->Xdescriptions->push_back(prefix + "virtual electric propellant");
                this->Xindex_virtual_electric_propellant_tank = this->Xdescriptions->size() - 1;

                //tell the spacecraft about this variable
                //global constraint
                this->mySpacecraft->appendGlobalElectricPropellantTank_Xindices(this->Xindex_virtual_electric_propellant_tank);
                this->mySpacecraft->appendGlobalElectricPropellantTank_Xscaleranges(this->Xupperbounds->back(), this->Xlowerbounds->back());
                //stage constraint
                this->mySpacecraft->appendElectricPropellantTank_Xindices(this->stageIndex, this->Xindex_virtual_electric_propellant_tank);
                this->mySpacecraft->appendElectricPropellantTank_Xscaleranges(this->stageIndex, this->Xupperbounds->back(), this->Xlowerbounds->back());
            }//end virtual electric propellant tank
        }//end calcbounds_virtual_propellant_tanks

        void ParallelShootingPhase::calcbounds_phase_main()
        {

            this->calcbounds_initial_coast();

            this->calcbounds_terminal_coast();

            //steps
            for (size_t stepIndex = 0; stepIndex < this->num_steps; ++stepIndex)
                this->mySteps[stepIndex].calcbounds_step();
        }//end calcbounds_phase_main

        void ParallelShootingPhase::calcbounds_initial_coast()
        {
            //StateAfterInitialCoast has derivatives with respect to EVERYTHING that affects the left boundary condition
            std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateBeforePhase = this->myDepartureEvent->get_Derivatives_of_StateAfterEvent();//Xindex, stateIndex, derivative value
            std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateBeforePhase_wrt_Time = this->myDepartureEvent->get_Derivatives_of_StateAfterEvent_wrt_Time();//Xindex, stateIndex, derivative value
                

            //Step 1: what does the left boundary have derivatives with respect to?
            {
                //Step 1.1: state variables
                for (size_t dIndex = 0; dIndex < Derivatives_of_StateBeforePhase.size(); ++dIndex)
                {
                    bool alreadyListed = false;
                    size_t Xindex = std::get<0>(Derivatives_of_StateBeforePhase[dIndex]);
                    for (size_t listIndex : this->ListOfVariablesAffectingLeftBoundary)
                    {
                        if (listIndex == Xindex)
                        {
                            alreadyListed = true;
                            break;
                        }
                    }
                    if (!alreadyListed)
                    {
                        this->ListOfVariablesAffectingLeftBoundary.push_back(Xindex);
                    }
                }

                this->DerivativesOfLeftBoundaryByVariable.resize(ListOfVariablesAffectingLeftBoundary.size());

                for (size_t listIndex = 0; listIndex < this->ListOfVariablesAffectingLeftBoundary.size(); ++listIndex)
                {
                    for (size_t dIndex = 0; dIndex < Derivatives_of_StateBeforePhase.size(); ++dIndex)
                    {
                        size_t Xindex = std::get<0>(Derivatives_of_StateBeforePhase[dIndex]);

                        if (Xindex == this->ListOfVariablesAffectingLeftBoundary[listIndex])
                        {
                            this->DerivativesOfLeftBoundaryByVariable[listIndex].push_back(dIndex);
                        }
                    }
                }

                //Step 1.2: time variables
                for (size_t dIndex = 0; dIndex < Derivatives_of_StateBeforePhase_wrt_Time.size(); ++dIndex)
                {
                    bool alreadyListed = false;
                    size_t Xindex = std::get<0>(Derivatives_of_StateBeforePhase_wrt_Time[dIndex]);
                    for (size_t listIndex : this->ListOfTimeVariablesAffectingLeftBoundary)
                    {
                        if (listIndex == Xindex)
                        {
                            alreadyListed = true;
                            break;
                        }
                    }
                    if (!alreadyListed)
                    {
                        this->ListOfTimeVariablesAffectingLeftBoundary.push_back(Xindex);
                    }
                }

                this->DerivativesOfLeftBoundaryByTimeVariable.resize(ListOfTimeVariablesAffectingLeftBoundary.size());

                for (size_t listIndex = 0; listIndex < this->ListOfTimeVariablesAffectingLeftBoundary.size(); ++listIndex)
                {
                    for (size_t dIndex = 0; dIndex < Derivatives_of_StateBeforePhase_wrt_Time.size(); ++dIndex)
                    {
                        size_t Xindex = std::get<0>(Derivatives_of_StateBeforePhase_wrt_Time[dIndex]);

                        if (Xindex == this->ListOfTimeVariablesAffectingLeftBoundary[listIndex])
                        {
                            this->DerivativesOfLeftBoundaryByTimeVariable[listIndex].push_back(dIndex);
                        }
                    }
                }
            }

            //Step 2: make the derivative skeleton of StateAfterInitialCoast
            //Step 2.1: 7-state w.r.t. things that affect 7-state
            //state
            for (std::vector<size_t> DerivativesOfLeftBoundaryThisVariable : this->DerivativesOfLeftBoundaryByVariable)
            {
                std::vector<size_t> dIndex_temp;

                for (size_t dIndex : DerivativesOfLeftBoundaryThisVariable)
                {
                    size_t stateIndex = std::get<1>(Derivatives_of_StateBeforePhase[dIndex]);

                    if (stateIndex < 7)
                    {
                        size_t Xindex = std::get<0>(Derivatives_of_StateBeforePhase[dIndex]);

                        for (size_t PropagatedStateIndex : {0, 1, 2, 3, 4, 5, 6, 8, 9})
                        {
                            this->Derivatives_of_state_after_initial_coast.push_back({ Xindex, PropagatedStateIndex, 0.0 });
                            dIndex_temp.push_back(this->Derivatives_of_state_after_initial_coast.size() - 1);
                        }

                        break;
                    }
                }

                if (dIndex_temp.size() > 0)
                    this->DerivativesOfStateAfterInitialCoastByVariable.push_back(dIndex_temp);
            }
            //time
            for (std::vector<size_t> DerivativesOfLeftBoundaryThisTimeVariable : this->DerivativesOfLeftBoundaryByTimeVariable)
            {
                std::vector<size_t> dIndex_temp;

                for (size_t dIndex : DerivativesOfLeftBoundaryThisTimeVariable)
                {
                    size_t stateIndex = std::get<1>(Derivatives_of_StateBeforePhase_wrt_Time[dIndex]);

                    if (stateIndex < 7)
                    {
                        size_t Xindex = std::get<0>(Derivatives_of_StateBeforePhase_wrt_Time[dIndex]);

                        for (size_t PropagatedStateIndex : {0, 1, 2, 3, 4, 5, 6, 8, 9})
                        {
                            this->Derivatives_of_state_after_initial_coast_wrt_Time.push_back({ Xindex, PropagatedStateIndex, 0.0 });
                            dIndex_temp.push_back(this->Derivatives_of_state_after_initial_coast_wrt_Time.size() - 1);
                        }

                        break;
                    }
                }

                if (dIndex_temp.size() > 0)
                    this->DerivativesOfStateAfterInitialCoastByTimeVariable.push_back(dIndex_temp);
            }
            

            //Step 2.2: epoch w.r.t. things that affect epoch
            //state
            for (std::vector<size_t> DerivativesOfLeftBoundaryThisVariable : this->DerivativesOfLeftBoundaryByVariable)
            {
                std::vector<size_t> dIndex_temp;

                for (size_t dIndex : DerivativesOfLeftBoundaryThisVariable)
                {
                    size_t stateIndex = std::get<1>(Derivatives_of_StateBeforePhase[dIndex]);

                    if (stateIndex == 7)
                    {
                        size_t Xindex = std::get<0>(Derivatives_of_StateBeforePhase[dIndex]);
                        this->Derivatives_of_state_after_initial_coast.push_back({ Xindex, 7, 0.0 });
                        dIndex_temp.push_back(this->Derivatives_of_state_after_initial_coast.size() - 1);

                        break;
                    }
                }

                if (dIndex_temp.size() > 0)
                    this->DerivativesOfStateAfterInitialCoastByVariable.push_back(dIndex_temp);
            }
            //time
            for (std::vector<size_t> DerivativesOfLeftBoundaryThisTimeVariable : this->DerivativesOfLeftBoundaryByTimeVariable)
            {
                std::vector<size_t> dIndex_temp;

                for (size_t dIndex : DerivativesOfLeftBoundaryThisTimeVariable)
                {
                    size_t stateIndex = std::get<1>(Derivatives_of_StateBeforePhase_wrt_Time[dIndex]);

                    if (stateIndex == 7)
                    {
                        size_t Xindex = std::get<0>(Derivatives_of_StateBeforePhase_wrt_Time[dIndex]);
                        this->Derivatives_of_state_after_initial_coast_wrt_Time.push_back({ Xindex, 7, 0.0 });
                        dIndex_temp.push_back(this->Derivatives_of_state_after_initial_coast_wrt_Time.size() - 1);

                        break;
                    }
                }

                if (dIndex_temp.size() > 0)
                    this->DerivativesOfStateAfterInitialCoastByTimeVariable.push_back(dIndex_temp);
            }
            //the tanks do not have derivatives with respect to any decision variables for a forced coast!
        }//end calcbounds_initial_coast

        void ParallelShootingPhase::calcbounds_terminal_coast()
        {
            //StateBeforeTerminalCoast has derivatives with respect to EVERYTHING that affects the left boundary condition
            std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateAfterPhase = this->myArrivalEvent->get_Derivatives_of_StateBeforeEvent();//Xindex, stateIndex, derivative value
            std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateAfterPhase_wrt_Time = this->myArrivalEvent->get_Derivatives_of_StateBeforeEvent_wrt_Time();//Xindex, stateIndex, derivative value

            //Step 1: what does the Right boundary have derivatives with respect to?
            {
                //Step 1.1: state variables
                for (size_t dIndex = 0; dIndex < Derivatives_of_StateAfterPhase.size(); ++dIndex)
                {
                    bool alreadyListed = false;
                    size_t Xindex = std::get<0>(Derivatives_of_StateAfterPhase[dIndex]);
                    for (size_t listIndex : this->ListOfVariablesAffectingRightBoundary)
                    {
                        if (listIndex == Xindex)
                        {
                            alreadyListed = true;
                            break;
                        }
                    }
                    if (!alreadyListed)
                    {
                        this->ListOfVariablesAffectingRightBoundary.push_back(Xindex);
                    }
                }

                this->DerivativesOfRightBoundaryByVariable.resize(ListOfVariablesAffectingRightBoundary.size());

                for (size_t listIndex = 0; listIndex < this->ListOfVariablesAffectingRightBoundary.size(); ++listIndex)
                {
                    for (size_t dIndex = 0; dIndex < Derivatives_of_StateAfterPhase.size(); ++dIndex)
                    {
                        size_t Xindex = std::get<0>(Derivatives_of_StateAfterPhase[dIndex]);

                        if (Xindex == this->ListOfVariablesAffectingRightBoundary[listIndex])
                        {
                            this->DerivativesOfRightBoundaryByVariable[listIndex].push_back(dIndex);
                        }
                    }
                }

                //Step 1.2: time variables
                for (size_t dIndex = 0; dIndex < Derivatives_of_StateAfterPhase_wrt_Time.size(); ++dIndex)
                {
                    bool alreadyListed = false;
                    size_t Xindex = std::get<0>(Derivatives_of_StateAfterPhase_wrt_Time[dIndex]);
                    for (size_t listIndex : this->ListOfTimeVariablesAffectingRightBoundary)
                    {
                        if (listIndex == Xindex)
                        {
                            alreadyListed = true;
                            break;
                        }
                    }
                    if (!alreadyListed)
                    {
                        this->ListOfTimeVariablesAffectingRightBoundary.push_back(Xindex);
                    }
                }

                this->DerivativesOfRightBoundaryByTimeVariable.resize(ListOfTimeVariablesAffectingRightBoundary.size());

                for (size_t listIndex = 0; listIndex < this->ListOfTimeVariablesAffectingRightBoundary.size(); ++listIndex)
                {
                    for (size_t dIndex = 0; dIndex < Derivatives_of_StateAfterPhase_wrt_Time.size(); ++dIndex)
                    {
                        size_t Xindex = std::get<0>(Derivatives_of_StateAfterPhase_wrt_Time[dIndex]);

                        if (Xindex == this->ListOfTimeVariablesAffectingRightBoundary[listIndex])
                        {
                            this->DerivativesOfRightBoundaryByTimeVariable[listIndex].push_back(dIndex);
                        }
                    }
                }
            }

            //Step 2: make the derivative skeleton of StateBeforeTerminalCoast
            //Step 2.1: 7-state w.r.t. things that affect 7-state
            //state
            for (std::vector<size_t> DerivativesOfRightBoundaryThisVariable : this->DerivativesOfRightBoundaryByVariable)
            {
                std::vector<size_t> dIndex_temp;

                for (size_t dIndex : DerivativesOfRightBoundaryThisVariable)
                {
                    size_t stateIndex = std::get<1>(Derivatives_of_StateAfterPhase[dIndex]);

                    if (stateIndex < 7)
                    {
                        size_t Xindex = std::get<0>(Derivatives_of_StateAfterPhase[dIndex]);

                        for (size_t PropagatedStateIndex : {0, 1, 2, 3, 4, 5, 6, 8, 9})
                        {
                            this->Derivatives_of_state_before_terminal_coast.push_back({ Xindex, PropagatedStateIndex, 0.0 });
                            dIndex_temp.push_back(this->Derivatives_of_state_before_terminal_coast.size() - 1);
                        }

                        break;
                    }
                }

                if (dIndex_temp.size() > 0)
                    this->DerivativesOfStateBeforeTerminalCoastByVariable.push_back(dIndex_temp);
            }
            //time
            for (std::vector<size_t> DerivativesOfRightBoundaryThisTimeVariable : this->DerivativesOfRightBoundaryByTimeVariable)
            {
                std::vector<size_t> dIndex_temp;

                for (size_t dIndex : DerivativesOfRightBoundaryThisTimeVariable)
                {
                    size_t stateIndex = std::get<1>(Derivatives_of_StateAfterPhase_wrt_Time[dIndex]);

                    if (stateIndex < 7)
                    {
                        size_t Xindex = std::get<0>(Derivatives_of_StateAfterPhase_wrt_Time[dIndex]);

                        for (size_t PropagatedStateIndex : {0, 1, 2, 3, 4, 5, 6, 8, 9})
                        {
                            this->Derivatives_of_state_before_terminal_coast_wrt_Time.push_back({ Xindex, PropagatedStateIndex, 0.0 });
                            dIndex_temp.push_back(this->Derivatives_of_state_before_terminal_coast_wrt_Time.size() - 1);
                        }

                        break;
                    }
                }

                if (dIndex_temp.size() > 0)
                    this->DerivativesOfStateBeforeTerminalCoastByTimeVariable.push_back(dIndex_temp);
            }

            //Step 2.2: epoch w.r.t. things that affect epoch
            //state
            for (std::vector<size_t> DerivativesOfRightBoundaryThisVariable : this->DerivativesOfRightBoundaryByVariable)
            {
                std::vector<size_t> dIndex_temp;

                for (size_t dIndex : DerivativesOfRightBoundaryThisVariable)
                {
                    size_t stateIndex = std::get<1>(Derivatives_of_StateAfterPhase[dIndex]);

                    if (stateIndex == 7)
                    {
                        size_t Xindex = std::get<0>(Derivatives_of_StateAfterPhase[dIndex]);
                        this->Derivatives_of_state_before_terminal_coast.push_back({ Xindex, 7, 0.0 });
                        dIndex_temp.push_back(this->Derivatives_of_state_before_terminal_coast.size() - 1);

                        break;
                    }
                }

                if (dIndex_temp.size() > 0)
                    this->DerivativesOfStateBeforeTerminalCoastByVariable.push_back(dIndex_temp);
            }
            //time
            for (std::vector<size_t> DerivativesOfRightBoundaryThisTimeVariable : this->DerivativesOfRightBoundaryByTimeVariable)
            {
                std::vector<size_t> dIndex_temp;

                for (size_t dIndex : DerivativesOfRightBoundaryThisTimeVariable)
                {
                    size_t stateIndex = std::get<1>(Derivatives_of_StateAfterPhase_wrt_Time[dIndex]);

                    if (stateIndex == 7)
                    {
                        size_t Xindex = std::get<0>(Derivatives_of_StateAfterPhase_wrt_Time[dIndex]);
                        this->Derivatives_of_state_before_terminal_coast_wrt_Time.push_back({ Xindex, 7, 0.0 });
                        dIndex_temp.push_back(this->Derivatives_of_state_before_terminal_coast_wrt_Time.size() - 1);

                        break;
                    }
                }

                if (dIndex_temp.size() > 0)
                    this->DerivativesOfStateBeforeTerminalCoastByTimeVariable.push_back(dIndex_temp);
            }

            //the tanks DO have derivatives with respect to their own virtual tank variables!
            this->Derivatives_of_state_before_terminal_coast.push_back({ this->Xindex_virtual_chemical_fuel_tank, 8, 1.0 });
            this->Derivatives_of_state_before_terminal_coast.push_back({ this->Xindex_virtual_electric_propellant_tank, 9, 1.0 });
            //the tanks do not have derivatives with respect to any decision variables for a forced coast!
        }//end calcbounds_terminal_coast

        void ParallelShootingPhase::calcbounds_fixed_inertial_control()
        {
            if (this->myJourneyOptions->num_interior_control_points > 1)
            {
                throw std::invalid_argument(this->prefix + ": fixed intertial control is incompatible with any derived class of ParallelShootingPhase if num_interior_control_points > 1. Set num_interior_control_points = 1 and the constraint will work.");
            }

            std::vector<std::string> controlNames({ "x", "y", "z", "command" });
            std::vector<size_t> Xindex_control_first_step = this->mySteps[0].getXindices_control(0);

            for (size_t stepIndex = 1; stepIndex < this->num_steps; ++stepIndex)
            {
                std::vector<size_t> step_Gindices_fixed_inertial_control_wrt_this_step_control;
                std::vector<size_t> step_Gindices_fixed_inertial_control_wrt_first_step_control;

                std::vector<size_t> Xindex_control_this_step = this->mySteps[stepIndex].getXindices_control(0);
                

                for (size_t controlIndex = 0; controlIndex < Xindex_control_this_step.size(); ++controlIndex)
                {
                    //u_i - u_0 = 0.0
                    //implement as nonlinear constraint so that SNOPT can temporarily violate it during the solve
                    this->Flowerbounds->push_back(-math::SMALL);
                    this->Fupperbounds->push_back(math::SMALL);
                    this->Fdescriptions->push_back(this->prefix + "Step " + std::to_string(stepIndex) + " u_" + controlNames[controlIndex] + " must match step 0");
                    this->create_sparsity_entry(this->Fdescriptions->size() - 1, Xindex_control_this_step[controlIndex], step_Gindices_fixed_inertial_control_wrt_this_step_control);
                    this->create_sparsity_entry(this->Fdescriptions->size() - 1, Xindex_control_first_step[controlIndex], step_Gindices_fixed_inertial_control_wrt_first_step_control);
                }

                this->Gindices_fixed_inertial_control_wrt_this_step_control.push_back(step_Gindices_fixed_inertial_control_wrt_this_step_control);
                this->Gindices_fixed_inertial_control_wrt_first_step_control.push_back(step_Gindices_fixed_inertial_control_wrt_first_step_control);
            }//end loop over steps
        }//end calcbounds_fixed_inertial_control

        void ParallelShootingPhase::process_phase(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            this->process_phase_left_boundary(X, Xindex, F, Findex, G, needG);

            this->process_phase_flight_time(X, Xindex, F, Findex, G, needG);

            this->process_phase_right_boundary(X, Xindex, F, Findex, G, needG);

            this->process_virtual_propellant_tanks(X, Xindex, F, Findex, G, needG);

            this->process_phase_main(X, Xindex, F, Findex, G, needG);

            //do we want fixed inertial control?
            if (this->myJourneyOptions->force_fixed_inertial_control)
                this->process_fixed_inertial_control(X, Xindex, F, Findex, G, needG);
        }//end process_phase

        void ParallelShootingPhase::process_virtual_propellant_tanks(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG)
        {
            this->virtual_chemical_fuel_used = X[Xindex++];
            this->virtual_electric_propellant_used = X[Xindex++];
        }//end process_virtual_propellant_tanks()

        void ParallelShootingPhase::process_phase_main(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            this->process_phase_initial_coast(X, Xindex, F, Findex, G, needG);

            this->process_phase_terminal_coast(X, Xindex, F, Findex, G, needG);

            //steps
            for (size_t stepIndex = 0; stepIndex < this->num_steps; ++stepIndex)
                this->mySteps[stepIndex].process_step(X, Xindex, F, Findex, G, needG);
        }//end process_phase_main

        void ParallelShootingPhase::process_phase_initial_coast(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //Step 1: initial tank state
            //set the tank states
            this->state_after_initial_TCM(8) = this->chemical_fuel_used; //TCM propellant, usually zero
            this->state_after_initial_TCM(9) = 0.0; //electric propellant

            //Step 2: initial coast
            if (this->hasInitialCoast)
            {
                this->InitialCoast_dStatedIndependentVariable.assign_zeros();
                this->InitialCoastPropagatorObject->setCurrentEpoch(this->state_after_initial_TCM(7));
                this->InitialCoastPropagatorObject->setIndexOfEpochInStateVec(7);
                this->InitialCoastPropagatorObject->setCurrentIndependentVariable(this->state_after_initial_TCM(7));
                this->InitialCoastPropagatorObject->propagate(this->InitialCoastDuration, needG);
                this->state_after_initial_coast(7) = this->state_after_initial_TCM(7) + this->InitialCoastDuration; 
                
                //form augmented STM
                if (needG)
                {
                    for (size_t i = 0; i < 7; ++i)
                    {
                        for (size_t j = 0; j < 7; ++j)
                            this->STM_Augmented_initial_coast(i, j) = this->STM_initial_coast(i, j);

                        //Phi_t terms
                        this->STM_Augmented_initial_coast(i, 7) = this->STM_initial_coast(i, 7);
                        this->STM_Augmented_initial_coast(i, this->stateIndex_phase_propagation_variable) = this->STM_initial_coast(i, 13);
                    }

                    //tanks
                    for (size_t i = 8; i < 10; ++i)
                    {
                        for (size_t j = 0; j < 7; ++j)
                            this->STM_Augmented_initial_coast(i, j) = this->STM_initial_coast(i, j);

                        //Phi_t terms
                        this->STM_Augmented_initial_coast(i, 7) = this->STM_initial_coast(i, 7);
                        this->STM_Augmented_initial_coast(i, this->stateIndex_phase_propagation_variable) = this->STM_initial_coast(i, 13);
                    }
                }//end derivatives
            }
            else
                this->state_after_initial_coast = this->state_after_initial_TCM;

            //Step 3: populate derivative tuples
            if (needG)
            {
                //Step 3.1: get derivatives from the boundary event
                std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateBeforePhase = this->myDepartureEvent->get_Derivatives_of_StateAfterEvent();//Xindex, stateIndex, derivative value
                std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateBeforePhase_wrt_Time = this->myDepartureEvent->get_Derivatives_of_StateAfterEvent_wrt_Time();//Xindex, stateIndex, derivative value

                //Step 3.2: cycle through decision variables that affect the states, find the derivatives of the post-coast state, and then put it in a useful place
                math::Matrix<double> dBoundaryState_dDecisionVariable(14, 1, 0.0);
                math::Matrix<double> dStateAfterInitialCoast_dDecisionVariable(14, 1, 0.0);

                //Step 3.2.1: non-time
                for (std::vector<size_t> DerivativesOfLeftBoundaryThisVariable : this->DerivativesOfLeftBoundaryByVariable)
                {
                    dBoundaryState_dDecisionVariable.assign_zeros();
                    size_t Xindex = std::get<0>(Derivatives_of_StateBeforePhase[DerivativesOfLeftBoundaryThisVariable[0]]);

                    //Step 3.2.1.1: populate dBoundaryState_dDecisionVariable for this variable
                    for (size_t dIndex : DerivativesOfLeftBoundaryThisVariable)
                    {
                        size_t stateIndex = std::get<1>(Derivatives_of_StateBeforePhase[dIndex]);
                        dBoundaryState_dDecisionVariable(stateIndex) = std::get<2>(Derivatives_of_StateBeforePhase[dIndex]);
                    }

                    //Step 3.2.1.2: chain!
                    if (this->hasInitialCoast)
                        dStateAfterInitialCoast_dDecisionVariable = this->STM_Augmented_initial_coast * dBoundaryState_dDecisionVariable;
                    else
                        dStateAfterInitialCoast_dDecisionVariable = dBoundaryState_dDecisionVariable;

                    //Step 3.2.1.3: put the stuff in the place
                    for (std::vector<size_t> DerivativesOfStateAfterInitialCoast : this->DerivativesOfStateAfterInitialCoastByVariable)
                    {
                        size_t cXindex = std::get<0>(this->Derivatives_of_state_after_initial_coast[DerivativesOfStateAfterInitialCoast[0]]);

                        if (cXindex == Xindex)
                        {
                            for (size_t dIndex : DerivativesOfStateAfterInitialCoast)
                            {
                                size_t stateIndex = std::get<1>(this->Derivatives_of_state_after_initial_coast[dIndex]);
                                std::get<2>(this->Derivatives_of_state_after_initial_coast[dIndex]) = dStateAfterInitialCoast_dDecisionVariable(stateIndex);
                            }
                        }
                    }
                }//end non-time

                //Step 3.2.2: time
                for (std::vector<size_t> DerivativesOfLeftBoundaryThisTimeVariable : this->DerivativesOfLeftBoundaryByTimeVariable)
                {
                    dBoundaryState_dDecisionVariable.assign_zeros();
                    size_t Xindex = std::get<0>(Derivatives_of_StateBeforePhase_wrt_Time[DerivativesOfLeftBoundaryThisTimeVariable[0]]);

                    //Step 3.2.2.1: populate dBoundaryState_dDecisionVariable for this variable
                    for (size_t dIndex : DerivativesOfLeftBoundaryThisTimeVariable)
                    {
                        size_t stateIndex = std::get<1>(Derivatives_of_StateBeforePhase_wrt_Time[dIndex]);
                        dBoundaryState_dDecisionVariable(stateIndex) = std::get<2>(Derivatives_of_StateBeforePhase_wrt_Time[dIndex]);
                    }

                    //Step 3.2.2.2: chain!
                    if (this->hasInitialCoast)
                        dStateAfterInitialCoast_dDecisionVariable = this->STM_Augmented_initial_coast * dBoundaryState_dDecisionVariable;
                    else
                        dStateAfterInitialCoast_dDecisionVariable = dBoundaryState_dDecisionVariable;

                    //Step 3.2.2.3: put the stuff in the place
                    for (std::vector<size_t> DerivativesOfStateAfterInitialCoast : this->DerivativesOfStateAfterInitialCoastByTimeVariable)
                    {
                        size_t cXindex = std::get<0>(this->Derivatives_of_state_after_initial_coast_wrt_Time[DerivativesOfStateAfterInitialCoast[0]]);

                        if (cXindex == Xindex)
                        {
                            for (size_t dIndex : DerivativesOfStateAfterInitialCoast)
                            {
                                size_t stateIndex = std::get<1>(this->Derivatives_of_state_after_initial_coast_wrt_Time[dIndex]);
                                std::get<2>(this->Derivatives_of_state_after_initial_coast_wrt_Time[dIndex]) = dStateAfterInitialCoast_dDecisionVariable(stateIndex);
                            }
                        }
                    }
                }//end non-time
            }//end derivatives
        }//end process_phase_initial_coast

        void ParallelShootingPhase::process_phase_terminal_coast(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //Step 1: tank states
            this->state_at_end_of_phase(8) = this->virtual_chemical_fuel_used;
            this->state_at_end_of_phase(9) = this->virtual_electric_propellant_used;

            //Step 2: terminal coast
            if (this->hasTerminalCoast)
            {
                this->TerminalCoast_dStatedIndependentVariable.assign_zeros();
                this->TerminalCoastPropagatorObject->setCurrentEpoch(this->state_at_end_of_phase(7));
                this->TerminalCoastPropagatorObject->setIndexOfEpochInStateVec(7);
                this->TerminalCoastPropagatorObject->setCurrentIndependentVariable(this->state_at_end_of_phase(7));
                this->TerminalCoastPropagatorObject->propagate(-this->TerminalCoastDuration, needG);
                this->state_before_terminal_coast(7) = this->state_at_end_of_phase(7) - (this->TerminalCoastDuration);

                //form augmented STM
                if (needG)
                {
                    for (size_t i = 0; i < 7; ++i)
                    {
                        for (size_t j = 0; j < 7; ++j)
                            this->STM_Augmented_terminal_coast(i, j) = this->STM_terminal_coast(i, j);

                        //Phi_t terms
                        this->STM_Augmented_terminal_coast(i, 7) = this->STM_terminal_coast(i, 7);
                        this->STM_Augmented_terminal_coast(i, this->stateIndex_phase_propagation_variable) = this->STM_terminal_coast(i, 13);
                    }

                    //tanks
                    for (size_t i = 8; i < 10; ++i)
                    {
                        for (size_t j = 0; j < 7; ++j)
                            this->STM_Augmented_terminal_coast(i, j) = this->STM_terminal_coast(i, j);

                        //Phi_t terms
                        this->STM_Augmented_terminal_coast(i, 7) = this->STM_terminal_coast(i, 7);
                        this->STM_Augmented_terminal_coast(i, this->stateIndex_phase_propagation_variable) = this->STM_terminal_coast(i, 13);
                    }
                }
            }
            else
                this->state_before_terminal_coast = this->state_at_end_of_phase;


            //Step 3: populate derivative tuples
            if (needG)
            {
                //Step 3.1: get derivatives from the boundary event
                std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateAfterPhase = this->myArrivalEvent->get_Derivatives_of_StateBeforeEvent();//Xindex, stateIndex, derivative value
                std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateAfterPhase_wrt_Time = this->myArrivalEvent->get_Derivatives_of_StateBeforeEvent_wrt_Time();//Xindex, stateIndex, derivative value

                                                                                                                                                                                    //Step 3.2: cycle through decision variables that affect the states, find the derivatives of the post-coast state, and then put it in a useful place
                math::Matrix<double> dBoundaryState_dDecisionVariable(14, 1, 0.0);
                math::Matrix<double> dStateBeforeTerminalCoast_dDecisionVariable(14, 1, 0.0);

                //Step 3.2.1: non-time
                for (std::vector<size_t> DerivativesOfRightBoundaryThisVariable : this->DerivativesOfRightBoundaryByVariable)
                {
                    dBoundaryState_dDecisionVariable.assign_zeros();
                    size_t Xindex = std::get<0>(Derivatives_of_StateAfterPhase[DerivativesOfRightBoundaryThisVariable[0]]);

                    //Step 3.2.1.1: populate dBoundaryState_dDecisionVariable for this variable
                    for (size_t dIndex : DerivativesOfRightBoundaryThisVariable)
                    {
                        size_t stateIndex = std::get<1>(Derivatives_of_StateAfterPhase[dIndex]);
                        dBoundaryState_dDecisionVariable(stateIndex) = std::get<2>(Derivatives_of_StateAfterPhase[dIndex]);
                    }

                    //Step 3.2.1.2: chain!
                    if (this->hasTerminalCoast)
                        dStateBeforeTerminalCoast_dDecisionVariable = this->STM_Augmented_terminal_coast * dBoundaryState_dDecisionVariable;
                    else
                        dStateBeforeTerminalCoast_dDecisionVariable = dBoundaryState_dDecisionVariable;

                    //Step 3.2.1.3: put the stuff in the place
                    for (std::vector<size_t> DerivativesOfStateBeforeTerminalCoast : this->DerivativesOfStateBeforeTerminalCoastByVariable)
                    {
                        size_t cXindex = std::get<0>(this->Derivatives_of_state_before_terminal_coast[DerivativesOfStateBeforeTerminalCoast[0]]);

                        if (cXindex == Xindex)
                        {
                            for (size_t dIndex : DerivativesOfStateBeforeTerminalCoast)
                            {
                                size_t stateIndex = std::get<1>(this->Derivatives_of_state_before_terminal_coast[dIndex]);
                                std::get<2>(this->Derivatives_of_state_before_terminal_coast[dIndex]) = dStateBeforeTerminalCoast_dDecisionVariable(stateIndex);
                            }
                        }
                    }
                }//end non-time

                 //Step 3.2.2: time
                for (std::vector<size_t> DerivativesOfRightBoundaryThisTimeVariable : this->DerivativesOfRightBoundaryByTimeVariable)
                {
                    dBoundaryState_dDecisionVariable.assign_zeros();
                    size_t Xindex = std::get<0>(Derivatives_of_StateAfterPhase_wrt_Time[DerivativesOfRightBoundaryThisTimeVariable[0]]);

                    //Step 3.2.2.1: populate dBoundaryState_dDecisionVariable for this variable
                    for (size_t dIndex : DerivativesOfRightBoundaryThisTimeVariable)
                    {
                        size_t stateIndex = std::get<1>(Derivatives_of_StateAfterPhase_wrt_Time[dIndex]);
                        dBoundaryState_dDecisionVariable(stateIndex) = std::get<2>(Derivatives_of_StateAfterPhase_wrt_Time[dIndex]);
                    }
                    
                    //Step 3.2.2.2: chain!
                    if (this->hasTerminalCoast)
                        dStateBeforeTerminalCoast_dDecisionVariable = this->STM_Augmented_terminal_coast * dBoundaryState_dDecisionVariable;
                    else
                        dStateBeforeTerminalCoast_dDecisionVariable = dBoundaryState_dDecisionVariable;

                    //Step 3.2.2.3: put the stuff in the place
                    for (std::vector<size_t> DerivativesOfStateBeforeTerminalCoast : this->DerivativesOfStateBeforeTerminalCoastByTimeVariable)
                    {
                        size_t cXindex = std::get<0>(this->Derivatives_of_state_before_terminal_coast_wrt_Time[DerivativesOfStateBeforeTerminalCoast[0]]);

                        if (cXindex == Xindex)
                        {
                            for (size_t dIndex : DerivativesOfStateBeforeTerminalCoast)
                            {
                                size_t stateIndex = std::get<1>(this->Derivatives_of_state_before_terminal_coast_wrt_Time[dIndex]);
                                std::get<2>(this->Derivatives_of_state_before_terminal_coast_wrt_Time[dIndex]) = dStateBeforeTerminalCoast_dDecisionVariable(stateIndex);
                            }
                        }
                    }
                }//end non-time
            }//end derivatives
        }//end process_phase_terminal_coast

        void ParallelShootingPhase::process_deltav_contribution(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //compute the total delta-v
            this->PhaseTotalDeterministicDeltav = this->myDepartureEvent->get_DeterministicDeltav()
                + this->myArrivalEvent->get_DeterministicDeltav();

            for (size_t step = 0; step < this->num_steps; ++step)
                this->PhaseTotalDeterministicDeltav += this->mySteps[step].getstepDeltav();

            this->PhaseTotalStatisticalDeltav = this->myDepartureEvent->get_StatisticalDeltav()
                + this->myArrivalEvent->get_StatisticalDeltav()
                + this->initial_TCM_magnitude;

            //TODO: derivatives
            if (this->myOptions->objective_type == ObjectiveFunctionType::MINIMIZE_DELTAV
                && needG)
            {
                //do stuff
            }
        }//end process_deltav_constribution

        void ParallelShootingPhase::process_fixed_inertial_control(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {

            std::vector<size_t> Xindex_control_first_step = this->mySteps[0].getXindices_control(0);

            for (size_t stepIndex = 1; stepIndex < this->num_steps; ++stepIndex)
            {
                std::vector<size_t> Xindex_control_this_step = this->mySteps[stepIndex].getXindices_control(0);

                for (size_t controlIndex = 0; controlIndex < Xindex_control_this_step.size(); ++controlIndex)
                {
                    //u_i - u_0 = 0.0
                    //implement as nonlinear constraint so that SNOPT can temporarily violate it during the solve
                    F[Findex++] = X[Xindex_control_this_step[controlIndex]] - X[Xindex_control_first_step[controlIndex]];

                    G[this->Gindices_fixed_inertial_control_wrt_this_step_control[stepIndex - 1][controlIndex]] = 1.0;
                    G[this->Gindices_fixed_inertial_control_wrt_first_step_control[stepIndex - 1][controlIndex]] = -1.0;
                }
            }//end loop over steps
        }//end process_fixed_inertial_control()

        void ParallelShootingPhase::output(std::ofstream& outputfile,
            size_t& eventcount)
        {
            //Step 1: output the departure event
            this->myDepartureEvent->output(outputfile, this->LaunchDate _GETVALUE, eventcount);
            this->mySpacecraft->setActiveStage(this->stageIndex);

            //Step 2: output the initial TCM if applicable
            this->output_initial_TCM(outputfile, eventcount);

            //Step 3: output the state at the middle of the initial coast, if applicable
            if (this->hasInitialCoast)
            {
                //Step 3.1: create a temporary output state and assign the initial coast propagator to it
                this->InitialCoastPropagatorObject->setStateRight(output_state);

                //Step 3.2: propagate
                this->InitialCoast_dStatedIndependentVariable.assign_zeros();
				this->InitialCoastPropagatorObject->setCurrentEpoch(this->state_after_initial_TCM(7));
				this->InitialCoastPropagatorObject->setIndexOfEpochInStateVec(7);
				this->InitialCoastPropagatorObject->setCurrentIndependentVariable(this->state_after_initial_TCM(7));
                this->InitialCoastPropagatorObject->propagate(this->InitialCoastDuration / 2.0, false);
                output_state(6) = this->state_after_initial_TCM(6);
                output_state(7) = this->state_after_initial_TCM(7) + this->InitialCoastDuration / 2.0;

                //Step 3.3: assign the initial coast propagator back to where it is supposed to go
                this->InitialCoastPropagatorObject->setStateRight(this->state_after_initial_coast);

                //Step 3.4: figure out spacecrafty things

                //Step 3.4.1: where am I relative to the sun?
                math::Matrix<doubleType> R_sc_Sun(3, 1, 0.0);
                if (this->myUniverse->central_body_SPICE_ID == 10)
                {
                    R_sc_Sun = output_state.getSubMatrix1D(0, 2);
                }
                else
                {
                    //where is the central body relative to the sun?
                    doubleType central_body_state_and_derivatives[12];
                    this->myUniverse->locate_central_body(output_state(7),
                        central_body_state_and_derivatives,
                        *this->myOptions,
                        false);

                    math::Matrix<doubleType> R_CB_Sun(3, 1, 0.0);
                    for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
                    {
                        R_CB_Sun(stateIndex) = central_body_state_and_derivatives[stateIndex];
                    }

                    R_sc_Sun = output_state.getSubMatrix1D(0, 2) + R_CB_Sun;
                }
                doubleType r_sc_sun_AU = R_sc_Sun.norm() / this->myOptions->AU;
                this->mySpacecraft->computePowerState(r_sc_sun_AU, output_state(7));

                doubleType power = this->mySpacecraft->getAvailablePower();
                doubleType active_power = this->mySpacecraft->getEPActivePower();

                //Step 3.5: print
                phase::write_output_line(outputfile,//outputfile
                    eventcount,//eventcount
                    "force-coast",//event_type
                    "deep-space",//event_location
                    this->InitialCoastDuration / 86400.0,// timestep_size,
                    -1,//flyby_altitude,
                    0,//BdotR
                    0,//BdotT
                    0,//angle1
                    0,//angle2
                    0,//C3
                    output_state,//state
                    math::Matrix<doubleType>(3, 1, 0.0),//dV
                    math::Matrix<doubleType>(3, 1, 0.0),//ThrustVector
                    0.0,//dVmag
                    0.0,//Thrust
                    0.0,//Isp
                    power,//AvailPower
                    0.0,//mdot
                    0,//number_of_active_engines
                    active_power,//active power
                    "none");//throttle level
            }

            //Step 4: output the steps
            for (size_t step = 0; step < this->num_steps; ++step)
                this->mySteps[step].output(outputfile, eventcount);

            //Step 5: output the state at the middle of the terminal coast, if applicable
            if (this->hasTerminalCoast)
            {
                //Step 5.1: create a temporary output state and assign the Terminal coast propagator to it
                this->TerminalCoastPropagatorObject->setStateRight(output_state);

                //Step 5.2: propagate
                this->TerminalCoast_dStatedIndependentVariable.assign_zeros();
				this->TerminalCoastPropagatorObject->setCurrentEpoch(this->state_at_end_of_phase(7));
				this->TerminalCoastPropagatorObject->setIndexOfEpochInStateVec(7);
				this->TerminalCoastPropagatorObject->setCurrentIndependentVariable(this->state_at_end_of_phase(7));
                this->TerminalCoastPropagatorObject->propagate(-this->TerminalCoastDuration / 2.0, false);
                output_state(6) = this->state_at_end_of_phase(6);
                output_state(7) = this->state_at_end_of_phase(7) - this->TerminalCoastDuration / 2.0;

                //Step 5.3: assign the Terminal coast propagator back to where it is supposed to go
                this->TerminalCoastPropagatorObject->setStateRight(this->state_before_terminal_coast);

                //Step 5.4: figure out spacecrafty things

                //Step 5.4.1: where am I relative to the sun?
                math::Matrix<doubleType> R_sc_Sun(3, 1, 0.0);
                if (this->myUniverse->central_body_SPICE_ID == 10)
                {
                    R_sc_Sun = output_state.getSubMatrix1D(0, 2);
                }
                else
                {
                    //where is the central body relative to the sun?
                    doubleType central_body_state_and_derivatives[12];
                    this->myUniverse->locate_central_body(output_state(7),
                        central_body_state_and_derivatives,
                        *this->myOptions,
                        false);

                    math::Matrix<doubleType> R_CB_Sun(3, 1, 0.0);
                    for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
                    {
                        R_CB_Sun(stateIndex) = central_body_state_and_derivatives[stateIndex];
                    }

                    R_sc_Sun = output_state.getSubMatrix1D(0, 2) + R_CB_Sun;
                }
                doubleType r_sc_sun_AU = R_sc_Sun.norm() / this->myOptions->AU;
                this->mySpacecraft->computePowerState(r_sc_sun_AU, output_state(7));

                doubleType power = this->mySpacecraft->getAvailablePower();
                doubleType active_power = this->mySpacecraft->getEPActivePower();

                //Step 5.5: print
                phase::write_output_line(outputfile,//outputfile
                    eventcount,//eventcount
                    "force-coast",//event_type
                    "deep-space",//event_location
                    this->TerminalCoastDuration / 86400.0,// timestep_size,
                    -1,//flyby_altitude,
                    0,//BdotR
                    0,//BdotT
                    0,//angle1
                    0,//angle2
                    0,//C3
                    output_state,//state
                    math::Matrix<doubleType>(3, 1, 0.0),//dV
                    math::Matrix<doubleType>(3, 1, 0.0),//ThrustVector
                    0.0,//dVmag
                    0.0,//Thrust
                    0.0,//Isp
                    power,//AvailPower
                    0.0,//mdot
                    0,//number_of_active_engines
                    active_power,//active_power
                    "none");//throttle level
            }

            //Step 6: output the arrival event
            this->myArrivalEvent->output(outputfile, this->LaunchDate _GETVALUE, eventcount);
        }//end output

        void ParallelShootingPhase::output_ephemeris(std::ofstream& outputfile, std::ofstream& acceleration_model_file)
        {
			// Reset the spacecraft
            this->mySpacecraft->setActiveStage(this->stageIndex);
			
            this->EphemerisOutputResolution = this->myOptions->integration_time_step_size;

            //Step 1: output the departure event
            this->myDepartureEvent->output_ephemeris(outputfile);

            //Step 2: output the initial coast
            if (this->hasInitialCoast)
            {
                //Step 2.1: temporarily assign the initial coast propagator to the output state
                this->InitialCoastPropagatorObject->setStateRight(output_state);

                //Step 2.2: propagate and print, skipping the first and last entry
                double timeToPropagate = this->EphemerisOutputResolution;
                double totalPropagationTime = this->InitialCoastDuration;
                while (timeToPropagate < totalPropagationTime)
                {
                    //Step 2.2.1: propagate
                    this->InitialCoast_dStatedIndependentVariable.assign_zeros();
					this->InitialCoastPropagatorObject->setCurrentEpoch(this->state_after_initial_TCM(7));
					this->InitialCoastPropagatorObject->setIndexOfEpochInStateVec(7);
					this->InitialCoastPropagatorObject->setCurrentIndependentVariable(this->state_after_initial_TCM(7));
                    this->InitialCoastPropagatorObject->propagate(timeToPropagate, false);
                    output_state(7) = this->state_after_initial_TCM(7) + timeToPropagate;

                    //Step 2.2.2: convert to Sun-centered if necessary
                    if (this->myUniverse->central_body.spice_ID != this->myOptions->forward_integrated_ephemeris_central_body_SPICE_ID)
                    {
                        double LT_dump;
                        double bodyStateDouble[6];
                        spkez_c(this->myUniverse->central_body.spice_ID, output_state(7)_GETVALUE - (51544.5 * 86400.0), "J2000", "NONE", this->myOptions->forward_integrated_ephemeris_central_body_SPICE_ID, bodyStateDouble, &LT_dump);

                        for (size_t stateIndex = 0; stateIndex < 6; ++stateIndex)
                            output_state(stateIndex) += bodyStateDouble[stateIndex];
                    }

                    //Step 2.2.3: print
                    if ((totalPropagationTime - timeToPropagate) > this->EphemerisOutputResolution)
                        this->write_ephemeris_line(outputfile,
                            output_state,
                            math::Matrix<doubleType>(3, 1, 0.0),//control vector
                            0.0,
                            0.0,
                            0.0,
                            0,
                            0.0,
                            "none");

                    //Step 2.2.4: increment propagatedEpoch
                    timeToPropagate += (totalPropagationTime - timeToPropagate) > this->EphemerisOutputResolution
                        ? this->EphemerisOutputResolution
                        : (totalPropagationTime - timeToPropagate);
                }

                //Step 2.3: assign the initial coast propagator back to where it is supposed to go
                this->InitialCoastPropagatorObject->setStateRight(this->state_after_initial_coast);
            }

            //Step 3: output the thrust arcs
            //forward step
            for (size_t step = 0; step < this->num_steps; ++step)
                this->mySteps[step].output_ephemeris(outputfile, acceleration_model_file);

             //Step 4: output the terminal coast
            if (this->hasTerminalCoast)
            {
                //Step 4.1: temporarily assign the initial coast propagator to the output state
                this->TerminalCoastPropagatorObject->setStateLeft(this->state_before_terminal_coast);
                this->TerminalCoastPropagatorObject->setStateRight(output_state);

                //Step 4.2: propagate and print, skipping the first entry
                double timeToPropagate = this->EphemerisOutputResolution;
                double totalPropagationTime = this->TerminalCoastDuration;
                while (timeToPropagate < totalPropagationTime)
                {
                    //Step 4.2.1: propagate
                    this->TerminalCoast_dStatedIndependentVariable.assign_zeros();
					this->TerminalCoastPropagatorObject->setCurrentEpoch(this->state_at_end_of_phase(7));
					this->TerminalCoastPropagatorObject->setIndexOfEpochInStateVec(7);
					this->TerminalCoastPropagatorObject->setCurrentIndependentVariable(this->state_at_end_of_phase(7));
                    this->TerminalCoastPropagatorObject->propagate(timeToPropagate, false);
                    output_state(7) = this->state_before_terminal_coast(7) + timeToPropagate;

                    //Step 4.2.2: convert to Sun-centered if necessary
                    if (this->myUniverse->central_body.spice_ID != this->myOptions->forward_integrated_ephemeris_central_body_SPICE_ID)
                    {
                        double LT_dump;
                        double bodyStateDouble[6];
                        spkez_c(this->myUniverse->central_body.spice_ID, output_state(7)_GETVALUE - (51544.5 * 86400.0), "J2000", "NONE", this->myOptions->forward_integrated_ephemeris_central_body_SPICE_ID, bodyStateDouble, &LT_dump);

                        for (size_t stateIndex = 0; stateIndex < 6; ++stateIndex)
                            output_state(stateIndex) += bodyStateDouble[stateIndex];
                    }

                    //Step 4.2.3: print
                    if ((totalPropagationTime - timeToPropagate) > this->EphemerisOutputResolution)
                        this->write_ephemeris_line(outputfile,
                            output_state,
                            math::Matrix<doubleType>(3, 1, 0.0),//control vector
                            0.0,
                            0.0,
                            0.0,
                            0,
                            0.0,
                            "none");

                    //Step 4.2.4: increment propagatedEpoch
                    timeToPropagate += (totalPropagationTime - timeToPropagate) > this->EphemerisOutputResolution
                        ? this->EphemerisOutputResolution
                        : (totalPropagationTime - timeToPropagate);
                }

                //Step 4.3: assign the terminal coast propagator back to where it is supposed to go
                this->TerminalCoastPropagatorObject->setStateLeft(this->state_at_end_of_phase);
                this->TerminalCoastPropagatorObject->setStateRight(this->state_before_terminal_coast);
            }

            //Step 5: output the arrival event
            this->myArrivalEvent->output_ephemeris(outputfile);
        }//end output_ephemeris

        void ParallelShootingPhase::output_STMs()
        {
			// Reset the spacecraft
            this->mySpacecraft->setActiveStage(this->stageIndex);
			
            std::ofstream summaryfile(this->myOptions->working_directory + "//" + boost::replace_all_copy(this->prefix, ": ", "") + ".stm_description");
            for (size_t stateIndex = 0; stateIndex < this->stateVectorNames.size(); ++stateIndex)
                summaryfile << this->stateVectorNames[stateIndex] << std::endl;
            summaryfile.close();

            size_t STMindex = 0;

            this->STM_Augmented_initial_coast.print_to_file(this->myOptions->working_directory + "//" + boost::replace_all_copy(this->prefix, ": ", "") + "_InitialCoast_STM_" + std::to_string(STMindex++) + "_forward.stm");

            for (size_t stepIndex = 0; stepIndex < this->num_steps; ++stepIndex)
                this->mySteps[stepIndex].output_STM(STMindex);

            this->STM_Augmented_terminal_coast.print_to_file(this->myOptions->working_directory + "//" + boost::replace_all_copy(this->prefix, ": ", "") + "_TerminalCoast_STM_" + std::to_string(STMindex++) + "_backward.stm");
        }//end output_STMs

        void ParallelShootingPhase::output_maneuver_and_target_spec(std::ofstream& maneuver_spec_file, std::ofstream& target_spec_file, bool& haveManeuverNeedTarget)
        {
			// Reset the spacecraft
            this->mySpacecraft->setActiveStage(this->stageIndex);

            //should write out a maneuver/target spec for the departure maneuver if there is one
            this->myDepartureEvent->output_maneuver_and_target_spec(maneuver_spec_file, target_spec_file, haveManeuverNeedTarget);
			
            //write out the step maneuver and target spec lines
            for (size_t stepIndex = 0; stepIndex < this->num_steps; ++stepIndex)
                this->mySteps[stepIndex].output_maneuver_and_target_spec(maneuver_spec_file, target_spec_file, haveManeuverNeedTarget);
                       
            //should write out a maneuver/target spec for the arrival maneuver if there is one
            this->myArrivalEvent->output_maneuver_and_target_spec(maneuver_spec_file, target_spec_file, haveManeuverNeedTarget);
        }//end output_maneuver_and_target_spec()
    }//end namespace phases
}//end namespace EMTG