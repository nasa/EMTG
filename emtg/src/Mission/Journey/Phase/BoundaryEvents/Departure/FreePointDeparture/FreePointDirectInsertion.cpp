
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

#include "FreePointDirectInsertion.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        FreePointDirectInsertion::FreePointDirectInsertion(const std::string& name,
                                                           const size_t& journeyIndex,
                                                           const size_t& phaseIndex,
                                                           size_t& stageIndex,
                                                           Astrodynamics::universe* Universe,
                                                           HardwareModels::Spacecraft* mySpacecraft,
                                                           missionoptions* myOptions,
                                                           ArrivalEvent* PreviousPhaseArrivalEvent) :
            FreePointDeparture::FreePointDeparture(name,
                                                   journeyIndex,
                                                   phaseIndex,
                                                   stageIndex,
                                                   Universe,
                                                   mySpacecraft,
                                                   myOptions,
                                                   PreviousPhaseArrivalEvent)
        {
            this->hasBipropManeuver = true;

            this->hasManeuver = true;

            this->InitialMassMultiplier = 1.0;

            this->departureImpulseDirection.resize(3, 1, 0.0);

        }

        //******************************************calcbounds methods

        //calcbounds
        void FreePointDirectInsertion::calcbounds(std::vector<size_t> timeVariables)
        {
            this->calcbounds_event_left_side(timeVariables);

            this->calcbounds_event_main();

            this->calcbounds_event_right_side();

            this->calcbounds_virtual_propellant_constraints();

            //delta-v if applicable
            if (this->myOptions->objective_type == ObjectiveFunctionType::MINIMIZE_DELTAV &&
                !(this->isFirstEventInMission && !this->myOptions->include_initial_impulse_in_cost))
            {
                this->calcbounds_deltav_contribution();
            }

            this->calcbounds_specialized_constraints();
        }//end calcbounds()

        void FreePointDirectInsertion::calcbounds_event_main()
        {
            //Step 1: initial velocity magnitude
            {
                this->Xlowerbounds->push_back(this->myJourneyOptions->initial_impulse_bounds[0]);
                this->Xupperbounds->push_back(this->myJourneyOptions->initial_impulse_bounds[1]);
                this->X_scale_factors->push_back(this->myUniverse->LU / this->myUniverse->TU);
                this->Xdescriptions->push_back(prefix + "magnitude of outgoing velocity asymptote");
                this->Xindex_C3 = this->Xdescriptions->size() - 1;

                //C3 affects all three velocity components and mass
                this->Derivatives_of_StateBeforeEvent.push_back({ this->Xindex_C3, 3, 1.0 });
                this->dIndex_C3_vx = this->Derivatives_of_StateBeforeEvent.size() - 1;
                this->Derivatives_of_StateBeforeEvent.push_back({ this->Xindex_C3, 4, 1.0 });
                this->dIndex_C3_vy = this->Derivatives_of_StateBeforeEvent.size() - 1;
                this->Derivatives_of_StateBeforeEvent.push_back({ this->Xindex_C3, 5, 1.0 });
                this->dIndex_C3_vz = this->Derivatives_of_StateBeforeEvent.size() - 1;
                this->Derivatives_of_StateBeforeEvent.push_back({ this->Xindex_C3, 6, 1.0 });
                this->dIndex_C3_mass = this->Derivatives_of_StateBeforeEvent.size() - 1;
            }

            //Step 2: RA
            {
                this->Xlowerbounds->push_back(this->myOptions->RLA_bounds[0] * math::deg2rad);
                this->Xupperbounds->push_back(this->myOptions->RLA_bounds[1] * math::deg2rad);
                this->X_scale_factors->push_back(math::TwoPI);
                this->Xdescriptions->push_back(prefix + "RA of departure asymptote");
                this->Xindex_RLA = this->Xdescriptions->size() - 1;
                this->Derivatives_of_StateBeforeEvent.push_back({ this->Xindex_RLA, 3, 1.0 });
                this->dIndex_RLA_vx = this->Derivatives_of_StateBeforeEvent.size() - 1;
                this->Derivatives_of_StateBeforeEvent.push_back({ this->Xindex_RLA, 4, 1.0 });
                this->dIndex_RLA_vy = this->Derivatives_of_StateBeforeEvent.size() - 1;
            }
            
            //Step 3: DEC
            {
                this->Xlowerbounds->push_back(-math::PI / 2.0);
                this->Xupperbounds->push_back(math::PI / 2.0);
                this->X_scale_factors->push_back(math::TwoPI);
                this->Xdescriptions->push_back(prefix + "DEC of departure asymptote");
                this->Xindex_DLA = this->Xdescriptions->size() - 1;
                this->Derivatives_of_StateBeforeEvent.push_back({ this->Xindex_DLA, 3, 1.0 });
                this->dIndex_DLA_vx = this->Derivatives_of_StateBeforeEvent.size() - 1;
                this->Derivatives_of_StateBeforeEvent.push_back({ this->Xindex_DLA, 4, 1.0 });
                this->dIndex_DLA_vy = this->Derivatives_of_StateBeforeEvent.size() - 1;
                this->Derivatives_of_StateBeforeEvent.push_back({ this->Xindex_DLA, 5, 1.0 });
                this->dIndex_DLA_vz = this->Derivatives_of_StateBeforeEvent.size() - 1;
            }

            //Step 4: if desired, initial mass multiplier
            if (this->isFirstEventInMission && this->myOptions->allow_initial_mass_to_vary)
            {
                this->Xlowerbounds->push_back(1.0e-13);
                this->Xupperbounds->push_back(1.0);
                this->X_scale_factors->push_back(1.0);
                this->Xdescriptions->push_back(this->prefix + "initial mass multiplier");
                this->Xindex_InitialMassMultiplier = this->Xdescriptions->size() - 1;
                this->Derivatives_of_StateBeforeEvent.push_back({ this->Xindex_InitialMassMultiplier, 6, 1.0 });
                this->dIndex_InitialMassMultiplier_mass = this->Derivatives_of_StateBeforeEvent.size() - 1;
            }

            //Step 5: if desired, constrain the maneuver to lie along the velocity vector
            if (this->myJourneyOptions->force_free_point_direct_insertion_along_velocity_vector)
            {
                for (size_t maneuver_velocity_index : {0, 1, 2})
                {
                    //Step 5.1: create the constraint
                    this->Flowerbounds->push_back(0.0);
                    this->Fupperbounds->push_back(0.0);
                    this->Fdescriptions->push_back(this->prefix + "constraint maneuver to lie along velocity vector " + this->stateNames[maneuver_velocity_index + 3]);

                    //Step 5.2: derivative wrt RLA
                    this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                        this->Xindex_RLA,
                        this->Gindex_maneuver_direction_constraint_wrt_RLA);

                    //Step 5.3: derivative wrt DLA
                    this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                        this->Xindex_DLA,
                        this->Gindex_maneuver_direction_constraint_wrt_DLA);

                    //Step 5.4: derivatives wrt boundary variables
                    std::vector<size_t> state_Gindex_maneuver_direction_constraint_wrt_free_point_velocity;
                    std::vector<size_t> state_dIndex_velocity_before_maneuver;

                    //non-time variables
                    for (size_t dIndex = 0; dIndex < this->Derivatives_of_StateBeforeEvent.size(); ++dIndex)
                    {
                        for (size_t boundary_velocity_index : {0, 1, 2})
                        {
                            //does this derivative entry affect the velocity entry of interest?
                            if (std::get<1>(this->Derivatives_of_StateBeforeEvent[dIndex]) == boundary_velocity_index + 3)
                            {
                                size_t Xindex = std::get<0>(this->Derivatives_of_StateBeforeEvent[dIndex]);

                                //we don't want to pick up the variable that control the departure impulse
                                if (Xindex != this->Xindex_C3
                                    && Xindex != this->Xindex_RLA
                                    && Xindex != this->Xindex_DLA)
                                {
                                    //hooray! we this entry is important!

                                    //save the dIndex
                                    state_dIndex_velocity_before_maneuver.push_back(dIndex);

                                    //create a sparsity pattern entry
                                    this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                                        Xindex,
                                        state_Gindex_maneuver_direction_constraint_wrt_free_point_velocity);
                                }//end code block to add a sparsity entry
                            }//end code block to execute if this is a relevant derivative entry
                        }//end loop over boundary velocity indices
                    }//end loop over non-time derivative entries

                    this->Gindex_maneuver_direction_constraint_wrt_free_point_velocity.push_back(state_Gindex_maneuver_direction_constraint_wrt_free_point_velocity);
                    this->dIndex_velocity_before_maneuver.push_back(state_dIndex_velocity_before_maneuver);

                    //time variables
                    std::vector<size_t> state_Gindex_maneuver_direction_constraint_wrt_free_point_velocity_time_variables;
                    std::vector<size_t> state_dIndex_velocity_before_maneuver_time_variables;

                    for (size_t dIndex = 0; dIndex < this->Derivatives_of_StateBeforeEvent_wrt_Time.size(); ++dIndex)
                    {
                        for (size_t boundary_velocity_index : {0, 1, 2})
                        {
                            //does this derivative entry affect the velocity entry of interest?
                            if (std::get<1>(this->Derivatives_of_StateBeforeEvent_wrt_Time[dIndex]) == maneuver_velocity_index + 3)
                            {
                                size_t Xindex = std::get<0>(this->Derivatives_of_StateBeforeEvent_wrt_Time[dIndex]);

                                //save the dIndex
                                state_dIndex_velocity_before_maneuver_time_variables.push_back(dIndex);

                                //create a sparsity pattern entry
                                this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                                    Xindex,
                                    state_Gindex_maneuver_direction_constraint_wrt_free_point_velocity_time_variables);
                            }//end code block to execute if this is a relevant derivative entry
                        }//end loop over boundary velocity entries
                    }//end loop over non-time derivative entries

                    this->Gindex_maneuver_direction_constraint_wrt_free_point_velocity_time_variables.push_back(state_Gindex_maneuver_direction_constraint_wrt_free_point_velocity_time_variables);
                    this->dIndex_velocity_before_maneuver_time_variables.push_back(state_dIndex_velocity_before_maneuver_time_variables);
                }//end loop over maneuver velocity index
            }//end force_free_point_direct_insertion_along_velocity_vector
        }//end calcbounds_event_main()

        void FreePointDirectInsertion::calcbounds_virtual_propellant_constraints()
        {
            //call the base class propellant tank calcbounds
            BoundaryEventBase::calcbounds_virtual_propellant_constraints();

            //specialty derivative entries for this phase type - burn with departure stage
            //derivatives of virtual chemical fuel with respect to initial mass and v-infinity, in that order
            this->create_sparsity_entry(this->Findex_VirtualChemicalFuelConstraint,
                this->X_index_of_first_decision_variable_in_this_event,
                true,
                "state mass",
                this->Gindex_dVirtualChemicalFuel_dLeftMass);

            this->create_sparsity_entry(this->Findex_VirtualChemicalFuelConstraint,
                this->X_index_of_first_decision_variable_in_this_event,
                true,
                "magnitude of outgoing velocity asymptote",
                this->Gindex_dVirtualChemicalFuel_dVinfinity);

            //derivatives of virtual chemical oxidizer with respect to initial mass and v-infinity, in that order
            this->create_sparsity_entry(this->Findex_VirtualChemicalOxidizerConstraint,
                this->X_index_of_first_decision_variable_in_this_event,
                true,
                "state mass",
                this->Gindex_dVirtualChemicalOxidizer_dLeftMass);

            this->create_sparsity_entry(this->Findex_VirtualChemicalOxidizerConstraint,
                this->X_index_of_first_decision_variable_in_this_event,
                true,
                "magnitude of outgoing velocity asymptote",
                this->Gindex_dVirtualChemicalOxidizer_dVinfinity);

            if (this->isFirstEventInMission && this->myOptions->allow_initial_mass_to_vary)
            {
                this->create_sparsity_entry(this->Findex_VirtualChemicalFuelConstraint,
                    this->X_index_of_first_decision_variable_in_this_event,
                    true,
                    "initial mass multiplier",
                    this->Gindex_dVirtualChemicalFuel_dInitialMassMultiplier);

                this->create_sparsity_entry(this->Findex_VirtualChemicalOxidizerConstraint,
                    this->X_index_of_first_decision_variable_in_this_event,
                    true,
                    "initial mass multiplier",
                    this->Gindex_dVirtualChemicalOxidizer_dInitialMassMultiplier);

            }
            
            //derivatives of virtual tanks with respect to time, because mass has a derivative with respect to time
            //WHY???
            /*if (this->myPropagatorType == PropagatorType::IntegratedPropagator)
            {
                this->create_sparsity_entry(this->Findex_VirtualChemicalFuelConstraint,
                    this->X_index_of_first_decision_variable_in_this_event,
                    true,
                    "event left state epoch",
                    this->Gindex_dVirtualChemicalFuel_dPropagationTime);

                this->create_sparsity_entry(this->Findex_VirtualChemicalOxidizerConstraint,
                    this->X_index_of_first_decision_variable_in_this_event,
                    true,
                    "event left state epoch",
                    this->Gindex_dVirtualChemicalOxidizer_dPropagationTime);
            }*/
        }//end calcbounds_virtual_propellant_constraint

        void FreePointDirectInsertion::calcbounds_deltav_contribution()
        {
            //add derivative with respect to launch v-infinity
            this->create_sparsity_entry(0,
                this->X_index_of_first_decision_variable_in_this_event,
                true,
                "magnitude of outgoing velocity asymptote",
                this->Gindex_dDeltav_dVinfinity);
        }//end calcbounds_virtual_deltav_constraint

        //******************************************process methods

        void FreePointDirectInsertion::process_event(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            this->process_staging();

            this->process_event_left_side(X, Xindex, F, Findex, G, needG);

            this->process_event_main(X, Xindex, F, Findex, G, needG);

            this->process_event_right_side(X, Xindex, F, Findex, G, needG);
            
            if (this->hasBipropManeuver)
                this->process_virtual_propellant_constraints(X, Xindex, F, Findex, G, needG);

            //virtual delta-v variable and constraint exist only if this event has deterministic maneuvers
            if (!(this->isFirstEventInMission && !this->myOptions->include_initial_impulse_in_cost))
            {
                this->process_deltav_contribution(X, Xindex, F, Findex, G, needG);
            }

            BoundaryEventBase::process_specialized_constraints(X, Xindex, F, Findex, G, needG);
        }//end process_event

        void FreePointDirectInsertion::
            process_event_main(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG)
        {
            //Step 1: initial velocity magnitude
            doubleType vinf_out = X[Xindex++];
            this->C3 = vinf_out * vinf_out;

            //Step 2: RA
            this->RLA = X[Xindex++];

            //Step 3: DEC
            this->DLA = X[Xindex++];

            //Step 4: if desired, initial mass multiplier
            if (this->isFirstEventInMission && this->myOptions->allow_initial_mass_to_vary)
            {
                this->InitialMassMultiplier = X[Xindex++];
            }

            //Step 5: form the departure impulse direction vector
            doubleType cosRLA = cos(RLA);
            doubleType sinRLA = sin(RLA);
            doubleType cosDLA = cos(DLA);
            doubleType sinDLA = sin(DLA);
            this->departureImpulseDirection(0) = cosRLA * cosDLA;
            this->departureImpulseDirection(1) = sinRLA * cosDLA;
            this->departureImpulseDirection(2) = sinDLA;
            
            //Step 6: if desired, constrain the delta-v vector to be along the velocity vector
            if (this->myJourneyOptions->force_free_point_direct_insertion_along_velocity_vector)
            {
                math::Matrix<doubleType> velocity_before_maneuver = this->state_before_event.getSubMatrix1D(3, 5);
                math::Matrix<doubleType> unit_velocity_before_maneuver = velocity_before_maneuver.unitize();


                for (size_t maneuver_velocity_index : {0, 1, 2})
                {
                    //Step 6.1: evaluate the constraint
                    F[Findex++] = this->departureImpulseDirection(maneuver_velocity_index) - unit_velocity_before_maneuver(maneuver_velocity_index);

                    //Step 6.2: derivatives
                    if (needG)
                    {
                        math::Matrix<doubleType> dConstraint_dBoundaryVelocity = unit_velocity_before_maneuver.unitDerivative(velocity_before_maneuver);

                        //Step 6.2: derivative wrt RLA
                        {
                            size_t Gindex = this->Gindex_maneuver_direction_constraint_wrt_RLA[maneuver_velocity_index];
                            switch (maneuver_velocity_index)
                            {
                                case 0: //x direction
                                {
                                    G[Gindex] = this->X_scale_factors->operator[](this->Xindex_RLA) * (-sinRLA * cosDLA)_GETVALUE;
                                    break;
                                }
                                case 1: //y direction
                                {
                                    G[Gindex] = this->X_scale_factors->operator[](this->Xindex_RLA) * (cosRLA * cosDLA)_GETVALUE;
                                    break;
                                }
                                case 2: //z direction
                                {
                                    G[Gindex] = 0.0;
                                    break;
                                }
                                default:
                                {
                                    //by definition this cannot happen. Default case included for completeness
                                    break;
                                }
                            }
                        }//end derivatives wrt RLA

                        //Step 6.3: derivative wrt DLA
                        {
                            size_t Gindex = this->Gindex_maneuver_direction_constraint_wrt_DLA[maneuver_velocity_index];
                            switch (maneuver_velocity_index)
                            {
                                case 0: //x direction
                                {
                                    G[Gindex] = this->X_scale_factors->operator[](this->Xindex_DLA) * (cosRLA * -sinDLA)_GETVALUE;
                                    break;
                                }
                                case 1: //y direction
                                {
                                    G[Gindex] = this->X_scale_factors->operator[](this->Xindex_DLA) * (sinRLA * -sinDLA)_GETVALUE;
                                    break;
                                }
                                case 2: //z direction
                                {
                                    G[Gindex] = this->X_scale_factors->operator[](this->Xindex_DLA) * cosDLA _GETVALUE;
                                    break;
                                }
                                default:
                                {
                                    //by definition this cannot happen. Default case included for completeness
                                    break;
                                }
                            }
                        }//end derivatives wrt DLA

                        //Step 6.4: derivatives wrt boundary variables
                        {
                            //first we need to zero out all the entries

                            //non-time variables
                            for (size_t Gindex : this->Gindex_maneuver_direction_constraint_wrt_free_point_velocity[maneuver_velocity_index])
                            {
                                G[Gindex] = 0.0;
                            }
                            for (size_t entryIndex = 0; entryIndex < this->dIndex_velocity_before_maneuver[maneuver_velocity_index].size(); ++entryIndex)
                            {
                                size_t dIndex = this->dIndex_velocity_before_maneuver[maneuver_velocity_index][entryIndex];
                                size_t Gindex = this->Gindex_maneuver_direction_constraint_wrt_free_point_velocity[maneuver_velocity_index][entryIndex];
                                size_t Xindex = std::get<0>(this->Derivatives_of_StateBeforeEvent[dIndex]);
                                size_t stateIndex = std::get<1>(this->Derivatives_of_StateBeforeEvent[dIndex]);

                                double TheDerivative = -dConstraint_dBoundaryVelocity(maneuver_velocity_index, stateIndex - 3)_GETVALUE * std::get<2>(this->Derivatives_of_StateBeforeEvent[dIndex]);
                                
                                G[Gindex] += this->X_scale_factors->operator[](Xindex)
                                    * TheDerivative;
                            }//end loop over non-time variables

                            //time variables
                            for (size_t Gindex : this->Gindex_maneuver_direction_constraint_wrt_free_point_velocity_time_variables[maneuver_velocity_index])
                            {
                                G[Gindex] = 0.0;
                            }
                            for (size_t entryIndex = 0; entryIndex < this->dIndex_velocity_before_maneuver_time_variables[maneuver_velocity_index].size(); ++entryIndex)
                            {
                                size_t dIndex = this->dIndex_velocity_before_maneuver_time_variables[maneuver_velocity_index][entryIndex];
                                size_t Gindex = this->Gindex_maneuver_direction_constraint_wrt_free_point_velocity_time_variables[maneuver_velocity_index][entryIndex];
                                size_t Xindex = std::get<0>(this->Derivatives_of_StateBeforeEvent_wrt_Time[dIndex]);
                                size_t stateIndex = std::get<1>(this->Derivatives_of_StateBeforeEvent[dIndex]);

                                double TheDerivative = -dConstraint_dBoundaryVelocity(maneuver_velocity_index, stateIndex - 3)_GETVALUE * std::get<2>(this->Derivatives_of_StateBeforeEvent_wrt_Time[dIndex]);

                                G[Gindex] = this->X_scale_factors->operator[](Xindex)
                                    * TheDerivative;
                            }//end loop over time variables
                        }//end loop over boundary velocity index
                    }//end derivatives
                }//end loop over maneuver velocity index
            }//end force_free_point_direct_insertion_along_velocity_vector

            //Step 7: add the delta-v to the state vector
            this->state_before_event(3) += vinf_out * this->departureImpulseDirection(0);
            this->state_before_event(4) += vinf_out * this->departureImpulseDirection(1);
            this->state_before_event(5) += vinf_out * this->departureImpulseDirection(2);

            if (needG)
            {
                //derivatives with respect to C3
                std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_C3_vx]) = (cosRLA * cosDLA)_GETVALUE;
                std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_C3_vy]) = (sinRLA * cosDLA)_GETVALUE;
                std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_C3_vz]) = sinDLA _GETVALUE;

                //derivatives with respect to RLA
                std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_RLA_vx]) = (vinf_out * -sinRLA * cosDLA)_GETVALUE;
                std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_RLA_vy]) = (vinf_out *  cosRLA * cosDLA)_GETVALUE;

                //derivatives with respect to DLA
                std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_DLA_vx]) = (vinf_out * cosRLA * -sinDLA)_GETVALUE;
                std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_DLA_vy]) = (vinf_out * sinRLA * -sinDLA)_GETVALUE;
                std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_DLA_vz]) = (vinf_out * cosDLA)_GETVALUE;
            }

            //Step 6: construct the new mass
            {
                //mass multiplier
                this->MassBeforeMultiplier = this->state_before_event(6);
                this->state_before_event(6) = this->MassBeforeMultiplier * this->InitialMassMultiplier;

                if (this->hasBipropManeuver)
                {
                    //perform a maneuver with the current stage
                    this->mySpacecraft->computeChemicalPropulsionPerformance(vinf_out, this->state_before_event(6), true, this->ChemicalManeuverType);

                    this->chemical_fuel_used = this->mySpacecraft->getChemFuelConsumedThisManeuver();
                    this->chemical_oxidizer_used = this->mySpacecraft->getChemOxidizerConsumedThisManeuver();

                    this->dChemicalFuel_dVinfinity = this->mySpacecraft->getdFuelConsumedThisManeuver_ddeltav();
                    this->dChemicalOxidizer_dVinfinity = this->mySpacecraft->getdOxidizerConsumedThisManeuver_ddeltav();

                    this->state_before_event(6) -= this->chemical_fuel_used + this->chemical_oxidizer_used;

                    std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_C3_mass]) = this->mySpacecraft->getdSpacecraftMassThisManeuver_ddeltav();
                }

                if (this->isFirstEventInMission && this->myOptions->allow_initial_mass_to_vary)
                    std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_InitialMassMultiplier_mass]) = (this->MassBeforeMultiplier 
                        * this->state_before_event(6) / (this->state_before_event(6) + this->chemical_fuel_used + this->chemical_oxidizer_used) )_GETVALUE;

                //derivative of mass with respect to encoded mass
                std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_mass_wrt_encodedMass]) = (this->state_before_event(6) / this->MassBeforeMultiplier) _GETVALUE;

                //derivative of mass with respect to epoch (important for when free point is being propagated)
                if (this->AllowStateToPropagate)
                {
                    for (size_t dIndex : this->dIndex_StateBeforeEvent_wrt_Time.back())
                    {
                        std::get<2>(this->Derivatives_of_StateBeforeEvent_wrt_Time[dIndex]) *= this->mySpacecraft->get_dMassAfterManeuver_dMassAtManeuver();
                    }
                }
            }
        }//process_event_main()

        void FreePointDirectInsertion::process_virtual_propellant_constraints(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //call the base class propellant constraints
            BoundaryEventBase::process_virtual_propellant_constraints(X, Xindex, F, Findex, G, needG);

            //event-specific

            if (needG)
            {
                size_t Gindex, Xindex;
                double Gentry;

                //derivatives with respect to mass
                Gindex = this->Gindex_dVirtualChemicalFuel_dLeftMass;
                Xindex = this->jGvar->operator[](Gindex);
                G[Gindex] = this->X_scale_factors->operator[](Xindex)
                    * ( -this->chemical_fuel_used / this->MassBeforeMultiplier ) _GETVALUE
                    * this->myUniverse->continuity_constraint_scale_factors(6);

                Gindex = this->Gindex_dVirtualChemicalOxidizer_dLeftMass;
                Xindex = this->jGvar->operator[](Gindex);
                G[Gindex] = this->X_scale_factors->operator[](Xindex) 
                    * ( -this->chemical_oxidizer_used / this->MassBeforeMultiplier ) _GETVALUE
                    * this->myUniverse->continuity_constraint_scale_factors(6);


                //derivatives with respect to v-infinity
                Gindex = this->Gindex_dVirtualChemicalFuel_dVinfinity;
                Xindex = this->jGvar->operator[](Gindex);
                G[Gindex] = this->X_scale_factors->operator[](Xindex) * this->dChemicalFuel_dVinfinity
                    * this->myUniverse->continuity_constraint_scale_factors(6);

                Gindex = this->Gindex_dVirtualChemicalOxidizer_dVinfinity;
                Xindex = this->jGvar->operator[](Gindex);
                G[Gindex] = this->X_scale_factors->operator[](Xindex) * this->dChemicalOxidizer_dVinfinity
                    * this->myUniverse->continuity_constraint_scale_factors(6);

                //derivatives with respect to initial mass multiplier
                if (this->isFirstEventInMission && this->myOptions->allow_initial_mass_to_vary)
                {
                    Gindex = this->Gindex_dVirtualChemicalFuel_dInitialMassMultiplier;
                    Xindex = this->jGvar->operator[](Gindex);
                    
                    Gentry = (this->chemical_fuel_used * this->MassBeforeMultiplier / (this->state_before_event(6) + this->chemical_fuel_used + this->chemical_oxidizer_used)) _GETVALUE;

                    G[Gindex] = this->X_scale_factors->operator[](Xindex) * -Gentry
                        * this->myUniverse->continuity_constraint_scale_factors(6);

                    Gindex = this->Gindex_dVirtualChemicalOxidizer_dInitialMassMultiplier;
                    Xindex = this->jGvar->operator[](Gindex);
            
                    Gentry = (this->chemical_oxidizer_used * this->MassBeforeMultiplier / (this->state_before_event(6) + this->chemical_fuel_used + this->chemical_oxidizer_used)) _GETVALUE;
            
                    G[Gindex] = this->X_scale_factors->operator[](Xindex) * -Gentry
                        * this->myUniverse->continuity_constraint_scale_factors(6);
                }

                //derivatives with respect to time, because mass has derivatives with respect to time
                //Again, huh?
                /*if (this->myPropagatorType == PropagatorType::IntegratedPropagator)
                {
                    Gindex = this->Gindex_dVirtualChemicalFuel_dPropagationTime;
                    Xindex = this->jGvar->operator[](Gindex);
                    G[Gindex] = this->X_scale_factors->operator[](Xindex)
                        * -this->mySpacecraft->get_dFuelConsumedThisManeuver_dMassAtManeuver() * this->STM(6, 13)
                        * this->myUniverse->continuity_constraint_scale_factors(6);

                    Gindex = this->Gindex_dVirtualChemicalOxidizer_dPropagationTime;
                    Xindex = this->jGvar->operator[](Gindex);
                    G[Gindex] = this->X_scale_factors->operator[](Xindex)
                        * -this->mySpacecraft->get_dOxidizerConsumedThisManeuver_dMassAtManeuver() * this->STM(6, 13)
                        * this->myUniverse->continuity_constraint_scale_factors(6);
                }*/
            }
        }//end process_virtual_propellant_constraints

        void FreePointDirectInsertion::process_deltav_contribution(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //add the event delta-v
            this->EventDeterministicDeltav = sqrt(this->C3);

            //event-specific
            if (this->myOptions->objective_type == ObjectiveFunctionType::MINIMIZE_DELTAV
                && needG)
            {
                //derivatives with respect to v-infinity
                size_t Gindex = this->Gindex_dDeltav_dVinfinity;
                size_t Xindex = this->jGvar->operator[](Gindex);
                G[Gindex] = this->X_scale_factors->operator[](Xindex);
            }
        }//end process_virtual_deltav_constraint()

         //******************************************output methods
        void FreePointDirectInsertion::output(std::ofstream& outputfile,
            const double& launchdate,
            size_t& eventcount)
        {

            this->mySpacecraft->setActiveStage(this->stageIndex);

            //base class
            FreePointDeparture::output(outputfile, launchdate, eventcount);

            std::string event_type = "departure";

            std::string boundary_name = "free point";

            math::Matrix<doubleType> empty3vector(3, 1, 0.0);
            math::Matrix<doubleType> dV(3, 1, std::vector<doubleType>({sqrt(this->C3) * cos(this->RLA) * cos(this->DLA),
                                                                       sqrt(this->C3) * sin(this->RLA) * cos(this->DLA),
                                                                       sqrt(this->C3) * sin(this->DLA) }));

            //where is the Sun?

            math::Matrix<doubleType> R_sc_Sun(3, 1, 0.0);

            if (this->myUniverse->central_body_SPICE_ID == 10)
            {
                R_sc_Sun = this->state_after_event.getSubMatrix1D(0, 2);
            }
            else
            {
                //where is the central body relative to the sun?

                doubleType central_body_state_and_derivatives[12];
                this->myUniverse->locate_central_body(this->state_after_event(7),
                    central_body_state_and_derivatives,
                    *this->myOptions,
                    false);

                math::Matrix<doubleType> R_CB_Sun(3, 1, 0.0);
                for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
                {
                    R_CB_Sun(stateIndex) = central_body_state_and_derivatives[stateIndex];
                }

                R_sc_Sun = this->state_after_event.getSubMatrix1D(0, 2) + R_CB_Sun;
            }

            this->mySpacecraft->computePowerState(R_sc_Sun.getSubMatrix1D(0, 2).norm() / this->myOptions->AU, this->state_after_event(7));


            double Isp = this->ChemicalManeuverType == PropulsionSystemChoice::Biprop ? this->mySpacecraft->getBipropIsp() : this->mySpacecraft->getMonopropIsp();

            write_output_line(outputfile,
                eventcount,
                event_type,
                boundary_name,
                this->EventTimeWidth,
                0.0,
                0.0,
                0.0,
                this->RLA,
                this->DLA,
                this->C3,
                this->state_after_event,
                dV,
                empty3vector,
                sqrt(this->C3),
                0.0,
                Isp,
                this->mySpacecraft->getAvailablePower(),
                0.0,
                0,
                0.0,
                "none");
        }//end output()

        void FreePointDirectInsertion::output_maneuver_and_target_spec(std::ofstream& maneuver_spec_file, std::ofstream& target_spec_file, bool& haveManeuverNeedTarget)

        {

            //Step 1: maneuver spec

            //Step 1.1: instantiate and populate a maneuver spec object

            maneuver_spec_line myManeuverSpecLine(this->name);



            double Isp = this->ChemicalManeuverType == PropulsionSystemChoice::Biprop ? this->mySpacecraft->getBipropIsp() : this->mySpacecraft->getMonopropIsp();

            doubleType massFlowRate = this->mySpacecraft->getchemthrust() / Isp / this->myOptions->g0;

            doubleType mass_before_maneuver = this->state_before_event(6) + this->chemical_fuel_used + this->chemical_oxidizer_used;

            doubleType maneuverDuration = (mass_before_maneuver - this->state_before_event(6)) / massFlowRate;




            math::Matrix<doubleType> dV(3, 1, std::vector<doubleType>({ sqrt(this->C3) * cos(this->RLA) * cos(this->DLA),
                                                                        sqrt(this->C3) * sin(this->RLA) * cos(this->DLA),
                                                                        sqrt(this->C3) * sin(this->DLA) }));



            myManeuverSpecLine.append_maneuver_spec_item("EME2000",

                this->state_before_event(7),

                dV.unitize(),

                mass_before_maneuver,

                this->state_before_event(6),

                this->mySpacecraft->getchemthrust(),

                massFlowRate,

                maneuverDuration,

                1.0);



            //Step 1.2: write the maneuver spec

            myManeuverSpecLine.write(maneuver_spec_file);



            //Step 1.3: signal that we need a target spec

            haveManeuverNeedTarget = true;

        }//end output_maneuver_and_target_spec()
    }//end namespace BoundaryEvents
}//end namespace EMTG