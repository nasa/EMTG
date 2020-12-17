
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

#include "EphemerisPeggedLaunchDirectInsertion.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        EphemerisPeggedLaunchDirectInsertion::EphemerisPeggedLaunchDirectInsertion(const std::string& name,
                                                                                            const size_t& journeyIndex,
                                                                                            const size_t& phaseIndex,
                                                                                            size_t& stageIndex,
                                                                                            Astrodynamics::universe* Universe,
                                                                                            HardwareModels::Spacecraft* mySpacecraft,
                                                                                            HardwareModels::LaunchVehicle* myLaunchVehicle,
                                                                                            missionoptions* myOptions,
                                                                                            ArrivalEvent* PreviousPhaseArrivalEvent) :
            EphemerisPeggedDeparture::EphemerisPeggedDeparture(name,
                                                               journeyIndex,
                                                               phaseIndex,
                                                               stageIndex,
                                                               Universe,
                                                               mySpacecraft,
                                                               myOptions,
                                                               PreviousPhaseArrivalEvent)
        {
            this->myLaunchVehicle = myLaunchVehicle;
            this->hasManeuver = true;
            
            if (this->isFirstEventInMission)
            {
                if (this->myLaunchVehicle->getName() == "Fixed_Initial_Mass")
                {
                    this->hasFixedInitialMass = !this->myOptions->allow_initial_mass_to_vary;
                    this->useLV = false;
                    this->hasBipropManeuver = false;
                    this->encodeInitialMass = false;
                }
                else if (this->myLaunchVehicle->getName() == "Depart_With_Spacecraft_Thrusters")
                {
                    this->hasFixedInitialMass = !this->myOptions->allow_initial_mass_to_vary;
                    this->useLV = false;
                    this->hasBipropManeuver = true;
                    this->encodeInitialMass = true;
                }
                else //any LV model
                {
                    this->hasFixedInitialMass = false;
                    this->useLV = true;
                    this->hasBipropManeuver = false;
                    this->encodeInitialMass = false;
                }
            }
            else //successive events, always a biprop maneuver and the mass is constrained to match the previous events
            {
                this->hasBipropManeuver = true;
                this->hasFixedInitialMass = false;
                this->useLV = false;
                this->encodeInitialMass = true;
                this->hasWaitTime = true;
            }

            this->InitialMassMultiplier = 1.0;
        }

        //******************************************calcbounds methods

        //calcbounds
        void EphemerisPeggedLaunchDirectInsertion::calcbounds(std::vector<size_t> timeVariables)
        {
            //we have to have this line because LaunchOrDirectInsertion does not call the base calcbounds_event_left_side()
            this->X_index_of_first_decision_variable_in_this_event = this->Xdescriptions->size();

            this->calcbounds_event_left_side(timeVariables);

            //this->calcbounds_event_main(); nothing really happens here

            this->EphemerisPeggedDeparture::calcbounds_event_right_side();

            if (this->hasBipropManeuver)
                this->calcbounds_virtual_propellant_constraints();

            //delta-v if applicable
            if (this->myOptions->objective_type == ObjectiveFunctionType::MINIMIZE_DELTAV &&
                !(this->isFirstEventInMission && !this->myOptions->include_initial_impulse_in_cost))
            {
                this->calcbounds_deltav_contribution();
            }

            this->calcbounds_specialized_constraints();
        }//end calcbounds()

        void EphemerisPeggedLaunchDirectInsertion::calcbounds_event_left_side(std::vector<size_t> timeVariables)
        {
            //Step 1: encode epoch/wait time and mass
            //if this is the first phase of the mission, encode an epoch
            if (this->isFirstEventInMission)
            {
                this->Xlowerbounds->push_back(this->myOptions->launch_window_open_date + this->myOptions->Journeys.front().wait_time_bounds[0]);
                this->Xupperbounds->push_back(fmin(this->myOptions->launch_window_open_date + this->myOptions->Journeys.front().wait_time_bounds[1], this->myBody->getEphemerisWindowClose()));
                this->X_scale_factors->push_back(1.0 / this->myUniverse->continuity_constraint_scale_factors(7));
                this->Xdescriptions->push_back(prefix + "event left state epoch");
                timeVariables.insert(timeVariables.begin(), this->Xdescriptions->size() - 1);

                if (this->encodeInitialMass)
                {
                    //mass
                    std::vector<double> MassBounds(2);
                    if (this->hasFixedInitialMass)
                    {
                        MassBounds[0] = this->myJourneyOptions->maximum_mass - 1.0e-13;
                        MassBounds[1] = this->myJourneyOptions->maximum_mass;
                    }
                    else
                    {
                        MassBounds[0] = 1.0e-13;
                        MassBounds[1] = this->myJourneyOptions->maximum_mass;
                    }

                    this->Xlowerbounds->push_back(MassBounds[0]);
                    this->Xupperbounds->push_back(MassBounds[1]);
                    this->X_scale_factors->push_back(1.0 / this->myUniverse->continuity_constraint_scale_factors(6));
                    this->Xdescriptions->push_back(prefix + "event left state mass");
                    size_t Xindex_mass = this->Xdescriptions->size() - 1;
                    this->Derivatives_of_StateBeforeEvent.push_back(std::make_tuple(Xindex_mass, 6, 1.0));
                    this->dIndex_mass_wrt_encodedMass = this->Derivatives_of_StateBeforeEvent.size() - 1;
                }
            }
            //otherwise there is a wait time and we need to encode a mass and call DepartureEvent::calcbounds_event_left_side()
            else
            {
                //base class (wait time)
                this->DepartureEvent::calcbounds_event_left_side();

                //mass
                std::vector<double> MassBounds(2);
                if (this->hasFixedInitialMass)
                {
                    MassBounds[0] = this->myJourneyOptions->maximum_mass - 1.0e-13;
                    MassBounds[1] = this->myJourneyOptions->maximum_mass;
                }
                else
                {
                    MassBounds[0] = 1.0e-13;
                    MassBounds[1] = this->myJourneyOptions->maximum_mass;
                }

                this->Xlowerbounds->push_back(MassBounds[0]);
                this->Xupperbounds->push_back(MassBounds[1]);
                this->X_scale_factors->push_back(1.0 / this->myUniverse->continuity_constraint_scale_factors(6));
                this->Xdescriptions->push_back(prefix + "event left state mass");
                size_t Xindex_mass = this->Xdescriptions->size() - 1;
                this->Derivatives_of_StateBeforeEvent.push_back(std::make_tuple(Xindex_mass, 6, 1.0));
                this->dIndex_mass_wrt_encodedMass = this->Derivatives_of_StateBeforeEvent.size() - 1;

                //also, there was a previous event, so there has to be a mass continuity constraint
                this->calcbounds_left_mass_continuity_constraint();
            }

            //Step 2: call calcbounds_left_epoch() and set up time derivatives
            this->calculate_dependencies_left_epoch(timeVariables);

            //all state variables except mass in an EphemerisPeggedBoundary event have a derivative with respect to epoch
            //we'll put in a dummy derivative of 0.0 for now, and later, when the event is processed, we'll do it right
            for (size_t stateIndex = 0; stateIndex < 6; ++stateIndex)
            {
                for (size_t Xepoch_index = 0; Xepoch_index < this->Xindices_EventLeftEpoch.size(); ++Xepoch_index)
                {
                    this->Derivatives_of_StateBeforeEvent_wrt_Time.push_back(std::make_tuple(this->Xindices_EventLeftEpoch[Xepoch_index], stateIndex, 0.0));
                }
            }

            //Step 3: initial velocity magnitude
            {
                this->Xlowerbounds->push_back(this->myJourneyOptions->initial_impulse_bounds[0]);
                this->Xupperbounds->push_back(this->myJourneyOptions->initial_impulse_bounds[1]);

                //make sure that if we are launching on a rocket that the maximum C3 does not exceed that of the rocket
                if (this->useLV)
                {
                    double C3max = this->myLaunchVehicle->getC3max();

                    if (this->Xupperbounds->back() > sqrt(C3max))
                        this->Xupperbounds->back() = sqrt(C3max);

                    double C3min = this->myLaunchVehicle->getC3min();

                    if (this->Xlowerbounds->back() < sqrt(C3min))
                        this->Xlowerbounds->back() = sqrt(C3min);
                }
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

            //Step 4: RA
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
            
            //Step 5: DEC
            {
                if (this->journeyIndex == 0) //if this is the first journey
                {
                    this->Xlowerbounds->push_back(this->myOptions->DLA_bounds[0] * math::deg2rad);
                    this->Xupperbounds->push_back(this->myOptions->DLA_bounds[1] * math::deg2rad);
                }
                else
                {
                    this->Xlowerbounds->push_back(-math::PI / 2.0);
                    this->Xupperbounds->push_back(math::PI / 2.0);
                }
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

            //Step 6: if desired, initial mass multiplier
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
        }

        void EphemerisPeggedLaunchDirectInsertion::calcbounds_virtual_propellant_constraints()
        {
            //call the base class propellant tank calcbounds
            BoundaryEventBase::calcbounds_virtual_propellant_constraints();

            //specialty derivative entries for this phase type - burn with departure stage
            //derivatives of virtual chemical fuel with respect to initial mass and v-infinity, in that order
            this->create_sparsity_entry(this->Findex_VirtualChemicalFuelConstraint,
                this->X_index_of_first_decision_variable_in_this_event,
                true,
                "mass",
                this->Gindices_dVirtualChemicalFuel_dLeftMass);

            this->create_sparsity_entry(this->Findex_VirtualChemicalFuelConstraint,
                this->X_index_of_first_decision_variable_in_this_event,
                true,
                "magnitude of outgoing velocity asymptote",
                this->Gindices_dVirtualChemicalFuel_dVinfinity);

            //derivatives of virtual chemical oxidizer with respect to initial mass and v-infinity, in that order
            this->create_sparsity_entry(this->Findex_VirtualChemicalOxidizerConstraint,
                this->X_index_of_first_decision_variable_in_this_event,
                true,
                "mass",
                this->Gindices_dVirtualChemicalOxidizer_dLeftMass);

            this->create_sparsity_entry(this->Findex_VirtualChemicalOxidizerConstraint,
                this->X_index_of_first_decision_variable_in_this_event,
                true,
                "magnitude of outgoing velocity asymptote",
                this->Gindices_dVirtualChemicalOxidizer_dVinfinity);
        }//end calcbounds_virtual_propellant_constraint

        void EphemerisPeggedLaunchDirectInsertion::calcbounds_deltav_contribution()
        {
            //add derivative with respect to launch v-infinity
            this->create_sparsity_entry(0,
                this->X_index_of_first_decision_variable_in_this_event,
                true,
                "magnitude of outgoing velocity asymptote",
                this->Gindex_dDeltav_dVinfinity);
        }//end calcbounds_virtual_deltav_constraint

        //******************************************process methods

        void EphemerisPeggedLaunchDirectInsertion::process_event(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            this->process_staging();

            this->process_event_left_side(X, Xindex, F, Findex, G, needG);

            //this->process_event_main(X, Xindex, F, Findex, G, needG); nothing happens here

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

        void EphemerisPeggedLaunchDirectInsertion::
            process_event_left_side(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG)
        {
            //Step 1: epoch/wait time
            if (!this->isFirstEventInMission)
            {
                this->DepartureEvent::process_event_left_side(X, Xindex, F, Findex, G, needG);

            }

            this->BoundaryEventBase::process_left_epoch(X, Xindex, F, Findex, G, needG);
            this->boundary_state(7) = this->state_before_event(7);

            //Step 2: mass
            if (this->encodeInitialMass)
            {
                this->state_before_event(6) = X[Xindex++];
            }
            else
            {
                this->state_before_event(6) = this->myJourneyOptions->maximum_mass;
            }
            this->boundary_state(6) = this->state_before_event(6);

            if (!this->isFirstEventInMission)
            {
                this->process_left_mass_continuity_constraint(X, Xindex, F, Findex, G, needG);
            }

            //Step 2: ephemeris and derivatives
            {
                doubleType body_state[12];
                this->myBody->locate_body(this->EventLeftEpoch,
                    body_state,
                    needG,
                    *this->myOptions);

                for (size_t stateIndex = 0; stateIndex < 6; ++stateIndex)
                {
                    this->state_before_event(stateIndex) = body_state[stateIndex];
                    this->boundary_state(stateIndex) = body_state[stateIndex];
                }

                if (needG)
                {
                    size_t nTimeVariables = this->Xindices_EventLeftEpoch.size();

                    if (this->isFirstEventInMission)
                    {
                        for (size_t stateIndex = 0; stateIndex < 6; ++stateIndex)
                        {
                            for (size_t Xepoch_index = 1; Xepoch_index < nTimeVariables + 1; ++Xepoch_index)
                            {
                                std::get<2>(this->Derivatives_of_StateBeforeEvent_wrt_Time[stateIndex * nTimeVariables + Xepoch_index + nTimeVariables - 1]) = body_state[stateIndex + 6]_GETVALUE;
                            }
                        }
                    }
                    else
                    {
                        for (size_t stateIndex = 0; stateIndex < 6; ++stateIndex)
                        {
                            for (size_t Xepoch_index = 0; Xepoch_index < nTimeVariables + 0; ++Xepoch_index)
                            {
                                std::get<2>(this->Derivatives_of_StateBeforeEvent_wrt_Time[nTimeVariables + stateIndex * nTimeVariables + Xepoch_index - 1]) = body_state[stateIndex + 6]_GETVALUE;
                            }
                        }
                    }
                }
            }//end ephemeris and derivatives

            //Step 3: initial velocity magnitude
            doubleType vinf_out = X[Xindex++];
            this->C3 = vinf_out * vinf_out;

            //Step 4: RA
            this->RLA = X[Xindex++];

            //Step 5: DEC
            this->DLA = X[Xindex++];

            //Step 6: if desired, initial mass multiplier
            if (this->isFirstEventInMission && this->myOptions->allow_initial_mass_to_vary)
            {
                this->InitialMassMultiplier = X[Xindex++];
            }

            //Step 7: construct the new velocity vector
            doubleType cosDLA = cos(DLA);
            doubleType sinDLA = sin(DLA);
            doubleType cosRLA = cos(RLA);
            doubleType sinRLA = sin(RLA);
            this->state_before_event(3) += vinf_out * cosRLA * cosDLA;
            this->state_before_event(4) += vinf_out * sinRLA * cosDLA;
            this->state_before_event(5) += vinf_out * sinDLA;

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

            //Step 8: construct the new mass
            if (this->useLV)
            {
                this->myLaunchVehicle->computePerformance(this->C3, this->myOptions->LV_margin);
                doubleType launch_mass = this->myLaunchVehicle->getDeliveredMass();

                if (launch_mass > this->myJourneyOptions->maximum_mass)
                {
                    this->MassBeforeMultiplier = this->myJourneyOptions->maximum_mass;
                    this->dm_dVinfinityOut = 0.0;
                }
                else
                {
                    this->MassBeforeMultiplier = launch_mass;
                    this->dm_dVinfinityOut = (this->myLaunchVehicle->getdmdC3() * (2.0 * vinf_out)
                        * this->InitialMassMultiplier)_GETVALUE;
                }

                this->state_before_event(6) = this->MassBeforeMultiplier * this->InitialMassMultiplier;
                std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_C3_mass]) = this->dm_dVinfinityOut;

                if (this->isFirstEventInMission && this->myOptions->allow_initial_mass_to_vary)
                    std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_InitialMassMultiplier_mass]) = this->MassBeforeMultiplier _GETVALUE;

                //if there is a fixed initial mass increment, apply it. Note you can't add mass this way, so it won't work unless it's positive
                if (this->myJourneyOptions->fixed_starting_mass_increment < 0.0)
                {
                    this->state_before_event(6) += this->myJourneyOptions->fixed_starting_mass_increment;
                }
                else if (this->myJourneyOptions->fixed_starting_mass_increment > math::SMALL)
                {
                    throw std::invalid_argument(this->prefix + " fixed starting mass increment is a positive number. You can't add mass right after launch, so it needs to be a negative number.");
                }
            }
            else
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

                    std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_mass_wrt_encodedMass]) = (this->state_before_event(6)
                        / (this->state_before_event(6) + this->chemical_oxidizer_used + this->chemical_fuel_used)) _GETVALUE;
                }
                else
                {
                    std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_C3_mass]) = 0.0;
                }

                if (this->isFirstEventInMission && this->myOptions->allow_initial_mass_to_vary)
                    std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_InitialMassMultiplier_mass]) = (this->MassBeforeMultiplier 
                        * this->state_before_event(6) / (this->state_before_event(6) + this->chemical_fuel_used + this->chemical_oxidizer_used) )_GETVALUE;
            }
        }//process_event_left_side()

        void EphemerisPeggedLaunchDirectInsertion::process_virtual_propellant_constraints(const std::vector<doubleType>& X,
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
                Gindex = this->Gindices_dVirtualChemicalFuel_dLeftMass;
                Xindex = this->jGvar->operator[](Gindex);
                if (this->isFirstEventInMission && this->myOptions->allow_initial_mass_to_vary)
                    Gentry = (this->chemical_fuel_used * this->MassBeforeMultiplier / (this->state_before_event(6) + this->chemical_fuel_used + this->chemical_oxidizer_used)) _GETVALUE;
                else
                    Gentry = (this->chemical_fuel_used / (this->state_before_event(6) + this->chemical_oxidizer_used + this->chemical_fuel_used)) _GETVALUE;
                G[Gindex] = this->X_scale_factors->operator[](Xindex) * -Gentry
                    * this->myUniverse->continuity_constraint_scale_factors(6);

                Gindex = this->Gindices_dVirtualChemicalOxidizer_dLeftMass;
                Xindex = this->jGvar->operator[](Gindex);
                if (this->isFirstEventInMission && this->myOptions->allow_initial_mass_to_vary)
                    Gentry = (this->chemical_oxidizer_used * this->MassBeforeMultiplier / (this->state_before_event(6) + this->chemical_fuel_used + this->chemical_oxidizer_used)) _GETVALUE;
                else
                    Gentry = (this->chemical_oxidizer_used / (this->state_before_event(6) + this->chemical_oxidizer_used + this->chemical_fuel_used)) _GETVALUE;
                G[Gindex] = this->X_scale_factors->operator[](Xindex) * -Gentry
                    * this->myUniverse->continuity_constraint_scale_factors(6);

                //derivatives with respect to v-infinity
                Gindex = this->Gindices_dVirtualChemicalFuel_dVinfinity;
                Xindex = this->jGvar->operator[](Gindex);
                G[Gindex] = this->X_scale_factors->operator[](Xindex) * this->dChemicalFuel_dVinfinity
                    * this->myUniverse->continuity_constraint_scale_factors(6);

                Gindex = this->Gindices_dVirtualChemicalOxidizer_dVinfinity;
                Xindex = this->jGvar->operator[](Gindex);
                G[Gindex] = this->X_scale_factors->operator[](Xindex) * this->dChemicalOxidizer_dVinfinity
                    * this->myUniverse->continuity_constraint_scale_factors(6);
            }
        }//end process_virtual_propellant_constraints

        void EphemerisPeggedLaunchDirectInsertion::process_deltav_contribution(const std::vector<doubleType>& X,
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
        void EphemerisPeggedLaunchDirectInsertion::output(std::ofstream& outputfile,
            const double& launchdate,
            size_t& eventcount)
        {
            this->mySpacecraft->setActiveStage(this->stageIndex);

            //base class
            EphemerisPeggedDeparture::output(outputfile, launchdate, eventcount);

            std::string event_type;
            if (!this->hasBipropManeuver)
                event_type = "launch";
            else
                event_type = "departure";

            std::string boundary_name = this->myBody->name;

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
                this->hasFixedInitialMass || this->useLV ? Isp : -1.0,
                this->mySpacecraft->getAvailablePower(),
                0.0,
                0,
                0.0,
                "none");
        }//end output()

        void EphemerisPeggedLaunchDirectInsertion::output_maneuver_and_target_spec(std::ofstream& maneuver_spec_file, std::ofstream& target_spec_file, bool& haveManeuverNeedTarget)
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