
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

#include "PeriapseLaunchOrImpulsiveDeparture.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        PeriapseLaunchOrImpulsiveDeparture::PeriapseLaunchOrImpulsiveDeparture(const std::string& name,
                                                                               const size_t& journeyIndex,
                                                                               const size_t& phaseIndex,
                                                                               size_t& stageIndex,
                                                                               Astrodynamics::universe* Universe,
                                                                               HardwareModels::Spacecraft* mySpacecraft,
                                                                               HardwareModels::LaunchVehicle* myLaunchVehicle,
                                                                               missionoptions* myOptions,
                                                                               ArrivalEvent* PreviousPhaseArrivalEvent) :
            PeriapseDeparture::PeriapseDeparture(name,
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


            this->Vunit.resize(3, 1, 0.0);
            
            if (this->isFirstEventInMission)
            {
                if (this->myLaunchVehicle->getName() == "Fixed_Initial_Mass")
                {
                    this->hasFixedInitialMass = !this->myOptions->allow_initial_mass_to_vary;
                    this->useLV = false;
                    this->hasBipropManeuver = false;    
                }
                else if (this->myLaunchVehicle->getName() == "Depart_With_Spacecraft_Thrusters")
                {
                    this->hasFixedInitialMass = !this->myOptions->allow_initial_mass_to_vary;
                    this->useLV = false;
                    this->hasBipropManeuver = true;
                }
                else //any LV model
                {
                    this->hasFixedInitialMass = false;
                    this->useLV = true;
                    this->hasBipropManeuver = false;
                }
            }
            else //successive events, always a biprop maneuver and the mass is inherited from the previous events
            {
                this->hasBipropManeuver = true;
                this->hasFixedInitialMass = true;
                this->useLV = false;
            }
        }

        //******************************************calcbounds methods

        //calcbounds
        void PeriapseLaunchOrImpulsiveDeparture::calcbounds()
        {
            //we have to have this line because LaunchOrDirectInsertion does not call the base calcbounds_event_left_side()
            this->X_index_of_first_decision_variable_in_this_event = this->Xdescriptions->size();

            this->calcbounds_event_left_side();

            //this->calcbounds_event_main(); nothing really happens here

            this->PeriapseDeparture::calcbounds_event_right_side();

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

        void PeriapseLaunchOrImpulsiveDeparture::calcbounds_event_left_side()
        {
            //Step 1: encode mass bounds
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

            //Step 2: encode radius bounds
            std::vector<double> RadiusBounds({ this->myUniverse->central_body.radius + this->myJourneyOptions->PeriapseDeparture_altitude_bounds[0],
                                               this->myUniverse->central_body.radius + this->myJourneyOptions->PeriapseDeparture_altitude_bounds[1] });

            //Step 3: encode velocity bounds
            std::vector<double> VelocityMagnitudeBounds({ sqrt(this->myUniverse->mu / RadiusBounds[0]),
                                                          sqrt(this->myUniverse->mu / RadiusBounds[1]) });

            //Step 4: call the base class calcbounds_event_left_side
            this->PeriapseDeparture::calcbounds_event_left_side(RadiusBounds, VelocityMagnitudeBounds, MassBounds);

            //Step 5: create inclination constraint
            if (this->isFirstEventInMission)
            {
                this->Flowerbounds->push_back(cos(this->myOptions->DLA_bounds[1] * math::deg2rad));
                this->Fupperbounds->push_back(cos(this->myOptions->DLA_bounds[0] * math::deg2rad));
                this->Fdescriptions->push_back(this->prefix + "Departure orbit cos(inclination)");

                //dependencies...inclination is dependent on angular momentum, so it is dependent on the entire parking 6-state
                this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                    this->Xindex_rMag,
                    this->Gindex_cosINC_rMAG);
                this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                    this->Xindex_RA,
                    this->Gindex_cosINC_RA);
                this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                    this->Xindex_DEC,
                    this->Gindex_cosINC_DEC);
                this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                    this->Xindex_vMag,
                    this->Gindex_cosINC_vMAG);

                if (this->myOptions->PeriapseBoundaryStateRepresentation == StateRepresentation::SphericalAZFPA)
                {
                    this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                        this->Xindex_AZ,
                        this->Gindex_cosINC_AZ);
                    this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                        this->Xindex_FPA,
                        this->Gindex_cosINC_FPA);
                }
                else if (this->myOptions->PeriapseBoundaryStateRepresentation == StateRepresentation::SphericalRADEC)
                {
                    this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                        this->Xindex_vRA,
                        this->Gindex_cosINC_vRA);
                    this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                        this->Xindex_vDEC,
                        this->Gindex_cosINC_vDEC);
                }
            }

            //Step 6: periapse maneuver velocity magnitude and mass constraint
            {
                double deltav_min = this->myJourneyOptions->initial_impulse_bounds[0];
                double deltav_max = this->myJourneyOptions->initial_impulse_bounds[1];
                this->Xlowerbounds->push_back(deltav_min);
                this->Xupperbounds->push_back(deltav_max);
                this->X_scale_factors->push_back(this->myUniverse->LU / this->myUniverse->TU);
                this->Xdescriptions->push_back(prefix + "magnitude of periapse maneuver");
                this->Xindex_PeriapseManeuver = this->Xdescriptions->size() - 1;

                //C3 affects all three velocity components and mass, and therefore so does periapse impulse
                this->Derivatives_of_StateBeforeEvent.push_back({ this->Xindex_PeriapseManeuver, 3, 1.0 });
                this->dIndex_PeriapseManeuver_vx = this->Derivatives_of_StateBeforeEvent.size() - 1;
                this->Derivatives_of_StateBeforeEvent.push_back({ this->Xindex_PeriapseManeuver, 4, 1.0 });
                this->dIndex_PeriapseManeuver_vy = this->Derivatives_of_StateBeforeEvent.size() - 1;
                this->Derivatives_of_StateBeforeEvent.push_back({ this->Xindex_PeriapseManeuver, 5, 1.0 });
                this->dIndex_PeriapseManeuver_vz = this->Derivatives_of_StateBeforeEvent.size() - 1;

                //if we are not using an LV then we are doing a regular maneuver and this affects mass
                if (!this->useLV && this->hasBipropManeuver)
                {
                    this->Derivatives_of_StateBeforeEvent.push_back({ this->Xindex_PeriapseManeuver, 6, 1.0 });
                    this->dIndex_PeriapseManeuver_mass = this->Derivatives_of_StateBeforeEvent.size() - 1;
                }
            }//end velocity magnitude

            //Step 7: if we are using a launch vehicle, create the launch vehicle mass constraint
            if (this->useLV)
            {
                this->Flowerbounds->push_back(math::SMALL);
                this->Fupperbounds->push_back(1.0);
                this->Fdescriptions->push_back(this->prefix + "spacecraft mass is less than launch vehicle capacity");

                //encoded mass
                this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                    this->Xindex_mass,
                    this->Gindex_LVmassConstraint_encodedMass);

                //periapse maneuver magnitude
                this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                    this->Xindex_PeriapseManeuver,
                    this->Gindex_LVmassConstraint_periapseManeuver);

                //initial parking orbit velocity
                this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                    this->Xindex_vMag,
                    this->Gindex_LVmassConstraint_initialOrbitVelocity);

                //derivative w.r.t parking orbit radius
                this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                    this->Xindex_rMag,
                    this->Gindex_LVmassConstraint_initialOrbitRadius);
            }
        }//end calcbounds_event_left_side()

        void PeriapseLaunchOrImpulsiveDeparture::calcbounds_virtual_propellant_constraints()
        {
            //call the base class propellant tank calcbounds
            BoundaryEventBase::calcbounds_virtual_propellant_constraints();

            //specialty derivative entries for this phase type - burn with departure stage
            //derivatives of virtual chemical fuel with respect to initial mass and v-infinity, in that order
            this->create_sparsity_entry(this->Findex_VirtualChemicalFuelConstraint,
                this->Xindex_mass,
                this->Gindices_dVirtualChemicalFuel_dLeftMass);

            this->create_sparsity_entry(this->Findex_VirtualChemicalFuelConstraint,
                this->Xindex_PeriapseManeuver,
                this->Gindices_dVirtualChemicalFuel_dVinfinity);

            //derivatives of virtual chemical oxidizer with respect to initial mass and v-infinity, in that order
            this->create_sparsity_entry(this->Findex_VirtualChemicalOxidizerConstraint,
                this->Xindex_mass,
                this->Gindices_dVirtualChemicalOxidizer_dLeftMass);

            this->create_sparsity_entry(this->Findex_VirtualChemicalOxidizerConstraint,
                this->Xindex_PeriapseManeuver,
                this->Gindices_dVirtualChemicalOxidizer_dVinfinity);
        }//end calcbounds_virtual_propellant_constraint

        void PeriapseLaunchOrImpulsiveDeparture::calcbounds_deltav_contribution()
        {
            //add derivative with respect to launch v-infinity
            this->create_sparsity_entry(0,
                this->Xindex_PeriapseManeuver,
                this->Gindex_dDeltav_dVinfinity);
        }//end calcbounds_virtual_deltav_constraint

        //******************************************process methods

        void PeriapseLaunchOrImpulsiveDeparture::process_event(const std::vector<doubleType>& X,
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
            if (this->myOptions->objective_type == ObjectiveFunctionType::MINIMIZE_DELTAV &&
                !(this->isFirstEventInMission && !this->myOptions->include_initial_impulse_in_cost))
            {
                this->process_deltav_contribution(X, Xindex, F, Findex, G, needG);
            }

            BoundaryEventBase::process_specialized_constraints(X, Xindex, F, Findex, G, needG);
        }//end process_event

        void PeriapseLaunchOrImpulsiveDeparture::process_event_left_side(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG)
        {
            //Step 1: base class
            this->PeriapseDeparture::process_event_left_side(X, Xindex, F, Findex, G, needG);

            //Step 2: Inclination constraint
            math::Matrix<doubleType> Xpark = this->state_before_event;
            math::Matrix<doubleType> Rpark = this->state_before_event.getSubMatrix1D(0, 2);
            math::Matrix<doubleType> Vpark = this->state_before_event.getSubMatrix1D(3, 5);
            math::Matrix<doubleType> Hpark = Rpark.cross(Vpark);
            doubleType hpark = Hpark.norm();
            if (this->isFirstEventInMission)
            {
                F[Findex++] = Hpark(2) / hpark;
            }

            //Step 3: extract periapse maneuver from decision vector and compute the new periapse velocity vector in cartesian coordinates
            //Step 3.1 extract maneuver magnitude
            this->PeriapseManeuverMagnitude = X[Xindex++];

            //Step 3.2 create a unit vector in the direction of velocity
            doubleType vPark;
            if (this->myOptions->PeriapseBoundaryStateRepresentation == StateRepresentation::SphericalRADEC)
            {
                vPark = X[this->Xindex_vMag];
                doubleType RA = X[this->Xindex_RA];
                doubleType DEC = X[this->Xindex_DEC];
                doubleType vRA = X[this->Xindex_vRA];
                doubleType vDEC = X[this->Xindex_vDEC];
                doubleType cosRA = cos(RA);
                doubleType sinRA = sin(RA);
                doubleType cosDEC = cos(DEC);
                doubleType sinDEC = sin(DEC);
                doubleType cosvRA = cos(vRA);
                doubleType sinvRA = sin(vRA);
                doubleType cosvDEC = cos(vDEC);
                doubleType sinvDEC = sin(vDEC);

                this->Vunit(0) = cosvRA * cosvDEC;
                this->Vunit(1) = sinvRA * cosvDEC;
                this->Vunit(2) = sinvDEC;
            }
            else if (this->myOptions->PeriapseBoundaryStateRepresentation == StateRepresentation::SphericalAZFPA)
            {
                vPark = X[this->Xindex_vMag];
                doubleType RA = X[this->Xindex_RA];
                doubleType DEC = X[this->Xindex_DEC];
                doubleType AZ = X[this->Xindex_AZ];
                doubleType FPA = X[this->Xindex_FPA];
                doubleType cosRA = cos(RA);
                doubleType sinRA = sin(RA);
                doubleType cosDEC = cos(DEC);
                doubleType sinDEC = sin(DEC);
                doubleType cosAZ = cos(AZ);
                doubleType sinAZ = sin(AZ);
                doubleType cosFPA = cos(FPA);
                doubleType sinFPA = sin(FPA);

                this->Vunit(0) = -(sinFPA*(sinAZ*sinRA + cosAZ * cosRA*sinDEC) - cosFPA * cosDEC*cosRA);
                this->Vunit(1) = (sinFPA*(cosRA*sinAZ - cosAZ * sinDEC*sinRA) + cosFPA * cosDEC*sinRA);
                this->Vunit(2) = (cosFPA*sinDEC + cosDEC * cosAZ*sinFPA);
            }

            //Step 3.3 compute the new velocity vector
            this->state_before_event(3) += Vunit(0) * this->PeriapseManeuverMagnitude;
            this->state_before_event(4) += Vunit(1) * this->PeriapseManeuverMagnitude;
            this->state_before_event(5) += Vunit(2) * this->PeriapseManeuverMagnitude;


            //Step 3.4: compute the outgoing C3
            doubleType rPark = X[this->Xindex_rMag];
            doubleType v_p = vPark + this->PeriapseManeuverMagnitude;
            doubleType v_esc2 = 2.0 * this->myUniverse->mu / rPark;
            this->C3 = v_p * v_p - v_esc2;

            //Step 4: if using an LV, do LV things, otherwise apply a biprop burn
            if (this->useLV)
            {

                this->myLaunchVehicle->computePerformance(this->C3, this->myOptions->LV_margin);
                doubleType launch_mass = this->myLaunchVehicle->getDeliveredMass();

                if (launch_mass > this->myJourneyOptions->maximum_mass)
                {
                    launch_mass = this->myJourneyOptions->maximum_mass;
                    this->dm_dPeriapseVelocity = 0.0;
                    this->dm_dParkingOrbitRadius = 0.0;
                }
                else
                {
                    this->dm_dPeriapseVelocity = (this->myLaunchVehicle->getdmdC3() * 2.0 * v_p)_GETVALUE;
                    this->dm_dParkingOrbitRadius = (this->myLaunchVehicle->getdmdC3() * 2.0 * this->myUniverse->mu / rPark / rPark)_GETVALUE;
                }

                F[Findex++] = (launch_mass - this->state_before_event(6)) / this->myJourneyOptions->maximum_mass;
            }
            else //periapse maneuver is a biprop burn OR is free
            {
                if (this->hasBipropManeuver)
                {
                    //perform a maneuver with the current stage
                    this->mySpacecraft->computeChemicalPropulsionPerformance(this->PeriapseManeuverMagnitude,
                        this->state_before_event(6),
                        true,
                        this->ChemicalManeuverType);

                    this->chemical_fuel_used = this->mySpacecraft->getChemFuelConsumedThisManeuver();
                    this->chemical_oxidizer_used = this->mySpacecraft->getChemOxidizerConsumedThisManeuver();

                    this->dChemicalFuel_dVinfinity = this->mySpacecraft->getdFuelConsumedThisManeuver_ddeltav();
                    this->dChemicalOxidizer_dVinfinity = this->mySpacecraft->getdOxidizerConsumedThisManeuver_ddeltav();

                    this->state_before_event(6) -= this->chemical_fuel_used + this->chemical_oxidizer_used;
                }
            }

            //Step 5: derivatives
            if (needG)
            {
                //Step 5.1: Inclination constraint
                if (this->isFirstEventInMission)
                {
                    //Step 5.1.1 generate the derivatives of inclination with respect to inertial state
                    double x = Rpark(0)_GETVALUE;
                    double y = Rpark(1)_GETVALUE;
                    double z = Rpark(2)_GETVALUE;
                    double xdot = Vpark(0)_GETVALUE;
                    double ydot = Vpark(1)_GETVALUE;
                    double zdot = Vpark(2)_GETVALUE;
                    double h = hpark _GETVALUE;
                    double h2 = h * h;
                    double h3 = h2 * h;
                    double Hx = Hpark(0)_GETVALUE;
                    double Hy = Hpark(1)_GETVALUE;
                    double Hz = Hpark(2)_GETVALUE;

                    double dcosINC_dHx = -(Hx*Hz) / h3;
                    double dcosINC_dHy = -(Hy*Hz) / h3;
                    double dcosINC_dHz = (Hx*Hx + Hy*Hy) / h3;

                    math::Matrix<double> dcosINC_dInertialState(6, 1, 0.0);
                    math::Matrix<double> dH_dInertialState(3, 6, 0.0);

                    dH_dInertialState(0, 0) = 0.0;
                    dH_dInertialState(1, 0) = -Vpark(2)_GETVALUE;
                    dH_dInertialState(2, 0) = Vpark(1)_GETVALUE;
                    dH_dInertialState(0, 1) = Vpark(2)_GETVALUE;
                    dH_dInertialState(1, 1) = 0.0;
                    dH_dInertialState(2, 1) = -Vpark(0)_GETVALUE;
                    dH_dInertialState(0, 2) = -Vpark(1)_GETVALUE;
                    dH_dInertialState(1, 2) = Vpark(0)_GETVALUE;
                    dH_dInertialState(2, 2) = 0.0;
                    dH_dInertialState(0, 3) = 0.0;
                    dH_dInertialState(1, 3) = Rpark(2)_GETVALUE;
                    dH_dInertialState(2, 3) = -Rpark(1)_GETVALUE;
                    dH_dInertialState(0, 4) = -Rpark(2)_GETVALUE;
                    dH_dInertialState(1, 4) = 0.0;
                    dH_dInertialState(2, 4) = Rpark(0)_GETVALUE;
                    dH_dInertialState(0, 5) = Rpark(1)_GETVALUE;
                    dH_dInertialState(1, 5) = -Rpark(0)_GETVALUE;
                    dH_dInertialState(2, 5) = 0.0;

                    for (size_t stateIndex = 0; stateIndex < 6; ++stateIndex)
                    {
                        dcosINC_dInertialState(stateIndex) = dcosINC_dHx * dH_dInertialState(0, stateIndex)
                            + dcosINC_dHy * dH_dInertialState(1, stateIndex)
                            + dcosINC_dHz * dH_dInertialState(2, stateIndex);
                    }

                    if (this->myOptions->PeriapseBoundaryStateRepresentation == StateRepresentation::SphericalAZFPA)
                    {
                        //Step 5.1.2: derivative w.r.t. rMag
                        G[this->Gindex_cosINC_rMAG] = this->X_scale_factors->operator[](this->jGvar->operator[](this->Gindex_cosINC_rMAG))
                            * (dcosINC_dInertialState(0) * std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_r[0]])
                                + dcosINC_dInertialState(1) * std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_r[1]])
                                + dcosINC_dInertialState(2) * std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_r[2]]));
                        //Step 5.1.3: derivative w.r.t. RA
                        G[this->Gindex_cosINC_RA] = dcosINC_dInertialState(0) * std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_RA[0]])
                            + dcosINC_dInertialState(1) * std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_RA[1]])
                            + dcosINC_dInertialState(3) * std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_RA[2]])
                            + dcosINC_dInertialState(4) * std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_RA[3]]);
                        //Step 5.1.4: derivative w.r.t. DEC
                        G[this->Gindex_cosINC_DEC] = dcosINC_dInertialState(0) * std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_DEC[0]])
                            + dcosINC_dInertialState(1) * std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_DEC[1]])
                            + dcosINC_dInertialState(2) * std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_DEC[2]])
                            + dcosINC_dInertialState(3) * std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_DEC[3]])
                            + dcosINC_dInertialState(4) * std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_DEC[4]])
                            + dcosINC_dInertialState(5) * std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_DEC[5]]);
                        //Step 5.1.5: derivative w.r.t. vMag
                        G[this->Gindex_cosINC_vMAG] = dcosINC_dInertialState(3) * std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_v[0]])
                            + dcosINC_dInertialState(4) * std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_v[1]])
                            + dcosINC_dInertialState(5) * std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_v[2]]);
                        //Step 5.1.6: derivative w.r.t. AZ
                        G[this->Gindex_cosINC_AZ] = dcosINC_dInertialState(3) * std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_AZ[0]])
                            + dcosINC_dInertialState(4) * std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_AZ[1]])
                            + dcosINC_dInertialState(5) * std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_AZ[2]]);
                        //Step 5.1.7: derivative w.r.t. FPA
                        G[this->Gindex_cosINC_FPA] = dcosINC_dInertialState(3) * std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_FPA[0]])
                            + dcosINC_dInertialState(4) * std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_FPA[1]])
                            + dcosINC_dInertialState(5) * std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_FPA[2]]);
                    }
                    else if (this->myOptions->PeriapseBoundaryStateRepresentation == StateRepresentation::SphericalRADEC)
                    {
                        //Step 5.1.2: derivative w.r.t. rMag
                        G[this->Gindex_cosINC_rMAG] = this->X_scale_factors->operator[](this->jGvar->operator[](this->Gindex_cosINC_rMAG))
                            * (   dcosINC_dInertialState(0) * std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_r[0]])
                                + dcosINC_dInertialState(1) * std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_r[1]])
                                + dcosINC_dInertialState(2) * std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_r[2]]));
                        //Step 5.1.3: derivative w.r.t. RA
                        G[this->Gindex_cosINC_RA] = dcosINC_dInertialState(0) * std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_RA[0]])
                            + dcosINC_dInertialState(1) * std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_RA[1]]);
                        //Step 5.1.4: derivative w.r.t. DEC
                        G[this->Gindex_cosINC_DEC] = dcosINC_dInertialState(0) * std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_DEC[0]])
                            + dcosINC_dInertialState(1) * std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_DEC[1]])
                            + dcosINC_dInertialState(2) * std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_DEC[2]]);
                        //Step 5.1.5: derivative w.r.t. vMag
                        G[this->Gindex_cosINC_vMAG] = dcosINC_dInertialState(3) * std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_v[0]])
                            + dcosINC_dInertialState(4) * std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_v[1]])
                            + dcosINC_dInertialState(5) * std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_v[2]]);
                        //Step 5.1.6: derivative w.r.t. vRA
                        G[this->Gindex_cosINC_vRA] = dcosINC_dInertialState(3) * std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_vRA[0]])
                            + dcosINC_dInertialState(4) * std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_vRA[1]]);
                        //Step 5.1.7: derivative w.r.t. vDEC
                        G[this->Gindex_cosINC_vDEC] = dcosINC_dInertialState(3) * std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_vDEC[0]])
                            + dcosINC_dInertialState(4) * std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_vDEC[1]])
                            + dcosINC_dInertialState(5) * std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_vDEC[2]]);
                    }
                }
                
                //Step 5.2: initial velocity vector
                double Multiplier = ((vPark + this->PeriapseManeuverMagnitude) / vPark) _GETVALUE;
                //Step 5.2.1: w.r.t. periapse maneuver magnitude
                std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_PeriapseManeuver_vx]) = Vunit(0)_GETVALUE;
                std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_PeriapseManeuver_vy]) = Vunit(1)_GETVALUE;
                std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_PeriapseManeuver_vz]) = Vunit(2)_GETVALUE;
                //Step 5.2.2: w.r.t. parking orbit velocity magnitude
                //no changes
                if (this->myOptions->PeriapseBoundaryStateRepresentation == StateRepresentation::SphericalRADEC)
                {
                    //Step 5.2.2: w.r.t. RA
                    //no changes
                    //Step 5.2.3: w.r.t. DEC
                    //no changes
                    //Step 5.2.4: w.r.t. vRA
                    std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_vRA[0]]) *= Multiplier;
                    std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_vRA[1]]) *= Multiplier;
                    ////Step 5.2.5: w.r.t. vDEC
                    std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_vDEC[0]]) *= Multiplier;
                    std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_vDEC[1]]) *= Multiplier;
                    std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_vDEC[2]]) *= Multiplier;
                }
                else if (this->myOptions->PeriapseBoundaryStateRepresentation == StateRepresentation::SphericalAZFPA)
                {
                    //Step 5.2.2: w.r.t. RA
                    std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_RA[2]]) *= Multiplier;
                    std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_RA[3]]) *= Multiplier;
                    //Step 5.2.3: w.r.t. DEC 
                    std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_DEC[3]]) *= Multiplier;
                    std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_DEC[4]]) *= Multiplier;
                    std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_DEC[5]]) *= Multiplier;
                    //Step 5.2.4: w.r.t. AZ 
                    std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_AZ[0]]) *= Multiplier;
                    std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_AZ[1]]) *= Multiplier;
                    std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_AZ[2]]) *= Multiplier;
                    ////Step 5.2.5: w.r.t. FPA 
                    std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_FPA[0]]) *= Multiplier;
                    std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_FPA[1]]) *= Multiplier;
                    std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_FPA[2]]) *= Multiplier;
                }

                //Step 5.3: mass
                if (this->useLV)
                {
                    size_t Gindex, Xindex;
                    //Step 5.3.1: derivative w.r.t. left state mass
                    //state left mass is subtracted off the right-hand side
                    Gindex = this->Gindex_LVmassConstraint_encodedMass;
                    Xindex = this->jGvar->operator[](Gindex);
                    G[Gindex] = -this->X_scale_factors->operator[](Xindex) / this->myJourneyOptions->maximum_mass;

                    //Step 5.3.2: derivative w.r.t. periapse maneuver manitude
                    Gindex = this->Gindex_LVmassConstraint_periapseManeuver;
                    Xindex = this->jGvar->operator[](Gindex);
                    G[Gindex] = this->X_scale_factors->operator[](Xindex) * this->dm_dPeriapseVelocity / this->myJourneyOptions->maximum_mass;

                    //Step 5.3.3: derivative w.r.t. parking orbit velocity magnitude
                    Gindex = this->Gindex_LVmassConstraint_initialOrbitVelocity;
                    Xindex = this->jGvar->operator[](Gindex);
                    G[Gindex] = this->X_scale_factors->operator[](Xindex) * this->dm_dPeriapseVelocity / this->myJourneyOptions->maximum_mass;

                    //Step 5.3.4: derivative w.r.t parking orbit radius
                    Gindex = this->Gindex_LVmassConstraint_initialOrbitRadius;
                    Xindex = this->jGvar->operator[](Gindex);
                    G[Gindex] = this->X_scale_factors->operator[](Xindex) * this->dm_dParkingOrbitRadius / this->myJourneyOptions->maximum_mass;

                }
                else if (this->hasBipropManeuver)
                {
                    std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_PeriapseManeuver_mass]) = this->mySpacecraft->getdSpacecraftMassThisManeuver_ddeltav();
                }
            }
        }//process_event_left_side()

        void PeriapseLaunchOrImpulsiveDeparture::process_virtual_propellant_constraints(const std::vector<doubleType>& X,
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
                Gentry = (this->chemical_fuel_used / this->state_before_event(6)) _GETVALUE;
                G[Gindex] = this->X_scale_factors->operator[](Xindex) * -Gentry
                    * this->myUniverse->continuity_constraint_scale_factors(6);

                Gindex = this->Gindices_dVirtualChemicalOxidizer_dLeftMass;
                Xindex = this->jGvar->operator[](Gindex);
                Gentry = (this->chemical_oxidizer_used / this->state_before_event(6)) _GETVALUE;
                G[Gindex] = this->X_scale_factors->operator[](Xindex) * -Gentry
                    * this->myUniverse->continuity_constraint_scale_factors(6);

                //derivatives with respect to v-infinity
                Gindex = this->Gindices_dVirtualChemicalFuel_dVinfinity;
                Xindex = this->jGvar->operator[](Gindex);
                G[Gindex] = this->X_scale_factors->operator[](Xindex) * -this->dChemicalFuel_dVinfinity
                    * this->myUniverse->continuity_constraint_scale_factors(6);

                Gindex = this->Gindices_dVirtualChemicalOxidizer_dVinfinity;
                Xindex = this->jGvar->operator[](Gindex);
                G[Gindex] = this->X_scale_factors->operator[](Xindex) * -this->dChemicalOxidizer_dVinfinity
                    * this->myUniverse->continuity_constraint_scale_factors(6);
            }
        }//end process_virtual_propellant_constraints

        void PeriapseLaunchOrImpulsiveDeparture::process_deltav_contribution(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //add the event delta-v
            this->EventDeterministicDeltav = this->PeriapseManeuverMagnitude;

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
        void PeriapseLaunchOrImpulsiveDeparture::output(std::ofstream& outputfile,
            const double& launchdate,
            size_t& eventcount)
        {
            this->mySpacecraft->setActiveStage(this->stageIndex);

            //base class
            PeriapseDeparture::output(outputfile, launchdate, eventcount);

            std::string event_type;
            if (!this->hasBipropManeuver)
                event_type = "launch";
            else
                event_type = "departure";

            std::string boundary_name = "periapse";

            math::Matrix<doubleType> empty3vector(3, 1, 0.0);
            math::Matrix<doubleType> dV = this->Vunit * this->PeriapseManeuverMagnitude;

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

            //compute RLA and DLA
            //first compute angular momentum
            math::Matrix<doubleType> R = this->state_after_event.getSubMatrix1D(0, 2);
            math::Matrix<doubleType> V = this->state_after_event.getSubMatrix1D(3, 5);

            math::Matrix<doubleType> hVec = R.cross(V);
            double mu = this->myUniverse->central_body.mu;
            double r = R.norm()_GETVALUE;
            double v = V.norm()_GETVALUE;
            double h = hVec.norm()_GETVALUE;

            //eccentricity vector
            math::Matrix<doubleType> eVec = (R * (v * v - mu / r) - V * R.dot(V)) / mu;

            //outgoing velocity asymptote
            math::Matrix<doubleType> sHat = (hVec.cross(eVec) * sqrt(C3) / mu - eVec) * 1 / (1.0 + C3 * (h * h / mu / mu));

            //RLA and DLA
            this->RLA = atan2(sHat(1), sHat(0));
            this->DLA = asin(sHat(2));

            write_output_line(outputfile,
                eventcount,
                event_type,
                boundary_name,
                this->EventTimeWidth,
                this->state_after_event.getSubMatrix1D(0, 2).norm() - this->myUniverse->central_body.radius,
                0.0,
                0.0,
                this->RLA,
                this->DLA,
                this->C3,
                this->state_after_event,
                dV,
                empty3vector,
                this->PeriapseManeuverMagnitude,
                0.0,
                this->hasFixedInitialMass || this->useLV ? Isp : -1.0,
                this->mySpacecraft->getAvailablePower(),
                0.0,
                0,
                0.0,
                "none");
        }

        void PeriapseLaunchOrImpulsiveDeparture::output_maneuver_and_target_spec(std::ofstream& maneuver_spec_file, std::ofstream& target_spec_file, bool& haveManeuverNeedTarget)
        { 

            //Step 1: maneuver spec
            //Step 1.1: instantiate and populate a maneuver spec object
            maneuver_spec_line myManeuverSpecLine(this->name);

            double Isp = this->ChemicalManeuverType == PropulsionSystemChoice::Biprop ? this->mySpacecraft->getBipropIsp() : this->mySpacecraft->getMonopropIsp();
            doubleType massFlowRate = this->mySpacecraft->getchemthrust() / Isp / this->myOptions->g0;
            doubleType mass_before_maneuver = this->state_before_event(6) + this->chemical_fuel_used + this->chemical_oxidizer_used;
            doubleType maneuverDuration = (mass_before_maneuver - this->state_before_event(6)) / massFlowRate;

            myManeuverSpecLine.append_maneuver_spec_item("EME2000",
                this->state_before_event(7),
                this->Vunit,
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