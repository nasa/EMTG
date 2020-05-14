
// EMTG: Evolutionary Mission Trajectory Generator
// An open-source global optimization tool for preliminary mission design
// Provided by NASA Goddard Space Flight Center
//
// Copyright (c) 2014 - 2017 United States Government as represented by the
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

#include "orbit_departure_2D.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        OrbitDeparture2D::OrbitDeparture2D(const std::string& name,
            const size_t& journeyIndex,
            const size_t& phaseIndex,
            size_t& stageIndex,
            Astrodynamics::universe* Universe,
            HardwareModels::Spacecraft* mySpacecraft,
            missionoptions* myOptions,
            ArrivalEvent* PreviousPhaseArrivalEvent) :
            DepartureEvent::DepartureEvent(name,
                journeyIndex,
                phaseIndex,
                stageIndex,
                Universe,
                mySpacecraft,
                myOptions,
                PreviousPhaseArrivalEvent)
        {
            if (!this->isFirstEventInMission)
                this->hasWaitTime = true;
            
            //orbit departures ALWAYS have a chemical maneuver
            this->hasBipropManeuver = true;

            if (!this->LeftBoundaryIsABody)
            {
                std::cerr << "j" << this->journeyIndex << "p" << this->phaseIndex << ": The left boundary of an OrbitDeparture2D event MUST be a body." << std::endl;
                throw 8;
            }
        }

        //******************************************calcbounds methods

        //calcbounds
        void OrbitDeparture2D::calcbounds()
        {
            this->calcbounds_event_left_side();

            this->calcbounds_event_main();

            this->calcbounds_event_right_side();

            this->calcbounds_virtual_propellant_constraints();

            //virtual delta-v variable and constraint exist only if this event has deterministic maneuvers
            if (!(this->isFirstEventInMission && !this->myOptions->include_initial_impulse_in_cost))
            {
                this->calcbounds_deltav_contribution();
            }
        }//end calcbounds()

        void OrbitDeparture2D::calcbounds_event_main()
        {
            //calcbounds for launch or direct insertion
            //First, we have outgoing velocity
            Xlowerbounds->push_back(this->myOptions->Journeys[journeyIndex].journey_initial_impulse_bounds[0]);
            Xupperbounds->push_back(this->myOptions->Journeys[journeyIndex].journey_initial_impulse_bounds[1]);
            X_scale_factors->push_back(this->myUniverse->LU / this->myUniverse->TU);
            Xdescriptions->push_back(prefix + "magnitude of outgoing velocity asymptote");

            //RLA
            Xlowerbounds->push_back(this->myOptions->RLA_bounds[0] * math::PI / 180.0);
            Xupperbounds->push_back(this->myOptions->RLA_bounds[1] * math::PI / 180.0);
            X_scale_factors->push_back(math::TwoPI);
            Xdescriptions->push_back(prefix + "RA of departure asymptote");


            //DLA
            if (this->journeyIndex == 0 && this->LeftBoundaryIsABody) //if this is the first journey and we are leaving from a planet, i.e. if this is a launch
            {
                Xlowerbounds->push_back(this->myOptions->DLA_bounds[0] * math::PI / 180.0);
                Xupperbounds->push_back(this->myOptions->DLA_bounds[1] * math::PI / 180.0);
            }
            else
            {
                Xlowerbounds->push_back(-math::PI / 2.0);
                Xupperbounds->push_back(math::PI / 2.0);
            }
            X_scale_factors->push_back(math::TwoPI);
            Xdescriptions->push_back(prefix + "DEC of departure asymptote");


        }//end calcbounds_event_main

        void OrbitDeparture2D::calcbounds_event_right_side()
        {
            //derivatives of the right-side linkage constraints with respect to event main decision variables
            BoundaryEventBase::calcbounds_event_right_side();

            //v-infinity
            for (size_t stateIndex = 3; stateIndex < 7; ++stateIndex) //does not affect position or epoch
            {
                //first we need to find out the F index of this state's right-side linkage constraint
                size_t Findex = this->F_index_of_first_constraint_in_this_event;
                while (Findex < this->Fdescriptions->size())
                {
                    if (Fdescriptions->at(Findex).find(prefix + "event right " + stateNames[stateIndex] + " linkage") < 1024)
                        break;
                    ++Findex;
                }

                this->create_sparsity_entry(Findex,
                    this->X_index_of_first_decision_variable_in_this_event,
                    true,
                    prefix + "magnitude of outgoing velocity asymptote",
                    Gindices_dStateAfterEventLinkage_dVinfinity);
            }

            //RLA
            for (size_t stateIndex = 3; stateIndex < 5; ++stateIndex) //does not affect position, zdot, mass, or epoch
            {
                //first we need to find out the F index of this state's right-side linkage constraint
                size_t Findex = this->F_index_of_first_constraint_in_this_event;
                while (Findex < this->Fdescriptions->size())
                {
                    if (Fdescriptions->at(Findex).find(prefix + "event right " + stateNames[stateIndex] + " linkage") < 1024)
                        break;
                    ++Findex;
                }

                this->create_sparsity_entry(Findex,
                    this->X_index_of_first_decision_variable_in_this_event,
                    true,
                    prefix + "RA of departure asymptote",
                    Gindices_dStateAfterEventLinkage_dRLA);
            }

            //DLA
            for (size_t stateIndex = 3; stateIndex < 6; ++stateIndex) //does not affect position, mass, or epoch
            {
                //first we need to find out the F index of this state's right-side linkage constraint
                size_t Findex = this->F_index_of_first_constraint_in_this_event;
                while (Findex < this->Fdescriptions->size())
                {
                    if (Fdescriptions->at(Findex).find(prefix + "event right " + stateNames[stateIndex] + " linkage") < 1024)
                        break;
                    ++Findex;
                }

                this->create_sparsity_entry(Findex,
                    this->X_index_of_first_decision_variable_in_this_event,
                    true,
                    prefix + "DEC of departure asymptote",
                    Gindices_dStateAfterEventLinkage_dDLA);
            }
        }


        void OrbitDeparture2D::calcbounds_virtual_propellant_constraints()
        {
            //call the base class propellant tank calcbounds
            BoundaryEventBase::calcbounds_virtual_propellant_constraints();

            //specialty derivative entries for this phase type - burn with departure stage
            //derivatives of virtual chemical fuel with respect to initial mass and v-infinity, in that order
            this->create_sparsity_entry(this->Findex_VirtualChemicalFuelConstraint,
                this->X_index_of_first_decision_variable_in_this_event,
                true,
                "mass",
                Gindices_dVirtualChemicalFuel_dLeftMass);

            this->create_sparsity_entry(this->Findex_VirtualChemicalFuelConstraint,
                this->X_index_of_first_decision_variable_in_this_event,
                true,
                "magnitude of outgoing velocity asymptote",
                Gindices_dVirtualChemicalFuel_dVinfinity);

            //derivatives of virtual chemical oxidizer with respect to initial mass and v-infinity, in that order
            this->create_sparsity_entry(this->Findex_VirtualChemicalOxidizerConstraint,
                this->X_index_of_first_decision_variable_in_this_event,
                true,
                "mass",
                Gindices_dVirtualChemicalOxidizer_dLeftMass);

            this->create_sparsity_entry(this->Findex_VirtualChemicalOxidizerConstraint,
                this->X_index_of_first_decision_variable_in_this_event,
                true,
                "magnitude of outgoing velocity asymptote",
                Gindices_dVirtualChemicalOxidizer_dVinfinity);
        }//end calcbounds_virtual_propellant_constraints()

        void OrbitDeparture2D::calcbounds_deltav_contribution()
        {
            //add derivative with respect to launch v-infinity
            this->create_sparsity_entry(0,
                this->X_index_of_first_decision_variable_in_this_event,
                true,
                "magnitude of outgoing velocity asymptote",
                Gindex_dVirtualEventTotalDeltav_dVinfinity);
        }//end calcbounds_virtual_deltav_constraint

        //******************************************process methods

        void OrbitDeparture2D::process_event(const std::vector<doubleType>& X,
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

            this->process_virtual_propellant_constraints(X, Xindex, F, Findex, G, needG);

            //virtual delta-v variable and constraint exist only if this event has deterministic maneuvers
            if (!(this->isFirstEventInMission && !this->myOptions->include_initial_impulse_in_cost))
            {
                this->process_deltav_contribution(X, Xindex, F, Findex, G, needG);
            }
        }//end process_event


        void OrbitDeparture2D::
            process_event_main(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG)
        {
            //Step 1: read in the decision variables and apply launch mass constraint
            //C3
            doubleType vinf_out = X[Xindex++];
            this->C3 = vinf_out * vinf_out;

            //RLA
            this->RLA = X[Xindex++];

            //DLA
            this->DLA = X[Xindex++];

            //Step 2: construct state_after_event_propagated
            //position, mass, and time do not change by default, although mass may later be modified by the departure maneuver
            for (size_t stateIndex = 0; stateIndex < 8; ++stateIndex)
                this->state_after_event_propagated(stateIndex) = this->state_before_event(stateIndex);
            doubleType cosDLA = cos(DLA);
            doubleType sinDLA = sin(DLA);
            doubleType cosRLA = cos(RLA);
            doubleType sinRLA = sin(RLA);
            this->state_after_event_propagated(3) += vinf_out * cosRLA * cosDLA;
            this->state_after_event_propagated(4) += vinf_out * sinRLA * cosDLA;
            this->state_after_event_propagated(5) += vinf_out * sinDLA;

            //Step 3: perform a maneuver with the stage
            {
                //Step 3.1 compute the required delta-v magnitude
                //TODO later we might make SMA and ECC variables?
                double SMA = this->myOptions->Journeys[this->journeyIndex].journey_departure_elements[0];
                double ECC = this->myOptions->Journeys[this->journeyIndex].journey_departure_elements[1];

                double r_p = SMA * (1 - ECC);

                //find the velocity of the spacecraft at periapse of the inbound hyperbola
                doubleType v_p_hyperbola = sqrt(this->C3 + 2 * this->myBody->mu / r_p);

                //find the velocity of the spacecraft at periapse of the initial elliptical orbit
                doubleType v_p_ellipse = sqrt(this->myBody->mu * (2 / r_p - 1 / SMA));

                //calculate the final deltaV in km/s
                this->dvDeparture = v_p_hyperbola - v_p_ellipse;

                //calculate the derivative of the maneuver with respect to inbound v-infinity magnitude
                this->ddvDeparture_dVinfinityOut = ( (v_p_hyperbola - v_p_ellipse) / this->dvDeparture * vinf_out / v_p_hyperbola )_GETVALUE;

                //Step 3.2 execute the maneuver
                this->mySpacecraft->computeChemicalPropulsionPerformance(this->dvDeparture, this->state_before_event(6), true, PropulsionSystemChoice::Biprop);

                this->chemical_fuel_used = this->mySpacecraft->getChemFuelConsumedThisManeuver();
                this->chemical_oxidizer_used = this->mySpacecraft->getChemOxidizerConsumedThisManeuver();

                this->dChemicalFuel_dVinfinity = this->mySpacecraft->getdFuelConsumedThisManeuver_ddeltav() * this->ddvDeparture_dVinfinityOut;
                this->dChemicalOxidizer_dVinfinity = this->mySpacecraft->getdOxidizerConsumedThisManeuver_ddeltav() * this->ddvDeparture_dVinfinityOut;

                this->dm_dVinfinityOut = this->mySpacecraft->getdSpacecraftMassThisManeuver_ddeltav() * this->ddvDeparture_dVinfinityOut;

                this->state_after_event_propagated(6) -= this->chemical_fuel_used + this->chemical_oxidizer_used;
            }

            //Step 4: construct the ETM
            this->ETM(6, 6) = (this->state_after_event_propagated(6) / this->state_before_event(6)) _GETVALUE;

            //Step 5: populate the derivatives of the right-side linkage constraint
            if (needG)
            {
                math::Matrix<doubleType> dstate_dDecisionVariable(8, 1, 0.0);
                math::Matrix<doubleType> dPropagatedRightState_dDecisionVariable(8, 1, 0.0);

                //derivatives with respect to v-infinity-out - affects velocity components and (sometimes) mass
                dstate_dDecisionVariable(3) = cosRLA * cosDLA;
                dstate_dDecisionVariable(4) = sinRLA * cosDLA;
                dstate_dDecisionVariable(5) = sinDLA;
                dstate_dDecisionVariable(6) = this->dm_dVinfinityOut / this->ETM(6, 6); //have to take out the ETM's modification to the state

                dPropagatedRightState_dDecisionVariable = this->ETM * dstate_dDecisionVariable;

                for (size_t stateIndex = 3; stateIndex < 7; ++stateIndex)
                {
                    size_t Gindex = this->Gindices_dStateAfterEventLinkage_dVinfinity[stateIndex - 3];
                    size_t Xindex = this->myOptions->jGvar[Gindex];

                    G[Gindex] = this->X_scale_factors->operator[](Xindex) * -dPropagatedRightState_dDecisionVariable(stateIndex) _GETVALUE
                        * this->myUniverse->continuity_constraint_scale_factors(stateIndex);
                }

                //derivatives with respect to RLA - affects velocity components
                dstate_dDecisionVariable.assign_zeros();
                dstate_dDecisionVariable(3) = vinf_out * -sinRLA * cosDLA;
                dstate_dDecisionVariable(4) = vinf_out * cosRLA * cosDLA;

                dPropagatedRightState_dDecisionVariable = this->ETM * dstate_dDecisionVariable;

                for (size_t stateIndex = 3; stateIndex < 5; ++stateIndex)
                {
                    size_t Gindex = this->Gindices_dStateAfterEventLinkage_dRLA[stateIndex - 3];
                    size_t Xindex = this->myOptions->jGvar[Gindex];

                    G[Gindex] = this->X_scale_factors->operator[](Xindex) * -dPropagatedRightState_dDecisionVariable(stateIndex) _GETVALUE
                        * this->myUniverse->continuity_constraint_scale_factors(stateIndex);
                }

                //derivatives with respect to DLA - affects velocity components
                dstate_dDecisionVariable.assign_zeros();
                dstate_dDecisionVariable(3) = vinf_out * cosRLA * -sinDLA;
                dstate_dDecisionVariable(4) = vinf_out * sinRLA * -sinDLA;
                dstate_dDecisionVariable(5) = vinf_out * cosDLA;

                dPropagatedRightState_dDecisionVariable = this->ETM * dstate_dDecisionVariable;

                for (size_t stateIndex = 3; stateIndex < 6; ++stateIndex)
                {
                    size_t Gindex = this->Gindices_dStateAfterEventLinkage_dDLA[stateIndex - 3];
                    size_t Xindex = this->myOptions->jGvar[Gindex];

                    G[Gindex] = this->X_scale_factors->operator[](Xindex) * -dPropagatedRightState_dDecisionVariable(stateIndex) _GETVALUE
                        * this->myUniverse->continuity_constraint_scale_factors(stateIndex);
                }
            }//end derivatives of right-side linkage constraint

        }//end process_event_main

        void OrbitDeparture2D::process_virtual_propellant_constraints(const std::vector<doubleType>& X,
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
                Xindex = this->myOptions->jGvar[Gindex];
                Gentry = (this->chemical_fuel_used / this->state_before_event(6)) _GETVALUE;
                G[Gindex] = this->X_scale_factors->operator[](Xindex) * -Gentry
                    * this->myUniverse->continuity_constraint_scale_factors(6);

                Gindex = this->Gindices_dVirtualChemicalOxidizer_dLeftMass;
                Xindex = this->myOptions->jGvar[Gindex];
                Gentry = (this->chemical_oxidizer_used / this->state_before_event(6)) _GETVALUE;
                G[Gindex] = this->X_scale_factors->operator[](Xindex) * -Gentry
                    * this->myUniverse->continuity_constraint_scale_factors(6);

                //derivatives with respect to v-infinity
                Gindex = this->Gindices_dVirtualChemicalFuel_dVinfinity;
                Xindex = this->myOptions->jGvar[Gindex];
                G[Gindex] = this->X_scale_factors->operator[](Xindex) * -this->dChemicalFuel_dVinfinity
                    * this->myUniverse->continuity_constraint_scale_factors(6);

                Gindex = this->Gindices_dVirtualChemicalOxidizer_dVinfinity;
                Xindex = this->myOptions->jGvar[Gindex];
                G[Gindex] = this->X_scale_factors->operator[](Xindex) * -this->dChemicalOxidizer_dVinfinity
                    * this->myUniverse->continuity_constraint_scale_factors(6);
            }
        }//end process_virtual_propellant_constraints

        void OrbitDeparture2D::process_deltav_contribution(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //add the event delta-v
            this->EventDeterministicDeltav = this->dvDeparture;

            //event-specific
            if (needG)
            {
                //derivatives with respect to v-infinity
                size_t Gindex = this->Gindex_dVirtualEventTotalDeltav_dVinfinity;
                size_t Xindex = this->myOptions->jGvar[Gindex];
                G[Gindex] = this->X_scale_factors->operator[](Xindex) * this->ddvDeparture_dVinfinityOut;
            }
        }//end process_virtual_deltav_constraint()

        //******************************************output methods
        void OrbitDeparture2D::output(std::ofstream& outputfile,
            const double& launchdate,
            size_t& eventcount)
        {
            this->mySpacecraft->setActiveStage(this->stageIndex);

            DepartureEvent::output(outputfile, launchdate, eventcount);

            std::string event_type = "departure";

            std::string boundary_name = this->myBody->name;

            math::Matrix<doubleType> empty3vector(3, 1, 0.0);
            math::Matrix<doubleType> dV(3, 1, std::vector<doubleType>({ sqrt(this->C3) * cos(this->RLA) * cos(this->DLA),
                sqrt(this->C3) * sin(this->RLA) * cos(this->DLA),
                sqrt(this->C3) * sin(this->DLA) }));

            this->mySpacecraft->computePowerState(this->state_after_event.getSubMatrix1D(0, 3).norm() / this->myUniverse->LU,
                (this->EventRightEpoch - launchdate), this->myOptions->power_margin);

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
                this->dvDeparture,
                0.0,
                this->mySpacecraft->getBipropIsp(),
                this->mySpacecraft->getProducedPower(),
                0.0,
                0,
                0.0,
                -1);
        }
    }//end namespace BoundaryEvents
}//end namespace EMTG