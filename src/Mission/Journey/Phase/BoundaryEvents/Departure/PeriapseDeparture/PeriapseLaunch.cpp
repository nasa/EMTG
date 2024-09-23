
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

#include "PeriapseLaunch.h"
#include "StateRepresentationFactory.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        PeriapseLaunch::PeriapseLaunch(const std::string& name,
                                       const size_t& journeyIndex,
                                       const size_t& phaseIndex,
                                       size_t& stageIndex,
                                       Astrodynamics::universe* Universe,
                                       HardwareModels::Spacecraft* mySpacecraft,
                                       HardwareModels::LaunchVehicle* myLaunchVehicle,
                                       missionoptions* myOptions,
                                       ArrivalEvent* PreviousPhaseArrivalEvent) :
            PeriapseBoundary()
        {
            this->myLaunchVehicle = myLaunchVehicle;

            this->initialize(name,
                journeyIndex,
                phaseIndex,
                stageIndex,
                Universe,
                mySpacecraft,
                myOptions,
                PreviousPhaseArrivalEvent);
        }

        void PeriapseLaunch::initialize(const std::string& name,
            const size_t& journeyIndex,
            const size_t& phaseIndex,
            size_t& stageIndex,
            Astrodynamics::universe* Universe,
            HardwareModels::Spacecraft* mySpacecraft,
            missionoptions* myOptions,
            ArrivalEvent* PreviousPhaseArrivalEvent)
        {
            this->isFirstEventInMission = true;

            this->FreePointBoundary::initialize(name,
                journeyIndex,
                phaseIndex,
                stageIndex,
                Universe,
                mySpacecraft,
                myOptions);

            this->periapseDistanceBounds[0] = this->myJourneyOptions->PeriapseDeparture_altitude_bounds[0] + this->myUniverse->central_body_radius;
            this->periapseDistanceBounds[1] = this->myJourneyOptions->PeriapseDeparture_altitude_bounds[1] + this->myUniverse->central_body_radius;

            this->myStateRepresentationEnum = StateRepresentation::OutgoingBplane;

            this->myStateRepresentation = Astrodynamics::CreateStateRepresentation(this->myStateRepresentationEnum, this->myUniverse->mu);

            if (this->myLaunchVehicle->getName() == "Fixed_Initial_Mass")
            {
                this->hasFixedInitialMass = !this->myOptions->allow_initial_mass_to_vary;
                this->useLV = false;
            }
            else //any LV model
            {
                this->hasFixedInitialMass = false;
                this->useLV = true;
            }
           
            this->hasFixedInitialMass = !this->myOptions->allow_initial_mass_to_vary;

            //this flag exists so that MinimizeInitialImpulseObjective knows that there is an initial velocity variable to minimize. It will actually minimize launch v-infinity instead.
            this->hasManeuver = true;

            this->LeftBoundaryIsABody = false;

            //otherwise kaboom, because we're based on FreePointBoundary
            this->AllowStateToPropagate = false;
            this->myEncodedReferenceFrame = ReferenceFrame::ICRF;
            this->ReferenceEpoch = 51545.0 * 86400.0;


            //this check has to be at the end because otherwise the throw actually causes a delete of a non-allocated StateRepresentation object.
            if (!(this->journeyIndex == 0 && this->phaseIndex == 0))
            {
                throw std::invalid_argument("PeriapseLaunch may only be the first event in the mission. Halting.");
            }

            //construct boundary constraints
            this->construct_boundary_constraints();
        }//end initialize

        //******************************************calcbounds methods

        //calcbounds
        void PeriapseLaunch::calcbounds(std::vector<size_t> timeVariables)
        {
            //we have to have this line because LaunchOrDirectInsertion does not call the base calcbounds_event_left_side()
            this->X_index_of_first_decision_variable_in_this_event = this->Xdescriptions->size();

            this->calcbounds_event_left_side(timeVariables);

            //this->calcbounds_event_main(); nothing really happens here

            this->PeriapseBoundary::calcbounds_event_right_side();

            //delta-v if applicable
            if (this->myOptions->objective_type == ObjectiveFunctionType::MINIMIZE_DELTAV && this->myOptions->include_initial_impulse_in_cost)
            {
                this->calcbounds_deltav_contribution();
            }

            this->calcbounds_specialized_constraints();
        }//end calcbounds()

        void PeriapseLaunch::calcbounds_event_left_side(std::vector<size_t> timeVariables)
        {

            //Step 1: compute state bounds in degrees so that they can be converted properly by FreePointBoundary
            std::vector< std::tuple<double, double> > stateBounds;
            //VINF
            stateBounds.push_back({ this->myJourneyOptions->initial_impulse_bounds[0], this->myJourneyOptions->initial_impulse_bounds[1] });
            //RHA
            stateBounds.push_back({ this->myOptions->RLA_bounds[0], this->myOptions->RLA_bounds[1] });
            //DHA
            stateBounds.push_back({ this->myOptions->DLA_bounds[0], this->myOptions->DLA_bounds[1] });
            //BRADIUS
            stateBounds.push_back({ math::SMALL, this->myUniverse->r_SOI });
            //BTHETA
            stateBounds.push_back({ -720.0, 720.0 });
            //TA
            stateBounds.push_back({ -math::SMALL, math::SMALL });
            //mass
            if (this->hasFixedInitialMass)
            {
                stateBounds.push_back({ this->myJourneyOptions->maximum_mass - 1.0e-13, this->myJourneyOptions->maximum_mass });
            }
            else
            {
                stateBounds.push_back({ 1.0e-13, this->myJourneyOptions->maximum_mass });
            }
            //initial epoch
            stateBounds.push_back({ this->myOptions->launch_window_open_date + this->myJourneyOptions->wait_time_bounds[0], this->myOptions->launch_window_open_date + this->myJourneyOptions->wait_time_bounds[1] });

            //Step 2: call the base class calcbounds_event_left_side
            this->FreePointBoundary::calcbounds_event_left_side(stateBounds, timeVariables);

            //we need to track the Xindex of the v-infinity magnitude variable that we just assigned
            this->Xindex_vinfinity_magnitude = this->Xindex_encoded_state[0];

            //Step 3: periapse distance constraint
            this->Flowerbounds->push_back(this->periapseDistanceBounds[0] / this->myUniverse->LU);
            this->Fupperbounds->push_back(this->periapseDistanceBounds[1] / this->myUniverse->LU);
            this->Fdescriptions->push_back(this->prefix + "distance from central body");

            for (size_t dIndex = 0; dIndex < this->Derivatives_of_StateBeforeEvent.size(); ++dIndex)
            {
                std::tuple<size_t, size_t, double>& derivativeEntry = this->Derivatives_of_StateBeforeEvent[dIndex];

                size_t stateIndex = std::get<1>(derivativeEntry);

                if (stateIndex < 3) //position
                {
                    size_t Xindex = std::get<0>(derivativeEntry);
                    this->dIndices_distance_constraint.push_back(dIndex);

                    this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                        Xindex,
                        this->Gindices_distance_constraint);
                }
            }//end derivatives for distance constraint

            //Step 4: create the launch vehicle mass constraint
            if (this->useLV)
            {
                this->Flowerbounds->push_back(-math::SMALL);
                if (this->myOptions->allow_initial_mass_to_vary)
                {
                    this->Fupperbounds->push_back(1.0);//allow the initial mass to be UP TO the launcher's capacity
                }
                else
                {
                    this->Fupperbounds->push_back(math::SMALL);//force the initial mass to match the launcher's capacity
                }
                this->Fdescriptions->push_back(this->prefix + "spacecraft mass is less than launch vehicle capacity");

                //encoded mass
                this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                    this->Xindex_encoded_state[6],
                    this->Gindex_LVmassConstraint_encodedMass);

                //v-infinity magnitude
                this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                    this->Xindex_vinfinity_magnitude,
                    this->Gindex_LVmassConstraint_vinfinity_magnitude);
            }
            
            //Step 5: mass multipliers
            this->calcbounds_mass_multipliers();
        }//end calcbounds_event_left_side()

        void PeriapseLaunch::calcbounds_deltav_contribution()
        {
            //add derivative with respect to launch v-infinity
            this->create_sparsity_entry(0,
                this->Xindex_vinfinity_magnitude,
                this->Gindex_dDeltav_dVinfinity);
        }//end calcbounds_virtual_deltav_constraint

        //******************************************process methods

        void PeriapseLaunch::process_event(const std::vector<doubleType>& X,
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
            if (this->myOptions->objective_type == ObjectiveFunctionType::MINIMIZE_DELTAV && this->myOptions->include_initial_impulse_in_cost)
            {
                this->process_deltav_contribution(X, Xindex, F, Findex, G, needG);
            }

            BoundaryEventBase::process_specialized_constraints(X, Xindex, F, Findex, G, needG);
        }//end process_event

        void PeriapseLaunch::process_event_left_side(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG)
        {
            //Step 0: C3, RLA, DLA
            doubleType vinf = X[this->Xindex_vinfinity_magnitude];
            this->C3 = vinf * vinf;
            this->RLA = X[this->Xindex_vinfinity_magnitude + 1];
            this->DLA = X[this->Xindex_vinfinity_magnitude + 2];

            //Step 1: base free point boundary
            this->FreePointBoundary::process_event_left_side(X, Xindex, F, Findex, G, needG);

            //Step 2: periapse distance constraint
            {
                doubleType distance = this->state_before_event.getSubMatrix1D(0, 2).norm();

                F[Findex++] = distance / this->myUniverse->LU;

                if (needG)
                {
                    //first zero out the derivatives
                    for (size_t Gindex : this->Gindices_distance_constraint)
                    {
                        G[Gindex] = 0.0;
                    }

                    //now assign the entries
                    for (size_t entryIndex = 0; entryIndex < this->dIndices_distance_constraint.size(); ++entryIndex )
                    {
                        size_t dIndex = this->dIndices_distance_constraint[entryIndex];

                        std::tuple<size_t, size_t, double>& derivativeEntry = this->Derivatives_of_StateBeforeEvent[dIndex];

                        size_t stateIndex = std::get<1>(derivativeEntry);

                        size_t Xindex = std::get<0>(derivativeEntry);

                        size_t Gindex = this->Gindices_distance_constraint[entryIndex];

                        double TheDerivative = (this->state_before_event(stateIndex) / distance)_GETVALUE * std::get<2>(this->Derivatives_of_StateBeforeEvent[dIndex]);

                        G[Gindex] += this->X_scale_factors->operator[](Xindex)
                            * TheDerivative
                            / this->myUniverse->LU;
                    }//end loop over derivative indices
                }//end derivatives
            }//end distance constraint

            //Step 3: launch mass constraint
            if (this->useLV)
            {

                this->myLaunchVehicle->computePerformance(C3, this->myOptions->LV_margin);

                doubleType launch_mass = this->myLaunchVehicle->getDeliveredMass();

                if (launch_mass > this->myJourneyOptions->maximum_mass)
                {
                    launch_mass = this->myJourneyOptions->maximum_mass;
                    this->dm_dvinf = 0.0;
                }
                else
                {
                    this->dm_dvinf = this->myLaunchVehicle->getdmdC3() * 2.0 * vinf _GETVALUE;
                }

                F[Findex++] = (launch_mass - this->state_before_event(6)) / this->myJourneyOptions->maximum_mass;

                if (needG)
                {
                    //wrt encoded mass
                    G[this->Gindex_LVmassConstraint_encodedMass] = -this->X_scale_factors->operator[](this->Xindex_encoded_state[6])
                        / this->myJourneyOptions->maximum_mass;
                    
                    //wrt v-infinity magnitude
                    G[this->Gindex_LVmassConstraint_vinfinity_magnitude] = this->X_scale_factors->operator[](this->Xindex_vinfinity_magnitude)
                        * this->dm_dvinf
                        / this->myJourneyOptions->maximum_mass;
                }
            }//end launch mass constraint

            //Step 4: mass increment
            this->process_mass_multipliers(X, Xindex, F, Findex, G, needG);
        }//process_event_left_side()
        
        void PeriapseLaunch::process_deltav_contribution(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //add the event delta-v
            this->EventDeterministicDeltav = X[this->Xindex_vinfinity_magnitude];

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
        void PeriapseLaunch::output(std::ofstream& outputfile,
            const double& launchdate,
            size_t& eventcount)
        {
            this->mySpacecraft->setActiveStage(this->stageIndex);

            std::string event_type = "launch";

            std::string boundary_name = "periapse";

            math::Matrix<doubleType> empty3vector(3, 1, 0.0);
            math::Matrix<doubleType> dV(3, 1, 0.0);

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
            double Isp = -1;
            
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
                0.0,
                0.0,
                -1.0,
                this->mySpacecraft->getAvailablePower(),
                0.0,
                0,
                0.0,
                "none");
        }

        void PeriapseLaunch::output_maneuver_and_target_spec(std::ofstream& maneuver_spec_file, std::ofstream& target_spec_file, bool& haveManeuverNeedTarget)
        { 
            //there is no maneuver in this boundary type
        }//end output_maneuver_and_target_spec()
    }//end namespace BoundaryEvents
}//end namespace EMTG