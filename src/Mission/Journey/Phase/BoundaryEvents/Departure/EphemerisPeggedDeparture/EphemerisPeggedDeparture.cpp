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

#include "EphemerisPeggedDeparture.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        EphemerisPeggedDeparture::EphemerisPeggedDeparture(const std::string& name,
                                                           const size_t& journeyIndex,
                                                           const size_t& phaseIndex,
                                                           size_t& stageIndex,
                                                           Astrodynamics::universe* Universe,
                                                           HardwareModels::Spacecraft* mySpacecraft,
                                                           missionoptions* myOptions,
                                                           ArrivalEvent* PreviousPhaseArrivalEvent)
        {
            this->initialize(name,
                             journeyIndex,
                             phaseIndex,
                             stageIndex,
                             Universe,
                             mySpacecraft,
                             myOptions,
                             PreviousPhaseArrivalEvent);
        }//end constructor
        
        void EphemerisPeggedDeparture::initialize(const std::string& name,
            const size_t& journeyIndex,
            const size_t& phaseIndex,
            size_t& stageIndex,
            Astrodynamics::universe* Universe,
            HardwareModels::Spacecraft* mySpacecraft,
            missionoptions* myOptions,
            ArrivalEvent* PreviousPhaseArrivalEvent)
        {
            this->EphemerisPeggedBoundary::initialize(name,
                journeyIndex,
                phaseIndex,
                stageIndex,
                Universe,
                mySpacecraft,
                myOptions);

            this->DepartureEvent::initialize(name,
                journeyIndex,
                phaseIndex,
                stageIndex,
                Universe,
                mySpacecraft,
                myOptions,
                PreviousPhaseArrivalEvent);

            if (this->myOptions->Journeys[journeyIndex].sequence[phaseIndex] == 0)
            {
                throw std::invalid_argument(this->name + " body is set to Universe central body. You can't currently design a mission to an ephemeris-pegged central body destination. That is effectively the same thing as requesting that the spacecraft fly to the origin of the universe. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }
            else
                this->myBody = &myUniverse->bodies[this->myOptions->Journeys[journeyIndex].sequence[phaseIndex] - 1];


            //If the parent Journey uses the integrated propagator and this boundary event's body is turned on as a perturber, throw an error!
            //Otherwise there will be a singularity in the integration and the user will not got a valid answer.
            if (this->myOptions->perturb_thirdbody
                && ((!this->myJourneyOptions->override_PropagatorType && this->myOptions->propagatorType == PropagatorType::IntegratedPropagator)
                    || (this->myJourneyOptions->override_PropagatorType && this->myJourneyOptions->propagatorType == PropagatorType::IntegratedPropagator)
                    || (this->myJourneyOptions->phase_type == PhaseType::PSFB || this->myJourneyOptions->phase_type == PhaseType::FBLT || this->myJourneyOptions->phase_type == PhaseType::FBLTS)))
            {
                for (size_t bodyIndex : this->myJourneyOptions->perturbation_bodies)
                {
                    if (this->myUniverse->bodies[bodyIndex - 1].spice_ID == this->myBody->spice_ID)
                        throw std::invalid_argument(this->name + " is an ephemeris pegged boundary event that occurs at " + this->myBody->name + ". Unfortunately, you have also set this journey to integrated propagation and have chosen " + this->myBody->name + " as a third-body perturber. This would cause a singularity in the integrator and so EMTG will not let you do it. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
                }
            }
       }//end initialize()
        
        //******************************************calcbounds methods
        void EphemerisPeggedDeparture::calcbounds_event_left_side(std::vector<size_t> timeVariables)
        {
            this->X_index_of_first_decision_variable_in_this_event = this->Xdescriptions->size();

            //Step 1: mass and epoch bounds
            //mass bounds
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

            //epoch bounds
            std::vector<double> EpochBounds(2);
            if (this->myJourneyOptions->bounded_departure_date)//bounded arrival date
            {
                EpochBounds[0] = this->myJourneyOptions->departure_date_bounds[0];
                EpochBounds[1] = this->myJourneyOptions->departure_date_bounds[1];
            }
            else if (this->isFirstEventInJourney)
            {
                double epoch_max = fmin(this->myOptions->launch_window_open_date + this->myJourneyOptions->wait_time_bounds[1], this->myBody->getEphemerisWindowClose());

                EpochBounds[0] = this->myOptions->launch_window_open_date + this->myJourneyOptions->wait_time_bounds[0];
                EpochBounds[1] = epoch_max;
            }
            else
            {
                EpochBounds[0] = this->myOptions->launch_window_open_date;
                EpochBounds[1] = this->myBody->getEphemerisWindowClose();
            }

            //Step 2: base departure class - creates the wait time if applicable
            DepartureEvent::calcbounds_event_left_side();

            //Step 3: base ephemeris pegged boundary
            EphemerisPeggedBoundary::calcbounds_event_left_side(MassBounds, EpochBounds, timeVariables);

            //Step 4: mass continuity
            if (this->hasWaitTime)
                this->calcbounds_left_mass_continuity_constraint();

            //Step 5: mass multipliers
            this->calcbounds_mass_multipliers();
        }//end calcbounds_event_left_side()


        void EphemerisPeggedDeparture::calcbounds_event_right_side()
        {
            //base class
            EphemerisPeggedBoundary::calcbounds_event_right_side();
        }//end calcbounds_event_right_side

        void EphemerisPeggedDeparture::calcbounds_specialized_constraints()
        {
            this->DepartureEvent::calcbounds_specialized_constraints();
        }

        //**************************************process functions
        void EphemerisPeggedDeparture::
            process_event_left_side(const std::vector<doubleType>& X,
                                    size_t& Xindex,
                                    std::vector<doubleType>& F,
                                    size_t& Findex,
                                    std::vector<double>& G,
                                    const bool& needG)
        {
            //base departure class
            this->DepartureEvent::process_event_left_side(X, Xindex, F, Findex, G, needG);

            //base ephemeris pegged boundary
            this->EphemerisPeggedBoundary::process_event_left_side(X, Xindex, F, Findex, G, needG);

            //mass continuity
            if (this->hasWaitTime)
                this->process_left_mass_continuity_constraint(X, Xindex, F, Findex, G, needG);

            //mass increment
            this->process_mass_multipliers(X, Xindex, F, Findex, G, needG);
        }//end process_event_left_side()

        

        void EphemerisPeggedDeparture::
            process_event_right_side(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG)
        {
            //base ephemeris pegged boundary
            EphemerisPeggedBoundary::process_event_right_side(X, Xindex, F, Findex, G, needG);
        }//end process_event_right_side()

         //******************************************output methods
        void EphemerisPeggedDeparture::output(std::ofstream& outputfile,
            const double& launchdate,
            size_t& eventcount)
        {
            this->mySpacecraft->setActiveStage(this->stageIndex);

            if (this->myOptions->output_dormant_journeys && this->hasWaitTime)
            {
                std::string event_type = "waiting";

                std::string boundary_name = this->myBody->name;

                math::Matrix<doubleType> empty3vector(3, 1, 0.0);

                math::Matrix<doubleType> waitState(8, 1, 0.0);
                doubleType waitEpoch;

                for (size_t step = 0; step < this->myOptions->num_timesteps; ++step)
                {
                    waitEpoch = this->EventLeftEpoch - this->EventWaitTime * ((double)(this->myOptions->num_timesteps - step) / this->myOptions->num_timesteps);

                    doubleType body_state_and_derivatives[12];
                    this->myBody->locate_body(waitEpoch,
                        body_state_and_derivatives,
                        false,
                        *this->myOptions);

                    //position/velocity
                    for (size_t stateIndex = 0; stateIndex < 6; ++stateIndex)
                        waitState(stateIndex) = body_state_and_derivatives[stateIndex];
                        
                    //mass
                    waitState(6) = this->state_before_event(6);

                    //epoch
                    waitState(7) = waitEpoch;


                    //where is the Sun?
                    math::Matrix<doubleType> R_sc_Sun(3, 1, 0.0);
                    if (this->myUniverse->central_body_SPICE_ID == 10)
                    {
                        R_sc_Sun = waitState.getSubMatrix1D(0, 2);
                    }
                    else
                    {
                        //where is the central body relative to the sun?
                        doubleType central_body_state_and_derivatives[12];
                        this->myUniverse->locate_central_body(waitState(7),
                            central_body_state_and_derivatives,
                            *this->myOptions,
                            false);
                        math::Matrix<doubleType> R_CB_Sun(3, 1, 0.0);
                        for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
                        {
                            R_CB_Sun(stateIndex) = central_body_state_and_derivatives[stateIndex];
                        }
                        R_sc_Sun = waitState.getSubMatrix1D(0, 2) + R_CB_Sun;
                    }
                    this->mySpacecraft->computePowerState(R_sc_Sun.getSubMatrix1D(0, 2).norm() / this->myOptions->AU, waitState(7));

                    write_output_line(outputfile,
                        eventcount,
                        event_type,
                        boundary_name,
                        this->EventTimeWidth,
                        0.0,
                        0.0,
                        0.0,
                        0.0,
                        0.0,
                        0.0,
                        waitState,
                        empty3vector,
                        empty3vector,
                        0.0,
                        0.0,
                        0.0,
                        this->mySpacecraft->getAvailablePower(),
                        0.0,
                        0,
                        0.0,
                        "none");
                }
            }
        }//end output()
        
        void EphemerisPeggedDeparture::output_ephemeris(std::ofstream& outputfile)
        {
            if (this->isFirstEventInMission)
                this->BoundaryEventBase::output_ephemeris(outputfile);
            else if (this->hasWaitTime)
            {
                //Step 0: we'll need an output vector
                math::Matrix<doubleType> output_state(8, 1, 0.0);

                ////Step 1: output the wait time - we do this by looking up the position of the body and printing it
                double EphemerisOutputResolution = 2.0 * math::PI * sqrt(this->myUniverse->central_body.radius * this->myUniverse->central_body.radius * this->myUniverse->central_body.radius / this->myUniverse->central_body.mu);

                double lookBackTime = this->EventWaitTime _GETVALUE - EphemerisOutputResolution;
                while (lookBackTime > EphemerisOutputResolution)
                {
                    //calculate the epoch of interest
                    output_state(6) = this->state_before_event(6);
                    output_state(7) = this->state_before_event(7) - lookBackTime;

                    //look up the state of the body on that epoch
                    doubleType body_state[12];

                    this->myBody->locate_body(output_state(7),
                        body_state,
                        false, 
                        *this->myOptions);

                    for (size_t stateIndex = 0; stateIndex < 6; ++stateIndex)
                        output_state(stateIndex) = body_state[stateIndex];

                    //make sure the state is in the sun-centered J2000 Earth equatorial frame
                    if (this->myUniverse->central_body.spice_ID != this->myOptions->forward_integrated_ephemeris_central_body_SPICE_ID)
                    {
                        double LT_dump;
                        double bodyStateDouble[6];
                        spkez_c(this->myUniverse->central_body.spice_ID, output_state(7)_GETVALUE - (51544.5 * 86400.0), "J2000", "NONE", this->myOptions->forward_integrated_ephemeris_central_body_SPICE_ID, bodyStateDouble, &LT_dump);

                        for (size_t stateIndex = 0; stateIndex < 6; ++stateIndex)
                            output_state(stateIndex) += bodyStateDouble[stateIndex];
                    }

                    //print the state
                    this->write_ephemeris_line(outputfile,
                        output_state,
                        math::Matrix<doubleType>(3, 1, 0.0),//control vector
                        0.0,
                        0.0,
                        0.0,
                        0,
                        0.0,
                        "none");

                    //update the look back time
                    lookBackTime -= (this->EventWaitTime _GETVALUE - lookBackTime) > EphemerisOutputResolution
                        ? EphemerisOutputResolution
                        : (this->EventWaitTime _GETVALUE - lookBackTime);
                }
                //now output the actual departure state
                this->BoundaryEventBase::output_ephemeris(outputfile);
            }//end while loop over ephemeris lookback
        }//end output_ephemeris()

        void EphemerisPeggedDeparture::output_maneuver_and_target_spec(std::ofstream& maneuver_spec_file, std::ofstream& target_spec_file, bool& haveManeuverNeedTarget)
        {
            //The base-class version of this method assumes that there is NOT an arrival maneuver,
            //and therefore does not write any lines in the maneuver spec
            //it also assumes bplane target values of 0.0
            if (haveManeuverNeedTarget)
            {
                //reset the flag that tells us we need a target spec line
                haveManeuverNeedTarget = false;

                //Step 1: initialize a target spec object
                target_spec_line myTargetSpecLine(this->name,
                    "EME2000",
                    this->myUniverse->central_body.name,
                    this->state_before_event(7),
                    this->state_before_event);

                //Step 2: write target spec object
                myTargetSpecLine.write(target_spec_file);
            }
        }//end output_maneuver_and_target_spec()
    }//end namespace BoundaryEvents
}//end namespace EMTG