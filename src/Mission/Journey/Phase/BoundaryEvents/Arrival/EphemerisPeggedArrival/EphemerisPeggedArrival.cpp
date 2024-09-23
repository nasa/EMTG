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

#include "EphemerisPeggedArrival.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        EphemerisPeggedArrival::EphemerisPeggedArrival(const std::string& name,
                                                                              const size_t& journeyIndex,
                                                                              const size_t& phaseIndex,
                                                                              size_t& stageIndex,
                                                                              Astrodynamics::universe* Universe,
                                                                              HardwareModels::Spacecraft* mySpacecraft,
                                                                              missionoptions* myOptions)
        {
            this->initialize(name,
                            journeyIndex,
                            phaseIndex,
                            stageIndex,
                            Universe,
                            mySpacecraft,
                            myOptions);
        }//end constructor

        void EphemerisPeggedArrival::initialize(const std::string& name,
            const size_t& journeyIndex,
            const size_t& phaseIndex,
            size_t& stageIndex,
            Astrodynamics::universe* Universe,
            HardwareModels::Spacecraft* mySpacecraft,
            missionoptions* myOptions)
        {
            this->EphemerisPeggedBoundary::initialize(name,
                journeyIndex,
                phaseIndex,
                stageIndex,
                Universe,
                mySpacecraft,
                myOptions);

            this->ArrivalEvent::initialize(name,
                journeyIndex,
                phaseIndex,
                stageIndex,
                Universe,
                mySpacecraft,
                myOptions);

            this->myBody = &myUniverse->bodies[this->myOptions->Journeys[this->journeyIndex].sequence[phaseIndex + 1] - 1];

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
        void EphemerisPeggedArrival::calcbounds_event_left_side(std::vector<size_t> timeVariables)
        {
            this->X_index_of_first_decision_variable_in_this_event = this->Xdescriptions->size();

            //Step 1: mass and epoch bounds
            //mass bounds
            std::vector<double> MassBounds({1.0e-13, this->myJourneyOptions->maximum_mass});

            //epoch bounds
            std::vector<double> EpochBounds(2);
            if (this->isLastEventInJourney
                && this->myJourneyOptions->timebounded == 2)//bounded arrival date
            {
                EpochBounds[0] = this->myJourneyOptions->arrival_date_bounds[0];
                EpochBounds[1] = this->myJourneyOptions->arrival_date_bounds[1];
            }
            else
            {
                EpochBounds[0] = this->myOptions->launch_window_open_date;
                EpochBounds[1] = this->myBody->getEphemerisWindowClose();
            }

            //Step 2: base ephemeris pegged boundary
            EphemerisPeggedBoundary::calcbounds_event_left_side(MassBounds, EpochBounds, timeVariables);
        }//end calcbounds_event_left_side()

        void EphemerisPeggedArrival::calcbounds_event_right_side()
        {
            //base classes
            EphemerisPeggedBoundary::calcbounds_event_right_side();

            ArrivalEvent::calcbounds_event_right_side();
        }//end calcbounds_event_right_side

        //**************************************process functions
        void EphemerisPeggedArrival::
            process_event_left_side(const std::vector<doubleType>& X,
                                    size_t& Xindex,
                                    std::vector<doubleType>& F,
                                    size_t& Findex,
                                    std::vector<double>& G,
                                    const bool& needG)
        {
            //base ephemeris pegged boundary
            EphemerisPeggedBoundary::process_event_left_side(X, Xindex, F, Findex, G, needG);

            //save the raw ephemeris state, i.e. before any maneuvering happens
            //we are cheating here, taking advantage of the fact that we "know" ephemeris pegged arrivals have no time width
            this->state_after_event_raw = this->state_before_event;
        }//end process_event_left_side()

        

        void EphemerisPeggedArrival::
            process_event_right_side(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG)
        {
            //Step 1: process the state
            this->state_after_event = this->state_before_event;
            this->state_after_event(6) -= this->chemical_fuel_used + this->chemical_oxidizer_used + this->electric_propellant_used;

            //Step 2: copy the derivative entries, but recognize that we don't want to kill off the derivative entries that are special to StateAfterEvent and not part of StateBefore Event
            //recognize that some multiplicative things may have happened to some states
            for (size_t dIndex = 0; dIndex < this->Derivatives_of_StateBeforeEvent.size(); ++dIndex)
            {
                this->Derivatives_of_StateAfterEvent[dIndex] = this->Derivatives_of_StateBeforeEvent[dIndex];
                
                size_t stateIndex = std::get<1>(this->Derivatives_of_StateAfterEvent[dIndex]);
                std::get<2>(this->Derivatives_of_StateAfterEvent[dIndex]) *= this->ETM(stateIndex, stateIndex);
            }

            for (size_t dIndex = 0; dIndex < this->Derivatives_of_StateBeforeEvent_wrt_Time.size(); ++dIndex)
                this->Derivatives_of_StateAfterEvent_wrt_Time[dIndex] = this->Derivatives_of_StateBeforeEvent_wrt_Time[dIndex];

            this->EventRightEpoch = this->state_after_event(7);

            //Step 3: perform arrival mass increment
            this->process_post_arrival_mass_increment();

            //Step 4: perform post-arrival delta-v
            this->process_post_arrival_deltav();

            //Step 5: adjust the mass derivative as necessary
            std::get<2>(this->Derivatives_of_StateAfterEvent[this->dIndex_mass_wrt_encodedMass]) *= this->ETM(6, 6);
        }//end process_event_right_side()

        //******************************************output methods
        void EphemerisPeggedArrival::output(std::ofstream& outputfile,
            const double& launchdate,
            size_t& eventcount)
        {            //output the end of the mission
            if (this->isLastEventInMission && this->myOptions->output_dormant_journeys)
            {
                std::string event_type = "waiting";

                std::string boundary_name = this->myBody->name;

                math::Matrix<doubleType> empty3vector(3, 1, 0.0);

                math::Matrix<doubleType> waitState(8, 1, 0.0);
                doubleType waitEpoch;

                for (size_t step = 1; step < this->myOptions->num_timesteps + 1; ++step)
                {
                    waitEpoch = this->EventLeftEpoch + this->myOptions->post_mission_wait_time * 86400.0 * ((double)(step) / this->myOptions->num_timesteps);

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
                }//end time loop
            }//end output dormant journeys
        }//end output()

        void EphemerisPeggedArrival::output_maneuver_and_target_spec(std::ofstream& maneuver_spec_file, std::ofstream& target_spec_file, bool& haveManeuverNeedTarget)
        {
            //The base-class version of this method assumes that there is NOT an arrival maneuver,
            //and therefore does not write any lines in the maneuver spec

            //only do anything if this is the last event in the mission
            if (this->isLastEventInMission
                && haveManeuverNeedTarget)
            {
                //reset the flag that tells us we need a target spec line
                haveManeuverNeedTarget = false;

                //Step 1: initialize a target spec object
                target_spec_line myTargetSpecLine(this->name,
                    "EME2000",
                    this->myUniverse->central_body.name,
                    this->state_after_event(7),
                    this->state_after_event);

                //Step 2: write target spec object
                myTargetSpecLine.write(target_spec_file);
            }
        }//end output_maneuver_and_target_spec()
    }//end namespace BoundaryEvents
}//end namespace EMTG