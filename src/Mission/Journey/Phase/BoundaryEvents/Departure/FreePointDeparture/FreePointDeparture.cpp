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

#include "FreePointDeparture.h"
#include "StateRepresentationFactory.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        FreePointDeparture::FreePointDeparture(const std::string& name,
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

        void FreePointDeparture::initialize(const std::string& name,
                                const size_t& journeyIndex,
                                const size_t& phaseIndex,
                                size_t& stageIndex,
                                Astrodynamics::universe* Universe,
                                HardwareModels::Spacecraft* mySpacecraft,
                                missionoptions* myOptions,
                                ArrivalEvent* PreviousPhaseArrivalEvent)
        {
            //we have to do this up front, before delegating
            this->AllowStateToPropagate = myOptions->Journeys[journeyIndex].AllowJourneyFreePointDepartureToPropagate;

            //if we are allowing the state to propagate and this is a later journey, then we need a wait time
            if (this->AllowStateToPropagate && journeyIndex > 0)
                this->hasWaitTime = true;

            this->FreePointBoundary::initialize(name,
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

            this->myStateRepresentationEnum = this->myJourneyOptions->departure_elements_state_representation;

            //create the state representation
            this->myStateRepresentation = Astrodynamics::CreateStateRepresentation(this->myStateRepresentationEnum, this->myUniverse->mu);
                                             
            this->myEncodedReferenceFrame = this->myJourneyOptions->departure_elements_frame;

            this->ReferenceEpoch = this->myJourneyOptions->departure_elements_reference_epoch;

            //are we using an object-referenced frame? If so, let's make a body
            if (this->myEncodedReferenceFrame == ReferenceFrame::ObjectReferenced)
            {
                if (this->myOptions->Journeys[journeyIndex].destination_list[0] < 1)
                {
                    throw std::invalid_argument(this->name + " reference body must be set to a body in the universe. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
                }
                else
                    this->myBody = &myUniverse->bodies[this->myOptions->Journeys[journeyIndex].destination_list[0] - 1];
            }
        }//end initialize()
        
        //******************************************calcbounds methods
        void FreePointDeparture::calcbounds_event_left_side(std::vector<size_t> timeVariables)
        {
            this->X_index_of_first_decision_variable_in_this_event = this->Xdescriptions->size();

            std::vector< std::tuple<double, double> > StateBounds;
            //Step 1: state bounds
            //Step 1.1: position and velocity
            for (size_t stateIndex = 0; stateIndex < 6; ++stateIndex)
            {
                if (this->myJourneyOptions->departure_elements_vary_flag[stateIndex])
                    StateBounds.push_back({ this->myJourneyOptions->departure_elements_bounds[2 * stateIndex], this->myJourneyOptions->departure_elements_bounds[2 * stateIndex + 1] });
                else
                    StateBounds.push_back({ this->myJourneyOptions->departure_elements[stateIndex] - 1.0e-13, this->myJourneyOptions->departure_elements[stateIndex] + 1.0e-13});
            }

            //Step 1.2: mass 
            if (this->isFirstEventInMission && this->hasFixedInitialMass)
            {
                StateBounds.push_back({ this->myJourneyOptions->maximum_mass - 1.0e-13, this->myJourneyOptions->maximum_mass });
            }
            else
            {
                StateBounds.push_back({ 1.0e-13, this->myJourneyOptions->maximum_mass });
            }

            //Step 1.3: epoch
            if (this->isFirstEventInMission)
            {
                StateBounds.push_back({ this->myOptions->launch_window_open_date + this->myOptions->Journeys.front().wait_time_bounds[0],
                                        this->myOptions->launch_window_open_date + this->myOptions->Journeys.front().wait_time_bounds[1] });
            }

            //Step 2: base departure class
            DepartureEvent::calcbounds_event_left_side();

            //Step 3: base free point boundary
            if (this->isFirstEventInMission)
                FreePointBoundary::calcbounds_event_left_side(StateBounds, timeVariables);
            else
            {//if this is not the first event in the mission, then we need to extract the state from the previous event - all except mass
                //Step 3.1: set the current stage
                this->mySpacecraft->setActiveStage(this->stageIndex);

                //Step 3.2: get the derivatives of the state after the previous event
                std::vector< std::tuple<size_t, size_t, double> >& PreviousEvent_Derivatives_of_StateBeforeEvent
                    = this->PreviousPhaseArrivalEvent->get_Derivatives_of_StateAfterEvent();


                std::vector< std::tuple<size_t, size_t, double> >& PreviousEvent_Derivatives_of_StateBeforeEvent_wrt_Time
                    = this->PreviousPhaseArrivalEvent->get_Derivatives_of_StateAfterEvent_wrt_Time();

                //Step 3.3: assemble the derivative skeleton
                for (size_t stateIndex = 0; stateIndex < 6; ++stateIndex)
                {
                    //non-time dependencies
                    std::vector<size_t> state_dIndex_StateAfterPreviousEvent_wrt_DecisionVariables;
                    std::vector<size_t> state_dIndex_StateBeforeEvent_wrt_DecisionVariables;

                    for (size_t dIndex = 0; dIndex < PreviousEvent_Derivatives_of_StateBeforeEvent.size(); ++dIndex)
                        if (std::get<1>(PreviousEvent_Derivatives_of_StateBeforeEvent[dIndex]) == stateIndex)
                        {
                            state_dIndex_StateAfterPreviousEvent_wrt_DecisionVariables.push_back(dIndex);
                            this->Derivatives_of_StateBeforeEvent.push_back(PreviousEvent_Derivatives_of_StateBeforeEvent[dIndex]);
                            state_dIndex_StateBeforeEvent_wrt_DecisionVariables.push_back(this->Derivatives_of_StateBeforeEvent.size() - 1);
                        }

                    this->dIndex_StateAfterPreviousEvent_wrt_DecisionVariables.push_back(state_dIndex_StateAfterPreviousEvent_wrt_DecisionVariables);
                    this->dIndex_StateBeforeEvent_wrt_DecisionVariables.push_back(state_dIndex_StateBeforeEvent_wrt_DecisionVariables);

                    //same thing, with respect to time
                    std::vector<size_t> state_dIndex_StateAfterPreviousEvent_wrt_DecisionVariables_wrt_Time;
                    std::vector<size_t> state_dIndex_StateBeforeEvent_wrt_DecisionVariables_wrt_Time;

                    for (size_t dIndex = 0; dIndex < PreviousEvent_Derivatives_of_StateBeforeEvent_wrt_Time.size(); ++dIndex)
                    {
                        if (std::get<1>(PreviousEvent_Derivatives_of_StateBeforeEvent_wrt_Time[dIndex]) == stateIndex)
                        {
                            state_dIndex_StateAfterPreviousEvent_wrt_DecisionVariables_wrt_Time.push_back(dIndex);
                            this->Derivatives_of_StateBeforeEvent_wrt_Time.push_back(PreviousEvent_Derivatives_of_StateBeforeEvent_wrt_Time[dIndex]);
                            state_dIndex_StateBeforeEvent_wrt_DecisionVariables_wrt_Time.push_back(this->Derivatives_of_StateBeforeEvent_wrt_Time.size() - 1);
                        }
                    }

                    this->dIndex_StateAfterPreviousEvent_wrt_DecisionVariables_wrt_Time.push_back(state_dIndex_StateAfterPreviousEvent_wrt_DecisionVariables_wrt_Time);
                    this->dIndex_StateBeforeEvent_wrt_DecisionVariables_wrt_Time.push_back(state_dIndex_StateBeforeEvent_wrt_DecisionVariables_wrt_Time);
                }//end loop over states


                //Step 3.4: encode a mass variable
                this->Xlowerbounds->push_back(1.0e-13);
                this->Xupperbounds->push_back(this->myJourneyOptions->maximum_mass);
                this->X_scale_factors->push_back(1.0 / this->myUniverse->continuity_constraint_scale_factors(6));
                this->Xdescriptions->push_back(prefix + "event left state mass");
                size_t Xindex_mass = this->Xdescriptions->size() - 1;

                //if we propagate with an integrator, everything has a derivative wrt mass
                if (this->myPropagatorType == PropagatorType::IntegratedPropagator && this->AllowStateToPropagate)
                {
                    for (size_t stateIndex : {0, 1, 2, 3, 4, 5})
                    {
                        this->Derivatives_of_StateBeforeEvent.push_back(std::make_tuple(Xindex_mass, stateIndex, 1.0));
                        this->dIndex_StateBeforeEvent_wrt_encoded_mass.push_back(this->Derivatives_of_StateBeforeEvent.size() - 1);
                    }
                }
                this->Derivatives_of_StateBeforeEvent.push_back(std::make_tuple(Xindex_mass, 6, 1.0));
                this->dIndex_mass_wrt_encodedMass = this->Derivatives_of_StateBeforeEvent.size() - 1;
                this->Xindex_encoded_state.push_back(Xindex_mass);

                //Step 3.5:  left mass continuity constraint
                this->calcbounds_left_mass_continuity_constraint();

                //Step 3.6: epoch
                this->calculate_dependencies_left_epoch(timeVariables);
                                
                //Step 3.7: if we allow the initial state to propagate then we have an additional set of derivatives of all of the states with respect to wait time
                if (this->AllowStateToPropagate && !this->isFirstEventInMission)
                {
                    for (size_t stateIndex = 0; stateIndex < 6; ++stateIndex)
                    {
                        this->Derivatives_of_StateBeforeEvent_wrt_Time.push_back({ this->Xindices_EventLeftEpoch.back(), stateIndex, 1.0 });
                        this->dIndex_StateBeforeEvent_wrt_DecisionVariables_wrt_Time[stateIndex].push_back(this->Derivatives_of_StateBeforeEvent_wrt_Time.size() - 1);
                    }
                }
            }//end successive departure
                       
            //mass multipliers
            this->calcbounds_mass_multipliers();
        }//end calcbounds_event_left_side()


        void FreePointDeparture::calcbounds_event_right_side()
        {
            //base class
            FreePointBoundary::calcbounds_event_right_side();
        }//end calcbounds_event_right_side

        void FreePointDeparture::calcbounds_specialized_constraints()
        {
            this->DepartureEvent::calcbounds_specialized_constraints();
        }

        //**************************************process functions
        void FreePointDeparture::
            process_event_left_side(const std::vector<doubleType>& X,
                                    size_t& Xindex,
                                    std::vector<doubleType>& F,
                                    size_t& Findex,
                                    std::vector<double>& G,
                                    const bool& needG)
        {
            //Step 1: base departure class
            this->DepartureEvent::process_event_left_side(X, Xindex, F, Findex, G, needG);

            //Step 2: base free point boundary
            if (this->isFirstEventInMission)
                this->FreePointBoundary::process_event_left_side(X, Xindex, F, Findex, G, needG);
            else
            {//extract mass from decision vector and copy everything else from the previous phase arrival event
                //Step 2.1: set the active stage
                this->mySpacecraft->setActiveStage(this->stageIndex);

                //Step 2.2: get the position and velocity, and their derivatives, from the previous event
                math::Matrix<doubleType>& PreviousEvent_StateAfterEvent
                    = this->PreviousPhaseArrivalEvent->get_state_after_event();

                for (size_t stateIndex = 0; stateIndex < 6; ++stateIndex)
                    this->state_before_event(stateIndex) = PreviousEvent_StateAfterEvent(stateIndex);


                //Step 2.3: process the left epoch
                this->process_left_epoch(X, Xindex, F, Findex, G, needG);

                //Step 2.4: extract mass from the decision vector
                this->state_before_event(6) = X[Xindex++];

                //Step 2.5: propagate, if appropriate
                doubleType PropagationTime = this->EventWaitTime + 1.0e-13;
                if (this->AllowStateToPropagate)
                {
                    this->total_number_of_states_to_integrate = needG
                        ? 10 + 13 * 13
                        : 10;
                    //Step 2.5.1: copy the state to a temporary vector
                    this->StateBeforeEventBeforePropagation.shallow_copy(this->state_before_event, 8);
                    this->StateBeforeEventBeforePropagation(7) = 0.0;                    
                    //epoch needs to *not* include wait time
                    for (size_t listIndex = 0; listIndex < this->Xindices_EventLeftEpoch.size() - 1; ++listIndex)
                    {
                        this->StateBeforeEventBeforePropagation(7) += X[this->Xindices_EventLeftEpoch[listIndex]];
                    }
                    //Step 2.5.2: then we need to propagate
                    this->dPropagatedStatedIndependentVariable.assign_zeros();
                    this->myPropagator->setCurrentEpoch(this->StateBeforeEventBeforePropagation(7));
                    this->myPropagator->setIndexOfEpochInStateVec(7);
                    this->myPropagator->setCurrentIndependentVariable(this->StateBeforeEventBeforePropagation(7));
                    this->myPropagator->propagate(PropagationTime, needG);
                    this->state_before_event(7) = this->EventLeftEpoch;
                    if (this->myPropagatorType == PropagatorType::IntegratedPropagator)
                    {
                        this->state_before_event(8) = 0.0; //chemical fuel - we're not actually moving the spacecraft, we're redefining a boundary point. So tankage does not change.
                        this->state_before_event(9) = 0.0; //chemical oxidizer
                    }
                }

                //Step 2.6 derivatives of position and velocity
                if (needG)
                {
                    //get the derivatives of the state after the previous event
                    std::vector< std::tuple<size_t, size_t, double> >& PreviousEvent_Derivatives_of_StateAfterEvent
                        = this->PreviousPhaseArrivalEvent->get_Derivatives_of_StateAfterEvent();
                    std::vector< std::tuple<size_t, size_t, double> >& PreviousEvent_Derivatives_of_StateAfterEvent_wrt_Time
                        = this->PreviousPhaseArrivalEvent->get_Derivatives_of_StateAfterEvent_wrt_Time();


                    if (this->AllowStateToPropagate)
                    {
                        //If the state was allowed to propagate, then we have to multiply the derivatives of the states w.r.t. non-time decision variables by STM entries
                        for (size_t departureStateIndex = 0; departureStateIndex < 6; ++departureStateIndex)
                        {
                            //with respect to variables affecting previous phase arrival event
                            for (size_t varIndex = 0; varIndex < this->dIndex_StateBeforeEvent_wrt_DecisionVariables[departureStateIndex].size(); ++varIndex)
                            {
                                size_t dIndex = this->dIndex_StateBeforeEvent_wrt_DecisionVariables[departureStateIndex][varIndex];
                                std::get<2>(this->Derivatives_of_StateBeforeEvent[dIndex]) = 0.0;
                                size_t Xindex = std::get<0>(this->Derivatives_of_StateBeforeEvent[dIndex]);

                                //with respect to decision variables affecting the previous event's arrival state
                                for (size_t arrivalStateIndex = 0; arrivalStateIndex < 6; ++arrivalStateIndex)
                                {
                                    for (size_t arrivalVarIndex = 0; arrivalVarIndex < this->dIndex_StateAfterPreviousEvent_wrt_DecisionVariables[arrivalStateIndex].size(); ++arrivalVarIndex)
                                    {
                                        size_t dIndex_previousEvent = this->dIndex_StateAfterPreviousEvent_wrt_DecisionVariables[arrivalStateIndex][arrivalVarIndex];
                                        if (Xindex == std::get<0>(PreviousEvent_Derivatives_of_StateAfterEvent[dIndex_previousEvent]))
                                        {
                                            std::get<2>(this->Derivatives_of_StateBeforeEvent[dIndex]) += this->STM(departureStateIndex, arrivalStateIndex) * std::get<2>(PreviousEvent_Derivatives_of_StateAfterEvent[dIndex_previousEvent]);
                                        }
                                    }
                                }
                            }

                            //with respect to current encoded mass
                            if (this->myPropagatorType == PropagatorType::IntegratedPropagator)
                            {
                                size_t dIndex = this->dIndex_StateBeforeEvent_wrt_encoded_mass[departureStateIndex];

                                std::get<2>(this->Derivatives_of_StateBeforeEvent[dIndex]) = this->STM(departureStateIndex, 6);
                            }

                            //then we have to handle time derivatives
                            for (size_t varIndex = 0; varIndex < this->dIndex_StateAfterPreviousEvent_wrt_DecisionVariables_wrt_Time[departureStateIndex].size(); ++varIndex)
                            {
                                size_t dIndex = this->dIndex_StateBeforeEvent_wrt_DecisionVariables_wrt_Time[departureStateIndex][varIndex];
                                std::get<2>(this->Derivatives_of_StateBeforeEvent_wrt_Time[dIndex]) = 0.0;
                                size_t Xindex = std::get<0>(this->Derivatives_of_StateBeforeEvent_wrt_Time[dIndex]);

                                //with respect to decision variables affecting the previous event's arrival state
                                for (size_t arrivalStateIndex = 0; arrivalStateIndex < 6; ++arrivalStateIndex)
                                {
                                    for (size_t arrivalVarIndex = 0; arrivalVarIndex < this->dIndex_StateAfterPreviousEvent_wrt_DecisionVariables_wrt_Time[arrivalStateIndex].size(); ++arrivalVarIndex)
                                    {
                                        size_t dIndex_previousEvent = this->dIndex_StateAfterPreviousEvent_wrt_DecisionVariables_wrt_Time[arrivalStateIndex][arrivalVarIndex];
                                        if (Xindex == std::get<0>(PreviousEvent_Derivatives_of_StateAfterEvent[dIndex_previousEvent]))
                                        {
                                            std::get<2>(this->Derivatives_of_StateBeforeEvent_wrt_Time[dIndex]) += this->STM(departureStateIndex, arrivalStateIndex) * std::get<2>(PreviousEvent_Derivatives_of_StateAfterEvent_wrt_Time[dIndex_previousEvent]);
                                        }
                                    }
                                }
                            }

                            //wait time
                            if (this->myPropagatorType == PropagatorType::KeplerPropagator)
                                std::get<2>(this->Derivatives_of_StateBeforeEvent_wrt_Time[this->dIndex_StateBeforeEvent_wrt_DecisionVariables_wrt_Time[departureStateIndex].back()]) = (PropagationTime >= 0.0 ? 1.0 : -1.0) * this->dPropagatedStatedIndependentVariable(departureStateIndex);
                            else
                                std::get<2>(this->Derivatives_of_StateBeforeEvent_wrt_Time[this->dIndex_StateBeforeEvent_wrt_DecisionVariables_wrt_Time[departureStateIndex].back()]) = this->STM(departureStateIndex, 13);
                        }//end loop over current stateIndex
                    }
                    else
                    {
                        for (size_t stateIndex = 0; stateIndex < 6; ++stateIndex)
                        {
                            //derivatives with respect to non-time decision variables
                            for (size_t varIndex = 0; varIndex < this->dIndex_StateBeforeEvent_wrt_DecisionVariables[stateIndex].size(); ++varIndex)
                            {
                                size_t dIndex_previousEvent = this->dIndex_StateAfterPreviousEvent_wrt_DecisionVariables[stateIndex][varIndex];
                                size_t dIndex = this->dIndex_StateBeforeEvent_wrt_DecisionVariables[stateIndex][varIndex];
                                std::get<2>(this->Derivatives_of_StateBeforeEvent[dIndex]) = std::get<2>(PreviousEvent_Derivatives_of_StateAfterEvent[dIndex_previousEvent]);
                            }

                            //derivatives with respect to time
                            for (size_t varIndex = 0; varIndex < this->dIndex_StateBeforeEvent_wrt_DecisionVariables_wrt_Time[stateIndex].size(); ++varIndex)
                            {
                                size_t dIndex_previousEvent = this->dIndex_StateAfterPreviousEvent_wrt_DecisionVariables_wrt_Time[stateIndex][varIndex];
                                size_t dIndex = this->dIndex_StateBeforeEvent_wrt_DecisionVariables_wrt_Time[stateIndex][varIndex];
                                std::get<2>(this->Derivatives_of_StateBeforeEvent_wrt_Time[dIndex]) = std::get<2>(PreviousEvent_Derivatives_of_StateAfterEvent_wrt_Time[dIndex_previousEvent]);
                            }
                        }//end loop over states
                    }
                }//end derivatives

                //Step 2.7: left mass continuity constraint
                 this->process_left_mass_continuity_constraint(X, Xindex, F, Findex, G, needG);

                //Step 2.8: mass increment
                this->process_mass_multipliers(X, Xindex, F, Findex, G, needG);

            }//end successive journey free point departure
        }//end process_event_left_side()

        

        void FreePointDeparture::
            process_event_right_side(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG)
        {
            //base ephemeris pegged boundary
            FreePointBoundary::process_event_right_side(X, Xindex, F, Findex, G, needG);
        }//end process_event_right_side()

         //******************************************output methods
        void FreePointDeparture::output(std::ofstream& outputfile,
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

                    //TODO compute wait-state for non-body boundary conditions, involves propagation

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

        void FreePointDeparture::output_ephemeris(std::ofstream& outputfile)
        {
            if (this->isFirstEventInMission)
                this->BoundaryEventBase::output_ephemeris(outputfile);
            else if (this->AllowStateToPropagate)
            {
                //Step 0: we'll need an output vector
                math::Matrix<doubleType> output_state = this->StateBeforeEventBeforePropagation;

                //Step 1: set output resolution
                double EphemerisOutputResolution;
                if (this->myJourneyOptions->override_integration_step_size)
                    EphemerisOutputResolution = this->myJourneyOptions->integration_step_size;
                else
                    EphemerisOutputResolution = this->myOptions->integration_time_step_size;

                //Step 2: output the wait time - we'll do this by integrating forward in time and plotapussing
                //Step 2.1: temporarily assign the initial coast propagator to the output state
                this->myPropagator->setStateRight(output_state);
                this->myPropagator->setCurrentEpoch(this->StateBeforeEventBeforePropagation(7));
                this->myPropagator->setIndexOfEpochInStateVec(7);
                this->myPropagator->setCurrentIndependentVariable(this->StateBeforeEventBeforePropagation(7));

                //Step 2.2: propagate and print, skipping the first and last entry
                doubleType timeToPropagate = EphemerisOutputResolution;
                doubleType totalPropagationTime = this->EventWaitTime;
                while (timeToPropagate < totalPropagationTime)
                {
                    //Step 2.2.1: propagate
                    this->myPropagator->setCurrentEpoch(this->StateBeforeEventBeforePropagation(7));
                    this->myPropagator->setIndexOfEpochInStateVec(7);
                    this->myPropagator->setCurrentIndependentVariable(this->StateBeforeEventBeforePropagation(7));
                    this->myPropagator->propagate(timeToPropagate, false);
                    output_state(7) = this->StateBeforeEventBeforePropagation(7) + timeToPropagate;

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
                    if ((totalPropagationTime - timeToPropagate) > EphemerisOutputResolution)
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
                    timeToPropagate += (totalPropagationTime - timeToPropagate) > EphemerisOutputResolution
                        ? EphemerisOutputResolution
                        : (totalPropagationTime - timeToPropagate);
                }

                //return the propagator to what it's supposed to be doing
                this->myPropagator->setStateRight(this->state_before_event);
            }//end while loop over wait time
        }//end output_ephemeris()
    }//end namespace BoundaryEvents
}//end namespace EMTG