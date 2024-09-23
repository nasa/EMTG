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

//EMTGv9 PSFBphase class (parallel shooting with finite burn)
//Jacob Englander 2-26-2018

#include "PSFBphase.h"

#include "IntegrationSchemeFactory.h"
#include "PropagatorFactory.h"
#include "PSFBstep_factory.h"

namespace EMTG
{
    namespace Phases
    {
        PSFBphase::PSFBphase(const std::string& name,
                const size_t& journeyIndex,
                const size_t& phaseIndex,
                size_t& stageIndex,
                phase* previousPhase,
                Astrodynamics::universe* Universe,
                HardwareModels::Spacecraft* mySpacecraft,
                HardwareModels::LaunchVehicle* myLaunchVehicle,
                missionoptions* myOptions) :
            ParallelShootingPhase::ParallelShootingPhase(name,
                journeyIndex,
                phaseIndex,
                stageIndex,
                previousPhase,
                Universe,
                mySpacecraft,
                myLaunchVehicle,
                myOptions)
        {
            //acceleration model object
            this->mySpacecraftAccelerationModel = new Astrodynamics::SpacecraftAccelerationModel(this->myOptions,
                this->myJourneyOptions,
                this->myUniverse,
                this->Xdescriptions,
                this->mySpacecraft,
                14); // STM size
            this->mySpacecraftAccelerationModel->setDutyCycle(this->PhaseDutyCycle);

            //EOM
            this->myEOM.setSpacecraftAccelerationModel(this->mySpacecraftAccelerationModel);

            //integration scheme
            this->myIntegrationScheme = CreateIntegrationScheme(&this->myEOM, this->numStatesToPropagate, 14);

            //integrators for initial and terminal coasts
            if (this->hasInitialCoast)
            {
                this->InitialCoastPropagatorObject = CreatePropagator(this->myOptions,
                    this->myUniverse,                    
                    this->numStatesToPropagate,
                    14,
                    this->state_after_initial_TCM,
                    this->state_after_initial_coast,
                    this->STM_initial_coast,
                    this->InitialCoast_dStatedIndependentVariable,
                    (Integration::Integrand*) &this->myEOM,
                    this->myIntegrationScheme,
                    &this->ForcedCoast_dStepSize_dPropagationVariable,
                    this->myJourneyOptions->override_integration_step_size
                        ? this->myJourneyOptions->integration_step_size
                        : this->myOptions->integration_time_step_size);
            }

            if (this->hasTerminalCoast)
            {
                this->TerminalCoastPropagatorObject = CreatePropagator(this->myOptions,
                    this->myUniverse,                    
                    this->numStatesToPropagate,
                    14,
                    this->state_at_end_of_phase,
                    this->state_before_terminal_coast,
                    this->STM_terminal_coast,
                    this->TerminalCoast_dStatedIndependentVariable,
                    (Integration::Integrand*) &this->myEOM,
                    this->myIntegrationScheme,
                    &this->ForcedCoast_dStepSize_dPropagationVariable,
                    this->myJourneyOptions->override_integration_step_size
                        ? this->myJourneyOptions->integration_step_size
                        : this->myOptions->integration_time_step_size);
            }

            //steps
            ParallelShootingStep* previousStep = NULL;
            for (size_t stepIndex = 0; stepIndex < this->num_steps; ++stepIndex)
            {
                this->mySteps.push_back(create_PSFB_step(this->name + "_Step" + std::to_string(stepIndex),
                    this->journeyIndex,
                    this->phaseIndex,
                    stepIndex,
                    this->stageIndex,
                    previousStep,
                    this,
                    this->myUniverse,
                    this->mySpacecraft,
                    this->myOptions));

                previousStep = &this->mySteps.back();
            }
        }//end constructor

        void PSFBphase::output_ephemeris(std::ofstream& outputfile, std::ofstream& acceleration_model_file)
        {
            // Reset the spacecraft
            this->mySpacecraft->setActiveStage(this->stageIndex);

            this->EphemerisOutputResolution = this->myOptions->integration_time_step_size;

            //Step 1: output the departure event
            this->myDepartureEvent->output_ephemeris(outputfile);

            math::Matrix<doubleType> temp_state;
            if (this->myOptions->generate_acceleration_model_instrumentation_file)
            {
                temp_state = this->myDepartureEvent->get_state_before_event();
                this->mySpacecraftAccelerationModel->populateInstrumentationFile(acceleration_model_file, temp_state, temp_state(7));
            }

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
                    temp_state = output_state;
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
                    {
                        if (this->myOptions->generate_acceleration_model_instrumentation_file)
                        {
                            this->mySpacecraftAccelerationModel->populateInstrumentationFile(acceleration_model_file, temp_state, temp_state(7));
                        }

                        this->write_ephemeris_line(outputfile,
                            output_state,
                            math::Matrix<doubleType>(3, 1, 0.0),//control vector
                            0.0,
                            0.0,
                            0.0,
                            0,
                            0.0,
                            "none");
                    }

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
                    temp_state = output_state;
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
                    {
                        if (this->myOptions->generate_acceleration_model_instrumentation_file)
                        {
                            this->mySpacecraftAccelerationModel->populateInstrumentationFile(acceleration_model_file, temp_state, temp_state(7));
                        }

                        this->write_ephemeris_line(outputfile,
                            output_state,
                            math::Matrix<doubleType>(3, 1, 0.0),//control vector
                            0.0,
                            0.0,
                            0.0,
                            0,
                            0.0,
                            "none");
                    }

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
            if (this->myOptions->generate_acceleration_model_instrumentation_file)
            {
                temp_state = this->myArrivalEvent->get_state_before_event();
                this->mySpacecraftAccelerationModel->populateInstrumentationFile(acceleration_model_file, temp_state, temp_state(7));
            }
        }//end output_ephemeris

        PSFBphase::~PSFBphase()
        {            
            delete this->myIntegrationScheme;

            delete this->mySpacecraftAccelerationModel;

            //forced coast STMs and such
            if (this->hasInitialCoast)
            {
                delete this->InitialCoastPropagatorObject;
            }

            if (this->hasTerminalCoast)
            {
                delete this->TerminalCoastPropagatorObject;
            }
        }//end destructor
    }//close namespace Phases
}//close namespace EMTG