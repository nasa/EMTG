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

//forward subphase for EMTGv9 MGAnDSMs
//Jacob Englander 8-24-2017

#include "Forward_MGAnDSMs_subphase.h"
#include "MGAnDSMs_maneuver_constraint_factory.h"
#include "PropagatorFactory.h"
#include "IntegrationSchemeFactory.h"
#include "MGAnDSMs_phase.h"
#include <algorithm>

namespace EMTG
{
    namespace Phases
    {
        Forward_MGAnDSMs_subphase::Forward_MGAnDSMs_subphase(const std::string& name,
            const size_t& journeyIndex,
            const size_t& phaseIndex,
            const size_t& subphaseIndex,
            const size_t& stageIndex,
            MGAnDSMs_phase* myPhase,
            MGAnDSMs_subphase* previousSubPhase,
            Astrodynamics::universe* myUniverse,
            HardwareModels::Spacecraft* mySpacecraft,
            missionoptions* myOptions) :
            MGAnDSMs_subphase(name,
                journeyIndex,
                phaseIndex,
                subphaseIndex,
                stageIndex,
                myPhase,
                previousSubPhase,
                myUniverse,
                mySpacecraft,
                myOptions)
        {
            this->name = name + "ForwardSubPhase" + std::to_string(subphaseIndex);
            if (this->subphaseIndex == 0)
                this->isFirstSubPhase = true;
            else
                this->isFirstSubPhase = false;
            
            size_t numberOfDSMs = this->myJourneyOptions->impulses_per_phase;
            div_t divresult = div(numberOfDSMs, 2);
            size_t subphaseMatchPointIndex;
            if (divresult.rem == 0 && numberOfDSMs > 0)
                subphaseMatchPointIndex = divresult.quot - 1;
            else
                subphaseMatchPointIndex = divresult.quot;

            size_t numberOfForwardSubphases = subphaseMatchPointIndex + 1;

            if (this->subphaseIndex == numberOfForwardSubphases - 1)
                this->BordersMatchPoint = true;
            else
                this->BordersMatchPoint = false;

            this->dMassAfterTCM_dDSMcomponents.resize(3, 1, 0.0);

            this->dFuelConsumedDSM_dMassAtDSM = 0.0;
            this->dOxidizerConsumedDSM_dMassAtDSM = 0.0;
            this->dMassAfterDSM_dMassAtDSM = 1.0;
            this->dFuelConsumedTCM_dMassAtTCM = 0.0;
            this->dOxidizerConsumedTCM_dMassAtTCM = 0.0;
            this->dMassAfterTCM_dMassAtTCM = 1.0;
            this->dFuelConsumedDSM_ddeltav = 0.0;
            this->dOxidizerConsumedDSM_ddeltav = 0.0;
            this->dMassAfterDSM_ddeltav = 0.0;
            this->dFuelConsumedTCM_ddeltav = 0.0;
            this->dOxidizerConsumedTCM_ddeltav = 0.0;
            this->dMassAfterTCM_ddeltav = 0.0;

            //maneuver constraints
            //first we need to define the prefix defines this maneuver and can be compared to the maneuver constraint definitions
            std::string ManeuverTag = "p" + std::to_string(this->phaseIndex)
                + "b" + std::to_string(this->subphaseIndex);

            std::string ManeuverEndTag = "pEndb" + std::to_string(this->subphaseIndex);

            std::string ManeuverEndLastBurnTag = "pEndbEnd";

            for (std::string& constraint : this->myJourneyOptions->ManeuverConstraintDefinitions)
            {
                if (constraint.find("#") != 0) //don't create a constraint if it is commented out
                {
                    if (constraint.find(ManeuverTag) < 1024
                        || (constraint.find(ManeuverEndTag) < 1024 && this->myPhase->getIsLastPhaseInJourney())
                        || (constraint.find(ManeuverEndLastBurnTag) < 1024 && this->myJourneyOptions->impulses_per_phase == 1))
                    {
                        if (constraint.find("epoch") < 1024
                            || constraint.find("magnitude") < 1024)
                        {
                            this->myManeuverConstraints.push_back(create_MGAnDSMs_maneuver_constraint(name,
                                journeyIndex,
                                phaseIndex,
                                subphaseIndex,
                                stageIndex,
                                this,
                                myUniverse,
                                mySpacecraft,
                                myOptions,
                                constraint,
                                "Forward"));
                        }//end epoch and magnitude constraint setup
                        else if (constraint.find("biprop") < 1024)
                        {
                            this->ChemicalManeuverType = PropulsionSystemChoice::Biprop;
                        }
                        else if (constraint.find("monoprop") < 1024)
                        {
                            this->ChemicalManeuverType = PropulsionSystemChoice::Monoprop;
                        }
                    }
                }//end constraint comment check
            }//end maneuver constraint setup
        }//end constructor

        void Forward_MGAnDSMs_subphase::configure_propagator()
        {
            //which propagator do we want?
            if (this->isKeplerian)
            {
                this->myPropagator = Astrodynamics::CreatePropagator(this->myOptions,
                    this->myUniverse,
                    6,
                    *this->spacecraft_state_minus_pointer,
                    this->StateAfterPropagationBeforeDSM,
                    this->STM,
                    *this->dPropagatedStatedIndependentVariable_pointer,
                    &this->BurnIndexDouble);
            }
            else //integrated propagator
            {
                //acceleration model
                this->mySpacecraftAccelerationModel = new Astrodynamics::SpacecraftAccelerationModel(this->myOptions,
                    this->myJourneyOptions,
                    this->myUniverse,
                    this->Xdescriptions,
                    this->mySpacecraft,
                    11); // STM size
                this->mySpacecraftAccelerationModel->setDutyCycle(1.0);

                //EOM
                this->myEOM = Astrodynamics::TimeDomainSpacecraftEOM();
                this->myEOM.setSpacecraftAccelerationModel(this->mySpacecraftAccelerationModel);

                //integration scheme
                this->myIntegrationScheme = CreateIntegrationScheme(&this->myEOM, 10, 11);

                //propagator
                this->myPropagator = Astrodynamics::CreatePropagator(this->myOptions,
                    this->myUniverse,
                    10,
                    11,
                    *this->spacecraft_state_minus_pointer,
                    this->StateAfterPropagationBeforeDSM,
                    this->STM,
                    *this->dPropagatedStatedIndependentVariable_pointer,
                    (Integration::Integrand*) &this->myEOM,
                    this->myIntegrationScheme,
                    &this->BurnIndexDouble,
                    this->myJourneyOptions->override_integration_step_size
                        ? this->myJourneyOptions->integration_step_size
                        : this->myOptions->integration_time_step_size);
            }
        }//end configure_propagator()

        //************************************calcbounds methods
        void Forward_MGAnDSMs_subphase::calcbounds(std::vector<size_t>& timeVariables)
        {
            this->timeVariables = timeVariables;

            this->First_X_entry_this_subphase = this->Xdescriptions->size();

            this->configure_propagator();

            this->calcbounds_subphase_time();

            this->calcbounds_DSM_components();
        }//end calcbounds()

        void Forward_MGAnDSMs_subphase::calcbounds_subphase_time()
        {
            //burn index variable
            if (this->subphaseIndex > 0)
                this->Xindex_BurnIndex = this->previousSubPhase->getXindex_burnIndex();

            this->Xlowerbounds->push_back(1.0e-4);
            this->Xupperbounds->push_back(1.0 - 1.0e-4);
            this->X_scale_factors->push_back(1.0);
            this->Xdescriptions->push_back(prefix + "burn index");
            this->Xindex_BurnIndex.push_back(this->Xdescriptions->size() - 1);
        }//end calcbounds_subphase_time()

         //************************************process methods
        void Forward_MGAnDSMs_subphase::process_subphase(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //Step 1: assign references to our various important things
            math::Matrix<doubleType>& state_before_subphase = *this->spacecraft_state_minus_pointer;
            math::Matrix<doubleType>& state_after_subphase = *this->spacecraft_state_plus_pointer;
            doubleType& PhaseFlightTime = *this->PhaseFlightTime_pointer;
            math::Matrix<double>& SPTM = *this->SPTM_pointer;
            math::Matrix<double>& dPropagatedStatedIndependentVariable = *this->dPropagatedStatedIndependentVariable_pointer;

            //Step 2: do fun things with decision variables
            this->process_subphase_time(X, Xindex, F, Findex, G, needG);

            this->process_DSM_components(X, Xindex, F, Findex, G, needG);

            //now do somewhat less fun things with propagation

            //Step 3: propagate
            this->myPropagator->setCurrentEpoch(state_before_subphase(7));
            this->myPropagator->setIndexOfEpochInStateVec(7);
            this->myPropagator->setCurrentIndependentVariable(state_before_subphase(7));
            this->myPropagator->propagate(this->SubPhaseTime, needG);

            this->ACS_fuel_used = (this->myOptions->trackACS ? this->myOptions->ACS_kg_per_day * this->SubPhaseTime / 86400.0 : 0.0);
            this->chemical_fuel_used += this->ACS_fuel_used;

            if (this->isKeplerian)
            {
                this->StateAfterPropagationBeforeDSM(6) = state_before_subphase(6) - this->ACS_fuel_used;
                this->StateAfterPropagationBeforeDSM(7) = state_before_subphase(7) + this->SubPhaseTime;
                this->StateAfterPropagationBeforeDSM(8) = state_before_subphase(8) + this->ACS_fuel_used; //virtual fuel
                this->StateAfterPropagationBeforeDSM(9) = state_before_subphase(9); //virtual oxidizer
            }

            //Step 4: perform the DSM
            //position and velocity
            for (size_t Vindex = 0; Vindex < 3; ++Vindex)
            {
                this->StateAfterDSMBeforeTCM(Vindex) = this->StateAfterPropagationBeforeDSM(Vindex);
                this->StateAfterDSMBeforeTCM(3 + Vindex) = this->StateAfterPropagationBeforeDSM(3 + Vindex) + this->DSM(Vindex);
            }

            //mass
            this->chemical_fuel_used = math::SMALL; //leak!
            this->chemical_oxidizer_used = math::SMALL; //leak!

            this->mySpacecraft->computeChemicalPropulsionPerformance(this->DSM_magnitude, this->StateAfterPropagationBeforeDSM(6), true, this->ChemicalManeuverType);
            this->DSM_fuel_used = this->mySpacecraft->getChemFuelConsumedThisManeuver();
            this->chemical_fuel_used += this->DSM_fuel_used;
            this->chemical_oxidizer_used += this->mySpacecraft->getChemOxidizerConsumedThisManeuver();
            this->StateAfterDSMBeforeTCM(6) = this->StateAfterPropagationBeforeDSM(6) - this->DSM_fuel_used - this->chemical_oxidizer_used;

            this->dFuelConsumedDSM_dMassAtDSM = this->mySpacecraft->get_dFuelConsumedThisManeuver_dMassAtManeuver();
            this->dOxidizerConsumedDSM_dMassAtDSM = this->mySpacecraft->get_dOxidizerConsumedThisManeuver_dMassAtManeuver();
            this->dMassAfterDSM_dMassAtDSM = this->mySpacecraft->get_dMassAfterManeuver_dMassAtManeuver();
            this->dFuelConsumedDSM_ddeltav = this->mySpacecraft->getdFuelConsumedThisManeuver_ddeltav();
            this->dOxidizerConsumedDSM_ddeltav = this->mySpacecraft->getdOxidizerConsumedThisManeuver_ddeltav();
            this->dMassAfterDSM_ddeltav = this->mySpacecraft->getdSpacecraftMassThisManeuver_ddeltav();

            //time
            this->StateAfterDSMBeforeTCM(7) = this->StateAfterPropagationBeforeDSM(7);
            this->DSMepoch = this->StateAfterDSMBeforeTCM(7);

            //tanks
            this->StateAfterDSMBeforeTCM(8) = StateAfterPropagationBeforeDSM(8) + this->DSM_fuel_used; //virtual fuel
            this->StateAfterDSMBeforeTCM(9) = StateAfterPropagationBeforeDSM(9) + this->chemical_oxidizer_used; //virtual oxidizer

            //Step 5: perform the TCM if necessary
            if (this->hasTCM)
            {
                this->TCM_magnitude = this->DSM_magnitude * this->myOptions->TCM_maneuver_fraction;

                //position and velocity
                for (size_t stateIndex = 0; stateIndex < 6; ++stateIndex)
                    state_after_subphase(stateIndex) = this->StateAfterDSMBeforeTCM(stateIndex);

                //mass
                this->mySpacecraft->computeChemicalPropulsionPerformance(this->TCM_magnitude, this->StateAfterDSMBeforeTCM(6), true, PropulsionSystemChoice::Monoprop);
                this->TCM_fuel_used = this->mySpacecraft->getChemFuelConsumedThisManeuver();
                this->chemical_fuel_used += this->TCM_fuel_used;
                state_after_subphase(6) = this->StateAfterDSMBeforeTCM(6) - this->TCM_fuel_used;

                this->dFuelConsumedTCM_dMassAtTCM = this->mySpacecraft->get_dFuelConsumedThisManeuver_dMassAtManeuver();
                this->dOxidizerConsumedTCM_dMassAtTCM = this->mySpacecraft->get_dOxidizerConsumedThisManeuver_dMassAtManeuver();
                this->dMassAfterTCM_dMassAtTCM = this->mySpacecraft->get_dMassAfterManeuver_dMassAtManeuver();
                this->dFuelConsumedTCM_ddeltav = this->mySpacecraft->getdFuelConsumedThisManeuver_ddeltav() * this->myOptions->TCM_maneuver_fraction;
                this->dOxidizerConsumedTCM_ddeltav = this->mySpacecraft->getdOxidizerConsumedThisManeuver_ddeltav() * this->myOptions->TCM_maneuver_fraction;
                this->dMassAfterTCM_ddeltav = this->mySpacecraft->getdSpacecraftMassThisManeuver_ddeltav() * this->myOptions->TCM_maneuver_fraction;


                //time
                state_after_subphase(7) = this->StateAfterDSMBeforeTCM(7);

                //tanks
                state_after_subphase(8) = StateAfterDSMBeforeTCM(8) + this->TCM_fuel_used; //virtual fuel
                state_after_subphase(9) = StateAfterDSMBeforeTCM(9); //virtual oxidizer
            }
            else
            {

                for (size_t stateIndex = 0; stateIndex < 10; ++stateIndex)
                    state_after_subphase(stateIndex) = this->StateAfterDSMBeforeTCM(stateIndex);
            }

            if (needG)
            {
                //Step 6: form the SPTM

                //Step 6.0: clear the SPTM
                SPTM.construct_identity_matrix();

                //Step 6.1: upper left stateMax x stateMax is the regular STM
                size_t stateMax = this->isKeplerian ? 6 : 10;
                for (size_t i = 0; i < stateMax; ++i)
                    for (size_t j = 0; j < stateMax; ++j)
                        SPTM(i, j) = this->STM(i, j);
                //Step 6.2: turn the upper right stateMax x 2 into the Phi_t terms developed by Lantoine for time and for burn index
                for (size_t i = 0; i < stateMax; ++i)
                {
                    //this really ought to be cleaned up...
                    if (this->isKeplerian)
                    {
                        SPTM(i, 10) = dPropagatedStatedIndependentVariable(i);
                        SPTM(i, 11) = dPropagatedStatedIndependentVariable(i)
                            * (PhaseFlightTime / this->BurnIndex) _GETVALUE;
                    }
                    else //integrated propagator
                    {
                        SPTM(i, 10) = this->STM(i, 10);
                        SPTM(i, 11) = this->STM(i, 10)
                            * (PhaseFlightTime / this->BurnIndex) _GETVALUE;
                    }
                }
                SPTM(7, 10) = this->BurnIndex _GETVALUE;
                SPTM(7, 11) = PhaseFlightTime _GETVALUE;

                //Step 6.3: mass
                SPTM(6, 6) = (state_after_subphase(6) / StateAfterPropagationBeforeDSM(6)) _GETVALUE;
                SPTM(6, 10) *= SPTM(6, 6);
                SPTM(6, 11) *= SPTM(6, 6);
                if (this->myOptions->trackACS)
                {
                    SPTM(6, 10) += -this->BurnIndex _GETVALUE * this->myOptions->ACS_kg_per_day / 86400.0 * SPTM(6, 6) / 2.0;
                    SPTM(6, 11) += -PhaseFlightTime _GETVALUE * this->myOptions->ACS_kg_per_day / 86400.0 * SPTM(6, 6);
                }

                //Step 6.4: tanks
                //derivatives of virtual fuel
                SPTM(8, 6) = ((this->chemical_fuel_used - this->ACS_fuel_used - math::SMALL) / StateAfterPropagationBeforeDSM(6)) _GETVALUE;
                SPTM(8, 10) = SPTM(6, 10) * this->dFuelConsumedDSM_dMassAtDSM;
                SPTM(8, 11) = SPTM(6, 10) * this->dFuelConsumedDSM_dMassAtDSM * (PhaseFlightTime / this->BurnIndex) _GETVALUE;
                if (this->myOptions->trackACS)
                {
                    SPTM(8, 10) += (this->BurnIndex * this->myOptions->ACS_kg_per_day / 86400.0 *
                        (1.0
                            - this->dFuelConsumedDSM_dMassAtDSM
                            - this->dMassAfterDSM_dMassAtDSM * this->dFuelConsumedTCM_dMassAtTCM)
                        )_GETVALUE;
                    SPTM(8, 11) += (PhaseFlightTime * this->myOptions->ACS_kg_per_day / 86400.0 *
                        (1.0
                            - this->dFuelConsumedDSM_dMassAtDSM
                            - this->dMassAfterDSM_dMassAtDSM * this->dFuelConsumedTCM_dMassAtTCM)
                        )_GETVALUE;
                }
                //derivatives of virtual oxidizer
                SPTM(9, 6) = ((this->chemical_oxidizer_used - math::SMALL) / StateAfterPropagationBeforeDSM(6)) _GETVALUE;
                if (this->myPropagatorType == PropagatorType::IntegratedPropagator)
                {
                    SPTM(9, 10) = STM(9, 10) + SPTM(6, 10) * this->dOxidizerConsumedDSM_dMassAtDSM;
                    SPTM(9, 11) = (STM(9, 10) + SPTM(6, 10) * this->dOxidizerConsumedDSM_dMassAtDSM) * (PhaseFlightTime / this->BurnIndex) _GETVALUE;
                }
                else
                {
                    SPTM(9, 10) = SPTM(6, 10) * this->dOxidizerConsumedDSM_dMassAtDSM;
                    SPTM(9, 11) = SPTM(6, 10) * this->dOxidizerConsumedDSM_dMassAtDSM * (PhaseFlightTime / this->BurnIndex) _GETVALUE;

                    if (this->myOptions->trackACS)
                    {
                        SPTM(9, 10) += (this->BurnIndex * this->myOptions->ACS_kg_per_day / 86400.0 * (-this->dOxidizerConsumedDSM_dMassAtDSM))_GETVALUE;
                        SPTM(9, 11) += (PhaseFlightTime * this->myOptions->ACS_kg_per_day / 86400.0 * (-this->dOxidizerConsumedDSM_dMassAtDSM))_GETVALUE;
                    }
                }

                //Step 7: derivatives of mass with respect to DSM components
                for (size_t Vindex = 0; Vindex < 3; ++Vindex)
                {
                    this->dMassAfterTCM_dDSMcomponents(Vindex) = (this->dMassAfterDSM_ddeltav * this->dMassAfterTCM_dMassAtTCM
                        + this->dMassAfterTCM_ddeltav)
                        * (this->DSM(Vindex) / this->DSM_magnitude) _GETVALUE;

                    this->dFuel_dDSMcomponents(Vindex) = -(this->dFuelConsumedDSM_ddeltav
                        + this->dMassAfterDSM_ddeltav * (this->dMassAfterTCM_dMassAtTCM - 1.0) + this->dMassAfterTCM_ddeltav)
                        * (this->DSM(Vindex) / this->DSM_magnitude) _GETVALUE;

                    this->dOxidizer_dDSMcomponents(Vindex) = -this->dOxidizerConsumedDSM_ddeltav
                        * (this->DSM(Vindex) / this->DSM_magnitude) _GETVALUE;
                }
            }//end derivatives
        }//end process_subphase()

        void Forward_MGAnDSMs_subphase::process_subphase_time(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //Step 1: extract the burn index
            this->BurnIndex = X[Xindex++];
            this->BurnIndexDouble = this->BurnIndex _GETVALUE;

            //Step 2: compute the subphase time
            doubleType& PhaseFlightTime = *this->PhaseFlightTime_pointer;
            this->SubPhaseTime = this->BurnIndex * PhaseFlightTime;
        }//end process_subphase_time()

         //************************************output methods
        void Forward_MGAnDSMs_subphase::output(std::ofstream& outputfile,
            size_t& eventcount)
        {
            //first we will propagate forward and print coast states at regular intervals
            size_t steps_per_subphase = std::max(1, (int)(this->BurnIndex _GETVALUE * this->num_timesteps));
            doubleType output_timestep = this->SubPhaseTime / steps_per_subphase;
            math::Matrix<doubleType>& state_before_subphase = *this->spacecraft_state_minus_pointer;

            //set the propagator temporarily to dump into the output state
            this->myPropagator->setStateRight(this->output_state);

            for (size_t step = 0; step < steps_per_subphase; ++step)
            {
                //propagate to the halfway point of this step
                this->myPropagator->setCurrentEpoch(state_before_subphase(7));
                this->myPropagator->setIndexOfEpochInStateVec(7);
                this->myPropagator->setCurrentIndependentVariable(state_before_subphase(7));
                this->myPropagator->propagate(output_timestep * (step + 0.5), false);

                if (this->isKeplerian)
                {
                    doubleType ACS_fuel_used = (this->myOptions->trackACS ? this->myOptions->ACS_kg_per_day * output_timestep * (step + 0.5) / 86400.0 : 0.0);
                    this->output_state(6) = this->spacecraft_state_minus_pointer->getentry(6, 0) - ACS_fuel_used;
                    this->output_state(7) = this->spacecraft_state_minus_pointer->getentry(7, 0) + output_timestep * (step + 0.5);
                }

                //print
                this->write_output_line(outputfile,
                    eventcount,
                    "coast",
                    "deep-space",
                    output_timestep / 86400.0,
                    0.0,
                    0.0,
                    this->output_state,
                    this->empty3,
                    0.0,
                    0.0);
            }

            //reset the propagator to go where it is supposed to
            this->myPropagator->setStateRight(this->StateAfterPropagationBeforeDSM);

            //then we will print the DSM
            double myIsp = (this->ChemicalManeuverType == PropulsionSystemChoice::Biprop ? this->mySpacecraft->getBipropIsp() : this->mySpacecraft->getMonopropIsp());
            this->write_output_line(outputfile,
                eventcount,
                "chem_burn",
                "deep-space",
                0.0,
                atan2(this->DSM(1), this->DSM(0)),
                asin(this->DSM(2) / this->DSM_magnitude),
                this->StateAfterDSMBeforeTCM,
                this->DSM,
                this->DSM_magnitude,
                myIsp);

            //then, if applicable, we will print the TCM
            if (this->hasTCM)
            {
                this->write_output_line(outputfile,
                    eventcount,
                    "TCM",
                    "deep-space",
                    0.0,
                    0.0,
                    0.0,
                    *(this->spacecraft_state_plus_pointer),
                    this->empty3,
                    this->TCM_magnitude,
                    this->mySpacecraft->getMonopropIsp());
            }
        }//end output()

        void Forward_MGAnDSMs_subphase::output_ephemeris(std::ofstream& outputfile, std::ofstream & acceleration_model_file)
        {
            //Step 0: we'll need an output vector
            math::Matrix<doubleType> state_before_subphase = *this->spacecraft_state_minus_pointer;
            math::Matrix<doubleType> current_left_hand_state = state_before_subphase;

            //Step 1: temporarily assign the propagator to the output state
            this->myPropagator->setStateLeft(current_left_hand_state);
            this->myPropagator->setStateRight(this->output_state);

            //Step 3.2: propagate and print, skipping the first entry
            double cumulative_time_propagated = this->EphemerisOutputResolution;
            while (cumulative_time_propagated < this->SubPhaseTime)
            {
                //Step 3.2.1: propagate
                this->myPropagator->setCurrentEpoch(current_left_hand_state(7));
                this->myPropagator->setIndexOfEpochInStateVec(7);
                this->myPropagator->setCurrentIndependentVariable(current_left_hand_state(7));
                this->myPropagator->propagate(this->EphemerisOutputResolution, false);

                if (this->isKeplerian)
                {
                    this->output_state(6) = state_before_subphase(6);
                    this->output_state(7) = state_before_subphase(7) + cumulative_time_propagated;
                }
                // keep a temp state that is not converted to heliocentric
                this->temp_state = this->output_state;
                current_left_hand_state = this->output_state;

                //Step 3.2.2: convert to Sun-centered if necessary
                if (this->myUniverse->central_body.spice_ID != this->myOptions->forward_integrated_ephemeris_central_body_SPICE_ID)
                {
                    double LT_dump;
                    double bodyStateDouble[6];
                    spkez_c(this->myUniverse->central_body.spice_ID, output_state(7)_GETVALUE - (51544.5 * 86400.0), "J2000", "NONE", this->myOptions->forward_integrated_ephemeris_central_body_SPICE_ID, bodyStateDouble, &LT_dump);

                    for (size_t stateIndex = 0; stateIndex < 6; ++stateIndex)
                        output_state(stateIndex) += bodyStateDouble[stateIndex];
                }

                //Step 3.2.3: print
                if (this->myOptions->generate_acceleration_model_instrumentation_file && this->myPropagatorType == PropagatorType::IntegratedPropagator)
                {
                    this->mySpacecraftAccelerationModel->populateInstrumentationFile(acceleration_model_file, this->temp_state, this->temp_state(7));
                }

                this->write_ephemeris_line(outputfile,
                    this->output_state,
                    math::Matrix<doubleType>(3, 1, 0.0),//control vector
                    0.0,
                    0.0,
                    0.0,
                    0,
                    0.0,
                    "none");

                //Step 3.2.4: increment propagatedEpoch
                cumulative_time_propagated += (this->SubPhaseTime _GETVALUE - cumulative_time_propagated) > this->EphemerisOutputResolution
                    ? this->EphemerisOutputResolution
                    : (this->SubPhaseTime _GETVALUE - cumulative_time_propagated);
            }

            //Step 3.3: reset the propagator to its original state vectors
            this->configure_propagator();
        }//end output_ephemeris()

        void Forward_MGAnDSMs_subphase::output_maneuver_and_target_spec(std::ofstream& maneuver_spec_file, std::ofstream& target_spec_file, bool& haveManeuverNeedTarget)
        {
            //Step 1: write a target spec line if the phase's departure event has a maneuver or if this is a later subphase
            if (haveManeuverNeedTarget)
            {
                //reset the flag that tells us we need a target spec line
                haveManeuverNeedTarget = false;

                //Step 1.1: configure target spec
                target_spec_line myTargetSpecLine(this->name,
                    "EME2000",
                    this->myUniverse->central_body.name,
                    this->StateAfterPropagationBeforeDSM(7),
                    this->StateAfterPropagationBeforeDSM);

                //Step 1.2: write target spec object
                myTargetSpecLine.write(target_spec_file);
            }

            //Step 2: maneuver spec for the DSM
            //Step 2.1: initialize a maneuver spec object
            double Isp = this->ChemicalManeuverType == PropulsionSystemChoice::Biprop ? this->mySpacecraft->getBipropIsp() : this->mySpacecraft->getMonopropIsp();
            doubleType massFlowRate = this->mySpacecraft->getchemthrust() / Isp / this->myOptions->g0;
            doubleType maneuverDuration = (this->StateAfterPropagationBeforeDSM(6) - this->StateAfterDSMBeforeTCM(6)) / massFlowRate;
            maneuver_spec_line myManeuverSpecLine(this->name + "_DSM");
            myManeuverSpecLine.append_maneuver_spec_item("EME2000",
                this->StateAfterPropagationBeforeDSM(7),
                this->DSM.unitize(),
                this->StateAfterPropagationBeforeDSM(6),
                this->StateAfterDSMBeforeTCM(6),
                this->mySpacecraft->getchemthrust(),
                massFlowRate,
                maneuverDuration,
                1.0);

            //Step 2.2: write maneuver spec object
            myManeuverSpecLine.write(maneuver_spec_file);

            //Step 2.3: signal that we need a target spec
            haveManeuverNeedTarget = true;
        }

    }//close namespace Phases
}//close namespace EMTG