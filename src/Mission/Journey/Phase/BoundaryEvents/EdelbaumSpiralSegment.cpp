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

//Edelbaum spiral segment
//Jacob Englander 10-24-2017

#include "EdelbaumSpiralSegment.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        EdelbaumSpiralSegment::EdelbaumSpiralSegment() :
            myDutyCycle(1.0),
            deltav(0.0),
            SpiralSegmentTime(0.0),
            mass_before_segment(1.0),
            mass_after_segment(1.0),
            virtual_electric_propellant(0.0),
            virtual_chemical_fuel(0.0),
            StateRelativeToJourneyCentralBody(math::Matrix<doubleType>(8, 1, 0.0)),
            PositionRelativeToSun(math::Matrix<doubleType>(3, 1, 0.0)),
            dPositionRelativeToSun_dt(math::Matrix<double>(3, 1, 0.0)),
            PositionCB_Sun(math::Matrix<doubleType>(3, 1, 0.0)),
            dPositionCB_Sun_dt(math::Matrix<double>(3, 1, 0.0))
        {}

        EdelbaumSpiralSegment::EdelbaumSpiralSegment(const std::string& name,
            const size_t& journeyIndex,
            const size_t& phaseIndex,
            const size_t& stageIndex,
            const size_t& spiralSegmentIndex,
            Astrodynamics::universe* myUniverse,
            Astrodynamics::body* myBody,
            HardwareModels::Spacecraft* mySpacecraft,
            missionoptions* myOptions,
            EdelbaumSpiral* mySpiral) :
            EdelbaumSpiralSegment::EdelbaumSpiralSegment()
        {
            this->name = name;
            this->journeyIndex = journeyIndex;
            this->phaseIndex = phaseIndex;
            this->stageIndex = stageIndex;
            this->spiralSegmentIndex = spiralSegmentIndex;
            this->myUniverse = myUniverse;
            this->myBody = myBody;
            this->mySpacecraft = mySpacecraft;
            this->myOptions = myOptions;
            this->myJourneyOptions = &this->myOptions->Journeys[this->journeyIndex];
            this->mySpiral = mySpiral;

            this->isHeliocentric = this->myUniverse->central_body_SPICE_ID == 10 ? true : false;

            //initialize the writey thing
            this->writey_thing::initialize(this->myOptions, this->myUniverse);
        }//end constructor
        
        void EdelbaumSpiralSegment::setup_calcbounds(std::vector<double>* Xupperbounds,
            std::vector<double>* Xlowerbounds,
            std::vector<double>* X_scale_factors,
            std::vector<double>* Fupperbounds,
            std::vector<double>* Flowerbounds,
            std::vector<double>* F_scale_factors,
            std::vector<std::string>* Xdescriptions,
            std::vector<std::string>* Fdescriptions,
            std::vector<size_t>* iGfun,
            std::vector<size_t>* jGvar,
            std::vector<std::string>* Gdescriptions,
            std::vector<size_t>* iAfun,
            std::vector<size_t>* jAvar,
            std::vector<std::string>* Adescriptions,
            std::vector<double>* A)
        {
            this->prefix = this->name + ": ";

            this->sparsey_thing::setup_calcbounds(Xupperbounds,
                Xlowerbounds,
                X_scale_factors,
                Fupperbounds,
                Flowerbounds,
                F_scale_factors,
                Xdescriptions,
                Fdescriptions,
                iGfun,
                jGvar,
                Gdescriptions,
                iAfun,
                jAvar,
                Adescriptions,
                A);
        }//end setup_calcbounds()

        void EdelbaumSpiralSegment::calcbounds()
        {
            //left-side epoch
            this->calculate_dependencies_epoch_time();

            //Step 3: spiral segment time
            //Step 3.1: let's start by hallucinating an upper bound on the flight time.
            double segment_flight_time_upper_bound = 1000.0 * 86400.0; //why not? easier than doing something smart

            //Step 3.2: now let's make the variable
            this->Xlowerbounds->push_back(math::SMALL);
            this->Xupperbounds->push_back(segment_flight_time_upper_bound);
            this->X_scale_factors->push_back(this->myUniverse->TU);
            this->Xdescriptions->push_back(this->prefix + "segment flight time");
            this->Xindex_SpiralSegmentTime = this->Xdescriptions->size() - 1;

            //Step 4: mass after segment
            this->Xlowerbounds->push_back(math::SMALL);
            this->Xupperbounds->push_back(this->myOptions->maximum_mass);
            this->X_scale_factors->push_back(this->myOptions->maximum_mass);
            this->Xdescriptions->push_back(this->prefix + "mass after segment");
            this->Xindex_mass_after_segment = this->Xdescriptions->size() - 1;

            //Step 5: virtual electric propellant
            {
                this->Xlowerbounds->push_back(math::SMALL);
                this->Xupperbounds->push_back(this->myOptions->maximum_mass);
                this->X_scale_factors->push_back(this->myOptions->maximum_mass);
                this->Xdescriptions->push_back(this->prefix + "virtual electric propellant");
                this->Xindex_virtual_electric_propellant = this->Xdescriptions->size() - 1;

                //tell the spacecraft about this variable
                //global constraint
                this->mySpacecraft->appendGlobalElectricPropellantTank_Xindices(this->Xindex_virtual_electric_propellant);
                this->mySpacecraft->appendGlobalElectricPropellantTank_Xscaleranges(this->Xupperbounds->back(), this->Xlowerbounds->back());
                //stage constraint
                this->mySpacecraft->appendElectricPropellantTank_Xindices(this->stageIndex, this->Xindex_virtual_electric_propellant);
                this->mySpacecraft->appendElectricPropellantTank_Xscaleranges(this->stageIndex, this->Xupperbounds->back(), this->Xlowerbounds->back());
            }//end virtual electric propellant tank

            //Step 6: virtual chemical fuel
            {
                this->Xlowerbounds->push_back(math::SMALL);
                this->Xupperbounds->push_back(this->myOptions->maximum_mass);
                this->X_scale_factors->push_back(this->myOptions->maximum_mass);
                this->Xdescriptions->push_back(this->prefix + "virtual chemical fuel");
                this->Xindex_virtual_chemical_fuel = this->Xdescriptions->size() - 1;
                
                //tell the spacecraft about this variable
                //global constraint
                this->mySpacecraft->appendGlobalChemicalFuelTank_Xindices(this->Xindex_virtual_chemical_fuel);
                this->mySpacecraft->appendGlobalChemicalFuelTank_Xscaleranges(this->Xupperbounds->back(), this->Xlowerbounds->back());
                //stage constraint
                this->mySpacecraft->appendChemicalFuelTank_Xindices(this->stageIndex, this->Xindex_virtual_chemical_fuel);
                this->mySpacecraft->appendChemicalFuelTank_Xscaleranges(this->stageIndex, this->Xupperbounds->back(), this->Xlowerbounds->back());
            }//end virtual chemical propellant tank

            //Step 7: continuity constraints
            //Step 7.1: time
            this->Flowerbounds->push_back(-math::SMALL);
            this->Fupperbounds->push_back(math::SMALL);
            this->Fdescriptions->push_back(this->prefix + "time continuity");
            this->Findex_SpiralSegmentTime = this->Fdescriptions->size() - 1;

            //the time continuity constraint has derivatives with respect to all previous time variables, left-hand mass, and segment encoded time
            for (size_t Xindex_previous_time_variable : this->Xindices_EventLeftEpoch)
                this->create_sparsity_entry(this->Findex_SpiralSegmentTime,
                    Xindex_previous_time_variable,
                    this->Gindex_SpiralSegmentTime_wrt_previous_time_variables);

            //the Xindex of the left mass can be found by scanning backward through the decision vector
            for (size_t Xindex = this->Xindex_SpiralSegmentTime; Xindex >= 0; --Xindex)
            {
                if (this->Xdescriptions->operator[](Xindex).find("mass") < 1024)
                {
                    this->Xindex_left_mass = Xindex;
                    break;
                }
            }

            this->create_sparsity_entry(this->Findex_SpiralSegmentTime,
                this->Xindex_left_mass,
                this->Gindex_SpiralSegmentTime_wrt_left_hand_mass);


            this->create_sparsity_entry(this->Findex_SpiralSegmentTime,
                this->Xindex_SpiralSegmentTime,
                this->Gindex_SpiralSegmentTime_wrt_SpiralSegmentTime);

            //Step 7.2: mass
            this->Flowerbounds->push_back(-math::SMALL);
            this->Fupperbounds->push_back(math::SMALL);
            this->Fdescriptions->push_back(this->prefix + "mass continuity");
            this->Findex_mass_after_segment = this->Fdescriptions->size() - 1;

            //the mass continuity constraint has derivatives with respect to all previous time variables, left-hand mass, segment encoded mass, and IF ACS IS ON also segment encoded time
            for (size_t Xindex_previous_time_variable : this->Xindices_EventLeftEpoch)
                this->create_sparsity_entry(this->Findex_mass_after_segment,
                    Xindex_previous_time_variable,
                    this->Gindex_mass_after_segment_wrt_previous_time_variables);

            this->create_sparsity_entry(this->Findex_mass_after_segment,
                this->Xindex_left_mass,
                this->Gindex_mass_after_segment_wrt_left_hand_mass);

            this->create_sparsity_entry(this->Findex_mass_after_segment,
                this->Xindex_mass_after_segment,
                this->Gindex_mass_after_segment_wrt_mass_after_segment);

            if (this->myOptions->trackACS)
            {
                this->create_sparsity_entry(this->Findex_mass_after_segment,
                    this->Xindex_SpiralSegmentTime,
                    this->Gindex_mass_after_segment_wrt_SpiralSegmentTime);
            }

            //Step 7.3: electric propellant
            this->Flowerbounds->push_back(-math::SMALL);
            this->Fupperbounds->push_back(math::SMALL);
            this->Fdescriptions->push_back(this->prefix + "virtual electric propellant");
            this->Findex_virtual_electric_propellant = this->Fdescriptions->size() - 1;
            
            //the electric propellant continuity constraint has derivatives with respect to all previous time variables, left-hand mass, and segment encoded electric propellant
            for (size_t Xindex_previous_time_variable : this->Xindices_EventLeftEpoch)
                this->create_sparsity_entry(this->Findex_virtual_electric_propellant,
                    Xindex_previous_time_variable,
                    this->Gindex_virtual_electric_propellant_wrt_previous_time_variables);

            this->create_sparsity_entry(this->Findex_virtual_electric_propellant,
                this->Xindex_left_mass,
                this->Gindex_virtual_electric_propellant_wrt_left_hand_mass);
            
            this->create_sparsity_entry(this->Findex_virtual_electric_propellant,
                this->Xindex_virtual_electric_propellant,
                this->Gindex_virtual_electric_propellant_wrt_virtual_electric_propellant);

            //Step 7.4: chemical fuel
            this->Flowerbounds->push_back(-math::SMALL);
            this->Fupperbounds->push_back(math::SMALL);
            this->Fdescriptions->push_back(this->prefix + "virtual chemical fuel");
            this->Findex_virtual_chemical_fuel = this->Fdescriptions->size() - 1;

            //the chemical fuel continuity constraint has derivatives with respect to segment encoded time and segment encoded chemical fuel
            if (this->myOptions->trackACS)
            {
                this->create_sparsity_entry(this->Findex_virtual_chemical_fuel,
                    this->Xindex_SpiralSegmentTime,
                    this->Gindex_virtual_chemical_fuel_wrt_SpiralSegmentTime);
            }
            
            this->create_sparsity_entry(this->Findex_virtual_chemical_fuel,
                this->Xindex_virtual_chemical_fuel,
                this->Gindex_virtual_chemical_fuel_wrt_virtual_chemical_fuel);
        }//end calcbounds()

        void EdelbaumSpiralSegment::calculate_dependencies_epoch_time()
        {
            //the left epoch of this event depends on whatever time variables the parent spiral currently has stored
            this->Xindices_EventLeftEpoch = this->mySpiral->get_Xindices_SpiralEndEpoch();
        }//end calculate_dependencies_left_epoch

        void EdelbaumSpiralSegment::process(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //Step 1: when am I?
            this->process_epoch_time(X, Xindex, F, Findex, G, needG);

            //Step 2: extract the thingums
            this->SpiralSegmentTime = X[Xindex++];
            this->mass_after_segment = X[Xindex++];
            this->virtual_electric_propellant = X[Xindex++];
            this->virtual_chemical_fuel = X[Xindex++];

            //Step 3: get the previous mass
            if (this->spiralSegmentIndex == 0)
                this->mass_before_segment = this->mySpiral->getInitialMass();
            else
                this->mass_before_segment = X[this->Xindex_left_mass];



            double& mu = this->myBody->mu;

            //Step 4: where am I?
            doubleType temp_body_state[12];
            doubleType DistanceFromSunAU;

            this->myBody->locate_body(this->segmentLeftEpoch,
                temp_body_state,
                needG,
                *this->myOptions);

            for (size_t stateIndex = 0; stateIndex < 6; ++stateIndex)
                this->StateRelativeToJourneyCentralBody(stateIndex) = temp_body_state[stateIndex];
            this->StateRelativeToJourneyCentralBody(6) = this->mass_before_segment;
            this->StateRelativeToJourneyCentralBody(7) = this->segmentLeftEpoch;

            if (this->isHeliocentric)
            {
                // get the position vector of the sun
                doubleType temp_CB_state[12];
                myUniverse->locate_central_body(this->segmentLeftEpoch, temp_CB_state, *this->myOptions, needG);
                // locate_central_body gives us the state of the central body relative to the sun 
                // we want the state of the sun relative to the central body though (so its negation)
                for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
                {
                    this->PositionCB_Sun(stateIndex) = temp_CB_state[stateIndex];
                    this->dPositionCB_Sun_dt(stateIndex) = temp_CB_state[stateIndex + 6] _GETVALUE;
                    this->PositionRelativeToSun(stateIndex) = this->PositionCB_Sun(stateIndex) + this->StateRelativeToJourneyCentralBody(stateIndex);
                    this->dPositionRelativeToSun_dt(stateIndex) = this->dPositionCB_Sun_dt(stateIndex) + temp_body_state[stateIndex + 6] _GETVALUE;
                }
            }
            else
            {
                // if we are dealing with a heliocentric trajectory, 
                // then the spacecraft's position relative to the sun is the same as its position relative to the central body
                for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
                {
                    this->PositionRelativeToSun(stateIndex) = this->StateRelativeToJourneyCentralBody(stateIndex);
                    this->dPositionRelativeToSun_dt(stateIndex) = temp_body_state[stateIndex + 6] _GETVALUE;
                }
            } // end if isHeliocentric
            DistanceFromSunAU = this->PositionRelativeToSun.norm() / this->myOptions->AU;

            //Step 3: what do I know about my power and propulsion state?
            this->mySpacecraft->computePowerState(DistanceFromSunAU, this->segmentLeftEpoch);

            this->mySpacecraft->computeElectricPropulsionPerformance(this->myDutyCycle);

            //Step 4: compute Isp and mass flow rate
            doubleType Isp = this->mySpacecraft->getEPIsp();
            double& g0 = this->myOptions->g0;
            doubleType delta_mass_expfun = exp(-this->deltav * 1000.0 / (Isp * g0));
            doubleType mdot = this->mySpacecraft->getEPMassFlowRate();

            //Step 5: compute the propellant used
            doubleType mass_after_segment_computed = this->mass_before_segment * delta_mass_expfun;
            this->electric_propellant_used = this->mass_before_segment - mass_after_segment_computed;

            //Step 6: compute the actual flight time
            this->actual_SpiralSegmentTime = this->electric_propellant_used / mdot;

            //Step 7: compute ACS propellant mass drop if applicable
            if (this->myOptions->trackACS)
            {
                this->chemical_fuel_used = this->myOptions->ACS_kg_per_day * this->SpiralSegmentTime / 86400.0;
                mass_after_segment_computed -= this->chemical_fuel_used;
            }

            //Step 8: compute the radius at the end of the spiral (only for output, so not overloaded)
            double v0 = sqrt(this->myBody->mu / *this->radius_before_segment);
            double vf = v0 - this->mySpiral->getDirection() * this->deltav _GETVALUE;
            *this->radius_after_segment = this->myBody->mu / vf / vf;


            //Step 9: enforce continuity constraints on time, mass, electric propellant, and chemical fuel
            F[Findex++] = (this->SpiralSegmentTime - this->actual_SpiralSegmentTime) / 86400.0;
            F[Findex++] = (this->mass_after_segment - mass_after_segment_computed) / this->myOptions->maximum_mass;
            F[Findex++] = (this->virtual_electric_propellant - this->electric_propellant_used) / this->myOptions->maximum_mass;
            F[Findex++] = (this->virtual_chemical_fuel - this->chemical_fuel_used) / this->myOptions->maximum_mass;

            //Step 10: derivatives
            if (needG)
            {
                //Step 8.1: pre-requisites
                double drSundt = ((this->PositionRelativeToSun(0) * this->dPositionRelativeToSun_dt(0)
                                 + this->PositionRelativeToSun(1) * this->dPositionRelativeToSun_dt(1)
                                 + this->PositionRelativeToSun(2) * this->dPositionRelativeToSun_dt(2)) / this->PositionRelativeToSun.norm()) _GETVALUE / this->myOptions->AU;

                double dPdt = (this->mySpacecraft->getdPdr() * drSundt + this->mySpacecraft->getdPdt());
                double dIsp_dt = this->mySpacecraft->getEPdIspdP() * dPdt;
                double dmdot_dt = this->mySpacecraft->getEPdMassFlowRatedP() * dPdt;

                double d_deltamass_expfun_dt = (1000.0 * this->deltav * delta_mass_expfun * dIsp_dt / g0 / Isp / Isp)_GETVALUE;

                double d_mass_after_segment_computed_dt = this->mass_before_segment _GETVALUE * d_deltamass_expfun_dt;

                double d_electric_propellant_used_dt = -d_mass_after_segment_computed_dt;

                double d_electric_propellant_used_dmass_before_segment = (1.0 - delta_mass_expfun) _GETVALUE;
                
                //Step 10.2: derivatives of time constraint
                //wrt all previous time variables
                double dSpiralSegmentTime_dt = ((this->electric_propellant_used * dmdot_dt - mdot * d_electric_propellant_used_dt) / mdot / mdot) _GETVALUE;
                for (size_t dIndex = 0; dIndex < this->Xindices_EventLeftEpoch.size(); ++dIndex)
                {
                    size_t Gindex = this->Gindex_SpiralSegmentTime_wrt_previous_time_variables[dIndex];
                    size_t Xindex = this->Xindices_EventLeftEpoch[dIndex];

                    G[Gindex] = this->X_scale_factors->operator[](Xindex)
                        * dSpiralSegmentTime_dt
                        / 86400.0;
                }
                //wrt left-hand mass
                double dSpiralSegmentTime_dmass_before_segment = (d_electric_propellant_used_dmass_before_segment / mdot)_GETVALUE;
                G[this->Gindex_SpiralSegmentTime_wrt_left_hand_mass] = -this->X_scale_factors->operator[](this->Xindex_left_mass)
                    * dSpiralSegmentTime_dmass_before_segment
                    / 86400.0;

                //wrt segment encoded time
                G[this->Gindex_SpiralSegmentTime_wrt_SpiralSegmentTime] = this->X_scale_factors->operator[](this->Xindex_SpiralSegmentTime)
                    / 86400.0;

                //Step 10.3: derivatives of mass constraint
                //wrt all previous time variables
                for (size_t dIndex = 0; dIndex < this->Xindices_EventLeftEpoch.size(); ++dIndex)
                {
                    size_t Gindex = this->Gindex_mass_after_segment_wrt_previous_time_variables[dIndex];
                    size_t Xindex = this->Xindices_EventLeftEpoch[dIndex];

                    G[Gindex] = this->X_scale_factors->operator[](Xindex)
                        * -d_mass_after_segment_computed_dt
                        / this->myOptions->maximum_mass;
                }

                //wrt left-hand mass
                G[this->Gindex_mass_after_segment_wrt_left_hand_mass] = this->X_scale_factors->operator[](this->Xindex_left_mass)
                    * -delta_mass_expfun _GETVALUE
                    / this->myOptions->maximum_mass;

                //wrt segment encoded mass
                G[this->Gindex_mass_after_segment_wrt_mass_after_segment] = this->X_scale_factors->operator[](this->Xindex_mass_after_segment)
                    / this->myOptions->maximum_mass;

                //IF ACS IS ON wrt segment encoded time
                if (this->myOptions->trackACS)
                {
                    G[this->Gindex_mass_after_segment_wrt_SpiralSegmentTime] = this->X_scale_factors->operator[](this->Xindex_SpiralSegmentTime)
                        * this->myOptions->ACS_kg_per_day / 86400.0
                        / this->myOptions->maximum_mass;
                }

                //Step 10.4: derivatives of virtual electric propellant
                //wrt all previous time variables
                for (size_t dIndex = 0; dIndex < this->Xindices_EventLeftEpoch.size(); ++dIndex)
                {
                    size_t Gindex = this->Gindex_virtual_electric_propellant_wrt_previous_time_variables[dIndex];
                    size_t Xindex = this->Xindices_EventLeftEpoch[dIndex];

                    G[Gindex] = this->X_scale_factors->operator[](Xindex)
                        * -d_electric_propellant_used_dt
                        / this->myOptions->maximum_mass;
                }

                //wrt left-hand mass
                G[this->Gindex_virtual_electric_propellant_wrt_left_hand_mass] = this->X_scale_factors->operator[](this->Xindex_left_mass)
                    * -d_electric_propellant_used_dmass_before_segment
                    / this->myOptions->maximum_mass;

                //wrt segment virtual electric propellant
                G[this->Gindex_virtual_electric_propellant_wrt_virtual_electric_propellant] = this->X_scale_factors->operator[](this->Xindex_virtual_electric_propellant)
                    / this->myOptions->maximum_mass;

                //Step 10.5: derivatives of virtual chemical fuel
                if (this->myOptions->trackACS)
                {
                    G[this->Gindex_virtual_chemical_fuel_wrt_SpiralSegmentTime] = this->X_scale_factors->operator[](this->Xindex_SpiralSegmentTime)
                        * -this->myOptions->ACS_kg_per_day / 86400.0
                        / this->myOptions->maximum_mass;
                }

                //wrt segment virtual chemical fuel
                G[this->Gindex_virtual_chemical_fuel_wrt_virtual_chemical_fuel] = this->X_scale_factors->operator[](this->Xindex_virtual_chemical_fuel)
                    / this->myOptions->maximum_mass;
            }
        }

        void EdelbaumSpiralSegment::process_epoch_time(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            this->LaunchDate = X[this->Xindices_EventLeftEpoch[0]];

            this->segmentLeftEpoch = 0.0;
            for (size_t listIndex = 0; listIndex < this->Xindices_EventLeftEpoch.size(); ++listIndex)
            {
                this->segmentLeftEpoch += X[this->Xindices_EventLeftEpoch[listIndex]];
            }
        }//end process_epoch_time

        void EdelbaumSpiralSegment::output(std::ofstream& outputfile,
            size_t& eventcount)
        {
            std::string eventType;

            if (this->spiralSegmentIndex == 0)
                eventType = "begin_spiral";
            else
                eventType = "LT_spiral";

            doubleType temp_body_state[12];

            this->myBody->locate_body(this->segmentLeftEpoch,
                temp_body_state,
                false,
                *this->myOptions);

            for (size_t stateIndex = 0; stateIndex < 6; ++stateIndex)
                this->StateRelativeToJourneyCentralBody(stateIndex) = temp_body_state[stateIndex];
            this->StateRelativeToJourneyCentralBody(6) = this->mass_before_segment;
            this->StateRelativeToJourneyCentralBody(7) = this->segmentLeftEpoch;

            if (this->isHeliocentric)
            {
                // get the position vector of the sun
                doubleType temp_state[12];
                myUniverse->locate_central_body(this->segmentLeftEpoch, temp_state, *this->myOptions, false);
                // locate_central_body gives us the state of the central body relative to the sun 
                // we want the state of the sun relative to the central body though (so its negation)
                for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
                {
                    this->PositionCB_Sun(stateIndex) = temp_state[stateIndex];
                    this->PositionRelativeToSun(stateIndex) = this->PositionCB_Sun(stateIndex) + this->StateRelativeToJourneyCentralBody(stateIndex);
                }
            }
            else
            {
                // if we are dealing with a heliocentric trajectory, 
                // then the spacecraft's position relative to the sun is the same as its position relative to the central body
                for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
                {
                    this->PositionRelativeToSun(stateIndex) = this->StateRelativeToJourneyCentralBody(stateIndex);
                }
            } // end if isHeliocentric
            doubleType DistanceFromSunAU = this->PositionRelativeToSun.norm() / this->myOptions->AU;

            this->mySpacecraft->computePowerState(DistanceFromSunAU, this->segmentLeftEpoch);

            this->mySpacecraft->computeElectricPropulsionPerformance(this->myDutyCycle);

            this->write_output_line(outputfile,
                eventcount,
                eventType,
                this->myBody->name,
                this->SpiralSegmentTime / 86400.0,
                *this->radius_before_segment - this->myBody->radius,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                this->StateRelativeToJourneyCentralBody,
                math::Matrix<doubleType>(3, 1, 0.0),
                math::Matrix<doubleType>(3, 1, 0.0),
                this->deltav,
                this->mySpacecraft->getEPthrust(),
                this->mySpacecraft->getEPIsp(),
                this->mySpacecraft->getAvailablePower(),
                this->mySpacecraft->getEPMassFlowRate(),
                this->mySpacecraft->getEPNumberOfActiveThrusters(),
                this->mySpacecraft->getEPActivePower(),
                this->mySpacecraft->getEPThrottleLevelString());


            //if this is the last segment, we need to print the end state
            if (this->spiralSegmentIndex == this->myOptions->spiral_segments - 1)
            {
                eventType = "end_spiral";

                this->myBody->locate_body(this->segmentLeftEpoch + this->SpiralSegmentTime,
                    temp_body_state,
                    false,
                    *this->myOptions);

                for (size_t stateIndex = 0; stateIndex < 6; ++stateIndex)
                    this->StateRelativeToJourneyCentralBody(stateIndex) = temp_body_state[stateIndex];
                this->StateRelativeToJourneyCentralBody(6) = this->mass_after_segment;
                this->StateRelativeToJourneyCentralBody(7) = this->segmentLeftEpoch + this->SpiralSegmentTime;
                
                if (this->isHeliocentric)
                {
                    // get the position vector of the sun
                    doubleType temp_state[12];
                    myUniverse->locate_central_body(this->segmentLeftEpoch + this->SpiralSegmentTime, temp_state, *this->myOptions, false);
                    // locate_central_body gives us the state of the central body relative to the sun 
                    // we want the state of the sun relative to the central body though (so its negation)
                    for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
                    {
                        this->PositionCB_Sun(stateIndex) = temp_state[stateIndex];
                        this->PositionRelativeToSun(stateIndex) = this->PositionCB_Sun(stateIndex) + this->StateRelativeToJourneyCentralBody(stateIndex);
                    }
                }
                else
                {
                    // if we are dealing with a heliocentric trajectory, 
                    // then the spacecraft's position relative to the sun is the same as its position relative to the central body
                    for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
                    {
                        this->PositionRelativeToSun(stateIndex) = this->StateRelativeToJourneyCentralBody(stateIndex);
                    }
                } // end if isHeliocentric
                DistanceFromSunAU = this->PositionRelativeToSun.norm() / this->myOptions->AU;

                this->mySpacecraft->computePowerState(DistanceFromSunAU, this->segmentLeftEpoch);

                this->mySpacecraft->computeElectricPropulsionPerformance(this->myDutyCycle);

                this->write_output_line(outputfile,
                    eventcount,
                    eventType,
                    this->myBody->name,
                    0.0,
                    *this->radius_after_segment - this->myBody->radius,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    this->StateRelativeToJourneyCentralBody,
                    math::Matrix<doubleType>(3, 1, 0.0),
                    math::Matrix<doubleType>(3, 1, 0.0),
                    0.0,
                    this->mySpacecraft->getEPthrust(),
                    this->mySpacecraft->getEPIsp(),
                    this->mySpacecraft->getAvailablePower(),
                    this->mySpacecraft->getEPMassFlowRate(),
                    this->mySpacecraft->getEPNumberOfActiveThrusters(),
                    this->mySpacecraft->getEPActivePower(),
                    this->mySpacecraft->getEPThrottleLevelString());
            }
        }
    }//end namespace BoundaryEvents
}//end namespace EMTG