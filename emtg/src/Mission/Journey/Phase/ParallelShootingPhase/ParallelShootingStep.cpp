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

//base parallel shooting step for EMTGv9
//Jacob Englander 2-20-2018

#include "ParallelShootingStep.h"
#include "ParallelShootingStep_maneuver_constraint_factory.h"

namespace EMTG
{
    namespace Phases
    {
        ParallelShootingStep::ParallelShootingStep() : 
            StepDutyCycle(1.0),
            stepDeltav(0.0),
            num_interior_control_points(1)
        {
            this->StateStepLeftInertial.resize(10, 1, 0.0);
            this->StateStepRightInertial.resize(10, 1, 0.0);
        }//end default constructor

        ParallelShootingStep::ParallelShootingStep(const std::string& name,
            const size_t& journeyIndex,
            const size_t& phaseIndex,
            const size_t& stepIndex,
            const size_t& stageIndex,
            ParallelShootingStep* previousStep,
            ParallelShootingPhase* myPhase,
            Astrodynamics::universe* myUniverse,
            HardwareModels::Spacecraft* mySpacecraft,
            missionoptions* myOptions)
        {
            this->initialize(name,
                journeyIndex,
                phaseIndex,
                stepIndex,
                stageIndex,
                previousStep,
                myPhase,
                myUniverse,
                mySpacecraft,
                myOptions);
        }//end constructor

        void ParallelShootingStep::initialize(const std::string& name,
            const size_t& journeyIndex,
            const size_t& phaseIndex,
            const size_t& stepIndex,
            const size_t& stageIndex,
            ParallelShootingStep* previousStep,
            ParallelShootingPhase* myPhase,
            Astrodynamics::universe* myUniverse,
            HardwareModels::Spacecraft* mySpacecraft,
            missionoptions* myOptions)
        {
            this->name = name;
            this->journeyIndex = journeyIndex;
            this->phaseIndex = phaseIndex;
            this->stepIndex = stepIndex;
            this->previousStep = previousStep;
            this->myPhase = myPhase;
            this->mySpacecraft = mySpacecraft;
            this->myJourneyOptions = &myOptions->Journeys[this->journeyIndex];

            //initialize the writey thing
            this->writey_thing::initialize(myOptions, myUniverse);

            this->prefix = this->name + ": ";

            this->num_interior_control_points = this->myJourneyOptions->num_interior_control_points;

            if (this->myJourneyOptions->override_integration_step_size)
            {
                this->EphemerisOutputResolution = this->myJourneyOptions->integration_step_size;
            }
            else
            {
                this->EphemerisOutputResolution = this->myOptions->integration_time_step_size;
            }

            this->numMatchConstraints = this->myPhase->get_numMatchConstraints();

            //control
            SpacecraftThrusterMode myThrusterMode = this->mySpacecraft->getCurrentStageOptions().getElectricPropulsionSystemOptions().getThrusterMode();
            if (myThrusterMode == FixedEfficiencyVSI
                || myThrusterMode == Stepped2D
                || myThrusterMode == Poly2D)
            {
                this->isVSI = true;
                this->num_controls = 4;
            }
            else
            {
                this->isVSI = false;
                this->num_controls = 3;
            }

            //size state
            this->StateAfterSubStepInertial.resize(this->num_interior_control_points, math::Matrix<doubleType>(10, 1, 0.0));

            //size various things
            this->throttle.resize(this->num_interior_control_points);
            this->ControlVector.resize(this->num_interior_control_points, math::Matrix<doubleType>(this->num_controls, 1, 0.0));

            //distance contraints
            std::string shortprefix = "p" + std::to_string(this->phaseIndex);
            for (std::string& constraintDefinition : this->myJourneyOptions->PhaseDistanceConstraintDefinitions)
            {
                if (constraintDefinition.find("#") != 0 
                    && constraintDefinition.find(shortprefix) < 1024)
                {
                    this->myDistanceConstraints.push_back(ParallelShootingStepDistanceConstraint(constraintDefinition,
                        this->journeyIndex,
                        this->phaseIndex,
                        this->stepIndex,
                        this,
                        this->myOptions,
                        this->myUniverse));
                }
            }


            //duty cycle, may be overriden in "maneuver constraints" below
            this->StepDutyCycle = this->myPhase->getPhaseDutyCycle();

            //maneuver constraints
            std::string burnprefix = "b" + std::to_string(this->stepIndex);
            for (std::string& constraintDefinition : this->myJourneyOptions->ManeuverConstraintDefinitions)
            {
                if (constraintDefinition.find("#") != 0) //don't create a constraint if it is commented out
                {
                    if (constraintDefinition.find(shortprefix) < 1024)
                    {
                        if (boost::to_lower_copy(constraintDefinition).find("dutycycle") < 1024)                       
                        {
                            if (boost::to_lower_copy(constraintDefinition).find(burnprefix) < 1024)
                            {
                                std::vector<std::string> ConstraintDefinitionCell;
                                boost::split(ConstraintDefinitionCell, constraintDefinition, boost::is_any_of("_"), boost::token_compress_on);

                                this->StepDutyCycle = std::stod(ConstraintDefinitionCell[2]);
                            }
                        }
                        else
                        {
                            for (size_t subStepIndex = 0; subStepIndex < this->num_interior_control_points; ++subStepIndex)
                            {
                                this->myManeuverConstraints.push_back(create_ParallelShootingStep_maneuver_constraint(this->name,
                                    this->journeyIndex,
                                    this->phaseIndex,
                                    this->stepIndex,
                                    subStepIndex,
                                    this->stageIndex,
                                    this,
                                    this->myUniverse,
                                    this->mySpacecraft,
                                    this->myOptions,
                                    constraintDefinition));
                            }//end call to factory for this interior point
                        }//end call to factory
                    }//end this constraint is in this phase
                }//end comment check
            }//end maneuver constraints

            //time derivative
            this->dStepTime_dPhaseFlightTime = 1.0 / this->myPhase->get_num_steps();

            //truth tables
            this->TruthTable_StateStepRightInertial_wrt_StateStepLeftInertial.resize(10, std::vector<size_t>(10, true));
            this->TruthTable_StateStepRightInertial_wrt_Control.resize(10, true);
            this->TruthTable_StateStepRightInertial_wrt_PhaseFlightTime.resize(10, true);

            //epoch has derivatives with respect only to itself
            this->TruthTable_StateStepRightInertial_wrt_StateStepLeftInertial[7] = std::vector<size_t>(10, false);
            this->TruthTable_StateStepRightInertial_wrt_StateStepLeftInertial[7][7] = true;

            //position, velocity, and mass have no derivative with respect to either tank variable
            for (size_t stateIndex = 0; stateIndex < 7; ++stateIndex)
            {
                this->TruthTable_StateStepRightInertial_wrt_StateStepLeftInertial[stateIndex][8] = false;
                this->TruthTable_StateStepRightInertial_wrt_StateStepLeftInertial[stateIndex][9] = false;
            }

            //chemical fuel only has a derivative with respect to itself and epoch, and even then only if ACS is on
            this->TruthTable_StateStepRightInertial_wrt_StateStepLeftInertial[8] = std::vector<size_t>(10, false);
            this->TruthTable_StateStepRightInertial_wrt_StateStepLeftInertial[8][8] = true;
            if (this->myOptions->trackACS)
                this->TruthTable_StateStepRightInertial_wrt_StateStepLeftInertial[8][7] = true;

            //position, velocity, mass, and electric propellant have derivatives with respect to control but chemical fuel does NOT
            this->TruthTable_StateStepRightInertial_wrt_Control[8] = false;

            //electric propellant does not have a derivative with respect to chemical fuel
            this->TruthTable_StateStepRightInertial_wrt_StateStepLeftInertial[9][8] = false;
        }//end initialize

        void ParallelShootingStep::setup_calcbounds(
            std::vector<double>* Xupperbounds,
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
            
            for (size_t constraintIndex = 0; constraintIndex < this->myDistanceConstraints.size(); ++constraintIndex)
                this->myDistanceConstraints[constraintIndex].setup_calcbounds(Xupperbounds,
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

            for (size_t constraintIndex = 0; constraintIndex < this->myManeuverConstraints.size(); ++constraintIndex)
                this->myManeuverConstraints[constraintIndex].setup_calcbounds(Xupperbounds,
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

        void ParallelShootingStep::output_STM(size_t& STMindex)
        {
            for (size_t subStepIndex = 0; subStepIndex < this->num_interior_control_points; ++subStepIndex)
                this->AugmentedSTM[subStepIndex].print_to_file(this->myOptions->working_directory + "//" + boost::replace_all_copy(this->prefix, ": ", "") + "_STM_" + std::to_string(STMindex++) + "_forward.stm");
        }//end output_STM

        //calcbounds
        void ParallelShootingStep::calcbounds_step_left_state()
        {
            if (this->myOptions->ParallelShootingStateRepresentation == StateRepresentation::SphericalRADEC)
            {
                //Step 1: position variables
                //Step 1.1: radius
                Xlowerbounds->push_back(this->myUniverse->central_body.radius + this->myUniverse->central_body.minimum_safe_flyby_altitude);
                Xupperbounds->push_back(this->myUniverse->r_SOI);
                X_scale_factors->push_back(this->myUniverse->LU);
                Xdescriptions->push_back(prefix + "left state r");
                this->Xindex_rMag = this->Xdescriptions->size() - 1;
                for (size_t stateIndex : {0, 1, 2})
                {
                    this->Derivatives_of_StateStepLeftInertial.push_back(std::make_tuple(this->Xdescriptions->size() - 1, stateIndex, 1.0));
                    this->dIndex_StateStepLeftInertial_wrt_r.push_back(this->Derivatives_of_StateStepLeftInertial.size() - 1);
                }

                //Step 1.2: RA
                this->Xlowerbounds->push_back(-8.0 * math::PI);
                this->Xupperbounds->push_back(8.0 * math::PI);
                this->X_scale_factors->push_back(1.0);
                this->Xdescriptions->push_back(this->prefix + "left state RA");
                this->Xindex_RA = this->Xdescriptions->size() - 1;
                //affects x, y
                for (size_t stateIndex : {0, 1})
                {
                    this->Derivatives_of_StateStepLeftInertial.push_back(std::make_tuple(this->Xdescriptions->size() - 1, stateIndex, 1.0));
                    this->dIndex_StateStepLeftInertial_wrt_RA.push_back(this->Derivatives_of_StateStepLeftInertial.size() - 1);
                }

                //Step 1.3: DEC
                this->Xlowerbounds->push_back(-0.5 * math::PI);
                this->Xupperbounds->push_back(0.5 * math::PI);
                this->X_scale_factors->push_back(1.0);
                this->Xdescriptions->push_back(this->prefix + "left state DEC");
                this->Xindex_DEC = this->Xdescriptions->size() - 1;
                for (size_t stateIndex : {0, 1, 2})
                {
                    this->Derivatives_of_StateStepLeftInertial.push_back(std::make_tuple(this->Xdescriptions->size() - 1, stateIndex, 1.0));
                    this->dIndex_StateStepLeftInertial_wrt_DEC.push_back(this->Derivatives_of_StateStepLeftInertial.size() - 1);
                }

                //Step 2: velocity variables
                //Step 2.1: v
                this->Xlowerbounds->push_back(0.0);
                this->Xupperbounds->push_back(10.0 * this->myUniverse->LU / this->myUniverse->TU);
                this->X_scale_factors->push_back(this->myUniverse->LU / this->myUniverse->TU);
                this->Xdescriptions->push_back(this->prefix + "left state v");
                this->Xindex_vMag = this->Xdescriptions->size() - 1;
                for (size_t stateIndex : {3, 4, 5})
                {
                    this->Derivatives_of_StateStepLeftInertial.push_back(std::make_tuple(this->Xdescriptions->size() - 1, stateIndex, 1.0));
                    this->dIndex_StateStepLeftInertial_wrt_v.push_back(this->Derivatives_of_StateStepLeftInertial.size() - 1);
                }


                //Step 2.2: velocity RA 
                this->Xlowerbounds->push_back(-8.0 * math::PI);
                this->Xupperbounds->push_back(8.0 * math::PI);
                this->X_scale_factors->push_back(1.0);
                this->Xdescriptions->push_back(this->prefix + "left state vRA");
                this->Xindex_vRA = this->Xdescriptions->size() - 1;
                for (size_t stateIndex : {3, 4})
                {
                    this->Derivatives_of_StateStepLeftInertial.push_back(std::make_tuple(this->Xdescriptions->size() - 1, stateIndex, 1.0));
                    this->dIndex_StateStepLeftInertial_wrt_vRA.push_back(this->Derivatives_of_StateStepLeftInertial.size() - 1);
                }

                //Step 2.3: velocity DEC 
                this->Xlowerbounds->push_back(-8.0 * math::PI);
                this->Xupperbounds->push_back(8.0 * math::PI);
                this->X_scale_factors->push_back(1.0);
                this->Xdescriptions->push_back(this->prefix + "left state vDEC");
                this->Xindex_vDEC = this->Xdescriptions->size() - 1;
                for (size_t stateIndex : {3, 4, 5})
                {
                    this->Derivatives_of_StateStepLeftInertial.push_back(std::make_tuple(this->Xdescriptions->size() - 1, stateIndex, 1.0));
                    this->dIndex_StateStepLeftInertial_wrt_vDEC.push_back(this->Derivatives_of_StateStepLeftInertial.size() - 1);
                }
            }
            else if (this->myOptions->ParallelShootingStateRepresentation == StateRepresentation::SphericalAZFPA)
            {
                //Step 1: position variables
                //Step 1.1: radius
                Xlowerbounds->push_back(this->myUniverse->central_body.radius + this->myUniverse->central_body.minimum_safe_flyby_altitude);
                Xupperbounds->push_back(this->myUniverse->r_SOI);
                X_scale_factors->push_back(this->myUniverse->LU);
                Xdescriptions->push_back(prefix + "left state r");
                this->Xindex_rMag = this->Xdescriptions->size() - 1;
                for (size_t stateIndex : {0, 1, 2})
                {
                    this->Derivatives_of_StateStepLeftInertial.push_back(std::make_tuple(this->Xdescriptions->size() - 1, stateIndex, 1.0));
                    this->dIndex_StateStepLeftInertial_wrt_r.push_back(this->Derivatives_of_StateStepLeftInertial.size() - 1);
                }

                //Step 1.2: RA
                this->Xlowerbounds->push_back(-8.0 * math::PI);
                this->Xupperbounds->push_back(8.0 * math::PI);
                this->X_scale_factors->push_back(1.0);
                this->Xdescriptions->push_back(this->prefix + "left state RA");
                this->Xindex_RA = this->Xdescriptions->size() - 1;
                //affects x, y
                for (size_t stateIndex : {0, 1, 3, 4})
                {
                    this->Derivatives_of_StateStepLeftInertial.push_back(std::make_tuple(this->Xdescriptions->size() - 1, stateIndex, 1.0));
                    this->dIndex_StateStepLeftInertial_wrt_RA.push_back(this->Derivatives_of_StateStepLeftInertial.size() - 1);
                }

                //Step 1.3: DEC
                this->Xlowerbounds->push_back(-0.5 * math::PI);
                this->Xupperbounds->push_back(0.5 * math::PI);
                this->X_scale_factors->push_back(1.0);
                this->Xdescriptions->push_back(this->prefix + "left state DEC");
                this->Xindex_DEC = this->Xdescriptions->size() - 1;
                for (size_t stateIndex : {0, 1, 2, 3, 4, 5})
                {
                    this->Derivatives_of_StateStepLeftInertial.push_back(std::make_tuple(this->Xdescriptions->size() - 1, stateIndex, 1.0));
                    this->dIndex_StateStepLeftInertial_wrt_DEC.push_back(this->Derivatives_of_StateStepLeftInertial.size() - 1);
                }

                //Step 2: velocity variables
                //Step 2.1: v
                this->Xlowerbounds->push_back(0.0);
                this->Xupperbounds->push_back(10.0 * this->myUniverse->LU / this->myUniverse->TU);
                this->X_scale_factors->push_back(this->myUniverse->LU / this->myUniverse->TU);
                this->Xdescriptions->push_back(this->prefix + "left state v");
                this->Xindex_vMag = this->Xdescriptions->size() - 1;
                for (size_t stateIndex : {3, 4, 5})
                {
                    this->Derivatives_of_StateStepLeftInertial.push_back(std::make_tuple(this->Xdescriptions->size() - 1, stateIndex, 1.0));
                    this->dIndex_StateStepLeftInertial_wrt_v.push_back(this->Derivatives_of_StateStepLeftInertial.size() - 1);
                }


                //Step 2.2: velocity AZ 
                this->Xlowerbounds->push_back(-8.0 * math::PI);
                this->Xupperbounds->push_back(8.0 * math::PI);
                this->X_scale_factors->push_back(1.0);
                this->Xdescriptions->push_back(this->prefix + "left state AZ");
                this->Xindex_AZ = this->Xdescriptions->size() - 1;
                for (size_t stateIndex : {3, 4, 5})
                {
                    this->Derivatives_of_StateStepLeftInertial.push_back(std::make_tuple(this->Xdescriptions->size() - 1, stateIndex, 1.0));
                    this->dIndex_StateStepLeftInertial_wrt_AZ.push_back(this->Derivatives_of_StateStepLeftInertial.size() - 1);
                }

                //Step 2.3: velocity DEC 
                this->Xlowerbounds->push_back(-8.0 * math::PI);
                this->Xupperbounds->push_back(8.0 * math::PI);
                this->X_scale_factors->push_back(1.0);
                this->Xdescriptions->push_back(this->prefix + "left state FPA");
                this->Xindex_FPA = this->Xdescriptions->size() - 1;
                for (size_t stateIndex : {3, 4, 5})
                {
                    this->Derivatives_of_StateStepLeftInertial.push_back(std::make_tuple(this->Xdescriptions->size() - 1, stateIndex, 1.0));
                    this->dIndex_StateStepLeftInertial_wrt_FPA.push_back(this->Derivatives_of_StateStepLeftInertial.size() - 1);
                }
            }

            //Step 3: mass variable
            this->Xlowerbounds->push_back(math::SMALL);
            this->Xupperbounds->push_back(this->myJourneyOptions->maximum_mass);
            this->X_scale_factors->push_back(1.0 / this->myUniverse->continuity_constraint_scale_factors(6));
            this->Xdescriptions->push_back(this->prefix + "left state mass");
            this->Xindex_mass = this->Xdescriptions->size() - 1;
            this->Derivatives_of_StateStepLeftInertial.push_back(std::make_tuple(this->Xindex_mass, 6, 1.0));
            this->dIndex_mass_wrt_mass = this->Derivatives_of_StateStepLeftInertial.size() - 1;

            //Step 4: epoch
            this->calculate_dependencies_epoch_time();

            //Step 5: virtual chemical fuel
            this->Xlowerbounds->push_back(0.0);
            this->Xupperbounds->push_back(this->myJourneyOptions->maximum_mass);
            this->X_scale_factors->push_back(this->myJourneyOptions->maximum_mass);
            this->Xdescriptions->push_back(prefix + "virtual chemical fuel");
            this->Xindex_chemical_fuel = this->Xdescriptions->size() - 1;
            this->Derivatives_of_StateStepLeftInertial.push_back(std::make_tuple(this->Xindex_chemical_fuel, 8, 1.0));
            this->dIndex_chemicalFuel_wrt_chemicalFuel = this->Derivatives_of_StateStepLeftInertial.size() - 1;

            //Step 6: virtual electric propellant
            this->Xlowerbounds->push_back(0.0);
            this->Xupperbounds->push_back(this->myJourneyOptions->maximum_mass);
            this->X_scale_factors->push_back(this->myJourneyOptions->maximum_mass);
            this->Xdescriptions->push_back(prefix + "virtual electric propellant");
            this->Xindex_electric_propellant = this->Xdescriptions->size() - 1;
            this->Derivatives_of_StateStepLeftInertial.push_back(std::make_tuple(this->Xindex_electric_propellant, 9, 1.0));
            this->dIndex_electricPropellant_wrt_electricPropellant = this->Derivatives_of_StateStepLeftInertial.size() - 1;
        }//end calcbounds_step_left_state()
        
        void ParallelShootingStep::calculate_dependencies_epoch_time()
        {
            //the left epoch of this event depends on all epoch and time variables before it
            for (size_t Xindex = 0; Xindex < this->Xdescriptions->size(); ++Xindex)
            {
                if (this->Xdescriptions->at(Xindex).find("epoch") < 1024
                    || this->Xdescriptions->at(Xindex).find("time") < 1024)
                {
                    this->Xindices_EventLeftEpoch.push_back(Xindex);
                    this->Derivatives_of_StateStepLeftInertial_wrt_Time.push_back({ Xindex, 7, 1.0 });
                }
            }

            //phase flight time is always the last entry and everything has a derivative with respect to it
            this->Xindex_PhaseFlightTime = this->Xindices_EventLeftEpoch.back();
        }//end calculate_dependencies_left_epoch

        void ParallelShootingStep::calcbounds_step_control()
        {
            std::vector<std::string> controlNames({ "x", "y", "z", "command" });
            this->dIndex_right_state_wrt_control.resize(this->num_interior_control_points, std::vector< std::vector<size_t> >(10));
            this->Xindices_control.resize(this->num_interior_control_points, std::vector<size_t>(this->num_controls));
            this->G_indices_control_magnitude_constraint_wrt_Control.resize(this->num_interior_control_points);

            //loop over substeps

            for (size_t subStepIndex = 0; subStepIndex < this->num_interior_control_points; ++subStepIndex)
            {
                for (size_t controlIndex = 0; controlIndex < 3; ++controlIndex)
                {
                    this->Xlowerbounds->push_back(-1.0);
                    this->Xupperbounds->push_back(1.0);
                    this->X_scale_factors->push_back(1.0);
                    this->Xdescriptions->push_back(prefix + "substep" + std::to_string(subStepIndex) + " u_" + controlNames[controlIndex]);
                    this->Xindices_control[subStepIndex][controlIndex] = this->Xdescriptions->size() - 1;

                    for (size_t stateIndex : {0, 1, 2, 3, 4, 5, 6, 9})//no affect on epoch or chemical fuel
                    {
                        this->Derivatives_of_StateStepRightInertial.push_back({ this->Xdescriptions->size() - 1, stateIndex, 1.0 });
                        this->dIndex_right_state_wrt_control[subStepIndex][stateIndex].push_back(this->Derivatives_of_StateStepRightInertial.size() - 1);
                    }
                }

                //u_command
                if (this->isVSI)
                {
                    this->Xlowerbounds->push_back(0.0);
                    this->Xupperbounds->push_back(1.0);
                    this->X_scale_factors->push_back(1.0);
                    this->Xdescriptions->push_back(prefix + "u_" + controlNames[3]);
                    this->Xindices_control[subStepIndex][3] = this->Xdescriptions->size() - 1;

                    for (size_t stateIndex : {0, 1, 2, 3, 4, 5, 6, 9})//no affect on epoch or chemical fuel
                    {
                        this->Derivatives_of_StateStepRightInertial.push_back({ this->Xdescriptions->size() - 1, stateIndex, 1.0 });
                        this->dIndex_right_state_wrt_control[subStepIndex][stateIndex].push_back(this->Derivatives_of_StateStepRightInertial.size() - 1);
                    }
                }

                //control magnitude constraint
                if (this->myJourneyOptions->force_unit_magnitude_control == ControlMagnitudeType::UnitMagnitude)
                {
                    this->Flowerbounds->push_back(1.0);
                }
                else
                {
                    this->Flowerbounds->push_back(0.0);
                }
                if (this->myJourneyOptions->force_unit_magnitude_control == ControlMagnitudeType::ZeroMagnitude)
                {
                    this->Fupperbounds->push_back(0.0);
                }
                else
                {
                    this->Fupperbounds->push_back(1.0);
                }
                this->Fdescriptions->push_back(prefix + "substep " + std::to_string(subStepIndex) + " control magnitude");

                //sparsity pattern
                for (size_t controlIndex = 0; controlIndex < this->num_controls; ++controlIndex)
                {
                    this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                        this->Xindices_control[subStepIndex][controlIndex],
                        this->G_indices_control_magnitude_constraint_wrt_Control[subStepIndex]);
                }
            }
        }//end calcbounds_step_control
        
        void ParallelShootingStep::calcbounds_step_main()
        {
            //ListOfVariablesAffectingCurrentStepLeftState
            //DerivativesOfStateLeftInertialByVariable[varIndex][stateIndex] contains a dIndex associated with that stateIndex's row in Derivatives_of_StateStepLeftInertial
            //Step 1: construct the right-side derivative tuples
            this->dIndex_StateStepRightInertial_wrt_LeftStateVariables.resize(10);
            this->dIndex_StateStepRightInertial_wrt_PreviousTimeVariables.resize(10);

            //for each varIndex that affects the left-hand state
            //  for each leftStateIndex that is affected by varIndex
            //    for each rightStateIndex
            //      if rightStateIndex is affected by leftStateIndex
            //        create a derivative entry of rightStateIndex wrt varIndex
            //        store the dIndex
            //Step 1: non-time
            this->DerivativesOfCurrentStepRightStateByVariable.resize(this->ListOfVariablesAffectingCurrentStepLeftState.size());
            for (size_t varIndex = 0; varIndex < this->ListOfVariablesAffectingCurrentStepLeftState.size(); ++varIndex)
            {
                size_t Xindex = this->ListOfVariablesAffectingCurrentStepLeftState[varIndex];

                for (size_t entryIndex = 0; entryIndex < this->DerivativesOfCurrentStepLeftStateByVariable[varIndex].size(); ++entryIndex)
                {
                    size_t leftStateIndex = std::get<0>(this->DerivativesOfCurrentStepLeftStateByVariable[varIndex][entryIndex]);

                    for (size_t rightStateIndex = 0; rightStateIndex < 10; ++rightStateIndex)
                    {
                        if (this->TruthTable_StateStepRightInertial_wrt_StateStepLeftInertial[rightStateIndex][leftStateIndex])
                        {
                            bool duplicateEntry = false;
                            for (std::tuple<size_t, size_t> OtherDerivativeIndexEntry : this->DerivativesOfCurrentStepRightStateByVariable[varIndex])
                            {
                                if (std::get<0>(OtherDerivativeIndexEntry) == rightStateIndex)
                                {
                                    duplicateEntry = true;
                                    break;
                                }
                            }
                            if (!duplicateEntry)
                            {
                                this->Derivatives_of_StateStepRightInertial.push_back({ Xindex, rightStateIndex, 0.0 });
                                this->DerivativesOfCurrentStepRightStateByVariable[varIndex].push_back({ rightStateIndex, this->Derivatives_of_StateStepRightInertial.size() - 1 });
                            }
                        }
                    }
                }
            }

            //Step 2: previous time (i.e. not propagation time)
            this->DerivativesOfCurrentStepRightStateByTimeVariable.resize(this->ListOfTimeVariablesAffectingCurrentStepLeftState.size());
            for (size_t varIndex = 0; varIndex < this->ListOfTimeVariablesAffectingCurrentStepLeftState.size(); ++varIndex)
            {
                size_t Xindex = this->ListOfTimeVariablesAffectingCurrentStepLeftState[varIndex];

                for (size_t entryIndex = 0; entryIndex < this->DerivativesOfCurrentStepLeftStateByTimeVariable[varIndex].size(); ++entryIndex)
                {
                    size_t leftStateIndex = std::get<0>(this->DerivativesOfCurrentStepLeftStateByTimeVariable[varIndex][entryIndex]);

                    for (size_t rightStateIndex = 0; rightStateIndex < 10; ++rightStateIndex)
                    {
                        if (this->TruthTable_StateStepRightInertial_wrt_StateStepLeftInertial[rightStateIndex][leftStateIndex])
                        {
                            bool duplicateEntry = false;
                            for (std::tuple<size_t, size_t> OtherDerivativeIndexEntry : this->DerivativesOfCurrentStepRightStateByTimeVariable[varIndex])
                            {
                                if (std::get<0>(OtherDerivativeIndexEntry) == rightStateIndex)
                                {
                                    duplicateEntry = true;
                                    break;
                                }
                            }
                            if (!duplicateEntry)
                            {
                                this->Derivatives_of_StateStepRightInertial_wrt_Time.push_back({ Xindex, rightStateIndex, 0.0 });
                                this->DerivativesOfCurrentStepRightStateByTimeVariable[varIndex].push_back({ rightStateIndex, this->Derivatives_of_StateStepRightInertial_wrt_Time.size() - 1 });
                            }
                        }
                    }
                }
            }
        }//end calcbounds_step_main

        void ParallelShootingStep::calcbounds_step_left_match_point_constraints()
        {
            //the left match point constraint is defined as current_left_state - previous_right_state
            //this is the generic version for an intermediate step - the first step in the phase overrides this method

            //Step 1: derivatives with respect to previous step right state
            std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_PreviousStepRightState = this->previousStep->get_Derivatives_of_StateStepRightInertial(); //[stateIndex][dIndex][Xindex, derivative value]
            std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_PreviousStepRightState_wrt_Time = this->previousStep->get_Derivatives_of_StateStepRightInertial_wrt_Time();//[stateIndex][dIndex][Xindex, derivative value]

            //Step 1.1: non-time                                                                                                                                                                                  //Step 1.1: non-time                                                                                                                                                                        ////Step 1.1: state variables
            {
                for (size_t stateIndex = 0; stateIndex < 10; ++stateIndex)
                {
                    for (size_t dIndex = 0; dIndex < Derivatives_of_PreviousStepRightState.size(); ++dIndex) //loop over derivative entries for this state
                    {
                        size_t previousStateIndex = std::get<1>(Derivatives_of_PreviousStepRightState[dIndex]);

                        if (stateIndex == previousStateIndex)
                        {
                            size_t Xindex = std::get<0>(Derivatives_of_PreviousStepRightState[dIndex]);

                            //if this variable is not yet in the list of variables that affect the previous step's right state, add it to the list
                            if (std::find(this->ListOfVariablesAffectingPreviousStepRightState.begin(), this->ListOfVariablesAffectingPreviousStepRightState.end(), Xindex)
                                == this->ListOfVariablesAffectingPreviousStepRightState.end())
                                this->ListOfVariablesAffectingPreviousStepRightState.push_back(Xindex);
                        }
                    }
                }

                this->DerivativesOfPreviousStepRightStateByVariable.resize(ListOfVariablesAffectingPreviousStepRightState.size());

                for (size_t listIndex = 0; listIndex < this->ListOfVariablesAffectingPreviousStepRightState.size(); ++listIndex)
                {
                    for (size_t stateIndex = 0; stateIndex < 10; ++stateIndex)
                    {
                        for (size_t dIndex = 0; dIndex < Derivatives_of_PreviousStepRightState.size(); ++dIndex) //loop over derivative entries for this state
                        {
                            size_t previousStateIndex = std::get<1>(Derivatives_of_PreviousStepRightState[dIndex]);

                            if (stateIndex == previousStateIndex)
                            {
                                size_t Xindex = std::get<0>(Derivatives_of_PreviousStepRightState[dIndex]);

                                if (Xindex == this->ListOfVariablesAffectingPreviousStepRightState[listIndex])
                                {
                                    this->DerivativesOfPreviousStepRightStateByVariable[listIndex].push_back({ stateIndex, dIndex });
                                }
                            }
                        }
                    }
                }
            }//end non-time

            //Step 1.2: time variables
            {
                for (size_t stateIndex = 0; stateIndex < 10; ++stateIndex)
                {
                    for (size_t dIndex = 0; dIndex < Derivatives_of_PreviousStepRightState_wrt_Time.size(); ++dIndex) //loop over derivative entries for this state
                    {
                        size_t previousStateIndex = std::get<1>(Derivatives_of_PreviousStepRightState_wrt_Time[dIndex]);

                        if (stateIndex == previousStateIndex)
                        {
                            size_t Xindex = std::get<0>(Derivatives_of_PreviousStepRightState_wrt_Time[dIndex]);

                            //if this variable is not yet in the list of variables that affect the previous step's right state, add it to the list
                            if (std::find(this->ListOfTimeVariablesAffectingPreviousStepRightState.begin(), this->ListOfTimeVariablesAffectingPreviousStepRightState.end(), Xindex)
                                == this->ListOfTimeVariablesAffectingPreviousStepRightState.end())
                                this->ListOfTimeVariablesAffectingPreviousStepRightState.push_back(Xindex);
                        }
                    }
                }

                this->DerivativesOfPreviousStepRightStateByTimeVariable.resize(ListOfTimeVariablesAffectingPreviousStepRightState.size());

                for (size_t listIndex = 0; listIndex < this->ListOfTimeVariablesAffectingPreviousStepRightState.size(); ++listIndex)
                {
                    for (size_t stateIndex = 0; stateIndex < 10; ++stateIndex)
                    {
                        for (size_t dIndex = 0; dIndex < Derivatives_of_PreviousStepRightState_wrt_Time.size(); ++dIndex) //loop over derivative entries for this state
                        {
                            size_t previousStateIndex = std::get<1>(Derivatives_of_PreviousStepRightState_wrt_Time[dIndex]);

                            if (stateIndex == previousStateIndex)
                            {
                                size_t Xindex = std::get<0>(Derivatives_of_PreviousStepRightState_wrt_Time[dIndex]);

                                if (Xindex == this->ListOfTimeVariablesAffectingPreviousStepRightState[listIndex])
                                {
                                    this->DerivativesOfPreviousStepRightStateByTimeVariable[listIndex].push_back({ stateIndex, dIndex });
                                }
                            }
                        }
                    }
                }
            }//end time variables

            //Step 2: derivatives with respect to current step left state
            //Step 2.1: non-time variables
            {
                for (size_t stateIndex = 0; stateIndex < 10; ++stateIndex)
                {
                    for (size_t dIndex = 0; dIndex < this->Derivatives_of_StateStepLeftInertial.size(); ++dIndex) //loop over derivative entries for this state
                    {
                        size_t previousStateIndex = std::get<1>(this->Derivatives_of_StateStepLeftInertial[dIndex]);

                        if (stateIndex == previousStateIndex)
                        {
                            size_t Xindex = std::get<0>(this->Derivatives_of_StateStepLeftInertial[dIndex]);

                            //if this variable is not yet in the list of variables that affect the previous step's right state, add it to the list
                            if (std::find(this->ListOfVariablesAffectingCurrentStepLeftState.begin(), this->ListOfVariablesAffectingCurrentStepLeftState.end(), Xindex)
                                == this->ListOfVariablesAffectingCurrentStepLeftState.end())
                            {
                                this->ListOfVariablesAffectingCurrentStepLeftState.push_back(Xindex);
                                if (stateIndex < 3)
                                    this->ListOfVariablesAffectingCurrentStepLeftPosition.push_back(Xindex);
                                else if (stateIndex < 6)
                                    this->ListOfVariablesAffectingCurrentStepLeftVelocity.push_back(Xindex);
                            }
                        }
                    }
                }

                this->DerivativesOfCurrentStepLeftStateByVariable.resize(this->ListOfVariablesAffectingCurrentStepLeftState.size());

                for (size_t listIndex = 0; listIndex < this->ListOfVariablesAffectingCurrentStepLeftState.size(); ++listIndex)
                {
                    for (size_t stateIndex = 0; stateIndex < 10; ++stateIndex)
                    {
                        for (size_t dIndex = 0; dIndex < this->Derivatives_of_StateStepLeftInertial.size(); ++dIndex) //loop over derivative entries for this state
                        {
                            size_t previousStateIndex = std::get<1>(this->Derivatives_of_StateStepLeftInertial[dIndex]);

                            if (stateIndex == previousStateIndex)
                            {
                                size_t Xindex = std::get<0>(this->Derivatives_of_StateStepLeftInertial[dIndex]);

                                if (Xindex == this->ListOfVariablesAffectingCurrentStepLeftState[listIndex])
                                {
                                    this->DerivativesOfCurrentStepLeftStateByVariable[listIndex].push_back({ stateIndex, dIndex });
                                }
                            }
                        }
                    }
                }
            }//end non-time

            //Step 2.2: time variables
            {
                for (size_t stateIndex = 0; stateIndex < 10; ++stateIndex)
                {
                    for (size_t dIndex = 0; dIndex < this->Derivatives_of_StateStepLeftInertial_wrt_Time.size(); ++dIndex) //loop over derivative entries for this state
                    {
                        size_t previousStateIndex = std::get<1>(this->Derivatives_of_StateStepLeftInertial_wrt_Time[dIndex]);

                        if (stateIndex == previousStateIndex)
                        {
                            size_t Xindex = std::get<0>(this->Derivatives_of_StateStepLeftInertial_wrt_Time[dIndex]);

                            //if this variable is not yet in the list of variables that affect the previous step's right state, add it to the list
                            if (std::find(this->ListOfTimeVariablesAffectingCurrentStepLeftState.begin(), this->ListOfTimeVariablesAffectingCurrentStepLeftState.end(), Xindex)
                                == this->ListOfTimeVariablesAffectingCurrentStepLeftState.end())
                                this->ListOfTimeVariablesAffectingCurrentStepLeftState.push_back(Xindex);
                        }
                    }
                }

                this->DerivativesOfCurrentStepLeftStateByTimeVariable.resize(this->ListOfTimeVariablesAffectingCurrentStepLeftState.size());

                for (size_t listIndex = 0; listIndex < this->ListOfTimeVariablesAffectingCurrentStepLeftState.size(); ++listIndex)
                {
                    for (size_t stateIndex = 0; stateIndex < 10; ++stateIndex)
                    {
                        for (size_t dIndex = 0; dIndex < this->Derivatives_of_StateStepLeftInertial_wrt_Time.size(); ++dIndex) //loop over derivative entries for this state
                        {
                            size_t previousStateIndex = std::get<1>(this->Derivatives_of_StateStepLeftInertial_wrt_Time[dIndex]);

                            if (stateIndex == previousStateIndex)
                            {
                                size_t Xindex = std::get<0>(this->Derivatives_of_StateStepLeftInertial_wrt_Time[dIndex]);

                                if (Xindex == this->ListOfTimeVariablesAffectingCurrentStepLeftState[listIndex])
                                {
                                    this->DerivativesOfCurrentStepLeftStateByTimeVariable[listIndex].push_back({ stateIndex, dIndex });
                                }
                            }
                        }
                    }
                }
            }//end time

            //Step 3: construct the match point constraints
            std::vector<std::string>& matchPointConstraintNames = this->myPhase->get_matchPointConstraintNames();
            for (size_t constraintIndex = 0; constraintIndex < this->numMatchConstraints; ++constraintIndex)
            {
                Flowerbounds->push_back(-math::SMALL);
                Fupperbounds->push_back(math::SMALL);
                Fdescriptions->push_back(prefix + "left match point " + this->myPhase->get_matchPointConstraintNames()[constraintIndex]);
                this->Findices_left_match_point_constraints.push_back(Fdescriptions->size() - 1);
            }//end construction of match point constraints

            //Step 4: derivatives with respect to previous step right boundary
            {
                //Step 4.1: non-time
                this->Gindices_LeftMatchPoint_wrt_PreviousStepRightState.resize(this->ListOfVariablesAffectingPreviousStepRightState.size());
                for (size_t varIndex = 0; varIndex < this->ListOfVariablesAffectingPreviousStepRightState.size(); ++varIndex)
                {
                    size_t Xindex = this->ListOfVariablesAffectingPreviousStepRightState[varIndex];

                    for (size_t entryIndex = 0; entryIndex < this->DerivativesOfPreviousStepRightStateByVariable[varIndex].size(); ++entryIndex)
                    {
                        size_t stateIndex = std::get<0>(this->DerivativesOfPreviousStepRightStateByVariable[varIndex][entryIndex]);

                        size_t constraintIndex = stateIndex < 7 ? stateIndex : stateIndex - 1;

                        this->create_sparsity_entry(this->Findices_left_match_point_constraints[constraintIndex],
                            Xindex,
                            this->Gindices_LeftMatchPoint_wrt_PreviousStepRightState[varIndex]);
                    }
                }

                //Step 4.2: time
                this->Gindices_LeftMatchPoint_wrt_PreviousStepRightStateTime.resize(this->ListOfTimeVariablesAffectingPreviousStepRightState.size());
                for (size_t varIndex = 0; varIndex < this->ListOfTimeVariablesAffectingPreviousStepRightState.size(); ++varIndex)
                {
                    size_t Xindex = this->ListOfTimeVariablesAffectingPreviousStepRightState[varIndex];

                    for (size_t entryIndex = 0; entryIndex < this->DerivativesOfPreviousStepRightStateByTimeVariable[varIndex].size(); ++entryIndex)
                    {
                        size_t stateIndex = std::get<0>(this->DerivativesOfPreviousStepRightStateByTimeVariable[varIndex][entryIndex]);

                        size_t constraintIndex = stateIndex < 7 ? stateIndex : stateIndex - 1;

                        this->create_sparsity_entry(this->Findices_left_match_point_constraints[constraintIndex],
                            Xindex,
                            this->Gindices_LeftMatchPoint_wrt_PreviousStepRightStateTime[varIndex]);
                    }
                }
            }//end derivatives with respect to previous step right boundary

            //Step 5: derivatives with respect to current step left state
            if (this->myOptions->ParallelShootingStateRepresentation == StateRepresentation::SphericalRADEC)
            {
                //r
                {
                    const size_t stateIndices[] = { 0, 1, 2 };
                    for (size_t entryIndex = 0; entryIndex < 3; ++entryIndex)
                    {
                        size_t stateIndex = stateIndices[entryIndex];
                        size_t dIndex = this->dIndex_StateStepLeftInertial_wrt_r[entryIndex];

                        size_t Xindex = std::get<0>(this->Derivatives_of_StateStepLeftInertial[dIndex]);

                        this->create_sparsity_entry(this->Findices_left_match_point_constraints[stateIndex],
                            Xindex,
                            this->Gindices_StepLeftMatchPoint_wrt_r);
                    }
                }
                //RA
                {
                    const size_t stateIndices[] = { 0, 1 };
                    for (size_t entryIndex = 0; entryIndex < 2; ++entryIndex)
                    {
                        size_t stateIndex = stateIndices[entryIndex];
                        size_t dIndex = this->dIndex_StateStepLeftInertial_wrt_RA[entryIndex];

                        size_t Xindex = std::get<0>(this->Derivatives_of_StateStepLeftInertial[dIndex]);

                        this->create_sparsity_entry(this->Findices_left_match_point_constraints[stateIndex],
                            Xindex,
                            this->Gindices_StepLeftMatchPoint_wrt_RA);
                    }
                }
                //DEC
                {
                    const size_t stateIndices[] = { 0, 1, 2 };
                    for (size_t entryIndex = 0; entryIndex < 3; ++entryIndex)
                    {
                        size_t stateIndex = stateIndices[entryIndex];
                        size_t dIndex = this->dIndex_StateStepLeftInertial_wrt_DEC[entryIndex];

                        size_t Xindex = std::get<0>(this->Derivatives_of_StateStepLeftInertial[dIndex]);

                        this->create_sparsity_entry(this->Findices_left_match_point_constraints[stateIndex],
                            Xindex,
                            this->Gindices_StepLeftMatchPoint_wrt_DEC);
                    }
                }
                //v
                {
                    const size_t stateIndices[] = { 3, 4, 5 };
                    for (size_t entryIndex = 0; entryIndex < 3; ++entryIndex)
                    {
                        size_t stateIndex = stateIndices[entryIndex];
                        size_t dIndex = this->dIndex_StateStepLeftInertial_wrt_v[entryIndex];

                        size_t Xindex = std::get<0>(this->Derivatives_of_StateStepLeftInertial[dIndex]);

                        this->create_sparsity_entry(this->Findices_left_match_point_constraints[stateIndex],
                            Xindex,
                            this->Gindices_StepLeftMatchPoint_wrt_v);
                    }
                }
                //vRA
                {
                    const size_t stateIndices[] = { 3, 4 };
                    for (size_t entryIndex = 0; entryIndex < 2; ++entryIndex)
                    {
                        size_t stateIndex = stateIndices[entryIndex];
                        size_t dIndex = this->dIndex_StateStepLeftInertial_wrt_vRA[entryIndex];

                        size_t Xindex = std::get<0>(this->Derivatives_of_StateStepLeftInertial[dIndex]);

                        this->create_sparsity_entry(this->Findices_left_match_point_constraints[stateIndex],
                            Xindex,
                            this->Gindices_StepLeftMatchPoint_wrt_vRA);
                    }
                }
                //vDEC
                {
                    const size_t stateIndices[] = { 3, 4, 5 };
                    for (size_t entryIndex = 0; entryIndex < 3; ++entryIndex)
                    {
                        size_t stateIndex = stateIndices[entryIndex];
                        size_t dIndex = this->dIndex_StateStepLeftInertial_wrt_vDEC[entryIndex];

                        size_t Xindex = std::get<0>(this->Derivatives_of_StateStepLeftInertial[dIndex]);

                        this->create_sparsity_entry(this->Findices_left_match_point_constraints[stateIndex],
                            Xindex,
                            this->Gindices_StepLeftMatchPoint_wrt_vDEC);
                    }
                }
            }//end SphericalRADEC
            else if (this->myOptions->ParallelShootingStateRepresentation == StateRepresentation::SphericalAZFPA)
            {
                //r
                {
                    const size_t stateIndices[] = { 0, 1, 2 };
                    for (size_t entryIndex = 0; entryIndex < 3; ++entryIndex)
                    {
                        size_t stateIndex = stateIndices[entryIndex];
                        size_t dIndex = this->dIndex_StateStepLeftInertial_wrt_r[entryIndex];

                        size_t Xindex = std::get<0>(this->Derivatives_of_StateStepLeftInertial[dIndex]);

                        this->create_sparsity_entry(this->Findices_left_match_point_constraints[stateIndex],
                            Xindex,
                            this->Gindices_StepLeftMatchPoint_wrt_r);
                    }
                }
                //RA
                {
                    const size_t stateIndices[] = { 0, 1, 3, 4 };
                    for (size_t entryIndex = 0; entryIndex < 4; ++entryIndex)
                    {
                        size_t stateIndex = stateIndices[entryIndex];
                        size_t dIndex = this->dIndex_StateStepLeftInertial_wrt_RA[entryIndex];

                        size_t Xindex = std::get<0>(this->Derivatives_of_StateStepLeftInertial[dIndex]);

                        this->create_sparsity_entry(this->Findices_left_match_point_constraints[stateIndex],
                            Xindex,
                            this->Gindices_StepLeftMatchPoint_wrt_RA);
                    }
                }
                //DEC
                {
                    const size_t stateIndices[] = { 0, 1, 2, 3, 4, 5 };
                    for (size_t entryIndex = 0; entryIndex < 6; ++entryIndex)
                    {
                        size_t stateIndex = stateIndices[entryIndex];
                        size_t dIndex = this->dIndex_StateStepLeftInertial_wrt_DEC[entryIndex];

                        size_t Xindex = std::get<0>(this->Derivatives_of_StateStepLeftInertial[dIndex]);

                        this->create_sparsity_entry(this->Findices_left_match_point_constraints[stateIndex],
                            Xindex,
                            this->Gindices_StepLeftMatchPoint_wrt_DEC);
                    }
                }
                //v
                {
                    const size_t stateIndices[] = { 3, 4, 5 };
                    for (size_t entryIndex = 0; entryIndex < 3; ++entryIndex)
                    {
                        size_t stateIndex = stateIndices[entryIndex];
                        size_t dIndex = this->dIndex_StateStepLeftInertial_wrt_v[entryIndex];

                        size_t Xindex = std::get<0>(this->Derivatives_of_StateStepLeftInertial[dIndex]);

                        this->create_sparsity_entry(this->Findices_left_match_point_constraints[stateIndex],
                            Xindex,
                            this->Gindices_StepLeftMatchPoint_wrt_v);
                    }
                }
                //AZ
                {
                    const size_t stateIndices[] = { 3, 4, 5 };
                    for (size_t entryIndex = 0; entryIndex < 3; ++entryIndex)
                    {
                        size_t stateIndex = stateIndices[entryIndex];
                        size_t dIndex = this->dIndex_StateStepLeftInertial_wrt_AZ[entryIndex];

                        size_t Xindex = std::get<0>(this->Derivatives_of_StateStepLeftInertial[dIndex]);

                        this->create_sparsity_entry(this->Findices_left_match_point_constraints[stateIndex],
                            Xindex,
                            this->Gindices_StepLeftMatchPoint_wrt_AZ);
                    }
                }
                //vDEC
                {
                    const size_t stateIndices[] = { 3, 4, 5 };
                    for (size_t entryIndex = 0; entryIndex < 3; ++entryIndex)
                    {
                        size_t stateIndex = stateIndices[entryIndex];
                        size_t dIndex = this->dIndex_StateStepLeftInertial_wrt_FPA[entryIndex];

                        size_t Xindex = std::get<0>(this->Derivatives_of_StateStepLeftInertial[dIndex]);

                        this->create_sparsity_entry(this->Findices_left_match_point_constraints[stateIndex],
                            Xindex,
                            this->Gindices_StepLeftMatchPoint_wrt_FPA);
                    }
                }
            }//end SphericalAZFPA
            //mass
            this->create_sparsity_entry(this->Findices_left_match_point_constraints[6],
                this->Xindex_mass,
                this->Gindex_StepLeftMatchPoint_wrt_Mass);
            //virtual chemical fuel
            this->create_sparsity_entry(this->Findices_left_match_point_constraints[7],
                this->Xindex_chemical_fuel,
                this->Gindex_StepLeftMatchPoint_wrt_ChemicalFuel);
            //virtual electric propellant
            this->create_sparsity_entry(this->Findices_left_match_point_constraints[8],
                this->Xindex_electric_propellant,
                this->Gindex_StepLeftMatchPoint_wrt_ElectricPropellant);
        }//end calcbounds_step_left_match_point_constraints()

        void ParallelShootingStep::calcbounds_distance_constraints()
        {
            for (size_t constraintIndex = 0; constraintIndex < this->myDistanceConstraints.size(); ++constraintIndex)
                this->myDistanceConstraints[constraintIndex].calcbounds();
        }//end calcbounds_distance_constraints()

        void ParallelShootingStep::calcbounds_maneuver_constraints()
        {
            for (size_t constraintIndex = 0; constraintIndex < this->myManeuverConstraints.size(); ++constraintIndex)
                this->myManeuverConstraints[constraintIndex].calcbounds();
        }//end calcbounds_maneuver_constraints()

        void ParallelShootingStep::calcbounds_waypoint_tracking()
        {
            //Step 1: create the virtual waypoint error variable
            this->Xlowerbounds->push_back(0.0);
            this->Xupperbounds->push_back(this->myUniverse->LU);
            this->X_scale_factors->push_back(1.0);
            this->Xdescriptions->push_back(this->prefix + "virtual mahalanobis distance");

            //Step 2: create the virtual waypoint error constraint
            this->Flowerbounds->push_back(-math::SMALL);
            this->Fupperbounds->push_back(math::SMALL);
            this->Fdescriptions->push_back(this->prefix + "virtual mahalanobis distance");
            this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                this->Xdescriptions->size() - 1,
                this->Gindex_virtual_mahalanobis_distance_wrt_virtual_mahalanobis_distance);
            
            //Step 3: derivatives with respect to step state position, velocity, and mass
            this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                this->Xindex_rMag,
                this->Gindex_virtual_mahalanobis_distance_wrt_rMag);

            this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                this->Xindex_RA,
                this->Gindex_virtual_mahalanobis_distance_wrt_RA);

            this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                this->Xindex_DEC,
                this->Gindex_virtual_mahalanobis_distance_wrt_DEC);

            this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                this->Xindex_vMag,
                this->Gindex_virtual_mahalanobis_distance_wrt_vMag);

            if (this->myOptions->ParallelShootingStateRepresentation == StateRepresentation::SphericalRADEC)
            {
                this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                    this->Xindex_vRA,
                    this->Gindex_virtual_mahalanobis_distance_wrt_vRA);

                this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                    this->Xindex_vDEC,
                    this->Gindex_virtual_mahalanobis_distance_wrt_vDEC);
            }
            else if (this->myOptions->ParallelShootingStateRepresentation == StateRepresentation::SphericalAZFPA)
            {
                this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                    this->Xindex_AZ,
                    this->Gindex_virtual_mahalanobis_distance_wrt_AZ);

                this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                    this->Xindex_FPA,
                    this->Gindex_virtual_mahalanobis_distance_wrt_FPA);
            }

            this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                this->Xindex_mass,
                this->Gindex_virtual_mahalanobis_distance_wrt_mass);

            //Step 4: derivatives with respect to time
            for (size_t Xindex : this->Xindices_EventLeftEpoch)
                this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                    Xindex,
                    this->Gindex_virtual_mahalanobis_distance_wrt_Time);

            //Step 5: resize the error matrix and covariance inverse
            this->WaypointError = math::Matrix<doubleType>(7, 1, 0.0);
            this->CovarianceInverse = math::Matrix<doubleType>(7, 7, 0.0);
            this->CovarianceInverseDouble = math::Matrix<double>(7, 7, 0.0);
        }//end calcbounds_waypoint_tracking()

        void ParallelShootingStep::process_step_left_state(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            if (this->myOptions->ParallelShootingStateRepresentation == StateRepresentation::SphericalRADEC)
            {
                //Step 1: extract the thingies
                doubleType r = X[Xindex++];
                doubleType RA = X[Xindex++];
                doubleType DEC = X[Xindex++];
                doubleType v = X[Xindex++];
                doubleType vRA = X[Xindex++];
                doubleType vDEC = X[Xindex++];
                this->StateStepLeftInertial(6) = X[Xindex++];//mass
                this->StateStepLeftInertial(8) = X[Xindex++];//chemical fuel
                this->StateStepLeftInertial(9) = X[Xindex++];//electric propellant

                //Step 2: convert to cartesian
                doubleType cosRA = cos(RA);
                doubleType sinRA = sin(RA);
                doubleType cosDEC = cos(DEC);
                doubleType sinDEC = sin(DEC);
                doubleType cosvRA = cos(vRA);
                doubleType sinvRA = sin(vRA);
                doubleType cosvDEC = cos(vDEC);
                doubleType sinvDEC = sin(vDEC);

                this->StateStepLeftInertial(0) = r * cosRA * cosDEC;
                this->StateStepLeftInertial(1) = r * sinRA * cosDEC;
                this->StateStepLeftInertial(2) = r * sinDEC;
                this->StateStepLeftInertial(3) = v * cosvRA * cosvDEC;
                this->StateStepLeftInertial(4) = v * sinvRA * cosvDEC;
                this->StateStepLeftInertial(5) = v * sinvDEC;

                //Step 3: epoch
                this->StateStepLeftInertial(7) = this->StepLeftEpoch;

                //Step 4: derivatives
                if (needG)
                {
                    double dx_dr = (cosRA * cosDEC)_GETVALUE;
                    double dy_dr = (sinRA * cosDEC)_GETVALUE;
                    double dz_dr = sinDEC _GETVALUE;

                    double dx_dRA = (-r * cosDEC * sinRA) _GETVALUE;
                    double dy_dRA = (r * cosDEC * cosRA) _GETVALUE;

                    double dx_dDEC = (-r * cosRA*sinDEC)_GETVALUE;
                    double dy_dDEC = (-r * sinDEC*sinRA)_GETVALUE;
                    double dz_dDEC = (r*cosDEC)_GETVALUE;

                    double dxdot_dv = (cosvRA * cosvDEC)_GETVALUE;
                    double dydot_dv = (sinvRA * cosvDEC)_GETVALUE;
                    double dzdot_dv = sinvDEC _GETVALUE;

                    double dxdot_dvRA = (-v * cosvDEC * sinvRA) _GETVALUE;
                    double dydot_dvRA = (v * cosvDEC * cosvRA) _GETVALUE;

                    double dxdot_dvDEC = (-v * cosvRA*sinvDEC)_GETVALUE;
                    double dydot_dvDEC = (-v * sinvDEC*sinvRA)_GETVALUE;
                    double dzdot_dvDEC = (v*cosvDEC)_GETVALUE;

                    std::get<2>(this->Derivatives_of_StateStepLeftInertial[this->dIndex_StateStepLeftInertial_wrt_r[0]]) = dx_dr;
                    std::get<2>(this->Derivatives_of_StateStepLeftInertial[this->dIndex_StateStepLeftInertial_wrt_r[1]]) = dy_dr;
                    std::get<2>(this->Derivatives_of_StateStepLeftInertial[this->dIndex_StateStepLeftInertial_wrt_r[2]]) = dz_dr;

                    std::get<2>(this->Derivatives_of_StateStepLeftInertial[this->dIndex_StateStepLeftInertial_wrt_RA[0]]) = dx_dRA;
                    std::get<2>(this->Derivatives_of_StateStepLeftInertial[this->dIndex_StateStepLeftInertial_wrt_RA[1]]) = dy_dRA;

                    std::get<2>(this->Derivatives_of_StateStepLeftInertial[this->dIndex_StateStepLeftInertial_wrt_DEC[0]]) = dx_dDEC;
                    std::get<2>(this->Derivatives_of_StateStepLeftInertial[this->dIndex_StateStepLeftInertial_wrt_DEC[1]]) = dy_dDEC;
                    std::get<2>(this->Derivatives_of_StateStepLeftInertial[this->dIndex_StateStepLeftInertial_wrt_DEC[2]]) = dz_dDEC;

                    std::get<2>(this->Derivatives_of_StateStepLeftInertial[this->dIndex_StateStepLeftInertial_wrt_v[0]]) = dxdot_dv;
                    std::get<2>(this->Derivatives_of_StateStepLeftInertial[this->dIndex_StateStepLeftInertial_wrt_v[1]]) = dydot_dv;
                    std::get<2>(this->Derivatives_of_StateStepLeftInertial[this->dIndex_StateStepLeftInertial_wrt_v[2]]) = dzdot_dv;

                    std::get<2>(this->Derivatives_of_StateStepLeftInertial[this->dIndex_StateStepLeftInertial_wrt_vRA[0]]) = dxdot_dvRA;
                    std::get<2>(this->Derivatives_of_StateStepLeftInertial[this->dIndex_StateStepLeftInertial_wrt_vRA[1]]) = dydot_dvRA;

                    std::get<2>(this->Derivatives_of_StateStepLeftInertial[this->dIndex_StateStepLeftInertial_wrt_vDEC[0]]) = dxdot_dvDEC;
                    std::get<2>(this->Derivatives_of_StateStepLeftInertial[this->dIndex_StateStepLeftInertial_wrt_vDEC[1]]) = dydot_dvDEC;
                    std::get<2>(this->Derivatives_of_StateStepLeftInertial[this->dIndex_StateStepLeftInertial_wrt_vDEC[2]]) = dzdot_dvDEC;
                }//end derivatives
            }//end SphericalRADEC
            else if (this->myOptions->ParallelShootingStateRepresentation == StateRepresentation::SphericalAZFPA)
            {
                //Step 1: extract the thingies
                doubleType r = X[Xindex++];
                doubleType RA = X[Xindex++];
                doubleType DEC = X[Xindex++];
                doubleType v = X[Xindex++];
                doubleType AZ = X[Xindex++];
                doubleType FPA = X[Xindex++];
                this->StateStepLeftInertial(6) = X[Xindex++];//mass
                this->StateStepLeftInertial(8) = X[Xindex++];//chemical fuel
                this->StateStepLeftInertial(9) = X[Xindex++];//electric propellant

                //Step 2: convert to cartesian
                doubleType cosRA = cos(RA);
                doubleType sinRA = sin(RA);
                doubleType cosDEC = cos(DEC);
                doubleType sinDEC = sin(DEC);
                doubleType cosAZ = cos(AZ);
                doubleType sinAZ = sin(AZ);
                doubleType cosFPA = cos(FPA);
                doubleType sinFPA = sin(FPA);

                this->StateStepLeftInertial(0) = r * cosRA * cosDEC;
                this->StateStepLeftInertial(1) = r * sinRA * cosDEC;
                this->StateStepLeftInertial(2) = r * sinDEC;
                this->StateStepLeftInertial(3) = -v * (sinFPA*(sinAZ*sinRA + cosAZ * cosRA*sinDEC) - cosFPA * cosDEC*cosRA);
                this->StateStepLeftInertial(4) = v * (sinFPA*(cosRA*sinAZ - cosAZ * sinDEC*sinRA) + cosFPA * cosDEC*sinRA);
                this->StateStepLeftInertial(5) = v * (cosFPA*sinDEC + cosDEC * cosAZ*sinFPA);

                //Step 3: epoch
                this->StateStepLeftInertial(7) = this->StepLeftEpoch;

                //Step 4: derivatives
                if (needG)
                {
                    double dx_dr = (cosRA * cosDEC)_GETVALUE;
                    double dy_dr = (sinRA * cosDEC)_GETVALUE;
                    double dz_dr = sinDEC _GETVALUE;

                    double dx_dRA = (-r * cosDEC * sinRA) _GETVALUE;
                    double dy_dRA = (r * cosDEC * cosRA) _GETVALUE;

                    double dx_dDEC = (-r * cosRA*sinDEC)_GETVALUE;
                    double dy_dDEC = (-r * sinDEC*sinRA)_GETVALUE;
                    double dz_dDEC = (r*cosDEC)_GETVALUE;

                    double dxdot_dRA = (-v * (sinFPA*(cosRA*sinAZ - cosAZ * sinDEC*sinRA) + cosFPA * cosDEC*sinRA)) _GETVALUE;
                    double dydot_dRA = (-v * (sinFPA*(sinAZ*sinRA + cosAZ * cosRA*sinDEC) - cosFPA * cosDEC*cosRA)) _GETVALUE;

                    double dxdot_dDEC = (-v * cosRA*(cosFPA*sinDEC + cosDEC * cosAZ*sinFPA)) _GETVALUE;
                    double dydot_dDEC = (-v * sinRA*(cosFPA*sinDEC + cosDEC * cosAZ*sinFPA)) _GETVALUE;
                    double dzdot_dDEC = (v*(cosFPA*cosDEC - cosAZ * sinFPA*sinDEC)) _GETVALUE;

                    double dxdot_dv = (-(sinFPA*(sinAZ*sinRA + cosAZ * cosRA*sinDEC) - cosFPA * cosDEC*cosRA))_GETVALUE;
                    double dydot_dv = ((sinFPA*(cosRA*sinAZ - cosAZ * sinDEC*sinRA) + cosFPA * cosDEC*sinRA))_GETVALUE;
                    double dzdot_dv = ((cosFPA*sinDEC + cosDEC * cosAZ*sinFPA))_GETVALUE;

                    double dxdot_dFPA = (-v * (cosFPA*(sinAZ*sinRA + cosAZ * cosRA*sinDEC) + cosDEC * sinFPA*cosRA)) _GETVALUE;
                    double dydot_dFPA = (v*(cosFPA*(cosRA*sinAZ - cosAZ * sinDEC*sinRA) - cosDEC * sinFPA*sinRA)) _GETVALUE;
                    double dzdot_dFPA = (v*cosFPA*cosDEC*cosAZ - v * sinFPA*sinDEC) _GETVALUE;

                    double dxdot_dAZ = (-v * sinFPA*(cosAZ*sinRA - cosRA * sinDEC*sinAZ)) _GETVALUE;
                    double dydot_dAZ = (v*sinFPA*(cosAZ*cosRA + sinDEC * sinAZ*sinRA)) _GETVALUE;
                    double dzdot_dAZ = (-v * cosDEC*sinFPA*sinAZ) _GETVALUE;

                    std::get<2>(this->Derivatives_of_StateStepLeftInertial[this->dIndex_StateStepLeftInertial_wrt_r[0]]) = dx_dr;
                    std::get<2>(this->Derivatives_of_StateStepLeftInertial[this->dIndex_StateStepLeftInertial_wrt_r[1]]) = dy_dr;
                    std::get<2>(this->Derivatives_of_StateStepLeftInertial[this->dIndex_StateStepLeftInertial_wrt_r[2]]) = dz_dr;

                    std::get<2>(this->Derivatives_of_StateStepLeftInertial[this->dIndex_StateStepLeftInertial_wrt_RA[0]]) = dx_dRA;
                    std::get<2>(this->Derivatives_of_StateStepLeftInertial[this->dIndex_StateStepLeftInertial_wrt_RA[1]]) = dy_dRA;
                    std::get<2>(this->Derivatives_of_StateStepLeftInertial[this->dIndex_StateStepLeftInertial_wrt_RA[2]]) = dxdot_dRA;
                    std::get<2>(this->Derivatives_of_StateStepLeftInertial[this->dIndex_StateStepLeftInertial_wrt_RA[3]]) = dydot_dRA;

                    std::get<2>(this->Derivatives_of_StateStepLeftInertial[this->dIndex_StateStepLeftInertial_wrt_DEC[0]]) = dx_dDEC;
                    std::get<2>(this->Derivatives_of_StateStepLeftInertial[this->dIndex_StateStepLeftInertial_wrt_DEC[1]]) = dy_dDEC;
                    std::get<2>(this->Derivatives_of_StateStepLeftInertial[this->dIndex_StateStepLeftInertial_wrt_DEC[2]]) = dz_dDEC;
                    std::get<2>(this->Derivatives_of_StateStepLeftInertial[this->dIndex_StateStepLeftInertial_wrt_DEC[3]]) = dxdot_dDEC;
                    std::get<2>(this->Derivatives_of_StateStepLeftInertial[this->dIndex_StateStepLeftInertial_wrt_DEC[4]]) = dydot_dDEC;
                    std::get<2>(this->Derivatives_of_StateStepLeftInertial[this->dIndex_StateStepLeftInertial_wrt_DEC[5]]) = dzdot_dDEC;

                    std::get<2>(this->Derivatives_of_StateStepLeftInertial[this->dIndex_StateStepLeftInertial_wrt_v[0]]) = dxdot_dv;
                    std::get<2>(this->Derivatives_of_StateStepLeftInertial[this->dIndex_StateStepLeftInertial_wrt_v[1]]) = dydot_dv;
                    std::get<2>(this->Derivatives_of_StateStepLeftInertial[this->dIndex_StateStepLeftInertial_wrt_v[2]]) = dzdot_dv;

                    std::get<2>(this->Derivatives_of_StateStepLeftInertial[this->dIndex_StateStepLeftInertial_wrt_AZ[0]]) = dxdot_dAZ;
                    std::get<2>(this->Derivatives_of_StateStepLeftInertial[this->dIndex_StateStepLeftInertial_wrt_AZ[1]]) = dydot_dAZ;
                    std::get<2>(this->Derivatives_of_StateStepLeftInertial[this->dIndex_StateStepLeftInertial_wrt_AZ[2]]) = dzdot_dAZ;

                    std::get<2>(this->Derivatives_of_StateStepLeftInertial[this->dIndex_StateStepLeftInertial_wrt_FPA[0]]) = dxdot_dFPA;
                    std::get<2>(this->Derivatives_of_StateStepLeftInertial[this->dIndex_StateStepLeftInertial_wrt_FPA[1]]) = dydot_dFPA;
                    std::get<2>(this->Derivatives_of_StateStepLeftInertial[this->dIndex_StateStepLeftInertial_wrt_FPA[2]]) = dzdot_dFPA;
                }//end derivatives
            }//end SphericalAZFPA

             //Step 5: epoch
            this->process_epoch_time(X, Xindex, F, Findex, G, needG);
        }//end process_step_left_state

        void ParallelShootingStep::process_epoch_time(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            this->PhaseFlightTime = X[this->Xindex_PhaseFlightTime];
            this->StepFlightTime = (this->PhaseFlightTime - this->myPhase->getInitialCoastDuration() - this->myPhase->getTerminalCoastDuration()) * this->dStepTime_dPhaseFlightTime;
            this->LaunchDate = X[this->Xindices_EventLeftEpoch[0]];

            if (this->stepIndex == 0)
                this->StepLeftEpoch = this->myPhase->get_StateAfterInitialCoast()(7);
            else
                this->StepLeftEpoch = this->previousStep->get_StateStepRightInertial()(7);

            this->StateStepLeftInertial(7) = this->StepLeftEpoch;
        }//end process_epoch_time
        
        void ParallelShootingStep::process_step_control(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            for (size_t subStepIndex = 0; subStepIndex < this->num_interior_control_points; ++subStepIndex)
            {
                //extract all of the control variables - yay!
                for (size_t controlIndex = 0; controlIndex < this->num_controls; ++controlIndex)
                {
                    this->ControlVector[subStepIndex](controlIndex) = X[Xindex++];
                }//end loop over controls

                    //control magnitude constraint
                this->throttle[subStepIndex] = sqrt(this->ControlVector[subStepIndex](0) * this->ControlVector[subStepIndex](0)
                                                  + this->ControlVector[subStepIndex](1) * this->ControlVector[subStepIndex](1)
                                                  + this->ControlVector[subStepIndex](2) * this->ControlVector[subStepIndex](2)
                                                  + 1.0e-25);
                F[Findex++] = this->throttle[subStepIndex];

                //derivatives
                if (needG)
                {
                    for (size_t controlIndex = 0; controlIndex < this->num_controls; ++controlIndex)
                    {
                        size_t Gindex = this->G_indices_control_magnitude_constraint_wrt_Control[subStepIndex][controlIndex];
                        size_t Xindex = this->jGvar->operator[](Gindex);

                        G[Gindex] = this->X_scale_factors->operator[](Xindex)
                                  * (this->ControlVector[subStepIndex](controlIndex) / this->throttle[subStepIndex]) _GETVALUE;
                    }//end loop over controls
                }//end derivatives
            }
        }//end process_step_control

        void ParallelShootingStep::process_step_left_match_point_constraints(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //Step 1: compute the match point constraints
            std::vector<size_t>& matchPointConstraintStateIndex = this->myPhase->get_matchPointConstraintStateIndex();
            math::Matrix<double>& continuity_constraint_scale_factors = this->myPhase->get_continuity_constraint_scale_factors();
            math::Matrix<doubleType>& PreviousStepRightInertial = this->previousStep->get_StateStepRightInertial();
            for (size_t constraintIndex = 0; constraintIndex < this->numMatchConstraints; ++constraintIndex)
            {
                size_t stateIndex = matchPointConstraintStateIndex[constraintIndex];
                F[Findex++] = (this->StateStepLeftInertial(stateIndex) - PreviousStepRightInertial(stateIndex))
                    * continuity_constraint_scale_factors(constraintIndex);
            }

            //Step 2: derivatives of the match point constraints
            if (needG)
            {
                math::Matrix<double> dMatchState_dDecisionVariable(this->StateStepLeftInertial.get_n(), 1, 0.0);

                if (this->myOptions->ParallelShootingStateRepresentation == StateRepresentation::SphericalRADEC)
                {
                    //Step 2.1: with respect to the current step
                    //r
                    {
                        const size_t stateIndices[] = { 0, 1, 2 };
                        for (size_t entryIndex = 0; entryIndex < 3; ++entryIndex)
                        {
                            size_t stateIndex = stateIndices[entryIndex];
                            size_t dIndex = this->dIndex_StateStepLeftInertial_wrt_r[entryIndex];

                            size_t Xindex = std::get<0>(this->Derivatives_of_StateStepLeftInertial[dIndex]);
                            size_t Gindex = this->Gindices_StepLeftMatchPoint_wrt_r[entryIndex];

                            G[Gindex] = this->X_scale_factors->operator[](Xindex)
                                * std::get<2>(this->Derivatives_of_StateStepLeftInertial[dIndex])
                                * continuity_constraint_scale_factors(stateIndex);
                        }
                    }
                    //RA
                    {
                        const size_t stateIndices[] = { 0, 1 };
                        for (size_t entryIndex = 0; entryIndex < 2; ++entryIndex)
                        {
                            size_t stateIndex = stateIndices[entryIndex];
                            size_t dIndex = this->dIndex_StateStepLeftInertial_wrt_RA[entryIndex];

                            size_t Xindex = std::get<0>(this->Derivatives_of_StateStepLeftInertial[dIndex]);
                            size_t Gindex = this->Gindices_StepLeftMatchPoint_wrt_RA[entryIndex];

                            G[Gindex] = this->X_scale_factors->operator[](Xindex)
                                * std::get<2>(this->Derivatives_of_StateStepLeftInertial[dIndex])
                                * continuity_constraint_scale_factors(stateIndex);
                        }
                    }
                    //DEC
                    {
                        const size_t stateIndices[] = { 0, 1, 2 };
                        for (size_t entryIndex = 0; entryIndex < 3; ++entryIndex)
                        {
                            size_t stateIndex = stateIndices[entryIndex];
                            size_t dIndex = this->dIndex_StateStepLeftInertial_wrt_DEC[entryIndex];

                            size_t Xindex = std::get<0>(this->Derivatives_of_StateStepLeftInertial[dIndex]);
                            size_t Gindex = this->Gindices_StepLeftMatchPoint_wrt_DEC[entryIndex];

                            G[Gindex] = this->X_scale_factors->operator[](Xindex)
                                * std::get<2>(this->Derivatives_of_StateStepLeftInertial[dIndex])
                                * continuity_constraint_scale_factors(stateIndex);
                        }
                    }
                    //v
                    {
                        const size_t stateIndices[] = { 3, 4, 5 };
                        for (size_t entryIndex = 0; entryIndex < 3; ++entryIndex)
                        {
                            size_t stateIndex = stateIndices[entryIndex];
                            size_t dIndex = this->dIndex_StateStepLeftInertial_wrt_v[entryIndex];

                            size_t Xindex = std::get<0>(this->Derivatives_of_StateStepLeftInertial[dIndex]);
                            size_t Gindex = this->Gindices_StepLeftMatchPoint_wrt_v[entryIndex];

                            G[Gindex] = this->X_scale_factors->operator[](Xindex)
                                * std::get<2>(this->Derivatives_of_StateStepLeftInertial[dIndex])
                                * continuity_constraint_scale_factors(stateIndex);
                        }
                    }
                    //vRA
                    {
                        const size_t stateIndices[] = { 3, 4 };
                        for (size_t entryIndex = 0; entryIndex < 2; ++entryIndex)
                        {
                            size_t stateIndex = stateIndices[entryIndex];
                            size_t dIndex = this->dIndex_StateStepLeftInertial_wrt_vRA[entryIndex];

                            size_t Xindex = std::get<0>(this->Derivatives_of_StateStepLeftInertial[dIndex]);
                            size_t Gindex = this->Gindices_StepLeftMatchPoint_wrt_vRA[entryIndex];

                            G[Gindex] = this->X_scale_factors->operator[](Xindex)
                                * std::get<2>(this->Derivatives_of_StateStepLeftInertial[dIndex])
                                * continuity_constraint_scale_factors(stateIndex);
                        }
                    }
                    //vDEC
                    {
                        const size_t stateIndices[] = { 3, 4, 5 };
                        for (size_t entryIndex = 0; entryIndex < 3; ++entryIndex)
                        {
                            size_t stateIndex = stateIndices[entryIndex];
                            size_t dIndex = this->dIndex_StateStepLeftInertial_wrt_vDEC[entryIndex];

                            size_t Xindex = std::get<0>(this->Derivatives_of_StateStepLeftInertial[dIndex]);
                            size_t Gindex = this->Gindices_StepLeftMatchPoint_wrt_vDEC[entryIndex];

                            G[Gindex] = this->X_scale_factors->operator[](Xindex)
                                * std::get<2>(this->Derivatives_of_StateStepLeftInertial[dIndex])
                                * continuity_constraint_scale_factors(stateIndex);
                        }
                    }
                }//end SphericalRADEC
                else if (this->myOptions->ParallelShootingStateRepresentation == StateRepresentation::SphericalAZFPA)
                {
                    //Step 2.1: with respect to the current step
                    //r
                    {
                        const size_t stateIndices[] = { 0, 1, 2 };
                        for (size_t entryIndex = 0; entryIndex < 3; ++entryIndex)
                        {
                            size_t stateIndex = stateIndices[entryIndex];
                            size_t dIndex = this->dIndex_StateStepLeftInertial_wrt_r[entryIndex];

                            size_t Xindex = std::get<0>(this->Derivatives_of_StateStepLeftInertial[dIndex]);
                            size_t Gindex = this->Gindices_StepLeftMatchPoint_wrt_r[entryIndex];

                            G[Gindex] = this->X_scale_factors->operator[](Xindex)
                                * std::get<2>(this->Derivatives_of_StateStepLeftInertial[dIndex])
                                * continuity_constraint_scale_factors(stateIndex);
                        }
                    }
                    //RA
                    {
                        const size_t stateIndices[] = { 0, 1, 3, 4 };
                        for (size_t entryIndex = 0; entryIndex < 4; ++entryIndex)
                        {
                            size_t stateIndex = stateIndices[entryIndex];
                            size_t dIndex = this->dIndex_StateStepLeftInertial_wrt_RA[entryIndex];

                            size_t Xindex = std::get<0>(this->Derivatives_of_StateStepLeftInertial[dIndex]);
                            size_t Gindex = this->Gindices_StepLeftMatchPoint_wrt_RA[entryIndex];

                            G[Gindex] = this->X_scale_factors->operator[](Xindex)
                                * std::get<2>(this->Derivatives_of_StateStepLeftInertial[dIndex])
                                * continuity_constraint_scale_factors(stateIndex);
                        }
                    }
                    //DEC
                    {
                        const size_t stateIndices[] = { 0, 1, 2, 3, 4, 5 };
                        for (size_t entryIndex = 0; entryIndex < 6; ++entryIndex)
                        {
                            size_t stateIndex = stateIndices[entryIndex];
                            size_t dIndex = this->dIndex_StateStepLeftInertial_wrt_DEC[entryIndex];

                            size_t Xindex = std::get<0>(this->Derivatives_of_StateStepLeftInertial[dIndex]);
                            size_t Gindex = this->Gindices_StepLeftMatchPoint_wrt_DEC[entryIndex];

                            G[Gindex] = this->X_scale_factors->operator[](Xindex)
                                * std::get<2>(this->Derivatives_of_StateStepLeftInertial[dIndex])
                                * continuity_constraint_scale_factors(stateIndex);
                        }
                    }
                    //v
                    {
                        const size_t stateIndices[] = { 3, 4, 5 };
                        for (size_t entryIndex = 0; entryIndex < 3; ++entryIndex)
                        {
                            size_t stateIndex = stateIndices[entryIndex];
                            size_t dIndex = this->dIndex_StateStepLeftInertial_wrt_v[entryIndex];

                            size_t Xindex = std::get<0>(this->Derivatives_of_StateStepLeftInertial[dIndex]);
                            size_t Gindex = this->Gindices_StepLeftMatchPoint_wrt_v[entryIndex];

                            G[Gindex] = this->X_scale_factors->operator[](Xindex)
                                * std::get<2>(this->Derivatives_of_StateStepLeftInertial[dIndex])
                                * continuity_constraint_scale_factors(stateIndex);
                        }
                    }
                    //AZ
                    {
                        const size_t stateIndices[] = { 3, 4, 5 };
                        for (size_t entryIndex = 0; entryIndex < 3; ++entryIndex)
                        {
                            size_t stateIndex = stateIndices[entryIndex];
                            size_t dIndex = this->dIndex_StateStepLeftInertial_wrt_AZ[entryIndex];

                            size_t Xindex = std::get<0>(this->Derivatives_of_StateStepLeftInertial[dIndex]);
                            size_t Gindex = this->Gindices_StepLeftMatchPoint_wrt_AZ[entryIndex];

                            G[Gindex] = this->X_scale_factors->operator[](Xindex)
                                * std::get<2>(this->Derivatives_of_StateStepLeftInertial[dIndex])
                                * continuity_constraint_scale_factors(stateIndex);
                        }
                    }
                    //FPA
                    {
                        const size_t stateIndices[] = { 3, 4, 5 };
                        for (size_t entryIndex = 0; entryIndex < 3; ++entryIndex)
                        {
                            size_t stateIndex = stateIndices[entryIndex];
                            size_t dIndex = this->dIndex_StateStepLeftInertial_wrt_FPA[entryIndex];

                            size_t Xindex = std::get<0>(this->Derivatives_of_StateStepLeftInertial[dIndex]);
                            size_t Gindex = this->Gindices_StepLeftMatchPoint_wrt_FPA[entryIndex];

                            G[Gindex] = this->X_scale_factors->operator[](Xindex)
                                * std::get<2>(this->Derivatives_of_StateStepLeftInertial[dIndex])
                                * continuity_constraint_scale_factors(stateIndex);
                        }
                    }
                }//end SphericalAZFPA
                //mass
                {
                    size_t Gindex = this->Gindex_StepLeftMatchPoint_wrt_Mass;
                    size_t Xindex = this->Xindex_mass;

                    G[Gindex] = this->X_scale_factors->operator[](Xindex)
                        * continuity_constraint_scale_factors(6);
                }
                //virtual chemical fuel
                {
                    size_t Gindex = this->Gindex_StepLeftMatchPoint_wrt_ChemicalFuel;
                    size_t Xindex = this->Xindex_chemical_fuel;

                    G[Gindex] = this->X_scale_factors->operator[](Xindex)
                        * continuity_constraint_scale_factors(7);
                }
                //virtual electric propellant
                {
                    size_t Gindex = this->Gindex_StepLeftMatchPoint_wrt_ElectricPropellant;
                    size_t Xindex = this->Xindex_electric_propellant;

                    G[Gindex] = this->X_scale_factors->operator[](Xindex)
                        * continuity_constraint_scale_factors(8);
                }


                //Step 2.2: with respect to the previous step
                std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_PreviousStepRightState = this->previousStep->get_Derivatives_of_StateStepRightInertial();//Xindex, stateIndex, derivative value
                std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_PreviousStepRightState_wrt_Time = this->previousStep->get_Derivatives_of_StateStepRightInertial_wrt_Time();//Xindex, stateIndex, derivative value

                //Step 2.2.1: non-time
                for (size_t varIndex = 0; varIndex < this->ListOfVariablesAffectingPreviousStepRightState.size(); ++varIndex)
                {
                    size_t Xindex = this->ListOfVariablesAffectingPreviousStepRightState[varIndex];

                    for (size_t entryIndex = 0; entryIndex < this->DerivativesOfPreviousStepRightStateByVariable[varIndex].size(); ++entryIndex)
                    {
                        size_t stateIndex = std::get<0>(this->DerivativesOfPreviousStepRightStateByVariable[varIndex][entryIndex]);

                        if (stateIndex != 7) //there is no constraint on epoch
                        {
                            size_t dIndex = std::get<1>(this->DerivativesOfPreviousStepRightStateByVariable[varIndex][entryIndex]);

                            size_t constraintIndex = stateIndex < 7 ? stateIndex : stateIndex - 1;

                            double TheDerivative = std::get<2>(Derivatives_of_PreviousStepRightState[dIndex]);

                            size_t Gindex = this->Gindices_LeftMatchPoint_wrt_PreviousStepRightState[varIndex][entryIndex];
                            size_t Xindex = this->jGvar->operator[](Gindex);

                            G[Gindex] = this->X_scale_factors->operator[](Xindex)
                                * -TheDerivative
                                * continuity_constraint_scale_factors(constraintIndex);
                        }
                    }
                }

                //Step 2.2.1: time
                for (size_t varIndex = 0; varIndex < this->ListOfTimeVariablesAffectingPreviousStepRightState.size(); ++varIndex)
                {
                    size_t Xindex = this->ListOfTimeVariablesAffectingPreviousStepRightState[varIndex];

                    for (size_t entryIndex = 0; entryIndex < this->DerivativesOfPreviousStepRightStateByTimeVariable[varIndex].size(); ++entryIndex)
                    {
                        size_t stateIndex = std::get<0>(this->DerivativesOfPreviousStepRightStateByTimeVariable[varIndex][entryIndex]);

                        if (stateIndex != 7) //there is no constraint on epoch
                        {
                            size_t dIndex = std::get<1>(this->DerivativesOfPreviousStepRightStateByTimeVariable[varIndex][entryIndex]);

                            size_t constraintIndex = stateIndex < 7 ? stateIndex : stateIndex - 1;

                            double TheDerivative = std::get<2>(Derivatives_of_PreviousStepRightState_wrt_Time[dIndex]);

                            size_t Gindex = this->Gindices_LeftMatchPoint_wrt_PreviousStepRightStateTime[varIndex][entryIndex];
                            size_t Xindex = this->jGvar->operator[](Gindex);

                            G[Gindex] = this->X_scale_factors->operator[](Xindex)
                                * -TheDerivative
                                * continuity_constraint_scale_factors(constraintIndex);
                        }
                    }
                }
            }//end match point constraint derivatives
        }//end process_step_left_match_point_constraints

        void ParallelShootingStep::process_distance_constraints(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            for (size_t constraintIndex = 0; constraintIndex < this->myDistanceConstraints.size(); ++constraintIndex)
                this->myDistanceConstraints[constraintIndex].process(X, Xindex, F, Findex, G, needG);
        }//end process_distance_constraints
        
        void ParallelShootingStep::process_maneuver_constraints(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            for (size_t constraintIndex = 0; constraintIndex < this->myManeuverConstraints.size(); ++constraintIndex)
                this->myManeuverConstraints[constraintIndex].process_constraint(X, Xindex, F, Findex, G, needG);
        }//end process_maneuver_constraints



        void ParallelShootingStep::process_waypoint_tracking(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //Step 1: extract the virtual waypoint error
            this->virtual_mahalanobis_distance = X[Xindex++];

            //Step 2: compute the actual waypoint error
            double waypoint_state[7];
            try
            {
                //this call will fail if you run off the end of a spline
                
                //get the waypoint state
                //this->myOptions->myWaypointTracker.getPosition(this->StateStepLeftInertial(7)_GETVALUE, waypoint_state);

                math::Matrix<doubleType> waypoint_state_vector(7, 1, std::vector<doubleType>(waypoint_state, waypoint_state + 7));

                //get the waypoint covariance (or rather, the inverse of said covariance
                //this->myOptions->myCovarianceReader.getPinv(this->StateStepLeftInertial(7)_GETVALUE, this->CovarianceInverseDouble);

                //this step is necessary to make it possible to check algorithmic derivatives
                for (size_t i = 0; i < 7; ++i)
                {
                    for (size_t j = 0; j < 7; ++j)
                        this->CovarianceInverse(i, j) = this->CovarianceInverseDouble(i, j);
                }
                

#ifdef AD_INSTRUMENTATION
                //state derivatives
                double waypoint_state_derivatives[7];
                //this->myOptions->myWaypointTracker.get7StateDerivative(this->StateStepLeftInertial(7)_GETVALUE, waypoint_state_derivatives);

                for (size_t stateIndex = 0; stateIndex < 7; ++stateIndex)
                {
                    for (size_t tIndex = 0; tIndex < this->Xindices_EventLeftEpoch.size(); ++tIndex)
                    {
                        size_t tXindex = this->Xindices_EventLeftEpoch[tIndex];
                        waypoint_state_vector(stateIndex).setDerivative(tIndex, waypoint_state_derivatives[stateIndex] * this->X_scale_factors->operator[](tXindex));
                    }
                    waypoint_state_vector(stateIndex).setDerivative(this->Xindex_PhaseFlightTime, waypoint_state_derivatives[stateIndex] * this->dStepTime_dPhaseFlightTime * this->X_scale_factors->operator[](this->Xindex_PhaseFlightTime));
                }

                //covariance derivatives
                //this->myOptions->myCovarianceReader.getPinvDerivative(this->StateStepLeftInertial(7)_GETVALUE, this->dCovarianceInverse_depoch);
                
                for (size_t i = 0; i < 7; ++i)
                {
                    for (size_t j = 0; j < 7; ++j)
                    {
                        for (size_t tIndex = 0; tIndex < this->Xindices_EventLeftEpoch.size(); ++tIndex)
                        {
                            size_t tXindex = this->Xindices_EventLeftEpoch[tIndex];
                            this->CovarianceInverse(i, j).setDerivative(tIndex, this->dCovarianceInverse_depoch(i, j) * this->X_scale_factors->operator[](tXindex));
                        }
                        this->CovarianceInverse(i, j).setDerivative(this->Xindex_PhaseFlightTime, this->dCovarianceInverse_depoch(i, j) * this->dStepTime_dPhaseFlightTime * this->X_scale_factors->operator[](this->Xindex_PhaseFlightTime));
                    }
                }
#endif
                
                //find the actual distance from the centroid
                this->WaypointError = this->StateStepLeftInertial.getSubMatrix1D(0, 6) - waypoint_state_vector;
                                
                //Step 3: compute the weighted waypoint error, where the weights themselves are functions of time
                //note that the .norm() is a bit of a misnomer - math::Matrix thinks that a 1x1 matrix is a matrix and .norm() turns it into a scalar
                this->mahalanobis_distance = sqrt((this->WaypointError.transpose() * this->CovarianceInverse * this->WaypointError).norm())
                    / this->myUniverse->LU;
            }
            catch(std::exception &error)
            {
                //if you do run off the end of a spline, then there isn't a waypoint to track so you can do whatever you want and the error is zero anyway
                this->mahalanobis_distance = 0.0;
            }
            
            //Step 4: compute the waypoint error constraint
            F[Findex++] = this->virtual_mahalanobis_distance - this->mahalanobis_distance;

            //Step 5: derivatives
            if (needG)
            {
                try
                {
                    //Step 5.1: components
                    //state derivatives
                    math::Matrix<double> waypoint_state_derivatives(7, 1, 0.0);
                    double temp[7];
                    //this->myOptions->myWaypointTracker.get7StateDerivative(this->StateStepLeftInertial(7)_GETVALUE, temp);
                    waypoint_state_derivatives.assign_all(temp);
                    
                    //covariance derivatives
                    //this->myOptions->myCovarianceReader.getPinvDerivative(this->StateStepLeftInertial(7)_GETVALUE, this->dCovarianceInverse_depoch);

                    //Step 5.2: derivative with respect to virtual waypoint error variable
                    {
                        size_t Gindex = this->Gindex_virtual_mahalanobis_distance_wrt_virtual_mahalanobis_distance;
                        size_t Xindex = this->jGvar->operator[](Gindex);

                        G[Gindex] = this->X_scale_factors->operator[](Xindex);
                    }

                    //Step 5.3: derivatives with respect to state
                    double M = this->mahalanobis_distance _GETVALUE;
                    math::Matrix<doubleType> dM_dX_adouble = (this->WaypointError.transpose() * (this->CovarianceInverse.transpose() + this->CovarianceInverse)) / (2.0 * M);

                    math::Matrix<double> dM_dX(7, 1, 0.0);
                    for (size_t stateIndex = 0; stateIndex < 7; ++stateIndex)
                        dM_dX(stateIndex) = dM_dX_adouble(stateIndex) _GETVALUE;

                    //now we can chain these with the derivatives of state with respect to decision variables

                    if (this->myOptions->ParallelShootingStateRepresentation == StateRepresentation::SphericalRADEC)
                    {
                        //rMag
                        {
                            size_t Gindex = this->Gindex_virtual_mahalanobis_distance_wrt_rMag;
                            G[Gindex] = 0.0;

                            for (size_t dIndex : this->dIndex_StateStepLeftInertial_wrt_r)
                            {
                                size_t stateIndex = std::get<1>(this->Derivatives_of_StateStepLeftInertial[dIndex]);
                                size_t Xindex = std::get<0>(this->Derivatives_of_StateStepLeftInertial[dIndex]);

                                double Gentry = dM_dX(stateIndex) * std::get<2>(this->Derivatives_of_StateStepLeftInertial[dIndex]);

                                G[Gindex] += this->X_scale_factors->operator[](Xindex)
                                    * Gentry
                                    / this->myUniverse->LU;
                            }
                        }
                        //RA
                        {
                            size_t Gindex = this->Gindex_virtual_mahalanobis_distance_wrt_RA;
                            G[Gindex] = 0.0;

                            for (size_t dIndex : this->dIndex_StateStepLeftInertial_wrt_RA)
                            {
                                size_t stateIndex = std::get<1>(this->Derivatives_of_StateStepLeftInertial[dIndex]);
                                size_t Xindex = std::get<0>(this->Derivatives_of_StateStepLeftInertial[dIndex]);

                                double Gentry = dM_dX(stateIndex) * std::get<2>(this->Derivatives_of_StateStepLeftInertial[dIndex]);

                                G[Gindex] += this->X_scale_factors->operator[](Xindex)
                                    * Gentry
                                    / this->myUniverse->LU;
                            }
                        }
                        //DEC
                        {
                            size_t Gindex = this->Gindex_virtual_mahalanobis_distance_wrt_DEC;
                            G[Gindex] = 0.0;

                            for (size_t dIndex : this->dIndex_StateStepLeftInertial_wrt_DEC)
                            {
                                size_t stateIndex = std::get<1>(this->Derivatives_of_StateStepLeftInertial[dIndex]);
                                size_t Xindex = std::get<0>(this->Derivatives_of_StateStepLeftInertial[dIndex]);

                                double Gentry = dM_dX(stateIndex) * std::get<2>(this->Derivatives_of_StateStepLeftInertial[dIndex]);

                                G[Gindex] += this->X_scale_factors->operator[](Xindex)
                                    * Gentry
                                    / this->myUniverse->LU;
                            }
                        }
                        //vMag
                        {
                            size_t Gindex = this->Gindex_virtual_mahalanobis_distance_wrt_vMag;
                            G[Gindex] = 0.0;

                            for (size_t dIndex : this->dIndex_StateStepLeftInertial_wrt_v)
                            {
                                size_t stateIndex = std::get<1>(this->Derivatives_of_StateStepLeftInertial[dIndex]);
                                size_t Xindex = std::get<0>(this->Derivatives_of_StateStepLeftInertial[dIndex]);

                                double Gentry = dM_dX(stateIndex) * std::get<2>(this->Derivatives_of_StateStepLeftInertial[dIndex]);

                                G[Gindex] += this->X_scale_factors->operator[](Xindex)
                                    * Gentry
                                    / this->myUniverse->LU;
                            }
                        }
                        //AZ
                        {
                            size_t Gindex = this->Gindex_virtual_mahalanobis_distance_wrt_vRA;
                            G[Gindex] = 0.0;

                            for (size_t dIndex : this->dIndex_StateStepLeftInertial_wrt_vRA)
                            {
                                size_t stateIndex = std::get<1>(this->Derivatives_of_StateStepLeftInertial[dIndex]);
                                size_t Xindex = std::get<0>(this->Derivatives_of_StateStepLeftInertial[dIndex]);

                                double Gentry = dM_dX(stateIndex) * std::get<2>(this->Derivatives_of_StateStepLeftInertial[dIndex]);

                                G[Gindex] += this->X_scale_factors->operator[](Xindex)
                                    * Gentry
                                    / this->myUniverse->LU;
                            }
                        }
                        //FPA
                        {
                            size_t Gindex = this->Gindex_virtual_mahalanobis_distance_wrt_vDEC;
                            G[Gindex] = 0.0;

                            for (size_t dIndex : this->dIndex_StateStepLeftInertial_wrt_vDEC)
                            {
                                size_t stateIndex = std::get<1>(this->Derivatives_of_StateStepLeftInertial[dIndex]);
                                size_t Xindex = std::get<0>(this->Derivatives_of_StateStepLeftInertial[dIndex]);

                                double Gentry = dM_dX(stateIndex) * std::get<2>(this->Derivatives_of_StateStepLeftInertial[dIndex]);

                                G[Gindex] += this->X_scale_factors->operator[](Xindex)
                                    * Gentry
                                    / this->myUniverse->LU;
                            }
                        }
                    }//end SphericalRADEC
                    else if (this->myOptions->ParallelShootingStateRepresentation == StateRepresentation::SphericalAZFPA)
                    {
                        //rMag
                        {
                            size_t Gindex = this->Gindex_virtual_mahalanobis_distance_wrt_rMag;
                            G[Gindex] = 0.0;

                            for (size_t dIndex : this->dIndex_StateStepLeftInertial_wrt_r)
                            {
                                size_t stateIndex = std::get<1>(this->Derivatives_of_StateStepLeftInertial[dIndex]);
                                size_t Xindex = std::get<0>(this->Derivatives_of_StateStepLeftInertial[dIndex]);

                                double Gentry = dM_dX(stateIndex) * std::get<2>(this->Derivatives_of_StateStepLeftInertial[dIndex]);

                                G[Gindex] += this->X_scale_factors->operator[](Xindex)
                                    * Gentry
                                    / this->myUniverse->LU;
                            }
                        }
                        //RA
                        {
                            size_t Gindex = this->Gindex_virtual_mahalanobis_distance_wrt_RA;
                            G[Gindex] = 0.0;

                            for (size_t dIndex : this->dIndex_StateStepLeftInertial_wrt_RA)
                            {
                                size_t stateIndex = std::get<1>(this->Derivatives_of_StateStepLeftInertial[dIndex]);
                                size_t Xindex = std::get<0>(this->Derivatives_of_StateStepLeftInertial[dIndex]);

                                double Gentry = dM_dX(stateIndex) * std::get<2>(this->Derivatives_of_StateStepLeftInertial[dIndex]);

                                G[Gindex] += this->X_scale_factors->operator[](Xindex)
                                    * Gentry
                                    / this->myUniverse->LU;
                            }
                        }
                        //DEC
                        {
                            size_t Gindex = this->Gindex_virtual_mahalanobis_distance_wrt_DEC;
                            G[Gindex] = 0.0;

                            for (size_t dIndex : this->dIndex_StateStepLeftInertial_wrt_DEC)
                            {
                                size_t stateIndex = std::get<1>(this->Derivatives_of_StateStepLeftInertial[dIndex]);
                                size_t Xindex = std::get<0>(this->Derivatives_of_StateStepLeftInertial[dIndex]);

                                double Gentry = dM_dX(stateIndex) * std::get<2>(this->Derivatives_of_StateStepLeftInertial[dIndex]);

                                G[Gindex] += this->X_scale_factors->operator[](Xindex)
                                    * Gentry
                                    / this->myUniverse->LU;
                            }
                        }
                        //vMag
                        {
                            size_t Gindex = this->Gindex_virtual_mahalanobis_distance_wrt_vMag;
                            G[Gindex] = 0.0;

                            for (size_t dIndex : this->dIndex_StateStepLeftInertial_wrt_v)
                            {
                                size_t stateIndex = std::get<1>(this->Derivatives_of_StateStepLeftInertial[dIndex]);
                                size_t Xindex = std::get<0>(this->Derivatives_of_StateStepLeftInertial[dIndex]);

                                double Gentry = dM_dX(stateIndex) * std::get<2>(this->Derivatives_of_StateStepLeftInertial[dIndex]);

                                G[Gindex] += this->X_scale_factors->operator[](Xindex)
                                    * Gentry
                                    / this->myUniverse->LU;
                            }
                        }
                        //AZ
                        {
                            size_t Gindex = this->Gindex_virtual_mahalanobis_distance_wrt_AZ;
                            G[Gindex] = 0.0;

                            for (size_t dIndex : this->dIndex_StateStepLeftInertial_wrt_AZ)
                            {
                                size_t stateIndex = std::get<1>(this->Derivatives_of_StateStepLeftInertial[dIndex]);
                                size_t Xindex = std::get<0>(this->Derivatives_of_StateStepLeftInertial[dIndex]);

                                double Gentry = dM_dX(stateIndex) * std::get<2>(this->Derivatives_of_StateStepLeftInertial[dIndex]);

                                G[Gindex] += this->X_scale_factors->operator[](Xindex)
                                    * Gentry
                                    / this->myUniverse->LU;
                            }
                        }
                        //FPA
                        {
                            size_t Gindex = this->Gindex_virtual_mahalanobis_distance_wrt_FPA;
                            G[Gindex] = 0.0;

                            for (size_t dIndex : this->dIndex_StateStepLeftInertial_wrt_FPA)
                            {
                                size_t stateIndex = std::get<1>(this->Derivatives_of_StateStepLeftInertial[dIndex]);
                                size_t Xindex = std::get<0>(this->Derivatives_of_StateStepLeftInertial[dIndex]);

                                double Gentry = dM_dX(stateIndex) * std::get<2>(this->Derivatives_of_StateStepLeftInertial[dIndex]);

                                G[Gindex] += this->X_scale_factors->operator[](Xindex)
                                    * Gentry
                                    / this->myUniverse->LU;
                            }
                        }
                    }//end SphericalAZFPA

                    //mass
                    {
                        size_t Gindex = this->Gindex_virtual_mahalanobis_distance_wrt_mass;
                        G[Gindex] = 0.0;

                        size_t dIndex = this->dIndex_mass_wrt_mass;

                        size_t stateIndex = std::get<1>(this->Derivatives_of_StateStepLeftInertial[dIndex]);
                        size_t Xindex = std::get<0>(this->Derivatives_of_StateStepLeftInertial[dIndex]);

                        double Gentry = dM_dX(stateIndex) * std::get<2>(this->Derivatives_of_StateStepLeftInertial[dIndex]);

                        G[Gindex] += this->X_scale_factors->operator[](Xindex)
                            * Gentry
                            / this->myUniverse->LU;
                    }

                    //Step 5.4: derivatives with respect to time variables affecting step left epoch
                    {
                        double dMdt = 0.0;

                        //contribution from x_c = .norm() is used to convert from 1x1 matrix to scalar
                        dMdt += -(dM_dX * waypoint_state_derivatives).norm();

                        //contribution from CovarianceInverse
                        for (size_t i = 0; i < 7; ++i)
                        {
                            for (size_t j = 0; j < 7; ++j)
                            {
                                dMdt += 1.0 / (2.0 * M) * (WaypointError(i) * WaypointError(j)) _GETVALUE * dCovarianceInverse_depoch(i, j);
                            }
                        }

                        for (size_t tIndex = 0; tIndex < this->Xindices_EventLeftEpoch.size(); ++tIndex)
                        {
                            size_t Gindex = this->Gindex_virtual_mahalanobis_distance_wrt_Time[tIndex];
                            size_t Xindex = this->jGvar->operator[](Gindex);

                            G[Gindex] = this->X_scale_factors->operator[](Xindex)
                                * dMdt;
                        }
                        //last time entry is phase flight time
                        G[this->Gindex_virtual_mahalanobis_distance_wrt_Time.back()] *= this->dStepTime_dPhaseFlightTime;
                    }
                }
                catch (std::exception &error)
                {
                    //if you ran off the end of a spline, the derivatives are all zero
                    //rMAG
                    {
                        size_t Gindex = this->Gindex_virtual_mahalanobis_distance_wrt_rMag;
                        size_t Xindex = this->jGvar->operator[](Gindex);

                        G[Gindex] = 0.0;
                    }
                    //RA
                    {
                        size_t Gindex = this->Gindex_virtual_mahalanobis_distance_wrt_RA;
                        size_t Xindex = this->jGvar->operator[](Gindex);

                        G[Gindex] = 0.0;
                    }
                    //DEC
                    {
                        size_t Gindex = this->Gindex_virtual_mahalanobis_distance_wrt_DEC;
                        size_t Xindex = this->jGvar->operator[](Gindex);

                        G[Gindex] = 0.0;
                    }

                    //Step 5.4: derivatives with respect to time variables affecting step left epoch
                    for (size_t tIndex = 0; tIndex < this->Xindices_EventLeftEpoch.size(); ++tIndex)
                    {
                        size_t Gindex = this->Gindex_virtual_mahalanobis_distance_wrt_Time[tIndex];
                        size_t Xindex = this->jGvar->operator[](tIndex);

                        G[Gindex] = 0.0;
                    }
                }//end if you ran off the spline
            }//end derivatives
        }//end process_waypoint_tracking()

        void ParallelShootingStep::process_derivative_tuples(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //for each variable affecting left-hand state
            //  form a perturbation vector
            //  multiply the perturbation vector by the STM to get the affect on the right-hand state
            //  apply that perturbation to any constraint whose stateIndex depends on a left state variable affected by this variable
            math::Matrix<double> dLeftState_dDecisionVariable(this->CumulativeAugmentedSTM[0].get_n(), 1, 0.0);
            math::Matrix<double> dRightState_dDecisionVariable(this->CumulativeAugmentedSTM[0].get_n(), 1, 0.0);


            //Step 1: non-time affect on the left state
            for (size_t varIndex = 0; varIndex < this->ListOfVariablesAffectingCurrentStepLeftState.size(); ++varIndex)
            {
                dLeftState_dDecisionVariable.assign_zeros();
                size_t Xindex = this->ListOfVariablesAffectingCurrentStepLeftState[varIndex];

                for (size_t entryIndex = 0; entryIndex < this->DerivativesOfCurrentStepLeftStateByVariable[varIndex].size(); ++entryIndex)
                {
                    size_t leftStateIndex = std::get<0>(this->DerivativesOfCurrentStepLeftStateByVariable[varIndex][entryIndex]);
                    size_t dIndex = std::get<1>(this->DerivativesOfCurrentStepLeftStateByVariable[varIndex][entryIndex]);

                    dLeftState_dDecisionVariable(leftStateIndex) = std::get<2>(this->Derivatives_of_StateStepLeftInertial[dIndex]);
                }

                dRightState_dDecisionVariable = this->CumulativeAugmentedSTM[0] * dLeftState_dDecisionVariable;

                for (size_t entryIndex = 0; entryIndex < this->DerivativesOfCurrentStepRightStateByVariable[varIndex].size(); ++entryIndex)
                {
                    size_t rightStateIndex = std::get<0>(this->DerivativesOfCurrentStepRightStateByVariable[varIndex][entryIndex]);
                    size_t dIndex = std::get<1>(this->DerivativesOfCurrentStepRightStateByVariable[varIndex][entryIndex]);

                    std::get<2>(this->Derivatives_of_StateStepRightInertial[dIndex]) = dRightState_dDecisionVariable(rightStateIndex);
                }
            }

            //Step 2: time affect on the left state
            for (size_t varIndex = 0; varIndex < this->ListOfTimeVariablesAffectingCurrentStepLeftState.size(); ++varIndex)
            {
                dLeftState_dDecisionVariable.assign_zeros();
                size_t Xindex = this->ListOfTimeVariablesAffectingCurrentStepLeftState[varIndex];

                for (size_t entryIndex = 0; entryIndex < this->DerivativesOfCurrentStepLeftStateByTimeVariable[varIndex].size(); ++entryIndex)
                {
                    size_t leftStateIndex = std::get<0>(this->DerivativesOfCurrentStepLeftStateByTimeVariable[varIndex][entryIndex]);
                    size_t dIndex = std::get<1>(this->DerivativesOfCurrentStepLeftStateByTimeVariable[varIndex][entryIndex]);

                    dLeftState_dDecisionVariable(leftStateIndex) = std::get<2>(this->Derivatives_of_StateStepLeftInertial_wrt_Time[dIndex]);
                }

                if (Xindex == this->Xindex_PhaseFlightTime)
                {
                    dLeftState_dDecisionVariable(7) = 0.0;
                    dLeftState_dDecisionVariable(13) = 1.0;
                }

                dRightState_dDecisionVariable = this->CumulativeAugmentedSTM[0] * dLeftState_dDecisionVariable;

                for (size_t entryIndex = 0; entryIndex < this->DerivativesOfCurrentStepRightStateByTimeVariable[varIndex].size(); ++entryIndex)
                {
                    size_t rightStateIndex = std::get<0>(this->DerivativesOfCurrentStepRightStateByTimeVariable[varIndex][entryIndex]);
                    size_t dIndex = std::get<1>(this->DerivativesOfCurrentStepRightStateByTimeVariable[varIndex][entryIndex]);

                    std::get<2>(this->Derivatives_of_StateStepRightInertial_wrt_Time[dIndex]) = dRightState_dDecisionVariable(rightStateIndex);
                }
            }

            //Step 4: control
            for (size_t subStepIndex = 0; subStepIndex < this->num_interior_control_points; ++subStepIndex)
            {
                for (size_t controlIndex = 0; controlIndex < this->num_controls; ++controlIndex)
                {
                    for (size_t stateIndex : {0, 1, 2, 3, 4, 5, 6, 9})//no affect on epoch or chemical fuel
                    {
                        size_t dIndex = this->dIndex_right_state_wrt_control[subStepIndex][stateIndex][controlIndex];

                        std::get<2>(this->Derivatives_of_StateStepRightInertial[dIndex]) = this->CumulativeAugmentedSTM[subStepIndex](stateIndex, 10 + controlIndex);
                    }
                }
            }
        }//end process_derivative_tuples
    }//close namespace Phases
}//close namespace EMTG