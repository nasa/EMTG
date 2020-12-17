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
#include "StateRepresentationFactory.h"


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

        ParallelShootingStep::~ParallelShootingStep()
        {
            delete this->myStateRepresentation;
        }//end destructor

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

            //create the state representation
            this->myStateRepresentation = Astrodynamics::CreateStateRepresentation(this->myOptions->ParallelShootingStateRepresentation, this->myUniverse->mu);

            //integrator stuff
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
            this->num_interior_control_points = this->myJourneyOptions->num_interior_control_points;

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
            this->StateStepLeftEncoded.resize(6, 1, 0.0);

            //size various things
            this->throttle.resize(this->num_interior_control_points);
            this->ControlVector.resize(this->num_interior_control_points, math::Matrix<doubleType>(this->num_controls, 1, 0.0));

            //distance contraints
            std::string shortprefix = "p" + std::to_string(this->phaseIndex);
            for (std::string& constraintDefinition : this->myJourneyOptions->PhaseDistanceConstraintDefinitions)
            {
                if (constraintDefinition.find("#") != 0
                    && (constraintDefinition.find(shortprefix) < 1024 || (constraintDefinition.find("pEnd") < 1024 && this->myPhase->getIsLastPhaseInJourney())))
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
                    if (constraintDefinition.find(shortprefix) < 1024 
                        || (constraintDefinition.find("pEnd") < 1024 && this->myPhase->getIsLastPhaseInJourney()))
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

            //left-hand side truth tables (default values, 6-state gets done in calcbounds_step_left_state())
            //almost everything is zero
            this->TruthTable_StateStepLeftCartesian_wrt_StateStepLeftEncoded.resize(10, std::vector<bool>(10, false));
            //mass and tanks are pass throughs
            this->TruthTable_StateStepLeftCartesian_wrt_StateStepLeftEncoded[6][6] = true;
            this->TruthTable_StateStepLeftCartesian_wrt_StateStepLeftEncoded[8][8] = true;
            this->TruthTable_StateStepLeftCartesian_wrt_StateStepLeftEncoded[9][9] = true;

            //right-hand side truth tables
            this->TruthTable_StateStepRightInertial_wrt_StateStepLeftInertial.resize(10, std::vector<bool>(10, true));
            this->TruthTable_StateStepRightInertial_wrt_Control.resize(10, true);
            this->TruthTable_StateStepRightInertial_wrt_PhaseFlightTime.resize(10, true);

            //epoch has derivatives with respect only to itself
            this->TruthTable_StateStepRightInertial_wrt_StateStepLeftInertial[7] = std::vector<bool>(10, false);
            this->TruthTable_StateStepRightInertial_wrt_StateStepLeftInertial[7][7] = true;

            //position, velocity, and mass have no derivative with respect to either tank variable
            for (size_t stateIndex = 0; stateIndex < 7; ++stateIndex)
            {
                this->TruthTable_StateStepRightInertial_wrt_StateStepLeftInertial[stateIndex][8] = false;
                this->TruthTable_StateStepRightInertial_wrt_StateStepLeftInertial[stateIndex][9] = false;
            }

            //chemical fuel only has a derivative with respect to itself and epoch, and even then only if ACS is on
            this->TruthTable_StateStepRightInertial_wrt_StateStepLeftInertial[8] = std::vector<bool>(10, false);
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
            //Step 1: create a vector of state names, bounds, and scale factors for the state representation of choice
            
            switch (this->myOptions->ParallelShootingStateRepresentation)
            {
            case StateRepresentation::SphericalRADEC:
            {
                this->statesToRepresent.push_back({ "r",  this->myUniverse->central_body.radius + this->myUniverse->central_body.minimum_safe_flyby_altitude, this->myUniverse->r_SOI, this->myUniverse->LU, false });
                this->statesToRepresent.push_back({ "RA",  -8.0 * math::PI, 8.0 * math::PI, 1.0, true });
                this->statesToRepresent.push_back({ "DEC",  -math::PIover2, math::PIover2, 1.0, true });
                this->statesToRepresent.push_back({ "v",  0.0, 10.0 * this->myUniverse->LU / this->myUniverse->TU, this->myUniverse->LU / this->myUniverse->TU, false });
                this->statesToRepresent.push_back({ "vRA",  -8.0 * math::PI, 8.0 * math::PI, 1.0, true });
                this->statesToRepresent.push_back({ "vDEC", -8.0 * math::PI, 8.0 * math::PI, 1.0, true });

                //the first three variables affect all position states
                for (size_t stateIndex : { 0, 1, 2 })
                {
                    for (size_t varIndex : { 0, 1, 2 })
                    {
                        this->TruthTable_StateStepLeftCartesian_wrt_StateStepLeftEncoded[stateIndex][varIndex] = true;
                    }
                }
                //the last three variables affect all velocity states
                for (size_t stateIndex : { 3, 4, 5 })
                {
                    for (size_t varIndex : { 3, 4, 5 })
                    {
                        this->TruthTable_StateStepLeftCartesian_wrt_StateStepLeftEncoded[stateIndex][varIndex] = true;
                    }
                }

                break;
            }
            case StateRepresentation::SphericalAZFPA:
            {
                this->statesToRepresent.push_back({ "r",  this->myUniverse->central_body.radius + this->myUniverse->central_body.minimum_safe_flyby_altitude, this->myUniverse->r_SOI, this->myUniverse->LU, false });
                this->statesToRepresent.push_back({ "RA",  -8.0 * math::PI, 8.0 * math::PI, 1.0, true });
                this->statesToRepresent.push_back({ "DEC",  -math::PIover2, math::PIover2, 1.0, true });
                this->statesToRepresent.push_back({ "v",  0.0, 10.0 * this->myUniverse->LU / this->myUniverse->TU, this->myUniverse->LU / this->myUniverse->TU, false });
                this->statesToRepresent.push_back({ "AZ",  -8.0 * math::PI, 8.0 * math::PI, 1.0, true });
                this->statesToRepresent.push_back({ "FPA",  -8.0 * math::PI, 8.0 * math::PI, 1.0, true });

                //the first three variables affect all position states
                for (size_t stateIndex : { 0, 1, 2 })
                {
                    for (size_t varIndex : { 0, 1, 2 })
                    {
                        this->TruthTable_StateStepLeftCartesian_wrt_StateStepLeftEncoded[stateIndex][varIndex] = true;
                    }
                }
                //the last three variables affect all velocity states
                for (size_t stateIndex : { 3, 4, 5 })
                {
                    for (size_t varIndex : { 3, 4, 5 })
                    {
                        this->TruthTable_StateStepLeftCartesian_wrt_StateStepLeftEncoded[stateIndex][varIndex] = true;
                    }
                }

                break;
            }
            case StateRepresentation::Cartesian:
            {
                this->statesToRepresent.push_back({ "x",  -this->myUniverse->r_SOI, this->myUniverse->r_SOI, this->myUniverse->LU, false });
                this->statesToRepresent.push_back({ "y",  -this->myUniverse->r_SOI, this->myUniverse->r_SOI, this->myUniverse->LU, false });
                this->statesToRepresent.push_back({ "z",  -this->myUniverse->r_SOI, this->myUniverse->r_SOI, this->myUniverse->LU, false });
                this->statesToRepresent.push_back({ "vx", -10.0 * this->myUniverse->LU / this->myUniverse->TU, 10.0 * this->myUniverse->LU / this->myUniverse->TU, this->myUniverse->LU / this->myUniverse->TU, false });
                this->statesToRepresent.push_back({ "vy", -10.0 * this->myUniverse->LU / this->myUniverse->TU, 10.0 * this->myUniverse->LU / this->myUniverse->TU, this->myUniverse->LU / this->myUniverse->TU, false });
                this->statesToRepresent.push_back({ "vz", -10.0 * this->myUniverse->LU / this->myUniverse->TU, 10.0 * this->myUniverse->LU / this->myUniverse->TU, this->myUniverse->LU / this->myUniverse->TU, false });

                //encoded states match 1:1 with cartesian states
                for (size_t stateIndex : { 0, 1, 2, 3, 4, 5 })
                {
                    size_t varIndex = stateIndex;
                    
                    this->TruthTable_StateStepLeftCartesian_wrt_StateStepLeftEncoded[stateIndex][varIndex] = true;
                }
                break;
            }
            case StateRepresentation::COE:
            {
                this->statesToRepresent.push_back({ "SMA",  -this->myUniverse->r_SOI, this->myUniverse->r_SOI, this->myUniverse->LU, false });
                this->statesToRepresent.push_back({ "ECC",  math::SMALL, 10.0, 1.0, false });//perfectly zero eccentricity can cause kablooey
                this->statesToRepresent.push_back({ "INC",   -math::PI, math::PI, 1.0, true });
                this->statesToRepresent.push_back({ "RAAN",  -8.0 * math::PI, 8.0 * math::PI, 1.0, true });
                this->statesToRepresent.push_back({ "AOP",  -8.0 * math::PI, 8.0 * math::PI, 1.0, true });
                this->statesToRepresent.push_back({ "TA",  -8.0 * math::PI, 8.0 * math::PI, 1.0, true });

                //all states affect everything
                for (size_t stateIndex : { 0, 1, 2, 3, 4, 5 })
                {
                    for (size_t varIndex : { 0, 1, 2, 3, 4, 5 })
                    {
                        this->TruthTable_StateStepLeftCartesian_wrt_StateStepLeftEncoded[stateIndex][varIndex] = true;
                    }
                }

                break;
            }
            case StateRepresentation::MEE:
            {
                this->statesToRepresent.push_back({ "P",  math::SMALL, this->myUniverse->r_SOI, this->myUniverse->LU, false });
                this->statesToRepresent.push_back({ "F",  -10.0, 10.0, 1.0, false });
                this->statesToRepresent.push_back({ "G",  -10.0, 10.0, 1.0, false });
                this->statesToRepresent.push_back({ "H",  -10.0, 10.0, 1.0, false });
                this->statesToRepresent.push_back({ "K",  -10.0, 10.0, 1.0, false });
                this->statesToRepresent.push_back({ "L",  -8.0 * math::PI, 8.0 * math::PI, 1.0, true });

                //all states affect everything
                for (size_t stateIndex : { 0, 1, 2, 3, 4, 5 })
                {
                    for (size_t varIndex : { 0, 1, 2, 3, 4, 5 })
                    {
                        this->TruthTable_StateStepLeftCartesian_wrt_StateStepLeftEncoded[stateIndex][varIndex] = true;
                    }
                }

                break;
            }
            default:
                throw std::invalid_argument("ParallelShootingStep does not recognize state representation " + std::to_string(this->myOptions->ParallelShootingStateRepresentation) + ". Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + ".");
            }

            //Step 2: create the variables and appropriate indices!
            for (std::tuple<std::string, double, double, double, bool>& State :  this->statesToRepresent)
            {

                this->Xdescriptions->push_back(prefix + "left state " + std::get<0>(State));
                this->Xlowerbounds->push_back(std::get<1>(State));
                this->Xupperbounds->push_back(std::get<2>(State));
                this->X_scale_factors->push_back(std::get<3>(State));
                this->encodedStateIsAngle.push_back(std::get<4>(State));
                this->Xindex_state_elements.push_back(this->Xdescriptions->size() - 1);

                std::vector<size_t> dIndex_StateStepLeftInertial_wrt_thisStateElement;

                size_t varIndex = this->Xindex_state_elements.size() - 1;

                //create a derivative entry if one exists
                for (size_t stateIndex : {0, 1, 2, 3, 4, 5})
                {
                    if (this->TruthTable_StateStepLeftCartesian_wrt_StateStepLeftEncoded[stateIndex][varIndex])
                    {
                        this->Derivatives_of_StateStepLeftInertial.push_back(std::make_tuple(this->Xdescriptions->size() - 1, stateIndex, 1.0));
                        dIndex_StateStepLeftInertial_wrt_thisStateElement.push_back(this->Derivatives_of_StateStepLeftInertial.size() - 1);
                    }
                    else //we have to keep the size of the dIndex matrix constant
                    {
                        dIndex_StateStepLeftInertial_wrt_thisStateElement.push_back(32767);
                    }
                }

                this->dIndex_StateStepLeftInertial_wrt_StateElements.push_back(dIndex_StateStepLeftInertial_wrt_thisStateElement);
            }//end loop over state variable creation

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

            //Step 7: create the continuity constraint scale factors, which need to be done here because we know the decison variable types
            this->continuity_constraint_scale_factors = this->myPhase->get_continuity_constraint_scale_factors();

            if (this->myOptions->ParallelShootingConstraintStateRepresentation == ConstraintStateRepresentation::Native)
            {
                for (size_t stateIndex : {0, 1, 2, 3, 4, 5})
                {
                    this->continuity_constraint_scale_factors(stateIndex) = (1.0 / std::get<3>(this->statesToRepresent[stateIndex]));
                }
            }


        }//end calcbounds_step_left_state()
        
        void ParallelShootingStep::calculate_dependencies_epoch_time()
        {
            //the left epoch of this event depends on all epoch and time variables before it
            std::vector<size_t> timeVariables = this->myPhase->getArrivalEvent()->get_Xindices_EventRightEpoch();
            for (size_t Xindex : timeVariables)
            {
                this->Xindices_EventLeftEpoch.push_back(Xindex);
                this->Derivatives_of_StateStepLeftInertial_wrt_Time.push_back({ Xindex, 7, 1.0 });
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
            std::vector<std::string> matchPointConstraintNames = this->myPhase->get_matchPointConstraintNames();
            
            if (this->myOptions->ParallelShootingConstraintStateRepresentation == ConstraintStateRepresentation::Native)
            {
                //override the 6-state names with whatever 6-state we're using
                for (size_t stateIndex : {0, 1, 2, 3, 4, 5})
                {
                    matchPointConstraintNames[stateIndex] = std::get<0>(this->statesToRepresent[stateIndex]);
                }
            }

            for (size_t constraintIndex = 0; constraintIndex < this->numMatchConstraints; ++constraintIndex)
            {
                Flowerbounds->push_back(-math::SMALL);
                Fupperbounds->push_back(math::SMALL);
                Fdescriptions->push_back(prefix + "left match point " + matchPointConstraintNames[constraintIndex]);
                this->Findices_left_match_point_constraints.push_back(Fdescriptions->size() - 1);
            }//end construction of match point constraints

            //Step 4: derivatives with respect to previous step right boundary
            {
                //Step 4.1: non-time
                this->Gindices_LeftMatchPoint6state_wrt_PreviousStepRightState.resize(this->ListOfVariablesAffectingPreviousStepRightState.size());
                this->Gindices_LeftMatchPointMassConstraints_wrt_PreviousStepRightState.resize(this->ListOfVariablesAffectingPreviousStepRightState.size());
                for (size_t varIndex = 0; varIndex < this->ListOfVariablesAffectingPreviousStepRightState.size(); ++varIndex)
                {
                    size_t Xindex = this->ListOfVariablesAffectingPreviousStepRightState[varIndex];

                    for (size_t entryIndex = 0; entryIndex < this->DerivativesOfPreviousStepRightStateByVariable[varIndex].size(); ++entryIndex)
                    {
                        size_t stateIndex = std::get<0>(this->DerivativesOfPreviousStepRightStateByVariable[varIndex][entryIndex]);

                        if (stateIndex < 6) //6-state constraint
                        {
                            //create entries in G for all six states wrt this variable and store their Gindex
                            for (size_t constraintIndex : {0, 1, 2, 3, 4, 5})
                            {
                                this->create_sparsity_entry(this->Findices_left_match_point_constraints[constraintIndex],
                                    Xindex,
                                    this->Gindices_LeftMatchPoint6state_wrt_PreviousStepRightState[varIndex]);
                            }
                        }
                        else //other constraints, such as the mass and tank constraints
                        {
                            //create entry in G for this specific state wrt this variable and store its Gindex
                            size_t constraintIndex = stateIndex < 7 ? stateIndex : stateIndex - 1;

                            this->create_sparsity_entry(this->Findices_left_match_point_constraints[constraintIndex],
                                Xindex,
                                this->Gindices_LeftMatchPointMassConstraints_wrt_PreviousStepRightState[varIndex]);
                        }
                    }
                }//end loop over non-time variables

                //Step 4.2: time
                this->Gindices_LeftMatchPoint6state_wrt_PreviousStepRightStateTime.resize(this->ListOfTimeVariablesAffectingPreviousStepRightState.size());
                this->Gindices_LeftMatchPointMassConstraints_wrt_PreviousStepRightStateTime.resize(this->ListOfTimeVariablesAffectingPreviousStepRightState.size());
                for (size_t varIndex = 0; varIndex < this->ListOfTimeVariablesAffectingPreviousStepRightState.size(); ++varIndex)
                {
                    size_t Xindex = this->ListOfTimeVariablesAffectingPreviousStepRightState[varIndex];

                    for (size_t entryIndex = 0; entryIndex < this->DerivativesOfPreviousStepRightStateByTimeVariable[varIndex].size(); ++entryIndex)
                    {
                        size_t stateIndex = std::get<0>(this->DerivativesOfPreviousStepRightStateByTimeVariable[varIndex][entryIndex]);

                        if (stateIndex < 6) //6-state constraint
                        {
                            //create entries in G for all six states wrt this variable and store their Gindex
                            for (size_t constraintIndex : {0, 1, 2, 3, 4, 5})
                            {
                                this->create_sparsity_entry(this->Findices_left_match_point_constraints[constraintIndex],
                                    Xindex,
                                    this->Gindices_LeftMatchPoint6state_wrt_PreviousStepRightStateTime[varIndex]);
                            }
                        }
                        else if (stateIndex != 7) //other constraints, such as the mass and tank constraints - note that there is no time constraint
                        {
                            //create entry in G for this specific state wrt this variable and store its Gindex
                            size_t constraintIndex = stateIndex < 7 ? stateIndex : stateIndex - 1;

                            this->create_sparsity_entry(this->Findices_left_match_point_constraints[constraintIndex],
                                Xindex,
                                this->Gindices_LeftMatchPointMassConstraints_wrt_PreviousStepRightStateTime[varIndex]);
                        }
                    }
                }//end loop over time variables
            }//end derivatives with respect to previous step right boundary

            //Step 5: derivatives with respect to current step left state.
            if (this->myOptions->ParallelShootingConstraintStateRepresentation == ConstraintStateRepresentation::CartesianConstraint)
            {
                for (size_t encodedStateIndex : {0, 1, 2, 3, 4, 5})
                {
                    for (size_t cartesianStateIndex : {0, 1, 2, 3, 4, 5})
                    {
                        if (this->TruthTable_StateStepLeftCartesian_wrt_StateStepLeftEncoded[cartesianStateIndex][encodedStateIndex])
                        {
                            size_t dIndex = this->dIndex_StateStepLeftInertial_wrt_StateElements[encodedStateIndex][cartesianStateIndex];

                            size_t Xindex = std::get<0>(this->Derivatives_of_StateStepLeftInertial[dIndex]);

                            this->create_sparsity_entry(this->Findices_left_match_point_constraints[cartesianStateIndex],
                                Xindex,
                                this->Gindices_StepLeftMatchPoint_wrt_StepLeftEncodedState);
                        }
                    }
                }
            }
            else if (this->myOptions->ParallelShootingConstraintStateRepresentation == ConstraintStateRepresentation::Native)
            {
                //We can assume that the encoded state is in the same state representation as the constraint, so these are linear
                for (size_t stateIndex : {0, 1, 2, 3, 4, 5})
                {
                    size_t Xindex = this->Xindex_state_elements[stateIndex];

                    this->create_sparsity_entry(this->Findices_left_match_point_constraints[stateIndex],
                        Xindex,
                        this->Gindices_StepLeftMatchPoint_wrt_StepLeftEncodedState);
                }
            }
            else
            {
                throw std::invalid_argument("Parallel shooting constraint state representation not recognized: " + std::to_string(this->myOptions->ParallelShootingConstraintStateRepresentation));
            }

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

        void ParallelShootingStep::process_step_left_state(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //Step 1: extract the states            
            //6-state
            for (size_t stateIndex : {0, 1, 2, 3, 4, 5})
            {
                this->StateStepLeftEncoded(stateIndex) = X[Xindex++];
            }
            
            //masses
            this->StateStepLeftInertial(6) = X[Xindex++];//mass
            this->StateStepLeftInertial(8) = X[Xindex++];//chemical fuel
            this->StateStepLeftInertial(9) = X[Xindex++];//electric propellant

            //Step 2: convert to cartesian
            math::Matrix<doubleType> stateLeftCartesianRepresentation = this->myStateRepresentation->convertFromRepresentationToCartesian(this->StateStepLeftEncoded, needG);

            for (size_t stateIndex : {0, 1, 2, 3, 4, 5})
            {
                this->StateStepLeftInertial(stateIndex) = stateLeftCartesianRepresentation(stateIndex);
            }

            //Step 3: derivatives
            if (needG)
            {
                //Step 3.1: extract the transformation matrix
                math::Matrix<doubleType> TransformationMatrix = this->myStateRepresentation->getRepresentationToCartesianTransitionMatrix();

                //Step 3.2: insert into the derivative tuples
                for (size_t nativeStateIndex : {0, 1, 2, 3, 4, 5})
                {
                    for (size_t cartesianStateIndex : {0, 1, 2, 3, 4, 5})
                    {
                        if (this->TruthTable_StateStepLeftCartesian_wrt_StateStepLeftEncoded[cartesianStateIndex][nativeStateIndex])
                        {
                            size_t dIndex = this->dIndex_StateStepLeftInertial_wrt_StateElements[nativeStateIndex][cartesianStateIndex];

                            std::get<2>(this->Derivatives_of_StateStepLeftInertial[dIndex]) = TransformationMatrix(cartesianStateIndex, nativeStateIndex)_GETVALUE;
                        }
                    }
                }
            }

             //Step 4: epoch
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
            //Step 1.1: get the necessary data
            std::vector<size_t>& matchPointConstraintStateIndex = this->myPhase->get_matchPointConstraintStateIndex();
            math::Matrix<doubleType>& PreviousStepRightInertial = this->previousStep->get_StateStepRightInertial();
            math::Matrix<doubleType> StateTransformationMatrix;
            math::Matrix<doubleType> PreviousStepRightNative;

            //Step 1.2: process the constraint
            //are we doing this in cartesian or the ParallelShootingPhase's "native" state representation
            //if cartesian, we'll use the "inertial" representations
            //if "native", we'll use the "encoded" representation for the current step and convert the previous step's inertial state back into native state representation
            if (this->myOptions->ParallelShootingConstraintStateRepresentation == ConstraintStateRepresentation::CartesianConstraint)
            {
                StateTransformationMatrix = math::Matrix<doubleType>(6, math::identity);

                for (size_t constraintIndex = 0; constraintIndex < this->numMatchConstraints; ++constraintIndex)
                {
                    size_t stateIndex = matchPointConstraintStateIndex[constraintIndex];
                    F[Findex++] = (this->StateStepLeftInertial(stateIndex) - PreviousStepRightInertial(stateIndex))
                        * continuity_constraint_scale_factors(constraintIndex);
                }
            }
            else if (this->myOptions->ParallelShootingConstraintStateRepresentation == ConstraintStateRepresentation::Native)
            {
                PreviousStepRightNative = this->myStateRepresentation->convertFromCartesianToRepresentation(PreviousStepRightInertial.getSubMatrix1D(0, 5), needG);
                StateTransformationMatrix = this->myStateRepresentation->getCartesianToRepresentationTransitionMatrix();

                //6-state
                for (size_t constraintIndex = 0; constraintIndex < 6; ++constraintIndex)
                {
                    size_t stateIndex = matchPointConstraintStateIndex[constraintIndex];

                    if (this->encodedStateIsAngle[stateIndex])
                    {
                        F[Findex++] = (math::fmod(this->StateStepLeftEncoded(stateIndex), math::TwoPI) - math::fmod(PreviousStepRightNative(stateIndex), math::TwoPI))
                            * continuity_constraint_scale_factors(constraintIndex);
                    }
                    else
                    {
                        F[Findex++] = (this->StateStepLeftEncoded(stateIndex) - PreviousStepRightNative(stateIndex))
                            * continuity_constraint_scale_factors(constraintIndex);
                    }
                }

                for (size_t constraintIndex = 6; constraintIndex < this->numMatchConstraints; ++constraintIndex)
                {
                    size_t stateIndex = matchPointConstraintStateIndex[constraintIndex];
                    F[Findex++] = (this->StateStepLeftInertial(stateIndex) - PreviousStepRightInertial(stateIndex))
                        * continuity_constraint_scale_factors(constraintIndex);
                }
            }
            else
            {
                throw std::invalid_argument("Parallel shooting constraint state representation not recognized: " + std::to_string(this->myOptions->ParallelShootingConstraintStateRepresentation));
            }

            //Step 2: derivatives of the match point constraints
            if (needG)
            {
                math::Matrix<double> dMatchState_dDecisionVariable(this->StateStepLeftInertial.get_n(), 1, 0.0);

                //Step 2.1: with respect to the current step 6-state
                if (this->myOptions->ParallelShootingConstraintStateRepresentation == ConstraintStateRepresentation::CartesianConstraint)
                {
                    //first clear the entries
                    for (size_t& Gindex : this->Gindices_StepLeftMatchPoint_wrt_StepLeftEncodedState)
                    {
                        G[Gindex] = 0.0;
                    }

                    //now populate them
                    for (size_t encodedStateIndex : {0, 1, 2, 3, 4, 5})
                    {
                        for (size_t cartesianStateIndex : {0, 1, 2, 3, 4, 5})
                        {
                            if (this->TruthTable_StateStepLeftCartesian_wrt_StateStepLeftEncoded[cartesianStateIndex][encodedStateIndex])
                            {
                                size_t dIndex = this->dIndex_StateStepLeftInertial_wrt_StateElements[encodedStateIndex][cartesianStateIndex];

                                size_t Gindex = this->Gindices_StepLeftMatchPoint_wrt_StepLeftEncodedState[dIndex];

                                size_t Xindex = std::get<0>(this->Derivatives_of_StateStepLeftInertial[dIndex]);

                                G[Gindex] = +this->X_scale_factors->operator[](Xindex)
                                    * std::get<2>(this->Derivatives_of_StateStepLeftInertial[dIndex])
                                    * continuity_constraint_scale_factors(cartesianStateIndex);
                            }
                        }
                    }
                }
                else if (this->myOptions->ParallelShootingConstraintStateRepresentation == ConstraintStateRepresentation::Native)
                {
                    //We can assume that the encoded state is in the same state representation as the constraint, so these are linear
                    for (size_t encodedStateIndex : {0, 1, 2, 3, 4, 5})
                    {
                        size_t dIndex = encodedStateIndex;

                        size_t Gindex = this->Gindices_StepLeftMatchPoint_wrt_StepLeftEncodedState[dIndex];

                        size_t Xindex = this->jGvar->operator[](Gindex);

                        G[Gindex] = this->X_scale_factors->operator[](Xindex)
                            * continuity_constraint_scale_factors(encodedStateIndex);
                    }
                }
                else
                {
                    throw std::invalid_argument("Parallel shooting constraint state representation not recognized: " + std::to_string(this->myOptions->ParallelShootingConstraintStateRepresentation));
                }

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
                    size_t MassEntryCounter = 0;
                    size_t SixStateEntryCounter = 0;

                    size_t Xindex = this->ListOfVariablesAffectingPreviousStepRightState[varIndex];
                    
                    //clear the existing 6-state entries, because we're going to build them up additively
                    for (size_t Gindex : this->Gindices_LeftMatchPoint6state_wrt_PreviousStepRightState[varIndex])
                    {
                        G[Gindex] = 0.0;
                    }

                    for (size_t entryIndex = 0; entryIndex < this->DerivativesOfPreviousStepRightStateByVariable[varIndex].size(); ++entryIndex)
                    {
                        size_t stateIndex = std::get<0>(this->DerivativesOfPreviousStepRightStateByVariable[varIndex][entryIndex]);
                        if (stateIndex < 6) //6-state
                        {
                            for (size_t constraintIndex : {0, 1, 2, 3, 4, 5})
                            {
                                //derivative entry is P(newStateIndex, oldStateIndex) * derivative
                                size_t dIndex = std::get<1>(this->DerivativesOfPreviousStepRightStateByVariable[varIndex][entryIndex]);
                                
                                double TheDerivative = std::get<2>(Derivatives_of_PreviousStepRightState[dIndex]) * StateTransformationMatrix(constraintIndex, stateIndex)_GETVALUE;

                                size_t Gindex = this->Gindices_LeftMatchPoint6state_wrt_PreviousStepRightState[varIndex][SixStateEntryCounter++];
                                size_t Xindex = this->jGvar->operator[](Gindex);

                                G[Gindex] += this->X_scale_factors->operator[](Xindex)
                                    * -TheDerivative
                                    * continuity_constraint_scale_factors(constraintIndex);
                            }
                        }//end 6-state entries
                        else //other constraints, such as mass and tanks
                        {
                            //directly populate the entry in G with the derivative entry
                            size_t dIndex = std::get<1>(this->DerivativesOfPreviousStepRightStateByVariable[varIndex][entryIndex]);

                            size_t constraintIndex = stateIndex < 7 ? stateIndex : stateIndex - 1;

                            double TheDerivative = std::get<2>(Derivatives_of_PreviousStepRightState[dIndex]);

                            size_t Gindex = this->Gindices_LeftMatchPointMassConstraints_wrt_PreviousStepRightState[varIndex][MassEntryCounter++];
                            size_t Xindex = this->jGvar->operator[](Gindex);

                            G[Gindex] = this->X_scale_factors->operator[](Xindex)
                                * -TheDerivative
                                * continuity_constraint_scale_factors(constraintIndex);
                        }//end non-6-state entries
                    }//end loop over derivative entries of the previous step's right-hand state
                }//end loop over non-time variables

                //Step 2.2.1: time
                for (size_t varIndex = 0; varIndex < this->ListOfTimeVariablesAffectingPreviousStepRightState.size(); ++varIndex)
                {
                    size_t Xindex = this->ListOfTimeVariablesAffectingPreviousStepRightState[varIndex];
                    size_t SixStateEntryCounter = 0;
                    size_t MassEntryCounter = 0;

                    //clear the existing 6-state entries, because we're going to build them up additively
                    for (size_t Gindex : this->Gindices_LeftMatchPoint6state_wrt_PreviousStepRightStateTime[varIndex])
                    {
                        G[Gindex] = 0.0;
                    }

                    for (size_t entryIndex = 0; entryIndex < this->DerivativesOfPreviousStepRightStateByTimeVariable[varIndex].size(); ++entryIndex)
                    {
                        size_t stateIndex = std::get<0>(this->DerivativesOfPreviousStepRightStateByTimeVariable[varIndex][entryIndex]);

                        if (stateIndex < 6) //6-state
                        {
                            for (size_t constraintIndex : {0, 1, 2, 3, 4, 5})
                            {
                                //derivative entry is P(newStateIndex, oldStateIndex) * derivative
                                size_t dIndex = std::get<1>(this->DerivativesOfPreviousStepRightStateByTimeVariable[varIndex][entryIndex]);

                                double TheDerivative = std::get<2>(Derivatives_of_PreviousStepRightState_wrt_Time[dIndex]) * StateTransformationMatrix(constraintIndex, stateIndex)_GETVALUE;

                                size_t Gindex = this->Gindices_LeftMatchPoint6state_wrt_PreviousStepRightStateTime[varIndex][SixStateEntryCounter++];
                                size_t Xindex = this->jGvar->operator[](Gindex);

                                G[Gindex] += this->X_scale_factors->operator[](Xindex)
                                    * -TheDerivative
                                    * continuity_constraint_scale_factors(constraintIndex);
                            }
                        }//end 6-state entries
                        else if (stateIndex != 7) //other constraints, such as mass and tanks - note that time does not have a constraint
                        {
                            //directly populate the entry in G with the derivative entry
                            size_t dIndex = std::get<1>(this->DerivativesOfPreviousStepRightStateByTimeVariable[varIndex][entryIndex]);

                            size_t constraintIndex = stateIndex < 7 ? stateIndex : stateIndex - 1;

                            double TheDerivative = std::get<2>(Derivatives_of_PreviousStepRightState_wrt_Time[dIndex]);

                            size_t Gindex = this->Gindices_LeftMatchPointMassConstraints_wrt_PreviousStepRightStateTime[varIndex][MassEntryCounter++];
                            size_t Xindex = this->jGvar->operator[](Gindex);

                            G[Gindex] = this->X_scale_factors->operator[](Xindex)
                                * -TheDerivative
                                * continuity_constraint_scale_factors(constraintIndex);
                        }//end non-6-state entries
                    }//end loop over derivative entries of the previous step's right-hand state
                }//end loop over non-time variables
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
                    dLeftState_dDecisionVariable(7) = (float)stepIndex / this->myPhase->get_num_steps();
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