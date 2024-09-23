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

#include "SpacecraftOptionsFactory.h"
#include "PowerSystemOptions.h"
#include "PropulsionSystemOptions.h"
#include "StageOptions.h"
#include "SpacecraftOptions.h"

#ifdef HAS_BUILT_IN_THRUSTERS
    #include "engine_model.h"
#endif


namespace EMTG
{
    namespace HardwareModels
    {
        SpacecraftOptions CreateSpacecraftOptions(const missionoptions& options)
        {
            //create an options structure that we will then populate
            SpacecraftOptions mySpacecraftOptions;

            //configure launch vehicle and spacecraft options objects. Right now this just mimics what is EMTGv9.cpp, but eventually we probably want a menu-driven outer-loop transcription
            //configure the LaunchVehicleOptions and SpacecraftOptions objects
            if (options.SpacecraftModelInput == EMTG::SpacecraftModelInputType::AssembleFromLibraries)
            {
                //creates a single-stage spacecraft with the appropriate hardware
                EMTG::HardwareModels::PowerSystemOptionsLibrary myPowerSystemOptionsLibrary(options.HardwarePath + "/" + options.PowerSystemsLibraryFile);
                EMTG::HardwareModels::PropulsionSystemOptionsLibrary myPropulsionSystemOptionsLibrary(options.HardwarePath + "/" + options.PropulsionSystemsLibraryFile);

                EMTG::HardwareModels::PowerSystemOptions myPowerSystemOptions(myPowerSystemOptionsLibrary.getPowerSystem(options.PowerSystemKey));
                myPowerSystemOptions.setPowerMargin(options.power_margin);
                myPowerSystemOptions.setPowerSystemDecayRefEpoch(options.power_system_decay_reference_epoch);
                EMTG::HardwareModels::PropulsionSystemOptions myChemicalPropulsionSystemOptions(myPropulsionSystemOptionsLibrary.getPropulsionSystem(options.ChemicalPropulsionSystemKey));
                EMTG::HardwareModels::PropulsionSystemOptions myElectricPropulsionSystemOptions(myPropulsionSystemOptionsLibrary.getPropulsionSystem(options.ElectricPropulsionSystemKey));

                myElectricPropulsionSystemOptions.setNumberOfStrings(options.number_of_electric_propulsion_systems);
                myElectricPropulsionSystemOptions.setSharpness(options.throttle_sharpness);

                EMTG::HardwareModels::StageOptions myStageOptions(myChemicalPropulsionSystemOptions, myElectricPropulsionSystemOptions, myPowerSystemOptions);

                myStageOptions.setThrottleLogic(options.throttle_logic_mode);
                myStageOptions.setThrottleSharpness(options.throttle_sharpness);

                mySpacecraftOptions = EMTG::HardwareModels::SpacecraftOptions(std::vector<EMTG::HardwareModels::StageOptions>(1, myStageOptions));
                mySpacecraftOptions.setName(options.mission_name);

                //tanks
                if (options.enable_electric_propellant_tank_constraint)
                {
                    mySpacecraftOptions.setEnableGlobalElectricPropellantTankConstraint(true);
                    mySpacecraftOptions.setGlobalElectricPropellantTankCapacity(options.maximum_electric_propellant);
                }
                if (options.enable_chemical_propellant_tank_constraint)
                {
                    mySpacecraftOptions.setEnableGlobalChemicalPropellantTankConstraint(true);
                    mySpacecraftOptions.setGlobalFuelTankCapacity(options.maximum_chemical_fuel);
                    mySpacecraftOptions.setGlobalOxidizerTankCapacity(options.maximum_chemical_oxidizer);
                }
                if (options.constrain_dry_mass)
                {
                    mySpacecraftOptions.setEnableGlobalDryMassConstraint(true);
                    mySpacecraftOptions.setGlobalDryMassBounds(std::vector<double>({ options.final_mass_constraint_bounds[0], options.final_mass_constraint_bounds[1] }));
                }
            }
            else if (options.SpacecraftModelInput == EMTG::SpacecraftModelInputType::ReadSpacecraftFile)
            {
                //ingests a .emtg_spacecraftoptions file with any number of stages
                mySpacecraftOptions.parse_input_file(options.HardwarePath + "/" + options.SpacecraftOptionsFile);

                for (size_t stageIndex = 0; stageIndex < mySpacecraftOptions.getNumberOfStages(); ++stageIndex)
				{
                    //have to make a new StageOptions and assign it
                    //make a new power system options
                    StageOptions myStage(mySpacecraftOptions.getStageOptions(stageIndex));
                    PowerSystemOptions myPowerSystemOptions(myStage.getPowerSystemOptions());
                                        
                    myPowerSystemOptions.setPowerMargin(options.power_margin);
                    myPowerSystemOptions.setPowerSystemDecayRefEpoch(options.power_system_decay_reference_epoch);

                    //put it on a copy of the stage
                    myStage.setPowerSystemOptions(myPowerSystemOptions);

                    //overwrite the old stage
                    mySpacecraftOptions.setStageOptions(stageIndex, myStage);
				}

                if (options.constrain_dry_mass)
                {
                    mySpacecraftOptions.setEnableGlobalDryMassConstraint(true);
                    mySpacecraftOptions.setGlobalDryMassBounds(std::vector<double>({ options.final_mass_constraint_bounds[0], options.final_mass_constraint_bounds[1] }));
                }
            }
            else// (options.SpacecraftModelInput == EMTG::SpacecraftModelInputType::AssembleFromMissionOptions)
            {
                //creates a single-stage spacecraft with the appropriate hardware

                //then configure the power system
                PowerSystemOptions myPowerSystemOptions;
                myPowerSystemOptions.setName("PowerFromMissionOptions");
                myPowerSystemOptions.setP0(options.power_at_1_AU);
                myPowerSystemOptions.setBusPowerType(options.spacecraft_power_model_type);
                myPowerSystemOptions.setPowerSupplyType(options.power_source_type);
                myPowerSystemOptions.setPowerSupplyCurveType(options.solar_power_model_type);
                myPowerSystemOptions.setDecayRate(options.power_decay_rate);
                myPowerSystemOptions.setPowerMargin(options.power_margin);
                myPowerSystemOptions.setPowerSystemDecayRefEpoch(options.power_system_decay_reference_epoch);
                myPowerSystemOptions.setGammaVector(options.solar_power_gamma);
                myPowerSystemOptions.setBusCoefficientVector(options.spacecraft_power_coefficients);

                //then the chemical propulsion system
                PropulsionSystemOptions myChemicalPropulsionSystemOptions;
                myChemicalPropulsionSystemOptions.setConstantThrust(options.Thrust);
                myChemicalPropulsionSystemOptions.setName("ChemicalPropulsionFromMissionOptions");
                myChemicalPropulsionSystemOptions.setConstantIsp(options.IspChem);
                myChemicalPropulsionSystemOptions.setMinimumOrMonopropIsp(options.TCM_Isp);
                myChemicalPropulsionSystemOptions.setMixtureRatio(options.bipropellant_mixture_ratio);
                myChemicalPropulsionSystemOptions.setg0(options.g0);

                //then the electric propulsion system
                PropulsionSystemOptions myElectricPropulsionSystemOptions;
                myElectricPropulsionSystemOptions.setName("ElectricPropulsionSystemFromMissionOptions");
                myElectricPropulsionSystemOptions.setNumberOfStrings(options.number_of_electric_propulsion_systems);
                myElectricPropulsionSystemOptions.setThrustScaleFactor(options.thrust_scale_factor);
                myElectricPropulsionSystemOptions.setConstantThrust(options.Thrust);
                myElectricPropulsionSystemOptions.setConstantIsp(options.IspLT);
                myElectricPropulsionSystemOptions.setMinimumOrMonopropIsp(options.IspLT_minimum);
                myElectricPropulsionSystemOptions.setFixedEfficiency(options.user_defined_engine_efficiency);
                myElectricPropulsionSystemOptions.setThrottleTableFile(options.HardwarePath + "/" + options.ThrottleTableFile);
                myElectricPropulsionSystemOptions.setMassPerString(0.0);
                myElectricPropulsionSystemOptions.setg0(options.g0);
                myElectricPropulsionSystemOptions.setPmin(options.engine_input_power_bounds[0]);
                myElectricPropulsionSystemOptions.setPmax(options.engine_input_power_bounds[1]);
				myElectricPropulsionSystemOptions.setSharpness(options.throttle_sharpness);
                
                switch (options.engine_type)
                {
                    case 0:
                        myElectricPropulsionSystemOptions.setThrusterMode(EMTG::SpacecraftThrusterMode::ConstantThrustIsp);
                        break;
                    case 1:
                        throw std::invalid_argument("User has selected engine_type 1: \"constant Isp, efficiency, EMTG chooses input power\". This engine type is not currently implemented. Halting program.");
                        break;
                    case 2:
                        myElectricPropulsionSystemOptions.setThrusterMode(EMTG::SpacecraftThrusterMode::FixedEfficiencyVSI);
                        throw std::invalid_argument("User has selected engine_type 2: \"choice of power model, constant efficiency, EMTG chooses Isp\". This engine type is not currently implemented. Halting program.");
                        break;
                    case 3:
                        myElectricPropulsionSystemOptions.setThrusterMode(EMTG::SpacecraftThrusterMode::FixedEfficiencyCSI);
                        break;
                    case 4:
                        myElectricPropulsionSystemOptions.setThrusterMode(EMTG::SpacecraftThrusterMode::FixedEfficiencyVSI);
                        throw std::invalid_argument("User has selected engine_type 4: \"continuously-varying specific impulse\". This engine type is not currently implemented in the phase transcriptions, although the spacecraft model does have it. Halting program.");
                        break;
                    case 5:
                        myElectricPropulsionSystemOptions.setThrusterMode(EMTG::SpacecraftThrusterMode::Poly1D);
                        myElectricPropulsionSystemOptions.setThrustCoefficients(options.engine_input_thrust_coefficients);
                        myElectricPropulsionSystemOptions.setMassFlowCoefficients(options.engine_input_mass_flow_rate_coefficients);

                        break;
                    case 29:
                        myElectricPropulsionSystemOptions.setThrusterMode(EMTG::SpacecraftThrusterMode::Stepped2D);
                        throw std::invalid_argument("User has selected engine_type 29: \"2D stepped model\". This engine type is not currently implemented in the phase transcriptions, although the spacecraft model does have it. Halting program.");
                        break;
                    case 30:
                        myElectricPropulsionSystemOptions.setThrusterMode(EMTG::SpacecraftThrusterMode::SteppedHThrust1D);
                        break;
                    case 31:
                        myElectricPropulsionSystemOptions.setThrusterMode(EMTG::SpacecraftThrusterMode::SteppedLMdot1D);
                        break;
                    case 32:
                        myElectricPropulsionSystemOptions.setThrusterMode(EMTG::SpacecraftThrusterMode::Poly2D);
                        throw std::invalid_argument("User has selected engine_type 32: \"2D polynomial\". This engine type is not currently implemented in the phase transcriptions, although the spacecraft model does have it. Halting program.");
                        break;
                    default: //built-in thruster model

#ifdef HAS_BUILT_IN_THRUSTERS
                        //recover coefficients from hard-coded models, but remember that the coefficients were in the reverse order back then
                        std::vector<double> temp_thrust_coefficients(7, 0.0);
                        std::vector<double> temp_mass_flow_coefficients(7, 0.0);
                        double temp_minP, temp_maxP;

                        myElectricPropulsionSystemOptions.setThrusterMode(EMTG::SpacecraftThrusterMode::Poly1D);

                        EMTG::Astrodynamics::get_thruster_coefficients_from_library(options,
                            temp_minP,
                            temp_maxP,
                            temp_thrust_coefficients[0], temp_thrust_coefficients[1], temp_thrust_coefficients[2], temp_thrust_coefficients[3], temp_thrust_coefficients[4], temp_thrust_coefficients[5], temp_thrust_coefficients[6],
                            temp_mass_flow_coefficients[0], temp_mass_flow_coefficients[1], temp_mass_flow_coefficients[2], temp_mass_flow_coefficients[3], temp_mass_flow_coefficients[4], temp_mass_flow_coefficients[5], temp_mass_flow_coefficients[6]);

                        myElectricPropulsionSystemOptions.setThrustCoefficients(temp_thrust_coefficients);
                        myElectricPropulsionSystemOptions.setMassFlowCoefficients(temp_mass_flow_coefficients);

                        myElectricPropulsionSystemOptions.setPmin(temp_minP);
                        myElectricPropulsionSystemOptions.setPmax(temp_maxP);
#else //HAS_BUILT_IN_THRUSTERS
                        throw std::invalid_argument("User has selected engine_type " + std::to_string(options.engine_type) + ". This thruster model is not included with your version of EMTG. You will have to construct a propulsion model file for it based on publicly available data. Halting program.");
#endif //HAS_BUILT_IN_THRUSTERS

                        break;
                }

                //build a stage
                EMTG::HardwareModels::StageOptions myStageOptions(myChemicalPropulsionSystemOptions, myElectricPropulsionSystemOptions, myPowerSystemOptions);

                myStageOptions.setThrottleLogic(options.throttle_logic_mode);
                myStageOptions.setmyChemicalPropulsionSystemOptionsName("ChemicalPropulsionFromMissionOptions");
                myStageOptions.setmyElectricPropulsionSystemOptionsName("ElectricPropulsionFromMissionOptions");
                myStageOptions.setmyPowerSystemOptionsName("PowerFromMissionOptions");

                //and attach it to the spacecraft!
                mySpacecraftOptions = EMTG::HardwareModels::SpacecraftOptions(std::vector<EMTG::HardwareModels::StageOptions>(1, myStageOptions));
                mySpacecraftOptions.setName(options.mission_name);

                //tanks
                if (options.enable_electric_propellant_tank_constraint)
                {
                    mySpacecraftOptions.setEnableGlobalElectricPropellantTankConstraint(true);
                    mySpacecraftOptions.setGlobalElectricPropellantTankCapacity(options.maximum_electric_propellant);
                }
                if (options.enable_chemical_propellant_tank_constraint)
                {
                    mySpacecraftOptions.setEnableGlobalChemicalPropellantTankConstraint(true);
                    mySpacecraftOptions.setGlobalFuelTankCapacity(options.maximum_chemical_fuel);
                    mySpacecraftOptions.setGlobalOxidizerTankCapacity(options.maximum_chemical_oxidizer);
                }
                if (options.constrain_dry_mass)
                {
                    mySpacecraftOptions.setEnableGlobalDryMassConstraint(true);
                    mySpacecraftOptions.setGlobalDryMassBounds(std::vector<double>({ options.final_mass_constraint_bounds[0], options.final_mass_constraint_bounds[1] }));
                }
            }

            return mySpacecraftOptions;
        }//end CreateSpacecraftOptions()
    }//end namespace HardwareModels
}//end namespace EMTG