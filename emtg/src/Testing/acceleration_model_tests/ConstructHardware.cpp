

#include "ConstructHardware.h"


void ConstructHardware(EMTG::missionoptions & options, EMTG::HardwareModels::LaunchVehicleOptions & myLaunchVehicleOptions, EMTG::HardwareModels::SpacecraftOptions & mySpacecraftOptions)
{

    //first configure the launch vehicle    
    std::vector<std::string> LV_names = { "Fixed_Initial_Mass",
                                         "Atlas_V_401",
                                         "Atlas_V_411",
                                         "Atlas_V_421",
                                         "Atlas_V_431",
                                         "Atlas_V_501",
                                         "Atlas_V_511",
                                         "Atlas_V_521",
                                         "Atlas_V_531",
                                         "Atlas_V_541",
                                         "Atlas_V_551",
                                         "F910",
                                         "F911",
                                         "Atlas_V_551_Star48",
                                         "Falcon_Heavy_(Expendable)",
                                         "Delta_IV_Heavy",
                                         "SLSb1" };
    if (options.LV_type == 0)
    {
        //fixed mass
        myLaunchVehicleOptions.setName("Fixed_Initial_Mass");
        myLaunchVehicleOptions.setCoefficients(std::vector<double>(1, options.maximum_mass));
        myLaunchVehicleOptions.setC3_lowerbound(options.Journeys.front().journey_initial_impulse_bounds[0] * options.Journeys.front().journey_initial_impulse_bounds[0]);
        myLaunchVehicleOptions.setC3_upperbound(options.Journeys.front().journey_initial_impulse_bounds[1] * options.Journeys.front().journey_initial_impulse_bounds[1]);
        myLaunchVehicleOptions.setDLA_lowerbound(options.DLA_bounds[0]);
        myLaunchVehicleOptions.setDLA_upperbound(options.DLA_bounds[1]);
    }
    else if (options.LV_type > 0)
    {
        //LV from library
        EMTG::HardwareModels::LaunchVehicleOptionsLibrary myLaunchVehicleOptionsLibrary(options.LaunchVehicleLibraryFile);
        myLaunchVehicleOptions = myLaunchVehicleOptionsLibrary.getLaunchVehicle(LV_names[options.LV_type]);
    }

    //then configure the power system
    EMTG::HardwareModels::PowerSystemOptions myPowerSystemOptions;
    myPowerSystemOptions.setName("PowerFromMissionOptions");
    myPowerSystemOptions.setP0(options.power_at_1_AU);
    myPowerSystemOptions.setBusPowerType(options.spacecraft_power_model_type);
    myPowerSystemOptions.setPowerSupplyType(options.power_source_type);
    myPowerSystemOptions.setPowerSupplyCurveType(options.solar_power_model_type);
    myPowerSystemOptions.setDecayRate(options.power_decay_rate);
    myPowerSystemOptions.setPowerMargin(options.power_margin);
    std::vector<double> temp(7);
    temp.assign(options.solar_power_gamma, options.solar_power_gamma + 7);
    myPowerSystemOptions.setGammaVector(temp);
    temp.clear();
    temp.assign(options.spacecraft_power_coefficients, options.spacecraft_power_coefficients + 3);
    myPowerSystemOptions.setBusCoefficientVector(temp);

    //then the chemical propulsion system
    EMTG::HardwareModels::PropulsionSystemOptions myChemicalPropulsionSystemOptions;
    myChemicalPropulsionSystemOptions.setName("ChemicalPropulsionFromMissionOptions");
    myChemicalPropulsionSystemOptions.setConstantIsp(options.IspChem);
    myChemicalPropulsionSystemOptions.setMinimumOrMonopropIsp(options.TCM_Isp);
    myChemicalPropulsionSystemOptions.setMixtureRatio(options.bipropellant_mixture_ratio);
    myChemicalPropulsionSystemOptions.setg0(options.g0);

    //then the electric propulsion system
    EMTG::HardwareModels::PropulsionSystemOptions myElectricPropulsionSystemOptions;
    myElectricPropulsionSystemOptions.setName("ElectricPropulsionSystemFromMissionOptions");
    myElectricPropulsionSystemOptions.setNumberOfStrings(options.number_of_electric_propulsion_systems);
    myElectricPropulsionSystemOptions.setThrustScaleFactor(options.thrust_scale_factor);
    myElectricPropulsionSystemOptions.setConstantThrust(options.Thrust);
    myElectricPropulsionSystemOptions.setConstantIsp(options.IspLT);
    myElectricPropulsionSystemOptions.setMinimumOrMonopropIsp(options.IspLT_minimum);
    myElectricPropulsionSystemOptions.setFixedEfficiency(options.user_defined_engine_efficiency);
    myElectricPropulsionSystemOptions.setThrottleTableFile(options.HardwarePath + options.ThrottleTableFile);
    myElectricPropulsionSystemOptions.setMassPerString(0.0);
    myElectricPropulsionSystemOptions.setg0(options.g0);
    myElectricPropulsionSystemOptions.setPmin(options.engine_input_power_bounds[0]);
    myElectricPropulsionSystemOptions.setPmax(options.engine_input_power_bounds[1]);

    //the following settings are set differently depending on which thruster model is used
    std::vector<double> temp_thrust_coefficients(7, 0.0);
    std::vector<double> temp_mass_flow_coefficients(7, 0.0);
    double temp_minP, temp_maxP;

    switch (options.engine_type)
    {
    case 0:
        myElectricPropulsionSystemOptions.setThrusterMode(EMTG::SpacecraftThrusterMode::ConstantThrustIsp);
        break;
    case 1:
        myElectricPropulsionSystemOptions.setThrusterMode(EMTG::SpacecraftThrusterMode::FixedEfficiencyCSI);
        break;
    case 2:
        myElectricPropulsionSystemOptions.setThrusterMode(EMTG::SpacecraftThrusterMode::FixedEfficiencyVSI);
        break;
    case 3:
        myElectricPropulsionSystemOptions.setThrusterMode(EMTG::SpacecraftThrusterMode::FixedEfficiencyCSI);
        break;
    case 4:
        myElectricPropulsionSystemOptions.setThrusterMode(EMTG::SpacecraftThrusterMode::FixedEfficiencyVSI);
        break;
    case 5:
        myElectricPropulsionSystemOptions.setThrusterMode(EMTG::SpacecraftThrusterMode::Poly1D);

        temp_thrust_coefficients.assign(options.engine_input_thrust_coefficients, options.engine_input_thrust_coefficients + 7);
        temp_mass_flow_coefficients.assign(options.engine_input_mass_flow_rate_coefficients, options.engine_input_mass_flow_rate_coefficients + 7);

        myElectricPropulsionSystemOptions.setThrustCoefficients(temp_thrust_coefficients);
        myElectricPropulsionSystemOptions.setMassFlowCoefficients(temp_mass_flow_coefficients);

        break;
    case 29:
        myElectricPropulsionSystemOptions.setThrusterMode(EMTG::SpacecraftThrusterMode::Stepped2D);
        break;
    case 30:
        myElectricPropulsionSystemOptions.setThrusterMode(EMTG::SpacecraftThrusterMode::SteppedHThrust1D);
        break;
    case 31:
        myElectricPropulsionSystemOptions.setThrusterMode(EMTG::SpacecraftThrusterMode::SteppedLMdot1D);
        break;
    case 32:
        myElectricPropulsionSystemOptions.setThrusterMode(EMTG::SpacecraftThrusterMode::Poly2D);
        break;
    default:
        //recover coefficients from hard-coded models, but remember that the coefficients were in the reverse order back then
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

        break;
    }

    //build a stage
    EMTG::HardwareModels::StageOptions myStageOptions(myChemicalPropulsionSystemOptions, myElectricPropulsionSystemOptions, myPowerSystemOptions);

    //and attach it to the spacecraft!
    mySpacecraftOptions = EMTG::HardwareModels::SpacecraftOptions(std::vector<EMTG::HardwareModels::StageOptions>(1, myStageOptions));
    mySpacecraftOptions.setName(options.mission_name);
}