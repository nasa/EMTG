#define BOOST_PYTHON_STATIC_LIB

#include "boost/python/module.hpp"
#include "boost/python/class.hpp"
#include "boost/python/dict.hpp"
#include "boost/python/import.hpp"
#include "boost/python/manage_new_object.hpp"

#include "LaunchVehicleOptions.h"
#include "LaunchVehicle.h"

#include "PowerSystemOptions.h"
#include "PowerSystem.h"

#include "PropulsionSystemOptions.h"
#include "ElectricPropulsionSystem.h"
#include "ChemicalPropulsionSystem.h"

#include "StageOptions.h"
#include "Stage.h"

#include "SpacecraftOptions.h"
#include "Spacecraft.h"

namespace EMTG
{
     namespace HardwareModels
     {
	    template<class T>
	    inline PyObject * managingPyObject(T *p)
	    {
	        return typename boost::python::manage_new_object::apply<T *>::type()(p);
	    }

	    template<class Copyable>
	    boost::python::object
	    generic__deepcopy__(boost::python::object copyable, boost::python::dict memo)
	    {
	        boost::python::object copyMod = boost::python::import("copy");
	        boost::python::object deepcopy = copyMod.attr("deepcopy");

	        Copyable *newCopyable(new Copyable(boost::python::extract<const Copyable&>(copyable)));
	        boost::python::object result(boost::python::detail::new_reference(managingPyObject(newCopyable)));

	        // HACK: copyableId shall be the same as the result of id(copyable) in Python -
	        // please tell me that there is a better way! (and which ;-p)
	        // int copyableId = (int)(copyable.ptr());
	        // memo[copyableId] = result;

	        boost::python::extract<boost::python::dict>(result.attr("__dict__"))().update(deepcopy(boost::python::extract<boost::python::dict>(copyable.attr("__dict__"))(),memo));

	        return result;
	    }

        BOOST_PYTHON_MODULE(PyHardware)
        {
			boost::python::class_<SpacecraftOptions>("SpacecraftOptions",boost::python::init<const std::string&>())
                .def("__deepcopy__",&generic__deepcopy__<SpacecraftOptions>)
                .def("parse_input_file",&SpacecraftOptions::parse_input_file)
                .def("write_output_file",&SpacecraftOptions::write_output_file)
                .def("remove_stage",&SpacecraftOptions::remove_stage)
                .def("setName",&SpacecraftOptions::setName)
                .def("setGlobalElectricPropellantTankCapacity",&SpacecraftOptions::setGlobalElectricPropellantTankCapacity)
                .def("setGlobalFuelTankCapacity",&SpacecraftOptions::setGlobalFuelTankCapacity)
                .def("setGlobalOxidizerTankCapacity",&SpacecraftOptions::setGlobalOxidizerTankCapacity)
                .def("setGlobalDryMassBounds",&SpacecraftOptions::setGlobalDryMassBounds)
                .def("setEnableGlobalElectricPropellantTankConstraint",&SpacecraftOptions::setEnableGlobalElectricPropellantTankConstraint)
                .def("setEnableGlobalChemicalPropellantTankConstraint",&SpacecraftOptions::setEnableGlobalChemicalPropellantTankConstraint)
                .def("setEnableGlobalDryMassConstraint",&SpacecraftOptions::setEnableGlobalDryMassConstraint)
                .def("setStageOptions",&SpacecraftOptions::setStageOptions)
                .def("getName",&SpacecraftOptions::getName)
                .def("getNumberOfStages",&SpacecraftOptions::getNumberOfStages)
                .def("getGlobalElectricPropellantTankCapacity",&SpacecraftOptions::getGlobalElectricPropellantTankCapacity)
                .def("getGlobalChemicalFuelTankCapacity",&SpacecraftOptions::getGlobalChemicalFuelTankCapacity)
                .def("getGlobalChemicalOxidizerTankCapacity",&SpacecraftOptions::getGlobalChemicalOxidizerTankCapacity)
                .def("getGlobalDryMassBounds",&SpacecraftOptions::getGlobalDryMassBounds)
                .def("getEnableGlobalElectricPropellantTankConstraint",&SpacecraftOptions::getEnableGlobalElectricPropellantTankConstraint)
                .def("getEnableGlobalChemicalPropellantTankConstraint",&SpacecraftOptions::getEnableGlobalChemicalPropellantTankConstraint)
                .def("getEnableGlobalDryMassConstraint",&SpacecraftOptions::getEnableGlobalDryMassConstraint)
                .def("getStageOptions",&SpacecraftOptions::getStageOptions)
				;

			boost::python::class_<StageOptions>("StageOptions",boost::python::init<const std::string&>())
                .def("__deepcopy__",&generic__deepcopy__<StageOptions>)
                .def("setName",&StageOptions::setName)
                .def("setBaseDryMass",&StageOptions::setBaseDryMass)
                .def("setAdapterMass",&StageOptions::setAdapterMass)
                .def("setElectricPropellantTankCapacity",&StageOptions::setElectricPropellantTankCapacity)
                .def("setChemicalFuelTankCapacity",&StageOptions::setChemicalFuelTankCapacity)
                .def("setChemicalOxidizerTankCapacity",&StageOptions::setChemicalOxidizerTankCapacity)
                .def("setEnableElectricPropellantTankConstraint",&StageOptions::setEnableElectricPropellantTankConstraint)
                .def("setEnableChemicalPropellantTankConstraint",&StageOptions::setEnableChemicalPropellantTankConstraint)
                .def("setEnableDryMassConstraint",&StageOptions::setEnableDryMassConstraint)
                .def("setmyPowerSystemOptionsName",&StageOptions::setmyPowerSystemOptionsName)
                .def("setmyElectricPropulsionSystemOptionsName",&StageOptions::setmyElectricPropulsionSystemOptionsName)
                .def("setmyChemicalPropulsionSystemOptionsName",&StageOptions::setmyChemicalPropulsionSystemOptionsName)
                .def("setThrottleLogic",&StageOptions::setThrottleLogic)
                .def("setThrottleSharpness",&StageOptions::setThrottleSharpness)
                .def("setElectricPropulsionSystemOptions",&StageOptions::setElectricPropulsionSystemOptions)
                .def("getName",&StageOptions::getName)
                .def("getBaseDryMass",&StageOptions::getBaseDryMass)
                .def("getAdapterMass",&StageOptions::getAdapterMass)
                .def("getElectricPropellantTankCapacity",&StageOptions::getElectricPropellantTankCapacity)
                .def("getChemicalFuelTankCapacity",&StageOptions::getChemicalFuelTankCapacity)
                .def("getChemicalOxidizerTankCapacity",&StageOptions::getChemicalOxidizerTankCapacity)
                .def("getEnableElectricPropellantTankConstraint",&StageOptions::getEnableElectricPropellantTankConstraint)
                .def("getEnableChemicalPropellantTankConstraint",&StageOptions::getEnableChemicalPropellantTankConstraint)
                .def("getEnableDryMassConstraint",&StageOptions::getEnableDryMassConstraint)
                .def("getmyChemicalPropulsionSystemOptionsName",&StageOptions::getmyChemicalPropulsionSystemOptionsName)
                .def("getmyElectricPropulsionSystemOptionsName",&StageOptions::getmyElectricPropulsionSystemOptionsName)
                .def("getmyPowerSystemOptionsName",&StageOptions::getmyPowerSystemOptionsName)
                .def("getElectricPropulsionSystemOptions",&StageOptions::getElectricPropulsionSystemOptions)
                .def("getChemicalPropulsionSystemOptions",&StageOptions::getChemicalPropulsionSystemOptions)
                .def("getPowerSystemOptions",&StageOptions::getPowerSystemOptions)
                .def("getThrottleLogic",&StageOptions::getThrottleLogic)
                .def("getThrottleSharpness",&StageOptions::getThrottleSharpness)
                .def("write_output_file",&StageOptions::write_output_file)
                .def("parse_input_file",static_cast<void (StageOptions::*)(const std::string& filename)>(&StageOptions::parse_input_file))
                .def("parse_input_file",static_cast<void (StageOptions::*)(std::ifstream& inputfile)>(&StageOptions::parse_input_file))
				;

			boost::python::class_<LaunchVehicleOptions>("LaunchVehicleOptions",boost::python::init<const std::string&>())
                .def("__deepcopy__",&generic__deepcopy__<LaunchVehicleOptions>)
                .def("parse_input_line",&LaunchVehicleOptions::parse_input_line)
                .def("write_output_line",&LaunchVehicleOptions::write_output_line)
                .def("getName",&LaunchVehicleOptions::getName)
                .def("getModelType",&LaunchVehicleOptions::getModelType)
                .def("getDLA_lowerbound",&LaunchVehicleOptions::getDLA_lowerbound)
                .def("getDLA_upperbound",&LaunchVehicleOptions::getDLA_upperbound)
                .def("getC3_lowerbound",&LaunchVehicleOptions::getC3_lowerbound)
                .def("getC3_upperbound",&LaunchVehicleOptions::getC3_upperbound)
                .def("getAdapterMass",&LaunchVehicleOptions::getAdapterMass)
                .def("getCoefficients",&LaunchVehicleOptions::getCoefficients)
                .def("getNumberOfCoefficients",&LaunchVehicleOptions::getNumberOfCoefficients)
                .def("getCoefficient",&LaunchVehicleOptions::getCoefficient)
                .def("setName",&LaunchVehicleOptions::setName)
                .def("setModelType",&LaunchVehicleOptions::setModelType)
                .def("setDLA_lowerbound",&LaunchVehicleOptions::setDLA_lowerbound)
                .def("setDLA_upperbound",&LaunchVehicleOptions::setDLA_upperbound)
                .def("setC3_lowerbound",&LaunchVehicleOptions::setC3_lowerbound)
                .def("setC3_upperbound",&LaunchVehicleOptions::setC3_upperbound)
                .def("setCoefficients",&LaunchVehicleOptions::setCoefficients)
                .def("setCoefficient",&LaunchVehicleOptions::setCoefficient)
                .def("setAdapterMass",&LaunchVehicleOptions::setAdapterMass)
				;

			boost::python::class_<LaunchVehicleOptionsLibrary>("LaunchVehicleOptionsLibrary",boost::python::init<const std::string&>())
                .def("__deepcopy__",&generic__deepcopy__<LaunchVehicleOptionsLibrary>)
                .def("write_output_file",&LaunchVehicleOptionsLibrary::write_output_file)
                .def("getLaunchVehicle",&LaunchVehicleOptionsLibrary::getLaunchVehicle)
                .def("getModelType",&LaunchVehicleOptionsLibrary::getModelType)
                .def("getDLA_lowerbound",&LaunchVehicleOptionsLibrary::getDLA_lowerbound)
                .def("getDLA_upperbound",&LaunchVehicleOptionsLibrary::getDLA_upperbound)
                .def("getC3_lowerbound",&LaunchVehicleOptionsLibrary::getC3_lowerbound)
                .def("getC3_upperbound",&LaunchVehicleOptionsLibrary::getC3_upperbound)
                .def("getAdapterMass",&LaunchVehicleOptionsLibrary::getAdapterMass)
                .def("getCoefficients",&LaunchVehicleOptionsLibrary::getCoefficients)
                .def("parse_input_file",static_cast<void (LaunchVehicleOptionsLibrary::*)(const std::string& filename)>(&LaunchVehicleOptionsLibrary::parse_input_file))
                .def("parse_input_file",static_cast<void (LaunchVehicleOptionsLibrary::*)(std::ifstream& inputfile)>(&LaunchVehicleOptionsLibrary::parse_input_file))
				;

			boost::python::class_<PowerSystemOptions>("PowerSystemOptions",boost::python::init<const std::string&>())
                .def("__deepcopy__",&generic__deepcopy__<PowerSystemOptions>)
                .def("parse_input_line",&PowerSystemOptions::parse_input_line)
                .def("write_output_line",&PowerSystemOptions::write_output_line)
                .def("getName",&PowerSystemOptions::getName)
                .def("getBusPowerType",&PowerSystemOptions::getBusPowerType)
                .def("getPowerSupplyType",&PowerSystemOptions::getPowerSupplyType)
                .def("getPowerSupplyCurveType",&PowerSystemOptions::getPowerSupplyCurveType)
                .def("getDecayRate",&PowerSystemOptions::getDecayRate)
                .def("getP0",&PowerSystemOptions::getP0)
                .def("getPowerMargin",&PowerSystemOptions::getPowerMargin)
                .def("getGammaVector",&PowerSystemOptions::getGammaVector)
                .def("getBusCoefficientVector",&PowerSystemOptions::getBusCoefficientVector)
                .def("getmass_per_kw",&PowerSystemOptions::getmass_per_kw)
                .def("getGamma",&PowerSystemOptions::getGamma)
                .def("getBusCoefficient",&PowerSystemOptions::getBusCoefficient)
                .def("setName",&PowerSystemOptions::setName)
                .def("setP0",&PowerSystemOptions::setP0)
                .def("setBusPowerType",&PowerSystemOptions::setBusPowerType)
                .def("setPowerSupplyType",&PowerSystemOptions::setPowerSupplyType)
                .def("setPowerSupplyCurveType",&PowerSystemOptions::setPowerSupplyCurveType)
                .def("setDecayRate",&PowerSystemOptions::setDecayRate)
                .def("setPowerMargin",&PowerSystemOptions::setPowerMargin)
                .def("setGammaVector",&PowerSystemOptions::setGammaVector)
                .def("setBusCoefficientVector",&PowerSystemOptions::setBusCoefficientVector)
                .def("setmass_per_kw",&PowerSystemOptions::setmass_per_kw)
                .def("setGamma",&PowerSystemOptions::setGamma)
                .def("setBusCoefficient",&PowerSystemOptions::setBusCoefficient)
				;

			boost::python::class_<PowerSystemOptionsLibrary>("PowerSystemOptionsLibrary",boost::python::init<const std::string&>())
                .def("__deepcopy__",&generic__deepcopy__<PowerSystemOptionsLibrary>)
                .def("write_output_file",&PowerSystemOptionsLibrary::write_output_file)
                .def("getPowerSystem",&PowerSystemOptionsLibrary::getPowerSystem)
                .def("getBusPowerType",&PowerSystemOptionsLibrary::getBusPowerType)
                .def("getPowerSupplyType",&PowerSystemOptionsLibrary::getPowerSupplyType)
                .def("getPowerSupplyCurveType",&PowerSystemOptionsLibrary::getPowerSupplyCurveType)
                .def("getDecayRate",&PowerSystemOptionsLibrary::getDecayRate)
                .def("getP0",&PowerSystemOptionsLibrary::getP0)
                .def("getPowerMargin",&PowerSystemOptionsLibrary::getPowerMargin)
                .def("getGamma",&PowerSystemOptionsLibrary::getGamma)
                .def("getBusCoefficients",&PowerSystemOptionsLibrary::getBusCoefficients)
                .def("getmass_per_kw",&PowerSystemOptionsLibrary::getmass_per_kw)
                .def("parse_input_file",static_cast<void (PowerSystemOptionsLibrary::*)(const std::string& filename)>(&PowerSystemOptionsLibrary::parse_input_file))
                .def("parse_input_file",static_cast<void (PowerSystemOptionsLibrary::*)(std::ifstream& inputfile)>(&PowerSystemOptionsLibrary::parse_input_file))
				;

			boost::python::class_<PropulsionSystemOptions>("PropulsionSystemOptions",boost::python::init<const std::string&>())
                .def("__deepcopy__",&generic__deepcopy__<PropulsionSystemOptions>)
                .def("parse_input_line",&PropulsionSystemOptions::parse_input_line)
                .def("write_output_line",&PropulsionSystemOptions::write_output_line)
                .def("getName",&PropulsionSystemOptions::getName)
                .def("getThrusterMode",&PropulsionSystemOptions::getThrusterMode)
                .def("getConstantThrust",&PropulsionSystemOptions::getConstantThrust)
                .def("getConstantIsp",&PropulsionSystemOptions::getConstantIsp)
                .def("getMinimumOrMonopropIsp",&PropulsionSystemOptions::getMinimumOrMonopropIsp)
                .def("getFixedEfficiency",&PropulsionSystemOptions::getFixedEfficiency)
                .def("getPmin",&PropulsionSystemOptions::getPmin)
                .def("getPmax",&PropulsionSystemOptions::getPmax)
                .def("getThrottleTableFile",&PropulsionSystemOptions::getThrottleTableFile)
                .def("getThrustCoefficients",&PropulsionSystemOptions::getThrustCoefficients)
                .def("getMassFlowCoefficients",&PropulsionSystemOptions::getMassFlowCoefficients)
                .def("getMassPerString",&PropulsionSystemOptions::getMassPerString)
                .def("getNumberOfStrings",&PropulsionSystemOptions::getNumberOfStrings)
                .def("getMixtureRatio",&PropulsionSystemOptions::getMixtureRatio)
                .def("getg0",&PropulsionSystemOptions::getg0)
                .def("getThrustScaleFactor",&PropulsionSystemOptions::getThrustScaleFactor)
                .def("getThrustCoefficient",&PropulsionSystemOptions::getThrustCoefficient)
                .def("getMassFlowCoefficient",&PropulsionSystemOptions::getMassFlowCoefficient)
                .def("setName",&PropulsionSystemOptions::setName)
                .def("setNumberOfStrings",&PropulsionSystemOptions::setNumberOfStrings)
                .def("setThrustScaleFactor",&PropulsionSystemOptions::setThrustScaleFactor)
                .def("setThrusterMode",&PropulsionSystemOptions::setThrusterMode)
                .def("setConstantThrust",&PropulsionSystemOptions::setConstantThrust)
                .def("setConstantIsp",&PropulsionSystemOptions::setConstantIsp)
                .def("setMinimumOrMonopropIsp",&PropulsionSystemOptions::setMinimumOrMonopropIsp)
                .def("setFixedEfficiency",&PropulsionSystemOptions::setFixedEfficiency)
                .def("setPmin",&PropulsionSystemOptions::setPmin)
                .def("setPmax",&PropulsionSystemOptions::setPmax)
                .def("setThrottleTableFile",&PropulsionSystemOptions::setThrottleTableFile)
                .def("setThrustCoefficients",&PropulsionSystemOptions::setThrustCoefficients)
                .def("setMassFlowCoefficients",&PropulsionSystemOptions::setMassFlowCoefficients)
                .def("setMassPerString",&PropulsionSystemOptions::setMassPerString)
                .def("setMixtureRatio",&PropulsionSystemOptions::setMixtureRatio)
                .def("setg0",&PropulsionSystemOptions::setg0)
                .def("setThrustCoefficient",&PropulsionSystemOptions::setThrustCoefficient)
                .def("setMassFlowCoefficient",&PropulsionSystemOptions::setMassFlowCoefficient)
				;

			boost::python::class_<PropulsionSystemOptionsLibrary>("PropulsionSystemOptionsLibrary",boost::python::init<const std::string&>())
                .def("__deepcopy__",&generic__deepcopy__<PropulsionSystemOptionsLibrary>)
                .def("write_output_file",&PropulsionSystemOptionsLibrary::write_output_file)
                .def("getPropulsionSystem",&PropulsionSystemOptionsLibrary::getPropulsionSystem)
                .def("getThrusterMode",&PropulsionSystemOptionsLibrary::getThrusterMode)
                .def("getConstantThrust",&PropulsionSystemOptionsLibrary::getConstantThrust)
                .def("getConstantIsp",&PropulsionSystemOptionsLibrary::getConstantIsp)
                .def("getMinimumOrMonopropIsp",&PropulsionSystemOptionsLibrary::getMinimumOrMonopropIsp)
                .def("getFixedEfficiency",&PropulsionSystemOptionsLibrary::getFixedEfficiency)
                .def("getPmin",&PropulsionSystemOptionsLibrary::getPmin)
                .def("getPmax",&PropulsionSystemOptionsLibrary::getPmax)
                .def("getThrottleTableFile",&PropulsionSystemOptionsLibrary::getThrottleTableFile)
                .def("getThrustCoefficients",&PropulsionSystemOptionsLibrary::getThrustCoefficients)
                .def("getMassFlowCoefficients",&PropulsionSystemOptionsLibrary::getMassFlowCoefficients)
                .def("getMassPerString",&PropulsionSystemOptionsLibrary::getMassPerString)
                .def("getNumberOfStrings",&PropulsionSystemOptionsLibrary::getNumberOfStrings)
                .def("getMixtureRatio",&PropulsionSystemOptionsLibrary::getMixtureRatio)
                .def("getThrustScaleFactor",&PropulsionSystemOptionsLibrary::getThrustScaleFactor)
                .def("parse_input_file",static_cast<void (PropulsionSystemOptionsLibrary::*)(const std::string& filename)>(&PropulsionSystemOptionsLibrary::parse_input_file))
                .def("parse_input_file",static_cast<void (PropulsionSystemOptionsLibrary::*)(std::ifstream& inputfile)>(&PropulsionSystemOptionsLibrary::parse_input_file))
				;

        }
    }
}
