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

//EMTGv9 mission class
//Jacob Englander 6-22-2017

#pragma once

#include "doubleType.h"

#include <vector>

#include "problem.h"
#include "journey.h"
#include "Spacecraft.h"
#include "LaunchVehicle.h"
#include "universe.h"
#include "missionoptions.h"
#include "ObjectiveFunctionBase.h"

#include "EMTG_solver_utilities.h"
#include "boost/ptr_container/ptr_vector.hpp"

namespace EMTG 
{
    //forward declaration of journey, allows "pointer to child" to work
    class Journey;

    //forward declaration of ObjectiveFunctionBase (pointer to child)
    namespace ObjectiveFunctions
    {
        class ObjectiveFunctionBase;
    }

    class Mission: public EMTG::problem
    {
    public:
        //constructor
        Mission();

        Mission(const missionoptions& options_in, 
            const std::vector<Astrodynamics::universe >& TheUniverse_in,
            const HardwareModels::LaunchVehicle& LaunchVehicle_in,
            const HardwareModels::Spacecraft& Spacecraft_in);

        //destructor
        ~Mission();

        //methods
        //evaluate function
        //return 0 if successful, 1 if failure
        virtual void evaluate(const std::vector<doubleType>& X,
                              std::vector<doubleType>& F,
                              std::vector<double>& G,
                              const bool& needG);

        //output function
        void output(const std::string& outputfilename);

        //function to write SPICE files
        void output_ephemeris();

        //function to write STMs
        void output_STMs();

        //method to write maneuver and target spec files
        void output_maneuver_and_target_spec();

        //bounds calculation function
        void calcbounds();

        //function to find constraints that are considered critical to the filament walker
        void locate_filament_critical_inequality_constraints();

        //function to construct an initial guess from inputs
        virtual std::vector<double> construct_initial_guess();

        //get
        Journey* getJourney(const size_t& journeyIndex) { return &this->Journeys[journeyIndex]; }
		
		doubleType getUnscaledObjective();

    private:    
        //utilities
        size_t create_sparsity_entry(const size_t& Findex,
            const size_t& Xstart,
            const bool& ForwardPass,
            const std::string& variable_name,
            size_t& sparsity_index_container);

        size_t create_sparsity_entry(const size_t& Findex,
            const size_t& Xstart,
            const bool& ForwardPass,
            const std::string& variable_name,
            std::vector<size_t>& sparsity_index_container);

        size_t create_sparsity_entry(const size_t& Findex,
            const size_t& Xindex,
            size_t& sparsity_index_container);

        size_t create_sparsity_entry(const size_t& Findex,
            const size_t& Xindex,
            std::vector<size_t>&  sparsity_index_container);

        void create_sparsity_entry_vector(const size_t& Findex,
            const size_t& Xstart,
            const bool& ForwardPass,
            const int& number_of_entries,
            const std::string& variable_name,
            std::vector<size_t>& sparsity_index_container);

        //vector of journeys
        boost::ptr_vector< Journey > Journeys;

        //vector of universes
        std::vector<Astrodynamics::universe > TheUniverse;

        //fields
        size_t number_of_journeys;
        size_t total_number_of_phases;
        doubleType total_deterministic_deltav;
        doubleType total_statistical_deltav;

        //spacecraft object
        HardwareModels::Spacecraft mySpacecraft;

        //launch vehicle object
        HardwareModels::LaunchVehicle myLaunchVehicle;

        //tank constraints
        std::vector<size_t> Gindex_dElectricPropellantConstraint_dVirtualElectricPropellant;
        std::vector<size_t> Gindex_dChemicalFuelConstraint_dVirtualChemicalFuel;
        std::vector<size_t> Gindex_dChemicalOxidizerConstraint_dVirtualChemicalOxidizer;
        std::vector< std::vector<size_t> > Gindex_PerStage_dElectricPropellantConstraint_dVirtualElectricPropellant;
        std::vector< std::vector<size_t> > Gindex_PerStage_dChemicalFuelConstraint_dVirtualChemicalFuel;
        std::vector< std::vector<size_t> > Gindex_PerStage_dChemicalOxidizerConstraint_dVirtualChemicalOxidizer;

        //final mass constraint
        std::vector<size_t> Gindex_derivatives_of_final_mass_constraint;
        std::vector<size_t> dIndex_final_mass_constraint;
        std::vector<size_t> dIndex_final_mass_constraint_wrt_Time;

        //global dry mass constraint
        std::vector<size_t> Gindex_derivatives_of_global_dry_mass_constraint;
        std::vector<size_t> dIndex_global_dry_mass_constraint;
        std::vector<size_t> dIndex_global_dry_mass_constraint_wrt_Time;
        std::vector<size_t> Gindex_dGlobalDryMassConstraint_dVirtualElectricPropellant;
        std::vector<size_t> Gindex_dGlobalDryMassConstraint_dVirtualChemicalFuel;
        std::vector<size_t> Gindex_dGlobalDryMassConstraint_dVirtualChemicalOxidizer;

        //time constraint
        doubleType TotalFlightTime;
        std::vector<size_t> timeconstraints_G_indices;
        std::vector<size_t> timeconstraints_X_indices;

        //pointer to objective function object
        bool ObjectiveSet;
        ObjectiveFunctions::ObjectiveFunctionBase* myObjectiveFunction;
    };

} /* namespace EMTG */
