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

/*
 * mission.cpp
 *
 *  Created on: Jul 17, 2012
 *      Author: Jacob
 */

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <random>

#include "mission.h"
#include "interpolator.h"

#include "ObjectiveFunctionFactory.h"

#include "boost/algorithm/string.hpp"
#include "boost/algorithm/string/predicate.hpp"
#include "boost/date_time.hpp"
#include "boost/date_time/local_time/local_date_time.hpp"

namespace EMTG 
{

    Mission::Mission()
		: problem::problem()
    {
		this->number_of_journeys = 0;
		this->ObjectiveSet = false;
        this->myObjectiveFunction = nullptr;
    }

    Mission::Mission(const missionoptions& options_in,
        const std::vector<Astrodynamics::universe >& TheUniverse_in,
        const HardwareModels::LaunchVehicle& LaunchVehicle_in,
        const HardwareModels::Spacecraft& Spacecraft_in)
		: Mission()
    {
        //first make a local copy of the options structure
        this->options = options_in;

        //make a local copy of the Universe vector
        this->TheUniverse = TheUniverse_in;

        //shorten the journeyoptions vector so that it matches the journeys we actually want
        this->options.number_of_journeys = std::min(this->options.number_of_journeys, this->options.stop_after_journey + 1);
        while (this->options.Journeys.size() > this->options.number_of_journeys)
            this->options.Journeys.pop_back();

        //configure the hardware objects
        this->myLaunchVehicle = LaunchVehicle_in;
        this->mySpacecraft = Spacecraft_in;
        this->mySpacecraft.resetStaging();

        //next, create the journeys
        this->total_number_of_phases = 0;
        this->number_of_journeys = this->options.number_of_journeys;
        size_t stageIndex = 0;
        try
        {
            for (size_t journeyIndex = 0; journeyIndex < this->number_of_journeys; ++journeyIndex)
            {
                Journey* PreviousJourney;
                if (journeyIndex == 0)
                    PreviousJourney = NULL;
                else
                    PreviousJourney = &this->Journeys.back();

                this->Journeys.push_back(new Journey(journeyIndex,
                    stageIndex,
                    &options,
                    &this->TheUniverse[journeyIndex],
                    PreviousJourney,
                    this,
                    &myLaunchVehicle,
                    &mySpacecraft));

                this->total_number_of_phases += this->Journeys.back().getNumberOfPhases();
            }
        }
        catch (std::exception &myError)
        {
            std::cout << "Failure in journey construction. Error message is:" << std::endl;
            std::cout << myError.what() << std::endl;
            throw;
        }

        //objective function object
        this->ObjectiveSet = true;
        this->myObjectiveFunction = ObjectiveFunctions::create_objective_function(&this->TheUniverse.back(),
            &this->mySpacecraft,
            &this->options,
            this);

        //set up calcbounds
        this->myObjectiveFunction->setup_calcbounds(&this->Xupperbounds,
                                                    &this->Xlowerbounds,
                                                    &this->X_scale_factors,
                                                    &this->Fupperbounds, 
                                                    &this->Flowerbounds,
													&this->F_scale_factors,
                                                    &this->Xdescriptions, 
                                                    &this->Fdescriptions, 
                                                    &this->iGfun,
                                                    &this->jGvar,
                                                    &this->Gdescriptions,
                                                    &this->iAfun,
                                                    &this->jAvar,
                                                    &this->Adescriptions,
                                                    &this->A);
        
        //calculate the upper and lower bounds on the decision variables and the constraints
        try
        {
            this->calcbounds();
        }
        catch (std::exception &myError)
        {
            std::cout << "Failure in mission::calcbounds(). Error message is:" << std::endl;
            std::cout << myError.what() << std::endl;
            throw;
        }


        //size the local "G storage" vector
        this->G.resize(this->Gdescriptions.size());

        //make sure that we have enough scale ranges
        assert(this->X_scale_factors.size() == this->Xupperbounds.size());

		// TODO: make these throws once we figure stuff out
		for (size_t k = 0; k < this->F.size(); ++k)
		{
			if (this->Fupperbounds[k] < this->Flowerbounds[k])
			{
				std::cout << std::endl << std::endl << "Upper and lower bounds are incompatible for constraint: " + this->Fdescriptions[k] + ", lower bound: " + std::to_string(this->Flowerbounds[k]) + ", upper bound: " + std::to_string(this->Fupperbounds[k]) << std::endl;
			}
		}

		for (size_t k = 0; k < this->X.size(); ++k)
		{
			if (this->Xupperbounds[k] < this->Xlowerbounds[k])
			{
				std::cout << std::endl << std::endl << "Upper and lower bounds are incompatible for decision variable: " + this->Xdescriptions[k] + ", lower bound: " + std::to_string(this->Xlowerbounds[k]) + ", upper bound: " + std::to_string(this->Xupperbounds[k]) << std::endl;
			}
		}

        //pass scale factors pointer to Universe and reset pointer to next universe
        for (size_t journeyIndex = 0; journeyIndex < this->number_of_journeys; ++journeyIndex)
        {
            TheUniverse[journeyIndex].set_X_scale_factors(&this->X_scale_factors);

            if (journeyIndex < this->number_of_journeys - 1)
                TheUniverse[journeyIndex].set_nextUniverse(TheUniverse[journeyIndex + 1]);
        }

    }//end constructor

    //destructor
    Mission::~Mission()
    {
        delete this->myObjectiveFunction;
    }//end destructor

	doubleType Mission::getUnscaledObjective()
	{
		return this->myObjectiveFunction->getUnscaledObjective();
	}

    //bounds calculation function
    //return 0 for success, 1 for failure
    void Mission::calcbounds()
    {
        //bounds on the objective function
        this->Flowerbounds.push_back(-math::LARGE);
        this->Fupperbounds.push_back(math::LARGE);
        this->Fdescriptions.push_back("objective function");

        //call the calcbounds() function for each journey
        for (size_t journeyIndex = 0; journeyIndex < number_of_journeys; ++journeyIndex)
        {
            this->Journeys[journeyIndex].calcbounds(Xupperbounds,
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
                A,
                synodic_periods);
        }

        //one final constraint, if applicable, for total time bounds. Plus we need to collect all the time variables so that we can write out the mission flight time
        std::vector<size_t> timeVariables = this->getJourney(this->options.number_of_journeys - 1)->getLastPhase()->getArrivalEvent()->get_Xindices_EventRightEpoch();
        for (size_t Xindex : timeVariables)
        {
            if (!(Xdescriptions[Xindex].find("epoch") < 1024))//we want to exclude launch epoch, which is the only thing called "epoch"
            {
                this->timeconstraints_X_indices.push_back(Xindex);
            }
        }

        if (this->options.global_timebounded)
        {
            this->Flowerbounds.push_back(options.total_flight_time_bounds[0] / options.total_flight_time_bounds[1]);
            this->Fupperbounds.push_back(1.0);
            this->Fdescriptions.push_back("Mission flight time bounds");

            //has derivatives with respect to all "time" variables in the mission. Note that we've already removed the launch epoch in the previous step
            for (size_t Xindex : this->timeconstraints_X_indices)
            {
                iGfun.push_back(Fdescriptions.size() - 1);
                jGvar.push_back(Xindex);
                Gdescriptions.push_back("Derivative of " + Fdescriptions.back()
                    + " F[" + std::to_string(iGfun.back()) + "] with respect to X["
                    + std::to_string(jGvar.back()) + "]: " + Xdescriptions[Xindex]);
                    
                this->timeconstraints_G_indices.push_back(Gdescriptions.size() - 1);
            }
        }//end total flight time constraint
        
        //objective function
        this->myObjectiveFunction->calcbounds();

        //propellant tank constraints
        //first for the spacecraft
        if (this->mySpacecraft.getSpacecraftOptions().getEnableGlobalElectricPropellantTankConstraint())
        {
            std::vector<size_t> Xindices_tank = this->mySpacecraft.getGlobalElectricPropellantTank_Xindices();
            size_t nX = Xindices_tank.size();

            if (nX > 0)
            {
                this->Flowerbounds.push_back(0.0);
                this->Fupperbounds.push_back(1.0 / (1.0 + this->options.electric_propellant_margin));
                this->Fdescriptions.push_back("global electric propellant tank constraint");

                this->mySpacecraft.setGlobalElectricPropellantTank_Findex(this->Fdescriptions.size() - 1);

                for (size_t virtualTankIndex = 0; virtualTankIndex < nX; ++virtualTankIndex)
                {
                    size_t Xindex = Xindices_tank[virtualTankIndex];
                    
                    this->create_sparsity_entry(this->Fdescriptions.size() - 1,
                        Xindex,
                        this->Gindex_dElectricPropellantConstraint_dVirtualElectricPropellant);
                }
            }
            else
                this->mySpacecraft.getSpacecraftOptions().setEnableGlobalElectricPropellantTankConstraint(false);
        }//end global electric propellant tank

        if (this->mySpacecraft.getSpacecraftOptions().getEnableGlobalChemicalPropellantTankConstraint())
        {
            std::vector<size_t> Xindices_tank = this->mySpacecraft.getGlobalChemicalFuelTank_Xindices();
            size_t nX = Xindices_tank.size();
            if (nX > 0)
            {
                //fuel
                this->Flowerbounds.push_back(0.0);
                this->Fupperbounds.push_back(1.0 / (1.0 + this->options.chemical_propellant_margin));
                this->Fdescriptions.push_back("global chemical fuel tank constraint");

                this->mySpacecraft.setGlobalChemicalFuelTank_Findex(this->Fdescriptions.size() - 1);

                for (size_t virtualTankIndex = 0; virtualTankIndex < nX; ++virtualTankIndex)
                {
                    size_t Xindex = Xindices_tank[virtualTankIndex];

                    this->create_sparsity_entry(this->Fdescriptions.size() - 1,
                        Xindex,
                        this->Gindex_dChemicalFuelConstraint_dVirtualChemicalFuel);
                }

                //oxidizer
                Xindices_tank = this->mySpacecraft.getGlobalChemicalOxidizerTank_Xindices();
                nX = Xindices_tank.size();

                if (nX > 0)
                {
                    this->Flowerbounds.push_back(0.0);
                    this->Fupperbounds.push_back(1.0 / (1.0 + this->options.chemical_propellant_margin));
                    this->Fdescriptions.push_back("global chemical oxidizer tank constraint");

                    for (size_t virtualTankIndex = 0; virtualTankIndex < nX; ++virtualTankIndex)
                    {
                        size_t Xindex = Xindices_tank[virtualTankIndex];

                        this->create_sparsity_entry(this->Fdescriptions.size() - 1,
                            Xindex,
                            this->Gindex_dChemicalOxidizerConstraint_dVirtualChemicalOxidizer);
                    }
                }
            }
            else
                this->mySpacecraft.getSpacecraftOptions().setEnableGlobalChemicalPropellantTankConstraint(false);
        }//end global chemical propellant tanks

        //global final mass constraint
        if (this->options.constrain_final_mass)
        {
            this->Flowerbounds.push_back(this->options.final_mass_constraint_bounds[0] / this->options.final_mass_constraint_bounds[1] - 1.0);
            this->Fupperbounds.push_back(0.0);
            this->Fdescriptions.push_back("mission final mass constraint");

            //this constraint has derivatives with respect to ALL variables that affect final arrival event mass
            size_t number_of_phases_in_final_journey = this->Journeys.back().getNumberOfPhases();

            BoundaryEvents::ArrivalEvent* myArrivalEvent = this->Journeys.back().getPhase(number_of_phases_in_final_journey - 1)->getArrivalEvent();
            std::vector< std::tuple<size_t, size_t, double> > Derivatives_of_StateAfterArrival = myArrivalEvent->get_Derivatives_of_StateAfterEvent();
            std::vector< std::tuple<size_t, size_t, double> > Derivatives_of_StateAfterArrival_wrt_Time = myArrivalEvent->get_Derivatives_of_StateAfterEvent_wrt_Time();

            //non-time variables
            for (size_t dIndex = 0; dIndex < Derivatives_of_StateAfterArrival.size(); ++dIndex)
            {
                if (std::get<1>(Derivatives_of_StateAfterArrival[dIndex]) == 6)
                {
                    this->dIndex_final_mass_constraint.push_back(dIndex);
                    this->create_sparsity_entry(this->Fdescriptions.size() - 1,
                        std::get<0>(Derivatives_of_StateAfterArrival[dIndex]),
                        this->Gindex_derivatives_of_final_mass_constraint);
                }
            }

            //time variables
            for (size_t dIndex = 0; dIndex < Derivatives_of_StateAfterArrival_wrt_Time.size(); ++dIndex)
            {
                if (std::get<1>(Derivatives_of_StateAfterArrival_wrt_Time[dIndex]) == 6)
                {
                    this->dIndex_final_mass_constraint_wrt_Time.push_back(dIndex);
                    this->create_sparsity_entry(this->Fdescriptions.size() - 1,
                        std::get<0>(Derivatives_of_StateAfterArrival_wrt_Time[dIndex]),
                        this->Gindex_derivatives_of_final_mass_constraint);
                }
            }
        }//end global final mass constraint

        //global dry mass constraint
        if (this->options.constrain_dry_mass)
        {
            this->Flowerbounds.push_back(this->options.final_mass_constraint_bounds[0] / this->options.final_mass_constraint_bounds[1] - 1.0);
            this->Fupperbounds.push_back(0.0);
            this->Fdescriptions.push_back("Global dry mass constraint");
            this->mySpacecraft.setGlobalDryMassConstraint_Findex(this->Fdescriptions.size() - 1);
            
            //this constraint has derivatives with respect to ALL variables that affect global_dry arrival event mass
            size_t number_of_phases_in_global_dry_journey = this->Journeys.back().getNumberOfPhases();

            BoundaryEvents::ArrivalEvent* myArrivalEvent = this->Journeys.back().getPhase(number_of_phases_in_global_dry_journey - 1)->getArrivalEvent();
            std::vector< std::tuple<size_t, size_t, double> > Derivatives_of_StateAfterArrival = myArrivalEvent->get_Derivatives_of_StateAfterEvent();
            std::vector< std::tuple<size_t, size_t, double> > Derivatives_of_StateAfterArrival_wrt_Time = myArrivalEvent->get_Derivatives_of_StateAfterEvent_wrt_Time();

            //non-time variables
            for (size_t dIndex = 0; dIndex < Derivatives_of_StateAfterArrival.size(); ++dIndex)
            {
                if (std::get<1>(Derivatives_of_StateAfterArrival[dIndex]) == 6)
                {
                    this->dIndex_global_dry_mass_constraint.push_back(dIndex);
                    this->create_sparsity_entry(this->Fdescriptions.size() - 1,
                        std::get<0>(Derivatives_of_StateAfterArrival[dIndex]),
                        this->Gindex_derivatives_of_global_dry_mass_constraint);
                }
            }

            //time variables
            for (size_t dIndex = 0; dIndex < Derivatives_of_StateAfterArrival_wrt_Time.size(); ++dIndex)
            {
                if (std::get<1>(Derivatives_of_StateAfterArrival_wrt_Time[dIndex]) == 6)
                {
                    this->dIndex_global_dry_mass_constraint_wrt_Time.push_back(dIndex);
                    this->create_sparsity_entry(this->Fdescriptions.size() - 1,
                        std::get<0>(Derivatives_of_StateAfterArrival_wrt_Time[dIndex]),
                        this->Gindex_derivatives_of_global_dry_mass_constraint);
                }
            }


            //has derivatives with respect to all virtual tank variables in the mission

            //electric tank
            std::vector<size_t> Xindices_tank = this->mySpacecraft.getGlobalElectricPropellantTank_Xindices();
            size_t nX = Xindices_tank.size();

            if (nX > 0)
            {
                for (size_t virtualTankIndex = 0; virtualTankIndex < nX; ++virtualTankIndex)
                {
                    size_t Xindex = Xindices_tank[virtualTankIndex];

                    this->create_sparsity_entry(this->Fdescriptions.size() - 1,
                        Xindex,
                        this->Gindex_dGlobalDryMassConstraint_dVirtualElectricPropellant);
                }
            }

            //chemical fuel tank
            Xindices_tank = this->mySpacecraft.getGlobalChemicalFuelTank_Xindices();
            nX = Xindices_tank.size();

            if (nX > 0)
            {
                for (size_t virtualTankIndex = 0; virtualTankIndex < nX; ++virtualTankIndex)
                {
                    size_t Xindex = Xindices_tank[virtualTankIndex];

                    this->create_sparsity_entry(this->Fdescriptions.size() - 1,
                        Xindex,
                        this->Gindex_dGlobalDryMassConstraint_dVirtualChemicalFuel);
                }
            }

            //chemical oxidizer tank
            Xindices_tank = this->mySpacecraft.getGlobalChemicalOxidizerTank_Xindices();
            nX = Xindices_tank.size();

            if (nX > 0)
            {
                for (size_t virtualTankIndex = 0; virtualTankIndex < nX; ++virtualTankIndex)
                {
                    size_t Xindex = Xindices_tank[virtualTankIndex];

                    this->create_sparsity_entry(this->Fdescriptions.size() - 1,
                        Xindex,
                        this->Gindex_dGlobalDryMassConstraint_dVirtualChemicalOxidizer);
                }
            }
        }//end global dry mass constraint

        //then for individual stages
        for (size_t stageIndex = 0; stageIndex < this->mySpacecraft.getSpacecraftOptions().getNumberOfStages(); ++stageIndex)
        {
            if (this->mySpacecraft.getStageOptions(stageIndex).getEnableElectricPropellantTankConstraint())
            {
                std::vector<size_t> Xindices_tank = this->mySpacecraft.getElectricPropellantTank_Xindices(stageIndex);
                size_t nX = Xindices_tank.size();
                
                if (nX > 0)
                {
                    this->Flowerbounds.push_back(0.0);
                    this->Fupperbounds.push_back(1.0 / (1.0 + this->options.electric_propellant_margin));
                    this->Fdescriptions.push_back("Stage " + std::to_string(stageIndex) + " electric propellant tank constraint");

                    this->mySpacecraft.setElectricPropellantTank_Findex(stageIndex, this->Fdescriptions.size() - 1);

                    std::vector<size_t> stage_Gindex_PerStage_dElectricPropellantConstraint_dVirtualElectricPropellant;
                    for (size_t virtualTankIndex = 0; virtualTankIndex < nX; ++virtualTankIndex)
                    {
                        size_t Xindex = Xindices_tank[virtualTankIndex];

                        this->create_sparsity_entry(this->Fdescriptions.size() - 1,
                            Xindex,
                            stage_Gindex_PerStage_dElectricPropellantConstraint_dVirtualElectricPropellant);
                    }
                    this->Gindex_PerStage_dElectricPropellantConstraint_dVirtualElectricPropellant.push_back(stage_Gindex_PerStage_dElectricPropellantConstraint_dVirtualElectricPropellant);
                }
                else
                    this->mySpacecraft.getStageOptions(stageIndex).setEnableElectricPropellantTankConstraint(false);
            }//end stage electric propellant tank

            if (this->mySpacecraft.getStageOptions(stageIndex).getEnableChemicalPropellantTankConstraint())
            {
                std::vector<size_t> Xindices_tank = this->mySpacecraft.getChemicalFuelTank_Xindices(stageIndex);
                size_t nX = Xindices_tank.size();

                if (nX > 0)
                {
                    //fuel
                    this->Flowerbounds.push_back(0.0);
                    this->Fupperbounds.push_back(1.0 / (1.0 + this->options.chemical_propellant_margin));
                    this->Fdescriptions.push_back("Stage " + std::to_string(stageIndex) + " chemical fuel tank constraint");

                    this->mySpacecraft.setChemicalFuelTank_Findex(stageIndex, this->Fdescriptions.size() - 1);

                    std::vector<size_t> stage_Gindex_PerStage_dChemicalFuelConstraint_dVirtualChemicalFuel;
                    for (size_t virtualTankIndex = 0; virtualTankIndex < nX; ++virtualTankIndex)
                    {
                        size_t Xindex = Xindices_tank[virtualTankIndex];

                        this->create_sparsity_entry(this->Fdescriptions.size() - 1,
                            Xindex,
                            stage_Gindex_PerStage_dChemicalFuelConstraint_dVirtualChemicalFuel);
                    }
                    this->Gindex_PerStage_dChemicalFuelConstraint_dVirtualChemicalFuel.push_back(stage_Gindex_PerStage_dChemicalFuelConstraint_dVirtualChemicalFuel);

                    //oxidizer
                    Xindices_tank = this->mySpacecraft.getChemicalOxidizerTank_Xindices(stageIndex);
                    nX = Xindices_tank.size();

                    if (nX > 0)
                    {
                        this->Flowerbounds.push_back(0.0);
                        this->Fupperbounds.push_back(1.0 / (1.0 + this->options.chemical_propellant_margin));
                        this->Fdescriptions.push_back("Stage " + std::to_string(stageIndex) + " chemical oxidizer tank constraint");

                        this->mySpacecraft.setChemicalOxidizerTank_Findex(stageIndex, this->Fdescriptions.size() - 1);

                        std::vector<size_t> stage_Gindex_PerStage_dChemicalOxidizerConstraint_dVirtualChemicalOxidizer;
                        for (size_t virtualTankIndex = 0; virtualTankIndex < nX; ++virtualTankIndex)
                        {
                            size_t Xindex = Xindices_tank[virtualTankIndex];

                            this->create_sparsity_entry(this->Fdescriptions.size() - 1,
                                Xindex,
                                stage_Gindex_PerStage_dChemicalOxidizerConstraint_dVirtualChemicalOxidizer);
                        }
                        this->Gindex_PerStage_dChemicalOxidizerConstraint_dVirtualChemicalOxidizer.push_back(stage_Gindex_PerStage_dChemicalOxidizerConstraint_dVirtualChemicalOxidizer);
                    }//end oxidizer tank
                }
                else
                    this->mySpacecraft.getStageOptions(stageIndex).setEnableChemicalPropellantTankConstraint(false);
            }//end stage chemical tank
        }//end stage propellant tanks

        //stage dry mass constraint
        for (size_t stageIndex = 0; stageIndex < this->mySpacecraft.getSpacecraftOptions().getNumberOfStages(); ++stageIndex)
        {
            if (this->mySpacecraft.getStageOptions(stageIndex).getEnableDryMassConstraint())
            {
                Flowerbounds.push_back(-math::LARGE);
                Fupperbounds.push_back(0.0);
                this->Fdescriptions.push_back("Stage " + std::to_string(stageIndex) + " dry mass constraint");


                //with respect to virtual electric tanks
                std::vector<size_t> Xindices_tank = this->mySpacecraft.getElectricPropellantTank_Xindices(stageIndex);
                size_t nX = Xindices_tank.size();

                for (size_t virtualTankIndex = 0; virtualTankIndex < nX; ++virtualTankIndex)
                {
                    size_t Gindex;
                    this->create_sparsity_entry(Fdescriptions.size() - 1,
                        Xindices_tank[virtualTankIndex],
                        Gindex);

                    this->mySpacecraft.appendDryMass_Gindices(stageIndex, Gindex);
                }

                //with respect to virtual chemical fuel tanks
                Xindices_tank = this->mySpacecraft.getChemicalFuelTank_Xindices(stageIndex);
                nX = Xindices_tank.size();

                for (size_t virtualTankIndex = 0; virtualTankIndex < nX; ++virtualTankIndex)
                {
                    size_t Gindex;
                    this->create_sparsity_entry(Fdescriptions.size() - 1,
                        Xindices_tank[virtualTankIndex],
                        Gindex);

                    this->mySpacecraft.appendDryMass_Gindices(stageIndex, Gindex);
                }

                //with respect to virtual chemical oxidizer tanks
                Xindices_tank = this->mySpacecraft.getChemicalOxidizerTank_Xindices(stageIndex);
                nX = Xindices_tank.size();

                for (size_t virtualTankIndex = 0; virtualTankIndex < nX; ++virtualTankIndex)
                {
                    size_t Gindex;
                    this->create_sparsity_entry(Fdescriptions.size() - 1,
                        Xindices_tank[virtualTankIndex],
                        Gindex);

                    this->mySpacecraft.appendDryMass_Gindices(stageIndex, Gindex);
                }

                //derivative entry with respect to stage final mass variable
                size_t Gindex;

                this->create_sparsity_entry(Fdescriptions.size() - 1,
                    this->mySpacecraft.getXindex_StageFinalMass(stageIndex),
                    Gindex);

                this->mySpacecraft.appendDryMass_Gindices(stageIndex, Gindex);
            }
        }


        this->total_number_of_NLP_parameters = Xupperbounds.size();
        this->total_number_of_constraints = Fupperbounds.size();

        this->X.resize(total_number_of_NLP_parameters, 0.0);
        this->X0.resize(total_number_of_NLP_parameters, 0.0);
        this->F.resize(total_number_of_constraints, 1.0e+100);
        this->G.resize(this->Gdescriptions.size());
        this->locate_equality_constraints();
    }//end calcbounds()

    void Mission::locate_filament_critical_inequality_constraints()
    {
        ////mission flight time constraint
        //for (size_t Findex = Fdescriptions.size() - 1; Findex >= 0; --Findex)
        //{
        //    if (Fdescriptions[Findex].find("Mission flight time bounds") < 1024)
        //    {
        //        this->F_indices_of_filament_critical_inequality_constraints.push_back(Findex);
        //        break;
        //    }
        //}
    }//end locate_filament_critical_inequality_constraints()


    //evaluate function
    void Mission::evaluate(const std::vector<doubleType>& X,
                          std::vector<doubleType>& F,
                          std::vector<double>& G,
                          const bool& needG)
    {
        size_t Xindex = 0;
        size_t Findex = 1; //F[0] is reserved for the objective function

        //reset delta-v
        this->total_deterministic_deltav = 0.0;
        this->total_statistical_deltav = 0.0;

        //reset staging
        this->mySpacecraft.resetStaging();

        //reset the Jacobian
        for (double& Gentry : G)
            Gentry = 0.0;

        //process all of the journeys
        for (size_t journeyIndex = 0; journeyIndex < this->number_of_journeys; ++journeyIndex)
        {
            this->Journeys[journeyIndex].process_journey(X,
                Xindex,
                F,
                Findex,
                G,
                needG);
        }

        // evaluate global time bounds
        this->TotalFlightTime = 0.0;
        for (size_t timeIndex = 0; timeIndex < this->timeconstraints_X_indices.size(); ++timeIndex)
        {
            this->TotalFlightTime += X[this->timeconstraints_X_indices[timeIndex]];
        }

        if (options.global_timebounded)
        {
            F[Findex++] = this->TotalFlightTime / this->options.total_flight_time_bounds[1];

            if (needG)
            {
                for (size_t timeIndex = 0; timeIndex < this->timeconstraints_X_indices.size(); ++timeIndex)
                {
                    G[this->timeconstraints_G_indices[timeIndex]] = this->X_scale_factors[this->timeconstraints_X_indices[timeIndex]] / this->options.total_flight_time_bounds[1];
                }
            }
        }

        //objective function
        this->myObjectiveFunction->process(X, Xindex, F, Findex, G, needG);

        //propellant tank constraints
        this->mySpacecraft.computePropellantState(X, this->options.electric_propellant_margin, this->options.chemical_propellant_margin);

        //first for the spacecraft
        if (this->mySpacecraft.getSpacecraftOptions().getEnableGlobalElectricPropellantTankConstraint())
        {
            F[Findex++] = this->mySpacecraft.getGlobalElectricPropellantUsed()
                / this->mySpacecraft.getSpacecraftOptions().getGlobalElectricPropellantTankCapacity();

            if (needG)
            {
                for (size_t tankIndex = 0; tankIndex < this->Gindex_dElectricPropellantConstraint_dVirtualElectricPropellant.size(); ++tankIndex)
                {
                    size_t Gindex = this->Gindex_dElectricPropellantConstraint_dVirtualElectricPropellant[tankIndex];
                    size_t Xindex = this->jGvar[Gindex];

                    G[Gindex] = this->X_scale_factors[Xindex]
                        / this->mySpacecraft.getSpacecraftOptions().getGlobalElectricPropellantTankCapacity();
                }
            }
        }//end global electric propellant tank

        if (this->mySpacecraft.getSpacecraftOptions().getEnableGlobalChemicalPropellantTankConstraint())
        {
            //fuel
            F[Findex++] = this->mySpacecraft.getGlobalChemicalFuelUsed()
                / this->mySpacecraft.getSpacecraftOptions().getGlobalChemicalFuelTankCapacity();

            if (needG)
            {
                for (size_t tankIndex = 0; tankIndex < this->Gindex_dChemicalFuelConstraint_dVirtualChemicalFuel.size(); ++tankIndex)
                {
                    size_t Gindex = this->Gindex_dChemicalFuelConstraint_dVirtualChemicalFuel[tankIndex];
                    size_t Xindex = this->jGvar[Gindex];

                    G[Gindex] = this->X_scale_factors[Xindex]
                        / this->mySpacecraft.getSpacecraftOptions().getGlobalChemicalFuelTankCapacity();
                }
            }
            
            //oxidizer
            if (this->mySpacecraft.getGlobalChemicalOxidizerTank_Xindices().size() > 0)
            {
                F[Findex++] = this->mySpacecraft.getGlobalChemicalOxidizerUsed()
                    / this->mySpacecraft.getSpacecraftOptions().getGlobalChemicalOxidizerTankCapacity();

                if (needG)
                {
                    for (size_t tankIndex = 0; tankIndex < this->Gindex_dChemicalOxidizerConstraint_dVirtualChemicalOxidizer.size(); ++tankIndex)
                    {
                        size_t Gindex = this->Gindex_dChemicalOxidizerConstraint_dVirtualChemicalOxidizer[tankIndex];
                        size_t Xindex = this->jGvar[Gindex];

                        G[Gindex] = this->X_scale_factors[Xindex]
                            / this->mySpacecraft.getSpacecraftOptions().getGlobalChemicalOxidizerTankCapacity();
                    }
                }
            }

        }//end global chemical propellant tanks

        //global dry mass constraint
        if (this->options.constrain_dry_mass)
        {
            size_t number_of_phases_in_final_journey = this->Journeys.back().getNumberOfPhases();
            BoundaryEvents::ArrivalEvent* myArrivalEvent = this->Journeys.back().getPhase(number_of_phases_in_final_journey - 1)->getArrivalEvent();
            math::Matrix<doubleType> FinalState = myArrivalEvent->get_state_after_event();

            doubleType dry_mass = FinalState(6)
                - this->mySpacecraft.getElectricPropellantMargin()
                - this->mySpacecraft.getChemicalFuelMargin()
                - this->mySpacecraft.getChemicalOxidizerMargin();

            F[Findex++] = dry_mass / this->options.final_mass_constraint_bounds[1] - 1.0;

            if (needG)
            {
                //wrt final state
                {
                    std::vector< std::tuple<size_t, size_t, double> > Derivatives_of_StateAfterArrival = myArrivalEvent->get_Derivatives_of_StateAfterEvent();
                    std::vector< std::tuple<size_t, size_t, double> > Derivatives_of_StateAfterArrival_wrt_Time = myArrivalEvent->get_Derivatives_of_StateAfterEvent_wrt_Time();

                    //non-time variables
                    for (size_t entryIndex = 0; entryIndex < this->dIndex_global_dry_mass_constraint.size(); ++entryIndex)
                    {
                        size_t dIndex = this->dIndex_global_dry_mass_constraint[entryIndex];
                        double TheDerivative = std::get<2>(Derivatives_of_StateAfterArrival[dIndex]);

                        size_t Gindex = this->Gindex_derivatives_of_global_dry_mass_constraint[entryIndex];
                        size_t Xindex = this->jGvar[Gindex];

                        G[Gindex] = this->X_scale_factors[Xindex]
                            * TheDerivative
                            / this->options.final_mass_constraint_bounds[1];
                    }

                    //time variables
                    for (size_t entryIndex = 0; entryIndex < this->dIndex_global_dry_mass_constraint_wrt_Time.size(); ++entryIndex)
                    {
                        size_t dIndex = this->dIndex_global_dry_mass_constraint_wrt_Time[entryIndex];
                        double TheDerivative = std::get<2>(Derivatives_of_StateAfterArrival_wrt_Time[dIndex]);

                        size_t Gindex = this->Gindex_derivatives_of_global_dry_mass_constraint[entryIndex];
                        size_t Xindex = this->jGvar[Gindex];

                        G[Gindex] = this->X_scale_factors[Xindex]
                            * TheDerivative
                            / this->options.final_mass_constraint_bounds[1];
                    }
                }

                //wrt electric propellant
                for (size_t Gindex : this->Gindex_dGlobalDryMassConstraint_dVirtualElectricPropellant)
                {
                    size_t Xindex = this->jGvar[Gindex];

                    G[Gindex] = -this->X_scale_factors[Xindex]
                        * this->options.electric_propellant_margin
                        / this->options.final_mass_constraint_bounds[1];
                }
                //wrt chemical fuel
                for (size_t Gindex : this->Gindex_dGlobalDryMassConstraint_dVirtualChemicalFuel)
                {
                    size_t Xindex = this->jGvar[Gindex];

                    G[Gindex] = -this->X_scale_factors[Xindex]
                        * this->options.chemical_propellant_margin
                        / this->options.final_mass_constraint_bounds[1];
                }

                //wrt chemical oxidizer
                for (size_t Gindex : this->Gindex_dGlobalDryMassConstraint_dVirtualChemicalOxidizer)
                {
                    size_t Xindex = this->jGvar[Gindex];

                    G[Gindex] = -this->X_scale_factors[Xindex]
                        * this->options.chemical_propellant_margin
                        / this->options.final_mass_constraint_bounds[1];
                }
            }//end derivatives
        }//end global dry mass constraint

        //global final mass constraint
        if (this->options.constrain_final_mass)
        {
            size_t number_of_phases_in_final_journey = this->Journeys.back().getNumberOfPhases();
            BoundaryEvents::ArrivalEvent* myArrivalEvent = this->Journeys.back().getPhase(number_of_phases_in_final_journey - 1)->getArrivalEvent();
            math::Matrix<doubleType> FinalState = myArrivalEvent->get_state_after_event();

            F[Findex++] = FinalState(6) / this->options.final_mass_constraint_bounds[1] - 1.0;

            if (needG)
            {
                std::vector< std::tuple<size_t, size_t, double> > Derivatives_of_StateAfterArrival = myArrivalEvent->get_Derivatives_of_StateAfterEvent();
                std::vector< std::tuple<size_t, size_t, double> > Derivatives_of_StateAfterArrival_wrt_Time = myArrivalEvent->get_Derivatives_of_StateAfterEvent_wrt_Time();

                //non-time variables
                for (size_t entryIndex = 0; entryIndex < this->dIndex_final_mass_constraint.size(); ++entryIndex)
                {
                    size_t dIndex = this->dIndex_final_mass_constraint[entryIndex];
                    double TheDerivative = std::get<2>(Derivatives_of_StateAfterArrival[dIndex]);

                    size_t Gindex = this->Gindex_derivatives_of_final_mass_constraint[entryIndex];
                    size_t Xindex = this->jGvar[Gindex];

                    G[Gindex] = this->X_scale_factors[Xindex]
                        * TheDerivative
                        / this->options.final_mass_constraint_bounds[1];
                }

                //time variables
                for (size_t entryIndex = 0; entryIndex < this->dIndex_final_mass_constraint_wrt_Time.size(); ++entryIndex)
                {
                    size_t dIndex = this->dIndex_final_mass_constraint_wrt_Time[entryIndex];
                    double TheDerivative = std::get<2>(Derivatives_of_StateAfterArrival_wrt_Time[dIndex]);

                    size_t Gindex = this->Gindex_derivatives_of_final_mass_constraint[entryIndex];
                    size_t Xindex = this->jGvar[Gindex];

                    G[Gindex] = this->X_scale_factors[Xindex]
                        * TheDerivative
                        / this->options.final_mass_constraint_bounds[1];
                }
            }
        }//end mission final mass constraint

        //then for individual stages
        for (size_t stageIndex = 0; stageIndex < this->mySpacecraft.getSpacecraftOptions().getNumberOfStages(); ++stageIndex)
        {
            if (this->mySpacecraft.getStageOptions(stageIndex).getEnableElectricPropellantTankConstraint())
            {
                F[Findex++] = this->mySpacecraft.getElectricPropellantUsed(stageIndex)
                    / this->mySpacecraft.getStageOptions(stageIndex).getElectricPropellantTankCapacity();

                if (needG)
                {
                    for (size_t tankIndex = 0; tankIndex < this->Gindex_PerStage_dElectricPropellantConstraint_dVirtualElectricPropellant[stageIndex].size(); ++tankIndex)
                    {
                        size_t Gindex = this->Gindex_PerStage_dElectricPropellantConstraint_dVirtualElectricPropellant[stageIndex][tankIndex];
                        size_t Xindex = this->jGvar[Gindex];

                        G[Gindex] = this->X_scale_factors[Xindex]
                            / this->mySpacecraft.getStageOptions(stageIndex).getElectricPropellantTankCapacity();
                    }
                }
            }//end stage electric propellant tank

            if (this->mySpacecraft.getStageOptions(stageIndex).getEnableChemicalPropellantTankConstraint())
            {
                //fuel
                F[Findex++] = this->mySpacecraft.getChemicalFuelUsed(stageIndex)
                    / this->mySpacecraft.getStageOptions(stageIndex).getChemicalFuelTankCapacity();

                if (needG)
                {
                    for (size_t tankIndex = 0; tankIndex < this->Gindex_PerStage_dChemicalFuelConstraint_dVirtualChemicalFuel[stageIndex].size(); ++tankIndex)
                    {
                        size_t Gindex = this->Gindex_PerStage_dChemicalFuelConstraint_dVirtualChemicalFuel[stageIndex][tankIndex];
                        size_t Xindex = this->jGvar[Gindex];

                        G[Gindex] = this->X_scale_factors[Xindex]
                            / this->mySpacecraft.getStageOptions(stageIndex).getChemicalFuelTankCapacity();
                    }
                }

                //oxidizer
                if (this->mySpacecraft.getChemicalOxidizerTank_Xindices(stageIndex).size() > 0)
                {
                    F[Findex++] = this->mySpacecraft.getChemicalOxidizerUsed(stageIndex)
                        / this->mySpacecraft.getStageOptions(stageIndex).getChemicalOxidizerTankCapacity();

                    if (needG)
                    {
                        for (size_t tankIndex = 0; tankIndex < this->Gindex_PerStage_dChemicalOxidizerConstraint_dVirtualChemicalOxidizer[stageIndex].size(); ++tankIndex)
                        {
                            size_t Gindex = this->Gindex_PerStage_dChemicalOxidizerConstraint_dVirtualChemicalOxidizer[stageIndex][tankIndex];
                            size_t Xindex = this->jGvar[Gindex];

                            G[Gindex] = this->X_scale_factors[Xindex]
                                / this->mySpacecraft.getStageOptions(stageIndex).getChemicalOxidizerTankCapacity();
                        }
                    }
                }
            }//end chemical tank
        }//end stage propellant

        //propellant tank derivatives are taken care of by SNOPT's linear constraint mode

        //stage dry mass constraint
        for (size_t stageIndex = 0; stageIndex < this->mySpacecraft.getSpacecraftOptions().getNumberOfStages(); ++stageIndex)
        {
            if (this->mySpacecraft.getStageOptions(stageIndex).getEnableDryMassConstraint())
            {
                this->mySpacecraft.setActiveStage(stageIndex);
                this->mySpacecraft.computeStageRequiredFinalMass(stageIndex);
                F[Findex++] = (this->mySpacecraft.getStageRequiredFinalMass(stageIndex) - X[this->mySpacecraft.getXindex_StageFinalMass(stageIndex)]) / this->mySpacecraft.getCurrentDryMass();

                if (needG)
                    this->mySpacecraft.populateDryMassDerivatives(stageIndex, G, this->options.electric_propellant_margin, this->options.chemical_propellant_margin);
            }
        }

        //compute total deterministic and statistical delta-v
        for (size_t journeyIndex = 0; journeyIndex < this->number_of_journeys; ++journeyIndex)
        {
            this->total_deterministic_deltav += this->Journeys[journeyIndex].getDeterministicDeltav();
            
            this->total_statistical_deltav += this->Journeys[journeyIndex].getStatisticalDeltav();
        }
    
        //test for errors
        if (std::isnan(F[0] _GETVALUE))
        {
            throw std::runtime_error("NaN in fitness value! Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
        }
    }//end evaluate()

    //output function
    void Mission::output(const std::string& outputfilename)
    {
        std::ofstream outputfile(outputfilename, std::ios::out | std::ios::trunc);
        outputfile << "Mission: " << options.mission_name << std::endl;
        outputfile << "Written by EMTG_v9 core program compiled " << __DATE__<< " " << __TIME__ << std::endl;

        //next, output summary lines describing each event in the mission
        int errcode = 0;
        size_t eventcount = 1;
        size_t jprint = 0;
        try
        {
            for (size_t journeyIndex = 0; journeyIndex < number_of_journeys; ++journeyIndex) 
            {
                this->Journeys[journeyIndex].output(outputfile, jprint, eventcount);
            }
        }
        catch (std::exception &error)
        {
            std::cout << "An error occurred while writing " << options.outputfile << std::endl;
            std::cout << error.what() << std::endl;
        }
    
        //skip 3 lines
        for (int k = 0; k < 2; ++k)
            outputfile << std::endl;

        outputfile.precision(20);
        //output total deltaV
        outputfile << "Total deterministic deltav (km/s): " << this->total_deterministic_deltav _GETVALUE << std::endl;
        outputfile << "Total statistical deltav (km/s): " << this->total_statistical_deltav _GETVALUE << std::endl;

        //output flight time in years
        outputfile << "Flight time (y): " << this->TotalFlightTime _GETVALUE / 86400.0 / 365.25 << std::endl;
        outputfile << std::endl;

        //output spacecraft mass and propellant
        this->mySpacecraft.output_mass_information(outputfile);
        outputfile << "Spacecraft: Final mass including propellant margin (kg): " << this->Journeys.back().getFinalMass() _GETVALUE << std::endl;
        outputfile << "Spacecraft: Dry mass (kg): " << (this->Journeys.back().getFinalMass()
            - this->mySpacecraft.getChemicalFuelMargin()
            - this->mySpacecraft.getChemicalOxidizerMargin()
            - this->mySpacecraft.getElectricPropellantMargin()) _GETVALUE << std::endl;
        outputfile << std::endl;

        //objective-function specific information
        this->myObjectiveFunction->output(outputfile);

        //worst constraint violation
        outputfile << std::endl;

        double feasibility, normalized_feasibility, distance_from_equality_filament, decision_variable_infeasibility;
        size_t worst_constraint, worst_decision_variable;

        this->check_feasibility(this->Xopt,
            this->F,
            worst_decision_variable,
            worst_constraint,
            feasibility,
            normalized_feasibility,
            distance_from_equality_filament,
            decision_variable_infeasibility);

        outputfile << "Worst constraint is F[" << worst_constraint << "]: " << this->Fdescriptions[worst_constraint] << std::endl;
        outputfile << "with violation " << feasibility << std::endl;
        outputfile << "Worst decision variable is X[" << worst_decision_variable << "]: " << this->Xdescriptions[worst_decision_variable] << std::endl;
        outputfile << "with unhappiness " << decision_variable_infeasibility << std::endl;

        //skip 3 lines
        for (int k = 0; k < 3; ++k)
            outputfile << std::endl;
    

        //finally, output the decision vector
        //upper bounds
        outputfile << "Xupperbounds";
        for (int k = 0; k < total_number_of_NLP_parameters; ++k)
        {
            if (this->Xdescriptions[k].find("epoch") < 1024 || this->Xdescriptions[k].find("time") < 1024)
            {
                outputfile << "," << this->Xupperbounds[k] / 86400.0;
            }
            else if (this->Xdescriptions[k].find("Sundman independent variable") < 1024)
            {
                //we need the journey index so that we can scale it by that journey's LU
                //prefix format is jJpP, so split the string by p and then remove the first character, convert to int
                std::vector<std::string> descriptionComponents;
                
                boost::split(descriptionComponents, this->Xdescriptions[k], boost::is_any_of("p"), boost::token_compress_on);

                std::string journeyPrefix = descriptionComponents[0].substr(1, descriptionComponents[0].size() - 1);

                size_t journeyIndex = std::stoi(journeyPrefix);

                outputfile << "," << this->Xupperbounds[k] / this->TheUniverse[journeyIndex].LU;
            }
            else
                outputfile << "," << this->Xupperbounds[k];
        }
        outputfile << std::endl;
        //decision vector
        outputfile << "Decision Vector:";
        //outputfile.unsetf(std::ios::fixed);
        //outputfile.precision(20);
        for (int k = 0; k < total_number_of_NLP_parameters; ++k)
        {
            if (this->Xdescriptions[k].find("epoch") < 1024 || this->Xdescriptions[k].find("time") < 1024)
            {
                outputfile << "," << this->Xopt[k]_GETVALUE / 86400.0;
            }
            else if (this->Xdescriptions[k].find("Sundman independent variable") < 1024)
            {
                //we need the journey index so that we can scale it by that journey's LU
                //prefix format is jJpP, so split the string by p and then remove the first character, convert to int
                std::vector<std::string> descriptionComponents;

                boost::split(descriptionComponents, this->Xdescriptions[k], boost::is_any_of("p"), boost::token_compress_on);

                std::string journeyPrefix = descriptionComponents[0].substr(1, descriptionComponents[0].size() - 1);

                size_t journeyIndex = std::stoi(journeyPrefix);

                outputfile << "," << this->Xopt[k] / this->TheUniverse[journeyIndex].LU;
            }
            else
                outputfile << "," << this->Xopt[k]_GETVALUE;
        }
        outputfile << std::endl;
        //lower bounds
        outputfile << "Xlowerbounds";
        for (int k = 0; k < total_number_of_NLP_parameters; ++k)
        {
            if (this->Xdescriptions[k].find("epoch") < 1024 || this->Xdescriptions[k].find("time") < 1024)
            {
                outputfile << "," << this->Xlowerbounds[k] / 86400.0;
            }
            else if (this->Xdescriptions[k].find("Sundman independent variable") < 1024)
            {
                //we need the journey index so that we can scale it by that journey's LU
                //prefix format is jJpP, so split the string by p and then remove the first character, convert to int
                std::vector<std::string> descriptionComponents;

                boost::split(descriptionComponents, this->Xdescriptions[k], boost::is_any_of("p"), boost::token_compress_on);

                std::string journeyPrefix = descriptionComponents[0].substr(1, descriptionComponents[0].size() - 1);

                size_t journeyIndex = std::stoi(journeyPrefix);

                outputfile << "," << this->Xlowerbounds[k] / this->TheUniverse[journeyIndex].LU;
            }
            else
                outputfile << "," << this->Xlowerbounds[k];
        }
        outputfile << std::endl;
        //descriptions
        outputfile << "Xdescriptions";
        for (int k = 0; k < total_number_of_NLP_parameters; ++k)
            outputfile << "," << Xdescriptions[k];
        outputfile << std::endl;

        outputfile << std::endl;
        outputfile << std::endl;

        //output the constraints
        //upper bounds
        outputfile << "Fupperbounds.";
        for (int k = 0; k < total_number_of_constraints; ++k)
            outputfile << "," << Fupperbounds[k];
        outputfile << std::endl;
        //constraint values
        outputfile << "Constraint_Vector";
        for (int k = 0; k < total_number_of_constraints; ++k)
            outputfile << "," << F[k]_GETVALUE;
        outputfile << std::endl;
        //lower bounds
        outputfile << "Flowerbounds";
        for (int k = 0; k < total_number_of_constraints; ++k)
            outputfile << "," << Flowerbounds[k];
        outputfile << std::endl;
        //descriptions
        outputfile << "Fdescriptions";
        for (int k = 0; k < total_number_of_constraints; ++k)
            outputfile << "," << Fdescriptions[k];
        outputfile << std::endl;
        
        outputfile << std::endl;
        outputfile << "user_data: " << this->options.user_data << std::endl;
        outputfile << std::endl;

        outputfile.close();
    }

    void Mission::output_ephemeris()
    {
        std::ofstream acceleration_model_file;
        const size_t bufsize = 256*1024;
        char buf1[bufsize];
        char buf2[bufsize];
        if (this->options.generate_acceleration_model_instrumentation_file)
        {
            std::string acceleration_model_string = this->options.working_directory + "//" + options.mission_name + ".accelerations";
            acceleration_model_file.open(acceleration_model_string, std::ios::out | std::ios::trunc);
            
            acceleration_model_file.rdbuf()->pubsetbuf(buf1, bufsize);

            acceleration_model_file.setf(std::ios::scientific, std::ios::floatfield);
            acceleration_model_file.precision(16);
        }

        //Step 1: write ephemeris file
        std::string ephemeris_string = this->options.working_directory + "//" + options.mission_name + ".ephemeris";
        std::ofstream ephemeris_file(ephemeris_string, std::ios::out | std::ios::trunc);
        ephemeris_file.rdbuf()->pubsetbuf(buf2, bufsize);
        ephemeris_file << "#epoch, x(km), y(km), z(km), vx(km/s), vy(km/s), vz(km/s)";

        if (this->options.append_mass_to_ephemeris_output)
            ephemeris_file << ", mass(kg)";

        if (this->options.append_control_to_ephemeris_output)
            ephemeris_file << ", ux, uy, uz";

        if (this->options.append_thrust_to_ephemeris_output)
            ephemeris_file << ", ThrustMagnitude(N)";

        if (this->options.append_mdot_to_ephemeris_output)
            ephemeris_file << ", MassFlowRate(kg/s)";

        if (this->options.append_Isp_to_ephemeris_output)
            ephemeris_file << ", Isp(s)";

        if (this->options.append_number_of_active_engines_to_ephemeris_output)
            ephemeris_file << ", NumberOfActiveThrusters";

        if (this->options.append_active_power_to_ephemeris_output)
            ephemeris_file << ", ActivePower(kW)";

        if (this->options.append_throttle_level_to_ephemeris_output)
            ephemeris_file << ", ThrottleLevel";

        //Step 9: newline
        ephemeris_file << std::endl;

        for (size_t journeyIndex = 0; journeyIndex < this->number_of_journeys; ++journeyIndex)
        {
            if (this->options.generate_acceleration_model_instrumentation_file)
            {
                if (journeyIndex > 0)
                {
                    acceleration_model_file << "\n";
                }
                acceleration_model_file << "#epoch (JD)";
                acceleration_model_file << ", sc_distance_from_cb" << ", x" << ", y" << ", z";
                acceleration_model_file << ", velocity_norm" << ", vx" << ", vy" << ", vz";
                acceleration_model_file << ", mass";
                acceleration_model_file << ", acceleration_norm" << ", ax" << ", ay" << ", az";;
                if (this->options.perturb_SRP)
                {
                    acceleration_model_file << ", SRP_norm, SRP_x, SRP_y, SRP_z";
                }
				for (size_t k = 0; k < this->options.Journeys.size(); ++k)
				{
					if (this->options.Journeys[k].perturb_drag)
					{
						acceleration_model_file << ", Drag_norm, Drag_x, Drag_y, Drag_z";
						break;
					}
				}
                acceleration_model_file << ", cb_point_mass_gravity_norm";
                acceleration_model_file << ", cb_point_mass_gravity_x";
                acceleration_model_file << ", cb_point_mass_gravity_y";
                acceleration_model_file << ", cb_point_mass_gravity_z";
                if (this->options.perturb_J2)
                {
                    acceleration_model_file << ", cb_J2_gravity_norm";
                    acceleration_model_file << ", cb_J2_gravity_x";
                    acceleration_model_file << ", cb_J2_gravity_y";
                    acceleration_model_file << ", cb_J2_gravity_z";
                }
                if (this->options.perturb_thirdbody)
                {
                    for (size_t k = 0; k < this->options.Journeys[journeyIndex].perturbation_bodies.size(); ++k)
                    {
                        acceleration_model_file << ", " + this->TheUniverse[journeyIndex].bodies[k].name + "_distance_from_cb";
                        acceleration_model_file << ", " + this->TheUniverse[journeyIndex].bodies[k].name + "_x";
                        acceleration_model_file << ", " + this->TheUniverse[journeyIndex].bodies[k].name + "_y";
                        acceleration_model_file << ", " + this->TheUniverse[journeyIndex].bodies[k].name + "_z";
                        acceleration_model_file << ", " + this->TheUniverse[journeyIndex].bodies[k].name + "_speed_wrt_cb";
                        acceleration_model_file << ", " + this->TheUniverse[journeyIndex].bodies[k].name + "_vx";
                        acceleration_model_file << ", " + this->TheUniverse[journeyIndex].bodies[k].name + "_vy";
                        acceleration_model_file << ", " + this->TheUniverse[journeyIndex].bodies[k].name + "_vz";
                        acceleration_model_file << ", sc_distance_from_" + this->TheUniverse[journeyIndex].bodies[k].name;
                        acceleration_model_file << ", " + this->TheUniverse[journeyIndex].bodies[k].name + "_2_sc_x";
                        acceleration_model_file << ", " + this->TheUniverse[journeyIndex].bodies[k].name + "_2_sc_y";
                        acceleration_model_file << ", " + this->TheUniverse[journeyIndex].bodies[k].name + "_2_sc_z";
                        acceleration_model_file << ", sc_speed_wrt_" + this->TheUniverse[journeyIndex].bodies[k].name;
                        acceleration_model_file << ", " + this->TheUniverse[journeyIndex].bodies[k].name + "_2_sc_vx";
                        acceleration_model_file << ", " + this->TheUniverse[journeyIndex].bodies[k].name + "_2_sc_vy";
                        acceleration_model_file << ", " + this->TheUniverse[journeyIndex].bodies[k].name + "_2_sc_vz";
                        acceleration_model_file << ", " + this->TheUniverse[journeyIndex].bodies[k].name + "_point_mass_gravity_norm";
                        acceleration_model_file << ", " + this->TheUniverse[journeyIndex].bodies[k].name + "_point_mass_gravity_x";
                        acceleration_model_file << ", " + this->TheUniverse[journeyIndex].bodies[k].name + "_point_mass_gravity_y";
                        acceleration_model_file << ", " + this->TheUniverse[journeyIndex].bodies[k].name + "_point_mass_gravity_z";
                    }
                }
                if (this->options.Journeys[journeyIndex].phase_type == 1 ||
                    this->options.Journeys[journeyIndex].phase_type == 3 ||
                    this->options.Journeys[journeyIndex].phase_type == 5)
                {
                    acceleration_model_file << ", thruster_norm, thruster_x, thruster_y, thruster_z";
                }
                acceleration_model_file.flush();
            }
            this->Journeys[journeyIndex].output_ephemeris(ephemeris_file, acceleration_model_file);
            
        }

        ephemeris_file.close();
        acceleration_model_file.close();

        //Step 2: write .cmd file
        std::string cmdstring = this->options.working_directory + "//" + options.mission_name + ".cmd";
        std::ofstream cmd_file(cmdstring, std::ios::out | std::ios::trunc);

        cmd_file << "\\begindata" << std::endl;
        cmd_file << "INPUT_DATA_TYPE = 'STATES'" << std::endl;
        cmd_file << "OUTPUT_SPK_TYPE = 9" << std::endl;
        cmd_file << "OBJECT_ID = " << this->options.spacecraft_SPICE_ID << std::endl;
        cmd_file << "OBJECT_NAME = '" << this->options.mission_name << "'" << std::endl;
        cmd_file << "CENTER_ID = " << this->options.forward_integrated_ephemeris_central_body_SPICE_ID << std::endl;
        cmd_file << "CENTER_NAME = 'CENTER'" << std::endl;
        cmd_file << "REF_FRAME_NAME = 'J2000'" << std::endl;
        cmd_file << "PRODUCER_ID = 'EMTGv9'" << std::endl;
        cmd_file << "DATA_ORDER = 'EPOCH X Y Z VX VY VZ'" << std::endl;
        cmd_file << "INPUT_DATA_UNITS = ('ANGLES=DEGREES' 'DISTANCES=km')" << std::endl;
        cmd_file << "TIME_WRAPPER = '# TDB'" << std::endl;
        cmd_file << "DATA_DELIMITER = ','" << std::endl;
        cmd_file << "LINES_PER_RECORD = 1" << std::endl;
        cmd_file << "IGNORE_FIRST_LINE = 0" << std::endl;
        cmd_file << "LEAPSECONDS_FILE = '" << boost::replace_all_copy(this->options.universe_folder, "\\", "/") + "/ephemeris_files/" + this->options.SPICE_leap_seconds_kernel << "'" << std::endl;
        cmd_file << "POLYNOM_DEGREE = 3" << std::endl;
        cmd_file << "SEGMENT_ID = 'SPK_STATES_09'" << std::endl;
        cmd_file << "\\begintext" << std::endl;

        cmd_file.close();

        //Step 3: write a Python file that makes and checks the .bsp
        std::ofstream pyfile(this->options.working_directory + "//bspwriter.py", std::ios::out | std::ios::trunc);

        pyfile << "import os" << std::endl;
		if (strcmp(this->options.spice_utility_extension.c_str(),"\"\""))
		{
        	pyfile << "mkspk_path = '" << boost::replace_all_copy(this->options.spice_utilities_path, "\\", "/") << "/mkspk'" << std::endl;
        	pyfile << "brief_path = '" << boost::replace_all_copy(this->options.spice_utilities_path, "\\", "/") << "/brief'" << std::endl;
		}
		else
		{
        	pyfile << "mkspk_path = '" << boost::replace_all_copy(this->options.spice_utilities_path, "\\", "/") << "/mkspk" << this->options.spice_utility_extension << "'" << std::endl;
        	pyfile << "brief_path = '" << boost::replace_all_copy(this->options.spice_utilities_path, "\\", "/") << "/brief" << this->options.spice_utility_extension << "'" << std::endl;
		}
		
		pyfile << "import sys " << std::endl;
		pyfile << "sys.path.append('" << boost::replace_all_copy(this->options.pyemtg_path, "\\", "/") << "')" << std::endl;
		pyfile << "sys.path.append('" << boost::replace_all_copy(this->options.pyemtg_path, "\\", "/") << "/SimpleMonteCarlo')" << std::endl;
		pyfile << "sys.path.append('" << boost::replace_all_copy(this->options.pyemtg_path, "\\", "/") << "/SpiceyPy_Utilities') " << std::endl;

		pyfile << "import clean_spiceicles " << std::endl;

		pyfile << "clean_spiceicles.do_the_stuff(['" << boost::replace_all_copy(this->options.working_directory, "\\", "/") << "','" << options.mission_name << ".ephemeris','" << boost::replace_all_copy(this->options.universe_folder, "\\", "/") << "/ephemeris_files/']) " << std::endl;
		
		
        pyfile << "if os.path.exists('" << options.mission_name << ".bsp') :" << std::endl;
        pyfile << "    os.remove('" << options.mission_name << ".bsp')" << std::endl;
        pyfile << "os.system(mkspk_path + ' -setup " << options.mission_name << ".cmd -input " << options.mission_name << "_clean.ephemeris -output " << options.mission_name << ".bsp')" << std::endl;
        pyfile << "os.system(brief_path + ' " << options.mission_name << ".bsp')" << std::endl;

        pyfile.close();
		
		if (this->options.call_system_to_generate_bsp)
		{
			std::string cmd = "cd " + this->options.working_directory + "; python bspwriter.py > bsp.log";
			system(cmd.c_str());
		}
    }//end output_ephemeris()

    void Mission::output_STMs()
    {
        for (size_t journeyIndex = 0; journeyIndex < this->number_of_journeys; ++journeyIndex)
            this->Journeys[journeyIndex].output_STMs();
    }

    void Mission::output_maneuver_and_target_spec()
    {
        std::ofstream maneuver_spec_file(this->options.working_directory + "//" + this->options.mission_name + ".mission_maneuver_spec");
        maneuver_spec_file << "<EVENTNAME>,<NUMBER_OF_MANEUVERS>,<FRAME>,<EPOCH(ET)>,<THRX>,<THRY>,<THRZ>,<THRMAG[N]>,<SMASS[kg]>,<MDOT[kg/s]>,<DUTY>,<FMASS[kg]>,<DV[km/s]>,<DUR[s]>,<EPOCH(ET seconds)>, repeat..." << std::endl;

        std::ofstream target_spec_file(this->options.working_directory + "//" + this->options.mission_name + ".mission_target_spec");
        target_spec_file << "<EVENTNAME>,<SAMPLE>,<CBODY>,<FRAME>,<EPOCH(ET)>,<X[km]>,<Y[km]>,<Z[km]>,<VX[km/s]>,<VY[km/s]>,<VZ[km/s]>,<MASS[kg]>,<B.R[km]>,<B.T[km]>,<EPOCH(ET seconds)>" << std::endl;

        bool haveManeuverNeedTarget = false;
        for (size_t journeyIndex = 0; journeyIndex < this->number_of_journeys; ++journeyIndex)
            this->Journeys[journeyIndex].output_maneuver_and_target_spec(maneuver_spec_file, target_spec_file, haveManeuverNeedTarget);
    }

    std::vector<double> Mission::construct_initial_guess()
    {
        //Step 1: make a random number generator and a distribution
        std::mt19937 RNG(time(0));
        std::uniform_real_distribution<> DoubleDistribution(0.0, 1.0);

        //Step 2: interpolate the guess
        std::vector<std::tuple<std::string, double> > processed_initial_guess;
        size_t Xindex = 0;

        while (Xindex < this->options.trialX.size())
        {
            //Step 2.1 is this a control variable or parallel shooting block? if so go to control interpolation logic
            if (std::get<0>(this->options.trialX[Xindex]).find("PSFB_Step0") < 1024
                || std::get<0>(this->options.trialX[Xindex]).find("PSBI_Step0") < 1024) //parallel shooting block
            {
                //Step 2.1.1: extract the step block
                bool reached_end_of_step_block = false;
                std::vector<std::tuple<std::string, double> > step_block_initial_guess;
                while (!reached_end_of_step_block)
                {
                    if (std::get<0>(this->options.trialX[Xindex]).find("PSFB_Step") < 1024
                        || std::get<0>(this->options.trialX[Xindex]).find("PSBI_Step") < 1024)
                        step_block_initial_guess.push_back(this->options.trialX[Xindex++]);
                    else
                        reached_end_of_step_block = true;

                    if (Xindex == this->options.trialX.size())//hey, we reached the end of the mission
                        reached_end_of_step_block = true;
                }//end while (!reached_end_of_step_block)

                //Step 2.1.2: determine if the control block needs to be interpolated, and if so do it
                //Step 2.1.2.1: first we need to know the number of variables in a block
                size_t size_of_step_block = 0;
                size_t number_of_interior_control_points_in_guess = 0;
                size_t num_controls = 3;

                for (size_t stepIndex = 0; stepIndex < step_block_initial_guess.size(); ++stepIndex)
                {
                    if (std::get<0>(step_block_initial_guess[stepIndex]).find("Step0") < 1024)
                    {
                        ++size_of_step_block;
                    }
                    if (std::get<0>(step_block_initial_guess[stepIndex]).find("u_x") < 1024)
                    {
                        ++number_of_interior_control_points_in_guess;
                    }
                    if (std::get<0>(step_block_initial_guess[stepIndex]).find("u_command") < 1024)
                    {
                        num_controls = 4;
                    }
                    else if (std::get<0>(step_block_initial_guess[stepIndex]).find("Step1") < 1024)
                    {
                        break;
                    }
                }//end loop to determine number of controls
				
                 //Step 2.1.2.2: now, knowing the number of controls, we can find the number of control steps
                size_t number_of_guess_steps = step_block_initial_guess.size() / size_of_step_block;

                //Step 2.1.2.3: if the number of initial guess steps does not equal the number of steps in the problem, interpolate
                //this needs to be right-sized to match the journey of interest, which may have overriden number of steps
                //for this we will need the journeyIndex and then we look up the right number of steps

                //first get the prefix for the current entry, which tells us which phase and journey we are in
                std::vector<std::string> parsed_description;
                boost::split(parsed_description, std::get<0>(step_block_initial_guess.front()), boost::is_any_of(":"));
                std::string temp = parsed_description[0];
                std::string prefix = temp.substr(0, temp.find("_"));

                size_t journeyIndex = std::stoi(prefix.substr(1, prefix.find_first_of("p")));
                if (journeyIndex > this->number_of_journeys - 1)
                    break;

                size_t num_steps = this->options.Journeys[journeyIndex].override_num_steps
                    ? this->options.Journeys[journeyIndex].number_of_steps
                    : this->options.num_timesteps;

                if (number_of_guess_steps != num_steps 
                    || number_of_interior_control_points_in_guess != this->options.Journeys[journeyIndex].num_interior_control_points)
                {
					
                    std::vector<std::tuple<std::string, double> > new_step_block_initial_guess;

                    //construct a data table of all control values
                    //we will do this by looping through the old vector and extracting the control values
                    //if this is a parallel shooting mission then we will also extract the state values - maybe we want to worry about that later, or even have the individual phases do it?
                    std::vector< std::vector< std::pair<double, double> > > Old_6state(6);
                    std::vector< std::pair<double, double> > Old_mass;
                    std::vector< std::pair<double, double> > Old_chemfuel;
                    std::vector< std::pair<double, double> > Old_electricpropellant;
                    std::vector< std::pair<double, double> > Oldu_x;
                    std::vector< std::pair<double, double> > Oldu_y;
                    std::vector< std::pair<double, double> > Oldu_z;
                    std::vector< std::pair<double, double> > Oldu_command;



                    size_t oldStepBlockIndex = 0;
                    //then each step
                    for (int step = 0; step < number_of_guess_steps; ++step)
                    {
                        size_t interior_control_substep = 0;

                        for (size_t stateIndex : {0, 1, 2, 3, 4, 5})
                        {
                            Old_6state[stateIndex].push_back({ (double)step / number_of_guess_steps, std::get<1>(step_block_initial_guess[oldStepBlockIndex++]) });
                        }

                        Old_mass.push_back(std::make_pair((double)step / number_of_guess_steps, std::get<1>(step_block_initial_guess[oldStepBlockIndex++])));
                        Old_chemfuel.push_back(std::make_pair((double)step / number_of_guess_steps, std::get<1>(step_block_initial_guess[oldStepBlockIndex++])));
                        Old_electricpropellant.push_back(std::make_pair((double)step / number_of_guess_steps, std::get<1>(step_block_initial_guess[oldStepBlockIndex++])));

                        for (size_t subStepIndex = 0; subStepIndex < number_of_interior_control_points_in_guess; ++subStepIndex)
                        {
                            Oldu_x.push_back(std::make_pair(((double)step + (double)subStepIndex / number_of_interior_control_points_in_guess) / number_of_guess_steps, std::get<1>(step_block_initial_guess[oldStepBlockIndex++])));
                            Oldu_y.push_back(std::make_pair(((double)step + (double)subStepIndex / number_of_interior_control_points_in_guess) / number_of_guess_steps, std::get<1>(step_block_initial_guess[oldStepBlockIndex++])));
                            Oldu_z.push_back(std::make_pair(((double)step + (double)subStepIndex / number_of_interior_control_points_in_guess) / number_of_guess_steps, std::get<1>(step_block_initial_guess[oldStepBlockIndex++])));

                            if (num_controls == 4)
                                Oldu_command.push_back(std::make_pair(((double)step + (double)subStepIndex / number_of_interior_control_points_in_guess) / number_of_guess_steps, std::get<1>(step_block_initial_guess[oldStepBlockIndex++])));
                        }
                    }

                    //and finally the right-hand side
                    size_t offset = size_of_step_block;

                    std::vector<std::string> stateNames(6);
                    for (size_t stateIndex : {0, 1, 2, 3, 4, 5})
                    {
                        std::vector<std::string> stateNameCell;
                        boost::split(stateNameCell, std::get<0>(step_block_initial_guess[oldStepBlockIndex - offset]), boost::is_any_of(":"), boost::token_compress_on);

                        stateNames[stateIndex] = stateNameCell.back();

                        Old_6state[stateIndex].push_back({ 1.0, std::get<1>(step_block_initial_guess[oldStepBlockIndex - offset--]) });
                    }

                    Old_mass.push_back(std::make_pair(1.0, std::get<1>(step_block_initial_guess[oldStepBlockIndex - offset--])));
                    Old_chemfuel.push_back(std::make_pair(1.0, std::get<1>(step_block_initial_guess[oldStepBlockIndex - offset--])));
                    Old_electricpropellant.push_back(std::make_pair(1.0, std::get<1>(step_block_initial_guess[oldStepBlockIndex - offset--])));
                    Oldu_x.push_back(std::make_pair(1.0, std::get<1>(step_block_initial_guess[oldStepBlockIndex - offset--])));
                    Oldu_y.push_back(std::make_pair(1.0, std::get<1>(step_block_initial_guess[oldStepBlockIndex - offset--])));
                    Oldu_z.push_back(std::make_pair(1.0, std::get<1>(step_block_initial_guess[oldStepBlockIndex - offset--])));


                    if (num_controls == 4)
                        Oldu_command.push_back(std::make_pair(1.0, std::get<1>(step_block_initial_guess[oldStepBlockIndex - offset--])));
                    
                    //which variables need to be wrapped around the clock?
                    std::vector<size_t> wrapIndices;
                    switch (this->options.ParallelShootingStateRepresentation)
                    {
                        case StateRepresentation::COE:
                            {
                                wrapIndices = { 3, 4, 5};
                                break;
                            }
                        case StateRepresentation::SphericalAZFPA:
                        {
                            wrapIndices = { 1, 4 };
                            break;
                        }
                        case StateRepresentation::SphericalRADEC:
                        {
                            wrapIndices = { 1, 4 };
                            break;
                        }
                        case StateRepresentation::MEE:
                        {
                            wrapIndices = { 5 };
                            break;
                        }
                    }

                    //now UNWRAP all wrapped elements so that the interpolator works properly
                    for (size_t wrapIndex : wrapIndices) //loop over wrapped variables
                    {
                        for (size_t timeStep = 1; timeStep < Old_6state[wrapIndex].size(); ++timeStep)
                        {
                            //if the state variable at this time step is more than 2*pi larger than the previous time step
                            //checking for 2*pi is not sufficient because the difference can become less than that due to propagation. So we use a 10% fudge factor.
                            while (std::get<1>(Old_6state[wrapIndex][timeStep]) > std::get<1>(Old_6state[wrapIndex][timeStep - 1]) + math::TwoPI * 0.9)
                            {
                                std::get<1>(Old_6state[wrapIndex][timeStep]) -= math::TwoPI;
                            }

                            //if the state variable at this itme step is more than 2*pi less than the previous time step
                            //checking for 2*pi is not sufficient because the difference can become less than that due to propagation. So we use a 10% fudge factor.
                            while (std::get<1>(Old_6state[wrapIndex][timeStep]) < std::get<1>(Old_6state[wrapIndex][timeStep - 1]) - math::TwoPI * 0.9)
                            {
                                std::get<1>(Old_6state[wrapIndex][timeStep]) += math::TwoPI;
                            }
                        }
                    }

                    //interpolate the data tables to fill in the new control vector
                    math::interpolator state0_interpolator(Old_6state[0]);
                    math::interpolator state1_interpolator(Old_6state[1]);
                    math::interpolator state2_interpolator(Old_6state[2]);
                    math::interpolator state3_interpolator(Old_6state[3]);
                    math::interpolator state4_interpolator(Old_6state[4]);
                    math::interpolator state5_interpolator(Old_6state[5]);
                    math::interpolator mass_interpolator(Old_mass);
                    math::interpolator chemfuel_interpolator(Old_chemfuel);
                    math::interpolator electricpropellant_interpolator(Old_electricpropellant);
                    math::interpolator ux_interpolator(Oldu_x);
                    math::interpolator uy_interpolator(Oldu_y);
                    math::interpolator uz_interpolator(Oldu_z);
                    math::interpolator ucommand_interpolator(Oldu_command);

                    for (int step = 0; step < num_steps; ++step)
                    {
                        double current_normalized_integration_variable = (step + 0.5) / num_steps;
                        std::stringstream stepstream;
                        stepstream << step;
                        
                        //interpolate and control the state vector
                        new_step_block_initial_guess.push_back({ prefix + "_Step" + std::to_string(step) + ":" + stateNames[0], state0_interpolator.interpolate(current_normalized_integration_variable) });
                        new_step_block_initial_guess.push_back({ prefix + "_Step" + std::to_string(step) + ":" + stateNames[1], state1_interpolator.interpolate(current_normalized_integration_variable) });
                        new_step_block_initial_guess.push_back({ prefix + "_Step" + std::to_string(step) + ":" + stateNames[2], state2_interpolator.interpolate(current_normalized_integration_variable) });
                        new_step_block_initial_guess.push_back({ prefix + "_Step" + std::to_string(step) + ":" + stateNames[3], state3_interpolator.interpolate(current_normalized_integration_variable) });
                        new_step_block_initial_guess.push_back({ prefix + "_Step" + std::to_string(step) + ":" + stateNames[4], state4_interpolator.interpolate(current_normalized_integration_variable) });
                        new_step_block_initial_guess.push_back({ prefix + "_Step" + std::to_string(step) + ":" + stateNames[5], state5_interpolator.interpolate(current_normalized_integration_variable) });

                        //now re-wrap to (0, 2*pi)
                        for (size_t wrapIndex : wrapIndices)
                        {
                            double value = std::get<1>(new_step_block_initial_guess[new_step_block_initial_guess.size() - 6 + wrapIndex]);

                            std::get<1>(new_step_block_initial_guess[new_step_block_initial_guess.size() - 6 + wrapIndex]) = fmod(value, math::TwoPI);
                        }

                        new_step_block_initial_guess.push_back({ prefix + "_Step" + std::to_string(step) + ": left state mass", mass_interpolator.interpolate(current_normalized_integration_variable) });
                        new_step_block_initial_guess.push_back({ prefix + "_Step" + std::to_string(step) + ": virtual chemical fuel", chemfuel_interpolator.interpolate(current_normalized_integration_variable) });
                        new_step_block_initial_guess.push_back({ prefix + "_Step" + std::to_string(step) + ": virtual electric propellant", electricpropellant_interpolator.interpolate(current_normalized_integration_variable) });

                        //interpolate and encode the control vector
                        for (size_t subStepIndex = 0; subStepIndex < this->options.Journeys[journeyIndex].num_interior_control_points; ++subStepIndex)
                        {
                            double current_substep_normalized_integration_variable = (step + ((double)subStepIndex + 0.5 / this->options.Journeys[journeyIndex].num_interior_control_points)) / num_steps;

                            double ux = ux_interpolator.interpolate(current_substep_normalized_integration_variable);
                            double uy = uy_interpolator.interpolate(current_substep_normalized_integration_variable);
                            double uz = uz_interpolator.interpolate(current_substep_normalized_integration_variable);
                            double umag = sqrt(ux*ux + uy * uy + uz * uz);

                            /*if (umag > 1.0)
                            {
                                ux /= umag;
                                uy /= umag;
                                uz /= umag;
                            }*/

                            new_step_block_initial_guess.push_back({ prefix + "_Step" + std::to_string(step) + ": substep" + std::to_string(subStepIndex) + " u_x", ux });
                            new_step_block_initial_guess.push_back({ prefix + "_Step" + std::to_string(step) + ": substep" + std::to_string(subStepIndex) + " u_y", uy });
                            new_step_block_initial_guess.push_back({ prefix + "_Step" + std::to_string(step) + ": substep" + std::to_string(subStepIndex) + " u_z", uz });

                            if (num_controls == 4)
                                new_step_block_initial_guess.push_back({ prefix + "_Step" + std::to_string(step) + ": substep" + std::to_string(subStepIndex) + " u_command", ucommand_interpolator.interpolate(current_substep_normalized_integration_variable) });
                        }
                    }

                    step_block_initial_guess = new_step_block_initial_guess;
                    
                    /*std::ofstream dumpfile("./tests/initial_guess.csv", std::ios::trunc);
                    dumpfile << "Variable name, value" << std::endl;
                    for (std::tuple<std::string, double> guessline : step_block_initial_guess)
                        dumpfile << std::get<0>(guessline) << "," << std::get<1>(guessline) << std::endl;
                    dumpfile.close();*/
                }

                for (size_t controlIndex = 0; controlIndex < step_block_initial_guess.size(); ++controlIndex)
                {
                    processed_initial_guess.push_back(step_block_initial_guess[controlIndex]);
                }//end loop to to copy controls
            }
            else if (std::get<0>(this->options.trialX[Xindex]).find("u_") < 1024) //two-point shooting control
            {
                //Step 2.1.1: extract the control block
                bool reached_end_of_control_block = false;
                std::vector<std::tuple<std::string, double> > control_block_initial_guess;
                //first we have to know how many steps there were in the initial guess
                while (!reached_end_of_control_block)
                {
                    if (std::get<0>(this->options.trialX[Xindex]).find("u_") < 1024)
                        control_block_initial_guess.push_back(this->options.trialX[Xindex++]);
                    else
                        reached_end_of_control_block = true;

                    if (Xindex == this->options.trialX.size())//hey, we reached the end of the mission
                        reached_end_of_control_block = true;
                }//end while (!reached_end_of_control_block)

                //Step 2.1.2: determine if the control block needs to be interpolated, and if so do it
                //Step 2.1.2.1: first we need to know the number of controls
                int number_of_controls = 0;
                for (size_t controlIndex = 0; controlIndex < control_block_initial_guess.size(); ++controlIndex)
                {
                    if (std::get<0>(control_block_initial_guess[controlIndex]).find("step 0") < 1024)
                    {
                        number_of_controls++;
                    }
                    else if (std::get<0>(control_block_initial_guess[controlIndex]).find("step 1") < 1024)
                    {
                        break;
                    }
                }//end loop to determine number of controls

                //Step 2.1.2.2: now, knowing the number of controls, we can find the number of control steps
                size_t number_of_guess_control_steps = control_block_initial_guess.size() / number_of_controls;

                //Step 2.1.2.3: if the number of initial guess steps does not equal the number of steps in the problem, interpolate
                //this needs to be right-sized to match the journey of interest, which may have overriden number of steps
                //for this we will need the journeyIndex and then we look up the right number of steps

                //first get the prefix for the current entry, which tells us which phase and journey we are in
                std::vector<std::string> parsed_description;
                boost::split(parsed_description, std::get<0>(control_block_initial_guess.front()), boost::is_any_of(":"));
                std::string temp = parsed_description[0];
                std::string prefix = temp.substr(0, temp.find("_"));

                size_t journeyIndex = std::stoi(prefix.substr(1, prefix.find_first_of("p")));

                if (journeyIndex > this->number_of_journeys - 1)
                    break;

                size_t num_steps = this->options.Journeys[journeyIndex].override_num_steps
                    ? this->options.Journeys[journeyIndex].number_of_steps
                    : this->options.num_timesteps;

                if (number_of_guess_control_steps != num_steps)
                {
                    std::vector<std::tuple<std::string, double> > new_control_block_initial_guess;

                    //construct a data table of all control values
                    //we will do this by looping through the old vector and extracting the control values
                    //if this is a parallel shooting mission then we will also extract the state values - maybe we want to worry about that later, or even have the individual phases do it?
                    std::vector< std::pair<double, double> > Oldu_x;
                    std::vector< std::pair<double, double> > Oldu_y;
                    std::vector< std::pair<double, double> > Oldu_z;
                    std::vector< std::pair<double, double> > Oldu_command;

                    //first encode the left hand side of the phase
                    size_t oldControlBlockIndex = 0;
                    Oldu_x.push_back(std::make_pair(0.0, std::get<1>(control_block_initial_guess[oldControlBlockIndex])));
                    Oldu_y.push_back(std::make_pair(0.0, std::get<1>(control_block_initial_guess[oldControlBlockIndex + 1])));
                    Oldu_z.push_back(std::make_pair(0.0, std::get<1>(control_block_initial_guess[oldControlBlockIndex + 2])));

                    if (number_of_controls == 4)
                        Oldu_command.push_back(std::make_pair(0.0, std::get<1>(control_block_initial_guess[oldControlBlockIndex + 3])));

                    //then each step
                    for (int step = 0; step < number_of_guess_control_steps; ++step)
                    {
                        Oldu_x.push_back(std::make_pair((step + 0.5) / number_of_guess_control_steps, std::get<1>(control_block_initial_guess[oldControlBlockIndex++])));
                        Oldu_y.push_back(std::make_pair((step + 0.5) / number_of_guess_control_steps, std::get<1>(control_block_initial_guess[oldControlBlockIndex++])));
                        Oldu_z.push_back(std::make_pair((step + 0.5) / number_of_guess_control_steps, std::get<1>(control_block_initial_guess[oldControlBlockIndex++])));

                        if (number_of_controls == 4)
                            Oldu_command.push_back(std::make_pair((step + 0.5) / number_of_guess_control_steps, std::get<1>(control_block_initial_guess[oldControlBlockIndex++])));
                    }

                    //and finally the right-hand side
                    Oldu_x.push_back(std::make_pair(1.0, std::get<1>(control_block_initial_guess[oldControlBlockIndex - number_of_controls])));
                    Oldu_y.push_back(std::make_pair(1.0, std::get<1>(control_block_initial_guess[oldControlBlockIndex - number_of_controls + 1])));
                    Oldu_z.push_back(std::make_pair(1.0, std::get<1>(control_block_initial_guess[oldControlBlockIndex - number_of_controls + 2])));

                    if (number_of_controls == 4)
                        Oldu_command.push_back(std::make_pair(1.0, std::get<1>(control_block_initial_guess[-number_of_controls + 3])));

                    //interpolate the data tables to fill in the new control vector
                    math::interpolator ux_interpolator(Oldu_x);
                    math::interpolator uy_interpolator(Oldu_y);
                    math::interpolator uz_interpolator(Oldu_z);
                    math::interpolator ucommand_interpolator(Oldu_command);

                    for (int step = 0; step < num_steps; ++step)
                    {
                        double current_normalized_integration_variable = (step + 0.5) / num_steps;
                        std::stringstream stepstream;
                        stepstream << step;
                        //interpolate and encode the control vector
                        new_control_block_initial_guess.push_back({ prefix + ": step " + std::to_string(step) + " u_x", ux_interpolator.interpolate(current_normalized_integration_variable) });
                        new_control_block_initial_guess.push_back({ prefix + ": step " + std::to_string(step) + " u_y", uy_interpolator.interpolate(current_normalized_integration_variable) });
                        new_control_block_initial_guess.push_back({ prefix + ": step " + std::to_string(step) + " u_z", uz_interpolator.interpolate(current_normalized_integration_variable) });

                        if (number_of_controls == 4)
                            new_control_block_initial_guess.push_back({ prefix + ": step " + std::to_string(step) + " u_command", ucommand_interpolator.interpolate(current_normalized_integration_variable) });
                    }

                    control_block_initial_guess = new_control_block_initial_guess;
                }

                for (size_t controlIndex = 0; controlIndex < control_block_initial_guess.size(); ++controlIndex)
                {
                    processed_initial_guess.push_back(control_block_initial_guess[controlIndex]);
                }//end loop to to copy controls
            }//end control handling
            else
            {
                processed_initial_guess.push_back(this->options.trialX[Xindex++]);
            }
        }//end while (!reached_end_of_initial_guess)

         //Step 3: step through decision variables
        std::vector<double> newX(this->total_number_of_NLP_parameters);
        for (size_t Xindex2 = 0; Xindex2 < this->total_number_of_NLP_parameters; ++Xindex2)
        {
            bool initializedFlag = false;
            size_t journeyIndex;
            //Step 3.1: Does this decision variable appear in the initial guess?
            for (size_t trialXindex2 = 0; trialXindex2 < processed_initial_guess.size(); ++trialXindex2)
            {
                bool match = false;
                //grab the initial guess entry if it matches
                if (this->Xdescriptions[Xindex2] == std::get<0>(processed_initial_guess[trialXindex2]))
                    match = true;

                if (match)
                {
                    initializedFlag = true;

                    //Step 3.1.1
                    //Time variables need to be converted from days to seconds. Otherwise copy directly.
                    if (this->Xdescriptions[Xindex2].find("epoch") < 1024
                        || this->Xdescriptions[Xindex2].find("time") < 1024)
                    {
                        newX[Xindex2] = std::get<1>(processed_initial_guess[trialXindex2]) * 86400.0;
                    }
                    else if (this->Xdescriptions[Xindex2].find("Sundman independent variable") < 1024)
                    {
                        //we need the journey index so that we can scale it by that journey's LU
                        //prefix format is jJpP, so split the string by p and then remove the first character, convert to int
                        std::vector<std::string> descriptionComponents;

                        boost::split(descriptionComponents, this->Xdescriptions[Xindex2], boost::is_any_of("p"), boost::token_compress_on);

                        std::string journeyPrefix = descriptionComponents[0].substr(1, descriptionComponents[0].size() - 1);

                        journeyIndex = std::stoi(journeyPrefix);
                        
                        newX[Xindex2] = std::get<1>(processed_initial_guess[trialXindex2]) *  this->TheUniverse[journeyIndex].LU;
                    }                    
                    else
                        newX[Xindex2] = std::get<1>(processed_initial_guess[trialXindex2]);

                    //freeze
                    std::vector<std::string> DescriptionCell;
                    boost::split(DescriptionCell, this->Xdescriptions[Xindex2], boost::is_any_of("p"), boost::token_compress_on);
                    size_t thisJourney = std::stoi(boost::replace_all_copy(DescriptionCell[0], "j", ""));
                    if (this->options.Journeys[thisJourney].freeze_decision_variables)
                    {
                        this->Xlowerbounds[trialXindex2] = newX[Xindex2] - math::SMALL;
                        this->Xupperbounds[trialXindex2] = newX[Xindex2] + math::SMALL;
                    }
                }
                else if (this->Xdescriptions[Xindex2].find("virtual waypoint error") < 1024)
                {
                    initializedFlag = true;

                    newX[Xindex2] = 0.0;
                }
            }

            //Step 3.2: warn if you need to
            if ((newX[Xindex2] - this->Xupperbounds[Xindex2]) / this->X_scale_factors[Xindex2] > this->options.snopt_feasibility_tolerance)
            {

                if (this->Xdescriptions[Xindex2].find("epoch") < 1024
                    || this->Xdescriptions[Xindex2].find("time") < 1024)
                {
                    std::cout << "WARNING: X[" << Xindex2 << "]: " << this->Xdescriptions[Xindex2] << " violates the upper bound by " << (newX[Xindex2] - this->Xupperbounds[Xindex2]) / 86400.0 << ", which exceeds the feasibility tolerance" << std::endl;
                }
                else if (this->Xdescriptions[Xindex2].find("Sundman independent variable") < 1024)
                {
                    std::cout << "WARNING: X[" << Xindex2 << "]: " << this->Xdescriptions[Xindex2] << " violates the upper bound by " << (newX[Xindex2] - this->Xupperbounds[Xindex2]) / this->TheUniverse[journeyIndex].LU << ", which exceeds the feasibility tolerance" << std::endl;
                }
                else
                    std::cout << "WARNING: X[" << Xindex2 << "]: " << this->Xdescriptions[Xindex2] << " violates the upper bound by " << newX[Xindex2] - this->Xupperbounds[Xindex2] << ", which exceeds the feasibility tolerance" << std::endl;
            }
            else if ((this->Xlowerbounds[Xindex2] - newX[Xindex2]) / this->X_scale_factors[Xindex2] > this->options.snopt_feasibility_tolerance)
            {

                if (this->Xdescriptions[Xindex2].find("epoch") < 1024
                    || this->Xdescriptions[Xindex2].find("time") < 1024)
                {
                    std::cout << "WARNING: X[" << Xindex2 << "]: " << this->Xdescriptions[Xindex2] << " violates the lower bound by " << (newX[Xindex2] - this->Xlowerbounds[Xindex2]) / 86400.0 << ", which exceeds the feasibility tolerance" << std::endl;
                }
                else if (this->Xdescriptions[Xindex2].find("Sundman independent variable") < 1024)
                {
                    std::cout << "WARNING: X[" << Xindex2 << "]: " << this->Xdescriptions[Xindex2] << " violates the lower bound by " << (newX[Xindex2] - this->Xlowerbounds[Xindex2]) / this->TheUniverse[journeyIndex].LU << ", which exceeds the feasibility tolerance" << std::endl;
                }
                else
                    std::cout << "WARNING: X[" << Xindex2 << "]: " << this->Xdescriptions[Xindex2] << " violates the lower bound by " << newX[Xindex2] - this->Xlowerbounds[Xindex2] << ", which exceeds the feasibility tolerance" << std::endl; 
            }

            //Step 3.3: If it didn't appear, put in a random value between the bounds
            if (!initializedFlag)
			{
				std::cout<<"Initial guess missing value for: " << this->Xdescriptions[Xindex2]<<std::endl;
                newX[Xindex2] = DoubleDistribution(RNG) * (this->Xupperbounds[Xindex2] - this->Xlowerbounds[Xindex2]) + this->Xlowerbounds[Xindex2];
			}
		}

        //Step 4: pass the guess back
        return newX;
    }
    

    //******************************************sparsity pattern stuff
    //we keep these in Mission so that we can have a shorter call interface
    size_t
        Mission::create_sparsity_entry(const size_t& Findex,
            const size_t& Xstart,
            const bool& ForwardPass,
            const std::string& variable_name,
            size_t& sparsity_index_container)
    {
        return solver_utilities::create_sparsity_entry(Findex,
            Xstart,
            ForwardPass,
            variable_name,
            Fdescriptions,
            Xdescriptions,
            Gdescriptions,
            iGfun,
            jGvar,
            sparsity_index_container);
    }

    size_t
        Mission::create_sparsity_entry(const size_t& Findex,
            const size_t& Xstart,
            const bool& ForwardPass,
            const std::string& variable_name,
            std::vector<size_t>& sparsity_index_container)
    {
        return solver_utilities::create_sparsity_entry(Findex,
            Xstart,
            ForwardPass,
            variable_name,
            Fdescriptions,
            Xdescriptions,
            Gdescriptions,
            iGfun,
            jGvar,
            sparsity_index_container);
    }

    size_t
        Mission::create_sparsity_entry(const size_t& Findex,
            const size_t& Xindex,
            size_t& sparsity_index_container)
    {
        return solver_utilities::create_sparsity_entry(Findex,
            Xindex,
            Fdescriptions,
            Xdescriptions,
            Gdescriptions,
            iGfun,
            jGvar,
            sparsity_index_container);
    }

    size_t
        Mission::create_sparsity_entry(const size_t& Findex,
            const size_t& Xindex,
            std::vector<size_t>& sparsity_index_container)
    {
        return solver_utilities::create_sparsity_entry(Findex,
            Xindex,
            Fdescriptions,
            Xdescriptions,
            Gdescriptions,
            iGfun,
            jGvar,
            sparsity_index_container);
    }

    void
        Mission::create_sparsity_entry_vector(const size_t& Findex,
            const size_t& Xstart,
            const bool& ForwardPass,
            const int& number_of_entries,
            const std::string& variable_name,
            std::vector<size_t>& sparsity_index_container)
    {
        solver_utilities::create_sparsity_entry_vector(Findex,
            Xstart,
            ForwardPass,
            number_of_entries,
            variable_name,
            Fdescriptions,
            Xdescriptions,
            Gdescriptions,
            iGfun,
            jGvar,
            sparsity_index_container);
    }
} /* namespace EMTG */
