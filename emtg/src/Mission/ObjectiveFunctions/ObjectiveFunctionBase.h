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

#pragma once

#include "doubleType.h"

#include <string>

#include "universe.h"
#include "missionoptions.h"
#include "Spacecraft.h"
#include "mission.h"

namespace EMTG
{
    //forward declaration of Mission (pointer to parent)
    class Mission;

    namespace ObjectiveFunctions
    {
        class ObjectiveFunctionBase : virtual public sparsey_thing
        {
        public:
            //constructors
            ObjectiveFunctionBase() {};

            ObjectiveFunctionBase(Astrodynamics::universe* Universe,
                HardwareModels::Spacecraft* mySpacecraft,
                missionoptions* myOptions,
                Mission* myMission);

            //destructor
            virtual ~ObjectiveFunctionBase() {};


            //public methods
            doubleType getObjectiveValue() const { return this->ObjectiveValue; }

            virtual void output(std::ofstream& outputfile) = 0;

            //calcbounds
            void setup_calcbounds(
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
                std::vector<double>* A);

            //calcbounds goes in the specialized event
            virtual void calcbounds() = 0;

            //process
            virtual void process(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG) = 0;
				
			doubleType getUnscaledObjective() {return UnscaledObjectiveValue;};

        protected:

            //fields
            std::string name;
            Astrodynamics::universe* myUniverse;
            HardwareModels::Spacecraft* mySpacecraft;
            missionoptions* myOptions;            

            //calcbounds things
            std::string prefix;

            //objective value
            doubleType ObjectiveValue;
            doubleType UnscaledObjectiveValue;
            double ObjectiveScale;

            //pointer to mission
            Mission* myMission;

            //derivative containers
            std::vector< std::tuple<size_t, double> > Derivatives_of_ObjectiveFunction;//Xindex, derivative value
            std::vector<size_t> Gindex_derivatives_of_objective_function;
                
        };
    }//end namespace ObjectiveFunctions
}//end namespace EMTG