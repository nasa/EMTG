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

//EphemerisPeggedSpiralArrival
//Jacob Englander 10-24-2017

#include "EdelbaumSpiral.h"
#include "EphemerisPeggedArrival.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        class EphemerisPeggedSpiralArrival : virtual public EphemerisPeggedArrival
        {
        public:
            EphemerisPeggedSpiralArrival(const std::string& name,
                const size_t& journeyIndex,
                const size_t& phaseIndex,
                size_t& stageIndex,
                Astrodynamics::universe* Universe,
                HardwareModels::Spacecraft* mySpacecraft,
                missionoptions* myOptions);

            ~EphemerisPeggedSpiralArrival(); 
            
            //we have to override setup-calcbounds so that we can set up the owned EdelbaumSpiral class and its owned classes
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

            void calcbounds(std::vector<size_t> timeVariables);

            void process_event(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            void output(std::ofstream& outputfile,
                const double& launchdate,
                size_t& eventcount);

        protected:
            void calcbounds_event_main();
            void calcbounds_event_right_side();
            void calcbounds_virtual_propellant_constraints() {}; //stub, propellant constraints live in the owned spiral segment objects

            void process_event_main(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            void process_event_right_side(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            void process_virtual_propellant_constraints(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG) {}; //stub, propellant constraints live in the owned spiral segment objects

            //fields
            EdelbaumSpiral* mySpiral;
            math::Matrix<doubleType> StateAfterSpiral;
            
            //derivative indices
            size_t Xindex_SpiralFlightTime;
            size_t Gindex_dSpiralFlightTimeConstraint_dSpiralFlightTime;
            std::vector<size_t> Gindex_dSpiralFlightTimeConstraint_dPreviousTimeVariable; //first entry is launch epoch
            size_t Gindex_dSpiralFlightTimeConstraint_dMassAtBeginningOfSpiral;

            size_t dIndex_EpochAfterSpiral_wrt_SpiralFlightTime;
            std::vector<size_t> dIndex_MassAfterSprial_wrt_PreviousTimeVariables;
            std::vector<size_t> dIndex_6state_wrt_SpiralFlightTime;

            std::vector<size_t> Gindex_VirtualElectricPropellantConstraint_wrt_PreviousTimeVariables;
            size_t Gindex_VirtualElectricPropellantConstraint_wrt_MassBeforeSpiral;
        };//end class EdelbaumSpiralDeparture
    }//close namespace BoundaryEvents
}//close namespace EMTG