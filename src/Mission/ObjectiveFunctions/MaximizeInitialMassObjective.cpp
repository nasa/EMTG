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

#include "MaximizeInitialMassObjective.h"

namespace EMTG
{
    namespace ObjectiveFunctions
{

        MaximizeInitialMassObjective::MaximizeInitialMassObjective(Astrodynamics::universe* Universe,
            HardwareModels::Spacecraft* mySpacecraft,
            missionoptions* myOptions,
            Mission* myMission) :
            ObjectiveFunctionBase(Universe,
                mySpacecraft,
                myOptions,
                myMission)
        {
            this->name = "MaximizeInitialMassObjectiveFunction";

            this->myDepartureEvent = this->myMission->getJourney(0)->getFirstPhase()->getDepartureEvent();

            this->ObjectiveScale = 1.0;
                
        }//end constructor

        void MaximizeInitialMassObjective::calcbounds()
        {
            //this objective function has derivatives with respect to ALL variables that affect final Departure event mass
            std::vector< std::tuple<size_t, size_t, double> > Derivatives_of_StateAfterDeparture = this->myDepartureEvent->get_Derivatives_of_StateAfterEvent();
            std::vector< std::tuple<size_t, size_t, double> > Derivatives_of_StateAfterDeparture_wrt_Time = this->myDepartureEvent->get_Derivatives_of_StateAfterEvent_wrt_Time();

            //non-time variables
            for (size_t dIndex = 0; dIndex < Derivatives_of_StateAfterDeparture.size(); ++dIndex)
            {
                if (std::get<1>(Derivatives_of_StateAfterDeparture[dIndex]) == 6)
                {
                    this->dIndex_initial_mass.push_back(dIndex);
                    this->Derivatives_of_ObjectiveFunction.push_back({ std::get<0>(Derivatives_of_StateAfterDeparture[dIndex]), std::get<2>(Derivatives_of_StateAfterDeparture[dIndex]) });
                    this->create_sparsity_entry(0,
                                                std::get<0>(Derivatives_of_StateAfterDeparture[dIndex]),
                                                this->Gindex_derivatives_of_objective_function);
                }
            }
            this->number_of_non_time_derivatives = this->Derivatives_of_ObjectiveFunction.size();

            //time variables
            this->number_of_derivatives_wrt_time_variables = 0;
            for (size_t dIndex = 0; dIndex < Derivatives_of_StateAfterDeparture_wrt_Time.size(); ++dIndex)
            {
                if (std::get<1>(Derivatives_of_StateAfterDeparture_wrt_Time[dIndex]) == 6)
                {
                    this->dIndex_initial_mass_wrt_Time.push_back(dIndex);
                    this->Derivatives_of_ObjectiveFunction.push_back({ std::get<0>(Derivatives_of_StateAfterDeparture_wrt_Time[dIndex]), std::get<2>(Derivatives_of_StateAfterDeparture_wrt_Time[dIndex]) });
                    this->create_sparsity_entry(0,
                                                std::get<0>(Derivatives_of_StateAfterDeparture_wrt_Time[dIndex]),
                                                this->Gindex_derivatives_of_objective_function);
                    ++this->number_of_derivatives_wrt_time_variables;
                }
            }

        }//end calcbounds()

        void MaximizeInitialMassObjective::process(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            math::Matrix<doubleType> InitialState = this->myDepartureEvent->get_state_after_event();
            this->ObjectiveValue = -this->ObjectiveScale * InitialState(6) * this->myUniverse->continuity_constraint_scale_factors(6);
            F[0] = this->ObjectiveValue;

            if (needG)
            {
                std::vector< std::tuple<size_t, size_t, double> > Derivatives_of_StateAfterDeparture = this->myDepartureEvent->get_Derivatives_of_StateAfterEvent();
                std::vector< std::tuple<size_t, size_t, double> > Derivatives_of_StateAfterDeparture_wrt_Time = this->myDepartureEvent->get_Derivatives_of_StateAfterEvent_wrt_Time();

                //non-time variables
                for (size_t dIndex = 0; dIndex < this->number_of_non_time_derivatives; ++dIndex)
                {
                    std::get<1>(this->Derivatives_of_ObjectiveFunction[dIndex]) = std::get<2>(Derivatives_of_StateAfterDeparture[this->dIndex_initial_mass[dIndex]]);

                    size_t Gindex = this->Gindex_derivatives_of_objective_function[dIndex];
                    size_t Xindex = this->jGvar->operator[](Gindex);

                    G[Gindex] = this->ObjectiveScale * this->X_scale_factors->operator[](Xindex)
                        * -std::get<1>(this->Derivatives_of_ObjectiveFunction[dIndex])
                        * this->myUniverse->continuity_constraint_scale_factors(6);
                }

                //time variables
                for (size_t dtIndex = 0; dtIndex < this->number_of_derivatives_wrt_time_variables; ++dtIndex)
                {
                    size_t dIndex = dtIndex + this->number_of_non_time_derivatives;

                    std::get<1>(this->Derivatives_of_ObjectiveFunction[dIndex]) = std::get<2>(Derivatives_of_StateAfterDeparture_wrt_Time[this->dIndex_initial_mass_wrt_Time[dtIndex]]);

                    size_t Gindex = this->Gindex_derivatives_of_objective_function[dIndex];
                    size_t Xindex = this->jGvar->operator[](Gindex);

                    G[Gindex] = this->ObjectiveScale * this->X_scale_factors->operator[](Xindex)
                        * -std::get<1>(this->Derivatives_of_ObjectiveFunction[dIndex])
                        * this->myUniverse->continuity_constraint_scale_factors(6);
                }
            }
        }//end process()

        void MaximizeInitialMassObjective::output(std::ofstream& outputfile)
        {
            outputfile << "Objective function is \"maximize initial mass\"" << std::endl;
            outputfile << "J = " << this->ObjectiveValue _GETVALUE << std::endl;
        }
    }//end namespace ObjectiveFunctions
}//end namespace EMTG