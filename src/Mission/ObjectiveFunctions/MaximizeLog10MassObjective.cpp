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

#include "MaximizeLog10MassObjective.h"

namespace EMTG
{
    namespace ObjectiveFunctions
{

        MaximizeLog10MassObjective::MaximizeLog10MassObjective(Astrodynamics::universe* Universe,
            HardwareModels::Spacecraft* mySpacecraft,
            missionoptions* myOptions,
            Mission* myMission) :
            ObjectiveFunctionBase(Universe,
                mySpacecraft,
                myOptions,
                myMission)
        {
            this->name = "MaximizeLog10MassObjectiveFunction";

            size_t number_of_phases_in_final_journey = this->myMission->getJourney(this->myOptions->number_of_journeys - 1)->getNumberOfPhases();

            this->myArrivalEvent = this->myMission->getJourney(this->myOptions->number_of_journeys - 1)->getPhase(number_of_phases_in_final_journey - 1)->getArrivalEvent();
                
        }//end constructor

        void MaximizeLog10MassObjective::calcbounds()
        {
            //this objective function has derivatives with respect to ALL variables that affect final arrival event mass
            std::vector< std::tuple<size_t, size_t, double> > Derivatives_of_StateAfterArrival = this->myArrivalEvent->get_Derivatives_of_StateAfterEvent();
            std::vector< std::tuple<size_t, size_t, double> > Derivatives_of_StateAfterArrival_wrt_Time = this->myArrivalEvent->get_Derivatives_of_StateAfterEvent_wrt_Time();

            //non-time variables
            for (size_t dIndex = 0; dIndex < Derivatives_of_StateAfterArrival.size(); ++dIndex)
            {
                if (std::get<1>(Derivatives_of_StateAfterArrival[dIndex]) == 6)
                {
                    this->dIndex_final_mass.push_back(dIndex);
                    this->Derivatives_of_ObjectiveFunction.push_back({ std::get<0>(Derivatives_of_StateAfterArrival[dIndex]), std::get<2>(Derivatives_of_StateAfterArrival[dIndex]) });
                    this->create_sparsity_entry(0,
                                                std::get<0>(Derivatives_of_StateAfterArrival[dIndex]),
                                                this->Gindex_derivatives_of_objective_function);
                }
            }
            this->number_of_non_time_derivatives = this->Derivatives_of_ObjectiveFunction.size();

            //time variables
            this->number_of_derivatives_wrt_time_variables = 0;
            for (size_t dIndex = 0; dIndex < Derivatives_of_StateAfterArrival_wrt_Time.size(); ++dIndex)
            {
                if (std::get<1>(Derivatives_of_StateAfterArrival_wrt_Time[dIndex]) == 6)
                {
                    this->dIndex_final_mass_wrt_Time.push_back(dIndex);
                    this->Derivatives_of_ObjectiveFunction.push_back({ std::get<0>(Derivatives_of_StateAfterArrival_wrt_Time[dIndex]), std::get<2>(Derivatives_of_StateAfterArrival_wrt_Time[dIndex]) });
                    this->create_sparsity_entry(0,
                                                std::get<0>(Derivatives_of_StateAfterArrival_wrt_Time[dIndex]),
                                                this->Gindex_derivatives_of_objective_function);
                    ++this->number_of_derivatives_wrt_time_variables;
                }
            }

        }//end calcbounds()

        void MaximizeLog10MassObjective::process(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            math::Matrix<doubleType> FinalState = this->myArrivalEvent->get_state_after_event();
            this->ObjectiveValue = -log10((FinalState(6) + 1.0) * this->myOptions->Journeys.back().maximum_mass);
            F[0] = this->ObjectiveValue;

            if (needG)
            {
                std::vector< std::tuple<size_t, size_t, double> > Derivatives_of_StateAfterArrival = this->myArrivalEvent->get_Derivatives_of_StateAfterEvent();
                std::vector< std::tuple<size_t, size_t, double> > Derivatives_of_StateAfterArrival_wrt_Time = this->myArrivalEvent->get_Derivatives_of_StateAfterEvent_wrt_Time();

                //non-time variables
                for (size_t dIndex = 0; dIndex < this->number_of_non_time_derivatives; ++dIndex)
                {
                    std::get<1>(this->Derivatives_of_ObjectiveFunction[dIndex]) = std::get<2>(Derivatives_of_StateAfterArrival[this->dIndex_final_mass[dIndex]]);

                    size_t Gindex = this->Gindex_derivatives_of_objective_function[dIndex];
                    size_t Xindex = this->jGvar->operator[](Gindex);

                    G[Gindex] = this->X_scale_factors->operator[](Xindex)
                        * -1.0 / (FinalState(6) + 1.0)_GETVALUE / log(10.0)
                        * std::get<1>(this->Derivatives_of_ObjectiveFunction[dIndex]);
                }

                //time variables
                for (size_t dtIndex = 0; dtIndex < this->number_of_derivatives_wrt_time_variables; ++dtIndex)
                {
                    size_t dIndex = dtIndex + this->number_of_non_time_derivatives;

                    std::get<1>(this->Derivatives_of_ObjectiveFunction[dIndex]) = std::get<2>(Derivatives_of_StateAfterArrival_wrt_Time[this->dIndex_final_mass_wrt_Time[dtIndex]]);

                    size_t Gindex = this->Gindex_derivatives_of_objective_function[dIndex];
                    size_t Xindex = this->jGvar->operator[](Gindex);

                    G[Gindex] = this->X_scale_factors->operator[](Xindex)
                        * -1.0 / (FinalState(6) + 1.0)_GETVALUE / log(10.0)
                        * std::get<1>(this->Derivatives_of_ObjectiveFunction[dIndex]);
                }
            }
        }//end process()

        void MaximizeLog10MassObjective::output(std::ofstream& outputfile)
        {
            outputfile << "Objective function is \"maximize log10(final mass)\"" << std::endl;
            outputfile << "J = " << this->ObjectiveValue _GETVALUE << std::endl;
        }
    }//end namespace ObjectiveFunctions
}//end namespace EMTG