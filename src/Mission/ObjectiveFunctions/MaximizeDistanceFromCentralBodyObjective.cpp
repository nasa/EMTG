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

#include "MaximizeDistanceFromCentralBodyObjective.h"

namespace EMTG
{
    namespace ObjectiveFunctions
{

        MaximizeDistanceFromCentralBodyObjective::MaximizeDistanceFromCentralBodyObjective(Astrodynamics::universe* Universe,
            HardwareModels::Spacecraft* mySpacecraft,
            missionoptions* myOptions,
            Mission* myMission) :
            ObjectiveFunctionBase(Universe,
                mySpacecraft,
                myOptions,
                myMission)
        {
            this->name = "MaximizeFinalPeriapsisObjectiveFunction";

            size_t number_of_phases_in_final_journey = this->myMission->getJourney(this->myOptions->number_of_journeys - 1)->getNumberOfPhases();

            this->myArrivalEvent = this->myMission->getJourney(this->myOptions->objective_journey - 1)->getPhase(number_of_phases_in_final_journey - 1)->getArrivalEvent();

            this->ObjectiveScale = 1.0;
                
			this->Gindices_objectivefunction_decision_variables_that_affect_state.resize(3);
			this->Gindices_objectivefunction_time_decision_variables_that_affect_state.resize(3);
        }//end constructor

        void MaximizeDistanceFromCentralBodyObjective::calcbounds()
        {
            //this objective function has derivatives with respect to ALL variables that affect final arrival event periapsis
            std::vector< std::tuple<size_t, size_t, double> > Derivatives_of_StateAfterArrival = this->myArrivalEvent->get_Derivatives_of_StateAfterEvent();
            std::vector< std::tuple<size_t, size_t, double> > Derivatives_of_StateAfterArrival_wrt_Time = this->myArrivalEvent->get_Derivatives_of_StateAfterEvent_wrt_Time();

            //non-time variables
			for (size_t stateIndex : {0, 1, 2})
			{
				std::vector<size_t> dIndex_position_thisState;
				for (size_t dIndex = 0; dIndex < Derivatives_of_StateAfterArrival.size(); ++dIndex)
				{
					if (std::get<1>(Derivatives_of_StateAfterArrival[dIndex]) == stateIndex)
					{
                        dIndex_position_thisState.push_back(dIndex);
						this->Derivatives_of_ObjectiveFunction.push_back({ std::get<0>(Derivatives_of_StateAfterArrival[dIndex]), 0.0 });
						this->create_sparsity_entry(0,
							std::get<0>(Derivatives_of_StateAfterArrival[dIndex]),
							this->Gindices_objectivefunction_decision_variables_that_affect_state[stateIndex]);
					}
				}
				this->dIndex_position.push_back(dIndex_position_thisState);
			}

            //time variables
			for (size_t stateIndex : {0, 1, 2})
			{
				std::vector<size_t> dIndex_position_thisState;
				for (size_t dIndex = 0; dIndex < Derivatives_of_StateAfterArrival_wrt_Time.size(); ++dIndex)
				{
					if (std::get<1>(Derivatives_of_StateAfterArrival_wrt_Time[dIndex]) == stateIndex)
					{
                        dIndex_position_thisState.push_back(dIndex);
						this->Derivatives_of_ObjectiveFunction.push_back({ std::get<0>(Derivatives_of_StateAfterArrival_wrt_Time[dIndex]), 0.0 });
						this->create_sparsity_entry(0,
							std::get<0>(Derivatives_of_StateAfterArrival_wrt_Time[dIndex]),
							this->Gindices_objectivefunction_time_decision_variables_that_affect_state[stateIndex]);
					}
				}
				this->dIndex_position_wrt_Time.push_back(dIndex_position_thisState);
			}
        }//end calcbounds()

        void MaximizeDistanceFromCentralBodyObjective::process(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            math::Matrix<doubleType> FinalState = this->myArrivalEvent->get_state_after_event();
			doubleType radiusOfOrbit = sqrt(FinalState(0) * FinalState(0) + FinalState(1) * FinalState(1) + FinalState(2) * FinalState(2));

            this->ObjectiveValue = -this->ObjectiveScale * radiusOfOrbit / this->myUniverse->LU;
            F[0] = this->ObjectiveValue;

			if (needG)
            {
                std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateAfterArrival = this->myArrivalEvent->get_Derivatives_of_StateAfterEvent();
                std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateAfterArrival_wrt_Time = this->myArrivalEvent->get_Derivatives_of_StateAfterEvent_wrt_Time();

                //first reset the time components
                for (size_t Gindex : this->Gindices_objectivefunction_time_decision_variables_that_affect_state[0])
                    G[Gindex] = 0.0;

                for (size_t Gindex : this->Gindices_objectivefunction_decision_variables_that_affect_state[0])
                    G[Gindex] = 0.0;

                for (size_t stateIndex : {0, 1, 2})
                {
                    double dRadius_dthisState = (FinalState(stateIndex) / radiusOfOrbit) _GETVALUE;

                    //non-time partials
                    for (size_t ddIndex = 0; ddIndex < this->Gindices_objectivefunction_decision_variables_that_affect_state[stateIndex].size(); ++ddIndex)
                    {                        
                        size_t Gindex = this->Gindices_objectivefunction_decision_variables_that_affect_state[stateIndex][ddIndex];
                        size_t dIndex = this->dIndex_position[stateIndex][ddIndex];
                        size_t Xindex = std::get<0>(Derivatives_of_StateAfterArrival[dIndex]);

                        G[Gindex] += -ObjectiveScale / this->myUniverse->LU
                            * this->X_scale_factors->operator[](Xindex)
                            * dRadius_dthisState * std::get<2>(Derivatives_of_StateAfterArrival[dIndex]);
                    }

                    //time partials
                    for (size_t ddIndex = 0; ddIndex < this->Gindices_objectivefunction_time_decision_variables_that_affect_state[stateIndex].size(); ++ddIndex)
                    {
                        size_t Gindex = this->Gindices_objectivefunction_time_decision_variables_that_affect_state[stateIndex][ddIndex];
                        size_t dIndex = this->dIndex_position_wrt_Time[stateIndex][ddIndex];
                        size_t Xindex = std::get<0>(Derivatives_of_StateAfterArrival_wrt_Time[dIndex]);

                        G[Gindex] += -ObjectiveScale / this->myUniverse->LU
                            * this->X_scale_factors->operator[](Xindex)
                            * dRadius_dthisState * std::get<2>(Derivatives_of_StateAfterArrival_wrt_Time[dIndex]);
                    }
                }//end loop over state
            }//end derivatives
        }//end process()

        void MaximizeDistanceFromCentralBodyObjective::output(std::ofstream& outputfile)
        {
            outputfile << "Objective function is \"maximize final periapsis\"" << std::endl;
            outputfile << "J = " << this->ObjectiveValue _GETVALUE << std::endl;
        }
    }//end namespace ObjectiveFunctions
}//end namespace EMTG