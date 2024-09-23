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

//forward MGAnDSMs maneuver epoch constraint with respect to next event
//9-26-2017

#include "MGAnDSMs_forward_maneuver_epoch_constraint_wrt_next_event.h"

#include "EMTG_solver_utilities.h"

#include "MGAnDSMs_subphase.h"

#include "boost/algorithm/string/split.hpp"

namespace EMTG
{
    namespace Phases
    {
        MGAnDSMs_forward_maneuver_epoch_constraint_wrt_next_event::MGAnDSMs_forward_maneuver_epoch_constraint_wrt_next_event(const std::string& name,
            const size_t& journeyIndex,
            const size_t& phaseIndex,
            const size_t& subphaseIndex,
            const size_t& stageIndex,
            MGAnDSMs_subphase* mySubPhase,
            Astrodynamics::universe* Universe,
            HardwareModels::Spacecraft* mySpacecraft,
            missionoptions* myOptions,
            const std::string& ConstraintDefinition) :
            MGAnDSMs_maneuver_epoch_constraint::MGAnDSMs_maneuver_epoch_constraint(name,
                journeyIndex,
                phaseIndex,
                subphaseIndex,
                stageIndex,
                mySubPhase,
                Universe,
                mySpacecraft,
                myOptions,
                ConstraintDefinition)
        {
        }//end constructor
        //************************************calcbounds methods
        void MGAnDSMs_forward_maneuver_epoch_constraint_wrt_next_event::calcbounds()
        {
            //Step 1: parse the constraint definition
            std::vector<std::string> ConstraintDefinitionCell;
            boost::split(ConstraintDefinitionCell, ConstraintDefinition, boost::is_any_of("_"), boost::token_compress_on);

            //Step 2: create the constraint
            this->lowerBound = std::stod(ConstraintDefinitionCell[3]) / 100.0;
            this->upperBound = std::stod(ConstraintDefinitionCell[4]) / 100.0;
            this->Flowerbounds->push_back(this->lowerBound - this->upperBound);
            this->Fupperbounds->push_back(0.0);
            this->Fdescriptions->push_back(prefix + "maneuver epoch relative to next event constraint");

            //sparsity pattern - derivatives with respect to phase flight time
            this->Xindex_PhaseFlightTime = this->mySubPhase->get_timeVariables().back();

            this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                this->Xindex_PhaseFlightTime,
                this->Gindex_wrt_TimeVariables);

            //derivatives with respect to burn index - only the NEXT EVENT's burn index because that is the measurement between now and the next event
            //there are two cases:
            //1: there is only one maneuver in this phase, so there is only one burn index and the derivative is (1.0 - eta)
            //2: there is another maneuver, either forward or backward, and we can find its burn index in the decision vector
            //if (this->myJourneyOptions->impulses_per_phase[this->phaseIndex] == 1)
            //{
            //    this->create_sparsity_entry(this->Fdescriptions->size() - 1,
            //        this->mySubPhase->getXindex_burnIndex().back(),
            //        this->Gindex_wrt_BurnIndices);
            //}//end case of single-maneuver subphase
            //else
            {
                if (this->mySubPhase->getBordersMatchPoint())
                {
                    //the next event is the last backward event, so it will be the last burn index in the decision vector at this point
                    for (size_t Xindex = this->Xdescriptions->size() - 1; Xindex > 0; --Xindex)
                    {
                        if (Xdescriptions->at(Xindex).find("burn index") < 1024)
                        {
                            this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                                Xindex,
                                this->Gindex_wrt_BurnIndices);
                            break;
                        }
                    }
                }
                else
                {
                    //the next event is the next forward event
                    bool skip = true;
                    for (size_t Xindex = this->mySubPhase->getFirst_X_entry_this_subphase(); Xindex < this->Xdescriptions->size(); ++Xindex)
                    {
                        if (Xdescriptions->at(Xindex).find("burn index") < 1024)
                        {
                            if (skip)
                                skip = false;
                            else
                            {
                                this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                                    Xindex,
                                    this->Gindex_wrt_BurnIndices);
                                break;
                            }
                        }
                    }
                }
            }//end case of looking for the next event
        }//end calcbounds()

         //******************************************process methods
        void MGAnDSMs_forward_maneuver_epoch_constraint_wrt_next_event::process_constraint(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {            
            //Step 1: evaluate the constraint
            size_t Gindex_BurnIndex = this->Gindex_wrt_BurnIndices.back();
            size_t Xindex_BurnIndex = this->jGvar->operator[](Gindex_BurnIndex);
            size_t Gindex_PhaseFlightTime = this->Gindex_wrt_TimeVariables.back();
            size_t Xindex_PhaseFlightTime = this->jGvar->operator[](Gindex_PhaseFlightTime);

            doubleType PhaseFlightTime = X[Xindex_PhaseFlightTime];
            doubleType burnIndex = X[Xindex_BurnIndex];

            doubleType NextSubPhaseTime = PhaseFlightTime * burnIndex;
            F[Findex++] = NextSubPhaseTime / 86400.0 / 100.0 - this->upperBound;

            //Step 2: derivatives
            if (needG)
            {
                //with respect to burn index                        
                G[Gindex_BurnIndex] = this->X_scale_factors->operator[](Xindex_BurnIndex)
                    * PhaseFlightTime _GETVALUE
                    / 86400.0 / 100.0;

                //phase flight time
                G[Gindex_PhaseFlightTime] = this->X_scale_factors->operator[](Xindex_PhaseFlightTime)
                    * burnIndex _GETVALUE
                    / 86400.0 / 100.0;
            }//end derivatives
        }//end process_constraint()
    }//close namespace Phases
}//close namespace EMTG