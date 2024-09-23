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

#include "RBP.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        namespace SpecializedConstraints
        {

            RBPconstraint::RBPconstraint(const std::string& name,
                const size_t& journeyIndex,
                const size_t& phaseIndex,
                const size_t& stageIndex,
                Astrodynamics::universe* Universe,
                HardwareModels::Spacecraft* mySpacecraft,
                missionoptions* myOptions,
                BoundaryEventBase* myBoundaryEvent,
                const std::string& constraintDefinition) :
                SpecializedBoundaryConstraintBase::SpecializedBoundaryConstraintBase(name,
                    journeyIndex,
                    phaseIndex,
                    stageIndex,
                    Universe,
                    mySpacecraft,
                    myOptions,
                    myBoundaryEvent,
                    constraintDefinition)
            {
                this->R_body_reference = math::Matrix<doubleType>(3, 1, 0.0);
                this->V_body_probe = math::Matrix<doubleType>(3, 1, 0.0);
                this->dcosRBP_dR_body_reference = math::Matrix<double>(3, 1, 0.0);
                this->dcosRBP_dV_body_probe = math::Matrix<double>(3, 1, 0.0);
                this->dR_body_probe_dt = math::Matrix<double>(3, 1, 0.0);
                this->dR_reference_cb_dt = math::Matrix<double>(3, 1, 0.0);

                //derived class version of myBoundaryEvent has a different type than SpecializedBoundaryConstraintBase uses
                this->myBoundaryEvent = myBoundaryEvent;
            }

            void RBPconstraint::calcbounds()
            {
                // Step 1: parse the constraint definition
				// constraintDefinition is like p0_arrival_RBP_cb_0.0_10.0 (numbers are in degrees)
                std::vector<std::string> ConstraintDefinitionCell;
                std::string splitty_thing = boost::algorithm::to_lower_copy(this->constraintDefinition);
                boost::split(ConstraintDefinitionCell, splitty_thing, boost::is_any_of("_"), boost::token_compress_on);

                //throw an error if the user is trying to apply this constraint to the wrong kind of boundary event
                if (ConstraintDefinitionCell[1] != "arrival")
                    throw std::invalid_argument("RBP constraint may only be applied to boundary events of class EphemerisPeggedArrivalWithVinfinity. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));

                //Step 2: which body are we doing this with respect to?
                if (ConstraintDefinitionCell[3] == "cb")
                {
                    this->isCentralBodyConstraint = true;
                    this->referenceName = this->myUniverse->central_body_name;
                }
                else
                {
                    this->isCentralBodyConstraint = false;

                    int bodyIndex = std::stoi(ConstraintDefinitionCell[3]) - 1;
                    this->myReference = &this->myUniverse->bodies[bodyIndex];
                    this->referenceName = this->myReference->name;
                }

                //Step 4: create the constraint
                //figure out the lower and upper bounds - convert from degrees to radians
                //the bigger number has the smaller cosine
                this->Flowerbounds->push_back(cos(std::stod(ConstraintDefinitionCell[5]) * math::PI / 180.0));
                
                this->Fupperbounds->push_back(cos(std::stod(ConstraintDefinitionCell[4]) * math::PI / 180.0));
                
                this->Fdescriptions->push_back(prefix + this->referenceName + "-Body-Probe constraint");

                //Step 5: sparsity pattern
                //derivatives with respect to anything influencing the boundary event's right-hand position and velocity vectors
                std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateAfterEvent = this->myBoundaryEvent->get_Derivatives_of_StateAfterEvent();
                std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateAfterEvent_wrt_Time = this->myBoundaryEvent->get_Derivatives_of_StateAfterEvent_wrt_Time();
                std::vector<size_t> dIndex_VbeforeEvent_dVinfinity_in = this->myBoundaryEvent->getdIndex_VbeforeEvent_dVinfinity_in();

                for (size_t stateIndex : {0, 1, 2})
                {
                    //non-time variables
                    std::vector<size_t> state_dIndex_with_respect_to_PositionAfterEvent;
                    std::vector<size_t> state_Gindex_constraint_wrt_PositionAfterEvent_variables;
                    for (size_t dIndex = 0; dIndex < Derivatives_of_StateAfterEvent.size(); ++dIndex)
                    {
                        if (std::get<1>(Derivatives_of_StateAfterEvent[dIndex]) == stateIndex)
                        {
                            state_dIndex_with_respect_to_PositionAfterEvent.push_back(dIndex);

                            this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                                std::get<0>(Derivatives_of_StateAfterEvent[dIndex]),
                                state_Gindex_constraint_wrt_PositionAfterEvent_variables);
                        }
                    }
                    this->dIndex_with_respect_to_PositionAfterEvent.push_back(state_dIndex_with_respect_to_PositionAfterEvent);
                    this->Gindex_constraint_wrt_PositionAfterEvent_variables.push_back(state_Gindex_constraint_wrt_PositionAfterEvent_variables);

                    //time variables
                    std::vector<size_t> state_dIndex_with_respect_to_PositionAfterEvent_wrt_Time;
                    std::vector<size_t> state_Gindex_constraint_wrt_PositionAfterEvent_time_variables;
                    for (size_t dIndex = 0; dIndex < Derivatives_of_StateAfterEvent_wrt_Time.size(); ++dIndex)
                    {
                        if (std::get<1>(Derivatives_of_StateAfterEvent_wrt_Time[dIndex]) == stateIndex)
                        {
                            state_dIndex_with_respect_to_PositionAfterEvent_wrt_Time.push_back(dIndex);

                            this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                                std::get<0>(Derivatives_of_StateAfterEvent_wrt_Time[dIndex]),
                                state_Gindex_constraint_wrt_PositionAfterEvent_time_variables);
                        }
                    }
                    this->dIndex_with_respect_to_PositionAfterEvent_wrt_Time.push_back(state_dIndex_with_respect_to_PositionAfterEvent_wrt_Time);
                    this->Gindex_constraint_wrt_PositionAfterEvent_time_variables.push_back(state_Gindex_constraint_wrt_PositionAfterEvent_time_variables);

                }

                //v-infinity
                for (size_t dIndex : dIndex_VbeforeEvent_dVinfinity_in)
                {
                    this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                        std::get<0>(Derivatives_of_StateAfterEvent[dIndex]),
                        Gindex_constraint_wrt_Vinfinity);
                }

                //derivatives with respect to time variables that affect reference body position
                if (!this->isCentralBodyConstraint)
                {
                    std::vector<size_t> timeVariables = this->myBoundaryEvent->get_Xindices_EventRightEpoch();
                    for (size_t Xindex : timeVariables)
                    {
                        this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                            Xindex,
                            this->Gindex_constraint_wrt_time_variables);
                    }
                }
            }//end calcbounds()

            void RBPconstraint::process_constraint(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG)
            {
                
                //Step 1: get the state relative to the central body
                math::Matrix<doubleType>& StateRelativeToCentralBody = this->myBoundaryEvent->get_state_after_event();

                //Step 2: get the vector from the body to the reference
                if (this->isCentralBodyConstraint) //is the reference body the central body?
                {
                    for (size_t stateIndex : {0, 1, 2})
                        this->R_body_reference(stateIndex) = -StateRelativeToCentralBody(stateIndex);
                }
                else //reference body is some other body in the Universe
                {
                    //Step 2.1: find the reference body relative to the central body
                    doubleType temp_body_state[12];
                    this->myReference->locate_body(StateRelativeToCentralBody(7),
                        temp_body_state,
                        needG,
                        *this->myOptions);

                    //Step 2.2: compute the vector from the body to the reference body
                    for (size_t stateIndex : {0, 1, 2})
                    {
                        this->R_body_reference(stateIndex) = -StateRelativeToCentralBody(stateIndex) + temp_body_state[stateIndex];
                        this->dR_reference_cb_dt(stateIndex) = temp_body_state[stateIndex + 6] _GETVALUE;
                    }
                }

                //Step 3: get the vector from the body to the probe (move backward one second along the v-infinity vector direction)
                math::Matrix<doubleType> Vinfinity = this->myBoundaryEvent->getVinfinityIn();

                for (size_t stateIndex : {0, 1, 2})
                {
                    this->V_body_probe(stateIndex) = -Vinfinity(stateIndex); //NEGATIVE DIRECTION!
                }

                //Step 4: compute the angle between the vectors
                this->cosRBP = this->R_body_reference.dot(this->V_body_probe) / (this->R_body_reference.norm() * this->V_body_probe.norm());
                this->Angle = math::safe_acos(this->cosRBP);

                //Step 5: apply the constraint - we will constrain the cos(Angle) instead of the Angle itself
                F[Findex++] = this->cosRBP;

                //Step 6: derivatives
                if (needG)
                {
                    //Step 6.1: we are going to be adding components on top of components, so let's start by zeroing out all of our G entries
                    for (size_t stateIndex : {0, 1, 2})
                    {
                        for (size_t Gindex : this->Gindex_constraint_wrt_PositionAfterEvent_variables[stateIndex])
                        {
                            G[Gindex] = 0.0;
                        }
                        for (size_t Gindex : this->Gindex_constraint_wrt_PositionAfterEvent_time_variables[stateIndex])
                        {
                            G[Gindex] = 0.0;
                        }
                    }
                    for (size_t Gindex : this->Gindex_constraint_wrt_time_variables)
                    {
                        G[Gindex] = 0.0;
                    }
					
                    //Step 6.2: derivatives with respect to non-time variables affecting state after boundary event
                    doubleType r_body_reference = this->R_body_reference.norm();
                    double dr_body_reference_drx = (this->R_body_reference(0) / r_body_reference)_GETVALUE;
                    double dr_body_reference_dry = (this->R_body_reference(1) / r_body_reference)_GETVALUE;
                    double dr_body_reference_drz = (this->R_body_reference(2) / r_body_reference)_GETVALUE;
                    doubleType v_body_probe = this->V_body_probe.norm();
                    double dv_body_probe_dvx = (this->V_body_probe(0) / v_body_probe)_GETVALUE;
                    double dv_body_probe_dvy = (this->V_body_probe(1) / v_body_probe)_GETVALUE;
                    double dv_body_probe_dvz = (this->V_body_probe(2) / v_body_probe)_GETVALUE;

                    double rx = this->R_body_reference(0)_GETVALUE;
                    double ry = this->R_body_reference(1)_GETVALUE;
                    double rz = this->R_body_reference(2)_GETVALUE;
                    double vx = this->V_body_probe(0)_GETVALUE;
                    double vy = this->V_body_probe(1)_GETVALUE;
                    double vz = this->V_body_probe(2)_GETVALUE;
                    double r = r_body_reference _GETVALUE;
                    double v = v_body_probe _GETVALUE;

                    this->dcosRBP_dR_body_reference(0) = vx / (r * v) - (rx * vx + ry * vy + rz * vz) * dr_body_reference_drx / (r * r * v);
                    this->dcosRBP_dR_body_reference(1) = vy / (r * v) - (rx * vx + ry * vy + rz * vz) * dr_body_reference_dry / (r * r * v);
                    this->dcosRBP_dR_body_reference(2) = vz / (r * v) - (rx * vx + ry * vy + rz * vz) * dr_body_reference_drz / (r * r * v);

                    this->dcosRBP_dV_body_probe(0) = rx / (r  *v) - (rx * vx + ry * vy + rz * vz) * dv_body_probe_dvx / (r * v * v);
                    this->dcosRBP_dV_body_probe(1) = ry / (r * v) - (rx * vx + ry * vy + rz * vz) * dv_body_probe_dvy / (r * v * v);
                    this->dcosRBP_dV_body_probe(2) = rz / (r * v) - (rx * vx + ry * vy + rz * vz) * dv_body_probe_dvz / (r * v * v);
                                        
                    //Step 6.3: derivatives with respect to non-time variables affecting position (there aren't any, this is ephemeris-pegged)

                    //Step 6.4: derivative with respect to time
                    std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateAfterEvent_wrt_Time = this->myBoundaryEvent->get_Derivatives_of_StateAfterEvent_wrt_Time();
                    
                    for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
                    {
                        for (size_t entryIndex = 0; entryIndex < this->dIndex_with_respect_to_PositionAfterEvent_wrt_Time[stateIndex].size(); ++entryIndex)
                        {
                            size_t dIndex = this->dIndex_with_respect_to_PositionAfterEvent_wrt_Time[stateIndex][entryIndex];
                            size_t Gindex = this->Gindex_constraint_wrt_PositionAfterEvent_time_variables[stateIndex][entryIndex];
                            size_t Xindex = this->jGvar->operator[](Gindex);

                            double Gentry = this->dcosRBP_dR_body_reference(stateIndex) * std::get<2>(Derivatives_of_StateAfterEvent_wrt_Time[dIndex]);

                            if (!this->isCentralBodyConstraint)
                            {
                                Gentry -= this->dcosRBP_dR_body_reference(stateIndex) * this->dR_reference_cb_dt(stateIndex);
                            }

                            G[Gindex] -= this->X_scale_factors->operator[](Xindex)
                                * Gentry;

                        }
                    }//end loop over states
                    

                    //Step 6.5: derivatives with respect to v-infinity
                    for (size_t velocityIndex : {0, 1, 2})
                    {
                        size_t Gindex = this->Gindex_constraint_wrt_Vinfinity[velocityIndex];
                        size_t Xindex = this->jGvar->operator[](Gindex);

                        double Gentry = this->dcosRBP_dV_body_probe(velocityIndex);

                        G[Gindex] += -this->X_scale_factors->operator[](Xindex)
                        	* Gentry;
                    }
					
                }//end derivatives
                
            }//end process_constraint()

            void RBPconstraint::output(std::ofstream& outputfile)
            {
                outputfile << this->myBoundaryEvent->getName() << " " << this->referenceName + "-Body-Probe angle (deg): " << this->Angle _GETVALUE * 180.0 / math::PI << std::endl;
            }//end output()

        }//end namespace SpecializedConstraints
    }//end namespace BoundaryEvents
}//end namespace EMTG