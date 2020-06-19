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

#include "angularMomentumReferenceAngle.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        namespace SpecializedConstraints
        {
            angularMomentumReferenceAngle::angularMomentumReferenceAngle(const std::string& name,
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
                this->hvec = math::Matrix<doubleType>(3, 1, 0.0);
                this->R_probe_ref = math::Matrix<doubleType>(3, 1, 0.0);
                this->R_cb_ref = math::Matrix<doubleType>(3, 1, 0.0);
                this->dR_cb_ref_dt = math::Matrix<double>(3, 1, 0.0);
                this->dR_cb_probe_dt = math::Matrix<double>(3, 1, 0.0);
                this->dV_cb_probe_dt = math::Matrix<double>(3, 1, 0.0);

                //derived class version of myBoundaryEvent has a different type than SpecializedBoundaryConstraintBase uses
                this->myBoundaryEvent = myBoundaryEvent;
            }

            void angularMomentumReferenceAngle::calcbounds()
            {
                // Step 1: parse the constraint definition
                // constraintDefinition is like p0_arrival_angularMomentumReferenceAngle_cb_0.0deg_15.0deg
                std::vector<std::string> ConstraintDefinitionCell;
                std::string splitty_thing = boost::algorithm::to_lower_copy(this->constraintDefinition);
                boost::split(ConstraintDefinitionCell, splitty_thing, boost::is_any_of("_"), boost::token_compress_on);

                //Step 2: what bodies are we doing this with respect to?
                if (ConstraintDefinitionCell[3] == "cb")
                {
                    this->refIsCentralBody = true;
                    this->refName = this->myUniverse->central_body_name;
                }
                else
                {
                    this->refIsCentralBody = false;

                    int bodyIndex = std::stoi(ConstraintDefinitionCell[3]) - 1;
                    this->myReference = &this->myUniverse->bodies[bodyIndex];
                    this->refName = this->myReference->name;
                }

                //Step 3: create the constraint
                // figure out the lower and upper bounds - depends on units
                this->Fdescriptions->push_back(prefix + " angle between angular momentum vector and " + this->refName);

                std::string lower_bound = ConstraintDefinitionCell[4];
                if (lower_bound.find("deg") < 1024)
                {
                    double lower_bound_float = std::stod(boost::erase_all_copy(lower_bound, "deg"));
                    this->Flowerbounds->push_back(lower_bound_float * math::deg2rad);
                }
                else if (lower_bound.find("rad") < 1024)
                {
                    double lower_bound_float = std::stod(boost::erase_all_copy(lower_bound, "rad"));
                    this->Flowerbounds->push_back(lower_bound_float);
                }
                else
                {
                    return throw std::logic_error("Please specify either deg or rad for angularMomentumReferenceAngle constraint bounds: " + this->Fdescriptions->back() + ". Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
                }

                std::string upper_bound = ConstraintDefinitionCell[5];
                if (upper_bound.find("deg") < 1024)
                {
                    double upper_bound_float = std::stod(boost::erase_all_copy(upper_bound, "deg"));
                    this->Fupperbounds->push_back(upper_bound_float * math::deg2rad);
                }
                else if (upper_bound.find("rad") < 1024)
                {
                    double upper_bound_float = std::stod(boost::erase_all_copy(upper_bound, "rad"));
                    this->Fupperbounds->push_back(upper_bound_float);
                }
                else
                {
                    return throw std::logic_error("Please specify either deg or rad for angularMomentumReferenceAngle constraint bounds: " + this->Fdescriptions->back() + ". Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
                }

                //Step 5: sparsity pattern
                //derivatives with respect to anything influencing the boundary event's right-hand position and velocity vectors
                std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateAfterEvent = this->myBoundaryEvent->get_Derivatives_of_StateAfterEvent();
                std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateAfterEvent_wrt_Time = this->myBoundaryEvent->get_Derivatives_of_StateAfterEvent_wrt_Time();
                std::vector<size_t> dIndex_VbeforeEvent_dVinfinity_in = this->myBoundaryEvent->getdIndex_VbeforeEvent_dVinfinity_in();

                for (size_t stateIndex : {0, 1, 2, 3, 4, 5})
                {
                    //non-time variables
                    std::vector<size_t> state_dIndex_with_respect_to_StateAfterEvent;
                    std::vector<size_t> state_Gindex_constraint_wrt_StateAfterEvent_variables;
                    for (size_t dIndex = 0; dIndex < Derivatives_of_StateAfterEvent.size(); ++dIndex)
                    {
                        if (std::get<1>(Derivatives_of_StateAfterEvent[dIndex]) == stateIndex)
                        {
                            state_dIndex_with_respect_to_StateAfterEvent.push_back(dIndex);

                            this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                                std::get<0>(Derivatives_of_StateAfterEvent[dIndex]),
                                state_Gindex_constraint_wrt_StateAfterEvent_variables);
                        }
                    }
                    this->dIndex_with_respect_to_StateAfterEvent.push_back(state_dIndex_with_respect_to_StateAfterEvent);
                    this->Gindex_constraint_wrt_StateAfterEvent_variables.push_back(state_Gindex_constraint_wrt_StateAfterEvent_variables);

                    //time variables
                    std::vector<size_t> state_dIndex_with_respect_to_StateAfterEvent_wrt_Time;
                    std::vector<size_t> state_Gindex_constraint_wrt_StateAfterEvent_time_variables;
                    for (size_t dIndex = 0; dIndex < Derivatives_of_StateAfterEvent_wrt_Time.size(); ++dIndex)
                    {
                        if (std::get<1>(Derivatives_of_StateAfterEvent_wrt_Time[dIndex]) == stateIndex)
                        {
                            state_dIndex_with_respect_to_StateAfterEvent_wrt_Time.push_back(dIndex);

                            this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                                std::get<0>(Derivatives_of_StateAfterEvent_wrt_Time[dIndex]),
                                state_Gindex_constraint_wrt_StateAfterEvent_time_variables);
                        }
                    }
                    this->dIndex_with_respect_to_StateAfterEvent_wrt_Time.push_back(state_dIndex_with_respect_to_StateAfterEvent_wrt_Time);
                    this->Gindex_constraint_wrt_StateAfterEvent_time_variables.push_back(state_Gindex_constraint_wrt_StateAfterEvent_time_variables);

                }

                //derivatives with respect to time variables that affect reference body position
                if (!this->refIsCentralBody)
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

            void angularMomentumReferenceAngle::process_constraint(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG)
            {

                //Step 1: get the state relative to the central body
                math::Matrix<doubleType>& SpacecraftStateRelativeToCentralBody = this->myBoundaryEvent->get_state_after_event();

                //Step 2: get the vector from the central body to the reference body
                if (this->refIsCentralBody) //is the reference body the central body?
                {
                    for (size_t stateIndex : {0, 1, 2})
                        this->R_cb_ref(stateIndex) = 0.0;
                }
                else //reference body is some other body in the Universe
                {
                    //Step 2.1: find the reference body relative to the central body
                    doubleType temp_body_state[12];
                    this->myReference->locate_body(SpacecraftStateRelativeToCentralBody(7),
                        temp_body_state,
                        needG,
                        *this->myOptions);

                    //Step 2.2: compute the vector from the body to the reference body
                    for (size_t stateIndex : {0, 1, 2})
                    {
                        this->R_cb_ref(stateIndex) = temp_body_state[stateIndex];
                        this->dR_cb_ref_dt(stateIndex) = temp_body_state[stateIndex + 6] _GETVALUE;
                    }
                }

                //Step 3: get the vector from the probe to the reference
                for (size_t stateIndex : {0, 1, 2})
                {
                    this->R_probe_ref(stateIndex) = this->R_cb_ref(stateIndex) - SpacecraftStateRelativeToCentralBody(stateIndex);
                }

                //Step 4: compute the angular momentum vector
                this->hvec = SpacecraftStateRelativeToCentralBody.getSubMatrix1D(0, 2).cross(SpacecraftStateRelativeToCentralBody.getSubMatrix1D(3, 5));

//#ifdef AD_INSTRUMENTATION
//                std::ofstream dhvec_drandv_algorithmic("dhvec_drandv_algorithmic.txt", std::ios::trunc);
//
//                for (size_t xxIndex : {13, 14, 15, 16, 17, 18})
//                {
//                    for (size_t hIndex : { 0, 1, 2 })
//                    {
//                        if (hIndex > 0)
//                            dhvec_drandv_algorithmic << ",";
//                        dhvec_drandv_algorithmic << this->hvec(hIndex).getDerivative(xxIndex) / this->X_scale_factors->operator[](xxIndex);
//                    }
//                    dhvec_drandv_algorithmic << std::endl;
//                }
//
//                dhvec_drandv_algorithmic.close();
//#endif

                //Step 5: compute the angle between the vectors
                math::Matrix<doubleType> zhat(3, 1, std::vector<doubleType>({ 0.0, 0.0, 1.0 }));
                doubleType A = (this->hvec.cross(this->R_probe_ref)).dot(zhat);
                doubleType B = this->hvec.dot(this->R_probe_ref);
                this->Angle = atan2(A, B);

                //Step 6: apply the constraint - we will constrain the cos(Angle) instead of the Angle itself
                F[Findex++] = this->Angle;

                //Step 7: derivatives
                if (needG)
                {
                    //Step 6.1: we are going to be adding components on top of components, so let's start by zeroing out all of our G entries
                    for (size_t stateIndex : {0, 1, 2, 3, 4, 5})
                    {
                        for (size_t Gindex : this->Gindex_constraint_wrt_StateAfterEvent_variables[stateIndex])
                        {
                            G[Gindex] = 0.0;
                        }
                        for (size_t Gindex : this->Gindex_constraint_wrt_StateAfterEvent_time_variables[stateIndex])
                        {
                            G[Gindex] = 0.0;
                        }
                    }
                    for (size_t Gindex : this->Gindex_constraint_wrt_time_variables)
                    {
                        G[Gindex] = 0.0;
                    }

                    //Step 6.2: derivatives with respect to non-time variables affecting state after boundary event
                    

                    math::Matrix<double> dA_dhvec = math::Matrix<double>(3, 1, std::vector<double>({ this->R_probe_ref(1)_GETVALUE, -this->R_probe_ref(0)_GETVALUE, 0.0 }));
                    
                    math::Matrix<double> dA_dR_probe_ref = math::Matrix<double>(3, 1, std::vector<double>({ -this->hvec(1)_GETVALUE, this->hvec(0)_GETVALUE, 0.0 }));

                    math::Matrix<double> dB_dhvec(3, 1, 0.0);
                    math::Matrix<double> dB_dR_probe_ref(3, 1, 0.0);

                    for (size_t index : {0, 1, 2})
                    {
                        dB_dhvec(index) = this->R_probe_ref(index)_GETVALUE;
                        dB_dR_probe_ref(index) = this->hvec(index)_GETVALUE;
                    }

                    //derivatives of angle with respect to vectors being angled
                    math::Matrix<double> dAngle_dhvec = dA_dhvec * (B / (A*A + B * B))_GETVALUE + dB_dhvec * (-A / (A*A + B * B))_GETVALUE;
                    math::Matrix<double> dAngle_dR_probe_ref = dA_dR_probe_ref * (B / (A*A + B * B))_GETVALUE + dB_dR_probe_ref * (-A / (A*A + B * B))_GETVALUE;

                    //derivatives of h_vec with respect to r and v
                    math::Matrix<double> R(3, 1, 0.0);
                    math::Matrix<double> V(3, 1, 0.0);
                    R(0) = SpacecraftStateRelativeToCentralBody(0)_GETVALUE;
                    R(1) = SpacecraftStateRelativeToCentralBody(1)_GETVALUE;
                    R(2) = SpacecraftStateRelativeToCentralBody(2)_GETVALUE;
                    V(0) = SpacecraftStateRelativeToCentralBody(3)_GETVALUE;
                    V(1) = SpacecraftStateRelativeToCentralBody(4)_GETVALUE;
                    V(2) = SpacecraftStateRelativeToCentralBody(5)_GETVALUE;
                    math::Matrix<double> dhvec_dRandV = R.crossDerivative(V); //makes a 6x3: upper 3x3 is dhvec_dR, lower 3x3 is dhvec_dV

                    //So now we have dAngle_dhvec and dhvec_dRandV. We need to chain them to create dAngle_dRandV
                    //each entry in dAngle_dRandV should be the derivative of Angle with respect to that state entry
                    math::Matrix<double> dAngle_dRandV = dhvec_dRandV * dAngle_dhvec;

                    //note that the position entries also have to carry a contribution from R_probe_ref, which is just dR_probe_ref_dR_cb_probe (which is -1) times dAngle_dR_probe_ref
                    for (size_t stateIndex : {0, 1, 2})
                        dAngle_dRandV(stateIndex) -= dAngle_dR_probe_ref(stateIndex);



//#ifdef AD_INSTRUMENTATION
//                    math::Matrix<double> dA_dRandV =  dhvec_dRandV * dA_dhvec;
//
//                    dA_dRandV.print_to_file("dA_dRandV_analytical.txt");
//                    std::ofstream dA_dRandV_algorithmic("dA_dRandV_algorithmic.txt", std::ios::trunc);
//
//                    for (size_t xxIndex : {13, 14, 15, 16, 17, 18})
//                    {
//                        dA_dRandV_algorithmic << A.getDerivative(xxIndex) / this->X_scale_factors->operator[](xxIndex) << std::endl;
//                    }
//
//                    dA_dRandV_algorithmic.close();
//#endif

                    //the time derivatives are dAngle_dt = dAngle_dhvec * dhvec_dt + dAngle_dR_probe_ref * dR_probe_ref_dt
                    //                                   = dAngle_dhvec * dAngle_dRandV * dRandV_dt + dAngle_dR_probe_ref * (dR_cb_ref_dt - dR_cb_probe_dt)

                    //Step 6.3: derivatives with respect to non-time variables affecting ONLY the probe's position
                    std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateAfterEvent = this->myBoundaryEvent->get_Derivatives_of_StateAfterEvent();
                    
                    for (size_t stateIndex = 0; stateIndex < 6; ++stateIndex)
                    {
                        for (size_t entryIndex = 0; entryIndex < this->dIndex_with_respect_to_StateAfterEvent[stateIndex].size(); ++entryIndex)
                        {
                            size_t dIndex = this->dIndex_with_respect_to_StateAfterEvent[stateIndex][entryIndex];
                            size_t Gindex = this->Gindex_constraint_wrt_StateAfterEvent_variables[stateIndex][entryIndex];
                            size_t Xindex = this->jGvar->operator[](Gindex);

                            double Gentry = dAngle_dRandV(stateIndex) * std::get<2>(Derivatives_of_StateAfterEvent[dIndex]);

                            G[Gindex] += this->X_scale_factors->operator[](Xindex)
                                * Gentry;

                        }
                    }//end loop over states

                    //Step 6.4: derivative with respect to time
                    std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StateAfterEvent_wrt_Time = this->myBoundaryEvent->get_Derivatives_of_StateAfterEvent_wrt_Time();
                    math::Matrix<double> dRandV_dt(6, 1, 0.0);
                    
                    //if we have time dependency in the boundary state, this loop will also scoop up time dependencies of the reference state
                    if (this->dIndex_with_respect_to_StateAfterEvent_wrt_Time[0].size() > 0)
                    {
                        for (size_t stateIndex = 0; stateIndex < 6; ++stateIndex)
                        {

                            for (size_t entryIndex = 0; entryIndex < this->dIndex_with_respect_to_StateAfterEvent_wrt_Time[stateIndex].size(); ++entryIndex)
                            {
                                size_t dIndex = this->dIndex_with_respect_to_StateAfterEvent_wrt_Time[stateIndex][entryIndex];
                                size_t Gindex = this->Gindex_constraint_wrt_StateAfterEvent_time_variables[stateIndex][entryIndex];
                                size_t Xindex = this->jGvar->operator[](Gindex);

                                dRandV_dt.assign_zeros();
                                dRandV_dt(stateIndex) = std::get<2>(Derivatives_of_StateAfterEvent_wrt_Time[dIndex]);

                                double dAngle_dt = dAngle_dRandV.dot(dRandV_dt) + dAngle_dR_probe_ref.dot(dR_cb_ref_dt - dR_cb_probe_dt);

                                G[Gindex] += this->X_scale_factors->operator[](Xindex)
                                    * dAngle_dt;
                            }
                        }//end loop over states
                    }
                    else //if only the reference is time dependent
                    {
                        double dAngle_dt = dAngle_dR_probe_ref.dot(dR_cb_ref_dt - dR_cb_probe_dt);

                        for (size_t Gindex : this->Gindex_constraint_wrt_time_variables)
                        {
                            size_t Xindex = this->jGvar->operator[](Gindex);


                            G[Gindex] += this->X_scale_factors->operator[](Xindex)
                                * dAngle_dt;
                        }
                    }
                }//end derivatives

            }//end process_constraint()

            void angularMomentumReferenceAngle::output(std::ofstream& outputfile)
            {
                outputfile << this->myBoundaryEvent->getName() << " angle between angular momentum vector and vector to " << this->refName + " (deg): " << this->Angle _GETVALUE * 180.0 / math::PI << std::endl;
            }//end output()

        }//end namespace SpecializedConstraints
    }//end namespace BoundaryEvents
}//end namespace EMTG