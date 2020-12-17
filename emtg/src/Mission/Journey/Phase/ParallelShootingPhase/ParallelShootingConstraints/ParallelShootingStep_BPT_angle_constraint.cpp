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

//ParallelShootingStep body-probe-thrust constraint
//Jacob Englander 4-13-2018

#include "ParallelShootingStep_BPT_angle_constraint.h"

#include "boost/algorithm/string/split.hpp"
#include "boost/algorithm/string/erase.hpp"

#include <string>
#include <tuple>

namespace EMTG
{
    namespace Phases
    {
        ParallelShootingStep_BPT_angle_constraint::ParallelShootingStep_BPT_angle_constraint(const std::string& name,
            const size_t& journeyIndex,
            const size_t& phaseIndex,
            const size_t& stepIndex,
            const size_t& subStepIndex,
            const size_t& stageIndex,
            ParallelShootingStep* myStep,
            Astrodynamics::universe* myUniverse,
            HardwareModels::Spacecraft* mySpacecraft,
            missionoptions* myOptions,
            const std::string& ConstraintDefinition) :
            ParallelShootingStep_maneuver_constraint(name,
                journeyIndex,
                phaseIndex,
                stepIndex,
                subStepIndex,
                stageIndex,
                myStep,
                myUniverse,
                mySpacecraft,
                myOptions,
                ConstraintDefinition)
        {

            std::vector<std::string> ConstraintDefinitionCell;
            boost::split(ConstraintDefinitionCell, this->ConstraintDefinition, boost::is_any_of("_"), boost::token_compress_on);

            if (boost::to_lower_copy(ConstraintDefinitionCell[2]) == "cb")
            {
                this->isDefinedRelativeToCentralBody = true;
                this->myBody = &this->myUniverse->central_body;
            }
            else
            {
                this->isDefinedRelativeToCentralBody = false;
                int bodyIndex = std::stoi(ConstraintDefinitionCell[2]);
                this->myBody = &this->myUniverse->bodies[bodyIndex];
            }

            this->name = this->name + " body-probe-thrust angle relative to " + this->myBody->name;

            //components of BPT minimum angle
            std::vector<std::string> fit_cell;
            boost::split(fit_cell, ConstraintDefinitionCell[3], boost::is_any_of(","), boost::token_compress_on);
            if (fit_cell.size() > 1)
            {
                this->MinimumAngle_a = std::stod(fit_cell[0]);
                this->MinimumAngle_b = std::stod(fit_cell[1]);
                this->MinimumAngle_c = std::stod(fit_cell[2]);
                this->MinimumAngle_d = std::stod(fit_cell[3]);
                this->ConstantMinimumAngle = false;
            }
            else
            {
                this->MinimumAngle_a = cos(std::stod(fit_cell[0]) * math::PI / 180.0);
                this->ConstantMinimumAngle = true;
            }

            //components of BPT maximum angle
            fit_cell.clear();
            boost::split(fit_cell, ConstraintDefinitionCell[4], boost::is_any_of(","), boost::token_compress_on);
            if (fit_cell.size() > 1)
            {
                this->MaximumAngle_a = std::stod(fit_cell[0]);
                this->MaximumAngle_b = std::stod(fit_cell[1]);
                this->MaximumAngle_c = std::stod(fit_cell[2]);
                this->MaximumAngle_d = std::stod(fit_cell[3]);
                this->ConstantMaximumAngle = false;
            }
            else
            {
                this->MaximumAngle_a = cos(std::stod(fit_cell[0]) * math::PI / 180.0);
                this->ConstantMaximumAngle = true;
            }

            this->spacecraft_to_reference = math::Matrix<doubleType>(3, 1, 0.0);
            this->dspacecraft_to_reference_depoch = math::Matrix<double>(3, 1, 0.0);
        }//end constructor

        void ParallelShootingStep_BPT_angle_constraint::calcbounds()
        {
            //there are two constraints
            //one for the minimum angle: -inf <= cos(theta) - cos(min_angle) <= 0.0
            //one for the maximum angle: 0.0 <= cos(theta) - cos(max_angle) <= inf

            //minimum angle
            {
                //Step 1: create the constraint
                this->Flowerbounds->push_back(-math::LARGE);
                this->Fupperbounds->push_back(0.0);
                this->Fdescriptions->push_back(this->name + " minimum");

                //Step 2: sparsity pattern

                //Step 2.1: non-time variables that affect position
                std::vector<size_t>& ListOfVariablesAffectingCurrentStepLeftPosition = this->myStep->getListOfVariablesAffectingCurrentStepLeftPosition();

                std::vector< std::vector< std::tuple<size_t, size_t> > >& DerivativesOfCurrentStepLeftStateByVariable = this->myStep->getDerivativesOfCurrentStepLeftStateByVariable();
                std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StepLeftState = this->myStep->get_Derivatives_of_StateStepLeftInertial();//Xindex, stateIndex, derivative value

                for (size_t varIndex = 0; varIndex < ListOfVariablesAffectingCurrentStepLeftPosition.size(); ++varIndex)
                {
                    size_t Xindex = ListOfVariablesAffectingCurrentStepLeftPosition[varIndex];
                    std::vector<size_t> dIndex_constraint_wrt_StepLeftPositionVariable;

                    std::vector< std::tuple<size_t, size_t> > DerivativesOfCurrentStepLeftStateThisVariable = DerivativesOfCurrentStepLeftStateByVariable[varIndex];

                    bool madeEntryFlag = false;

                    for (std::tuple<size_t, size_t> thisDerivativeEntry : DerivativesOfCurrentStepLeftStateThisVariable)
                    {
                        size_t stateIndex = std::get<0>(thisDerivativeEntry);

                        if (stateIndex < 3)
                        {
                            size_t dIndex = std::get<1>(thisDerivativeEntry);
                            dIndex_constraint_wrt_StepLeftPositionVariable.push_back(dIndex);

                            if (!madeEntryFlag)
                            {
                                this->create_sparsity_entry(Fdescriptions->size() - 1,
                                    Xindex,
                                    this->Gindices_BPT_constraintMinimumAngle_wrt_StepLeftPosition);

                                madeEntryFlag = true;
                            }
                        }
                    }

                    this->dIndex_MinimumAngle_constraint_wrt_StepLeftPosition.push_back(dIndex_constraint_wrt_StepLeftPositionVariable);
                }

                //Step 2.2: time variables
                if (!this->isDefinedRelativeToCentralBody)
                {
                    std::vector<size_t>& ListOfTimeVariablesAffectingCurrentStepLeftState = this->myStep->getListOfTimeVariablesAffectingCurrentStepLeftState();
                    for (size_t varIndex = 0; varIndex < ListOfTimeVariablesAffectingCurrentStepLeftState.size(); ++varIndex)
                    {
                        size_t Xindex = ListOfTimeVariablesAffectingCurrentStepLeftState[varIndex];

                        this->create_sparsity_entry(Fdescriptions->size() - 1,
                            Xindex,
                            this->Gindices_BPT_constraintMinimumAngle_wrt_StepLeftTime);
                    }
                }//end time variable dependencies, which only exist if the constraint is NOT defined relative to the central body

                //Step 2.3: control variables
                std::vector<size_t>& Xindices_control = this->myStep->getXindices_control(this->subStepIndex);
                for (size_t varIndex = 0; varIndex < 3; ++varIndex)
                {
                    size_t Xindex = Xindices_control[varIndex];

                    this->create_sparsity_entry(Fdescriptions->size() - 1,
                        Xindex,
                        this->Gindices_BPT_constraintMinimumAngle_wrt_control);
                }
            }

            //maximum angle
            {
                //Step 1: create the constraint
                this->Flowerbounds->push_back(0.0);
                this->Fupperbounds->push_back(math::LARGE);
                this->Fdescriptions->push_back(this->name + " maximum");

                //Step 2: sparsity pattern

                //Step 2.1: non-time variables that affect state
                std::vector<size_t>& ListOfVariablesAffectingCurrentStepLeftPosition = this->myStep->getListOfVariablesAffectingCurrentStepLeftPosition();

                std::vector< std::vector< std::tuple<size_t, size_t> > >& DerivativesOfCurrentStepLeftStateByVariable = this->myStep->getDerivativesOfCurrentStepLeftStateByVariable();
                std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StepLeftState = this->myStep->get_Derivatives_of_StateStepLeftInertial();//Xindex, stateIndex, derivative value

                for (size_t varIndex = 0; varIndex < ListOfVariablesAffectingCurrentStepLeftPosition.size(); ++varIndex)
                {
                    size_t Xindex = ListOfVariablesAffectingCurrentStepLeftPosition[varIndex];
                    std::vector<size_t> dIndex_constraint_wrt_StepLeftPositionVariable;

                    std::vector< std::tuple<size_t, size_t> > DerivativesOfCurrentStepLeftStateThisVariable = DerivativesOfCurrentStepLeftStateByVariable[varIndex];

                    bool madeEntryFlag = false;

                    for (std::tuple<size_t, size_t> thisDerivativeEntry : DerivativesOfCurrentStepLeftStateThisVariable)
                    {
                        size_t stateIndex = std::get<0>(thisDerivativeEntry);

                        if (stateIndex < 3)
                        {
                            size_t dIndex = std::get<1>(thisDerivativeEntry);
                            dIndex_constraint_wrt_StepLeftPositionVariable.push_back(dIndex);

                            if (!madeEntryFlag)
                            {
                                this->create_sparsity_entry(Fdescriptions->size() - 1,
                                    Xindex,
                                    this->Gindices_BPT_constraintMaximumAngle_wrt_StepLeftPosition);

                                madeEntryFlag = true;
                            }
                        }
                    }

                    this->dIndex_MaximumAngle_constraint_wrt_StepLeftPosition.push_back(dIndex_constraint_wrt_StepLeftPositionVariable);
                }

                //Step 2.2: time variables
                if (!this->isDefinedRelativeToCentralBody)
                {
                    std::vector<size_t>& ListOfTimeVariablesAffectingCurrentStepLeftState = this->myStep->getListOfTimeVariablesAffectingCurrentStepLeftState();
                    for (size_t varIndex = 0; varIndex < ListOfTimeVariablesAffectingCurrentStepLeftState.size(); ++varIndex)
                    {
                        size_t Xindex = ListOfTimeVariablesAffectingCurrentStepLeftState[varIndex];

                        this->create_sparsity_entry(Fdescriptions->size() - 1,
                            Xindex,
                            this->Gindices_BPT_constraintMaximumAngle_wrt_StepLeftTime);
                    }
                }//end time variable dependencies, which only exist if the constraint is NOT relative to the central body

                //Step 2.3: control variables
                std::vector<size_t>& Xindices_control = this->myStep->getXindices_control(this->subStepIndex);
                for (size_t varIndex = 0; varIndex < 3; ++varIndex)
                {
                    size_t Xindex = Xindices_control[varIndex];

                    this->create_sparsity_entry(Fdescriptions->size() - 1,
                        Xindex,
                        this->Gindices_BPT_constraintMaximumAngle_wrt_control);
                }
            }
        }//end calcbounds


        void ParallelShootingStep_BPT_angle_constraint::process_constraint(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //Step 1: evaluate the constraint

            //Step 1.1: get the control vector
            math::Matrix<doubleType>& ControlVector = this->myStep->getControlVector(this->subStepIndex);

            //Step 1.2: get the position vector of the spacecraft relative to the body
            math::Matrix<doubleType>& SpacecraftState = this->myStep->get_StateStepLeftInertial();
            doubleType DistanceToReferenceBody;
                        
            doubleType body_state[12];
            this->myBody->locate_body(
                SpacecraftState(7),
                body_state,
                (needG),
                *this->myOptions);
            for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
            {
                this->spacecraft_to_reference(stateIndex) = -(SpacecraftState(stateIndex) - body_state[stateIndex]);
                this->dspacecraft_to_reference_depoch(stateIndex) = -body_state[6 + stateIndex] _GETVALUE;
            }

            DistanceToReferenceBody = this->spacecraft_to_reference.norm();
            doubleType r_AU = DistanceToReferenceBody / this->myUniverse->LU;

            //Step 1.3: evaluate the angle between the vector to the reference and the control vector
            this->cosAngleThrustVectorToSun = this->spacecraft_to_reference.dot(ControlVector) / this->spacecraft_to_reference.norm() / ControlVector.norm();


            //Step 1.4: minimum angle constraint
            // -inf <= cos(theta) - cos(min_angle) <= 0.0
            doubleType MinBoundary, MaxBoundary;

            if (this->ConstantMinimumAngle)
            {
                MinBoundary = this->MinimumAngle_a;
            }
            else
            {
                MinBoundary = this->MinimumAngle_d + (this->MinimumAngle_a - this->MinimumAngle_d)
                    / (1.0 + pow(r_AU / this->MinimumAngle_c, this->MinimumAngle_b));
            }

            F[Findex++] = this->cosAngleThrustVectorToSun - MinBoundary;

            //Step 1.5: maximum angle constraint
            // one for the maximum angle: 0.0 <= cos(theta) - cos(max_angle) <= inf
            if (this->ConstantMaximumAngle)
            {
                MaxBoundary = this->MaximumAngle_a;
            }
            else
            {
                MaxBoundary = this->MaximumAngle_d + (this->MaximumAngle_a - this->MaximumAngle_d)
                    / (1.0 + pow(r_AU / this->MaximumAngle_c, this->MaximumAngle_b));
            }
            F[Findex++] = this->cosAngleThrustVectorToSun - MaxBoundary;
            
            //Step 2: derivatives
            if (needG)
            {
                math::Matrix<double> dConstraintPointState_dDecisionVariable(10, 1, 0.0);

                double rx = this->spacecraft_to_reference(0) _GETVALUE;
                double ry = this->spacecraft_to_reference(1) _GETVALUE;
                double rz = this->spacecraft_to_reference(2) _GETVALUE;
                double ux = ControlVector(0) _GETVALUE;
                double uy = ControlVector(1) _GETVALUE;
                double uz = ControlVector(2) _GETVALUE;

                double r = this->spacecraft_to_reference.norm() _GETVALUE;
                double u = ControlVector.norm() _GETVALUE;

                double dr_drx = rx / r;
                double dr_dry = ry / r;
                double dr_drz = rz / r;
                double du_dux = ux / u;
                double du_duy = uy / u;
                double du_duz = uz / u;


                //first let's do the derivative of the actual thrust angle component, which is the same for both constraints
                double dFdU[3];

                
                double dFdR[6];

                dFdR[0] = ux / (r*u) - (rx*ux + ry * uy + rz * uz)*dr_drx / (r * r * u);
                dFdR[1] = uy / (r*u) - (rx*ux + ry * uy + rz * uz)*dr_dry / (r * r * u);
                dFdR[2] = uz / (r*u) - (rx*ux + ry * uy + rz * uz)*dr_drz / (r * r * u);
                dFdR[3] = 0.0;
                dFdR[4] = 0.0;
                dFdR[5] = 0.0;

                dFdU[0] = rx / (r*u) - (rx*ux + ry * uy + rz * uz)*du_dux / (r * u * u);
                dFdU[1] = ry / (r*u) - (rx*ux + ry * uy + rz * uz)*du_duy / (r * u * u);
                dFdU[2] = rz / (r*u) - (rx*ux + ry * uy + rz * uz)*du_duz / (r * u * u);

                //Step 2.1: non-time variables affecting spacecraft state
                std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StepLeftState = this->myStep->get_Derivatives_of_StateStepLeftInertial();//Xindex, stateIndex, derivative value

                for (size_t varIndex = 0; varIndex < this->dIndex_MinimumAngle_constraint_wrt_StepLeftPosition.size(); ++varIndex)
                {
                    dConstraintPointState_dDecisionVariable.assign_zeros();

                    for (size_t dIndex : this->dIndex_MinimumAngle_constraint_wrt_StepLeftPosition[varIndex])
                    {
                        size_t stateIndex = std::get<1>(Derivatives_of_StepLeftState[dIndex]);

                        dConstraintPointState_dDecisionVariable(stateIndex) = std::get<2>(Derivatives_of_StepLeftState[dIndex]);
                    }

                    //minimum angle constraint
                    size_t Gindex = this->Gindices_BPT_constraintMinimumAngle_wrt_StepLeftPosition[varIndex];
                    size_t Xindex = this->jGvar->operator[](Gindex);

                    G[Gindex] = this->X_scale_factors->operator[](Xindex)
                        *   -(dFdR[0] * dConstraintPointState_dDecisionVariable(0)
                            + dFdR[1] * dConstraintPointState_dDecisionVariable(1)
                            + dFdR[2] * dConstraintPointState_dDecisionVariable(2));

                    //maximum angle constraint
                    Gindex = this->Gindices_BPT_constraintMaximumAngle_wrt_StepLeftPosition[varIndex];
                    Xindex = this->jGvar->operator[](Gindex);

                    G[Gindex] = this->X_scale_factors->operator[](Xindex)
                        *   -(dFdR[0] * dConstraintPointState_dDecisionVariable(0)
                            + dFdR[1] * dConstraintPointState_dDecisionVariable(1)
                            + dFdR[2] * dConstraintPointState_dDecisionVariable(2));
                }//end non-time derivatives

                //Step 2.2: time variables affecting spacecraft state
                if (!this->isDefinedRelativeToCentralBody)
                {
                    std::vector< std::vector< std::tuple<size_t, size_t> > >& DerivativesOfCurrentStepLeftStateByTimeVariable = this->myStep->getDerivativesOfCurrentStepLeftStateByTimeVariable();
                    std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StepLeftState_wrt_Time = this->myStep->get_Derivatives_of_StateStepLeftInertial_wrt_Time();//Xindex, stateIndex, derivative value

                    for (size_t varIndex = 0; varIndex < DerivativesOfCurrentStepLeftStateByTimeVariable.size(); ++varIndex)
                    {
                        dConstraintPointState_dDecisionVariable.assign_zeros();

                        for (size_t entryIndex = 0; entryIndex < DerivativesOfCurrentStepLeftStateByTimeVariable[varIndex].size(); ++entryIndex)
                        {
                            size_t stateIndex = std::get<0>(DerivativesOfCurrentStepLeftStateByTimeVariable[varIndex][entryIndex]);
                            size_t dIndex = std::get<1>(DerivativesOfCurrentStepLeftStateByTimeVariable[varIndex][entryIndex]);

                            dConstraintPointState_dDecisionVariable(stateIndex) = std::get<2>(Derivatives_of_StepLeftState_wrt_Time[dIndex]);
                        }

                        size_t Gindex = this->Gindices_BPT_constraintMinimumAngle_wrt_StepLeftTime[varIndex];
                        size_t Xindex = this->jGvar->operator[](Gindex);

                        G[Gindex] = this->X_scale_factors->operator[](Xindex)
                            * -(dFdR[0] * dConstraintPointState_dDecisionVariable(0) - this->dspacecraft_to_reference_depoch(0)
                                + dFdR[1] * dConstraintPointState_dDecisionVariable(1) - this->dspacecraft_to_reference_depoch(1)
                                + dFdR[2] * dConstraintPointState_dDecisionVariable(2) - this->dspacecraft_to_reference_depoch(2));

                        //maximum angle constraint
                        Gindex = this->Gindices_BPT_constraintMaximumAngle_wrt_StepLeftTime[varIndex];
                        Xindex = this->jGvar->operator[](Gindex);

                        G[Gindex] = this->X_scale_factors->operator[](Xindex)
                            * -(dFdR[0] * dConstraintPointState_dDecisionVariable(0) - this->dspacecraft_to_reference_depoch(0)
                                + dFdR[1] * dConstraintPointState_dDecisionVariable(1) - this->dspacecraft_to_reference_depoch(1)
                                + dFdR[2] * dConstraintPointState_dDecisionVariable(2) - this->dspacecraft_to_reference_depoch(2));
                    }//end time derivatives
                }//end switch for time variables, which only exist if the constraint is NOT defined relative to the central body

                //Step 2.3: control
                for (size_t controlIndex = 0; controlIndex < 3; ++controlIndex)
                {
                    //minimum angle constraint
                    size_t Gindex = this->Gindices_BPT_constraintMinimumAngle_wrt_control[controlIndex];
                    size_t Xindex = this->jGvar->operator[](Gindex);

                    G[Gindex] = this->X_scale_factors->operator[](Xindex)
                        * dFdU[controlIndex];

                    //maximum angle constraint
                    Gindex = this->Gindices_BPT_constraintMaximumAngle_wrt_control[controlIndex];
                    Xindex = this->jGvar->operator[](Gindex);

                    G[Gindex] = this->X_scale_factors->operator[](Xindex)
                        * dFdU[controlIndex];
                }

                //Step 3: bound for minimum angle constraint
                if (!this->ConstantMinimumAngle)
                {
                    double a = this->MinimumAngle_a;
                    double b = this->MinimumAngle_b;
                    double c = this->MinimumAngle_c;
                    double d = this->MinimumAngle_d;
                    doubleType dMinBoundary_dr = b * pow(r_AU / c, b) * (a - d) / r_AU
                        / (pow(r_AU / c, b) + 1.0) / (pow(r_AU / c, b) + 1.0) / this->myUniverse->LU;

                    //Step 3.1: position
                    std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StepLeftState = this->myStep->get_Derivatives_of_StateStepLeftInertial();//Xindex, stateIndex, derivative value

                    for (size_t varIndex = 0; varIndex < this->dIndex_MinimumAngle_constraint_wrt_StepLeftPosition.size(); ++varIndex)
                    {
                        dConstraintPointState_dDecisionVariable.assign_zeros();

                        for (size_t dIndex : this->dIndex_MinimumAngle_constraint_wrt_StepLeftPosition[varIndex])
                        {
                            size_t stateIndex = std::get<1>(Derivatives_of_StepLeftState[dIndex]);

                            dConstraintPointState_dDecisionVariable(stateIndex) = std::get<2>(Derivatives_of_StepLeftState[dIndex]);
                        }

                        size_t Gindex = this->Gindices_BPT_constraintMinimumAngle_wrt_StepLeftPosition[varIndex];
                        size_t Xindex = this->jGvar->operator[](Gindex);

                        G[Gindex] += this->X_scale_factors->operator[](Xindex) * -dMinBoundary_dr _GETVALUE
                            * (rx / r * dConstraintPointState_dDecisionVariable(0)
                                + ry / r * dConstraintPointState_dDecisionVariable(1)
                                + rz / r * dConstraintPointState_dDecisionVariable(2));
                    }//end non-time derivatives

                    //Step 3.2: time
                    if (!this->isDefinedRelativeToCentralBody)
                    {
                        std::vector< std::vector< std::tuple<size_t, size_t> > >& DerivativesOfCurrentStepLeftStateByTimeVariable = this->myStep->getDerivativesOfCurrentStepLeftStateByTimeVariable();
                        std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StepLeftState_wrt_Time = this->myStep->get_Derivatives_of_StateStepLeftInertial_wrt_Time();//Xindex, stateIndex, derivative value

                        for (size_t varIndex = 0; varIndex < DerivativesOfCurrentStepLeftStateByTimeVariable.size(); ++varIndex)
                        {
                            dConstraintPointState_dDecisionVariable.assign_zeros();

                            for (size_t entryIndex = 0; entryIndex < DerivativesOfCurrentStepLeftStateByTimeVariable[varIndex].size(); ++entryIndex)
                            {
                                size_t stateIndex = std::get<0>(DerivativesOfCurrentStepLeftStateByTimeVariable[varIndex][entryIndex]);
                                size_t dIndex = std::get<1>(DerivativesOfCurrentStepLeftStateByTimeVariable[varIndex][entryIndex]);

                                dConstraintPointState_dDecisionVariable(stateIndex) = std::get<2>(Derivatives_of_StepLeftState_wrt_Time[dIndex]);
                            }

                            size_t Gindex = this->Gindices_BPT_constraintMinimumAngle_wrt_StepLeftTime[varIndex];
                            size_t Xindex = this->jGvar->operator[](Gindex);

                            G[Gindex] += this->X_scale_factors->operator[](Xindex) * -dMinBoundary_dr _GETVALUE
                                * (rx / r * this->dspacecraft_to_reference_depoch(0)
                                    + ry / r * this->dspacecraft_to_reference_depoch(1)
                                    + rz / r * this->dspacecraft_to_reference_depoch(2));
                        }//end time derivatives
                    }//end switch for time variables, which only exist if the constraint is NOT defined relative to the central body
                }//end bound for minimum angle constraint
                
                //Step 4: bound for maximum angle constraint
                if (!this->ConstantMaximumAngle)
                {
                    double a = this->MaximumAngle_a;
                    double b = this->MaximumAngle_b;
                    double c = this->MaximumAngle_c;
                    double d = this->MaximumAngle_d;
                    doubleType dMaxBoundary_dr = b * pow(r_AU / c, b) * (a - d) / r_AU
                        / (pow(r_AU / c, b) + 1.0) / (pow(r_AU / c, b) + 1.0) / this->myUniverse->LU;

                    //Step 3.1: position
                    std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StepLeftState = this->myStep->get_Derivatives_of_StateStepLeftInertial();//Xindex, stateIndex, derivative value

                    for (size_t varIndex = 0; varIndex < this->dIndex_MaximumAngle_constraint_wrt_StepLeftPosition.size(); ++varIndex)
                    {
                        dConstraintPointState_dDecisionVariable.assign_zeros();

                        for (size_t dIndex : this->dIndex_MaximumAngle_constraint_wrt_StepLeftPosition[varIndex])
                        {
                            size_t stateIndex = std::get<1>(Derivatives_of_StepLeftState[dIndex]);

                            dConstraintPointState_dDecisionVariable(stateIndex) = std::get<2>(Derivatives_of_StepLeftState[dIndex]);
                        }

                        size_t Gindex = this->Gindices_BPT_constraintMaximumAngle_wrt_StepLeftPosition[varIndex];
                        size_t Xindex = this->jGvar->operator[](Gindex);

                        G[Gindex] += this->X_scale_factors->operator[](Xindex) * -dMaxBoundary_dr _GETVALUE
                            * (rx / r * dConstraintPointState_dDecisionVariable(0)
                                + ry / r * dConstraintPointState_dDecisionVariable(1)
                                + rz / r * dConstraintPointState_dDecisionVariable(2));
                    }//end non-time derivatives

                    //Step 3.2: time
                    if (!this->isDefinedRelativeToCentralBody)
                    {
                        std::vector< std::vector< std::tuple<size_t, size_t> > >& DerivativesOfCurrentStepLeftStateByTimeVariable = this->myStep->getDerivativesOfCurrentStepLeftStateByTimeVariable();
                        std::vector< std::tuple<size_t, size_t, double> >& Derivatives_of_StepLeftState_wrt_Time = this->myStep->get_Derivatives_of_StateStepLeftInertial_wrt_Time();//Xindex, stateIndex, derivative value

                        for (size_t varIndex = 0; varIndex < DerivativesOfCurrentStepLeftStateByTimeVariable.size(); ++varIndex)
                        {
                            dConstraintPointState_dDecisionVariable.assign_zeros();

                            for (size_t entryIndex = 0; entryIndex < DerivativesOfCurrentStepLeftStateByTimeVariable[varIndex].size(); ++entryIndex)
                            {
                                size_t stateIndex = std::get<0>(DerivativesOfCurrentStepLeftStateByTimeVariable[varIndex][entryIndex]);
                                size_t dIndex = std::get<1>(DerivativesOfCurrentStepLeftStateByTimeVariable[varIndex][entryIndex]);

                                dConstraintPointState_dDecisionVariable(stateIndex) = std::get<2>(Derivatives_of_StepLeftState_wrt_Time[dIndex]);
                            }

                            size_t Gindex = this->Gindices_BPT_constraintMaximumAngle_wrt_StepLeftTime[varIndex];
                            size_t Xindex = this->jGvar->operator[](Gindex);

                            G[Gindex] += this->X_scale_factors->operator[](Xindex) * -dMaxBoundary_dr _GETVALUE
                                * (rx / r * this->dspacecraft_to_reference_depoch(0)
                                    + ry / r * this->dspacecraft_to_reference_depoch(1)
                                    + rz / r * this->dspacecraft_to_reference_depoch(2));
                        }//end time derivatives
                    }//end switch for time variables, which only exist if the constraint is NOT defined relative to the central body
                }//end bound for Maximum angle constraint
            }//end derivatives
        }//end process_constraint()
    }//close namespace Phases
}//close namespace EMTG