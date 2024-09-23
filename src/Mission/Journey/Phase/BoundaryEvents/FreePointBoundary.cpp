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

#include "FreePointBoundary.h"
#include "PropagatorFactory.h"
#include "IntegrationSchemeFactory.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        //constructors
        FreePointBoundary::FreePointBoundary() :
            StateBeforeEventEncoded(math::Matrix<doubleType>(8, 1, 0.0)),
            dStateBeforeEventEncoded_dEncodedElements(math::Matrix<doubleType>(8, math::identity)),
            dStateBeforeEvent_dStateBeforeEventEncoded(math::Matrix<doubleType>(8, math::identity)),
            dStateBeforeEventICRF_dStateBeforeEventEncoded(math::Matrix<doubleType>(8, math::identity)),
            dStateBeforeEvent_dStateBeforeEventICRF(math::Matrix<doubleType>(8, math::identity)),
            R_from_local_to_ICRF(math::Matrix<doubleType>(8, math::identity)),
            dR_from_local_to_ICRF_dt(math::Matrix<doubleType>(8, math::identity)),
            dIndex_StateBeforeEvent_wrt_EncodedState(7)
        {
            this->LeftBoundaryIsABody = false;
        }

        FreePointBoundary::FreePointBoundary(const std::string& name,
            const size_t& journeyIndex,
            const size_t& phaseIndex,
            size_t& stageIndex,
            Astrodynamics::universe* Universe,
            HardwareModels::Spacecraft* mySpacecraft,
            missionoptions* myOptions) :
            FreePointBoundary::FreePointBoundary()
        {
            this->initialize(name,
                             journeyIndex,
                             phaseIndex,
                             stageIndex,
                             Universe,
                             mySpacecraft,
                             myOptions);
        }
        
        void FreePointBoundary::initialize(const std::string& name,
                const size_t& journeyIndex,
                const size_t& phaseIndex,
                size_t& stageIndex,
                Astrodynamics::universe* Universe,
                HardwareModels::Spacecraft* mySpacecraft,
                missionoptions* myOptions)
        {
            //base class
            this->BoundaryEventBase::initialize(name,
                journeyIndex, 
                phaseIndex,
                stageIndex,
                Universe, 
                mySpacecraft,
                myOptions);

            if (this->myPropagatorType == PropagatorType::KeplerPropagator)
            {
                this->STM = math::Matrix<double>(6, math::identity);
                this->dPropagatedStatedIndependentVariable = math::Matrix<double>(6, 1, 0.0);
                this->StateBeforeEventBeforePropagation = math::Matrix<doubleType>(8, 1, 0.0);
            }
            else
            {
                this->STM = math::Matrix<double>(11, math::identity);
                this->dPropagatedStatedIndependentVariable = math::Matrix<double>(10, 2, 0.0);
                this->StateBeforeEventBeforePropagation = math::Matrix<doubleType>(10, 1, 0.0);
                this->state_before_event = math::Matrix<doubleType>(10, 1, 0.0);
            }

            //set up propagation
            if (this->AllowStateToPropagate)
            {
                this->dPropagationTime_dIndependentVariable = 1.0;

                //instantiate a propagator
                if (this->myPropagatorType == PropagatorType::KeplerPropagator)
                {

                    this->myPropagator = Astrodynamics::CreatePropagator(myOptions,
                        Universe,
                        6,
                        this->StateBeforeEventBeforePropagation,
                        this->state_before_event,
                        this->STM,
                        this->dPropagatedStatedIndependentVariable,
                        &this->dPropagationTime_dIndependentVariable);
                }
                else
                {

                    //acceleration model
                    this->mySpacecraftAccelerationModel = new Astrodynamics::SpacecraftAccelerationModel(myOptions,
                        &myOptions->Journeys[journeyIndex],
                        Universe,
                        this->Xdescriptions,
                        mySpacecraft,
                        11); //STM size
                    this->mySpacecraftAccelerationModel->setDutyCycle(1.0);

                    //EOM
                    this->myEOM = Astrodynamics::TimeDomainSpacecraftEOM();
                    this->myEOM.setSpacecraftAccelerationModel(this->mySpacecraftAccelerationModel);

                    //integration scheme
                    this->myIntegrationScheme = CreateIntegrationScheme(&this->myEOM, 10, 11);

                    //propagator
                    this->myPropagator = Astrodynamics::CreatePropagator(myOptions,
                        Universe,
                        10,
                        11,
                        this->StateBeforeEventBeforePropagation,
                        this->state_before_event,
                        this->STM,
                        this->dPropagatedStatedIndependentVariable,
                        (Integration::Integrand*) &this->myEOM,
                        this->myIntegrationScheme,
                        &this->dPropagationTime_dIndependentVariable,
                        myOptions->Journeys[journeyIndex].override_integration_step_size //can't use class members because they don't get initialized until later
                            ? myOptions->Journeys[journeyIndex].integration_step_size
                            : myOptions->integration_time_step_size);
                }

            }
        }//end initialize()

        FreePointBoundary::~FreePointBoundary()
        {
            delete this->myStateRepresentation;

            if (this->AllowStateToPropagate)
            {
                delete this->myPropagator;

                if (this->myPropagatorType == PropagatorType::IntegratedPropagator)
                {
                    //delete force modely things
                    delete this->mySpacecraftAccelerationModel;
                    delete this->myIntegrationScheme;
                }
            }
        }

        //method to handle bounds for event left state
        //state is encoded in cartesian coordinates
        void FreePointBoundary::calcbounds_event_left_side(const std::vector< std::tuple<double, double> >& StateBounds, std::vector<size_t> timeVariables)
        {
            //Step 1: set the current stage
            this->mySpacecraft->setActiveStage(this->stageIndex);

            //Step 2: first six decision variables and sparsity entries
            std::vector< std::tuple<std::string, double, bool> > statesToRepresent; //state name, scale factor, isAngle
            switch (this->myStateRepresentationEnum)
            {
                case StateRepresentation::SphericalRADEC:
                {
                    statesToRepresent.push_back({ "r",  this->myUniverse->LU, false });
                    statesToRepresent.push_back({ "RA",  1.0, true });
                    statesToRepresent.push_back({ "DEC", 1.0, true });
                    statesToRepresent.push_back({ "v", this->myUniverse->LU / this->myUniverse->TU, false });
                    statesToRepresent.push_back({ "vRA",  1.0, true });
                    statesToRepresent.push_back({ "vDEC", 1.0, true });

                    break;
                }
                case StateRepresentation::SphericalAZFPA:
                {
                    statesToRepresent.push_back({ "r",  this->myUniverse->LU, false });
                    statesToRepresent.push_back({ "RA", 1.0, true });
                    statesToRepresent.push_back({ "DEC", 1.0, true });
                    statesToRepresent.push_back({ "v", this->myUniverse->LU / this->myUniverse->TU, false });
                    statesToRepresent.push_back({ "AZ", 1.0, true });
                    statesToRepresent.push_back({ "FPA", 1.0, true });

                    break;
                }
                case StateRepresentation::Cartesian:
                {
                    statesToRepresent.push_back({ "x",  this->myUniverse->LU, false });
                    statesToRepresent.push_back({ "y",  this->myUniverse->LU, false });
                    statesToRepresent.push_back({ "z",  this->myUniverse->LU, false });
                    statesToRepresent.push_back({ "vx", this->myUniverse->LU / this->myUniverse->TU, false });
                    statesToRepresent.push_back({ "vy", this->myUniverse->LU / this->myUniverse->TU, false });
                    statesToRepresent.push_back({ "vz", this->myUniverse->LU / this->myUniverse->TU, false });

                    break;
                }
                case StateRepresentation::COE:
                {
                    statesToRepresent.push_back({ "SMA",  this->myUniverse->LU, false });
                    statesToRepresent.push_back({ "ECC",  1.0, false });
                    statesToRepresent.push_back({ "INC",  1.0, true });
                    statesToRepresent.push_back({ "RAAN", 1.0, true });
                    statesToRepresent.push_back({ "AOP", 1.0, true });
                    statesToRepresent.push_back({ "TA",  1.0, true });

                    break;
                }
                case StateRepresentation::MEE:
                {
                    statesToRepresent.push_back({ "P",  this->myUniverse->LU, false });
                    statesToRepresent.push_back({ "F",  1.0, false });
                    statesToRepresent.push_back({ "G",  1.0, false });
                    statesToRepresent.push_back({ "H",  1.0, false });
                    statesToRepresent.push_back({ "K",  1.0, false });
                    statesToRepresent.push_back({ "L",  1.0, true });

                    break;
                }
                case StateRepresentation::IncomingBplane:
                {
                    statesToRepresent.push_back({ "VINFin",  this->myUniverse->LU / this->myUniverse->TU, false });
                    statesToRepresent.push_back({ "RHAin",  1.0, true });
                    statesToRepresent.push_back({ "DHAin",  1.0, true });
                    statesToRepresent.push_back({ "BRADIUSin", this->myUniverse->LU, false });
                    statesToRepresent.push_back({ "BTHETAin", 1.0, true });
                    statesToRepresent.push_back({ "TAin",  1.0, true });

                    break;
                }
                case StateRepresentation::OutgoingBplane:
                {
                    statesToRepresent.push_back({ "VINFout",  this->myUniverse->LU / this->myUniverse->TU, false });
                    statesToRepresent.push_back({ "RHAout",  1.0, true });
                    statesToRepresent.push_back({ "DHAout",  1.0, true });
                    statesToRepresent.push_back({ "BRADIUSout", this->myUniverse->LU, false });
                    statesToRepresent.push_back({ "BTHETAout", 1.0, true });
                    statesToRepresent.push_back({ "TAout",  1.0, true });

                    break;
                }
                default:
                    throw std::invalid_argument("FreePointBoundary does not recognize state representation " + std::to_string(this->myStateRepresentationEnum) + ". Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + ".");
            }

            for (size_t stateIndex = 0; stateIndex < 7; ++stateIndex)
            {
                //is this a mass?
                if (stateIndex == 6)
                {
                    if (this->isFirstEventInMission && !this->myOptions->allow_initial_mass_to_vary)
                        this->Xlowerbounds->push_back(std::get<1>(StateBounds[6]) - math::SMALL);
                    else
                        this->Xlowerbounds->push_back(std::get<0>(StateBounds[6]));
                    this->Xupperbounds->push_back(std::get<1>(StateBounds[6]));

                    X_scale_factors->push_back(1.0 / this->myUniverse->continuity_constraint_scale_factors(stateIndex));
                    Xdescriptions->push_back(prefix + "event left state " + stateNames[stateIndex]);
                }
                else
                {
                    //is this an angle? if so, transform the bounds from degrees to radians
                    if (std::get<2>(statesToRepresent[stateIndex]) == true) 
                    {
                        Xlowerbounds->push_back(std::get<0>(StateBounds[stateIndex]) * math::PI / 180.0);
                        Xupperbounds->push_back(std::get<1>(StateBounds[stateIndex]) * math::PI / 180.0);
                    }
                    else
                    {
                        Xlowerbounds->push_back(std::get<0>(StateBounds[stateIndex]));
                        Xupperbounds->push_back(std::get<1>(StateBounds[stateIndex]));
                    }
                    X_scale_factors->push_back(std::get<1>(statesToRepresent[stateIndex]));
                    Xdescriptions->push_back(prefix + "event left state " + std::get<0>(statesToRepresent[stateIndex]));
                }
                this->Xindex_encoded_state.push_back(this->Xdescriptions->size() - 1);

                //in a FreePointBoundary, all state have a dependence on all other states and we need to track them
                for (size_t dependentStateIndex = 0; dependentStateIndex < 7; ++dependentStateIndex)
                {
                    this->Derivatives_of_StateBeforeEvent.push_back({ Xdescriptions->size() - 1, dependentStateIndex, 1.0 });
                    this->dIndex_StateBeforeEvent_wrt_EncodedState[dependentStateIndex].push_back(this->Derivatives_of_StateBeforeEvent.size() - 1);
                }
            }//end loop over 6-state variables

            this->dIndex_mass_wrt_encodedMass = this->dIndex_StateBeforeEvent_wrt_EncodedState[6][6];

            //Step 3: left epoch
            if (this->isFirstEventInMission)
            {
                this->Xlowerbounds->push_back(std::get<0>(StateBounds[7]));
                this->Xupperbounds->push_back(std::get<1>(StateBounds[7]));
                this->X_scale_factors->push_back(1.0 / this->myUniverse->continuity_constraint_scale_factors(7));
                this->Xdescriptions->push_back(prefix + "event left state epoch");
                this->Xindex_encoded_state.push_back(this->Xdescriptions->size() - 1);
                timeVariables.insert(timeVariables.begin(), this->Xdescriptions->size() - 1);
            }

            this->calculate_dependencies_left_epoch(timeVariables);

            //Step 4: If we allow the state to propagate
            if (this->AllowStateToPropagate
                || this->myEncodedReferenceFrame == ReferenceFrame::J2000_BCF
                || this->myEncodedReferenceFrame == ReferenceFrame::TrueOfDate_BCI
                || this->myEncodedReferenceFrame == ReferenceFrame::TrueOfDate_BCF
                || this->myEncodedReferenceFrame == ReferenceFrame::SAM
                || this->myEncodedReferenceFrame == ReferenceFrame::ObjectReferenced)
            {
                //all state variables, including mass in an FreePointBoundary event have a derivative with respect to epoch
                //we'll put in a dummy derivative of 0.0 for now, and later, when the event is processed, we'll do it right
                for (size_t stateIndex = 0; stateIndex < 7; ++stateIndex)
                {
                    std::vector<size_t> dIndex_StateBeforeEvent_i_wrt_Time;
                    for (size_t Xepoch_index = 0; Xepoch_index < this->Xindices_EventLeftEpoch.size(); ++Xepoch_index)
                    {
                        this->Derivatives_of_StateBeforeEvent_wrt_Time.push_back({this->Xindices_EventLeftEpoch[Xepoch_index], stateIndex, 0.0} );
                        dIndex_StateBeforeEvent_i_wrt_Time.push_back(this->Derivatives_of_StateBeforeEvent_wrt_Time.size() - 1);
                    }
                    this->dIndex_StateBeforeEvent_wrt_Time.push_back(dIndex_StateBeforeEvent_i_wrt_Time);
                }//end time derivatives, only relevant if state is allowed to propagate
            }//end derivatives due to propagating the state
        }//end calcbounds_event_left_side()

        void FreePointBoundary::calcbounds_event_right_side()
        {
            this->Derivatives_of_StateAfterEvent = this->Derivatives_of_StateBeforeEvent;
            this->Derivatives_of_StateAfterEvent_wrt_Time = this->Derivatives_of_StateBeforeEvent_wrt_Time;

            this->Xindices_EventRightEpoch = this->Xindices_EventLeftEpoch;
        }//end calcbounds_event_right_side

        //process methods
        void FreePointBoundary::process_event_left_side(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG)
        {
            //Step 1: set the active stage
            this->mySpacecraft->setActiveStage(this->stageIndex);

            //Step 2: extract state - position, velocity, and mass
            for (size_t stateIndex = 0; stateIndex < 7; ++stateIndex)
            {
                this->StateBeforeEventEncoded(stateIndex) = X[Xindex++];
            }
            this->state_before_event(6) = this->StateBeforeEventEncoded(6);

            //Step 2.1: convert to Cartesian, store the derivatives
            //Step 2.1.1: convert
            math::Matrix<doubleType> stateBeforeEventCartesian = this->myStateRepresentation->convertFromRepresentationToCartesian(this->StateBeforeEventEncoded.getSubMatrix1D(0, 5), needG);

            //Step 2.1.2: put in state container - overwrite the original encoded state
            for (size_t stateIndex : {0, 1, 2, 3, 4, 5})
            {
                this->StateBeforeEventEncoded(stateIndex) = stateBeforeEventCartesian(stateIndex);
            }
            
            //Step 2.1.3: derivatives
            if (needG)
            {
                //Step 2.1.3.1: extract the transformation matrix
                math::Matrix<doubleType> TransformationMatrix = this->myStateRepresentation->getRepresentationToCartesianTransitionMatrix();

                //Step 2.1.3.2: insert into the derivative tuples
                for (size_t nativeStateIndex : {0, 1, 2, 3, 4, 5})
                {
                    for (size_t cartesianStateIndex : {0, 1, 2, 3, 4, 5})
                    {
                        this->dStateBeforeEventEncoded_dEncodedElements(cartesianStateIndex, nativeStateIndex) = TransformationMatrix(cartesianStateIndex, nativeStateIndex);
                    }
                }
            }//end derivatives

            //Step 3: figure out the current epoch
            this->process_left_epoch(X, Xindex, F, Findex, G, needG);
            this->StateBeforeEventEncoded(7) = this->state_before_event(7);

            //Step 4: frame rotation
            //we have to rotate frames if the state was not supplied in ICRF
            if (!(this->myEncodedReferenceFrame == ReferenceFrame::ICRF))
            {
                math::Matrix<doubleType> dstate_dt(3, 1, 0.0);
                math::Matrix<doubleType> referenceVector(3, 1, 0.0);
                math::Matrix<doubleType> dreferenceVector_dt(3, 1, 0.0);
                math::Matrix<doubleType> Rotated3Vector;
                math::Matrix<doubleType> dRotated3Vector_dt;
                math::Matrix<doubleType> dRotated3Vector_dOriginal3Vector;
                math::Matrix<doubleType> dRotated3Vector_dreferenceVector;

                if (this->myEncodedReferenceFrame == ReferenceFrame::SAM)
                {
                    //we need the vector to the sun and its time derivative
                    doubleType centralBodyState[12];
                    this->myUniverse->locate_central_body(this->EventLeftEpoch,
                        centralBodyState,
                        *this->myOptions,
                        needG);

                    //set vector *to* Sun
                    for (size_t stateIndex : {0, 1, 2})
                    {
                        referenceVector(stateIndex) = -centralBodyState[stateIndex];
                        dreferenceVector_dt(stateIndex) = -centralBodyState[stateIndex + 6];
                    }

                    //rotate the vectors into the body's J2000BCI frame
                    math::Matrix<doubleType> rot_in = referenceVector;
                    
                    this->myUniverse->LocalFrame.rotate_frame_to_frame(ReferenceFrame::ICRF,
                        rot_in,
                        ReferenceFrame::J2000_BCI,
                        referenceVector,
                        this->EventLeftEpoch);

                    rot_in = dreferenceVector_dt;

                    this->myUniverse->LocalFrame.rotate_frame_to_frame(ReferenceFrame::ICRF,
                        rot_in,
                        ReferenceFrame::J2000_BCI,
                        dreferenceVector_dt,
                        this->EventLeftEpoch);
                }
                else if (this->myEncodedReferenceFrame == ReferenceFrame::ObjectReferenced)
                {
                    //we need the vector to the refernce body and its time derivative
                    referenceVector.resize(6, 1, 0.0);
                    dreferenceVector_dt.resize(6, 1, 0.0);

                    doubleType bodyState[12];
                    this->myBody->locate_body(this->EventLeftEpoch,
                        bodyState,
                        needG,
                        *this->myOptions);

                    //set vector *to* reference body
                    for (size_t stateIndex : {0, 1, 2, 3, 4, 5})
                    {
                        referenceVector(stateIndex) = bodyState[stateIndex];
                        dreferenceVector_dt(stateIndex) = bodyState[stateIndex + 6];
                    }
                }

                this->myUniverse->LocalFrame.rotate_frame_to_frame(this->myEncodedReferenceFrame,//OriginalReferenceFrame
                    this->StateBeforeEventEncoded.getSubMatrix1D(0, 2),//state
                    dstate_dt,//dstate_dt - this state is ENCODED, so this is always zero
                    referenceVector,//referenceVector,
                    dreferenceVector_dt,//dreferenceVector_dt
                    ReferenceFrame::ICRF,//RotatedReferenceFrame
                    Rotated3Vector,//Rstate
                    dRotated3Vector_dOriginal3Vector,//dRstate_dstate
                    dRotated3Vector_dt,//dRstate_dt
                    dRotated3Vector_dreferenceVector,//dRstate_dreferenceVector
                    this->EventLeftEpoch,//ReferenceEpochJD
                    needG);//GenerateDerivatives

                for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
                    this->state_before_event(stateIndex) = Rotated3Vector(stateIndex);

                if (needG)
                {
                    for (size_t i = 0; i < 3; ++i)
                    {
                        for (size_t j = 0; j < 3; ++j)
                        {
                            this->dStateBeforeEventICRF_dStateBeforeEventEncoded(i, j) = dRotated3Vector_dOriginal3Vector(i, j);
                        }
                    }

                    //time derivative of position
                    for (size_t i = 0; i < 3; ++i)
                    {
                        this->dStateBeforeEventICRF_dStateBeforeEventEncoded(i, 7) = dRotated3Vector_dt(i) _GETVALUE;
                    }
                }

                this->myUniverse->LocalFrame.rotate_frame_to_frame(this->myEncodedReferenceFrame,//OriginalReferenceFrame
                    this->StateBeforeEventEncoded.getSubMatrix1D(3, 5),//state
                    dstate_dt,//dstate_dt - this state is ENCODED, so this is always zero
                    referenceVector,//referenceVector,
                    dreferenceVector_dt,//dreferenceVector_dt
                    ReferenceFrame::ICRF,//RotatedReferenceFrame
                    Rotated3Vector,//Rstate
                    dRotated3Vector_dOriginal3Vector,//dRstate_dstate
                    dRotated3Vector_dt,//dRstate_dt
                    dRotated3Vector_dreferenceVector,//dRstate_dreferenceVector
                    this->EventLeftEpoch,//ReferenceEpochJD
                    needG);//GenerateDerivatives

                for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
                    this->state_before_event(stateIndex + 3) = Rotated3Vector(stateIndex);
                               
                if (needG)
                {
                    //derivative of velocity in w.r.t. velocity out
                    for (size_t i = 0; i < 3; ++i)
                    {
                        for (size_t j = 0; j < 3; ++j)
                        {
                            this->dStateBeforeEventICRF_dStateBeforeEventEncoded(i + 3, j + 3) = dRotated3Vector_dOriginal3Vector(i, j);
                        }
                    }

                    //time derivative of velocity
                    for (size_t i = 0; i < 3; ++i)
                    {
                        this->dStateBeforeEventICRF_dStateBeforeEventEncoded(i + 3, 7) = dRotated3Vector_dt(i) _GETVALUE;
                    }
                }
            }//end rotation
            else
                this->state_before_event.shallow_copy(this->StateBeforeEventEncoded, 8);

            //Step 5: state propagation
            //Step 5.1: first we need to determine the propagation time
            doubleType PropagationTime = this->state_before_event(7) - this->ReferenceEpoch + 1.0e-13;
            if (this->AllowStateToPropagate)
            {
                //Step 5.2: copy the state to a temporary vector
                this->StateBeforeEventBeforePropagation.shallow_copy(this->state_before_event, 8);

                //Step 5.3: then we need to propagate
                this->total_number_of_states_to_integrate = needG
                    ? 11 + 11 * 11
                    : 11;
                this->dPropagatedStatedIndependentVariable.assign_zeros();
                this->myPropagator->setCurrentEpoch(this->ReferenceEpoch);
                this->myPropagator->setIndexOfEpochInStateVec(7);
                this->myPropagator->setCurrentIndependentVariable(this->ReferenceEpoch);
                this->myPropagator->propagate(PropagationTime, needG);
                this->state_before_event(7) = this->StateBeforeEventBeforePropagation(7);

                if (this->myPropagatorType == PropagatorType::IntegratedPropagator)
                {
                    this->state_before_event(8) = 0.0; //chemical fuel - we're not actually moving the spacecraft, we're redefining a boundary point. So tankage does not change.
                    this->state_before_event(9) = 0.0; //chemical oxidizer
                }
            }

            if (needG)
            {
                //Step 6: form the augmented STM - we only need to do this if propagating, otherwise the matrix is identity
                size_t stateMax = this->myPropagatorType == PropagatorType::IntegratedPropagator ? 7 : 6;
                for (size_t i = 0; i < stateMax; ++i)
                {
                    for (size_t j = 0; j < stateMax; ++j)
                        this->dStateBeforeEvent_dStateBeforeEventICRF(i, j) = this->STM(i, j);

                    if (this->myPropagatorType == PropagatorType::KeplerPropagator)
                        this->dStateBeforeEvent_dStateBeforeEventICRF(i, 7) = (PropagationTime >= 0.0 ? 1.0 : -1.0) * this->dPropagatedStatedIndependentVariable(i);
                    else
                        this->dStateBeforeEvent_dStateBeforeEventICRF(i, 7) = this->STM(i, 10); //remember, the initial epoch is the propagation time, too
                }

                //Step 7: assemble the derivatives of state_before_event with respect to the encoded state
                this->dStateBeforeEvent_dStateBeforeEventEncoded = this->dStateBeforeEvent_dStateBeforeEventICRF * this->dStateBeforeEventICRF_dStateBeforeEventEncoded * this->dStateBeforeEventEncoded_dEncodedElements;

                //Step 8: insert the derivatives into the appropriate tuples
                for (size_t encodedStateIndex = 0; encodedStateIndex < 7; ++encodedStateIndex)
                {
                    for (size_t stateIndex = 0; stateIndex < 7; ++stateIndex)
                    {
                        std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_StateBeforeEvent_wrt_EncodedState[stateIndex][encodedStateIndex]]) = this->dStateBeforeEvent_dStateBeforeEventEncoded(stateIndex, encodedStateIndex) _GETVALUE;
                    }
                }

                //Step 9: special things if the state is allowed to propagate or we are using a time-dependent frame
                if (this->AllowStateToPropagate
                    || this->myEncodedReferenceFrame == ReferenceFrame::J2000_BCF
                    || this->myEncodedReferenceFrame == ReferenceFrame::TrueOfDate_BCI
                    || this->myEncodedReferenceFrame == ReferenceFrame::TrueOfDate_BCF
                    || this->myEncodedReferenceFrame == ReferenceFrame::SAM
                    || this->myEncodedReferenceFrame == ReferenceFrame::ObjectReferenced)
                {
                    //put time derivatives into tuples
                    //all state variables except mass in an EphemerisPeggedBoundary event have a derivative with respect to epoch
                    for (size_t stateIndex = 0; stateIndex < 7; ++stateIndex)
                    {
                        for (size_t Xepoch_index = 0; Xepoch_index < this->Xindices_EventLeftEpoch.size(); ++Xepoch_index)
                        {
                            std::get<2>(this->Derivatives_of_StateBeforeEvent_wrt_Time[this->dIndex_StateBeforeEvent_wrt_Time[stateIndex][Xepoch_index]]) = this->dStateBeforeEvent_dStateBeforeEventEncoded(stateIndex, 7) _GETVALUE;
                        }
                    }
                }//end time derivatives, only relevant if state is allowed to propagate
            }//end derivative stuff

            //set the boundary state, which we need for printing
            this->boundary_state.shallow_copy(this->state_before_event, 8);
        }//end process_event_left_side()

        void FreePointBoundary::process_event_right_side(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            this->state_after_event.shallow_copy(this->state_before_event, 8);
            
            for (size_t dIndex = 0; dIndex < this->Derivatives_of_StateBeforeEvent.size(); ++dIndex)
                this->Derivatives_of_StateAfterEvent[dIndex] = this->Derivatives_of_StateBeforeEvent[dIndex];

            for (size_t dIndex = 0; dIndex < this->Derivatives_of_StateBeforeEvent_wrt_Time.size(); ++dIndex)
                this->Derivatives_of_StateAfterEvent_wrt_Time[dIndex] = this->Derivatives_of_StateBeforeEvent_wrt_Time[dIndex];

            this->EventRightEpoch = this->EventLeftEpoch;
        }//end process_event_right_side()
    }//close namespace events
}//close namespace EMTG