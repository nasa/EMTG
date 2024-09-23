
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

#include "PeriapseBoundary.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        PeriapseBoundary::PeriapseBoundary() :
            FreePointBoundary::FreePointBoundary()
        {}

        PeriapseBoundary::PeriapseBoundary(const std::string& name,
            const size_t& journeyIndex,
            const size_t& phaseIndex,
            size_t& stageIndex,
            Astrodynamics::universe* Universe,
            HardwareModels::Spacecraft* mySpacecraft,
            missionoptions* myOptions) :
            PeriapseBoundary()

        {
            this->initialize(name,
                journeyIndex,
                phaseIndex,
                stageIndex,
                Universe,
                mySpacecraft,
                myOptions);
        }//end constructor

        void PeriapseBoundary::initialize(const std::string& name,
            const size_t& journeyIndex,
            const size_t& phaseIndex,
            size_t& stageIndex,
            Astrodynamics::universe* Universe,
            HardwareModels::Spacecraft* mySpacecraft,
            missionoptions* myOptions)
        {
            //periapse states don't get propagated
            this->AllowStateToPropagate = false;

            //periapse boundary states are always encoded in ICRF
            this->myEncodedReferenceFrame = ReferenceFrame::ICRF;

            //do we need the rdotv constraint?
            if (this->myStateRepresentationEnum == StateRepresentation::IncomingBplane
                || this->myStateRepresentationEnum == StateRepresentation::OutgoingBplane
                || this->myStateRepresentationEnum == StateRepresentation::COE
                || this->myStateRepresentationEnum == StateRepresentation::SphericalAZFPA)
            {
                this->need_rdotv_constraint = false;
            }
            else
            {
                this->need_rdotv_constraint = true;
            }

            //do we need the distance constraint?
            if (this->myStateRepresentationEnum == StateRepresentation::SphericalAZFPA
                || this->myStateRepresentationEnum == StateRepresentation::SphericalRADEC)
            {
                this->need_distance_constraint = false;
            }
            else
            {
                this->need_distance_constraint = true;
            }

            //base class
            this->FreePointBoundary::initialize(name,
                journeyIndex,
                phaseIndex,
                stageIndex,
                Universe,
                mySpacecraft,
                myOptions);

            this->LeftBoundaryIsABody = false;
        }//end initialize


        //******************************************calcbounds methods
        void PeriapseBoundary::calcbounds_event_left_side(const std::vector<double>& RadiusBounds,
                                                          const std::vector<double>& VelocityMagnitudeBounds,
                                                          std::vector<size_t> timeVariables)
        {
            //Step 0: what should the state bounds be?
            std::vector< std::tuple<double, double> > stateBounds; //in degrees because FreePointBoundary converts
            switch (this->myStateRepresentationEnum)
            {
                case StateRepresentation::SphericalRADEC:
                {
                    //radius
                    stateBounds.push_back({ RadiusBounds[0], RadiusBounds[1] });
                    //RA
                    stateBounds.push_back({ -360.0 * 4, 360.0 * 4 });
                    //DEC
                    stateBounds.push_back({ -90, 90.0 });
                    //velocity magnitude
                    stateBounds.push_back({ VelocityMagnitudeBounds[0], VelocityMagnitudeBounds[1] });
                    //vRA
                    stateBounds.push_back({ -360.0 * 4, 360.0 * 4 });
                    //vDEC
                    stateBounds.push_back({ -90.0, 90.0 });

                    break;
                }
                case StateRepresentation::SphericalAZFPA:
                {
                    //radius
                    stateBounds.push_back({ RadiusBounds[0], RadiusBounds[1] });
                    //RA
                    stateBounds.push_back({ -360.0 * 4, 360.0 * 4 });
                    //DEC
                    stateBounds.push_back({ -90, 90.0 });
                    //velocity magnitude
                    stateBounds.push_back({ VelocityMagnitudeBounds[0], VelocityMagnitudeBounds[1] });
                    //AZ
                    stateBounds.push_back({ -360.0 * 4, 360.0 * 4 });
                    //FPA
                    stateBounds.push_back({ 90.0 - math::SMALL, 90.0 + math::SMALL });

                    break;
                }
                case StateRepresentation::Cartesian:
                {
                    //x
                    stateBounds.push_back({ -RadiusBounds[1], RadiusBounds[1] });
                    //y
                    stateBounds.push_back({ -RadiusBounds[1], RadiusBounds[1] });
                    //z
                    stateBounds.push_back({ -RadiusBounds[1], RadiusBounds[1] });
                    //vx
                    stateBounds.push_back({ -VelocityMagnitudeBounds[1], VelocityMagnitudeBounds[1] });
                    //vy
                    stateBounds.push_back({ -VelocityMagnitudeBounds[1], VelocityMagnitudeBounds[1] });
                    //vz
                    stateBounds.push_back({ -VelocityMagnitudeBounds[1], VelocityMagnitudeBounds[1] });

                    break;
                }
                case StateRepresentation::COE:
                {
                    //SMA
                    stateBounds.push_back({ -RadiusBounds[1], RadiusBounds[1] });
                    //ECC
                    stateBounds.push_back({ math::SMALL, 10.0 });
                    //INC
                    stateBounds.push_back({ -180.0, 180.0 });
                    //RAAN
                    stateBounds.push_back({ -360.0 * 4, 360.0 * 4 });
                    //AOP
                    stateBounds.push_back({ -360.0 * 4, 360.0 * 4 });
                    //TA
                    stateBounds.push_back({ -math::SMALL, math::SMALL });

                    break;
                }
                case StateRepresentation::MEE:
                {
                    //P
                    stateBounds.push_back({ 0.0, RadiusBounds[1] });
                    //F
                    stateBounds.push_back({ -10.0, 10.0 });
                    //G
                    stateBounds.push_back({ -10.0, 10.0 });
                    //H
                    stateBounds.push_back({ -10.0, 10.0 });
                    //K
                    stateBounds.push_back({ -10.0, 10.0 });
                    //L
                    stateBounds.push_back({ -360.0 * 4, 360.0 * 4 });
                                       
                    break;
                }
                case StateRepresentation::IncomingBplane:
                {
                    //VINF
                    stateBounds.push_back({ VelocityMagnitudeBounds[0], VelocityMagnitudeBounds[1] });
                    //RHA
                    stateBounds.push_back({ -360.0 * 4, 360.0 * 4 });
                    //DHA
                    stateBounds.push_back({ -90.0, 90.0 });
                    //BRADIUS
                    stateBounds.push_back({ math::SMALL, this->myUniverse->r_SOI });
                    //BTHETA
                    stateBounds.push_back({ -360.0 * 4, 360.0 * 4 });
                    //TA
                    stateBounds.push_back({ -math::SMALL, math::SMALL });

                    break;
                }
                case StateRepresentation::OutgoingBplane:
                {
                    //VINF
                    stateBounds.push_back({ VelocityMagnitudeBounds[0], VelocityMagnitudeBounds[1] });
                    //RHA
                    stateBounds.push_back({ -360.0 * 4, 360.0 * 4 });
                    //DHA
                    stateBounds.push_back({ -90.0, 90.0 });
                    //BRADIUS
                    stateBounds.push_back({ math::SMALL, this->myUniverse->r_SOI });
                    //BTHETA
                    stateBounds.push_back({ -360.0 * 4, 360.0 * 4 });
                    //TA
                    stateBounds.push_back({ -math::SMALL, math::SMALL });

                    break;
                }
                default:
                    throw std::invalid_argument("PeriapseBoundary does not recognize state representation " + std::to_string(this->myStateRepresentationEnum) + ". Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + ".");
            }
            //mass
            stateBounds.push_back({ math::SMALL, this->myJourneyOptions->maximum_mass });

            //Step 1: set state bounds by calling FreePointBoundary
            this->FreePointBoundary::calcbounds_event_left_side(stateBounds, timeVariables);

            //Step 2: r dot v = 0 constraint. This is not needed for state representations that encode true anomaly or flight path angle
            if (this->need_rdotv_constraint)
            {
                this->Flowerbounds->push_back(-1.0e-13);
                this->Fupperbounds->push_back(1.0e-13);
                this->Fdescriptions->push_back(this->prefix + "periapse is an apse");

                //for the purpose of solver robustness, we have several different versions of this constraint for different state representations
                if (this->myStateRepresentationEnum == StateRepresentation::SphericalRADEC)
                {
                    //in order, derivatives wrt RA, DEC, vRA, vDEC
                    for (size_t stateIndex : {1, 2, 4, 5})
                    {
                        this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                            this->Xindex_encoded_state[stateIndex],
                            this->Gindices_rdotv);
                    }
                }//end derivatives for rdotv constraint, SphericalRADEC
                else //the generic case
                {
                    for (size_t dIndex = 0; dIndex < this->Derivatives_of_StateBeforeEvent.size(); ++dIndex)
                    {
                        std::tuple<size_t, size_t, double>& derivativeEntry = this->Derivatives_of_StateBeforeEvent[dIndex];

                        size_t stateIndex = std::get<1>(derivativeEntry);

                        if (stateIndex < 6) //6-state
                        {
                            size_t Xindex = std::get<0>(derivativeEntry);
                            this->dIndices_rdotv.push_back(dIndex);

                            this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                                Xindex,
                                this->Gindices_rdotv);
                        }
                    }
                }//end derivatives for rdotv constraint, generic case
            }//end RdotV constraint

            //Step 3: distance from central body constraint. This is not needed for state representations that encode distance from central body
            if (this->need_distance_constraint)
            {
                this->Flowerbounds->push_back(this->periapseDistanceBounds[0] / this->myUniverse->LU);
                this->Fupperbounds->push_back(this->periapseDistanceBounds[1] / this->myUniverse->LU);
                this->Fdescriptions->push_back(this->prefix + "distance from central body");

                for (size_t dIndex = 0; dIndex < this->Derivatives_of_StateBeforeEvent.size(); ++dIndex)
                {
                    std::tuple<size_t, size_t, double>& derivativeEntry = this->Derivatives_of_StateBeforeEvent[dIndex];

                    size_t stateIndex = std::get<1>(derivativeEntry);

                    if (stateIndex < 3) //position
                    {
                        size_t Xindex = std::get<0>(derivativeEntry);
                        this->dIndices_distance_constraint.push_back(dIndex);

                        this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                            Xindex,
                            this->Gindices_distance_constraint);
                    }
                }//end derivatives for distance constraint
            }//end distance constraint
        }//end calcbounds_left_side()

        void PeriapseBoundary::calcbounds_event_right_side()
        {
            //base class
            this->FreePointBoundary::calcbounds_event_right_side();
        }//end calcbounds_right_side()

        //******************************************process methods

        void PeriapseBoundary::process_event_left_side(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //Step 1: the state vector itself is handled by the base class
            this->FreePointBoundary::process_event_left_side(X, Xindex, F, Findex, G, needG);

            //Step 2: rdotv constraint
            if (this->need_rdotv_constraint)
            {
                //for the purpose of solver robustness, we have several different versions of this constraint for different state representations
                if (this->myStateRepresentationEnum == StateRepresentation::SphericalRADEC)
                {
                    const doubleType& RA = X[this->Xindex_encoded_state[1]];
                    const doubleType& DEC = X[this->Xindex_encoded_state[2]];
                    const doubleType& vRA = X[this->Xindex_encoded_state[4]];
                    const doubleType& vDEC = X[this->Xindex_encoded_state[5]];

                    doubleType cosRA = cos(RA);
                    doubleType sinRA = sin(RA);
                    doubleType cosDEC = cos(DEC);
                    doubleType sinDEC = sin(DEC);
                    doubleType cosvRA = cos(vRA);
                    doubleType sinvRA = sin(vRA);
                    doubleType cosvDEC = cos(vDEC);
                    doubleType sinvDEC = sin(vDEC);

                    doubleType rdotv = sinDEC * sinvDEC + cosDEC * cosRA*cosvDEC*cosvRA + cosDEC * sinRA*cosvDEC*sinvRA;
                    F[Findex++] = rdotv;

                    if (needG)
                    {
                        //r dot v constraint
                        double drdotv_dRA = (cosDEC*cosRA*cosvDEC*sinvRA - cosDEC * sinRA*cosvDEC*cosvRA) _GETVALUE;
                        double drdotv_dDEC = -(cosRA*sinDEC*cosvDEC*cosvRA - cosDEC * sinvDEC + sinDEC * sinRA*cosvDEC*sinvRA)_GETVALUE;
                        double drdotv_dvRA = -(cosDEC*cosRA*cosvDEC*sinvRA - cosDEC * sinRA*cosvDEC*cosvRA) _GETVALUE;
                        double drdotv_dvDEC = -(cosDEC*cosRA*cosvRA*sinvDEC - sinDEC * cosvDEC + cosDEC * sinRA*sinvDEC*sinvRA) _GETVALUE;

                        G[this->Gindices_rdotv[0]] = this->X_scale_factors->operator[](this->Xindex_encoded_state[1]) * drdotv_dRA;
                        G[this->Gindices_rdotv[1]] = this->X_scale_factors->operator[](this->Xindex_encoded_state[2]) * drdotv_dDEC;
                        G[this->Gindices_rdotv[2]] = this->X_scale_factors->operator[](this->Xindex_encoded_state[4]) * drdotv_dvRA;
                        G[this->Gindices_rdotv[3]] = this->X_scale_factors->operator[](this->Xindex_encoded_state[5]) * drdotv_dvDEC;
                    }//end derivatives, SphericalRADEC
                }//end rdotv constraint, SphericalRADEC
                else
                {
                    doubleType rdotv = this->state_before_event.getSubMatrix1D(0, 2).dot(this->state_before_event.getSubMatrix1D(3, 5));

                    F[Findex++] = rdotv / this->myUniverse->LU / this->myUniverse->LU * this->myUniverse->TU;

                    if (needG)
                    {
                        //first zero out the derivatives
                        for (size_t Gindex : this->Gindices_rdotv)
                        {
                            G[Gindex] = 0.0;
                        }

                        //now assign the entries
                        for (size_t entryIndex = 0; entryIndex < this->dIndices_rdotv.size(); ++entryIndex)
                        {
                            size_t dIndex = this->dIndices_rdotv[entryIndex];

                            std::tuple<size_t, size_t, double>& derivativeEntry = this->Derivatives_of_StateBeforeEvent[dIndex];

                            size_t stateIndex = std::get<1>(derivativeEntry);

                            size_t Xindex = std::get<0>(derivativeEntry);

                            size_t Gindex = this->Gindices_rdotv[entryIndex];

                            if (stateIndex < 3) //position
                            {
                                //cartesian component
                                double drdotv_dCartesianState = this->state_before_event(stateIndex + 3) _GETVALUE;

                                //chain to the actual state representation
                                double TheDerivative = drdotv_dCartesianState * std::get<2>(this->Derivatives_of_StateBeforeEvent[dIndex]);

                                //assign
                                G[Gindex] += this->X_scale_factors->operator[](Xindex)
                                    * TheDerivative
                                    / this->myUniverse->LU / this->myUniverse->LU * this->myUniverse->TU;
                            }
                            else //velocity
                            {
                                double drdotv_dCartesianState = this->state_before_event(stateIndex - 3) _GETVALUE;

                                //chain to the actual state representation
                                double TheDerivative = drdotv_dCartesianState * std::get<2>(this->Derivatives_of_StateBeforeEvent[dIndex]);

                                //assign
                                G[Gindex] += this->X_scale_factors->operator[](Xindex)
                                    * TheDerivative
                                    / this->myUniverse->LU / this->myUniverse->LU * this->myUniverse->TU;
                            }
                        }
                    }//end derivatives
                }//end generic case
            }//end rdotv constraint

            //Step 3: periapse distance constraint
            if (this->need_distance_constraint)
            {
                doubleType distance = this->state_before_event.getSubMatrix1D(0, 2).norm();

                F[Findex++] = distance / this->myUniverse->LU;

                if (needG)
                {
                    //first zero out the derivatives
                    for (size_t Gindex : this->Gindices_distance_constraint)
                    {
                        G[Gindex] = 0.0;
                    }

                    //now assign the entries
                    for (size_t entryIndex = 0; entryIndex < this->dIndices_distance_constraint.size(); ++entryIndex)
                    {
                        size_t dIndex = this->dIndices_distance_constraint[entryIndex];

                        std::tuple<size_t, size_t, double>& derivativeEntry = this->Derivatives_of_StateBeforeEvent[dIndex];

                        size_t stateIndex = std::get<1>(derivativeEntry);

                        size_t Xindex = std::get<0>(derivativeEntry);

                        size_t Gindex = this->Gindices_distance_constraint[entryIndex];

                        double TheDerivative = (this->state_before_event(stateIndex) / distance)_GETVALUE * std::get<2>(this->Derivatives_of_StateBeforeEvent[dIndex]);

                        G[Gindex] += this->X_scale_factors->operator[](Xindex)
                            * TheDerivative
                            / this->myUniverse->LU;
                    }//end loop over derivative indices
                }//end derivatives
            }//end distance constraint

            //Step 4: set the boundary state, which we need for printing
            this->boundary_state.shallow_copy(this->state_before_event);
        }//end process_event_left_side()

        
        void PeriapseBoundary::process_event_right_side(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //handled by the base class
            this->FreePointBoundary::process_event_right_side(X, Xindex, F, Findex, G, needG);
        }//end process_event_right_side

    }//end namespace BoundaryEvents
}//end namespace EMTG