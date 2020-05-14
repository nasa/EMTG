
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
            BoundaryEventBase::BoundaryEventBase()
        {}

        PeriapseBoundary::PeriapseBoundary(const std::string& name,
            const size_t& journeyIndex,
            const size_t& phaseIndex,
            size_t& stageIndex,
            Astrodynamics::universe* Universe,
            HardwareModels::Spacecraft* mySpacecraft,
            missionoptions* myOptions)

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
            //base class
            this->BoundaryEventBase::initialize(name,
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
                                                          const std::vector<double>& MassBounds)
        {

            //Step 1: set the current stage
            this->mySpacecraft->setActiveStage(this->stageIndex);

            if (this->myOptions->PeriapseBoundaryStateRepresentation == StateRepresentation::SphericalRADEC)
            {
                //Step 2: position variables
                //Step 2.1: radius
                Xlowerbounds->push_back(RadiusBounds[0]);
                Xupperbounds->push_back(RadiusBounds[1]);
                X_scale_factors->push_back(this->myUniverse->central_body.radius);
                Xdescriptions->push_back(prefix + "event left state r");
                this->Xindex_rMag = this->Xdescriptions->size() - 1;
                for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
                {
                    this->Derivatives_of_StateBeforeEvent.push_back(std::make_tuple(this->Xdescriptions->size() - 1, stateIndex, 1.0));
                    this->dIndex_periapse_state_wrt_r.push_back(this->Derivatives_of_StateBeforeEvent.size() - 1);
                }

                //Step 2.2: RA
                this->Xlowerbounds->push_back(-8.0 * math::PI);
                this->Xupperbounds->push_back(8.0 * math::PI);
                this->X_scale_factors->push_back(1.0);
                this->Xdescriptions->push_back(this->prefix + "event left state RA");
                this->Xindex_RA = this->Xdescriptions->size() - 1;
                //affects x, y
                for (size_t stateIndex = 0; stateIndex < 2; ++stateIndex)
                {
                    this->Derivatives_of_StateBeforeEvent.push_back(std::make_tuple(this->Xdescriptions->size() - 1, stateIndex, 1.0));
                    this->dIndex_periapse_state_wrt_RA.push_back(this->Derivatives_of_StateBeforeEvent.size() - 1);
                }

                //Step 2.3: DEC
                this->Xlowerbounds->push_back(-0.5 * math::PI);
                this->Xupperbounds->push_back(0.5 * math::PI);
                this->X_scale_factors->push_back(1.0);
                this->Xdescriptions->push_back(this->prefix + "event left state DEC");
                this->Xindex_DEC = this->Xdescriptions->size() - 1;
                for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
                {
                    this->Derivatives_of_StateBeforeEvent.push_back(std::make_tuple(this->Xdescriptions->size() - 1, stateIndex, 1.0));
                    this->dIndex_periapse_state_wrt_DEC.push_back(this->Derivatives_of_StateBeforeEvent.size() - 1);
                }

                //Step 3: velocity variables
                //Step 3.1: v
                this->Xlowerbounds->push_back(VelocityMagnitudeBounds[0]);
                this->Xupperbounds->push_back(VelocityMagnitudeBounds[1]);
                this->X_scale_factors->push_back(this->myUniverse->LU / this->myUniverse->TU);
                this->Xdescriptions->push_back(this->prefix + "event left state v");
                this->Xindex_vMag = this->Xdescriptions->size() - 1;
                for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
                {
                    this->Derivatives_of_StateBeforeEvent.push_back(std::make_tuple(this->Xdescriptions->size() - 1, stateIndex + 3, 1.0));
                    this->dIndex_periapse_state_wrt_v.push_back(this->Derivatives_of_StateBeforeEvent.size() - 1);
                }

                //Step 3.2: velocity RA
                this->Xlowerbounds->push_back(-8.0 * math::PI);
                this->Xupperbounds->push_back(8.0 * math::PI);
                this->X_scale_factors->push_back(1.0);
                this->Xdescriptions->push_back(this->prefix + "event left state vRA");
                this->Xindex_vRA = this->Xdescriptions->size() - 1;
                for (size_t stateIndex = 0; stateIndex < 2; ++stateIndex)
                {
                    this->Derivatives_of_StateBeforeEvent.push_back(std::make_tuple(this->Xdescriptions->size() - 1, stateIndex + 3, 1.0));
                    this->dIndex_periapse_state_wrt_vRA.push_back(this->Derivatives_of_StateBeforeEvent.size() - 1);
                }

                //Step 3.3: velocity DEC
                this->Xlowerbounds->push_back(-math::PI / 2.0 - math::SMALL);
                this->Xupperbounds->push_back(math::PI / 2.0 + math::SMALL);
                this->X_scale_factors->push_back(1.0);
                this->Xdescriptions->push_back(this->prefix + "event left state vDEC");
                this->Xindex_vDEC = this->Xdescriptions->size() - 1;
                for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
                {
                    this->Derivatives_of_StateBeforeEvent.push_back(std::make_tuple(this->Xdescriptions->size() - 1, stateIndex + 3, 1.0));
                    this->dIndex_periapse_state_wrt_vDEC.push_back(this->Derivatives_of_StateBeforeEvent.size() - 1);
                }
            }//end SphericalRADEC
            else if (this->myOptions->PeriapseBoundaryStateRepresentation == StateRepresentation::SphericalAZFPA)
            {
                //Step 2: position variables
                //Step 2.1: radius
                Xlowerbounds->push_back(RadiusBounds[0]);
                Xupperbounds->push_back(RadiusBounds[1]);
                X_scale_factors->push_back(this->myUniverse->central_body.radius);
                Xdescriptions->push_back(prefix + "event left state r");
                this->Xindex_rMag = this->Xdescriptions->size() - 1;
                for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
                {
                    this->Derivatives_of_StateBeforeEvent.push_back(std::make_tuple(this->Xdescriptions->size() - 1, stateIndex, 1.0));
                    this->dIndex_periapse_state_wrt_r.push_back(this->Derivatives_of_StateBeforeEvent.size() - 1);
                }

                //Step 2.2: RA
                this->Xlowerbounds->push_back(-8.0 * math::PI);
                this->Xupperbounds->push_back(8.0 * math::PI);
                this->X_scale_factors->push_back(1.0);
                this->Xdescriptions->push_back(this->prefix + "event left state RA");
                this->Xindex_RA = this->Xdescriptions->size() - 1;
                //affects x, y, xdot, ydot but not z, zdot 
                for (size_t stateIndex = 0; stateIndex < 2; ++stateIndex)
                {
                    this->Derivatives_of_StateBeforeEvent.push_back(std::make_tuple(this->Xdescriptions->size() - 1, stateIndex, 1.0));
                    this->dIndex_periapse_state_wrt_RA.push_back(this->Derivatives_of_StateBeforeEvent.size() - 1);
                }
                for (size_t stateIndex = 3; stateIndex < 5; ++stateIndex)
                {
                    this->Derivatives_of_StateBeforeEvent.push_back(std::make_tuple(this->Xdescriptions->size() - 1, stateIndex, 1.0));
                    this->dIndex_periapse_state_wrt_RA.push_back(this->Derivatives_of_StateBeforeEvent.size() - 1);
                }

                //Step 2.3: DEC
                this->Xlowerbounds->push_back(-0.5 * math::PI);
                this->Xupperbounds->push_back(0.5 * math::PI);
                this->X_scale_factors->push_back(1.0);
                this->Xdescriptions->push_back(this->prefix + "event left state DEC");
                this->Xindex_DEC = this->Xdescriptions->size() - 1;
                for (size_t stateIndex = 0; stateIndex < 6; ++stateIndex)
                {
                    this->Derivatives_of_StateBeforeEvent.push_back(std::make_tuple(this->Xdescriptions->size() - 1, stateIndex, 1.0));
                    this->dIndex_periapse_state_wrt_DEC.push_back(this->Derivatives_of_StateBeforeEvent.size() - 1);
                }

                //Step 3: velocity variables
                //Step 3.1: v
                this->Xlowerbounds->push_back(VelocityMagnitudeBounds[0]);
                this->Xupperbounds->push_back(VelocityMagnitudeBounds[1]);
                this->X_scale_factors->push_back(this->myUniverse->LU / this->myUniverse->TU);
                this->Xdescriptions->push_back(this->prefix + "event left state v");
                this->Xindex_vMag = this->Xdescriptions->size() - 1;
                for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
                {
                    this->Derivatives_of_StateBeforeEvent.push_back(std::make_tuple(this->Xdescriptions->size() - 1, stateIndex + 3, 1.0));
                    this->dIndex_periapse_state_wrt_v.push_back(this->Derivatives_of_StateBeforeEvent.size() - 1);
                }

                //Step 3.2: velocity AZ
                this->Xlowerbounds->push_back(-8.0 * math::PI);
                this->Xupperbounds->push_back(8.0 * math::PI);
                this->X_scale_factors->push_back(1.0);
                this->Xdescriptions->push_back(this->prefix + "event left state AZ");
                this->Xindex_AZ = this->Xdescriptions->size() - 1;
                for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
                {
                    this->Derivatives_of_StateBeforeEvent.push_back(std::make_tuple(this->Xdescriptions->size() - 1, stateIndex + 3, 1.0));
                    this->dIndex_periapse_state_wrt_AZ.push_back(this->Derivatives_of_StateBeforeEvent.size() - 1);
                }

                //Step 3.3: velocity FPA
                this->Xlowerbounds->push_back(math::PI / 2.0 - math::SMALL);
                this->Xupperbounds->push_back(math::PI / 2.0 + math::SMALL);
                this->X_scale_factors->push_back(1.0);
                this->Xdescriptions->push_back(this->prefix + "event left state FPA");
                this->Xindex_FPA = this->Xdescriptions->size() - 1;
                for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
                {
                    this->Derivatives_of_StateBeforeEvent.push_back(std::make_tuple(this->Xdescriptions->size() - 1, stateIndex + 3, 1.0));
                    this->dIndex_periapse_state_wrt_FPA.push_back(this->Derivatives_of_StateBeforeEvent.size() - 1);
                }
            }//end SphericalAZFPA

            //Step 4: mass variable
            this->Xlowerbounds->push_back(MassBounds[0]);
            this->Xupperbounds->push_back(MassBounds[1]);
            this->X_scale_factors->push_back(1.0 / this->myUniverse->continuity_constraint_scale_factors(6));
            this->Xdescriptions->push_back(this->prefix + "event left state mass");
            this->Xindex_mass = this->Xdescriptions->size() - 1;
            this->Derivatives_of_StateBeforeEvent.push_back(std::make_tuple(this->Xindex_mass, 6, 1.0));
            this->dIndex_mass_wrt_encodedMass = this->Derivatives_of_StateBeforeEvent.size() - 1;

            //Step 5: r dot v = 0 constraint
            if (this->myOptions->PeriapseBoundaryStateRepresentation == StateRepresentation::SphericalRADEC)
            {
                this->Flowerbounds->push_back(-1.0e-13);
                this->Fupperbounds->push_back(1.0e-13);
                this->Fdescriptions->push_back(this->prefix + "periapse is an apse");

                this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                    this->Xindex_RA,
                    this->Gindex_rdotv_wrt_RA);
                this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                    this->Xindex_DEC,
                    this->Gindex_rdotv_wrt_DEC);
                this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                    this->Xindex_vRA,
                    this->Gindex_rdotv_wrt_vRA);
                this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                    this->Xindex_vDEC,
                    this->Gindex_rdotv_wrt_vDEC);
            }//end RdotV constraint

            this->calculate_dependencies_left_epoch();
        }//end calcbounds_left_side()

        void PeriapseBoundary::calcbounds_event_right_side()
        {
            this->Derivatives_of_StateAfterEvent = this->Derivatives_of_StateBeforeEvent;
            this->Derivatives_of_StateAfterEvent_wrt_Time = this->Derivatives_of_StateBeforeEvent_wrt_Time;
        }//end calcbounds_right_side()

        //******************************************process methods

        void PeriapseBoundary::process_event_left_side(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {

            if (this->myOptions->PeriapseBoundaryStateRepresentation == StateRepresentation::SphericalRADEC)
            {
                //Step 1: extract the thingies
                doubleType r = X[Xindex++];
                doubleType RA = X[Xindex++];
                doubleType DEC = X[Xindex++];
                doubleType v = X[Xindex++];
                doubleType vRA = X[Xindex++];
                doubleType vDEC = X[Xindex++];
                this->state_before_event(6) = X[Xindex++];

                //Step 2: convert to cartesian
                doubleType cosRA = cos(RA);
                doubleType sinRA = sin(RA);
                doubleType cosDEC = cos(DEC);
                doubleType sinDEC = sin(DEC);
                doubleType cosvRA = cos(vRA);
                doubleType sinvRA = sin(vRA);
                doubleType cosvDEC = cos(vDEC);
                doubleType sinvDEC = sin(vDEC);

                this->state_before_event(0) = r * cosRA * cosDEC;
                this->state_before_event(1) = r * sinRA * cosDEC;
                this->state_before_event(2) = r * sinDEC;
                this->state_before_event(3) = v * cosvRA * cosvDEC;
                this->state_before_event(4) = v * sinvRA * cosvDEC;
                this->state_before_event(5) = v * sinvDEC;

                //Step 3: apply r dot v = 0 constraint
                doubleType rdotv = sinDEC * sinvDEC + cosDEC * cosRA*cosvDEC*cosvRA + cosDEC * sinRA*cosvDEC*sinvRA;
                F[Findex++] = rdotv;

                //Step 4: derivatives
                if (needG)
                {
                    double dx_dr = (cosRA * cosDEC)_GETVALUE;
                    double dy_dr = (sinRA * cosDEC)_GETVALUE;
                    double dz_dr = sinDEC _GETVALUE;

                    double dx_dRA = (-r * cosDEC * sinRA) _GETVALUE;
                    double dy_dRA = (r * cosDEC * cosRA) _GETVALUE;

                    double dx_dDEC = (-r * cosRA*sinDEC)_GETVALUE;
                    double dy_dDEC = (-r*sinDEC*sinRA)_GETVALUE;
                    double dz_dDEC = (r*cosDEC)_GETVALUE;

                    double dxdot_dv = (cosvRA * cosvDEC)_GETVALUE;
                    double dydot_dv = (sinvRA * cosvDEC)_GETVALUE;
                    double dzdot_dv = sinvDEC _GETVALUE;

                    double dxdot_dvRA = (-v * cosvDEC * sinvRA) _GETVALUE;
                    double dydot_dvRA = (v * cosvDEC * cosvRA) _GETVALUE;

                    double dxdot_dvDEC = (-v * cosvRA*sinvDEC)_GETVALUE;
                    double dydot_dvDEC = (-v * sinvDEC*sinvRA)_GETVALUE;
                    double dzdot_dvDEC = (v*cosvDEC)_GETVALUE;

                    std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_r[0]]) = dx_dr;
                    std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_r[1]]) = dy_dr;
                    std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_r[2]]) = dz_dr;

                    std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_RA[0]]) = dx_dRA;
                    std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_RA[1]]) = dy_dRA;

                    std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_DEC[0]]) = dx_dDEC;
                    std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_DEC[1]]) = dy_dDEC;
                    std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_DEC[2]]) = dz_dDEC;

                    std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_v[0]]) = dxdot_dv;
                    std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_v[1]]) = dydot_dv;
                    std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_v[2]]) = dzdot_dv;

                    std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_vRA[0]]) = dxdot_dvRA;
                    std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_vRA[1]]) = dydot_dvRA;

                    std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_vDEC[0]]) = dxdot_dvDEC;
                    std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_vDEC[1]]) = dydot_dvDEC;
                    std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_vDEC[2]]) = dzdot_dvDEC;

                    //r dot v constraint
                    double drdotv_dRA = (cosDEC*cosRA*cosvDEC*sinvRA - cosDEC*sinRA*cosvDEC*cosvRA) _GETVALUE;
                    double drdotv_dDEC = -(cosRA*sinDEC*cosvDEC*cosvRA - cosDEC*sinvDEC + sinDEC*sinRA*cosvDEC*sinvRA)_GETVALUE;
                    double drdotv_dvRA = -(cosDEC*cosRA*cosvDEC*sinvRA - cosDEC*sinRA*cosvDEC*cosvRA) _GETVALUE;
                    double drdotv_dvDEC = -(cosDEC*cosRA*cosvRA*sinvDEC - sinDEC*cosvDEC + cosDEC*sinRA*sinvDEC*sinvRA) _GETVALUE;

                    G[this->Gindex_rdotv_wrt_RA] = this->X_scale_factors->operator[](this->Xindex_RA) * drdotv_dRA;
                    G[this->Gindex_rdotv_wrt_DEC] = this->X_scale_factors->operator[](this->Xindex_DEC) * drdotv_dDEC;
                    G[this->Gindex_rdotv_wrt_vRA] = this->X_scale_factors->operator[](this->Xindex_vRA) * drdotv_dvRA;
                    G[this->Gindex_rdotv_wrt_vDEC] = this->X_scale_factors->operator[](this->Xindex_vDEC) * drdotv_dvDEC;
                }
            }//end SphericalRADEC
            else if (this->myOptions->PeriapseBoundaryStateRepresentation == StateRepresentation::SphericalAZFPA)
            {
                //Step 1: extract the thingies
                doubleType r = X[Xindex++];
                doubleType RA = X[Xindex++];
                doubleType DEC = X[Xindex++];
                doubleType v = X[Xindex++];
                doubleType AZ = X[Xindex++];
                doubleType FPA = X[Xindex++];
                this->state_before_event(6) = X[Xindex++];

                //Step 2: convert to cartesian
                doubleType cosRA = cos(RA);
                doubleType sinRA = sin(RA);
                doubleType cosDEC = cos(DEC);
                doubleType sinDEC = sin(DEC);
                doubleType cosAZ = cos(AZ);
                doubleType sinAZ = sin(AZ);
                doubleType cosFPA = cos(FPA);
                doubleType sinFPA = sin(FPA);

                this->state_before_event(0) = r * cosRA * cosDEC;
                this->state_before_event(1) = r * sinRA * cosDEC;
                this->state_before_event(2) = r * sinDEC;
                this->state_before_event(3) = -v * (sinFPA*(sinAZ*sinRA + cosAZ * cosRA*sinDEC) - cosFPA * cosDEC*cosRA);
                this->state_before_event(4) = v * (sinFPA*(cosRA*sinAZ - cosAZ * sinDEC*sinRA) + cosFPA * cosDEC*sinRA);
                this->state_before_event(5) = v * (cosFPA*sinDEC + cosDEC * cosAZ*sinFPA);

                //Step 3: derivatives
                if (needG)
                {
                    double dx_dr = (cosRA * cosDEC)_GETVALUE;
                    double dy_dr = (sinRA * cosDEC)_GETVALUE;
                    double dz_dr = sinDEC _GETVALUE;

                    double dx_dRA = (-r * cosDEC * sinRA) _GETVALUE;
                    double dy_dRA = (r * cosDEC * cosRA) _GETVALUE;

                    double dx_dDEC = (-r * cosRA*sinDEC)_GETVALUE;
                    double dy_dDEC = (-r * sinDEC*sinRA)_GETVALUE;
                    double dz_dDEC = (r*cosDEC)_GETVALUE;

                    double dxdot_dRA = (-v * (sinFPA*(cosRA*sinAZ - cosAZ * sinDEC*sinRA) + cosFPA * cosDEC*sinRA)) _GETVALUE;
                    double dydot_dRA = (-v * (sinFPA*(sinAZ*sinRA + cosAZ * cosRA*sinDEC) - cosFPA * cosDEC*cosRA)) _GETVALUE;

                    double dxdot_dDEC = (-v * cosRA*(cosFPA*sinDEC + cosDEC * cosAZ*sinFPA)) _GETVALUE;
                    double dydot_dDEC = (-v * sinRA*(cosFPA*sinDEC + cosDEC * cosAZ*sinFPA)) _GETVALUE;
                    double dzdot_dDEC = (v*(cosFPA*cosDEC - cosAZ * sinFPA*sinDEC)) _GETVALUE;

                    double dxdot_dv = (-(sinFPA*(sinAZ*sinRA + cosAZ * cosRA*sinDEC) - cosFPA * cosDEC*cosRA))_GETVALUE;
                    double dydot_dv = ((sinFPA*(cosRA*sinAZ - cosAZ * sinDEC*sinRA) + cosFPA * cosDEC*sinRA))_GETVALUE;
                    double dzdot_dv = ((cosFPA*sinDEC + cosDEC * cosAZ*sinFPA))_GETVALUE;

                    double dxdot_dFPA = (-v * (cosFPA*(sinAZ*sinRA + cosAZ * cosRA*sinDEC) + cosDEC * sinFPA*cosRA)) _GETVALUE;
                    double dydot_dFPA = (v*(cosFPA*(cosRA*sinAZ - cosAZ * sinDEC*sinRA) - cosDEC * sinFPA*sinRA)) _GETVALUE;
                    double dzdot_dFPA = (v*cosFPA*cosDEC*cosAZ - v * sinFPA*sinDEC) _GETVALUE;

                    double dxdot_dAZ = (-v * sinFPA*(cosAZ*sinRA - cosRA * sinDEC*sinAZ)) _GETVALUE;
                    double dydot_dAZ = (v*sinFPA*(cosAZ*cosRA + sinDEC * sinAZ*sinRA)) _GETVALUE;
                    double dzdot_dAZ = (-v * cosDEC*sinFPA*sinAZ) _GETVALUE;

                    std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_r[0]]) = dx_dr;
                    std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_r[1]]) = dy_dr;
                    std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_r[2]]) = dz_dr;

                    std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_RA[0]]) = dx_dRA;
                    std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_RA[1]]) = dy_dRA;
                    std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_RA[2]]) = dxdot_dRA;
                    std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_RA[3]]) = dydot_dRA;

                    std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_DEC[0]]) = dx_dDEC;
                    std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_DEC[1]]) = dy_dDEC;
                    std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_DEC[2]]) = dz_dDEC;
                    std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_DEC[3]]) = dxdot_dDEC;
                    std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_DEC[4]]) = dydot_dDEC;
                    std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_DEC[5]]) = dzdot_dDEC;

                    std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_v[0]]) = dxdot_dv;
                    std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_v[1]]) = dydot_dv;
                    std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_v[2]]) = dzdot_dv;

                    std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_AZ[0]]) = dxdot_dAZ;
                    std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_AZ[1]]) = dydot_dAZ;
                    std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_AZ[2]]) = dzdot_dAZ;

                    std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_FPA[0]]) = dxdot_dFPA;
                    std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_FPA[1]]) = dydot_dFPA;
                    std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_periapse_state_wrt_FPA[2]]) = dzdot_dFPA;
                }
            }//end SphericalAZFPA


            //set the boundary state, which we need for printing
            this->boundary_state.shallow_copy(this->state_before_event);
        }//end process_event_left_side()

        
        void PeriapseBoundary::process_event_right_side(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            this->state_after_event.shallow_copy(this->state_before_event);
            
            for (size_t dIndex = 0; dIndex < this->Derivatives_of_StateBeforeEvent.size(); ++dIndex)
                this->Derivatives_of_StateAfterEvent[dIndex] = this->Derivatives_of_StateBeforeEvent[dIndex];

            for (size_t dIndex = 0; dIndex < this->Derivatives_of_StateBeforeEvent_wrt_Time.size(); ++dIndex)
                this->Derivatives_of_StateAfterEvent_wrt_Time[dIndex] = this->Derivatives_of_StateBeforeEvent_wrt_Time[dIndex];

            this->EventRightEpoch = this->EventLeftEpoch;
        }//end process_event_right_side

    }//end namespace BoundaryEvents
}//end namespace EMTG