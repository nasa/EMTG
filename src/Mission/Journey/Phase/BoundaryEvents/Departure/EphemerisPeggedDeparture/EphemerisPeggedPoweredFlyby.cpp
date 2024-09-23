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

#include "EphemerisPeggedPoweredFlyby.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        EphemerisPeggedPoweredFlyby::EphemerisPeggedPoweredFlyby(const std::string& name,
            const size_t& journeyIndex,
            const size_t& phaseIndex,
            size_t& stageIndex,
            Astrodynamics::universe* Universe,
            HardwareModels::Spacecraft* mySpacecraft,
            missionoptions* myOptions,
			EphemerisPeggedArrivalWithVinfinity* PreviousPhaseArrivalEvent) :
            EphemerisPeggedFlybyOut::EphemerisPeggedFlybyOut(name,
                journeyIndex,
                phaseIndex,
                stageIndex,
                Universe,
                mySpacecraft,
                myOptions,
                PreviousPhaseArrivalEvent)
        {
            this->hasBipropManeuver = true;
            this->hasManeuver = true;

            this->dDeltavMagnitude_dVinfinity_in.resize(3, 1, 0.0);
            this->dDeltavMagnitude_dVinfinity_out.resize(3, 1, 0.0);

            this->Gindices_TurnAngleConstraint_with_respect_to_Vinfinity.resize(3);
            this->Gindices_VirtualChemicalFuel_with_respect_to_Vinfinity.resize(3);
            this->Gindices_VirtualChemicalOxidizer_with_respect_to_Vinfinity.resize(3);
            this->Gindices_Deltav_with_respect_to_Vinfinity.resize(3);
        }

        //******************************************calcbounds methods

        //calcbounds
        void EphemerisPeggedPoweredFlyby::calcbounds(std::vector<size_t> timeVariables)
        {
            this->calcbounds_event_left_side(timeVariables);

            //create bounds for the v-infinity
			double vinf_max = fmax(this->myUniverse->LU / this->myUniverse->TU, this->PreviousPhaseArrivalEvent->get_Vinfinity_upperbound());

            std::vector< std::tuple<double, double> > vinfBounds = std::vector< std::tuple<double, double> >(3, std::make_tuple(-vinf_max, vinf_max));

            this->calcbounds_event_main(vinfBounds);

            this->calcbounds_event_right_side();

            this->calcbounds_virtual_propellant_constraints();

            if (this->myOptions->objective_type == ObjectiveFunctionType::MINIMIZE_DELTAV)
                this->calcbounds_deltav_contribution();

            this->calcbounds_specialized_constraints();
        }//end calcbounds()

        void EphemerisPeggedPoweredFlyby::
            calcbounds_event_main(const std::vector< std::tuple<double, double> >& vinfBounds)
        {
            //Step 1: base class
            this->EphemerisPeggedFlybyOut::calcbounds_event_main(vinfBounds);

            //Step 2: effect of v-infinity on mass
            for (size_t Vindex = 0; Vindex < 3; ++Vindex)
            {
                this->Derivatives_of_StateBeforeEvent.push_back( {this->Xindices_Vinfinity_in[Vindex], 6, 1.0} );
                this->dIndices_Mass_wrt_Vinfinity_in.push_back(this->Derivatives_of_StateBeforeEvent.size() - 1);

                this->Derivatives_of_StateBeforeEvent.push_back( {this->Xindices_Vinfinity_out[Vindex], 6, 1.0} );
                this->dIndices_Mass_wrt_Vinfinity_out.push_back(this->Derivatives_of_StateBeforeEvent.size() - 1);
            }

            //Step 3: encode periapse distance
            if (this->myJourneyOptions->override_flyby_altitude_bounds)
            {
                Xlowerbounds->push_back(this->myBody->radius + this->myJourneyOptions->flyby_altitude_bounds[0]);
                Xupperbounds->push_back(this->myBody->radius + this->myJourneyOptions->flyby_altitude_bounds[1]);
            }
            else
            {
                Xlowerbounds->push_back(this->myBody->radius + this->myBody->minimum_safe_flyby_altitude);
                if (this->myBody->mass < 1.0e+25)
                    Xupperbounds->push_back(10.0 * this->myBody->radius);
                else
                    Xupperbounds->push_back(300.0 * this->myBody->radius);
            }
            X_scale_factors->push_back(this->myBody->radius);
            Xdescriptions->push_back(prefix + "flyby periapse radius");
            //effect of flyby periapse distance on mass
            this->Derivatives_of_StateBeforeEvent.push_back(std::make_tuple<size_t, size_t, double>(Xdescriptions->size() - 1, 6, 1.0));
            this->dIndex_Mass_wrt_FlybyPeriapseDistance = this->Derivatives_of_StateBeforeEvent.size() - 1;


            //Step 4: turn angle feasibility constraint
            Flowerbounds->push_back(-math::SMALL);
            Fupperbounds->push_back(math::SMALL);
            Fdescriptions->push_back(prefix + "Powered flyby turn angle constraint");

            //derivatives with respect to periapse distance, v-infinity in, and v-infinity out
            //rp
            this->create_sparsity_entry(Fdescriptions->size() - 1,
                Xdescriptions->size() - 1,
                Gindices_TurnAngleConstraint_with_respect_to_FlybyPeriapseDistance);

            //V-infinity
            for (size_t Vindex = 0; Vindex < 3; ++Vindex)
            {
                this->create_sparsity_entry(Fdescriptions->size() - 1,
                    this->Xindices_Vinfinity_in[Vindex],
                    this->Gindices_TurnAngleConstraint_with_respect_to_Vinfinity[Vindex]);

                this->create_sparsity_entry(Fdescriptions->size() - 1,
                    this->Xindices_Vinfinity_out[Vindex],
                    this->Gindices_TurnAngleConstraint_with_respect_to_Vinfinity[Vindex]);
            }
        }//end calcbounds_event_main()

        void EphemerisPeggedPoweredFlyby::calcbounds_virtual_propellant_constraints()
        {
            //call the base class propellant tank calcbounds
            BoundaryEventBase::calcbounds_virtual_propellant_constraints();

            //specialty derivative entries for this phase type

            //derivatives of virtual chemical fuel with respect to initial mass, periapse distance, and v-infinity, in that order
            this->create_sparsity_entry(this->Findex_VirtualChemicalFuelConstraint,
                this->PreviousPhaseArrivalEvent->getX_index_of_first_decision_variable_in_this_event(),
                true,
                "event left state mass",
                Gindices_dVirtualChemicalFuel_dLeftMass);

            //periapse distance
            this->create_sparsity_entry(this->Findex_VirtualChemicalFuelConstraint,
                this->X_index_of_first_decision_variable_in_this_event,
                true,
                "flyby periapse radius",
                this->Gindices_VirtualChemicalFuel_with_respect_to_FlybyPeriapseDistance);

            //v-infinity
            for (size_t Vindex = 0; Vindex < 3; ++Vindex)
            {
                this->create_sparsity_entry(this->Findex_VirtualChemicalFuelConstraint,
                    this->Xindices_Vinfinity_in[Vindex],
                    this->Gindices_VirtualChemicalFuel_with_respect_to_Vinfinity[Vindex]);

                this->create_sparsity_entry(this->Findex_VirtualChemicalFuelConstraint,
                    this->Xindices_Vinfinity_out[Vindex],
                    this->Gindices_VirtualChemicalFuel_with_respect_to_Vinfinity[Vindex]);
            }

            //derivatives of virtual chemical oxidizer with respect to initial mass, periapse distance, and v-infinity, in that order
            this->create_sparsity_entry(this->Findex_VirtualChemicalOxidizerConstraint,
                this->PreviousPhaseArrivalEvent->getX_index_of_first_decision_variable_in_this_event(),
                true,
                "event left state mass",
                Gindices_dVirtualChemicalOxidizer_dLeftMass);

            //periapse distance
            this->create_sparsity_entry(this->Findex_VirtualChemicalOxidizerConstraint,
                this->X_index_of_first_decision_variable_in_this_event,
                true,
                "flyby periapse radius",
                this->Gindices_VirtualChemicalOxidizer_with_respect_to_FlybyPeriapseDistance);

            //v-infinity
            for (size_t Vindex = 0; Vindex < 3; ++Vindex)
            {
                this->create_sparsity_entry(this->Findex_VirtualChemicalOxidizerConstraint,
                    this->Xindices_Vinfinity_in[Vindex],
                    this->Gindices_VirtualChemicalOxidizer_with_respect_to_Vinfinity[Vindex]);

                this->create_sparsity_entry(this->Findex_VirtualChemicalOxidizerConstraint,
                    this->Xindices_Vinfinity_out[Vindex],
                    this->Gindices_VirtualChemicalOxidizer_with_respect_to_Vinfinity[Vindex]);
            }
        }//end calcbounds_virtual_propellant_constraints()

        void EphemerisPeggedPoweredFlyby::calcbounds_deltav_contribution()
        {
            //add derivatives with respect to periapse distance and v-infinity
            //periapse distance
            this->create_sparsity_entry(0,
                this->X_index_of_first_decision_variable_in_this_event,
                true,
                "flyby periapse radius",
                this->Gindices_Deltav_with_respect_to_FlybyPeriapseDistance);

            //v-infinity
            for (size_t Vindex = 0; Vindex < 3; ++Vindex)
            {
                this->create_sparsity_entry(0,
                    this->Xindices_Vinfinity_in[Vindex],
                    this->Gindices_Deltav_with_respect_to_Vinfinity[Vindex]);

                this->create_sparsity_entry(0,
                    this->Xindices_Vinfinity_out[Vindex],
                    this->Gindices_Deltav_with_respect_to_Vinfinity[Vindex]);
            }
        }//end calcbounds_deltav_contribution()

        //******************************************process methods

        void EphemerisPeggedPoweredFlyby::process_event(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            this->process_staging();

            this->process_event_left_side(X, Xindex, F, Findex, G, needG);

            this->process_event_main(X, Xindex, F, Findex, G, needG);

            this->process_event_right_side(X, Xindex, F, Findex, G, needG);

            this->process_virtual_propellant_constraints(X, Xindex, F, Findex, G, needG);

            this->process_deltav_contribution(X, Xindex, F, Findex, G, needG);

            this->process_specialized_constraints(X, Xindex, F, Findex, G, needG);
        }//end process_event

        void EphemerisPeggedPoweredFlyby::process_event_main(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //Step 1: base class
            this->EphemerisPeggedFlybyOut::process_event_main(X, Xindex, F, Findex, G, needG);

            //Step 2: extract flyby periapse distance
            this->FlybyPeriapseDistance = X[Xindex++];

            this->FlybyAltitude = this->FlybyPeriapseDistance - this->myBody->radius;

            //Step 3: compute turn angle constraint
            doubleType V0dotV0 = Vinfinity_in.dot(Vinfinity_in);
            doubleType VfdotVf = Vinfinity_out.dot(Vinfinity_out);
            doubleType vinf_in = sqrt(V0dotV0);
            doubleType vinf_out = sqrt(VfdotVf);
            double mu = this->myBody->mu;

            //compute turn angle
            doubleType denom = vinf_in * vinf_out;
            this->FlybyTurnAngle = math::safe_acos(this->Vinfinity_in.dot(this->Vinfinity_out) / ((denom > math::SMALL) ? denom : math::SMALL));

            //compute eccentricity of incoming and outgoing hyperbolas
            doubleType e_in = 1 + V0dotV0 * this->FlybyPeriapseDistance / mu;
            e_in = (e_in > math::SMALL) ? e_in : math::SMALL;

            doubleType e_out = 1 + VfdotVf * this->FlybyPeriapseDistance / mu;
            e_out = (e_out > math::SMALL) ? e_out : math::SMALL;


            //apply the flyby turn angle feasibility constraint
            F[Findex++] = math::safe_asin(1 / e_in) + math::safe_asin(1 / e_out) - this->FlybyTurnAngle;

            //Step 4: compute powered flyby delta-v
            doubleType A_dv = sqrt(V0dotV0 + 2.0 * mu / this->FlybyPeriapseDistance);
            doubleType B_dv = sqrt(VfdotVf + 2.0 * mu / this->FlybyPeriapseDistance);
            this->FlybyDeltavSigned = B_dv - A_dv;
            this->FlybyDeltavMagnitude = sqrt(FlybyDeltavSigned * FlybyDeltavSigned) + 1.0e-10;


            //Step 5: perform the powered flyby burn
            this->mySpacecraft->computeChemicalPropulsionPerformance(this->FlybyDeltavMagnitude, this->state_before_event(6), true, this->ChemicalManeuverType); //biprop

            this->chemical_fuel_used += this->mySpacecraft->getChemFuelConsumedThisManeuver();
            this->chemical_oxidizer_used += this->mySpacecraft->getChemOxidizerConsumedThisManeuver();

            this->dChemicalFuel_dDeltavMagnitude = this->mySpacecraft->getdFuelConsumedThisManeuver_ddeltav();
            this->dChemicalOxidizer_dDeltavMagnitude = this->mySpacecraft->getdOxidizerConsumedThisManeuver_ddeltav();

            this->dChemicalFuel_dInitialMass = (this->chemical_fuel_used / this->state_before_event(6)) _GETVALUE;
            this->dChemicalOxidizer_dInitialMass = (this->chemical_oxidizer_used / this->state_before_event(6)) _GETVALUE;

            this->state_before_event(6) -= (this->chemical_fuel_used + this->chemical_oxidizer_used);

            //derivatives of powered flyby turn angle constraint and delta-v
            //(thanks SymPy!)
            if (needG)
            {
                size_t Gindex, Xindex;

                double v0x = this->Vinfinity_in(0) _GETVALUE;
                double v0y = this->Vinfinity_in(1) _GETVALUE;
                double v0z = this->Vinfinity_in(2) _GETVALUE;
                double vfx = this->Vinfinity_out(0) _GETVALUE;
                double vfy = this->Vinfinity_out(1) _GETVALUE;
                double vfz = this->Vinfinity_out(2) _GETVALUE;
                double V0dotVf = Vinfinity_in.dot(Vinfinity_out) _GETVALUE;
                double e_in2 = (e_in * e_in) _GETVALUE;
                double e_out2 = (e_out * e_out) _GETVALUE;
                double mf = this->state_before_event(6) _GETVALUE;
                double rp = this->FlybyPeriapseDistance _GETVALUE;

                //with respect to rp
                Gindex = Gindices_TurnAngleConstraint_with_respect_to_FlybyPeriapseDistance;
                Xindex = this->jGvar->operator[](Gindex);
                G[Gindex] = this->X_scale_factors->operator[](Xindex)
                    * (-(VfdotVf) / (mu*e_out2*sqrt(1 - 1 / e_out2)) - (V0dotV0) / (mu*e_in2*sqrt(1 - 1 / e_in2))) _GETVALUE;

                //with respect to v_infinity_in and v_infinity_out
                double dF_dvfx = ( (-v0x*(V0dotVf) / (sqrt(V0dotV0*VfdotVf)*V0dotV0) + vfx / sqrt(V0dotV0*VfdotVf)) / sqrt(1 - V0dotVf*V0dotVf / (V0dotV0*VfdotVf)) - 2 * rp*v0x / (mu*e_in2*sqrt(1 - 1 / e_in2)) )_GETVALUE;
                double dF_dvfy = ( (-v0y*(V0dotVf) / (sqrt(V0dotV0*VfdotVf)*V0dotV0) + vfy / sqrt(V0dotV0*VfdotVf)) / sqrt(1 - V0dotVf*V0dotVf / (V0dotV0*VfdotVf)) - 2 * rp*v0y / (mu*e_in2*sqrt(1 - 1 / e_in2)) )_GETVALUE;
                double dF_dvfz = ( (-v0z*(V0dotVf) / (sqrt(V0dotV0*VfdotVf)*V0dotV0) + vfz / sqrt(V0dotV0*VfdotVf)) / sqrt(1 - V0dotVf*V0dotVf / (V0dotV0*VfdotVf)) - 2 * rp*v0z / (mu*e_in2*sqrt(1 - 1 / e_in2)) )_GETVALUE;


                double dF_dv0x = ( (v0x / sqrt(V0dotV0*VfdotVf) - vfx*(V0dotVf) / (sqrt(V0dotV0*VfdotVf)*VfdotVf)) / sqrt(1 - V0dotVf*V0dotVf / (V0dotV0*VfdotVf)) - 2 * rp*vfx / (mu*e_out2*sqrt(1 - 1 / e_out2)) )_GETVALUE;
                double dF_dv0y = ( (v0y / sqrt(V0dotV0*VfdotVf) - vfy*(V0dotVf) / (sqrt(V0dotV0*VfdotVf)*VfdotVf)) / sqrt(1 - V0dotVf*V0dotVf / (V0dotV0*VfdotVf)) - 2 * rp*vfy / (mu*e_out2*sqrt(1 - 1 / e_out2)) )_GETVALUE;
                double dF_dv0z = ( (v0z / sqrt(V0dotV0*VfdotVf) - vfz*(V0dotVf) / (sqrt(V0dotV0*VfdotVf)*VfdotVf)) / sqrt(1 - V0dotVf*V0dotVf / (V0dotV0*VfdotVf)) - 2 * rp*vfz / (mu*e_out2*sqrt(1 - 1 / e_out2)) )_GETVALUE;
                
                Gindex = Gindices_TurnAngleConstraint_with_respect_to_Vinfinity[0][1];
                Xindex = this->jGvar->operator[](Gindex);
                G[Gindex] = this->X_scale_factors->operator[](Xindex) * dF_dv0x;

                Gindex = Gindices_TurnAngleConstraint_with_respect_to_Vinfinity[1][1];
                Xindex = this->jGvar->operator[](Gindex);
                G[Gindex] = this->X_scale_factors->operator[](Xindex) * dF_dv0y;

                Gindex = Gindices_TurnAngleConstraint_with_respect_to_Vinfinity[2][1];
                Xindex = this->jGvar->operator[](Gindex);
                G[Gindex] = this->X_scale_factors->operator[](Xindex) * dF_dv0z;

                Gindex = Gindices_TurnAngleConstraint_with_respect_to_Vinfinity[0][0];
                Xindex = this->jGvar->operator[](Gindex);
                G[Gindex] = this->X_scale_factors->operator[](Xindex) * dF_dvfx;

                Gindex = Gindices_TurnAngleConstraint_with_respect_to_Vinfinity[1][0];
                Xindex = this->jGvar->operator[](Gindex);
                G[Gindex] = this->X_scale_factors->operator[](Xindex) * dF_dvfy;

                Gindex = Gindices_TurnAngleConstraint_with_respect_to_Vinfinity[2][0];
                Xindex = this->jGvar->operator[](Gindex);
                G[Gindex] = this->X_scale_factors->operator[](Xindex) * dF_dvfz;

                //for anti-velocity burns, change the signs of all the derivatives of deltav
                double sign = this->FlybyDeltavSigned > 0 ? -1.0 : 1.0;

                this->dDeltavMagnitude_dFlybyPeriapseDistance = sign * mu / (rp * rp) * (1 / B_dv - 1 / A_dv) _GETVALUE;
                this->dDeltavMagnitude_dVinfinity_in(0) = sign * v0x / A_dv _GETVALUE;
                this->dDeltavMagnitude_dVinfinity_in(1) = sign * v0y / A_dv _GETVALUE;
                this->dDeltavMagnitude_dVinfinity_in(2) = sign * v0z / A_dv _GETVALUE;
                this->dDeltavMagnitude_dVinfinity_out(0) = sign * -vfx / B_dv _GETVALUE;
                this->dDeltavMagnitude_dVinfinity_out(1) = sign * -vfy / B_dv _GETVALUE;
                this->dDeltavMagnitude_dVinfinity_out(2) = sign * -vfz / B_dv _GETVALUE;

                double dSpacecraftMass_dDeltav = -(this->dChemicalFuel_dDeltavMagnitude + dChemicalOxidizer_dDeltavMagnitude);

                //derivatives of mass with respect to periapse altitude
                std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_Mass_wrt_FlybyPeriapseDistance]) = -dSpacecraftMass_dDeltav
                    * this->dDeltavMagnitude_dFlybyPeriapseDistance;

                //derivatives of mass with respect to v-infinity
                for (size_t Vindex = 0; Vindex < 3; ++Vindex)
                {
                    //v-infinity in
                    std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndices_Mass_wrt_Vinfinity_in[Vindex]]) =
                        -dSpacecraftMass_dDeltav * this->dDeltavMagnitude_dVinfinity_in(Vindex);

                    //v-infinity out
                    std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndices_Mass_wrt_Vinfinity_out[Vindex]]) =
                        -dSpacecraftMass_dDeltav * this->dDeltavMagnitude_dVinfinity_out(Vindex);
                }

                //derivatives of mass with respect to mass
                std::get<2>(this->Derivatives_of_StateBeforeEvent[this->dIndex_CurrentEventMass_PreviousEventMass]) =
                    (this->state_before_event(6) / (this->state_before_event(6) + this->chemical_fuel_used + this->chemical_oxidizer_used))_GETVALUE;
            }
        }//end process_event_main()

        void EphemerisPeggedPoweredFlyby::process_virtual_propellant_constraints(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //Step 1: base class
            BoundaryEventBase::process_virtual_propellant_constraints(X, Xindex, F, Findex, G, needG);

            //Step 2: specialized class
            {
                size_t Gindex, Xindex;

                //Step 2.1: derivatives with respect to left encoded mass
                Gindex = this->Gindices_dVirtualChemicalFuel_dLeftMass;
                Xindex = this->jGvar->operator[](Gindex);
                G[Gindex] = this->X_scale_factors->operator[](Xindex)
                    * -this->dChemicalFuel_dInitialMass
                    * this->myUniverse->continuity_constraint_scale_factors(6);

                Gindex = this->Gindices_dVirtualChemicalOxidizer_dLeftMass;
                Xindex = this->jGvar->operator[](Gindex);
                G[Gindex] = this->X_scale_factors->operator[](Xindex)
                    * -this->dChemicalOxidizer_dInitialMass
                    * this->myUniverse->continuity_constraint_scale_factors(6);

                //Step 2.2: derivatives with respect to periapse distance
                Gindex = this->Gindices_VirtualChemicalFuel_with_respect_to_FlybyPeriapseDistance;
                Xindex = this->jGvar->operator[](Gindex);
                G[Gindex] = this->X_scale_factors->operator[](Xindex)
                    * this->dChemicalFuel_dDeltavMagnitude * this->dDeltavMagnitude_dFlybyPeriapseDistance
                    * this->myUniverse->continuity_constraint_scale_factors(6);

                Gindex = this->Gindices_VirtualChemicalOxidizer_with_respect_to_FlybyPeriapseDistance;
                Xindex = this->jGvar->operator[](Gindex);
                G[Gindex] = this->X_scale_factors->operator[](Xindex)
                    * this->dChemicalOxidizer_dDeltavMagnitude * this->dDeltavMagnitude_dFlybyPeriapseDistance
                    * this->myUniverse->continuity_constraint_scale_factors(6);

                //Step 2.3: derivatives with respect to v-infinity
                for (size_t Vindex = 0; Vindex < 3; ++Vindex)
                {
                    //v-infinity in
                    Gindex = this->Gindices_VirtualChemicalFuel_with_respect_to_Vinfinity[Vindex][0];
                    Xindex = this->jGvar->operator[](Gindex);
                    G[Gindex] = this->X_scale_factors->operator[](Xindex)
                        * this->dChemicalFuel_dDeltavMagnitude * this->dDeltavMagnitude_dVinfinity_in(Vindex)
                        * this->myUniverse->continuity_constraint_scale_factors(6);

                    Gindex = this->Gindices_VirtualChemicalOxidizer_with_respect_to_Vinfinity[Vindex][0];
                    Xindex = this->jGvar->operator[](Gindex);
                    G[Gindex] = this->X_scale_factors->operator[](Xindex)
                        * this->dChemicalOxidizer_dDeltavMagnitude * this->dDeltavMagnitude_dVinfinity_in(Vindex)
                        * this->myUniverse->continuity_constraint_scale_factors(6);

                    //v-infinity out
                    Gindex = this->Gindices_VirtualChemicalFuel_with_respect_to_Vinfinity[Vindex][1];
                    Xindex = this->jGvar->operator[](Gindex);
                    G[Gindex] = this->X_scale_factors->operator[](Xindex)
                        * this->dChemicalFuel_dDeltavMagnitude * this->dDeltavMagnitude_dVinfinity_out(Vindex)
                        * this->myUniverse->continuity_constraint_scale_factors(6);

                    Gindex = this->Gindices_VirtualChemicalOxidizer_with_respect_to_Vinfinity[Vindex][1];
                    Xindex = this->jGvar->operator[](Gindex);
                    G[Gindex] = this->X_scale_factors->operator[](Xindex)
                        * this->dChemicalOxidizer_dDeltavMagnitude * this->dDeltavMagnitude_dVinfinity_out(Vindex)
                        * this->myUniverse->continuity_constraint_scale_factors(6);
                }
            }
        }//end process_virtual_propellant_constraints()

        void EphemerisPeggedPoweredFlyby::process_deltav_contribution(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //Step 1: assign the event delta-v
            this->EventDeterministicDeltav = this->FlybyDeltavMagnitude;


            //Step 2: derivatives
            if (this->myOptions->objective_type == ObjectiveFunctionType::MINIMIZE_DELTAV
                && needG)
            {
                size_t Gindex, Xindex;

                //Step 3.1: derivatives with respect to periapse distance
                Gindex = Gindices_Deltav_with_respect_to_FlybyPeriapseDistance;
                Xindex = this->jGvar->operator[](Gindex);
                G[Gindex] = this->X_scale_factors->operator[](Xindex)
                    * this->dDeltavMagnitude_dFlybyPeriapseDistance;

                //Step 3.3: derivatives with respect to v-infinity
                for (size_t Vindex = 0; Vindex < 3; ++Vindex)
                {
                    //v-infinity in
                    Gindex = Gindices_Deltav_with_respect_to_Vinfinity[Vindex][0];
                    Xindex = this->jGvar->operator[](Gindex);
                    G[Gindex] = this->X_scale_factors->operator[](Xindex)
                        * this->dDeltavMagnitude_dVinfinity_in(Vindex);

                    //v-infinity out
                    Gindex = Gindices_Deltav_with_respect_to_Vinfinity[Vindex][1];
                    Xindex = this->jGvar->operator[](Gindex);
                    G[Gindex] = this->X_scale_factors->operator[](Xindex)
                        * this->dDeltavMagnitude_dVinfinity_out(Vindex);
                }
            }
        }//end process_deltav_contribution()

         //******************************************output methods
        void EphemerisPeggedPoweredFlyby::output(std::ofstream& outputfile,
            const double& launchdate,
            size_t& eventcount)
        {
            this->mySpacecraft->setActiveStage(this->stageIndex);

            std::string event_type;
            
            if (this->FlybyDeltavMagnitude > 1.0e-4)
                event_type = "pwr_flyby";
            else
                event_type = "upwr_flyby";

            std::string boundary_name = this->myBody->name;

            math::Matrix<doubleType> empty3vector(3, 1, 0.0);

            //where is the Sun?
            math::Matrix<doubleType> R_sc_Sun(3, 1, 0.0);
            if (this->myUniverse->central_body_SPICE_ID == 10)            {
                R_sc_Sun = this->state_after_event.getSubMatrix1D(0, 2);
            }
            else
            {
                //where is the central body relative to the sun?
                doubleType central_body_state_and_derivatives[12];
                this->myUniverse->locate_central_body(this->state_after_event(7),
                    central_body_state_and_derivatives,
                    *this->myOptions,
                    false);
                math::Matrix<doubleType> R_CB_Sun(3, 1, 0.0);
                for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
                {
                    R_CB_Sun(stateIndex) = central_body_state_and_derivatives[stateIndex];
                }
                R_sc_Sun = this->state_after_event.getSubMatrix1D(0, 2) + R_CB_Sun;
            }
            this->mySpacecraft->computePowerState(R_sc_Sun.getSubMatrix1D(0, 2).norm() / this->myOptions->AU, this->state_after_event(7));

            //rotate the v-infinity vector to the local frame
            this->myUniverse->LocalFrame.construct_rotation_matrices(this->EventLeftEpoch / 86400.0 + 2400000.5, false);

            math::Matrix<doubleType> R_ICRF_to_local = this->myUniverse->LocalFrame.get_R(EMTG::ReferenceFrame::ICRF,
                this->myJourneyOptions->arrival_elements_frame);

            math::Matrix<doubleType> VinfinityLocalFrame = R_ICRF_to_local * this->Vinfinity_out;

            //compute RA and DEC
            doubleType RA = atan2(VinfinityLocalFrame(1), VinfinityLocalFrame(0));
            doubleType DEC = asin(VinfinityLocalFrame(2) / VinfinityLocalFrame.norm());

            //compute BdotR and BdotT
            math::Matrix<doubleType> periapse_state = this->calculate_flyby_periapse_state();
            this->myBplane.define_bplane(this->get_periapse_state());
            this->BdotR = this->myBplane.getBdotR();
            this->BdotT = this->myBplane.getBdotT();

            double Isp = this->ChemicalManeuverType == PropulsionSystemChoice::Biprop ? this->mySpacecraft->getBipropIsp() : this->mySpacecraft->getMonopropIsp();

            write_output_line(outputfile,
                eventcount,
                event_type,
                boundary_name,
                this->EventTimeWidth,
                this->FlybyAltitude,
                this->BdotR,
                this->BdotT,
                RA,
                DEC,
                this->Vinfinity_out.dot(this->Vinfinity_out),
                this->state_after_event,
                this->Vinfinity_out,
                empty3vector,
                this->FlybyDeltavMagnitude,
                0.0,
                Isp,
                this->mySpacecraft->getAvailablePower(),
                0.0,
                0,
                0.0,
                "none");
        }

        math::Matrix<doubleType> EphemerisPeggedPoweredFlyby::get_periapse_state()
        {
            //calculate unit vector pointed towards periapse
            math::Matrix<doubleType> periapse_position_unit_vector = (this->Vinfinity_in.unitize() - this->Vinfinity_out.unitize()).unitize();

            //calculate angular momentum unit vector
            math::Matrix<doubleType> angular_momentum_vector = this->Vinfinity_in.cross(this->Vinfinity_out);

            //calculate velocity unit vector at periapse
            math::Matrix<doubleType> periapse_velocity_unit_vector = (angular_momentum_vector.cross(periapse_position_unit_vector)).unitize();

            //calculate velocity magnitude at periapse
            doubleType periapse_velocity_magnitude = sqrt(2 * this->myBody->mu / (this->myBody->radius + this->FlybyAltitude) + this->Vinfinity_in.dot(this->Vinfinity_in));

            //transform from unit vector space to state space
            math::Matrix<doubleType> periapse_position_vector = periapse_position_unit_vector * (this->myBody->radius + this->FlybyAltitude);
            math::Matrix<doubleType> periapse_velocity_vector = periapse_velocity_unit_vector * periapse_velocity_magnitude;

            math::Matrix<doubleType> periapse_state = periapse_position_vector.vert_cat(periapse_velocity_vector);

            return periapse_state;
        }//end get_periapse_state()

        void EphemerisPeggedPoweredFlyby::output_maneuver_and_target_spec(std::ofstream& maneuver_spec_file, std::ofstream& target_spec_file, bool& haveManeuverNeedTarget)
        {
            //in this case, we don't actually know the periapse state, so we'll use the approximated one that we already output in the .emtg file
            math::Matrix<doubleType> periapse_state = this->get_periapse_state();

            if (haveManeuverNeedTarget)
            {
                //reset the flag that tells us we need a target spec line
                haveManeuverNeedTarget = false;

                //Step 1: target spec

                //Step 1.1: initialize a target spec object
                this->myBplane.define_bplane(periapse_state);

                target_spec_line myTargetSpecLine(this->name,
                    "EME2000",
                    this->myUniverse->central_body.name,
                    this->state_after_event(7),
                    this->state_after_event,
                    this->myBplane.getBdotR(),
                    this->myBplane.getBdotT());

                //Step 1.2: write target spec object
                myTargetSpecLine.write(target_spec_file);
            }

            //Step 2: maneuver spec
            //Step 2.1: instantiate and populate a maneuver spec object
            maneuver_spec_line myManeuverSpecLine(this->name);

            double Isp = this->ChemicalManeuverType == PropulsionSystemChoice::Biprop ? this->mySpacecraft->getBipropIsp() : this->mySpacecraft->getMonopropIsp();
            doubleType massFlowRate = this->mySpacecraft->getchemthrust() / Isp / this->myOptions->g0;
            doubleType mass_before_maneuver = this->state_before_event(6) + this->chemical_fuel_used + this->chemical_oxidizer_used;
            doubleType maneuverDuration = (mass_before_maneuver - this->state_before_event(6)) / massFlowRate;
            
            //maneuver direction is along track or anti-track, depending on the flyby delta-v 
            math::Matrix<doubleType> periapse_maneuver_vector = periapse_state.getSubMatrix1D(3, 5).unitize() * this->FlybyDeltavSigned;

            myManeuverSpecLine.append_maneuver_spec_item("EME2000",
                this->state_before_event(7),
                periapse_maneuver_vector.unitize(),
                mass_before_maneuver,
                this->state_before_event(6),
                this->mySpacecraft->getchemthrust(),
                massFlowRate,
                maneuverDuration,
                1.0);

            //Step 2.2: write the maneuver spec
            myManeuverSpecLine.write(maneuver_spec_file);

            //Step 2.3: signal that we need a target spec
            haveManeuverNeedTarget = true;
        }//end output_maneuver_and_target_spec()
    }//end namespace BoundaryEvents
}//end namespace EMTG