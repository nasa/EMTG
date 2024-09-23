// EMTG: Evolutionary Mission Trajectory Generator
// An open-source global optimization tool for preliminary mission design
// Provided by NASA Goddard Space Flight Center
//
// Copyright (c) 2014 - 2017 United States Government as represented by the
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

#include "FreePointChemRendezvous.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        FreePointChemRendezvous::FreePointChemRendezvous(const std::string name,
            const size_t& journeyIndex,
            const size_t& phaseIndex,
            size_t& stageIndex,
            Astrodynamics::universe* Universe,
            HardwareModels::Spacecraft* mySpacecraft,
            missionoptions* myOptions) :
            FreePointArrivalWithVinfinity::FreePointArrivalWithVinfinity(name,
                journeyIndex,
                phaseIndex,
                stageIndex,
                Universe,
                mySpacecraft,
                myOptions)
        {
            this->hasBipropManeuver = true;

            if (this->hasTCM)
                this->hasMonopropManeuver = true;
        }

        //******************************************calcbounds methods

        //calcbounds
        void FreePointChemRendezvous::calcbounds()
        {
            //create bounds for the v-infinity
            double vinf_max = this->myUniverse->LU / this->myUniverse->TU;

            std::vector< std::tuple<double, double> > vinfBounds = std::vector< std::tuple<double, double> >(3, std::make_tuple(-vinf_max, vinf_max));

            this->calcbounds_event_left_side();

            this->calcbounds_event_main(vinfBounds);

            this->calcbounds_event_right_side();

            this->calcbounds_virtual_propellant_constraints();

            if (this->myOptions->objective_type == ObjectiveFunctionType::MINIMIZE_DELTAV)
                this->calcbounds_deltav_contribution();

            this->calcbounds_specialized_constraints();
        }//end calcbounds()

        void FreePointChemRendezvous::
            calcbounds_event_right_side()
        {
            //Step 1: base class
            this->FreePointArrivalWithVinfinity::calcbounds_event_right_side();

            //Step 2: derivative of mass with respect to v-infinity
            for (size_t Vindex = 0; Vindex < 3; ++Vindex)
            {
                size_t Xindex = std::get<0>(this->Derivatives_of_StateBeforeEvent[this->dIndex_VbeforeEvent_dVinfinity_in[Vindex]]);
                this->Derivatives_of_StateAfterEvent.push_back({ Xindex, 6, 1.0 });
                this->dIndex_mass_after_event_wrt_v_infinity.push_back(this->Derivatives_of_StateAfterEvent.size() - 1);
            }
        }

        void FreePointChemRendezvous::calcbounds_virtual_propellant_constraints()
        {
            //call the base class propellant tank calcbounds
            BoundaryEventBase::calcbounds_virtual_propellant_constraints();

            //specialty derivative entries for this phase type

            //derivatives of virtual chemical fuel with respect to initial mass and v-infinity, in that order
            this->create_sparsity_entry(this->Findex_VirtualChemicalFuelConstraint,
                this->X_index_of_first_decision_variable_in_this_event,
                true,
                "mass",
                this->Gindices_dVirtualChemicalFuel_dLeftMass);

            //v-infinity
            for (size_t Vindex = 0; Vindex < 3; ++Vindex)
            {
                this->create_sparsity_entry(this->Findex_VirtualChemicalFuelConstraint,
                    std::get<0>(this->Derivatives_of_StateBeforeEvent[dIndex_VbeforeEvent_dVinfinity_in[Vindex]]),
                    this->Gindices_dVirtualChemicalFuel_dVinfinity);
            }

            //derivatives of virtual chemical oxidizer with respect to initial mass and v-infinity, in that order
            this->create_sparsity_entry(this->Findex_VirtualChemicalOxidizerConstraint,
                this->X_index_of_first_decision_variable_in_this_event,
                true,
                "mass",
                this->Gindices_dVirtualChemicalOxidizer_dLeftMass);

            //v-infinity
            for (size_t Vindex = 0; Vindex < 3; ++Vindex)
            {
                this->create_sparsity_entry(this->Findex_VirtualChemicalOxidizerConstraint,
                    std::get<0>(this->Derivatives_of_StateBeforeEvent[dIndex_VbeforeEvent_dVinfinity_in[Vindex]]),
                    this->Gindices_dVirtualChemicalOxidizer_dVinfinity);
            }
        }//end calcbounds_virtual_propellant_constraints()

        void FreePointChemRendezvous::calcbounds_deltav_contribution()
        {
            //add derivatives with respect to incoming v-infinity
            for (size_t Vindex = 0; Vindex < 3; ++Vindex)
            {
                this->create_sparsity_entry(0,
                    std::get<0>(this->Derivatives_of_StateBeforeEvent[dIndex_VbeforeEvent_dVinfinity_in[Vindex]]),
                    this->Gindices_dDeltav_dVinfinity);
            }
        }//end calcbounds_virtual_deltav_constraint()

        //******************************************process methods
        void FreePointChemRendezvous::process_event(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            this->process_event_left_side(X, Xindex, F, Findex, G, needG);

            this->process_event_main(X, Xindex, F, Findex, G, needG);

            this->process_event_right_side(X, Xindex, F, Findex, G, needG);

            this->process_virtual_propellant_constraints(X, Xindex, F, Findex, G, needG);

            this->process_deltav_contribution(X, Xindex, F, Findex, G, needG);

            BoundaryEventBase::process_specialized_constraints(X, Xindex, F, Findex, G, needG);
        }//end process_event

        void FreePointChemRendezvous::process_event_main(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //Step 1: call the base class
            this->FreePointArrivalWithVinfinity::process_event_main(X, Xindex, F, Findex, G, needG);

            //Step 2: reset the propellant
            this->chemical_fuel_used = 0.0;
            this->chemical_oxidizer_used = 0.0;

            //Step 3: perform a TCM if applicable
            if (this->hasTCM)
            {
                //perform a TCM with the current stage monoprop system
                this->mySpacecraft->computeChemicalPropulsionPerformance(this->TCM_magnitude, this->state_before_event(6), true, PropulsionSystemChoice::Monoprop); //monoprop

                this->chemical_fuel_used += this->mySpacecraft->getChemFuelConsumedThisManeuver();
            }

            //Step 3: perform the arrival maneuver
            this->ArrivalDeltavMagnitude = this->Vinfinity_in.norm();

            this->mySpacecraft->computeChemicalPropulsionPerformance(this->ArrivalDeltavMagnitude, this->state_before_event(6) - this->chemical_fuel_used, true, this->ChemicalManeuverType); //biprop

            this->chemical_fuel_used += this->mySpacecraft->getChemFuelConsumedThisManeuver();
            this->chemical_oxidizer_used += this->mySpacecraft->getChemOxidizerConsumedThisManeuver();

            this->dChemicalFuel_ddeltavMagnitude = this->mySpacecraft->getdFuelConsumedThisManeuver_ddeltav();
            this->dChemicalOxidizer_ddeltavMagnitude = this->mySpacecraft->getdOxidizerConsumedThisManeuver_ddeltav();

           
            //derivatives
            if (needG)
            {
                //derivative of mass with respect to mass
                this->ETM(6, 6) = exp(-1000.0 * this->TCM_magnitude / (this->mySpacecraft->getMonopropIsp() * this->myOptions->g0) - 1000.0 * this->ArrivalDeltavMagnitude / (this->mySpacecraft->getBipropIsp() * this->myOptions->g0)) _GETVALUE;
                

                //derivative of mass with respect to v-infinity
                for (size_t Vindex = 0; Vindex < 3; ++Vindex)
                {
                    this->ETM(6, 3 + Vindex) = (this->Vinfinity_in(Vindex) / this->ArrivalDeltavMagnitude * this->mySpacecraft->getdSpacecraftMassThisManeuver_ddeltav()) _GETVALUE;
                }
            }
        }//end process_event_main()

        void FreePointChemRendezvous::process_event_right_side(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //Step 1: base class
            this->FreePointArrival::process_event_right_side(X, Xindex, F, Findex, G, needG);
            
            //Step 2: specialized entries
            std::get<2>(this->Derivatives_of_StateAfterEvent[this->dIndex_mass_wrt_encodedMass]) = this->ETM(6, 6);


            for (size_t Vindex = 0; Vindex < 3; ++Vindex)
            {
                std::get<2>(this->Derivatives_of_StateAfterEvent[this->dIndex_mass_after_event_wrt_v_infinity[Vindex]]) = this->ETM(6, 3 + Vindex);
            }
        }//end process_event_right_side()

        void FreePointChemRendezvous::process_virtual_propellant_constraints(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //call the base class propellant constraints
            BoundaryEventBase::process_virtual_propellant_constraints(X, Xindex, F, Findex, G, needG);

            //event-specific
            if (needG)
            {
                size_t Xindex, Gindex;
                for (size_t Vindex = 0; Vindex < 3; ++Vindex)
                {
                    //fuel
                    //wrt mass
                    Gindex = Gindices_dVirtualChemicalFuel_dLeftMass;
                    Xindex = this->myOptions->jGvar[Gindex];
                    double Gentry = (this->chemical_fuel_used / this->state_before_event(6)) _GETVALUE;
                    G[Gindex] = this->myOptions->X_scale_factors[Xindex] * -Gentry
                        * this->myUniverse->continuity_constraint_scale_factors(6);

                    //wrt v-infinity
                    Gindex = Gindices_dVirtualChemicalFuel_dVinfinity[Vindex];
                    Xindex = this->myOptions->jGvar[Gindex];

                    double ddvMag_dVinfinity = (this->Vinfinity_in(Vindex) / this->ArrivalDeltavMagnitude) _GETVALUE;
                    
                    G[Gindex] = this->myOptions->X_scale_factors[Xindex]
                        * this->dChemicalFuel_ddeltavMagnitude * -ddvMag_dVinfinity
                        * this->myUniverse->continuity_constraint_scale_factors(6);

                    //oxidizer
                    //wrt mass
                    Gindex = Gindices_dVirtualChemicalOxidizer_dLeftMass;
                    Xindex = this->myOptions->jGvar[Gindex];
                    Gentry = (this->chemical_oxidizer_used / this->state_before_event(6)) _GETVALUE;
                    G[Gindex] = this->myOptions->X_scale_factors[Xindex] * -Gentry
                        * this->myUniverse->continuity_constraint_scale_factors(6);

                    //wrt v-infinity
                    Gindex = Gindices_dVirtualChemicalOxidizer_dVinfinity[Vindex];
                    Xindex = this->myOptions->jGvar[Gindex];

                    G[Gindex] = this->myOptions->X_scale_factors[Xindex]
                        * this->dChemicalOxidizer_ddeltavMagnitude * -ddvMag_dVinfinity
                        * this->myUniverse->continuity_constraint_scale_factors(6);
                }
            }
        }//end process_virtual_propellant_constraints()

        void FreePointChemRendezvous::process_deltav_contribution(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //add the event delta-v
            this->EventDeterministicDeltav = this->ArrivalDeltavMagnitude;

            //event-specific
            if (this->myOptions->objective_type == ObjectiveFunctionType::MINIMIZE_DELTAV
                && needG)
            {
                size_t Xindex, Gindex;
                for (size_t Vindex = 0; Vindex < 3; ++Vindex)
                {
                    //derivatives with respect to v-infinity
                    double ddvMag_dVinfinity = (this->Vinfinity_in(Vindex) / this->ArrivalDeltavMagnitude) _GETVALUE;

                    Gindex = this->Gindices_dDeltav_dVinfinity[Vindex];
                    Xindex = this->myOptions->jGvar[Gindex];
                    G[Gindex] = this->myOptions->X_scale_factors[Xindex] * ddvMag_dVinfinity;
                }
            }
        }//end process_virtual_deltav_constraint()

         //******************************************output methods
        void FreePointChemRendezvous::output(std::ofstream& outputfile,
            const double& launchdate,
            size_t& eventcount)
        {
            std::string event_type = "rendezvous";

            std::string boundary_name = "free point";

            math::Matrix<doubleType> empty3vector(3, 1, 0.0);

            this->mySpacecraft->computePowerState(this->state_after_event.getSubMatrix1D(0, 3).norm() / this->myUniverse->LU,
                (this->EventRightEpoch - launchdate), this->myOptions->power_margin);

            //compute RA and DEC for incoming asymptote in the body's local frame
            doubleType RA, DEC;
            RA = 0.0;
            DEC = 0.0;

            write_output_line(outputfile,
                eventcount,
                event_type,
                boundary_name,
                this->EventTimeWidth,
                0.0,
                0.0,
                0.0,
                RA,
                DEC,
                this->Vinfinity_in.norm(),
                this->state_after_event,
                -this->Vinfinity_in,
                empty3vector,
                this->Vinfinity_in.norm(),
                0.0,
                this->mySpacecraft->getBipropIsp(),
                this->mySpacecraft->getProducedPower(),
                0.0,
                0,
                0.0);
        }
    }//end namespace BoundaryEvents
}//end namespace EMTG