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

#include "EphemerisPeggedChemRendezvous.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        EphemerisPeggedChemRendezvous::EphemerisPeggedChemRendezvous(const std::string& name,
            const size_t& journeyIndex,
            const size_t& phaseIndex,
            size_t& stageIndex,
            Astrodynamics::universe* Universe,
            HardwareModels::Spacecraft* mySpacecraft,
            missionoptions* myOptions) :
            EphemerisPeggedArrivalWithVinfinity::EphemerisPeggedArrivalWithVinfinity(name,
                journeyIndex,
                phaseIndex,
                stageIndex,
                Universe,
                mySpacecraft,
                myOptions)
        {
            this->hasBipropManeuver = true;
            this->hasManeuver = true;

            if (this->hasTCM)
                this->hasMonopropManeuver = true;
        }

        //******************************************calcbounds methods

        //calcbounds
        void EphemerisPeggedChemRendezvous::calcbounds(std::vector<size_t> timeVariables)
        {
            //create bounds for the v-infinity
            double vinf_max = this->myUniverse->LU / this->myUniverse->TU;

            std::vector< std::tuple<double, double> > vinfBounds = std::vector< std::tuple<double, double> >(3, std::make_tuple(-vinf_max, vinf_max));

            this->calcbounds_event_left_side(timeVariables);

            this->calcbounds_event_main(vinfBounds);

            this->calcbounds_event_right_side();

            this->calcbounds_virtual_propellant_constraints();

            if (this->myOptions->objective_type == ObjectiveFunctionType::MINIMIZE_DELTAV)
                this->calcbounds_deltav_contribution();

            this->calcbounds_specialized_constraints();
        }//end calcbounds()

        void EphemerisPeggedChemRendezvous::
            calcbounds_event_right_side()
        {
            //Step 1: base class
            this->EphemerisPeggedArrivalWithVinfinity::calcbounds_event_right_side();

            //Step 2: derivative of mass with respect to v-infinity
            for (size_t Vindex = 0; Vindex < 3; ++Vindex)
            {
                size_t Xindex = std::get<0>(this->Derivatives_of_StateBeforeEvent[this->dIndex_VbeforeEvent_dVinfinity_in[Vindex]]);
                this->Derivatives_of_StateAfterEvent.push_back({ Xindex, 6, 1.0 });
                this->dIndex_mass_after_event_wrt_v_infinity.push_back(this->Derivatives_of_StateAfterEvent.size() - 1);
            }
        }

        void EphemerisPeggedChemRendezvous::calcbounds_virtual_propellant_constraints()
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

        void EphemerisPeggedChemRendezvous::calcbounds_deltav_contribution()
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
        void EphemerisPeggedChemRendezvous::process_event(const std::vector<doubleType>& X,
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

            BoundaryEventBase::process_specialized_constraints(X, Xindex, F, Findex, G, needG);
        }//end process_event

        void EphemerisPeggedChemRendezvous::process_event_main(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //Step 1: call the base class
            this->EphemerisPeggedArrivalWithVinfinity::process_event_main(X, Xindex, F, Findex, G, needG);

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

            this->mySpacecraft->computeChemicalPropulsionPerformance(this->ArrivalDeltavMagnitude, this->state_before_event(6) - this->chemical_fuel_used, true, this->ChemicalManeuverType); //could be monoprop or biprop

            this->chemical_fuel_used += this->mySpacecraft->getChemFuelConsumedThisManeuver();
            this->chemical_oxidizer_used += this->mySpacecraft->getChemOxidizerConsumedThisManeuver();

            this->dChemicalFuel_ddeltavMagnitude = this->mySpacecraft->getdFuelConsumedThisManeuver_ddeltav();
            this->dChemicalOxidizer_ddeltavMagnitude = this->mySpacecraft->getdOxidizerConsumedThisManeuver_ddeltav();

           
            //derivatives
            if (needG)
            {
                //derivative of mass with respect to mass
                double Isp = this->ChemicalManeuverType == PropulsionSystemChoice::Biprop ? this->mySpacecraft->getBipropIsp() : this->mySpacecraft->getMonopropIsp();
                this->ETM(6, 6) = exp(-1000.0 * this->TCM_magnitude / (this->mySpacecraft->getMonopropIsp() * this->myOptions->g0) - 1000.0 * this->ArrivalDeltavMagnitude / (Isp * this->myOptions->g0)) _GETVALUE;
                

                //derivative of mass with respect to v-infinity
                for (size_t Vindex = 0; Vindex < 3; ++Vindex)
                {
                    this->ETM(6, 3 + Vindex) = (this->Vinfinity_in(Vindex) / this->ArrivalDeltavMagnitude * this->mySpacecraft->getdSpacecraftMassThisManeuver_ddeltav()) _GETVALUE;
                }
            }
        }//end process_event_main()

        void EphemerisPeggedChemRendezvous::process_event_right_side(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //Step 1: base class
            this->EphemerisPeggedArrival::process_event_right_side(X, Xindex, F, Findex, G, needG);
            
            //Step 2: specialized entries
            std::get<2>(this->Derivatives_of_StateAfterEvent[this->dIndex_mass_wrt_encodedMass]) = this->ETM(6, 6);


            for (size_t Vindex = 0; Vindex < 3; ++Vindex)
            {
                std::get<2>(this->Derivatives_of_StateAfterEvent[this->dIndex_mass_after_event_wrt_v_infinity[Vindex]]) = this->ETM(6, 3 + Vindex)
                    * (this->state_after_event(6) / (this->state_after_event(6) + this->journey_end_propellant_used))_GETVALUE;
            }
        }//end process_event_right_side()

        void EphemerisPeggedChemRendezvous::process_virtual_propellant_constraints(const std::vector<doubleType>& X,
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
                    Xindex = this->jGvar->operator[](Gindex);
                    double Gentry = (this->chemical_fuel_used / this->state_before_event(6)) _GETVALUE;
                    G[Gindex] = this->X_scale_factors->operator[](Xindex) * -Gentry
                        * this->myUniverse->continuity_constraint_scale_factors(6);

                    //wrt v-infinity
                    Gindex = Gindices_dVirtualChemicalFuel_dVinfinity[Vindex];
                    Xindex = this->jGvar->operator[](Gindex);

                    double ddvMag_dVinfinity = (this->Vinfinity_in(Vindex) / this->ArrivalDeltavMagnitude) _GETVALUE;
                    
                    G[Gindex] = this->X_scale_factors->operator[](Xindex)
                        * this->dChemicalFuel_ddeltavMagnitude * ddvMag_dVinfinity
                        * this->myUniverse->continuity_constraint_scale_factors(6);

                    //oxidizer
                    //wrt mass
                    Gindex = Gindices_dVirtualChemicalOxidizer_dLeftMass;
                    Xindex = this->jGvar->operator[](Gindex);
                    Gentry = (this->chemical_oxidizer_used / this->state_before_event(6)) _GETVALUE;
                    G[Gindex] = this->X_scale_factors->operator[](Xindex) * -Gentry
                        * this->myUniverse->continuity_constraint_scale_factors(6);

                    //wrt v-infinity
                    Gindex = Gindices_dVirtualChemicalOxidizer_dVinfinity[Vindex];
                    Xindex = this->jGvar->operator[](Gindex);

                    G[Gindex] = this->X_scale_factors->operator[](Xindex)
                        * this->dChemicalOxidizer_ddeltavMagnitude * ddvMag_dVinfinity
                        * this->myUniverse->continuity_constraint_scale_factors(6);
                }
            }
        }//end process_virtual_propellant_constraints()

        void EphemerisPeggedChemRendezvous::process_deltav_contribution(const std::vector<doubleType>& X,
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
                    Xindex = this->jGvar->operator[](Gindex);
                    G[Gindex] = this->X_scale_factors->operator[](Xindex) * ddvMag_dVinfinity;
                }
            }
        }//end process_virtual_deltav_constraint()

         //******************************************output methods
        void EphemerisPeggedChemRendezvous::output(std::ofstream& outputfile,
            const double& launchdate,
            size_t& eventcount)
        {
            this->mySpacecraft->setActiveStage(this->stageIndex);

            std::string event_type = "rendezvous";

            std::string boundary_name = this->myBody->name;

            math::Matrix<doubleType> empty3vector(3, 1, 0.0);

            //where is the Sun?
            math::Matrix<doubleType> R_sc_Sun(3, 1, 0.0);
            if (this->myUniverse->central_body_SPICE_ID == 10)
            {
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

            //compute RA and DEC for incoming asymptote in the body's local frame
            doubleType RA, DEC;
            //rotate the v-infinity vector to the local frame
            this->myUniverse->LocalFrame.construct_rotation_matrices(this->EventLeftEpoch / 86400.0 + 2400000.5, false);

            math::Matrix<doubleType> R_ICRF_to_local = this->myUniverse->LocalFrame.get_R(EMTG::ReferenceFrame::ICRF,
                this->myJourneyOptions->arrival_elements_frame);

            math::Matrix<doubleType> VinfinityLocalFrame = R_ICRF_to_local * this->Vinfinity_in;

            //compute RA and DEC
            RA = atan2(VinfinityLocalFrame(1), VinfinityLocalFrame(0));
            DEC = asin(VinfinityLocalFrame(2) / VinfinityLocalFrame.norm());

            double Isp = this->ChemicalManeuverType == PropulsionSystemChoice::Biprop ? this->mySpacecraft->getBipropIsp() : this->mySpacecraft->getMonopropIsp();

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
                this->Vinfinity_in.dot(this->Vinfinity_in),
                this->state_after_event,
                -this->Vinfinity_in,
                empty3vector,
                this->Vinfinity_in.norm(),
                0.0,
                Isp,
                this->mySpacecraft->getAvailablePower(),
                0.0,
                0,
                0.0,
                "none");

            //output end of mission if appropriate
            this->EphemerisPeggedArrival::output(outputfile, launchdate, eventcount);
        }//end output()

        void EphemerisPeggedChemRendezvous::output_maneuver_and_target_spec(std::ofstream& maneuver_spec_file, std::ofstream& target_spec_file, bool& haveManeuverNeedTarget)
        {
            //Here there IS an arrival maneuver, so we will do the following steps:
            //1. write out a target spec for whatever maneuvers previously occurred, i.e. the state before the maneuver
            //2. write out a maneuver spec for the arrival maneuver
            //3. write out a target spec for the arrival maneuver, i.e. state_after_event

            math::Matrix<doubleType> output_state = this->state_before_event;

            if (haveManeuverNeedTarget)
            {
                //Step 1: target for whatever happened before this boundary event
                //reset the flag that tells us we need a target spec line
                haveManeuverNeedTarget = false;
                //Step 1.1: initialize a target spec object
                for (size_t vIndex : {0, 1, 2})
                    output_state(3 + vIndex) -= this->Vinfinity_in(vIndex);

                target_spec_line myTargetSpecLine(this->name,
                    "EME2000",
                    this->myUniverse->central_body.name,
                    output_state(7),
                    output_state);

                //Step 1.2: write target spec object
                myTargetSpecLine.write(target_spec_file);
            }

            //Step 2: maneuver spec for the arrival maneuver
            //Step 2.1: initialize a maneuver spec object
            double Isp = this->ChemicalManeuverType == PropulsionSystemChoice::Biprop ? this->mySpacecraft->getBipropIsp() : this->mySpacecraft->getMonopropIsp();
            doubleType massFlowRate = this->mySpacecraft->getchemthrust() / Isp / this->myOptions->g0;
            double maneuverDuration = 0.0;//this is a bit of a hack - it prevents the target from being AFTER the maneuver. We know that in the real world we have to move the maneuver backward in time.
            maneuver_spec_line myManeuverSpecLine(this->name + "_arrival_maneuver");
            myManeuverSpecLine.append_maneuver_spec_item("EME2000",
                this->state_after_event(7),
                -this->Vinfinity_in.unitize(),
                output_state(6),
                this->state_after_event(6),
                this->mySpacecraft->getchemthrust(),
                massFlowRate,
                maneuverDuration,
                1.0);

            //Step 2.2: write maneuver spec object
            myManeuverSpecLine.write(maneuver_spec_file);

            //Step 3: target spec for arrival maneuver
            //Step 3.1: initialize a target spec object
            target_spec_line myTargetSpecLine(this->name + "_arrival_maneuver",
                "EME2000",
                this->myUniverse->central_body.name,
                this->state_after_event(7),
                this->state_after_event);

            //Step 3.2: write target spec object
            myTargetSpecLine.write(target_spec_file);
        }//end output_maneuver_and_target_spec()
    }//end namespace BoundaryEvents
}//end namespace EMTG