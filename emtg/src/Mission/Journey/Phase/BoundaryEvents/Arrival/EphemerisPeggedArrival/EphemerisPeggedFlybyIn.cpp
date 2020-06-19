
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

#include "EphemerisPeggedFlybyIn.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        EphemerisPeggedFlybyIn::EphemerisPeggedFlybyIn(const std::string& name,
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
        }

        //******************************************calcbounds methods

        //calcbounds
        void EphemerisPeggedFlybyIn::calcbounds(std::vector<size_t> timeVariables)
        {
            //create bounds for the v-infinity
            double vinf_max = this->myUniverse->LU / this->myUniverse->TU;

            std::vector< std::tuple<double, double> > vinfBounds = std::vector< std::tuple<double, double> >(3, std::make_tuple(-vinf_max, vinf_max));

            this->calcbounds_event_left_side(timeVariables);

            this->EphemerisPeggedArrivalWithVinfinity::calcbounds_event_main(vinfBounds);

            this->calcbounds_event_right_side();

            if (this->hasTCM)
                this->calcbounds_virtual_propellant_constraints();

            this->calcbounds_specialized_constraints();
        }//end calcbounds()
        
        //calcbounds for propellant constraints - TCM
        void EphemerisPeggedFlybyIn::calcbounds_virtual_propellant_constraints()
        {
            if (this->hasTCM)
                this->hasMonopropManeuver = true;

            BoundaryEventBase::calcbounds_virtual_propellant_constraints();

            if (this->hasTCM)
            {
                //derivative of virtual chemical fuel with respect to initial mass
                this->create_sparsity_entry(this->Findex_VirtualChemicalFuelConstraint,
                    this->X_index_of_first_decision_variable_in_this_event,
                    true,
                    "mass",
                    Gindices_dVirtualChemicalFuel_TCM_dLeftMass);
            }
        }//end calcbounds_virtual_propellant_constraints

        //******************************************process methods
        void EphemerisPeggedFlybyIn::process_event(const std::vector<doubleType>& X,
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

            if (this->hasTCM)
                this->process_virtual_propellant_constraints(X, Xindex, F, Findex, G, needG);

            this->BoundaryEventBase::process_specialized_constraints(X, Xindex, F, Findex, G, needG);
        }//end process_event

        void EphemerisPeggedFlybyIn::process_event_main(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //base class
            this->EphemerisPeggedArrivalWithVinfinity::process_event_main(X, Xindex, F, Findex, G, needG);

            //TCM
            if (this->hasTCM)
            {
                //perform a TCM with the current stage monoprop system
                this->mySpacecraft->computeChemicalPropulsionPerformance(this->TCM_magnitude, this->state_before_event(6), true, PropulsionSystemChoice::Monoprop); //monoprop

                this->chemical_fuel_used = this->mySpacecraft->getChemFuelConsumedThisManeuver();

                this->ETM(6, 6) = ( (this->state_before_event(6) - this->chemical_fuel_used) / this->state_before_event(6))_GETVALUE;
            }
        }//end process_event_main()

        void EphemerisPeggedFlybyIn::process_event_right_side(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //base class
            this->EphemerisPeggedArrival::process_event_right_side(X, Xindex, F, Findex, G, needG);

            std::get<2>(this->Derivatives_of_StateAfterEvent[this->dIndex_mass_wrt_encodedMass]) = this->ETM(6, 6);
        }//end process_event_right_side()

        void EphemerisPeggedFlybyIn::process_virtual_propellant_constraints(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //call the base class propellant constraints
            BoundaryEventBase::process_virtual_propellant_constraints(X, Xindex, F, Findex, G, needG);

            if (this->hasTCM && needG)
            {
                size_t Gindex = this->Gindices_dVirtualChemicalFuel_TCM_dLeftMass;
                size_t Xindex = this->jGvar->operator[](Gindex);
                double Gentry = (this->chemical_fuel_used / this->state_before_event(6)) _GETVALUE;
                G[Gindex] = this->X_scale_factors->operator[](Xindex) * -Gentry
                    * this->myUniverse->continuity_constraint_scale_factors(6);
            }
        }//end process_virtual_propellant_constraints()

        void EphemerisPeggedFlybyIn::output(std::ofstream& outputfile,
            const double& launchdate,
            size_t& eventcount)
        {
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


            this->mySpacecraft->computePowerState(R_sc_Sun.norm() / this->myOptions->AU, this->state_after_event(7));

            doubleType power = this->mySpacecraft->getAvailablePower();

            //TCM if applicable
            if (this->hasTCM)
            {
                this->write_output_line(outputfile,//outputfile
                    eventcount,//eventcount
                    "TCM",//event_type
                    "deep-space",//event_location
                    0.0,// timestep_size,
                    -1,//flyby_altitude,
                    0,//BdotR
                    0,//BdotT
                    0,//angle1
                    0,//angle2
                    0,//C3
                    this->state_after_event,//state
                    math::Matrix<doubleType>(3, 1, 0.0),//dV
                    math::Matrix<doubleType>(3, 1, 0.0),//ThrustVector
                    this->TCM_magnitude,//dVmag
                    0.0,//Thrust
                    this->mySpacecraft->getMonopropIsp(),//Isp
                    power,//AvailPower
                    0.0,//mdot
                    0,//number_of_active_engines
                    0.0,
                    "none");//active_power)
            }

            //output end of mission if appropriate
            this->EphemerisPeggedArrival::output(outputfile, launchdate, eventcount);
        }//end output()
    }//end namespace BoundaryEvents
}//end namespace EMTG