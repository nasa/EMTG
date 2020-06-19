
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

#include "FreePointIntercept.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        FreePointIntercept::FreePointIntercept(const std::string& name,
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
        }

        //******************************************calcbounds methods

        //calcbounds
        void FreePointIntercept::calcbounds(std::vector<size_t> timeVariables)
        {
            //create bounds for the v-infinity
            double vinf_max = this->myJourneyOptions->final_velocity[1];

            std::vector< std::tuple<double, double> > vinfBounds = std::vector< std::tuple<double, double> >(3, std::make_tuple(-vinf_max, vinf_max));

            this->calcbounds_event_left_side(timeVariables);

            this->calcbounds_event_main(vinfBounds);

            this->calcbounds_event_right_side();

            if (this->hasTCM)
                this->calcbounds_virtual_propellant_constraints();

            this->calcbounds_specialized_constraints();
        }//end calcbounds()
        
        //calcbounds for propellant constraints - TCM
        void FreePointIntercept::calcbounds_virtual_propellant_constraints()
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

        void FreePointIntercept::process_event(const std::vector<doubleType>& X,
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

        void FreePointIntercept::process_virtual_propellant_constraints(const std::vector<doubleType>& X,
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

         //******************************************output methods
        void FreePointIntercept::output(std::ofstream& outputfile,
            const double& launchdate,
            size_t& eventcount)
        {
            this->mySpacecraft->setActiveStage(this->stageIndex);

            std::string event_type = "intercept";

            std::string boundary_name = "free point";

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
                this->Vinfinity_in.dot(this->Vinfinity_in),
                this->state_after_event,
                -this->Vinfinity_in,
                empty3vector,
                this->Vinfinity_in.norm(),
                0.0,
                0.0,
                this->mySpacecraft->getAvailablePower(),
                0.0,
                0,
                0.0,
                "none");
        }
    }//end namespace BoundaryEvents
}//end namespace EMTG