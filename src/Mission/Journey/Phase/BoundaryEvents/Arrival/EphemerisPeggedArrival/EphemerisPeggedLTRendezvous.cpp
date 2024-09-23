
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

#include "EphemerisPeggedLTRendezvous.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        EphemerisPeggedLTRendezvous::EphemerisPeggedLTRendezvous(const std::string& name,
            const size_t& journeyIndex,
            const size_t& phaseIndex,
            size_t& stageIndex,
            Astrodynamics::universe* Universe,
            HardwareModels::Spacecraft* mySpacecraft,
            missionoptions* myOptions) :
                EphemerisPeggedArrival(name,
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
        void EphemerisPeggedLTRendezvous::calcbounds(std::vector<size_t> timeVariables)
        {
            this->calcbounds_event_left_side(timeVariables);

            //no calcbounds_event_main because LT_rendezvous is trivial

            this->calcbounds_event_right_side();

            this->calcbounds_specialized_constraints();
        }//end calcbounds()

        void EphemerisPeggedLTRendezvous::calcbounds_event_right_side()
        {
            //base class
            EphemerisPeggedArrival::calcbounds_event_right_side();

            //no specialized variables
        }

        //******************************************process methods

        void EphemerisPeggedLTRendezvous::process_event(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            this->process_staging();

            this->process_event_left_side(X, Xindex, F, Findex, G, needG);

            this->process_event_right_side(X, Xindex, F, Findex, G, needG);

            this->BoundaryEventBase::process_specialized_constraints(X, Xindex, F, Findex, G, needG);
        }//end process_event

        void EphemerisPeggedLTRendezvous::process_event_right_side(const std::vector<doubleType>& X,
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

        //******************************************output methods
        void EphemerisPeggedLTRendezvous::output(std::ofstream& outputfile,
            const double& launchdate,
            size_t& eventcount)
        {
            this->mySpacecraft->setActiveStage(this->stageIndex);

            std::string event_type = "LT_rndzvs";

            std::string boundary_name = this->myBody->name;

            math::Matrix<doubleType> empty3vector(3, 1, 0.0);

            //where is the Sun?
            math::Matrix<doubleType> R_sc_Sun(3, 1, 0.0);
            if (this->myUniverse->central_body_SPICE_ID == 10)
            {
                R_sc_Sun = this->state_after_event.getSubMatrix1D(0, 2);
            }
            else
            {                //where is the central body relative to the sun?
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

            write_output_line(outputfile,
                eventcount,
                event_type,
                boundary_name,
                this->EventTimeWidth,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                this->state_after_event,
                empty3vector,
                empty3vector,
                0.0,
                0.0,
                0.0,
                this->mySpacecraft->getAvailablePower(),
                0.0,
                0,
                0.0,
                "none");

            //output end of mission if appropriate
            this->EphemerisPeggedArrival::output(outputfile, launchdate, eventcount);
        }//end output()

        void EphemerisPeggedLTRendezvous::output_maneuver_and_target_spec(std::ofstream& maneuver_spec_file, std::ofstream& target_spec_file, bool& haveManeuverNeedTarget)
        {
            //The base-class version of this method assumes that there is NOT an arrival maneuver,
            //and therefore does not write any lines in the maneuver spec
            if (haveManeuverNeedTarget)
            {
                //reset the flag that tells us we need a target spec line
                haveManeuverNeedTarget = false;

                //Step 1: initialize a target spec object
                target_spec_line myTargetSpecLine(this->name,
                    "EME2000",
                    this->myUniverse->central_body.name,
                    this->state_after_event(7),
                    this->state_after_event);

                //Step 2: write target spec object
                myTargetSpecLine.write(target_spec_file);
            }
        }//end output_maneuver_and_target_spec()
    }//end namespace BoundaryEvents
}//end namespace EMTG