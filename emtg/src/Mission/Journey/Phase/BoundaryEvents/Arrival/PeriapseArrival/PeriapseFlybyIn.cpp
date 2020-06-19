
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

#include "PeriapseFlybyIn.h"
#include "bplane.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        PeriapseFlybyIn::PeriapseFlybyIn(const std::string& name,
            const size_t& journeyIndex,
            const size_t& phaseIndex,
            size_t& stageIndex,
            Astrodynamics::universe* Universe,
            HardwareModels::Spacecraft* mySpacecraft,
            missionoptions* myOptions)
        {
            this->PeriapseArrival::initialize(name,
                                              journeyIndex,
                                              phaseIndex,
                                              stageIndex,
                                              Universe,
                                              mySpacecraft,
                                              myOptions);
        }

        //******************************************calcbounds methods

        //calcbounds
        void PeriapseFlybyIn::calcbounds(std::vector<size_t> timeVariables)
        {
            this->calcbounds_event_left_side(timeVariables);

            //no calcbounds_event_main because PeriapseFlybyIn is trivial

            this->calcbounds_event_right_side();

            this->calcbounds_specialized_constraints();
        }//end calcbounds()

        void PeriapseFlybyIn::calcbounds_event_left_side(std::vector<size_t> timeVariables)
        {
            this->X_index_of_first_decision_variable_in_this_event = this->Xdescriptions->size();

            //Step 1: set the current stage
            this->mySpacecraft->setActiveStage(this->stageIndex);

            //Step 2: radius bounds
            double minRadius, maxRadius;
            if (this->myJourneyOptions->PeriapseArrival_override_altitude)
            {
                minRadius = this->myUniverse->central_body.radius + this->myJourneyOptions->PeriapseArrival_altitude_bounds[0];
                maxRadius = this->myUniverse->central_body.radius + this->myJourneyOptions->PeriapseArrival_altitude_bounds[1];
            }
            else
            {
                minRadius = this->myUniverse->central_body.radius + this->myUniverse->central_body.minimum_safe_flyby_altitude;
                if (this->myUniverse->central_body.mass < 1.0e+25)
                    maxRadius = 10.0 * this->myUniverse->central_body.radius;
                else
                    maxRadius = 300.0 * this->myUniverse->central_body.radius;
            }
            std::vector<double> RadiusBounds({ minRadius, maxRadius });

            //Step 3: v
            std::vector<double> VelocityMagnitudeBounds({ this->myJourneyOptions->final_velocity[0], this->myJourneyOptions->final_velocity[1] });
            
            //Step 4: base class
            this->PeriapseArrival::calcbounds_event_left_side(RadiusBounds, VelocityMagnitudeBounds, timeVariables);
        }//end calcbounds_left_side()

        void PeriapseFlybyIn::calcbounds_event_right_side()
        {
            //base class
            PeriapseArrival::calcbounds_event_right_side();

            //no specialized variables
        }//end calcbounds_right_side()

        //******************************************process methods

        void PeriapseFlybyIn::process_event(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            this->process_staging();

            this->process_staging();

            this->process_event_left_side(X, Xindex, F, Findex, G, needG);

            this->process_event_right_side(X, Xindex, F, Findex, G, needG);

            this->BoundaryEventBase::process_specialized_constraints(X, Xindex, F, Findex, G, needG);
        }//end process_event

        void PeriapseFlybyIn::process_event_left_side(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //base class - nothing else happens
            this->PeriapseArrival::process_event_left_side(X, Xindex, F, Findex, G, needG);
        }//end process_event_left_side()

        //******************************************output methods
        void PeriapseFlybyIn::output(std::ofstream& outputfile,
            const double& launchdate,
            size_t& eventcount)
        {
            this->mySpacecraft->setActiveStage(this->stageIndex);

            std::string event_type = "periapse";

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


            //find bplane parameters
            doubleType BdotR = 0.0;
            doubleType BdotT = 0.0;
            if (this->myUniverse->central_body.spice_ID != 10)
            {
                Astrodynamics::bplane myBplane(this->myUniverse->central_body.mu);
                myBplane.define_bplane(this->state_before_event);
                BdotR = myBplane.getBdotR();
                BdotT = myBplane.getBdotT();
            }

            write_output_line(outputfile,
                eventcount,
                event_type,
                boundary_name,
                this->EventTimeWidth,
                this->state_after_event.getSubMatrix1D(0, 2).norm() - this->myUniverse->central_body.radius,
                BdotR,
                BdotT,
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
        }
    }//end namespace BoundaryEvents
}//end namespace EMTG