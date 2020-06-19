
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

#include "EphemerisReferencedFreeDirectDepartureInterior.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        EphemerisReferencedFreeDirectDepartureInterior::EphemerisReferencedFreeDirectDepartureInterior(const std::string& name,
            const size_t& journeyIndex,
            const size_t& phaseIndex,
            size_t& stageIndex,
            Astrodynamics::universe* Universe,
            HardwareModels::Spacecraft* mySpacecraft,
            missionoptions* myOptions,
            ArrivalEvent* PreviousPhaseArrivalEvent) :
            EphemerisReferencedDepartureInterior::EphemerisReferencedDepartureInterior(name,
                    journeyIndex,
                    phaseIndex,
                    stageIndex,
                    Universe,
                    mySpacecraft,
                    myOptions,
                    PreviousPhaseArrivalEvent)

        {
            if (!this->isFirstEventInMission)
                this->hasWaitTime = true;
        }

        //******************************************calcbounds methods

        //calcbounds
        void EphemerisReferencedFreeDirectDepartureInterior::calcbounds(std::vector<size_t> timeVariables)
        {
            std::vector<double> RAbounds({ -2880.0, 2880.0 });
            std::vector<double> DECbounds({ -90.0, 90.0 });

            std::vector<double> MassBounds(2);
            if (this->hasFixedInitialMass)
            {
                MassBounds[0] = this->myOptions->maximum_mass - 1.0e-13;
                MassBounds[1] = this->myOptions->maximum_mass;
            }
            else
            {

                MassBounds[0] = 1.0e-13;
                MassBounds[1] = this->myOptions->maximum_mass;
            }
            
            //epoch bounds
            std::vector<double> EpochBounds(2);
            if (this->myJourneyOptions->bounded_departure_date)//bounded arrival date
            {
                EpochBounds[0] = this->myJourneyOptions->departure_date_bounds[0];
                EpochBounds[1] = this->myJourneyOptions->departure_date_bounds[1];
            }
            else if (this->isFirstEventInJourney)
            {
                EpochBounds[0] = this->myOptions->launch_window_open_date + this->myJourneyOptions->wait_time_bounds[0];
                EpochBounds[1] = this->myOptions->launch_window_open_date + this->myJourneyOptions->wait_time_bounds[1];
            }
            else
            {
                EpochBounds[0] = this->myOptions->launch_window_open_date;
                EpochBounds[1] = this->myBody->getEphemerisWindowClose();
            }

            this->calcbounds_event_interface_state(RAbounds, DECbounds, MassBounds, EpochBounds, timeVariables);

            this->calcbounds_event_left_side(timeVariables);

            //no calcbounds_event_main because free direct departure is trivial

            this->calcbounds_event_right_side();

            this->calcbounds_specialized_constraints();
        }//end calcbounds()

        void EphemerisReferencedFreeDirectDepartureInterior::calcbounds_event_interface_state(const std::vector<double>& RAbounds,
            const std::vector<double>& DECbounds,
            std::vector<double>& MassBounds,
            const std::vector<double>& EpochBounds,
            std::vector<size_t> timeVariables)
        {
            this->EphemerisReferencedDepartureInterior::calcbounds_event_interface_state(RAbounds, DECbounds, MassBounds, EpochBounds, timeVariables);
        }//end calcbounds_event_interface_state()


        void EphemerisReferencedFreeDirectDepartureInterior::calcbounds_event_left_side(std::vector<size_t> timeVariables)
        {
            //base class
            this->EphemerisReferencedDepartureInterior::calcbounds_event_left_side(timeVariables);

            //no specialized variables
        }//end calcbounds_event_left_side()

        void EphemerisReferencedFreeDirectDepartureInterior::calcbounds_event_right_side()
        {
            //base class
            this->EphemerisReferencedDepartureInterior::calcbounds_event_right_side();

            //no specialized variables
        }//end calcbounds_event_right_side()

        //******************************************process methods

        void EphemerisReferencedFreeDirectDepartureInterior::process_event(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            this->process_staging();

            this->process_event_interface_state(X, Xindex, F, Findex, G, needG);

            this->process_event_left_side(X, Xindex, F, Findex, G, needG);

            this->process_event_right_side(X, Xindex, F, Findex, G, needG);

            this->process_specialized_constraints(X, Xindex, F, Findex, G, needG);
        }//end process_event

        void EphemerisReferencedFreeDirectDepartureInterior::process_event_interface_state(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //velocity on the boundary of the ellipsoid is zero
            for (size_t stateIndex = 3; stateIndex < 6; ++stateIndex)
            {
                this->state_on_interface_spherical(stateIndex) = 0.0;
                this->state_on_interface_cartesian(stateIndex) = 0.0;
            }

            this->EphemerisReferencedDepartureInterior::process_event_interface_state(X, Xindex, F, Findex, G, needG);
        }//end process_event_interface_state()

        void EphemerisReferencedFreeDirectDepartureInterior::process_event_right_side(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //base class
            this->EphemerisReferencedDepartureInterior::process_event_right_side(X, Xindex, F, Findex, G, needG);

            std::get<2>(this->Derivatives_of_StateAfterEvent[this->dIndex_mass_wrt_encodedMass]) = this->ETM(6, 6);
        }//end process_event_right_side()

        //******************************************output methods
        void EphemerisReferencedFreeDirectDepartureInterior::output(std::ofstream& outputfile,
            const double& launchdate,
            size_t& eventcount)
        {
            this->mySpacecraft->setActiveStage(this->stageIndex);

            std::string event_type = "departure";

            std::string boundary_name = this->myUniverse->central_body_name + "_BE"; //"_BE" means "boundary ellipsoid"

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


            //we print state BEFORE event here, and we will print state AFTER event in the next journey
            //that way we get the right reference frame!
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
                this->state_before_event,
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
        }//end output()

    }//end namespace BoundaryEvents
}//end namespace EMTG