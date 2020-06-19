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

#include "PeriapseDeparture.h"
#include "StateRepresentationFactory.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        PeriapseDeparture::PeriapseDeparture(const std::string& name,
                                             const size_t& journeyIndex,
                                             const size_t& phaseIndex,
                                             size_t& stageIndex,
                                             Astrodynamics::universe* Universe,
                                             HardwareModels::Spacecraft* mySpacecraft,
                                             missionoptions* myOptions,
                                             ArrivalEvent* PreviousPhaseArrivalEvent)
        {
            this->initialize(name,
                             journeyIndex,
                             phaseIndex,
                             stageIndex,
                             Universe,
                             mySpacecraft,
                             myOptions,
                             PreviousPhaseArrivalEvent);
        }//end constructor

        void PeriapseDeparture::initialize(const std::string& name,
                                const size_t& journeyIndex,
                                const size_t& phaseIndex,
                                size_t& stageIndex,
                                Astrodynamics::universe* Universe,
                                HardwareModels::Spacecraft* mySpacecraft,
                                missionoptions* myOptions,
                                ArrivalEvent* PreviousPhaseArrivalEvent)
        {

            //set this boundary's state representation
            this->myStateRepresentationEnum = myOptions->PeriapseBoundaryStateRepresentation;

            //periapse arrival cannot use IncomingBplane
            if (this->myStateRepresentationEnum == StateRepresentation::IncomingBplane)
            {
                std::cout << "In Journey " << this->journeyIndex << "'s departure event, the state representation is set to "
                    << StateRepresentationStrings[this->myStateRepresentationEnum == StateRepresentation::IncomingBplane]
                    << ". PeriapseArrival automatically switches this to "
                    << StateRepresentationStrings[this->myStateRepresentationEnum == StateRepresentation::OutgoingBplane] << std::endl;
                this->myStateRepresentationEnum = StateRepresentation::OutgoingBplane;
            }

            this->myStateRepresentation = Astrodynamics::CreateStateRepresentation(this->myStateRepresentationEnum, Universe->mu);

            this->PeriapseBoundary::initialize(name,
                                                journeyIndex,
                                                phaseIndex,
                                                stageIndex,
                                                Universe,
                                                mySpacecraft,
                                                myOptions);

            this->DepartureEvent::initialize(name,
                                             journeyIndex,
                                             phaseIndex,
                                             stageIndex,
                                             Universe,
                                             mySpacecraft,
                                             myOptions,
                                             PreviousPhaseArrivalEvent);

            //set periapse distance bounds
            this->periapseDistanceBounds[0] = this->myJourneyOptions->PeriapseDeparture_altitude_bounds[0] + this->myUniverse->central_body_radius;
            this->periapseDistanceBounds[1] = this->myJourneyOptions->PeriapseDeparture_altitude_bounds[1] + this->myUniverse->central_body_radius;
        }//end initialize()
                               
        
        //******************************************calcbounds methods
        void PeriapseDeparture::calcbounds_event_left_side(const std::vector<double>& RadiusBounds,
                                                           const std::vector<double>& VelocityMagnitudeBounds,
                                                           std::vector<size_t> timeVariables)
        {
            //Step 1: mass bounds
            std::vector<double> MassBounds;            

            if (this->isFirstEventInMission && this->hasFixedInitialMass)
            {
                MassBounds = std::vector<double>({ this->myJourneyOptions->maximum_mass - 1.0e-13, this->myJourneyOptions->maximum_mass });
            }
            else
            {
                MassBounds = std::vector<double>({ 1.0e-13, this->myJourneyOptions->maximum_mass });
            }

            //Step 2: delegate
            this->calcbounds_event_left_side(RadiusBounds, 
                VelocityMagnitudeBounds,
                MassBounds,
                timeVariables);
        }//end calcbounds_event_left_side()

        void PeriapseDeparture::calcbounds_event_left_side(const std::vector<double>& RadiusBounds,
                                                           const std::vector<double>& VelocityMagnitudeBounds,
                                                           const std::vector<double>& MassBounds,
                                                           std::vector<size_t> timeVariables)
        {
            //Step 1: epoch bounds
            if (this->isFirstEventInMission)
            {
                this->Xlowerbounds->push_back(this->myOptions->launch_window_open_date + this->myOptions->Journeys.front().wait_time_bounds[0]);
                this->Xupperbounds->push_back(this->myOptions->launch_window_open_date + this->myOptions->Journeys.front().wait_time_bounds[1]);
                this->X_scale_factors->push_back(1.0 / this->myUniverse->continuity_constraint_scale_factors(7));
                this->Xdescriptions->push_back(prefix + "event left state epoch");
                timeVariables.insert(timeVariables.begin(), this->Xdescriptions->size() - 1);
            }

            //Step 2: base departure class
            this->DepartureEvent::calcbounds_event_left_side();

            //Step 3: base boundary class
            this->PeriapseBoundary::calcbounds_event_left_side(RadiusBounds,
                VelocityMagnitudeBounds,
                timeVariables);
                       
            //mass multipliers
            this->calcbounds_mass_multipliers();
        }//end calcbounds_event_left_side()


        void PeriapseDeparture::calcbounds_event_right_side()
        {
            //base class
            this->PeriapseBoundary::calcbounds_event_right_side();
        }//end calcbounds_event_right_side

        void PeriapseDeparture::calcbounds_specialized_constraints()
        {
            this->DepartureEvent::calcbounds_specialized_constraints();
        }

        //**************************************process functions
        void PeriapseDeparture::process_event_left_side(const std::vector<doubleType>& X,
                                    size_t& Xindex,
                                    std::vector<doubleType>& F,
                                    size_t& Findex,
                                    std::vector<double>& G,
                                    const bool& needG)
        {
            //Step 1: epoch
            this->BoundaryEventBase::process_left_epoch(X, Xindex, F, Findex, G, needG);

            //Step 2: base departure class
            this->DepartureEvent::process_event_left_side(X, Xindex, F, Findex, G, needG);

            //Step 3: base periapse boundary
            this->PeriapseBoundary::process_event_left_side(X, Xindex, F, Findex, G, needG);

            //Step 4: mass continuity
            if (this->hasWaitTime)
                this->process_left_mass_continuity_constraint(X, Xindex, F, Findex, G, needG);

            //Step 5: mass increment
            this->process_mass_multipliers(X, Xindex, F, Findex, G, needG);
        }//end process_event_left_side()

        

        void PeriapseDeparture::process_event_right_side(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG)
        {
            //base ephemeris pegged boundary
            this->PeriapseBoundary::process_event_right_side(X, Xindex, F, Findex, G, needG);
        }//end process_event_right_side()

         //******************************************output methods
        void PeriapseDeparture::output(std::ofstream& outputfile,
            const double& launchdate,
            size_t& eventcount)
        {
            this->mySpacecraft->setActiveStage(this->stageIndex);

            if (this->myOptions->output_dormant_journeys && this->hasWaitTime)
            {
                std::string event_type = "waiting";

                std::string boundary_name = this->myBody->name;

                math::Matrix<doubleType> empty3vector(3, 1, 0.0);

                math::Matrix<doubleType> waitState(8, 1, 0.0);
                doubleType waitEpoch;

                for (size_t step = 0; step < this->myOptions->num_timesteps; ++step)
                {
                    waitEpoch = this->EventLeftEpoch - this->EventWaitTime * ((double)(this->myOptions->num_timesteps - step) / this->myOptions->num_timesteps);

                    //TODO compute wait-state for non-body boundary conditions, involves propagation

                    //where is the Sun?
                    math::Matrix<doubleType> R_sc_Sun(3, 1, 0.0);
                    if (this->myUniverse->central_body_SPICE_ID == 10)
                    {
                        R_sc_Sun = waitState.getSubMatrix1D(0, 2);
                    }
                    else
                    {
                        //where is the central body relative to the sun?
                        doubleType central_body_state_and_derivatives[12];
                        this->myUniverse->locate_central_body(waitState(7),
                            central_body_state_and_derivatives,
                            *this->myOptions,
                            false);                                               
                        math::Matrix<doubleType> R_CB_Sun(3, 1, 0.0);
                        for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
                        {                            R_CB_Sun(stateIndex) = central_body_state_and_derivatives[stateIndex];
                        }                                               
                        R_sc_Sun = waitState.getSubMatrix1D(0, 2) + R_CB_Sun;
                    }
                    this->mySpacecraft->computePowerState(R_sc_Sun.getSubMatrix1D(0, 2).norm() / this->myOptions->AU, waitState(7));

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
                        waitState,
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
            }
        }//end output()

    }//end namespace BoundaryEvents
}//end namespace EMTG