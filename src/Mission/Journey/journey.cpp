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

//EMTGv9 Journey class
//7-2-2017

#include <iostream>
#include <fstream>
#include <algorithm>
#include <iterator>

//#include "MGALTS_phase.h"
#include "MGALTphase.h"
#include "FBLTphase.h"
#include "PSFBphase.h"
#include "PSBIphase.h"
#include "MGAnDSMs_phase.h"
#include "CoastPhase.h"
#include "SundmanCoastPhase.h"

#ifdef HAS_PROBEENTRYPHASE
#include "ProbeEntryPhase.h"
#endif

#ifdef HAS_CONTROLLAWTHRUSTPHASE
#include "ControlLawThrustPhase.h"
#endif

#include "journey.h"
#include "missionoptions.h"
#include "universe.h"
#include "EMTG_time_utilities.h"
#include "EMTG_solver_utilities.h"

namespace EMTG
{

    Journey::Journey() :
            number_of_phases(0)
    {
        // default constructor does nothing (is never used)
    }

    Journey::Journey(const size_t& journeyIndex,
        size_t& stageIndex,
        missionoptions* options,
        Astrodynamics::universe* Universe,
        Journey* PreviousJourney,
        Mission* myMission,
        HardwareModels::LaunchVehicle* launchvehicle,
        HardwareModels::Spacecraft* spacecraft)
    {
        //index
        this->journeyIndex = journeyIndex;

        //pointer to previous journey
        this->PreviousJourney = PreviousJourney;

        //pointers to options and universe
        this->myOptions = options;
        this->myJourneyOptions = &(this->myOptions->Journeys[this->journeyIndex]);
        this->myUniverse = Universe;

        //correct the universe's scale factors, specifically mass
        this->myUniverse->continuity_constraint_scale_factors(6) = 1.0 / this->myJourneyOptions->maximum_mass;

        //store pointers to hardware
        this->myLaunchVehicle = launchvehicle;
        this->mySpacecraft = spacecraft;

        //pointer to mission
        this->myMission = myMission;

        //create the phases

        this->phase_type_codes = std::vector<std::string>({ "MGALTS", "FBLTS", "MGALT", "FBLT", "PSBI",  "PSFB", "MGAnDSMs", "CoastPhase", "SundmanCoastPhase", "VARIABLE_PHASE_TYPE", "ProbeEntryPhase", "ControlLawThrustPhase" });

        this->number_of_phases = this->myJourneyOptions->number_of_phases;

        this->boundary_states.resize(this->number_of_phases + 1);

        for (size_t phaseIndex = 0; phaseIndex < this->number_of_phases; ++phaseIndex)
        {
            Phases::phase* previousPhase;
            if (phaseIndex == 0 && this->journeyIndex == 0)
                previousPhase = NULL;
            else if (phaseIndex > 0)
                previousPhase = &this->phases.back();
            else
                previousPhase = &(this->PreviousJourney->phases.back());

            PhaseType phase_type = this->myOptions->mission_type == PhaseType::VARIABLE_PHASE_TYPE ? this->myJourneyOptions->phase_type : this->myOptions->mission_type;

            std::string PhaseName = "j" + std::to_string(this->journeyIndex) + "p" + std::to_string(phaseIndex)
                + this->phase_type_codes[phase_type];

            switch (phase_type)
            {
                case EMTG::PhaseType::MGALTS:
                    {
                        //this phase is an MGA-LT-S phase
                        throw std::invalid_argument("The MGALT-S phase type is not yet implemented. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
                        //this->phases.push_back(new MGA_LT_S_phase(j, p, options, Universe, PreviousPhase, this, this->myLaunchVehicle, this->mySpacecraft));
                        break;
                    }
                case EMTG::PhaseType::FBLTS:
                    {
                        //this phase is an FBLT-S phase
                        throw std::invalid_argument("The FBLT-S phase type is not yet implemented. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
                        //this->phases.push_back(new FBLT_S_phase(j, p, options, Universe));
                        break;
                    }
                case EMTG::PhaseType::MGALT:
                    {
                        //this phase is an MGALT phase
                        this->phases.push_back(new Phases::MGALTphase(PhaseName,
                            this->journeyIndex,
                            phaseIndex,
                            stageIndex,
                            previousPhase,
                            this->myUniverse,
                            this->mySpacecraft,
                            this->myLaunchVehicle,
                            this->myOptions));
                        break;
                    }
                case EMTG::PhaseType::FBLT:
                    {
                        //this phase is an FBLT phase
                        this->phases.push_back(new Phases::FBLTphase(PhaseName,
                            this->journeyIndex,
                            phaseIndex,
                            stageIndex,
                            previousPhase,
                            this->myUniverse,
                            this->mySpacecraft,
                            this->myLaunchVehicle,
                            this->myOptions));
                        break;
                    }
                case EMTG::PhaseType::PSBI:
                {
                    //this phase is a PSBIphase
                    this->phases.push_back(new Phases::PSBIphase(PhaseName,
                        this->journeyIndex,
                        phaseIndex,
                        stageIndex,
                        previousPhase,
                        this->myUniverse,
                        this->mySpacecraft,
                        this->myLaunchVehicle,
                        this->myOptions));
                    break;
                }
                case EMTG::PhaseType::PSFB:
                {
                    //this phase is a PSFBphase
                    this->phases.push_back(new Phases::PSFBphase(PhaseName,
                        this->journeyIndex,
                        phaseIndex,
                        stageIndex,
                        previousPhase,
                        this->myUniverse,
                        this->mySpacecraft,
                        this->myLaunchVehicle,
                        this->myOptions));
                    break;
                }
                case EMTG::PhaseType::MGAnDSMs:
                    {
                        //this phase is an MGA-nDSM-s phase
                        this->phases.push_back(new Phases::MGAnDSMs_phase(PhaseName,
                            this->journeyIndex,
                            phaseIndex,
                            stageIndex,
                            previousPhase,
                            this->myUniverse,
                            this->mySpacecraft,
                            this->myLaunchVehicle,
                            this->myOptions));
                        break;
                    }
                case EMTG::PhaseType::CoastPhase:
                {
                    //this phase is an CoastPhase
                    this->phases.push_back(new Phases::CoastPhase(PhaseName,
                        this->journeyIndex,
                        phaseIndex,
                        stageIndex,
                        previousPhase,
                        this->myUniverse,
                        this->mySpacecraft,
                        this->myLaunchVehicle,
                        this->myOptions));
                    break;
                }
                case EMTG::PhaseType::SundmanCoastPhase:
                {
                    //this phase is a SundmanCoastPhase
                    this->phases.push_back(new Phases::SundmanCoastPhase(PhaseName,
                        this->journeyIndex,
                        phaseIndex,
                        stageIndex,
                        previousPhase,
                        this->myUniverse,
                        this->mySpacecraft,
                        this->myLaunchVehicle,
                        this->myOptions));
                    break;
                }
                case EMTG::PhaseType::ProbeEntryPhase:
                {
                    //this phase is a ProbeEntryPhase
#ifdef HAS_PROBEENTRYPHASE
                    this->phases.push_back(new Phases::ProbeEntryPhase(PhaseName,
                        this->journeyIndex,
                        phaseIndex,
                        stageIndex,
                        previousPhase,
                        this->myUniverse,
                        this->mySpacecraft,
                        this->myLaunchVehicle,
                        this->myOptions));
#else
                    throw std::invalid_argument("The ProbeEntryPhase phase type is not included in this version of EMTG. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
#endif
                    break;
                }
                case EMTG ::PhaseType::ControlLawThrustPhase:
                {
                    //this phase is a SundmanCoastPhase
#ifdef HAS_CONTROLLAWTHRUSTPHASE
                    this->phases.push_back(new Phases::ControlLawThrustPhase(PhaseName,
                        this->journeyIndex,
                        phaseIndex,
                        stageIndex,
                        previousPhase,
                        this->myUniverse,
                        this->mySpacecraft,
                        this->myLaunchVehicle,
                        this->myOptions));
#else
                    throw std::invalid_argument("The ControlLawThrustPhase phase type is not included in this version of EMTG. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
#endif
                    break;
                }
                case EMTG::PhaseType::VARIABLE_PHASE_TYPE:
                {
                    //Variable phase type
                    throw std::invalid_argument("Variable phase type is not a phase type, it's a metaphysical theme. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
                    break;
                }
            }
        }
    }//end Journey constructor

    void Journey::calcbounds(std::vector<double>& Xupperbounds,
        std::vector<double>& Xlowerbounds,
        std::vector<double>& X_scale_factors,
        std::vector<double>& Fupperbounds,
        std::vector<double>& Flowerbounds,
		std::vector<double>& F_scale_factors,
        std::vector<std::string>& Xdescriptions,
        std::vector<std::string>& Fdescriptions,
        std::vector<size_t>& iGfun,
        std::vector<size_t>& jGvar,
        std::vector<std::string>& Gdescriptions,
        std::vector<size_t>& iAfun,
        std::vector<size_t>& jAvar,
        std::vector<std::string>& Adescriptions,
        std::vector<double>& A,
        std::vector<double>& synodic_periods)
    {
        std::string prefix = "j" + std::to_string(this->journeyIndex) + ": ";
        this->first_X_entry_in_journey = Xdescriptions.size();
        this->first_F_entry_in_journey = Fdescriptions.size();
        this->first_G_entry_in_journey = Gdescriptions.size();

        this->setup_calcbounds(&Xupperbounds,
            &Xlowerbounds,
            &X_scale_factors,
            &Fupperbounds,
            &Flowerbounds,
			&F_scale_factors,
            &Xdescriptions,
            &Fdescriptions,
            &iGfun,
            &jGvar,
            &Gdescriptions,
            &iAfun,
            &jAvar,
            &Adescriptions,
            &A);

        //phase-specific bounds and constraints
        for (size_t phaseIndex = 0; phaseIndex < this->number_of_phases; ++phaseIndex)
        {
            this->phases[phaseIndex].setup_calcbounds(&Xupperbounds,
                &Xlowerbounds,
                &X_scale_factors,
                &Fupperbounds,
                &Flowerbounds,
				&F_scale_factors,
                &Xdescriptions,
                &Fdescriptions,
                &iGfun,
                &jGvar,
                &Gdescriptions,
                &iAfun,
                &jAvar,
                &Adescriptions,
                &A);

            
            this->phases[phaseIndex].calcbounds();

            synodic_periods.push_back(this->phases[phaseIndex].getsynodic_period());
        }

    
        //journey time constraints

        if (this->myJourneyOptions->bounded_departure_date)
        {
            this->Flowerbounds->push_back((this->myJourneyOptions->departure_date_bounds[0] - this->myOptions->launch_window_open_date) / (this->myJourneyOptions->departure_date_bounds[1] - this->myOptions->launch_window_open_date) - 1.0);
            this->Fupperbounds->push_back(0.0);
            this->Fdescriptions->push_back(prefix + "journey departure date bounds");

            if (this->journeyIndex == 0)//if this constraint is applied to the first journey, it's special because there are no preceding time entries to search
            {
                //launch epoch is the first entry in the first phase's departure event's time vector
                size_t Xindex = this->phases.front().getDepartureEvent()->get_Xindices_EventLeftEpoch().front();
                this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                    Xindex,
                    this->departure_date_constraint_G_indices);
            }
            else
            {
                //we need all time variables from the previous journeys
                std::vector<size_t> timeVariables = this->PreviousJourney->getLastPhase()->getArrivalEvent()->get_Xindices_EventRightEpoch();
                for (size_t Xindex : timeVariables)
                {
                    this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                        Xindex,
                        this->departure_date_constraint_G_indices);
                }

                //we also need the current journey wait time if applicable
                if (this->phases.front().getDepartureEvent()->getHasWaitTime())
                {
                    //wait time is the first entry in the first phase's departure event's time vector
                    size_t Xindex = this->phases.front().getDepartureEvent()->get_Xindices_EventLeftEpoch().front();

                    this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                        Xindex,
                        this->departure_date_constraint_G_indices);

                }
            }
        }//end bounded departure date

        if (this->myJourneyOptions->timebounded == 1 //bounded journey flight time
            || this->myJourneyOptions->timebounded == 3) //bounded aggregate flight time
        {
            this->Flowerbounds->push_back(this->myJourneyOptions->flight_time_bounds[0] / this->myJourneyOptions->flight_time_bounds[1] - 1);
            this->Fupperbounds->push_back(0.0);
            this->Fdescriptions->push_back(prefix + "journey flight time bounds");

            if (this->myJourneyOptions->timebounded == 3)
            {
                //all time variables to the end of this journey - the final phase's arrival event's right side time vector
                //but NOT the launch epoch. Which means not the first variable, what ever it is.
                std::vector<size_t> timeVariables = this->phases.back().getArrivalEvent()->get_Xindices_EventRightEpoch();
                timeVariables.erase(timeVariables.begin());
                for (size_t Xindex : timeVariables)
                {
                    this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                        Xindex,
                        this->timeconstraints_G_indices);
                }
            }
            else //this->myJourneyOptions->timebounded == 1
            {
                //this is trickier. We're going to take the first event's and the last event's time vectors and then keep only the things that are unique
                std::vector<size_t> left_timeVariables = this->phases.front().getDepartureEvent()->get_Xindices_EventLeftEpoch();
                std::vector<size_t> right_timeVariables = this->phases.back().getArrivalEvent()->get_Xindices_EventRightEpoch();

                std::vector<size_t> timeVariables;

                //using std::set_symmetric_difference because it's cool - have to sort everything first (should be pre-sorted but let's be safe)
                std::sort(left_timeVariables.begin(), left_timeVariables.end());
                std::sort(right_timeVariables.begin(), right_timeVariables.end());

                std::set_symmetric_difference(
                    left_timeVariables.begin(), left_timeVariables.end(),
                    right_timeVariables.begin(), right_timeVariables.end(),
                    std::back_inserter(timeVariables));

                for (size_t Xindex : timeVariables)
                {
                    this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                        Xindex,
                        this->timeconstraints_G_indices);
                }
            }
        }//end journey flight time bounds
        else if (this->myJourneyOptions->timebounded == 2)//bounded arrival date
        {
            this->Flowerbounds->push_back((this->myJourneyOptions->arrival_date_bounds[0] - this->myOptions->launch_window_open_date) / (this->myJourneyOptions->arrival_date_bounds[1] - this->myOptions->launch_window_open_date) - 1.0);
            this->Fupperbounds->push_back(0.0);
            this->Fdescriptions->push_back(prefix + "journey arrival date bounds");

            //all time variables to the end of this journey - the final phase's arrival event's right side time vector
            std::vector<size_t> timeVariables = this->phases.back().getArrivalEvent()->get_Xindices_EventRightEpoch();
            for (size_t Xindex : timeVariables)
            {
                this->create_sparsity_entry(this->Fdescriptions->size() - 1,
                    Xindex,
                    this->timeconstraints_G_indices);
            }
        }//end bounded arrival date
    }//end Journey::calcbounds()

    //evaluate function
    //return 0 if successful, 1 if failure
    void Journey::process_journey(const std::vector<doubleType>& X,
        size_t& Xindex,
        std::vector<doubleType>& F,
        size_t& Findex,
        std::vector<double>& G,
        const bool& needG)

    {
        //process all of the phases
        for (size_t phaseIndex = 0; phaseIndex < this->number_of_phases; ++phaseIndex)
        {
            //TODO: if this is NOT the first journey, transform the coordinate system into the central body frame of the current journey
            //maybe this should be done in the departure event?

            this->phases[phaseIndex].process_phase(X, Xindex, F, Findex, G, needG);

            //TODO: at the end of the journey, transform the coordinate system back into the central body frame of the Sun (EMTG global reference frame)
            //maybe this should be done in the arrival event?
        }

        //journey-level constraints

        //journey time constraints

        if (this->myJourneyOptions->bounded_departure_date)
        {
            doubleType departure_date = 0.0;

            for (size_t entry = 0; entry < this->departure_date_constraint_G_indices.size(); ++entry)
            {
                size_t Gindex = this->departure_date_constraint_G_indices[entry];
                size_t Xindex = this->jGvar->operator[](Gindex);

                departure_date += X[Xindex];
            }

            F[Findex++] = (departure_date - this->myOptions->launch_window_open_date) / (this->myJourneyOptions->departure_date_bounds[1] - this->myOptions->launch_window_open_date) - 1.0;

            if (needG)
            {
                //the first n-1 derivatives are with respect to flight times
                for (size_t entry = 0; entry < this->departure_date_constraint_G_indices.size(); ++entry)
                {
                    size_t Gindex = this->departure_date_constraint_G_indices[entry];
                    size_t Xindex = this->jGvar->operator[](Gindex);
                    G[Gindex] = this->X_scale_factors->operator[](Xindex)
                        / (this->myJourneyOptions->departure_date_bounds[1] - this->myOptions->launch_window_open_date);
                }
            }
        }//end bounded departure date

        if (this->myJourneyOptions->timebounded)
        {
            if (this->myJourneyOptions->timebounded == 1 //bounded flight time
                || this->myJourneyOptions->timebounded == 3)//bounded aggregate flight time
            {
                F[Findex] = 0.0;
                for (size_t timeIndex = 0; timeIndex < this->timeconstraints_G_indices.size(); ++timeIndex)
                {
                    size_t Gindex = this->timeconstraints_G_indices[timeIndex];
                    size_t Xindex = this->jGvar->operator[](Gindex);

                    F[Findex] += X[Xindex] / this->myJourneyOptions->flight_time_bounds[1];
                }
                F[Findex++] -= 1.0;

                //derivative with respect to journey flight time variables
                if (needG)
                {
                    for (size_t timeIndex = 0; timeIndex < this->timeconstraints_G_indices.size(); ++timeIndex)
                    {
                        size_t Gindex = this->timeconstraints_G_indices[timeIndex];
                        size_t Xindex = this->jGvar->operator[](Gindex);

                        G[Gindex] = this->X_scale_factors->operator[](Xindex)
                            / this->myJourneyOptions->flight_time_bounds[1];
                    }
                }
            }
            else if (this->myJourneyOptions->timebounded == 2) //bounded arrival date
            {
                doubleType arrival_date = 0.0;

                for (size_t entry = 0; entry < this->timeconstraints_G_indices.size(); ++entry)
                {
                    size_t Gindex = this->timeconstraints_G_indices[entry];
                    size_t Xindex = this->jGvar->operator[](Gindex);

                    arrival_date += X[Xindex];
                }

                F[Findex++] = ((arrival_date - this->myOptions->launch_window_open_date) / (this->myJourneyOptions->arrival_date_bounds[1] - this->myOptions->launch_window_open_date) - 1.0);

                if (needG)
                {
                    //the first n-1 derivatives are with respect to flight times
                    for (size_t entry = 0; entry < this->timeconstraints_G_indices.size(); ++entry)
                    {
                        size_t Gindex = this->timeconstraints_G_indices[entry];
                        size_t Xindex = this->jGvar->operator[](Gindex);
                        G[Gindex] = this->X_scale_factors->operator[](Xindex)
                            / (this->myJourneyOptions->arrival_date_bounds[1] - this->myOptions->launch_window_open_date);
                    }
                }
            }
        }//end journey flight time constraints
    }//end process_journey()

    //output functions
    void Journey::output_journey_header(std::ofstream& outputfile,
        size_t& jprint)
    {
        //first output a bunch of header stuff


        outputfile.precision(20);

        outputfile << std::endl;
        outputfile << "Journey: " << jprint << std::endl;
        outputfile << "Journey name: " << this->myJourneyOptions->journey_name << std::endl;
        outputfile << "Central Body: " << this->myUniverse->central_body_name << std::endl;
        outputfile << "Radius (km): " << this->myUniverse->central_body_radius << std::endl;
        outputfile << "mu (km^3/s^2): " << this->myUniverse->mu << std::endl;
        outputfile << "Characteristic length unit (km): " << this->myUniverse->LU << std::endl;
        std::vector<std::string> FrameDefinitions({ "ICRF","J2000_BCI","J2000_BCF","TrueOfDate_BCI","TrueOfDate_BCF" });
        outputfile << "Frame: " << FrameDefinitions[this->myOptions->output_file_frame] << std::endl;

        if (this->myOptions->output_file_frame > ReferenceFrame::ICRF)
        {
            outputfile << "alpha0: " << this->myUniverse->LocalFrame.get_alpha0() << std::endl;

            if (this->myOptions->output_file_frame == ReferenceFrame::TrueOfDate_BCI || this->myOptions->output_file_frame == ReferenceFrame::TrueOfDate_BCF)
                outputfile << "alphadot: " << this->myUniverse->LocalFrame.get_alphadot() << std::endl;

            outputfile << "delta0: " << this->myUniverse->LocalFrame.get_delta0() << std::endl;

            if (this->myOptions->output_file_frame == ReferenceFrame::TrueOfDate_BCI || this->myOptions->output_file_frame == ReferenceFrame::TrueOfDate_BCF)
                outputfile << "deltadot: " << this->myUniverse->LocalFrame.get_deltadot() << std::endl;
            
            if (this->myOptions->output_file_frame == ReferenceFrame::J2000_BCF || this->myOptions->output_file_frame == ReferenceFrame::TrueOfDate_BCF)
            {
                outputfile << "W0: " << this->myUniverse->LocalFrame.getW0() << std::endl;

                if (this->myOptions->output_file_frame == ReferenceFrame::TrueOfDate_BCF)
                    outputfile << "Wdot: " << this->myUniverse->LocalFrame.getWdot() << std::endl;
            }
        }

        if (this->myJourneyOptions->phase_type == PhaseType::MGALT
            || this->myJourneyOptions->phase_type == PhaseType::FBLT
            || this->myJourneyOptions->phase_type == PhaseType::MGALTS
            || this->myJourneyOptions->phase_type == PhaseType::FBLTS
            || this->myJourneyOptions->phase_type == PhaseType::PSFB
            || this->myJourneyOptions->phase_type == PhaseType::PSBI)
            outputfile << "Thruster duty cycle: " << (this->myJourneyOptions->override_duty_cycle ? this->myJourneyOptions->duty_cycle : this->myOptions->engine_duty_cycle) << std::endl;
        outputfile << std::endl;

        //next, column headers

        //column headers line 1
        outputfile.width(5); outputfile << "#";
        outputfile.width(3); outputfile << " | ";
        outputfile.width(16); outputfile << "JulianDate";
        outputfile.width(3); outputfile << " | ";
        outputfile.width(11); outputfile << "MM/DD/YYYY";
        outputfile.width(3); outputfile << " | ";
        outputfile.width(12); outputfile << "event type";
        outputfile.width(3); outputfile << " | ";
        outputfile.width(25); outputfile << "location";
        outputfile.width(3); outputfile << " | ";
        outputfile.width(15); outputfile << "step size";
        outputfile.width(3); outputfile << " | ";
        outputfile.width(19); outputfile << "altitude";
        outputfile.width(3); outputfile << " | ";
        outputfile.width(19); outputfile << "BdotR";
        outputfile.width(3); outputfile << " | ";
        outputfile.width(19); outputfile << "BdotT";
        outputfile.width(3); outputfile << " | ";
        outputfile.width(8); outputfile << "RA";
        outputfile.width(3); outputfile << " | ";
        outputfile.width(8); outputfile << "DEC";
        outputfile.width(3); outputfile << " | ";
        outputfile.width(14); outputfile << "C3";
        outputfile.width(3); outputfile << " | ";
        outputfile.width(19); outputfile << " x";
        outputfile.width(3); outputfile << " | ";
        outputfile.width(19); outputfile << " y";
        outputfile.width(3); outputfile << " | ";
        outputfile.width(19); outputfile << " z";
        outputfile.width(3); outputfile << " | ";
        outputfile.width(19); outputfile << " xdot";
        outputfile.width(3); outputfile << " | ";
        outputfile.width(19); outputfile << " ydot";
        outputfile.width(3); outputfile << " | ";
        outputfile.width(19); outputfile << " zdot";
        outputfile.width(3); outputfile << " | ";
        outputfile.width(19); outputfile << " dV_x";
        outputfile.width(3); outputfile << " | ";
        outputfile.width(19); outputfile << " dV_y";
        outputfile.width(3); outputfile << " | ";
        outputfile.width(19); outputfile << " dV_z";
        outputfile.width(3); outputfile << " | ";
        outputfile.width(19); outputfile << " T_x";
        outputfile.width(3); outputfile << " | ";
        outputfile.width(19); outputfile << " T_y";
        outputfile.width(3); outputfile << " | ";
        outputfile.width(19); outputfile << " T_z";
        outputfile.width(3); outputfile << " | ";
        outputfile.width(17); outputfile << "|dV| (km/s)";
        outputfile.width(3); outputfile << " | ";
        outputfile.width(14); outputfile << "Avail. Thrust";
        outputfile.width(3); outputfile << " | ";
        outputfile.width(14); outputfile << "Isp";
        outputfile.width(3); outputfile << " | ";
        outputfile.width(14); outputfile << "Avail. Power";
        outputfile.width(3); outputfile << " | ";
        outputfile.width(19); outputfile << "Mass Flow";
        outputfile.width(3); outputfile << " | ";
        outputfile.width(14); outputfile << "mass";
        outputfile.width(3); outputfile << " | ";
        outputfile.width(14); outputfile << "number of";
        outputfile.width(3); outputfile << " | ";
        outputfile.width(14); outputfile << "active power";
        outputfile.width(3); outputfile << " | ";
        outputfile.width(14); outputfile << "throttle level";
        outputfile.width(3); outputfile << " | ";
        outputfile << std::endl;

        //column headers line 2
        outputfile.width(5); outputfile << "";
        outputfile.width(3); outputfile << " | ";
        outputfile.width(16); outputfile << " (ET/TDB)";
        outputfile.width(3); outputfile << " | ";
        outputfile.width(11); outputfile << "";
        outputfile.width(3); outputfile << " | ";
        outputfile.width(12); outputfile << "";
        outputfile.width(3); outputfile << " | ";
        outputfile.width(25); outputfile << "";
        outputfile.width(3); outputfile << " | ";
        outputfile.width(15); outputfile << "(days)";
        outputfile.width(3); outputfile << " | ";
        outputfile.width(19); outputfile << "(km)";
        outputfile.width(3); outputfile << " | ";
        outputfile.width(19); outputfile << "(km)";
        outputfile.width(3); outputfile << " | ";
        outputfile.width(19); outputfile << "(km)";
        outputfile.width(3); outputfile << " | ";
        outputfile.width(8); outputfile << "degrees";
        outputfile.width(3); outputfile << " | ";
        outputfile.width(8); outputfile << "degrees";
        outputfile.width(3); outputfile << " | ";
        outputfile.width(14); outputfile << "(km^2/s^2)";
        outputfile.width(3); outputfile << " | ";
        outputfile.width(19); outputfile << "(km)";
        outputfile.width(3); outputfile << " | ";
        outputfile.width(19); outputfile << "(km)";
        outputfile.width(3); outputfile << " | ";
        outputfile.width(19); outputfile << "(km)";
        outputfile.width(3); outputfile << " | ";
        outputfile.width(19); outputfile << "(km/s)";
        outputfile.width(3); outputfile << " | ";
        outputfile.width(19); outputfile << "(km/s)";
        outputfile.width(3); outputfile << " | ";
        outputfile.width(19); outputfile << "(km/s)";
        outputfile.width(3); outputfile << " | ";
        outputfile.width(19); outputfile << "(km/s)";
        outputfile.width(3); outputfile << " | ";
        outputfile.width(19); outputfile << "(km/s)";
        outputfile.width(3); outputfile << " | ";
        outputfile.width(19); outputfile << "(km/s)";
        outputfile.width(3); outputfile << " | ";
        outputfile.width(19); outputfile << "(N)";
        outputfile.width(3); outputfile << " | ";
        outputfile.width(19); outputfile << "(N)";
        outputfile.width(3); outputfile << " | ";
        outputfile.width(19); outputfile << "(N)";
        outputfile.width(3); outputfile << " | ";
        outputfile.width(17); outputfile << "throttle (0-1)";
        outputfile.width(3); outputfile << " | ";
        outputfile.width(14); outputfile << "(N)";
        outputfile.width(3); outputfile << " | ";
        outputfile.width(14); outputfile << "(s)";
        outputfile.width(3); outputfile << " | ";
        outputfile.width(14); outputfile << "(kW)";
        outputfile.width(3); outputfile << " | ";
        outputfile.width(19); outputfile << "rate (kg/s)";
        outputfile.width(3); outputfile << " | ";
        outputfile.width(14); outputfile << "(kg)";
        outputfile.width(3); outputfile << " | ";
        outputfile.width(14); outputfile << "active engines";
        outputfile.width(3); outputfile << " | ";
        outputfile.width(14); outputfile << "(kW)";
        outputfile.width(3); outputfile << " | ";
        outputfile.width(14); outputfile << "";
        outputfile.width(3); outputfile << " | ";
        outputfile << std::endl;


        for (size_t k = 0; k < 627; ++k)
            outputfile << "-";
        outputfile << std::endl;
    }

    void Journey::output(std::ofstream& outputfile,
        size_t& jprint,
        size_t& eventcount)
    {
        //print this journey
        //start with the header
        this->output_journey_header(outputfile, jprint);
        ++jprint;

        for (int phaseIndex = 0; phaseIndex < this->number_of_phases; ++phaseIndex)
        {
            this->phases[phaseIndex].output(outputfile, eventcount);
        }

        //skip two lines
        outputfile << std::endl;
        outputfile << std::endl;

        outputfile.precision(20);

        //print the departure mass increment
        this->phases.front().getDepartureEvent()->output_mass_increment(outputfile);
        this->phases.back().getArrivalEvent()->output_mass_increment(outputfile);

        if (this->myJourneyOptions->journey_end_deltav > 0.0)
            this->phases.back().getArrivalEvent()->output_post_arrival_maneuver(outputfile);
    
        //**************************************************************************************************************
        //print the boundary states
        //create a 3-element storage vector that will be used every time something needs to be rotated to the local frame
        math::Matrix<doubleType> r(3,1);
        math::Matrix<doubleType> v(3,1);
        math::Matrix<doubleType> disp_r(3,1);
        math::Matrix<doubleType> disp_v(3,1);

        //construct the rotation matrix from ICRF to the Universe J2000 BCI frame
        outputfile << "Boundary states:" << std::endl;
        for (size_t i = 0; i < this->number_of_phases + 1; ++i)
        {
            this->boundary_states[i] = i < this->number_of_phases ?
                this->phases[i].getDepartureEvent()->get_boundary_state() //we want the state BEFORE the departure maneuver
                : this->phases[i - 1].getArrivalEvent()->get_boundary_state(); //we want the state AFTER the arrival maneuver
            outputfile << "Boundary: " << i + 1;
            outputfile.precision(8);
            r(0) = boundary_states[i](0);
            r(1) = boundary_states[i](1);
            r(2) = boundary_states[i](2);
            v(0) = boundary_states[i](3);
            v(1) = boundary_states[i](4);
            v(2) = boundary_states[i](5);
            
            this->myUniverse->LocalFrame.rotate_frame_to_frame(ReferenceFrame::ICRF, r, this->myOptions->output_file_frame, disp_r);
            this->myUniverse->LocalFrame.rotate_frame_to_frame(ReferenceFrame::ICRF, v, this->myOptions->output_file_frame, disp_v);

            for (int j = 0; j < 3; ++j)
                outputfile << " " << disp_r(j) _GETVALUE;
            for (int j = 0; j < 3; ++j)
                outputfile << " " << disp_v(j) _GETVALUE;
            outputfile << std::endl;
        }
        outputfile << std::endl;

        //**************************************************************************************************************
        //output the flyby periapse states
        size_t flybyIndex = 1;
        for (size_t phaseIndex = 0; phaseIndex < this->number_of_phases; ++phaseIndex)
        {
            this->phases[phaseIndex].getDepartureEvent()->output_periapse_state(flybyIndex, outputfile);
        }

        //print any information about specialized constraints

        for (size_t phaseIndex = 0; phaseIndex < this->number_of_phases; ++phaseIndex)
        {
            this->phases[phaseIndex].getDepartureEvent()->output_specialized_constraints(outputfile);
            this->phases[phaseIndex].getArrivalEvent()->output_specialized_constraints(outputfile);
        }

        //skip a line, print the flight time, and skip one more line
        outputfile << std::endl;
        doubleType flightTime = this->phases.back().getArrivalEvent()->get_state_after_event()(7) - this->phases.front().getDepartureEvent()->get_state_before_event()(7);
        outputfile << "Journey flight time (seconds): " << flightTime << std::endl;
        outputfile << "Journey flight time (days): " << flightTime / 86400.0 << std::endl;
        outputfile << std::endl;

        outputfile << "End journey" << std::endl;
    }


    void Journey::output_ephemeris(std::ofstream& outputfile, std::ofstream& acceleration_model_file)
    {
        for (size_t phaseIndex = 0; phaseIndex < this->number_of_phases; ++phaseIndex)
            this->phases[phaseIndex].output_ephemeris(outputfile, acceleration_model_file);
    }//end output_ephemeris()

    void Journey::output_STMs()
    {
        for (size_t phaseIndex = 0; phaseIndex < this->number_of_phases; ++phaseIndex)
            this->phases[phaseIndex].output_STMs();
    }//end output_STMs

    void Journey::output_maneuver_and_target_spec(std::ofstream& maneuver_spec_file, std::ofstream& target_spec_file, bool& haveManeuverNeedTarget)
    {
        for (size_t phaseIndex = 0; phaseIndex < this->number_of_phases; ++phaseIndex)
            this->phases[phaseIndex].output_maneuver_and_target_spec(maneuver_spec_file, target_spec_file, haveManeuverNeedTarget);
    }//end output_maneuver_and_target_spec()
       
    doubleType Journey::getDeterministicDeltav() const
    {
        doubleType temp_deltav = 0.0;

        for (size_t phaseIndex = 0; phaseIndex < this->number_of_phases; ++phaseIndex)
            temp_deltav += this->phases[phaseIndex].getDeterministicDeltav();
        
        return temp_deltav;
    }

    doubleType Journey::getStatisticalDeltav() const
    {
        doubleType temp_deltav = 0.0;

        for (size_t phaseIndex = 0; phaseIndex < this->number_of_phases; ++phaseIndex)
            temp_deltav += this->phases[phaseIndex].getStatisticalDeltav();

        return temp_deltav;
    }
} /* namespace EMTG */
