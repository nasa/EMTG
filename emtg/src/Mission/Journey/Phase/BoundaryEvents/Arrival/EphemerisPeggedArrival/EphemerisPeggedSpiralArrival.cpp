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

//Edelbaum spiral departure
//Jacob Englander 10-24-2017

#include "EphemerisPeggedSpiralArrival.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        
        EphemerisPeggedSpiralArrival::EphemerisPeggedSpiralArrival(const std::string& name,
            const size_t& journeyIndex,
            const size_t& phaseIndex,
            size_t& stageIndex,
            Astrodynamics::universe* Universe,
            HardwareModels::Spacecraft* mySpacecraft,
            missionoptions* myOptions) :
            EphemerisPeggedArrival::EphemerisPeggedArrival(name,
                journeyIndex,
                phaseIndex,
                stageIndex,
                Universe,
                mySpacecraft,
                myOptions)
        {
            //create EdelbaumSpiral
            this->mySpiral = new EdelbaumSpiral(name,
                this->journeyIndex,
                this->phaseIndex,
                this->stageIndex,
                this->myUniverse,
                this->myBody,
                this->mySpacecraft,
                this->myOptions);

            //set up spiral
            this->mySpiral->setSpiralStartRadius(this->myJourneyOptions->capture_spiral_starting_radius);
            this->mySpiral->setSpiralEndRadius(this->myJourneyOptions->capture_spiral_final_radius);
            this->mySpiral->setDirection(-1.0);

            if (this->myJourneyOptions->capture_spiral_final_radius > this->myJourneyOptions->capture_spiral_starting_radius)
                throw std::invalid_argument(this->name + " capture spiral final radius is greater than capture spiral starting radius. This is not physical. Please change your input file.");

            this->mySpiral->initialize_spiral_segments();

            this->EventDeterministicDeltav = this->mySpiral->getTotalDeltav();
        }//end constructor

        EphemerisPeggedSpiralArrival::~EphemerisPeggedSpiralArrival()
        {
            delete this->mySpiral;
        }//end ~EphemerisPeggedSpiralArrival()

        void EphemerisPeggedSpiralArrival::setup_calcbounds(
                std::vector<double>* Xupperbounds,
                std::vector<double>* Xlowerbounds,
                std::vector<double>* X_scale_factors,
                std::vector<double>* Fupperbounds,
                std::vector<double>* Flowerbounds,
                std::vector<double>* F_scale_factors,
                std::vector<std::string>* Xdescriptions,
                std::vector<std::string>* Fdescriptions,
                std::vector<size_t>* iGfun,
                std::vector<size_t>* jGvar,
                std::vector<std::string>* Gdescriptions,
                std::vector<size_t>* iAfun,
                std::vector<size_t>* jAvar,
                std::vector<std::string>* Adescriptions,
                std::vector<double>* A)
        {
            this->BoundaryEventBase::setup_calcbounds(
                Xupperbounds,
                Xlowerbounds,
                X_scale_factors,
                Fupperbounds,
                Flowerbounds,
                F_scale_factors,
                Xdescriptions,
                Fdescriptions,
                iGfun,
                jGvar,
                Gdescriptions,
                iAfun,
                jAvar,
                Adescriptions,
                A);

            this->mySpiral->setup_calcbounds(
                Xupperbounds,
                Xlowerbounds,
                X_scale_factors,
                Fupperbounds,
                Flowerbounds,
                F_scale_factors,
                Xdescriptions,
                Fdescriptions,
                iGfun,
                jGvar,
                Gdescriptions,
                iAfun,
                jAvar,
                Adescriptions,
                A);
        }

        void EphemerisPeggedSpiralArrival::calcbounds(std::vector<size_t> timeVariables)
        {
            this->calcbounds_event_left_side(timeVariables);

            this->calcbounds_event_main();

            this->calcbounds_event_right_side();

            this->calcbounds_specialized_constraints();
        }//end calcbounds()

        void EphemerisPeggedSpiralArrival::calcbounds_event_main()
        {
            this->mySpiral->calcbounds(this->Xindices_EventLeftEpoch);
        }//end calcbounds_event_main()

        void EphemerisPeggedSpiralArrival::calcbounds_event_right_side()
        {
            //base class
            this->EphemerisPeggedArrival::calcbounds_event_right_side();

            //but the derivatives of state after event actually come from the spiral
            this->Derivatives_of_StateAfterEvent = this->mySpiral->get_Derivatives_of_StateAfterSpiral();
            this->Derivatives_of_StateAfterEvent_wrt_Time = this->mySpiral->get_Derivatives_of_StateAfterSpiral_wrt_Time();

            //store the Xindices of the right epoch
            this->Xindices_EventRightEpoch = this->mySpiral->get_Xindices_SpiralEndEpoch();
        }//end calcbounds_event_right_side

        void EphemerisPeggedSpiralArrival::process_event(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            this->process_event_left_side(X, Xindex, F, Findex, G, needG);

            this->process_event_main(X, Xindex, F, Findex, G, needG);

            this->process_event_right_side(X, Xindex, F, Findex, G, needG);

            this->process_specialized_constraints(X, Xindex, F, Findex, G, needG);
        }//end process_event()

        void EphemerisPeggedSpiralArrival::process_event_main(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            this->mySpiral->setInitialMass(this->state_before_event(6));
            this->mySpiral->setInitialEpoch(this->state_before_event(7));
            this->mySpiral->process_spiral(X, Xindex, F, Findex, G, needG);
        }//end process_event_main()

        void EphemerisPeggedSpiralArrival::process_event_right_side(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            this->state_after_event = this->mySpiral->getStateAfterSpiral();
            this->EventRightEpoch = this->state_after_event(7);

            if (needG)
            {
                this->Derivatives_of_StateAfterEvent = this->mySpiral->get_Derivatives_of_StateAfterSpiral();
                this->Derivatives_of_StateAfterEvent_wrt_Time = this->mySpiral->get_Derivatives_of_StateAfterSpiral_wrt_Time();
            }
        }//end process_event_right_side()

        void EphemerisPeggedSpiralArrival::output(std::ofstream& outputfile,
            const double& launchdate,
            size_t& eventcount)
        {
            this->mySpiral->output(outputfile, launchdate, eventcount);

            //output end of mission if appropriate
            this->EphemerisPeggedArrival::output(outputfile, launchdate, eventcount);
        }
    }//close namespace BoundaryEvents
}//close namespace EMTG