
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

#include "EphemerisPeggedFlybyOut.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        EphemerisPeggedFlybyOut::EphemerisPeggedFlybyOut(const std::string& name,
            const size_t& journeyIndex,
            const size_t& phaseIndex,
            size_t& stageIndex,
            Astrodynamics::universe* Universe,
            HardwareModels::Spacecraft* mySpacecraft,
            missionoptions* myOptions,
			EphemerisPeggedArrivalWithVinfinity* PreviousPhaseArrivalEvent) :
            EphemerisPeggedDeparture::EphemerisPeggedDeparture(name,
                journeyIndex,
                phaseIndex,
                stageIndex,
                Universe,
                mySpacecraft,
                myOptions,
                PreviousPhaseArrivalEvent)
        {
            this->Vinfinity_out.resize(3, 1, 0.0);
            this->Vinfinity_in.resize(3, 1, 0.0);

			this->PreviousPhaseArrivalEvent = PreviousPhaseArrivalEvent;

            this->myBplane = Astrodynamics::bplane(this->myBody->mu);
        }

        //******************************************calcbounds methods
        void EphemerisPeggedFlybyOut::
            calcbounds_event_left_side(std::vector<size_t> timeVariables)
        {
            //Step 1: calcbounds on the left epoch
            this->calculate_dependencies_left_epoch(timeVariables);

            //Step 2: track the Xindex of the first variable in the event
            this->X_index_of_first_decision_variable_in_this_event = this->Xdescriptions->size();

            //Step 3: get derivative skeleton from the previous event
            std::vector< std::tuple<size_t, size_t, double> >& PreviousEvent_Derivatives_of_StateAfterEvent
                = this->PreviousPhaseArrivalEvent->get_Derivatives_of_StateAfterEvent();

            //Step 4: find the Xindex of the previous event's mass
            size_t Xindex_previous_event_end_mass = 0;
            for (int Xindex = this->X_index_of_first_decision_variable_in_this_event - 1; Xindex >= 0; --Xindex)
            {
                if (this->Xdescriptions->at(Xindex).find("mass") < 1024
                    && this->Xdescriptions->at(Xindex).find("event") < 1024)
                {
                    Xindex_previous_event_end_mass = Xindex;
                    break;
                }
            }

            //Step 5: find the dIndex of the previous event's state with respect to its mass
            for (size_t dIndex = 0; dIndex < PreviousEvent_Derivatives_of_StateAfterEvent.size(); ++dIndex)
            {
                size_t Xindex = std::get<0>(PreviousEvent_Derivatives_of_StateAfterEvent[dIndex]);
                if (Xindex == Xindex_previous_event_end_mass)
                    this->dIndex_PreviousEventMass_PreviousEventMass = dIndex;
            }

            //Step 6: construct this phase's derivative skeleton
            this->Derivatives_of_StateBeforeEvent_wrt_Time = this->PreviousPhaseArrivalEvent->get_Derivatives_of_StateAfterEvent_wrt_Time();

            this->Derivatives_of_StateBeforeEvent.push_back(PreviousEvent_Derivatives_of_StateAfterEvent[this->dIndex_PreviousEventMass_PreviousEventMass]);
            this->dIndex_CurrentEventMass_PreviousEventMass = this->Derivatives_of_StateBeforeEvent.size() - 1;
        }//end calcbounds_event_left_side()

        void EphemerisPeggedFlybyOut::
            calcbounds_event_main(const std::vector< std::tuple<double, double> >& vinfBounds)
        {
            //Step 1: encode V-infinity
            for (size_t Vindex = 0; Vindex < 3; ++Vindex)
            {
                Xlowerbounds->push_back(std::get<0>(vinfBounds[Vindex]));
                Xupperbounds->push_back(std::get<1>(vinfBounds[Vindex]));
                Xdescriptions->push_back(prefix + "V_infinity_" + this->stateNames[Vindex]);
                X_scale_factors->push_back(this->myUniverse->LU / this->myUniverse->TU);
                this->Derivatives_of_StateBeforeEvent.push_back({ Xdescriptions->size() - 1, Vindex + 3, 1.0 });
                this->dIndex_VbeforeEvent_dVinfinity_out.push_back(this->Derivatives_of_StateBeforeEvent.size() - 1);
                this->Xindices_Vinfinity_out.push_back(Xdescriptions->size() - 1);
            }

            //Step 2: locate the previous phase's V_infinity_in components
            for (size_t Vindex = 0; Vindex < 3; ++Vindex)
            {
                std::string variableName = "V_infinity_" + this->stateNames[Vindex];

                for (int Xindex = this->PreviousPhaseArrivalEvent->getX_index_of_first_decision_variable_in_this_event(); Xindex < this->X_index_of_first_decision_variable_in_this_event; ++Xindex)
                {
                    if (Xdescriptions->at(Xindex).find(variableName) < 1024)
                    {
                        this->Xindices_Vinfinity_in.push_back(Xindex);
                        break;
                    }
                }
            }
        }//end calcbounds_event_main()

         //******************************************process methods
        void EphemerisPeggedFlybyOut::process_event_left_side(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //Step 1: copy the state from the previous event
            this->state_before_event = this->PreviousPhaseArrivalEvent->get_state_after_event_raw();
            this->state_before_event(6) = this->PreviousPhaseArrivalEvent->get_state_after_event()(6);

            this->boundary_state = this->state_before_event;

            if (needG)
            {
                //extract the derivatives from the previous event and put them where we need them
                //Step 2: get derivatives from the previous event
                std::vector< std::tuple<size_t, size_t, double> >& PreviousEvent_Derivatives_of_StateBeforeEvent
                    = this->PreviousPhaseArrivalEvent->get_Derivatives_of_StateAfterEvent();

                //Step 3: put them where we need them
                //Step 3.1: wrt epoch
                this->Derivatives_of_StateBeforeEvent_wrt_Time = this->PreviousPhaseArrivalEvent->get_Derivatives_of_StateAfterEvent_wrt_Time();

                //Step 3.2: wrt mass
                this->Derivatives_of_StateBeforeEvent[this->dIndex_CurrentEventMass_PreviousEventMass]
                    = PreviousEvent_Derivatives_of_StateBeforeEvent[this->dIndex_PreviousEventMass_PreviousEventMass];
            }
        }

        void EphemerisPeggedFlybyOut::process_event_main(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //Step 1: extract v-infinity and add to state_before_event
            for (size_t Vindex = 0; Vindex < 3; ++Vindex)
            {
                this->Vinfinity_out(Vindex) = X[Xindex++];
                this->state_before_event(3 + Vindex) += this->Vinfinity_out(Vindex);
            }

            this->C3 = this->Vinfinity_out.dot(this->Vinfinity_out);

            //Step 2: extract incoming v-infinity
            for (size_t Vindex = 0; Vindex < 3; ++Vindex)
            {
                this->Vinfinity_in(Vindex) = X[this->Xindices_Vinfinity_in[Vindex]];
            }
        }//end process_event_left_side()

        math::Matrix<doubleType> EphemerisPeggedFlybyOut::calculate_flyby_periapse_state()
        {
            //calculate unit vector pointed towards periapse
            math::Matrix<doubleType> periapse_position_unit_vector = (this->Vinfinity_in.unitize() - this->Vinfinity_out.unitize()).unitize();

            //calculate angular momentum unit vector
            math::Matrix<doubleType> angular_momentum_vector = this->Vinfinity_in.cross(this->Vinfinity_out);

            //calculate velocity unit vector at periapse
            math::Matrix<doubleType> periapse_velocity_unit_vector = (angular_momentum_vector.cross(periapse_position_unit_vector)).unitize();

            //calculate velocity magnitude at periapse
            doubleType periapse_velocity_magnitude = sqrt(2 * this->myBody->mu / (this->myBody->radius + this->FlybyAltitude) + this->Vinfinity_in.dot(this->Vinfinity_out));

            //transform from unit vector space to state space
            math::Matrix<doubleType> periapse_position_vector = periapse_position_unit_vector * (this->myBody->radius + this->FlybyAltitude);
            math::Matrix<doubleType> periapse_velocity_vector = periapse_velocity_unit_vector * periapse_velocity_magnitude;

            return periapse_position_vector.vert_cat(periapse_velocity_vector);
        }//end calculate_flyby_periapse_state()


        void EphemerisPeggedFlybyOut::output_periapse_state(size_t& flybyIndex, std::ofstream& outputfile)
        {
            math::Matrix<doubleType> periapse_state = this->get_periapse_state();
            outputfile << "Flyby: " << flybyIndex++ << " (" << this->myBody->name << ")";
            for (int k = 0; k < 6; ++k)
                outputfile << " " << periapse_state(k) _GETVALUE;
            outputfile << std::endl;
        }
    }//end namespace BoundaryEvents
}//end namespace EMTG