
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

#include "EphemerisPeggedIntercept.h"
#include "bplane.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        EphemerisPeggedIntercept::EphemerisPeggedIntercept(const std::string& name,
            const size_t& journeyIndex,
            const size_t& phaseIndex,
            size_t& stageIndex,
            Astrodynamics::universe* Universe,
            HardwareModels::Spacecraft* mySpacecraft,
            missionoptions* myOptions) :
            EphemerisPeggedFlybyIn::EphemerisPeggedFlybyIn(name,
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
        void EphemerisPeggedIntercept::calcbounds(std::vector<size_t> timeVariables)
        {

            this->calcbounds_event_left_side(timeVariables);

            this->calcbounds_event_main();

            this->calcbounds_event_right_side();

            if (this->hasTCM)
                this->calcbounds_virtual_propellant_constraints();

            this->calcbounds_specialized_constraints();
        }//end calcbounds()

        void EphemerisPeggedIntercept::calcbounds_event_main()
        {
            //create bounds for the v-infinity
            double vinf_max = this->myOptions->Journeys[this->journeyIndex].final_velocity[1];

            std::vector< std::tuple<double, double> > vinfBounds = std::vector< std::tuple<double, double> >(3, std::make_tuple(-vinf_max, vinf_max));

            //base class
            this->EphemerisPeggedArrivalWithVinfinity::calcbounds_event_main(vinfBounds);

            //v-infinity magnitude constraint
            Flowerbounds->push_back(this->myJourneyOptions->final_velocity[0] - this->myJourneyOptions->final_velocity[1]);
            Fupperbounds->push_back(0.0);
            Fdescriptions->push_back(prefix + "V_infinity_in magnitude");

            for (size_t Vindex = 0; Vindex < 3; ++Vindex)
            {
                this->create_sparsity_entry(Fdescriptions->size() - 1,
                    std::get<0>(this->Derivatives_of_StateBeforeEvent[dIndex_VbeforeEvent_dVinfinity_in[Vindex]]),
                    this->Gindices_incoming_velocity_magnitude_wrt_Vinfinity_in);
            }
        }//end calcbounds_event_main()

        //******************************************process methods

        void EphemerisPeggedIntercept::process_event(const std::vector<doubleType>& X,
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

        void EphemerisPeggedIntercept::process_event_main(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //call the base class
            this->EphemerisPeggedFlybyIn::process_event_main(X, Xindex, F, Findex, G, needG);

            //v-infinity magnitude constraint
            doubleType v_infinity_in_magnitude = this->Vinfinity_in.norm();
            F[Findex++] = v_infinity_in_magnitude - this->myJourneyOptions->final_velocity[1];

            if (needG)
            {
                for (size_t Vindex = 0; Vindex < 3; ++Vindex)
                {
                    size_t Gindex = this->Gindices_incoming_velocity_magnitude_wrt_Vinfinity_in[Vindex];
                    size_t Xindex = this->jGvar->operator[](Gindex);

                    G[Gindex] = this->X_scale_factors->operator[](Xindex)
                        * (this->Vinfinity_in(Vindex) / v_infinity_in_magnitude)_GETVALUE;
                }
            }
        }//end process_event_main()

         //******************************************output methods
        void EphemerisPeggedIntercept::output(std::ofstream& outputfile,
            const double& launchdate,
            size_t& eventcount)
        {
            this->mySpacecraft->setActiveStage(this->stageIndex);

            std::string event_type = "intercept";

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
                    this->mySpacecraft->getAvailablePower(),//AvailPower
                    0.0,//mdot
                    0,//number_of_active_engines
                    0.0,
                    "none");//active_power)
            }

            //compute RA and DEC for incoming asymptote in the body's local frame
            doubleType RA, DEC;
            //rotate the v-infinity vector to the local frame
            this->myUniverse->LocalFrame.construct_rotation_matrices(this->EventLeftEpoch / 86400.0 + 2400000.5, false);

            math::Matrix<doubleType> R_ICRF_to_local = this->myUniverse->LocalFrame.get_R(EMTG::ReferenceFrame::ICRF,
                this->myJourneyOptions->arrival_elements_frame);

            math::Matrix<doubleType> VinfinityLocalFrame = R_ICRF_to_local * this->Vinfinity_in;

            //compute RA and DEC
            RA = atan2(VinfinityLocalFrame(1), VinfinityLocalFrame(0));
            DEC = math::safe_asin(VinfinityLocalFrame(2) / (VinfinityLocalFrame.norm() + math::SMALL));

            //b-plane coordinates if this is the last event in the mission
            doubleType BdotR = 0.0;
            doubleType BdotT = 0.0;
            doubleType periapse_distance = 0.0;
            if (this->isLastEventInMission)
            {
                Astrodynamics::bplane myBplane(this->myUniverse->central_body.mu);
                math::Matrix<doubleType> periapse_state = this->get_periapse_state().vert_cat(math::Matrix<doubleType>(1, 1, this->state_before_event(6)));
                myBplane.define_bplane_without_bend_angle(periapse_state, this->Vinfinity_in);
                periapse_distance = periapse_state.norm();
                BdotR = myBplane.getBdotR();
                BdotT = myBplane.getBdotT();
            }

            this->write_output_line(outputfile,
                eventcount,
                event_type,
                boundary_name,
                this->EventTimeWidth,
                0.0,
                BdotT,
                BdotR,
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

            //output end of mission if appropriate
            this->EphemerisPeggedArrival::output(outputfile, launchdate, eventcount);
        }//end output()

        math::Matrix<doubleType> EphemerisPeggedIntercept::get_periapse_state()
        {
            //This code will only be used for the last intercept in a mission, and assumes that said intercept is a zero-turn flyby.
            //We will proceed using the same assumptions as in the Lucy mission - that a zero-turn flyby always occurs along the Sun-body line
            //Jacob thinks this is a good place to start for most missions, and hey, that's how Lucy works and they're paying for this.
            //therefore this code is derived from the Python code that we used in Lucy Phase B

            //get the vector to the Sun
            math::Matrix<doubleType> vector_to_sun(3, 1, 0.0);
            if (this->myUniverse->central_body.spice_ID == 10)
            {
                //the central body is the sun, so the vector to the sun is just the negative of the target body's position vector
                vector_to_sun = -this->state_before_event.getSubMatrix1D(0, 2);
            }
            else
            {
                //the central body is not the sun, so we need to locate the central body relative to the sun and then add it to the target body's position vector
                //and then negate the whole thing
                doubleType central_body_state[12];
                this->myUniverse->locate_central_body(this->state_before_event(7),
                    central_body_state,
                    *this->myOptions,
                    false);

                for (size_t stateIndex : {0, 1, 2})
                    vector_to_sun(stateIndex) = -(central_body_state[stateIndex] + this->state_before_event(stateIndex));
            }

            math::Matrix<doubleType> h_hat = (vector_to_sun).unitcross(this->Vinfinity_in);

            math::Matrix<doubleType> v_hat = this->Vinfinity_in.unitize();

            math::Matrix<doubleType> Spacecraft_Position_ICRF = v_hat.unitcross(h_hat) * this->myJourneyOptions->zero_turn_flyby_distance;

            return Spacecraft_Position_ICRF.vert_cat(this->Vinfinity_in);
        }//end get_periapse_state()

        void EphemerisPeggedIntercept::output_maneuver_and_target_spec(std::ofstream& maneuver_spec_file, std::ofstream& target_spec_file, bool& haveManeuverNeedTarget)
        {
            //IFF this is the last event in the mission, write a target spec
            //assume that we are doing a Lucy-style flyby of the sub-solar point

            if (this->isLastEventInMission
                && haveManeuverNeedTarget)
            {
                //reset the flag that tells us we need a target spec line
                haveManeuverNeedTarget = false;

                //Step 1.1: initialize a target spec object
                Astrodynamics::bplane myBplane(this->myUniverse->central_body.mu);
                math::Matrix<doubleType> periapse_state = this->get_periapse_state().vert_cat(math::Matrix<doubleType>(1, 1, this->state_before_event(6)));
                myBplane.define_bplane_without_bend_angle(periapse_state, this->Vinfinity_in);

                target_spec_line myTargetSpecLine(this->name,
                    "EME2000",
                    this->myBody->name,
                    this->state_before_event(7),
                    periapse_state,
                    myBplane.getBdotR(),
                    myBplane.getBdotT());

                //Step 1.2: write target spec object
                myTargetSpecLine.write(target_spec_file);
            }
        }//end output_maneuver_and_target_spec()
    }//end namespace BoundaryEvents
}//end namespace EMTG