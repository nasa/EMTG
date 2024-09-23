
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

#include "EphemerisPeggedZeroTurnFlyby.h"
#include "EMTG_math.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        EphemerisPeggedZeroTurnFlyby::EphemerisPeggedZeroTurnFlyby(const std::string& name,
            const size_t& journeyIndex,
            const size_t& phaseIndex,
            size_t& stageIndex,
            Astrodynamics::universe* Universe,
            HardwareModels::Spacecraft* mySpacecraft,
            missionoptions* myOptions,
			EphemerisPeggedArrivalWithVinfinity* PreviousPhaseArrivalEvent) :
            EphemerisPeggedFlybyOut::EphemerisPeggedFlybyOut(name,
                journeyIndex,
                phaseIndex,
                stageIndex,
                Universe,
                mySpacecraft,
                myOptions,
                PreviousPhaseArrivalEvent)
        {
        }

        //******************************************calcbounds methods

        //calcbounds
        void EphemerisPeggedZeroTurnFlyby::calcbounds(std::vector<size_t> timeVariables)
        {
            this->calcbounds_event_left_side(timeVariables);

            //create bounds for the v-infinity
            double vinf_max = this->myUniverse->LU / this->myUniverse->TU;

            std::vector< std::tuple<double, double> > vinfBounds = std::vector< std::tuple<double, double> >(3, std::make_tuple(-vinf_max, vinf_max));

            this->calcbounds_event_main(vinfBounds);

            this->calcbounds_event_right_side();

            this->calcbounds_specialized_constraints();
        }//end calcbounds()

        void EphemerisPeggedZeroTurnFlyby
            ::calcbounds_event_main(const std::vector< std::tuple<double, double> >& vinfBounds)
        {
            //base class
            this->EphemerisPeggedFlybyOut::calcbounds_event_main(vinfBounds);

            //constraint: outgoing velocity V-infinity must match incoming V-infinity
            this->Gindices_Vinfinity_match.resize(3);
            for (size_t Vindex = 0; Vindex < 3; ++Vindex)
            {
                Flowerbounds->push_back(-math::SMALL);
                Fupperbounds->push_back(math::SMALL);
                Fdescriptions->push_back(prefix + "flyby v-infinity " + this->stateNames[3 + Vindex] + " match");
                
                this->create_sparsity_entry(Fdescriptions->size() - 1,
                    this->Xindices_Vinfinity_in[Vindex],
                    this->Gindices_Vinfinity_match[Vindex]);

                this->create_sparsity_entry(Fdescriptions->size() - 1,
                    this->Xindices_Vinfinity_out[Vindex],
                    this->Gindices_Vinfinity_match[Vindex]);
            }

        }//end calcbounds_event_main()

        //******************************************process methods

        void EphemerisPeggedZeroTurnFlyby::process_event(const std::vector<doubleType>& X,
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

            this->process_specialized_constraints(X, Xindex, F, Findex, G, needG);
        }//end process_event

        void EphemerisPeggedZeroTurnFlyby::process_event_main(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //Step 1: base class
            this->EphemerisPeggedFlybyOut::process_event_main(X, Xindex, F, Findex, G, needG);

            //Step 2: constraint: outgoing velocity V-infinity must match incoming V-infinity
            for (size_t Vindex = 0; Vindex < 3; ++Vindex)
            {
                F[Findex++] = ( this->Vinfinity_out(Vindex) - this->Vinfinity_in(Vindex) )
                    * this->myUniverse->continuity_constraint_scale_factors(3 + Vindex);
            }

            if (needG)
            {
                for (size_t Vindex = 0; Vindex < 3; ++Vindex)
                {
                    //with respect to v-infinity-in
                    size_t Gindex = this->Gindices_Vinfinity_match[Vindex][0];
                    size_t Xindex = this->jGvar->operator[](Gindex);
                    
                    G[Gindex] = this->X_scale_factors->operator[](Xindex) * -1.0 * this->myUniverse->continuity_constraint_scale_factors(3 + Vindex);

                    //with respect to v-infinity-out
                    Gindex = this->Gindices_Vinfinity_match[Vindex][1];
                    Xindex = this->jGvar->operator[](Gindex);

                    G[Gindex] = this->X_scale_factors->operator[](Xindex) * 1.0 * this->myUniverse->continuity_constraint_scale_factors(3 + Vindex);
                }
            }
        }//end process_event_main()

         //******************************************output methods
        void EphemerisPeggedZeroTurnFlyby::output(std::ofstream& outputfile,
            const double& launchdate,
            size_t& eventcount)
        {
            this->mySpacecraft->setActiveStage(this->stageIndex);

            std::string event_type = "zeroflyby";

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
                math::Matrix<doubleType> R_CB_Sun(3, 1, 0.0);                for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
                {
                    R_CB_Sun(stateIndex) = central_body_state_and_derivatives[stateIndex];
                }
                R_sc_Sun = this->state_after_event.getSubMatrix1D(0, 2) + R_CB_Sun;
            }            
            this->mySpacecraft->computePowerState(R_sc_Sun.getSubMatrix1D(0, 2).norm() / this->myOptions->AU, this->state_after_event(7));

            //compute RA and DEC for incoming asymptote in the body's local frame
            doubleType RA, DEC;
            //rotate the v-infinity vector to the local frame
            this->myUniverse->LocalFrame.construct_rotation_matrices(this->EventLeftEpoch / 86400.0 + 2400000.5, false);

            math::Matrix<doubleType> R_ICRF_to_local = this->myUniverse->LocalFrame.get_R(EMTG::ReferenceFrame::ICRF,
                this->myJourneyOptions->arrival_elements_frame);

            math::Matrix<doubleType> VinfinityLocalFrame = R_ICRF_to_local * this->Vinfinity_out;

            //compute RA and DEC
            RA = atan2(VinfinityLocalFrame(1), VinfinityLocalFrame(0));
            DEC = asin(VinfinityLocalFrame(2) / VinfinityLocalFrame.norm());

            //BdotT and BdotR
            Astrodynamics::bplane myBplane(this->myUniverse->central_body.mu);
            math::Matrix<doubleType> periapse_state = this->get_periapse_state().vert_cat(math::Matrix<doubleType>(1, 1, this->state_before_event(6)));
            this->myBplane.define_bplane_without_bend_angle(periapse_state, this->Vinfinity_in);
            this->BdotR = this->myBplane.getBdotR();
            this->BdotT = this->myBplane.getBdotT();

            write_output_line(outputfile,
                eventcount,
                event_type,
                boundary_name,
                this->EventTimeWidth,
                periapse_state.norm(),
                this->BdotR,
                this->BdotT,
                RA,
                DEC,
                this->Vinfinity_out.dot(this->Vinfinity_out),
                this->state_after_event,
                -this->Vinfinity_out,
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


        math::Matrix<doubleType> EphemerisPeggedZeroTurnFlyby::get_periapse_state()
        {
            //We will proceed using the same assumptions as in the Lucy mission - that a zero-turn flyby always occurs along the Sun-body line
            //Jacob thinks this is a good place to start for most missions, and hey, that's how Lucy works and they're paying for this.
            //therefore this code is derived from the Python code that we used in Lucy Phase B

            math::Matrix<doubleType> periapse_position(3, 1, 0.0);

            //get the vector to the Sun
            math::Matrix<doubleType> vector_to_sun(3, 1, 0.0);
            math::Matrix<doubleType> velocity_relative_to_sun(3, 1, 0.0);
            if (this->myUniverse->central_body.spice_ID == 10)
            {
                //the central body is the sun, so the vector to the sun is just the negative of the target body's position vector
                vector_to_sun = -this->state_before_event.getSubMatrix1D(0, 2);
                velocity_relative_to_sun = this->state_before_event.getSubMatrix1D(3, 5);
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
                {
                    vector_to_sun(stateIndex) = -(central_body_state[stateIndex] + this->state_before_event(stateIndex));
                    velocity_relative_to_sun(stateIndex) = central_body_state[stateIndex + 3] + this->state_before_event(stateIndex + 3);
                }
            }

            
            //mimic what we did in Phase B, for testing
            //compute the basis vectors of the target frame
            Astrodynamics::frame EclipticFrame;
            EclipticFrame.initialize_ecliptic();
            math::Matrix<doubleType> Vinf_ecliptic = EclipticFrame.get_R_from_ICRF_to_J2000BCI() * this->Vinfinity_in;
            math::Matrix<doubleType> RtoSun_ecliptic = EclipticFrame.get_R_from_ICRF_to_J2000BCI() * vector_to_sun;
            math::Matrix<doubleType> Xhat = Vinf_ecliptic.unitize();
            math::Matrix<doubleType> Yhat = RtoSun_ecliptic.unitcross(Xhat);
            math::Matrix<doubleType> Zhat = Xhat.unitcross(Yhat);

            math::Matrix<doubleType> R_from_Ecliptic_to_Target = Xhat.horz_cat(Yhat.horz_cat(Zhat));

            math::Matrix<doubleType> R_from_Target_to_Ecliptic = R_from_Ecliptic_to_Target.transpose();

            math::Matrix<doubleType> rp_ecliptic = R_from_Target_to_Ecliptic * math::Matrix<doubleType>(3, 1, std::vector<doubleType>({ 0.0, 0.0, this->myJourneyOptions->zero_turn_flyby_distance }));

            math::Matrix<doubleType> rp_ICRF = EclipticFrame.get_R_from_J2000BCI_to_ICRF() * rp_ecliptic;

            //********************************************now the real way

            //compute the target-centered angular momentum unit vector - the cross of the vector to the sun line with the v-infinity
            math::Matrix<doubleType> h_hat = (vector_to_sun).unitcross(this->Vinfinity_in);

            math::Matrix<doubleType> v_hat = this->Vinfinity_in.unitize();

            math::Matrix<doubleType> Spacecraft_Position_ICRF = v_hat.unitcross(h_hat) * this->myJourneyOptions->zero_turn_flyby_distance;

            return Spacecraft_Position_ICRF.vert_cat(this->Vinfinity_in);
        }//end get_periapse_state()

        //output
        void EphemerisPeggedZeroTurnFlyby::output_maneuver_and_target_spec(std::ofstream& maneuver_spec_file, std::ofstream& target_spec_file, bool& haveManeuverNeedTarget)
        {
            //We will proceed using the same assumptions as in the Lucy mission - that a zero-turn flyby always occurs along the Sun-body line
            //Jacob thinks this is a good place to start for most missions, and hey, that's how Lucy works and they're paying for this.

            if (haveManeuverNeedTarget)
            {
                //reset the flag that tells us we need a target spec line
                haveManeuverNeedTarget = false;

                //Step 1.1: initialize a target spec object
                math::Matrix<doubleType> periapse_state = this->get_periapse_state().vert_cat(math::Matrix<doubleType>(1, 1, this->state_before_event(6)));

                this->myBplane.define_bplane_without_bend_angle(periapse_state, this->Vinfinity_in);
                this->myBplane.define_bplane(periapse_state);
                target_spec_line myTargetSpecLine(this->name,
                    "EME2000",
                    this->myBody->name,
                    this->state_before_event(7),
                    periapse_state,
                    this->myBplane.getBdotR(),
                    this->myBplane.getBdotT());

                //Step 1.2: write target spec object
                myTargetSpecLine.write(target_spec_file);
            }
        }//end output_maneuver_and_target_spec()
    }//end namespace BoundaryEvents
}//end namespace EMTG