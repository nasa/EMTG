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

#include "EphemerisPeggedUnpoweredFlyby.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        EphemerisPeggedUnpoweredFlyby::EphemerisPeggedUnpoweredFlyby(const std::string& name,
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
        void EphemerisPeggedUnpoweredFlyby::calcbounds(std::vector<size_t> timeVariables)
        {
            this->calcbounds_event_left_side(timeVariables);

            //create bounds for the v-infinity
            double vinf_max = fmax(this->myUniverse->LU / this->myUniverse->TU, this->PreviousPhaseArrivalEvent->get_Vinfinity_upperbound());

            std::vector< std::tuple<double, double> > vinfBounds = std::vector< std::tuple<double, double> >(3, std::make_tuple(-vinf_max, vinf_max));
            
            this->calcbounds_event_main(vinfBounds);

            this->calcbounds_event_right_side();

            this->calcbounds_specialized_constraints();
        }//end calcbounds()

        void EphemerisPeggedUnpoweredFlyby::
            calcbounds_event_main(const std::vector< std::tuple<double, double> >& vinfBounds)
        {
            //Step 1: base class
            this->EphemerisPeggedFlybyOut::calcbounds_event_main(vinfBounds);

            //Step 2: constraint: outgoing velocity V-infinity magnitude must match incoming V-infinity magnitude
            this->Gindices_Vinfinity_magnitude_match.resize(3);

            Flowerbounds->push_back(-math::SMALL);
            Fupperbounds->push_back(math::SMALL);
            Fdescriptions->push_back(prefix + "flyby v-infinity magnitude match");

            for (size_t Vindex = 0; Vindex < 3; ++Vindex)
            {
                this->create_sparsity_entry(Fdescriptions->size() - 1,
                    this->Xindices_Vinfinity_in[Vindex],
                    this->Gindices_Vinfinity_magnitude_match[Vindex]);

                this->create_sparsity_entry(Fdescriptions->size() - 1,
                    this->Xindices_Vinfinity_out[Vindex],
                    this->Gindices_Vinfinity_magnitude_match[Vindex]);
            }

            //Step 3: constraint: flyby altitude must be safe
            if (this->myJourneyOptions->override_flyby_altitude_bounds)
            {
                //zero is the Universe file's minimum flyby altitude
                Flowerbounds->push_back((this->myBody->minimum_safe_flyby_altitude - this->myJourneyOptions->flyby_altitude_bounds[1]) / (this->myBody->minimum_safe_flyby_altitude + this->myBody->radius));
                Fupperbounds->push_back((this->myBody->minimum_safe_flyby_altitude - this->myJourneyOptions->flyby_altitude_bounds[0]) / (this->myBody->minimum_safe_flyby_altitude + this->myBody->radius));
            }
            else
            {
                if (this->myBody->mu / this->myUniverse->mu < 1.0e+5)
                    Flowerbounds->push_back(-20.0);
                else
                    Flowerbounds->push_back(-300.0);
                Fupperbounds->push_back(0.0);
            }
            Fdescriptions->push_back(prefix + "flyby altitude constraint (above minimum altitude but below [10x/300x] altitude for [rocky/gas] planets");

            this->Gindices_turn_angle.resize(3);

            for (size_t Vindex = 0; Vindex < 3; ++Vindex)
            {
                this->create_sparsity_entry(Fdescriptions->size() - 1,
                    this->Xindices_Vinfinity_in[Vindex],
                    this->Gindices_turn_angle[Vindex]);

                this->create_sparsity_entry(Fdescriptions->size() - 1,
                    this->Xindices_Vinfinity_out[Vindex],
                    this->Gindices_turn_angle[Vindex]);
            }

        }//end calcbounds_event_main()
        //******************************************process methods

        void EphemerisPeggedUnpoweredFlyby::process_event(const std::vector<doubleType>& X,
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

        void EphemerisPeggedUnpoweredFlyby::process_event_main(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //Step 1: base class
            this->EphemerisPeggedFlybyOut::process_event_main(X, Xindex, F, Findex, G, needG);

            //Step 2: constraint: outgoing velocity V-infinity magnitude must match incoming V-infinity magnitude
            doubleType VinfinityInMagnitude = this->Vinfinity_in.norm() + 1.0e-25;
            doubleType VinfinityOutMagnitude = this->Vinfinity_out.norm() + 1.0e-25;

            F[Findex++] = (VinfinityOutMagnitude - VinfinityInMagnitude)
                / this->myUniverse->LU * this->myUniverse->TU;

            if (needG)
            {
                for (size_t Vindex = 0; Vindex < 3; ++Vindex)
                {
                    //with respect to v-infinity-in
                    size_t Gindex = this->Gindices_Vinfinity_magnitude_match[Vindex][0];
                    size_t Xindex = this->jGvar->operator[](Gindex);
                    
                    G[Gindex] = this->X_scale_factors->operator[](Xindex)
                        * (-this->Vinfinity_in(Vindex) / VinfinityInMagnitude) _GETVALUE
                        / this->myUniverse->LU * this->myUniverse->TU;

                    //with respect to v-infinity-out
                    Gindex = this->Gindices_Vinfinity_magnitude_match[Vindex][1];
                    Xindex = this->jGvar->operator[](Gindex);

                    G[Gindex] = this->X_scale_factors->operator[](Xindex)
                        * (this->Vinfinity_out(Vindex) / VinfinityOutMagnitude) _GETVALUE
                        / this->myUniverse->LU * this->myUniverse->TU;
                }
            }

            //Step 3: constraint: flyby altitude must be safe
            //Step 3.1: compute the flyby turn angle
            doubleType cos_turn_angle = this->Vinfinity_in.dot(this->Vinfinity_out) / (VinfinityInMagnitude * VinfinityOutMagnitude);
            this->FlybyTurnAngle = math::safe_acos(cos_turn_angle);

            //Step 3.2: compute the required periapse altitude
            doubleType C3 = VinfinityOutMagnitude * VinfinityOutMagnitude;
            this->FlybyAltitude = this->FlybyTurnAngle == 0.0 ? 0.0 : this->myBody->mu / C3 * (1 / sin(this->FlybyTurnAngle / 2.0) - 1) - this->myBody->radius;

            //Step 3.3: constraint
            F[Findex++] = (this->myBody->minimum_safe_flyby_altitude - this->FlybyAltitude) / (this->myBody->minimum_safe_flyby_altitude + this->myBody->radius);

            //Step 3.4: derivatives
            if (needG)
            {
                double Vfx = this->Vinfinity_out(0) _GETVALUE;
                double Vfy = this->Vinfinity_out(1) _GETVALUE;
                double Vfz = this->Vinfinity_out(2) _GETVALUE;
                double V0x = this->Vinfinity_in(0) _GETVALUE;
                double V0y = this->Vinfinity_in(1) _GETVALUE;
                double V0z = this->Vinfinity_in(2) _GETVALUE;

                double VfdotVf = this->Vinfinity_out.dot(Vinfinity_out) _GETVALUE;
                double VfdotV0 = this->Vinfinity_in.dot(Vinfinity_out) _GETVALUE;
                double V0dotV0 = this->Vinfinity_in.dot(Vinfinity_in) _GETVALUE;

                double term1 = VfdotV0 / sqrt(VfdotVf * V0dotV0);
                double acos_term1_over2 = 0.5 * math::safe_acos(term1)_GETVALUE;

                double mu = this->myBody->mu;

                double termAcoeff = 2.0 * mu * (1 / sin(acos_term1_over2) - 1) / (VfdotVf*VfdotVf);
                double termBcoeff = mu * cos(acos_term1_over2) / ((term1 - 1) * sqrt(1 - VfdotV0*VfdotV0 / (V0dotV0*VfdotVf)) * sqrt(V0dotV0 * VfdotVf*VfdotVf*VfdotVf*VfdotVf*VfdotVf));
                double termA, termB;
                double Gentry;
                size_t Gindex, Xindex;

                //Vfx
                termA = termAcoeff * Vfx;
                termB = termBcoeff * (V0x*Vfy*Vfy - V0y*Vfx*Vfy + V0x*Vfz*Vfz - V0z*Vfx*Vfz);

                Gindex = Gindices_turn_angle[0][1];
                Xindex = this->jGvar->operator[](Gindex);

                Gentry = termA + termB;

                G[Gindex] = this->X_scale_factors->operator[](Xindex) * Gentry / (this->myBody->minimum_safe_flyby_altitude + this->myBody->radius);

                //Vfy
                termA = termAcoeff * Vfy;
                termB = termBcoeff * (V0y*Vfx*Vfx - V0x*Vfy*Vfx + V0y*Vfz*Vfz - V0z*Vfy*Vfz);

                Gindex = Gindices_turn_angle[1][1];
                Xindex = this->jGvar->operator[](Gindex);

                Gentry = termA + termB;

                G[Gindex] = this->X_scale_factors->operator[](Xindex) * Gentry / (this->myBody->minimum_safe_flyby_altitude + this->myBody->radius);

                //Vfz
                termA = termAcoeff * Vfz;
                termB = termBcoeff * (V0z*Vfx*Vfx - V0x*Vfz*Vfx + V0z*Vfy*Vfy - V0y*Vfz*Vfy);

                Gindex = Gindices_turn_angle[2][1];
                Xindex = this->jGvar->operator[](Gindex);

                Gentry = termA + termB;

                G[Gindex] = this->X_scale_factors->operator[](Xindex) * Gentry / (this->myBody->minimum_safe_flyby_altitude + this->myBody->radius);

                //derivatives with respect to V0
                double Coeff = mu * cos(acos_term1_over2) / ((term1 - 1) * sqrt((1 - term1*term1) * V0dotV0*V0dotV0*V0dotV0 * VfdotVf*VfdotVf*VfdotVf));

                //V0x
                Gindex = Gindices_turn_angle[0][0];
                Xindex = this->jGvar->operator[](Gindex);

                Gentry = Coeff * (Vfx*V0y*V0y - V0x*Vfy*V0y + Vfx*V0z*V0z - V0x*Vfz*V0z);

                G[Gindex] = this->X_scale_factors->operator[](Xindex) * Gentry / (this->myBody->minimum_safe_flyby_altitude + this->myBody->radius);

                //V0y
                Gindex = Gindices_turn_angle[1][0];
                Xindex = this->jGvar->operator[](Gindex);

                Gentry = Coeff * (Vfy*V0x*V0x - V0y*Vfx*V0x + Vfy*V0z*V0z - V0y*Vfz*V0z);

                G[Gindex] = this->X_scale_factors->operator[](Xindex) * Gentry / (this->myBody->minimum_safe_flyby_altitude + this->myBody->radius);

                //V0z
                Gindex = Gindices_turn_angle[2][0];
                Xindex = this->jGvar->operator[](Gindex);

                Gentry = Coeff * (Vfz*V0x*V0x - V0z*Vfx*V0x + Vfz*V0y*V0y - V0z*Vfy*V0y);

                G[Gindex] = this->X_scale_factors->operator[](Xindex) * Gentry / (this->myBody->minimum_safe_flyby_altitude + this->myBody->radius);

            }
        }//end process_event_main()

         //******************************************output methods
        void EphemerisPeggedUnpoweredFlyby::output(std::ofstream& outputfile,
            const double& launchdate,
            size_t& eventcount)
        {
            this->mySpacecraft->setActiveStage(this->stageIndex);

            std::string event_type = "upwr_flyby";

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

            //rotate the v-infinity vector to the local frame
            this->myUniverse->LocalFrame.construct_rotation_matrices(this->EventLeftEpoch / 86400.0 + 2400000.5, false);

            math::Matrix<doubleType> R_ICRF_to_local = this->myUniverse->LocalFrame.get_R(EMTG::ReferenceFrame::ICRF,
                this->myJourneyOptions->arrival_elements_frame);

            math::Matrix<doubleType> VinfinityLocalFrame = R_ICRF_to_local * this->Vinfinity_out;

            //compute RA and DEC
            doubleType RA = atan2(VinfinityLocalFrame(1), VinfinityLocalFrame(0));
            doubleType DEC = asin(VinfinityLocalFrame(2) / VinfinityLocalFrame.norm());

            //compute BdotR and BdotT
            math::Matrix<doubleType> periapse_state = this->calculate_flyby_periapse_state();
            this->myBplane.define_bplane(this->get_periapse_state());
            this->BdotR = this->myBplane.getBdotR();
            this->BdotT = this->myBplane.getBdotT();

            write_output_line(outputfile,
                eventcount,
                event_type,
                boundary_name,
                this->EventTimeWidth,
                this->FlybyAltitude,
                this->BdotR,
                this->BdotT,
                RA,
                DEC,
                this->Vinfinity_out.dot(this->Vinfinity_out),
                this->state_after_event,
                this->Vinfinity_out,
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

        math::Matrix<doubleType> EphemerisPeggedUnpoweredFlyby::get_periapse_state()
        {
            //calculate unit vector pointed towards periapse
            math::Matrix<doubleType> periapse_position_unit_vector = (this->Vinfinity_in.unitize() - this->Vinfinity_out.unitize()).unitize();

            //calculate angular momentum unit vector
            math::Matrix<doubleType> angular_momentum_vector = this->Vinfinity_in.cross(this->Vinfinity_out);

            //calculate velocity unit vector at periapse
            math::Matrix<doubleType> periapse_velocity_unit_vector = (angular_momentum_vector.cross(periapse_position_unit_vector)).unitize();

            //calculate velocity magnitude at periapse
            doubleType periapse_velocity_magnitude = sqrt(2 * this->myBody->mu / (this->myBody->radius + this->FlybyAltitude) + this->Vinfinity_in.dot(this->Vinfinity_in));

            //transform from unit vector space to state space
            math::Matrix<doubleType> periapse_position_vector = periapse_position_unit_vector * (this->myBody->radius + this->FlybyAltitude);
            math::Matrix<doubleType> periapse_velocity_vector = periapse_velocity_unit_vector * periapse_velocity_magnitude;

            math::Matrix<doubleType> periapse_state = periapse_position_vector.vert_cat(periapse_velocity_vector);

            return periapse_state;
        }//end get_periapse_state()

        void EphemerisPeggedUnpoweredFlyby::output_maneuver_and_target_spec(std::ofstream& maneuver_spec_file, std::ofstream& target_spec_file, bool& haveManeuverNeedTarget)
        {
            //in this case, we don't actually know the periapse state, so we'll use the approximated one that we already output in the .emtg file
            
            if (haveManeuverNeedTarget)
            {
                //reset the flag that tells us we need a target spec line
                haveManeuverNeedTarget = false;

                //Step 1: initialize a target spec object
                this->myBplane.define_bplane(this->get_periapse_state());

                target_spec_line myTargetSpecLine(this->name,
                    "EME2000",
                    this->myUniverse->central_body.name,
                    this->state_after_event(7),
                    this->state_after_event,
                    this->myBplane.getBdotR(),
                    this->myBplane.getBdotT());

                //Step 2: write target spec object
                myTargetSpecLine.write(target_spec_file);
            }
        }//end output_maneuver_and_target_spec()
    }//end namespace BoundaryEvents
}//end namespace EMTG