// EMTG: Evolutionary Mission Trajectory Generator
// An open-source global optimization tool for preliminary mission design
// Provided by NASA Goddard Space Flight Center
//
// Copyright (c) 2014 - 2017 United States Government as represented by the
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
        EphemerisPeggedUnpoweredFlyby::EphemerisPeggedUnpoweredFlyby(const std::string name,
            const size_t& journeyIndex,
            const size_t& phaseIndex,
            size_t& stageIndex,
            Astrodynamics::universe* Universe,
            HardwareModels::Spacecraft* mySpacecraft,
            missionoptions* myOptions,
            ArrivalEvent* PreviousPhaseArrivalEvent) :
            EphemerisPeggedFlybyOut::EphemerisPeggedFlybyOut(name,
                journeyIndex,
                phaseIndex,
                stageIndex,
                Universe,
                mySpacecraft,
                myOptions,
                PreviousPhaseArrivalEvent)
        {
            this->Gindices_TurnAngleConstraint_with_respect_to_Vinfinity.resize(3);
        }

        //******************************************calcbounds methods

        //calcbounds
        void EphemerisPeggedUnpoweredFlyby::calcbounds()
        {
            this->calcbounds_event_left_side();

            //create bounds for the v-infinity
            double vinf_max = this->myUniverse->LU / this->myUniverse->TU;

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
            
            //Step 3: encode periapse distance
            if (this->myOptions->Journeys[this->journeyIndex].journey_override_flyby_altitude_bounds)
            {
                Xlowerbounds->push_back(this->myBody->radius + this->myOptions->Journeys[this->journeyIndex].journey_flyby_altitude_bounds[0]);
                Xupperbounds->push_back(this->myBody->radius + this->myOptions->Journeys[this->journeyIndex].journey_flyby_altitude_bounds[1]);
            }
            else
            {
                Xlowerbounds->push_back(this->myBody->radius + this->myBody->minimum_safe_flyby_altitude);
                if (this->myBody->mass < 1.0e+25)
                    Xupperbounds->push_back(10.0 * this->myBody->radius);
                else
                    Xupperbounds->push_back(300.0 * this->myBody->radius);
            }
            X_scale_factors->push_back(this->myBody->radius);
            Xdescriptions->push_back(prefix + "flyby periapse radius");

            //Step 4: turn angle feasibility constraint
            Flowerbounds->push_back(-math::SMALL);
            Fupperbounds->push_back(math::SMALL);
            Fdescriptions->push_back(prefix + "Unpowered flyby turn angle constraint");

            //derivatives with respect to periapse distance, v-infinity in, and v-infinity out
            //rp
            this->create_sparsity_entry(Fdescriptions->size() - 1,
                Xdescriptions->size() - 1,
                Gindices_TurnAngleConstraint_with_respect_to_FlybyPeriapseDistance);

            //V-infinity
            for (size_t Vindex = 0; Vindex < 3; ++Vindex)
            {
                this->create_sparsity_entry(Fdescriptions->size() - 1,
                    this->Xindices_Vinfinity_in[Vindex],
                    this->Gindices_TurnAngleConstraint_with_respect_to_Vinfinity[Vindex]);

                this->create_sparsity_entry(Fdescriptions->size() - 1,
                    this->Xindices_Vinfinity_out[Vindex],
                    this->Gindices_TurnAngleConstraint_with_respect_to_Vinfinity[Vindex]);
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
                    size_t Xindex = this->myOptions->jGvar[Gindex];
                    
                    G[Gindex] = this->myOptions->X_scale_factors[Xindex]
                        * (-this->Vinfinity_in(Vindex) / VinfinityInMagnitude) _GETVALUE
                        / this->myUniverse->LU * this->myUniverse->TU;

                    //with respect to v-infinity-out
                    Gindex = this->Gindices_Vinfinity_magnitude_match[Vindex][1];
                    Xindex = this->myOptions->jGvar[Gindex];

                    G[Gindex] = this->myOptions->X_scale_factors[Xindex]
                        * (this->Vinfinity_out(Vindex) / VinfinityOutMagnitude) _GETVALUE
                        / this->myUniverse->LU * this->myUniverse->TU;
                }
            }

            //Step 3: extract flyby periapse distance and compute altitude
            this->FlybyPeriapseDistance = X[Xindex++];
            this->FlybyAltitude = this->FlybyPeriapseDistance - this->myBody->radius;

            //Step 3: compute turn angle constraint

            doubleType V0dotV0 = Vinfinity_in.dot(Vinfinity_in);
            doubleType VfdotVf = Vinfinity_out.dot(Vinfinity_out);
            doubleType vinf_in = sqrt(V0dotV0);
            doubleType vinf_out = sqrt(VfdotVf);
            double mu = this->myBody->mu;

            //compute turn angle
            doubleType denom = vinf_in * vinf_out;
            this->FlybyTurnAngle = math::safe_acos(this->Vinfinity_in.dot(this->Vinfinity_out) / ((denom > math::SMALL) ? denom : math::SMALL));

            //compute eccentricity of incoming and outgoing hyperbolas
            doubleType e_in = 1 + V0dotV0 * this->FlybyPeriapseDistance / mu;
            e_in = (e_in > math::SMALL) ? e_in : math::SMALL;

            doubleType e_out = 1 + VfdotVf * this->FlybyPeriapseDistance / mu;
            e_out = (e_out > math::SMALL) ? e_out : math::SMALL;

            //apply the flyby turn angle feasibility constraint
            F[Findex++] = math::safe_asin(1 / e_in) + math::safe_asin(1 / e_out) - this->FlybyTurnAngle;

            //derivatives of powered flyby turn angle constraint
            //(thanks SymPy!)
            if (needG)
            {
                size_t Gindex, Xindex;

                double v0x = this->Vinfinity_in(0) _GETVALUE;
                double v0y = this->Vinfinity_in(1) _GETVALUE;
                double v0z = this->Vinfinity_in(2) _GETVALUE;
                double vfx = this->Vinfinity_out(0) _GETVALUE;
                double vfy = this->Vinfinity_out(1) _GETVALUE;
                double vfz = this->Vinfinity_out(2) _GETVALUE;
                double V0dotVf = Vinfinity_in.dot(Vinfinity_out) _GETVALUE;
                double e_in2 = (e_in * e_in) _GETVALUE;
                double e_out2 = (e_out * e_out) _GETVALUE;
                double rp = this->FlybyPeriapseDistance _GETVALUE;

                //with respect to rp
                Gindex = Gindices_TurnAngleConstraint_with_respect_to_FlybyPeriapseDistance;
                Xindex = this->myOptions->jGvar[Gindex];
                G[Gindex] = this->myOptions->X_scale_factors[Xindex]
                    * (-(VfdotVf) / (mu*e_out2*sqrt(1 - 1 / e_out2)) - (V0dotV0) / (mu*e_in2*sqrt(1 - 1 / e_in2))) _GETVALUE;

                //with respect to v_infinity_in and v_infinity_out
                double dF_dvfx = ((-v0x*(V0dotVf) / (sqrt(V0dotV0*VfdotVf)*V0dotV0) + vfx / sqrt(V0dotV0*VfdotVf)) / sqrt(1 - V0dotVf*V0dotVf / (V0dotV0*VfdotVf)) - 2 * rp*v0x / (mu*e_in2*sqrt(1 - 1 / e_in2)))_GETVALUE;
                double dF_dvfy = ((-v0y*(V0dotVf) / (sqrt(V0dotV0*VfdotVf)*V0dotV0) + vfy / sqrt(V0dotV0*VfdotVf)) / sqrt(1 - V0dotVf*V0dotVf / (V0dotV0*VfdotVf)) - 2 * rp*v0y / (mu*e_in2*sqrt(1 - 1 / e_in2)))_GETVALUE;
                double dF_dvfz = ((-v0z*(V0dotVf) / (sqrt(V0dotV0*VfdotVf)*V0dotV0) + vfz / sqrt(V0dotV0*VfdotVf)) / sqrt(1 - V0dotVf*V0dotVf / (V0dotV0*VfdotVf)) - 2 * rp*v0z / (mu*e_in2*sqrt(1 - 1 / e_in2)))_GETVALUE;


                double dF_dv0x = ((v0x / sqrt(V0dotV0*VfdotVf) - vfx*(V0dotVf) / (sqrt(V0dotV0*VfdotVf)*VfdotVf)) / sqrt(1 - V0dotVf*V0dotVf / (V0dotV0*VfdotVf)) - 2 * rp*vfx / (mu*e_out2*sqrt(1 - 1 / e_out2)))_GETVALUE;
                double dF_dv0y = ((v0y / sqrt(V0dotV0*VfdotVf) - vfy*(V0dotVf) / (sqrt(V0dotV0*VfdotVf)*VfdotVf)) / sqrt(1 - V0dotVf*V0dotVf / (V0dotV0*VfdotVf)) - 2 * rp*vfy / (mu*e_out2*sqrt(1 - 1 / e_out2)))_GETVALUE;
                double dF_dv0z = ((v0z / sqrt(V0dotV0*VfdotVf) - vfz*(V0dotVf) / (sqrt(V0dotV0*VfdotVf)*VfdotVf)) / sqrt(1 - V0dotVf*V0dotVf / (V0dotV0*VfdotVf)) - 2 * rp*vfz / (mu*e_out2*sqrt(1 - 1 / e_out2)))_GETVALUE;

                Gindex = Gindices_TurnAngleConstraint_with_respect_to_Vinfinity[0][1];
                Xindex = this->myOptions->jGvar[Gindex];
                G[Gindex] = this->myOptions->X_scale_factors[Xindex] * dF_dv0x;

                Gindex = Gindices_TurnAngleConstraint_with_respect_to_Vinfinity[1][1];
                Xindex = this->myOptions->jGvar[Gindex];
                G[Gindex] = this->myOptions->X_scale_factors[Xindex] * dF_dv0y;

                Gindex = Gindices_TurnAngleConstraint_with_respect_to_Vinfinity[2][1];
                Xindex = this->myOptions->jGvar[Gindex];
                G[Gindex] = this->myOptions->X_scale_factors[Xindex] * dF_dv0z;

                Gindex = Gindices_TurnAngleConstraint_with_respect_to_Vinfinity[0][0];
                Xindex = this->myOptions->jGvar[Gindex];
                G[Gindex] = this->myOptions->X_scale_factors[Xindex] * dF_dvfx;

                Gindex = Gindices_TurnAngleConstraint_with_respect_to_Vinfinity[1][0];
                Xindex = this->myOptions->jGvar[Gindex];
                G[Gindex] = this->myOptions->X_scale_factors[Xindex] * dF_dvfy;

                Gindex = Gindices_TurnAngleConstraint_with_respect_to_Vinfinity[2][0];
                Xindex = this->myOptions->jGvar[Gindex];
                G[Gindex] = this->myOptions->X_scale_factors[Xindex] * dF_dvfz;
            }//end derivatives
        }//end process_event_main()

         //******************************************output methods
        void EphemerisPeggedUnpoweredFlyby::output(std::ofstream& outputfile,
            const double& launchdate,
            size_t& eventcount)
        {
            std::string event_type = "upwr_flyby";

            std::string boundary_name = this->myBody->name;

            math::Matrix<doubleType> empty3vector(3, 1, 0.0);

            this->mySpacecraft->computePowerState(this->state_after_event.getSubMatrix1D(0, 3).norm() / this->myUniverse->LU,
                (this->EventRightEpoch - launchdate), this->myOptions->power_margin);

            //rotate the v-infinity vector to the local frame
            this->myUniverse->LocalFrame.construct_rotation_matrices(this->EventLeftEpoch / 86400.0 + 2400000.5, false);

            math::Matrix<doubleType> R_ICRF_to_local = this->myUniverse->LocalFrame.get_R(EMTG::ReferenceFrame::ICRF,
                this->myOptions->Journeys[this->journeyIndex].journey_arrival_elements_frame);

            math::Matrix<doubleType> VinfinityLocalFrame = R_ICRF_to_local * this->Vinfinity_out;

            //compute RA and DEC
            doubleType RA = atan2(VinfinityLocalFrame(1), VinfinityLocalFrame(0));
            doubleType DEC = asin(VinfinityLocalFrame(2) / VinfinityLocalFrame.norm());

            //compute BdotR and BdotT
            math::Matrix<doubleType> periapse_state = this->calculate_flyby_periapse_state();
            math::Matrix<doubleType> periapse_R(3, 1);
            periapse_R(0) = periapse_state(0);
            periapse_R(1) = periapse_state(1);
            periapse_R(2) = periapse_state(2);
            this->myBplane.define_bplane(this->Vinfinity_in, this->state_before_event);
            this->myBplane.compute_BdotR_BdotT_from_periapse_position(this->myBody->mu, this->Vinfinity_in, periapse_R, this->BdotR, this->BdotT);

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
                this->mySpacecraft->getProducedPower(),
                0.0,
                0,
                0.0);
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
        }
    }//end namespace BoundaryEvents
}//end namespace EMTG