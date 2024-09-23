
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

#include "EphemerisPeggedMomentumTransfer.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        EphemerisPeggedMomentumTransfer::EphemerisPeggedMomentumTransfer(const std::string& name,
            const size_t& journeyIndex,
            const size_t& phaseIndex,
            size_t& stageIndex,
            Astrodynamics::universe* Universe,
            HardwareModels::Spacecraft* mySpacecraft,
            missionoptions* myOptions) :
            EphemerisPeggedIntercept::EphemerisPeggedIntercept(name,
                journeyIndex,
                phaseIndex,
                stageIndex,
                Universe,
                mySpacecraft,
                myOptions)
        {
			this->deltaV.resize(3, 1, 0.0);
            this->ddeltaV_dSpacecraftMass.resize(3, 1, 0.0);

            //this is the really scary part. We have to change the "maximum_mass" field in the journeyoptions for all FUTURE journeys
            //otherwise all encoded masses won't include the mass of the body that we are transferring momentum to
            for (size_t futureJourneyIndex = this->journeyIndex + 1; futureJourneyIndex < this->myOptions->number_of_journeys; ++futureJourneyIndex)
                this->myOptions->Journeys[futureJourneyIndex].maximum_mass = this->myBody->mass;
        }

        //******************************************calcbounds methods

        //calcbounds
        void EphemerisPeggedMomentumTransfer::calcbounds(std::vector<size_t> timeVariables)
        {
            this->calcbounds_event_left_side(timeVariables);

            this->calcbounds_event_main();

            this->calcbounds_event_right_side();

            if (this->hasTCM)
                this->calcbounds_virtual_propellant_constraints();

            this->calcbounds_specialized_constraints();
        }//end calcbounds()

        void EphemerisPeggedMomentumTransfer::calcbounds_event_main()
        {
            //base class
            this->EphemerisPeggedIntercept::calcbounds_event_main();

            //create derivatives of velocity after event with respect to spacecraft mass
            for (size_t Vindex = 0; Vindex < 3; ++Vindex)
            {
                this->Derivatives_of_StateAfterEvent.push_back({ this->Xindex_mass, Vindex + 3, 1.0 });
                this->dIndex_ddeltaV_dSpacecraftMass.push_back(this->Derivatives_of_StateBeforeEvent.size() - 1);
            }
        }//end calcbounds_event_main()

        //******************************************process methods

        void EphemerisPeggedMomentumTransfer::process_event(const std::vector<doubleType>& X,
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

        void EphemerisPeggedMomentumTransfer::process_event_main(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //call the base class
            this->EphemerisPeggedIntercept::process_event_main(X, Xindex, F, Findex, G, needG);

			double beta_impact = this->myJourneyOptions->impact_momentum_enhancement_factor;
			doubleType mass_SC = this->state_before_event(6);
			double mass_asteroid = this->myBody->mass;

			// Velocity imparted on the asteroid by the S/C's impact
			this->ddeltav_dVinfinity_in = (beta_impact * mass_SC / (mass_SC + mass_asteroid)) _GETVALUE;
            this->ddeltaV_dSpacecraftMass(0) = (this->Vinfinity_in(0) * beta_impact * (mass_asteroid / (mass_SC + mass_asteroid) / (mass_SC + mass_asteroid))) _GETVALUE;
            this->ddeltaV_dSpacecraftMass(1) = (this->Vinfinity_in(1) * beta_impact * (mass_asteroid / (mass_SC + mass_asteroid) / (mass_SC + mass_asteroid))) _GETVALUE;
            this->ddeltaV_dSpacecraftMass(2) = (this->Vinfinity_in(2) * beta_impact * (mass_asteroid / (mass_SC + mass_asteroid) / (mass_SC + mass_asteroid))) _GETVALUE;

			this->deltaV(0) = this->Vinfinity_in(0) * this->ddeltav_dVinfinity_in;		// Delta_Vx of impact on the asteroid by the S/C
			this->deltaV(1) = this->Vinfinity_in(1) * this->ddeltav_dVinfinity_in;		// Delta_Vy of impact on the asteroid by the S/C
			this->deltaV(2) = this->Vinfinity_in(2) * this->ddeltav_dVinfinity_in;		// Delta_Vz of impact on the asteroid by the S/C
        }//end process_event_main()


		void EphemerisPeggedMomentumTransfer::
			process_event_right_side(const std::vector<doubleType>& X,
				size_t& Xindex,
				std::vector<doubleType>& F,
				size_t& Findex,
				std::vector<double>& G,
				const bool& needG)
		{
			//Step 1: base class
			this->EphemerisPeggedArrival::process_event_right_side(X, Xindex, F, Findex, G, needG);

			//Step 2: modify the velocity and mass AFTER the momentum exchange event
			this->state_after_event(3) += (this->deltaV(0) - this->Vinfinity_in(0));
			this->state_after_event(4) += (this->deltaV(1) - this->Vinfinity_in(1));
			this->state_after_event(5) += (this->deltaV(2) - this->Vinfinity_in(2));

			//Step 3: modify the partials of state_after_event with respect to v-infinity and mass to properly reflect what happened
			//Step 3.1: partials with respect to v-infinity
            for (size_t dindex : this->dIndex_VbeforeEvent_dVinfinity_in)
			{
				std::get<2>(this->Derivatives_of_StateAfterEvent[dindex]) = this->ddeltav_dVinfinity_in;
			}
            //Step 3.2: partials with respect to mass
            for (size_t vIndex : {0, 1, 2})
            {
                size_t dIndex = this->dIndex_ddeltaV_dSpacecraftMass[vIndex];
                std::get<2>(this->Derivatives_of_StateAfterEvent[dIndex]) = this->ddeltaV_dSpacecraftMass(vIndex);
            }
            

			//Step 4: the mass is now the mass of the asteroid
			this->state_after_event(6) = this->myBody->mass; //have a conversation about ejecta later
			std::get<2>(this->Derivatives_of_StateAfterEvent[this->dIndex_mass_wrt_encodedMass]) = 0.0; //we're just using the mass of the asteroid
		}//end process_event_right_side()

         //******************************************output methods
        void EphemerisPeggedMomentumTransfer::output(std::ofstream& outputfile,
            const double& launchdate,
            size_t& eventcount)
        {
            this->mySpacecraft->setActiveStage(this->stageIndex);

            std::string event_type = "momtransfer";

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
                    power,//AvailPower
                    0.0,//mdot
                    0,//number_of_active_engines
                    0.0,
                    "none");
            }

            //compute RA and DEC for incoming asymptote in the body's local frame
            doubleType RA, DEC;
            //rotate the v-infinity vector to the local frame
            this->myUniverse->LocalFrame.construct_rotation_matrices(this->EventLeftEpoch / 86400.0 + 2400000.5, false);

            math::Matrix<doubleType> R_ICRF_to_local = this->myUniverse->LocalFrame.get_R(EMTG::ReferenceFrame::ICRF,
                this->myOptions->Journeys[this->journeyIndex].arrival_elements_frame);

            math::Matrix<doubleType> VinfinityLocalFrame = R_ICRF_to_local * this->Vinfinity_in;

            //compute RA and DEC
            RA = atan2(VinfinityLocalFrame(1), VinfinityLocalFrame(0));
            DEC = math::safe_asin(VinfinityLocalFrame(2) / (VinfinityLocalFrame.norm() + math::SMALL));

            this->write_output_line(outputfile,
                eventcount,
                event_type,
                boundary_name,
                this->EventTimeWidth,
                0.0,
                0.0,
                0.0,
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
    }//end namespace BoundaryEvents
}//end namespace EMTG