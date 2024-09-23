
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

#include "EphemerisPeggedMomentumExchange.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        EphemerisPeggedMomentumExchange::EphemerisPeggedMomentumExchange(const std::string& name,
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
        void EphemerisPeggedMomentumExchange::calcbounds()
        {

            this->calcbounds_event_left_side();

            this->calcbounds_event_main();

            this->calcbounds_event_right_side();

            //if (this->hasTCM)
            //    this->calcbounds_virtual_propellant_constraints();

            this->calcbounds_specialized_constraints(); //?????????????????????????????????????????????????????????????????????????????????????????????????????????????
        }//end calcbounds()

        void EphemerisPeggedIntercept::calcbounds_event_main()
        {
			this->myArrivalEvent = this->myMission->getJourney(this->myOptions->number_of_journeys - 1)->getPhase(number_of_phases_in_final_journey - 1)->getArrivalEvent();
			math::Matrix<doubleType> FinalState = this->myArrivalEvent->get_state_after_event();

			if (momentum_exchange_type == 'kinetic impactor')
			{
				double beta_impact   = myOptions.momentum_enhancement_factor;
				double mass_SC		 = FinalState(6);
				double number_of_SC	 = myOptions.number_of_SC;
				double mass_asteroid = Universe.bodies[options.Journeys[myOptions.number_of_journeys - 1].destination_list[1] - 1].mass;
				
				// Velocity imparted on the asteroid by the S/C's impact
				Delta_v[0] = v_impact[0] * beta_impact * mass_SC * number_of_SC / (mass_SC * number_of_SC + mass_asteroid);		// Delta_Vx of impact on the asteroid by the S/C
				Delta_v[1] = v_impact[1] * beta_impact * mass_SC * number_of_SC / (mass_SC * number_of_SC + mass_asteroid);		// Delta_Vy of impact on the asteroid by the S/C
				Delta_v[2] = v_impact[2] * beta_impact * mass_SC * number_of_SC / (mass_SC * number_of_SC + mass_asteroid);		// Delta_Vz of impact on the asteroid by the S/C				
			}
			else if (momentum_exchange_type == 'explosive device')
			{
				double explosive_magnitude = myOptions.explosive_device_magnitude;
				if (arrival_type == 'flyby')
				{
					math::Matrix<doubleType> position_wrt_target;
					double x = FinalState(1);
					double y = FinalState(2);
					double z = FinalState(3);
					
					double v_explosive_device_azimuth = math::safe_asin(y/x)/EMTG_math::safe_acos(y/x);
					double v_explosive_device_elevation = math::safe_asin(z/ sqrt(x. ^ 2 + y. ^ 2)))/ math::safe_acos(z/sqrt(x. ^ 2 + y. ^ 2)));


				}
				else if (arrival_type == 'rendezvous')
				{
					double v_explosive_device_azimuth;
					double v_explosive_device_elevation;

					for (int entry = Xdescriptions.size() - 1; entry >= 0; --entry)
					{
						if (Xdescriptions[entry].find(prefix + "v_explosive_device_azimuth") < 1024)
						{
							v_explosive_device_azimuth = X[entry];
						}
						else if (Xdescriptions[entry].find(prefix + "v_explosive_device_elevation") < 1024)
						{
							v_explosive_device_elevation = X[entry];
						}
					}
				}
				
				// Velocity imparted on the asteroid by the NED explosion
				Delta_v[0] = explosive_magnitude * cos(v_explosive_device_azimuth)*cos(v_explosive_device_elevation);		// Delta_Vx of explosion on the asteroid by the S/C
				Delta_v[1] = explosive_magnitude * sin(v_explosive_device_azimuth)*cos(v_explosive_device_elevation);		// Delta_Vy of explosion on the asteroid by the S/C
				Delta_v[2] = explosive_magnitude * sin(v_explosive_device_elevation);										// Delta_Vz of explosion on the asteroid by the S/C
			}
			else
			{
				cout << "Not yet implemented" << endl;
				exit(EXIT_FAILURE);
			}
            
			//create bounds for the v-infinity
            double vinf_max = this->myOptions->Journeys[this->journeyIndex].journey_final_velocity[1];

            std::vector< std::tuple<double, double> > vinfBounds = std::vector< std::tuple<double, double> >(3, std::make_tuple(-vinf_max, vinf_max));

            //base class
            this->EphemerisPeggedArrivalWithVinfinity::calcbounds_event_main(vinfBounds);

            //v-infinity magnitude constraint
            Flowerbounds->push_back(this->myJourneyOptions->journey_final_velocity[0] - this->myJourneyOptions->journey_final_velocity[1]);
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

        void EphemerisPeggedMomentumExchange::process_event(const std::vector<doubleType>& X,
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

            /*if (this->hasTCM)
                this->process_virtual_propellant_constraints(X, Xindex, F, Findex, G, needG);*/

            this->BoundaryEventBase::process_specialized_constraints(X, Xindex, F, Findex, G, needG);
        }//end process_event

        void EphemerisPeggedMomentumExchange::process_event_main(const std::vector<doubleType>& X,
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
            F[Findex++] = v_infinity_in_magnitude - this->myJourneyOptions->journey_final_velocity[1];

            if (needG)
            {
                for (size_t Vindex = 0; Vindex < 3; ++Vindex)
                {
                    size_t Gindex = this->Gindices_incoming_velocity_magnitude_wrt_Vinfinity_in[Vindex];
                    size_t Xindex = this->myOptions->jGvar[Gindex];

                    G[Gindex] = this->X_scale_factors->operator[](Xindex)
                        * (this->Vinfinity_in(Vindex) / v_infinity_in_magnitude)_GETVALUE;
                }
            }
        }//end process_event_main()

         //******************************************output methods
        void EphemerisPeggedMomentumExchange::output(std::ofstream& outputfile,
            const double& launchdate,
            size_t& eventcount)
        {
            this->mySpacecraft->setActiveStage(this->stageIndex);

            std::string event_type = "MomentumExchange";

            std::string boundary_name = this->myBody->name;

            math::Matrix<doubleType> empty3vector(3, 1, 0.0);

            this->mySpacecraft->computePowerState(this->state_after_event.getSubMatrix1D(0, 3).norm() / this->myUniverse->LU,
                (this->EventRightEpoch - launchdate), this->myOptions->power_margin);

            doubleType power = this->mySpacecraft->getAvailablePower();

            ////TCM if applicable
            //if (this->hasTCM)
            //{
            //    this->write_output_line(outputfile,//outputfile
            //        eventcount,//eventcount
            //        "TCM",//event_type
            //        "deep-space",//event_location
            //        0.0,// timestep_size,
            //        -1,//flyby_altitude,
            //        0,//BdotR
            //        0,//BdotT
            //        0,//angle1
            //        0,//angle2
            //        0,//C3
            //        this->state_after_event,//state
            //        math::Matrix<doubleType>(3, 1, 0.0),//dV
            //        math::Matrix<doubleType>(3, 1, 0.0),//ThrustVector
            //        this->TCM_magnitude,//dVmag
            //        0.0,//Thrust
            //        this->mySpacecraft->getMonopropIsp(),//Isp
            //        power,//AvailPower
            //        0.0,//mdot
            //        0,//number_of_active_engines
            //        0.0,
            //        -1);//active_power)
            //}

            //compute RA and DEC for incoming asymptote in the body's local frame
            doubleType RA, DEC;
            //rotate the v-infinity vector to the local frame
            this->myUniverse->LocalFrame.construct_rotation_matrices(this->EventLeftEpoch / 86400.0 + 2400000.5, false);

            math::Matrix<doubleType> R_ICRF_to_local = this->myUniverse->LocalFrame.get_R(EMTG::ReferenceFrame::ICRF,
                this->myOptions->Journeys[this->journeyIndex].journey_arrival_elements_frame);

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
                this->mySpacecraft->getProducedPower(),
                0.0,
                0,
                0.0,
                -1);
        }
    }//end namespace BoundaryEvents
}//end namespace EMTG