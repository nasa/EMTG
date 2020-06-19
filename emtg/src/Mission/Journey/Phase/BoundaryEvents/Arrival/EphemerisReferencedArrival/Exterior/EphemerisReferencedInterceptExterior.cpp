
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

#include "EphemerisReferencedInterceptExterior.h"
#include "bplane.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        EphemerisReferencedInterceptExterior::EphemerisReferencedInterceptExterior(const std::string& name,
            const size_t& journeyIndex,
            const size_t& phaseIndex,
            size_t& stageIndex,
            Astrodynamics::universe* Universe,
            HardwareModels::Spacecraft* mySpacecraft,
            missionoptions* myOptions) :
            EphemerisReferencedArrivalWithVinfinityExterior::EphemerisReferencedArrivalWithVinfinityExterior(name,
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
        void EphemerisReferencedInterceptExterior::calcbounds(std::vector<size_t> timeVariables)
        {
            std::vector<double> RAbounds({ -2880.0, 2880.0 });
            std::vector<double> DECbounds({ -90.0, 90.0 });
            std::vector<double> MassBounds({1.0e-13, this->myJourneyOptions->maximum_mass });

            //epoch bounds
            std::vector<double> EpochBounds(2);
            if (this->isLastEventInJourney
                && this->myJourneyOptions->timebounded == 2)//bounded arrival date
            {
                EpochBounds[0] = this->myJourneyOptions->arrival_date_bounds[0];
                EpochBounds[1] = this->myJourneyOptions->arrival_date_bounds[1];
            }
            else
            {
                EpochBounds[0] = this->myOptions->launch_window_open_date;
                EpochBounds[1] = this->myBody->getEphemerisWindowClose();
            }

            this->calcbounds_event_interface_state(this->myJourneyOptions->final_velocity,
                RAbounds, 
                DECbounds,
                MassBounds, 
                EpochBounds,
                timeVariables);

            this->calcbounds_event_left_side(timeVariables);

            this->calcbounds_event_right_side();



            if (this->hasTCM)

                this->calcbounds_virtual_propellant_constraints();

            this->calcbounds_specialized_constraints();
        }//end calcbounds()

        

        //calcbounds for propellant constraints - TCM

        void EphemerisReferencedInterceptExterior::calcbounds_virtual_propellant_constraints()

        {
            if (this->hasTCM)
            {
                this->hasMonopropManeuver = true;
            }

            BoundaryEventBase::calcbounds_virtual_propellant_constraints();

            if (this->hasTCM)
            {
                //derivative of virtual chemical fuel with respect to initial mass
                this->create_sparsity_entry(this->Findex_VirtualChemicalFuelConstraint,
                    this->X_index_of_first_decision_variable_in_this_event,
                    true,
                    "mass",
                    Gindices_dVirtualChemicalFuel_TCM_dLeftMass);
            }

        }//end calcbounds_virtual_propellant_constraints

         //******************************************process methods

         void EphemerisReferencedInterceptExterior::process_event(const std::vector<doubleType>& X,
             size_t& Xindex,
             std::vector<doubleType>& F,
             size_t& Findex,
             std::vector<double>& G,
             const bool& needG)
         {
             this->process_staging();
             this->process_event_interface_state(X, Xindex, F, Findex, G, needG);
             this->process_event_left_side(X, Xindex, F, Findex, G, needG);
             this->process_event_main(X, Xindex, F, Findex, G, needG);
             this->process_event_right_side(X, Xindex, F, Findex, G, needG);

             if (this->hasTCM)
             {
                 this->process_virtual_propellant_constraints(X, Xindex, F, Findex, G, needG);
             }

             this->process_specialized_constraints(X, Xindex, F, Findex, G, needG);
         }//end process_event()



         void EphemerisReferencedInterceptExterior::process_event_main(const std::vector<doubleType>& X,
                                                                       size_t& Xindex,
                                                                       std::vector<doubleType>& F,
                                                                       size_t& Findex,
                                                                       std::vector<double>& G,
                                                                       const bool& needG)

         {
             //TCM
             if (this->hasTCM)
             {
                 //perform a TCM with the current stage monoprop system
                 this->mySpacecraft->computeChemicalPropulsionPerformance(this->TCM_magnitude, this->state_before_event(6), true, PropulsionSystemChoice::Monoprop); //monoprop
                 this->chemical_fuel_used = this->mySpacecraft->getChemFuelConsumedThisManeuver();
                 this->ETM(6, 6) = ((this->state_before_event(6) - this->chemical_fuel_used) / this->state_before_event(6))_GETVALUE;
             }

         }//end process_event_main()

         void EphemerisReferencedInterceptExterior::process_virtual_propellant_constraints(const std::vector<doubleType>& X,
                                                                                           size_t& Xindex,
                                                                                           std::vector<doubleType>& F,
                                                                                           size_t& Findex,
                                                                                           std::vector<double>& G,
                                                                                           const bool& needG)
         {
             //call the base class propellant constraints
             BoundaryEventBase::process_virtual_propellant_constraints(X, Xindex, F, Findex, G, needG);

             if (this->hasTCM && needG)
             {
                 size_t Gindex = this->Gindices_dVirtualChemicalFuel_TCM_dLeftMass;
                 size_t Xindex = this->jGvar->operator[](Gindex);
                 double Gentry = (this->chemical_fuel_used / this->state_before_event(6)) _GETVALUE;
                 G[Gindex] = this->X_scale_factors->operator[](Xindex) * -Gentry
                     * this->myUniverse->continuity_constraint_scale_factors(6);
             }

         }//end process_virtual_propellant_constraints()

         //******************************************output methods
         void EphemerisReferencedInterceptExterior::output(std::ofstream& outputfile,
             const double& launchdate,
             size_t& eventcount)
         {
             this->mySpacecraft->setActiveStage(this->stageIndex);
             std::string event_type = "intercept";
             std::string boundary_name = this->myBody->name + "_BE"; //"_BE" means "boundary ellipsoid"
             math::Matrix<doubleType> empty3vector(3, 1, 0.0);
             
             //where is the Sun?
             math::Matrix<doubleType> R_sc_Sun(3, 1, 0.0);
             if (this->myBody->spice_ID == 10)
             {
                 R_sc_Sun = this->state_before_event.getSubMatrix1D(0, 2);
             }
             else
             {
                 //where is the central body relative to the sun?
                 doubleType central_body_state_and_derivatives[12];
                 this->myUniverse->locate_central_body(this->state_before_event(7),
                     central_body_state_and_derivatives,
                     *this->myOptions,
                     false);

                 math::Matrix<doubleType> R_CB_Sun(3, 1, 0.0);
                 for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
                 {
                     R_CB_Sun(stateIndex) = central_body_state_and_derivatives[stateIndex];
                 }

                 R_sc_Sun = this->state_before_event.getSubMatrix1D(0, 2) + R_CB_Sun;
             }

             this->mySpacecraft->computePowerState(R_sc_Sun.getSubMatrix1D(0, 2).norm() / this->myOptions->AU, this->state_before_event(7));

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
                     this->state_before_event,//state
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

             math::Matrix<doubleType> Vinfinity_in = this->state_on_interface_cartesian.getSubMatrix1D(3, 5);

             math::Matrix<doubleType> VinfinityLocalFrame = R_ICRF_to_local * Vinfinity_in;

             //make an output state in the frame of the left side of the event but with the mass after the event
             math::Matrix<doubleType> output_state = this->state_before_event;
             output_state(6) = this->state_after_event(6);

             //compute RA and DEC
             RA = atan2(VinfinityLocalFrame(1), VinfinityLocalFrame(0));
             DEC = asin(VinfinityLocalFrame(2) / VinfinityLocalFrame.norm());
             //we print state BEFORE event here, and we will print state AFTER event in the next journey
             //that way we get the right reference frame!
             write_output_line(outputfile,
                 eventcount,
                 event_type,
                 boundary_name,
                 this->EventTimeWidth,
                 0.0,
                 0.0,
                 0.0,
                 RA,
                 DEC,
                 this->state_on_interface_spherical(3) * this->state_on_interface_spherical(3),
                 output_state,
                 -Vinfinity_in,
                 empty3vector,
                 this->state_on_interface_spherical(3),
                 0.0,
                 0.0,
                 this->mySpacecraft->getAvailablePower(),
                 0.0,
                 0,
                 0.0,
                 "none");
         }//end output()

         void EphemerisReferencedInterceptExterior::output_maneuver_and_target_spec(std::ofstream& maneuver_spec_file, std::ofstream& target_spec_file, bool& haveManeuverNeedTarget)
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
                 myBplane.define_bplane(this->state_after_event);

                 target_spec_line myTargetSpecLine(this->name,
                     "EME2000",
                     this->myBody->name,
                     this->state_after_event(7),
                     this->state_after_event,
                     myBplane.getBdotR(),
                     myBplane.getBdotT());

                 //Step 1.2: write target spec object
                 myTargetSpecLine.write(target_spec_file);
             }
         }//end output_maneuver_and_target_spec()

    }//end namespace BoundaryEvents
}//end namespace EMTG