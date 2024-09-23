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

//base subphase for EMTGv9 MGAnDSMs
//Jacob Englander 8-23-2017

#include "MGAnDSMs_subphase.h"
#include "mjd_to_mdyhms.h"

#include "EMTG_solver_utilities.h"


#include <iostream>
#include <fstream>
#include <sstream>

namespace EMTG
{
    namespace Phases
    {
        MGAnDSMs_subphase::MGAnDSMs_subphase() :
            chemical_fuel_used(0.0),
            chemical_oxidizer_used(0.0),
            DSM_fuel_used(0.0),
            TCM_fuel_used(0.0),
            hasTCM(false),
            DSM_magnitude(0.0),
            TCM_magnitude(0.0),
            DSM(math::Matrix<doubleType>(3, 1, 0.0)),
            ChemicalManeuverType(PropulsionSystemChoice::Biprop),
            dFuel_dDSMcomponents(math::Matrix<double>(3, 1, 0.0)),
            dOxidizer_dDSMcomponents(math::Matrix<double>(3, 1, 0.0))
        {

            this->StateAfterPropagationBeforeDSM = math::Matrix<doubleType>(10, 1, 0.0);
            this->StateAfterDSMBeforeTCM = math::Matrix<doubleType>(10, 1, 0.0);
            this->output_state.resize(10 + 11 * 11, 1, 0.0);
            this->empty3.resize(3, 1, 0.0);
            this->myPropagator = nullptr;
            this->mySpacecraftAccelerationModel = nullptr;
            this->myIntegrationScheme = nullptr;
        }

        MGAnDSMs_subphase::MGAnDSMs_subphase(const std::string& name,
            const size_t& journeyIndex,
            const size_t& phaseIndex,
            const size_t& subphaseIndex,
            const size_t& stageIndex,
            MGAnDSMs_phase* myPhase,
            MGAnDSMs_subphase* previousSubPhase,
            Astrodynamics::universe* myUniverse,
            HardwareModels::Spacecraft* mySpacecraft,
            missionoptions* myOptions) :
            MGAnDSMs_subphase()
        {
            //initialize the writey thing
            this->writey_thing::initialize(myOptions, myUniverse);

            //initialize other things
            this->journeyIndex = journeyIndex;
            this->phaseIndex = phaseIndex;
            this->subphaseIndex = subphaseIndex;
            this->stageIndex = stageIndex;
            this->myUniverse = myUniverse;
            this->mySpacecraft = mySpacecraft;
            this->myJourneyOptions = &(this->myOptions->Journeys[this->journeyIndex]);
            this->myPhase = myPhase;
            this->previousSubPhase = previousSubPhase;

            this->num_timesteps = this->myJourneyOptions->override_num_steps
                ? this->myJourneyOptions->number_of_steps
                : this->myOptions->num_timesteps;

            if (this->myJourneyOptions->override_integration_step_size)
            {
                this->EphemerisOutputResolution = this->myJourneyOptions->integration_step_size;
            }
            else
            {
                this->EphemerisOutputResolution = this->myOptions->integration_time_step_size;
            }

            this->hasTCM = this->myOptions->TCM_maneuver_fraction > 1.0e-10 ? true : false;

            if (this->myJourneyOptions->override_PropagatorType)
            {
                this->myPropagatorType = this->myJourneyOptions->propagatorType;
            }
            else
            {
                this->myPropagatorType = this->myOptions->propagatorType;
            }

            if (this->myPropagatorType == PropagatorType::KeplerPropagator)
            {
                this->isKeplerian = true;
                this->STM = math::Matrix<double>(6, math::identity);
            }
            else //integrated propagator
            {
                this->isKeplerian = false;
                this->STM = math::Matrix<double>(11, math::identity);
            }
        }//end constructor

        MGAnDSMs_subphase::~MGAnDSMs_subphase()
        {
            delete this->myPropagator;
            
            if (this->myPropagatorType == PropagatorType::IntegratedPropagator)
            {
                //delete force modely things
                delete this->mySpacecraftAccelerationModel;
                delete this->myIntegrationScheme;
            }
        }

        void MGAnDSMs_subphase::setup_calcbounds(std::vector<double>* Xupperbounds,
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
            this->prefix = this->name + ": ";

            this->sparsey_thing::setup_calcbounds(Xupperbounds,
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
            
            for (size_t constraintIndex = 0; constraintIndex < this->myManeuverConstraints.size(); ++constraintIndex)
                this->myManeuverConstraints[constraintIndex].setup_calcbounds(Xupperbounds,
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
        }//end setup_calcbounds()

        //************************************calcbounds methods
        void MGAnDSMs_subphase::calcbounds_DSM_components()
        {
            //DSM components
            std::vector<std::string> DSM_component_names({ "x", "y", "z" });
            for (size_t Vindex = 0; Vindex < 3; ++Vindex)
            {
                Xlowerbounds->push_back(-10.0);
                Xupperbounds->push_back(10.0);
                X_scale_factors->push_back(10.0);
                Xdescriptions->push_back(prefix + "DSM " + DSM_component_names[Vindex] + " component");
                this->Xindex_DSM_components.push_back(Xdescriptions->size() - 1);
            }
        }//end calcbounds_DSM_components()

        void MGAnDSMs_subphase::calcbounds_maneuver_constraints()
        {
            for (size_t constraintIndex = 0; constraintIndex < this->myManeuverConstraints.size(); ++constraintIndex)
                this->myManeuverConstraints[constraintIndex].calcbounds();
        }//end calcbounds_maneuver_constraints()

         //******************************************process methods

        void MGAnDSMs_subphase::process_DSM_components(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //DSM components
            for (size_t Vindex = 0; Vindex < 3; ++Vindex)
            {
                this->DSM(Vindex) = X[Xindex++];
            }

            this->DSM_magnitude = this->DSM.norm() + 1.0e-25;
        }//end process_DSM_components()

        void MGAnDSMs_subphase::process_maneuver_constraints(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            for (size_t constraintIndex = 0; constraintIndex < this->myManeuverConstraints.size(); ++constraintIndex)
                this->myManeuverConstraints[constraintIndex].process_constraint(X, Xindex, F, Findex, G, needG);
        }//end process_maneuver_constraints

         //******************************************output methods
        void MGAnDSMs_subphase::write_output_line(std::ofstream& outputfile,
                size_t& eventcount,
                const std::string& event_type,
                const std::string& location,
                const doubleType& timestep_size,
                const doubleType& angle1,
                const doubleType& angle2,
                const math::Matrix<doubleType>& state,
                const math::Matrix<doubleType>& dV,
                const doubleType dVmag,
                const doubleType Isp)
        {
            //create a 3-element storage vector that will be used every time something needs to be rotated to the local frame
            math::Matrix<doubleType> display_vector(3, 1, 0.0);
            math::Matrix<doubleType> rot_in_vec(3, 1, 0.0);

            outputfile.width(5); outputfile << eventcount++;
            outputfile.width(3); outputfile << " | ";

            //output the event epoch in both MJD and MM/DD/YYYY
            double current_epoch_MJD = state(7) _GETVALUE / 86400.0;
            outputfile.width(16); outputfile.setf(std::ios::fixed, std::ios::floatfield); outputfile.precision(8); outputfile << current_epoch_MJD + 2400000.5;
            outputfile.width(3); outputfile << " | ";
            int month, day, year, hrs, mins;
            double secs;
            mjd_to_mdyhms(current_epoch_MJD, &month, &day, &year, &hrs, &mins, &secs);
            std::stringstream datestream;
            datestream << month << "/" << day << "/" << year;
            outputfile.width(11); outputfile << datestream.str();
            outputfile.width(3); outputfile << " | ";

            outputfile.width(12); outputfile << event_type;
            outputfile.width(3); outputfile << " | ";

            outputfile.width(25); outputfile << location;
            outputfile.width(3); outputfile << " | ";

            outputfile.width(15); outputfile << timestep_size _GETVALUE;
            outputfile.width(3); outputfile << " | ";

            //rp, BdotR, and BdotT
            outputfile.width(19); outputfile << "-";
            outputfile.width(3); outputfile << " | ";
            outputfile.width(19); outputfile << "-";
            outputfile.width(3); outputfile << " | ";
            outputfile.width(19); outputfile << "-";
            outputfile.width(3); outputfile << " | ";

            //control angles
            if (!(event_type == "coast" || event_type == "Scoast" || event_type == "TCM"))
            {
                outputfile.precision(3); outputfile.width(8); outputfile << fmod(angle1 _GETVALUE, 2.0*math::PI) * 180.0 / math::PI;
                outputfile.width(3); outputfile << " | ";
                outputfile.precision(3); outputfile.width(8); outputfile << fmod(angle2 _GETVALUE, 2.0*math::PI) * 180.0 / math::PI;
                outputfile.width(3); outputfile << " | ";
            }
            else
            {
                outputfile.precision(3); outputfile.width(8); outputfile << 0.0;
                outputfile.width(3); outputfile << " | ";
                outputfile.precision(3); outputfile.width(8); outputfile << 0.0;
                outputfile.width(3); outputfile << " | ";
            }

            //C3 - always zero
            outputfile.width(14); outputfile << "-";
            outputfile.width(3); outputfile << " | ";


            //state at event +
            outputfile.precision(8);
            for (size_t k = 0; k < 3; ++k)
                rot_in_vec(k) = state(k);
            this->myUniverse->LocalFrame.rotate_frame_to_frame(ReferenceFrame::ICRF, rot_in_vec, this->myOptions->output_file_frame, display_vector);
            for (size_t k = 0; k<3; ++k)
            {
                outputfile.width(19); outputfile << display_vector(k) _GETVALUE;
                outputfile.width(3); outputfile << " | ";
            }
            outputfile.precision(8);
            for (size_t k = 0; k < 3; ++k)
                rot_in_vec(k) = state(k + 3);
            this->myUniverse->LocalFrame.rotate_frame_to_frame(ReferenceFrame::ICRF, rot_in_vec, this->myOptions->output_file_frame, display_vector);
            for (size_t k = 0; k<3; ++k)
            {
                outputfile.width(19); outputfile << display_vector(k) _GETVALUE;
                outputfile.width(3); outputfile << " | ";
            }


            //deltaV - only for burns
            if (event_type == "chem_burn")

            {
                outputfile.precision(12);
                for (size_t k = 0; k < 3; ++k)
                    rot_in_vec(k) = dV(k);
                this->myUniverse->LocalFrame.rotate_frame_to_frame(ReferenceFrame::ICRF, rot_in_vec, this->myOptions->output_file_frame, display_vector);
                for (int k = 0; k < 3; ++k)
                {
                    outputfile.width(19); outputfile << display_vector(k) _GETVALUE;
                    outputfile.width(3); outputfile << " | ";
                }
            } 
            else
            {
                for (int k = 0; k < 3; ++k)
                {
                    outputfile.width(19); outputfile << "-";
                    outputfile.width(3); outputfile << " | ";
                }
            }

            //Thrust
            for (int k = 0; k < 3; ++k)
            {
                outputfile.width(19); outputfile << "-";
                outputfile.width(3); outputfile << " | ";
            }

            //dV magnitude
            outputfile.precision(5);
            outputfile.width(17); outputfile << dVmag _GETVALUE;
            outputfile.width(3); outputfile << " | ";

            //thrust, Isp, power
            outputfile.precision(5);
            if (event_type == "coast"
                || event_type == "Scoast"
                || event_type == "force-coast")
            {
                outputfile.width(14); outputfile << "-";
                outputfile.width(3); outputfile << " | ";
            }
            else
            {
                outputfile.width(14); outputfile << "impulse";
                outputfile.width(3); outputfile << " | ";
            }

            outputfile.precision(0);
            if (event_type == "coast"
                || event_type == "Scoast"
                || event_type == "force-coast")
            {
                outputfile.width(14); outputfile << "-";
                outputfile.width(3); outputfile << " | ";
            }
            else if (Isp > 0)
            {
                outputfile.width(14); outputfile << Isp _GETVALUE;
                outputfile.width(3); outputfile << " | ";
            }
            else
            {
                outputfile.width(14); outputfile << "UNHANDLED EVENT TYPE";
                outputfile.width(3); outputfile << " | ";
            }

            math::Matrix<doubleType> R_sc_Sun(3, 1, 0.0);
            if (this->myUniverse->central_body_SPICE_ID == 10)
            {
                R_sc_Sun = state.getSubMatrix1D(0, 2);
            }
            else
            {
                //where is the central body relative to the sun?
                doubleType central_body_state_and_derivatives[12];
                this->myUniverse->locate_central_body(state(7),
                    central_body_state_and_derivatives,
                    *this->myOptions,
                    false);

                math::Matrix<doubleType> R_CB_Sun(3, 1, 0.0);
                for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
                {
                    R_CB_Sun(stateIndex) = central_body_state_and_derivatives[stateIndex];
                }

                R_sc_Sun = state.getSubMatrix1D(0, 2) + R_CB_Sun;
            }
            doubleType r_sc_sun_AU = R_sc_Sun.norm() / this->myOptions->AU;
            this->mySpacecraft->computePowerState(r_sc_sun_AU, state(7));

            doubleType power = this->mySpacecraft->getAvailablePower();
            doubleType active_power = this->mySpacecraft->getEPActivePower();

            //power
            outputfile.precision(5);
            outputfile.width(14); outputfile << power _GETVALUE;
            outputfile.width(3); outputfile << " | ";
            
            //mass flow rate
            outputfile.width(19); outputfile << "-";
            outputfile.width(3); outputfile << " | ";

            //mass
            outputfile.width(14); outputfile.precision(4); outputfile << state(6) _GETVALUE;

            outputfile.width(3); outputfile << " | ";

            //number of active engines
            outputfile.width(14);
            outputfile << "-";

            //active power
            outputfile.width(3); outputfile << " | ";
            outputfile.precision(5);
            outputfile.width(14); outputfile << active_power _GETVALUE;
            outputfile.width(3); outputfile << " | ";

            //throttle level
            outputfile.width(14); outputfile << "-";
            outputfile.width(3); outputfile << " | ";

            outputfile << std::endl;
        }
    }//close namespace Phases
}//close namespace EMTG