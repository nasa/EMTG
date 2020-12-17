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

#include "BoundaryEventBase.h"
#include "EMTG_solver_utilities.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        BoundaryEventBase::BoundaryEventBase() :
            name("dummy_event"),
            stateNames(std::vector<std::string> ({ "x", "y", "z", "xdot", "ydot", "zdot", "mass", "epoch" })),
            COENames(std::vector<std::string>({ "SMA", "ECC", "INC", "RAAN", "AOP", "TA" })),
            isFirstEventInMission(false),
            EventHasTimeWidth(false),
            journeyIndex(0),
            phaseIndex(0),
            stageIndex(0),
            state_before_event(math::Matrix<doubleType>(8, 1, 0.0)),
            state_after_event(math::Matrix<doubleType>(8, 1, 0.0)),
            state_after_event_propagated(math::Matrix<doubleType>(8, 1, 0.0)),
            boundary_state(math::Matrix<doubleType>(8, 1, 0.0)),
            ETM(math::Matrix<double>(8, math::identity)),
            EventDeterministicDeltav(0.0),
            virtual_electric_propellant_used(0.0),
            virtual_chemical_fuel_used(0.0),
            virtual_chemical_oxidizer_used(0.0),
            electric_propellant_used(0.0),
            chemical_fuel_used(0.0),
            chemical_oxidizer_used(0.0),
            EventTimeWidth(0.0),
            hasElectricManeuver(false),
            hasMonopropManeuver(false),
            hasBipropManeuver(false),
            hasManeuver(false),
			compute_orbit_elements(false),
            TCM_magnitude(0.0),
            hasTCM(false),
            ChemicalManeuverType(PropulsionSystemChoice::Biprop),
            C3(0.0)
            {
                
            }//end default constructor

            //constructor with arguements
            BoundaryEventBase::BoundaryEventBase(const std::string& name,
                const size_t& journeyIndex,
                const size_t& phaseIndex,
                size_t& stageIndex,
                Astrodynamics::universe* Universe,
                HardwareModels::Spacecraft* mySpacecraft,
                missionoptions* myOptions) :
                BoundaryEventBase()
            {
                this->initialize(name,
                                 journeyIndex,
                                 phaseIndex,
                                 stageIndex,
                                 Universe,
                                 mySpacecraft,
                                 myOptions);
            }//end constructor with arguments

            void BoundaryEventBase::initialize(const std::string& name,
                const size_t& journeyIndex,
                const size_t& phaseIndex,
                size_t& stageIndex,
                Astrodynamics::universe* Universe,
                HardwareModels::Spacecraft* mySpacecraft,
                missionoptions* myOptions)
            {
                this->name = name;
                this->journeyIndex = journeyIndex;
                this->phaseIndex = phaseIndex;
                this->stageIndex = stageIndex;
                this->mySpacecraft = mySpacecraft;
                this->myJourneyOptions = &(myOptions->Journeys[this->journeyIndex]);

                //initialize the writey thing
                this->writey_thing::initialize(myOptions, Universe);

                //initialize the propagator type
                if (this->myJourneyOptions->override_PropagatorType)
                {
                    this->myPropagatorType = this->myJourneyOptions->propagatorType;
                }
                else
                {
                    this->myPropagatorType = this->myOptions->propagatorType;
                }
            }//end initialize()
            
            //******************************************calcbounds methods

            //function to set up calcbounds() methods
            void BoundaryEventBase::setup_calcbounds(
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


                //specialized constraints
                for (size_t constraintIndex = 0; constraintIndex < this->mySpecializedConstraints.size(); ++constraintIndex)
                    this->mySpecializedConstraints[constraintIndex].setup_calcbounds(Xupperbounds,
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
            
            void BoundaryEventBase::calculate_dependencies_left_epoch(std::vector<size_t> timeVariables)
            {
                for (size_t Xindex : timeVariables)
                {
                    this->Xindices_EventLeftEpoch.push_back(Xindex);
                    this->Derivatives_of_StateBeforeEvent_wrt_Time.push_back({ Xindex, 7, 1.0 });
                }
            }//end calculate_dependencies_left_epoch

            void BoundaryEventBase::calcbounds_specialized_constraints()
            {
				for (size_t constraintIndex = 0; constraintIndex < this->mySpecializedConstraints.size(); ++constraintIndex)
				{
					try
					{
						this->mySpecializedConstraints[constraintIndex].calcbounds();
					}
					catch (const std::logic_error & ex)
					{
						std::cerr << ex.what() << std::endl;
                        throw;
						// TODO: throw once we figure out some stuff in MGAnDSMs_subphase's destructor;
					}
				}
            }//end calcbounds_specialized_constraints

            //method to create bounds for the virtual propellant constraints
            void BoundaryEventBase::calcbounds_virtual_propellant_constraints()
            {
                //Step 1: virtual electric propellant
                if (this->hasElectricManeuver)
                {
                    this->Xlowerbounds->push_back(0.0);
                    this->Xupperbounds->push_back(this->myJourneyOptions->maximum_mass);
                    this->X_scale_factors->push_back(this->myJourneyOptions->maximum_mass);
                    this->Xdescriptions->push_back(prefix + "virtual electric propellant");
                    size_t Xindex = this->Xdescriptions->size() - 1;

                    this->Flowerbounds->push_back(-math::SMALL);
                    this->Fupperbounds->push_back(math::SMALL);
                    this->Fdescriptions->push_back(prefix + "virtual electric propellant constraint");
                    this->Findex_VirtualElectricPropellantConstraint = this->Fdescriptions->size() - 1;

                    this->create_sparsity_entry(this->Findex_VirtualElectricPropellantConstraint,
                        Xindex,
                        this->Gindex_dVirtualElectricPropellantConstraint_dVirtualElectricPropellant);

                    //tell the spacecraft about this variable
                    //global constraint
                    this->mySpacecraft->appendGlobalElectricPropellantTank_Xindices(Xindex);
                    this->mySpacecraft->appendGlobalElectricPropellantTank_Xscaleranges(Xupperbounds->back(), Xlowerbounds->back());
                    //stage constraint
                    this->mySpacecraft->appendElectricPropellantTank_Xindices(this->stageIndex, Xindex);
                    this->mySpacecraft->appendElectricPropellantTank_Xscaleranges(this->stageIndex, Xupperbounds->back(), Xlowerbounds->back());
                }//end virtual electric propellant constraint

                //Step 2: virtual chemical fuel
                if (this->hasBipropManeuver || this->hasMonopropManeuver)
                {
                    this->Xlowerbounds->push_back(0.0);
                    this->Xupperbounds->push_back(this->myJourneyOptions->maximum_mass);
                    this->X_scale_factors->push_back(this->myJourneyOptions->maximum_mass);
                    this->Xdescriptions->push_back(prefix + "virtual chemical fuel");
                    size_t Xindex = this->Xdescriptions->size() - 1;

                    this->Flowerbounds->push_back(-math::SMALL);
                    this->Fupperbounds->push_back(math::SMALL);
                    this->Fdescriptions->push_back(prefix + "virtual chemical fuel constraint");
                    this->Findex_VirtualChemicalFuelConstraint = this->Fdescriptions->size() - 1;

                    this->create_sparsity_entry(this->Findex_VirtualChemicalFuelConstraint,
                        Xindex,
                        this->Gindex_dVirtualChemicalFuelConstraint_dVirtualChemicalFuel);

                    //tell the spacecraft about this variable
                    //global constraint
                    this->mySpacecraft->appendGlobalChemicalFuelTank_Xindices(Xindex);
                    this->mySpacecraft->appendGlobalChemicalFuelTank_Xscaleranges(Xupperbounds->back(), Xlowerbounds->back());
                    //stage constraint
                    this->mySpacecraft->appendChemicalFuelTank_Xindices(this->stageIndex, Xindex);
                    this->mySpacecraft->appendChemicalFuelTank_Xscaleranges(this->stageIndex, Xupperbounds->back(), Xlowerbounds->back());
                }//end virtual chemical fuel constraint

                //Step 3: virtual chemical oxidizer
                if (this->hasBipropManeuver)
                {
                    this->Xlowerbounds->push_back(0.0);
                    this->Xupperbounds->push_back(this->myJourneyOptions->maximum_mass);
                    this->X_scale_factors->push_back(this->myJourneyOptions->maximum_mass);
                    this->Xdescriptions->push_back(prefix + "virtual chemical oxidizer");
                    size_t Xindex = this->Xdescriptions->size() - 1;

                    this->Flowerbounds->push_back(-math::SMALL);
                    this->Fupperbounds->push_back(math::SMALL);
                    this->Fdescriptions->push_back(prefix + "virtual chemical oxidizer constraint");
                    this->Findex_VirtualChemicalOxidizerConstraint = this->Fdescriptions->size() - 1;

                    this->create_sparsity_entry(this->Findex_VirtualChemicalOxidizerConstraint,
                        Xindex,
                        this->Gindex_dVirtualChemicalOxidizerConstraint_dVirtualChemicalOxidizer);

                    //tell the spacecraft about this variable
                    //global constraint
                    this->mySpacecraft->appendGlobalChemicalOxidizerTank_Xindices(Xindex);
                    this->mySpacecraft->appendGlobalChemicalOxidizerTank_Xscaleranges(Xupperbounds->back(), Xlowerbounds->back());
                    //stage constraint
                    this->mySpacecraft->appendChemicalOxidizerTank_Xindices(this->stageIndex, Xindex);
                    this->mySpacecraft->appendChemicalOxidizerTank_Xscaleranges(this->stageIndex, Xupperbounds->back(), Xlowerbounds->back());
                }//end virtual chemical oxidizer constraint
            }//end calcbounds_virtual_propellant_constraints

            //******************************************process methods
            void BoundaryEventBase::process_staging()
            {
                this->mySpacecraft->setActiveStage(this->stageIndex);
            }//end process_staging()

            void BoundaryEventBase::process_left_epoch(const std::vector<doubleType>& X,
                    size_t& Xindex,
                    std::vector<doubleType>& F,
                    size_t& Findex,
                    std::vector<double>& G,
                    const bool& needG)
            {
                //if this is the first event in the mission, we need to advance past the launch epoch!
                if (this->isFirstEventInMission)
                    ++Xindex;

                this->EventLeftEpoch = 0.0;
                for (size_t listIndex = 0; listIndex < this->Xindices_EventLeftEpoch.size(); ++listIndex)
                {
                    this->EventLeftEpoch += X[this->Xindices_EventLeftEpoch[listIndex]];
                }
                this->state_before_event(7) = this->EventLeftEpoch;

                //let's do a check to make sure that the central body exists wrt the Sun at this epoch. If it doesn't then we're going to have a kaboom so let's throw an error and make SNOPT walk back
                if (this->myUniverse->central_body.spice_ID != 10)//because the Sun always exists wrt itself{
                {
                    double latestAllowedEpoch = this->myUniverse->MySplineEphemUniverse->getEphemerisWindowClose(this->myUniverse->central_body.spice_ID, 10);
                    if (this->EventLeftEpoch > latestAllowedEpoch)
                    {
                        throw std::runtime_error(this->name + " left epoch is beyond the length of your ephemeris"
                            + " for the central body with respect to the Sun. EMTG is trying to poll the ephemeris at " + std::to_string(this->EventLeftEpoch _GETVALUE / 86400.0)
                            + " but the ephemeris ends at " + std::to_string(latestAllowedEpoch / 86400.0) + "."
                            + " Either don't worry about it because SNOPT and MBH will fix this for you, or get a longer ephemeris."
                            + " Should you choose to debug, place a breakpoint at " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
                    }
                }
            }//end process_left_epoch

            //function to process the virtual propellant constraint
            void
                BoundaryEventBase::process_virtual_propellant_constraints(const std::vector<doubleType>& X,
                                                                                        size_t& Xindex,
                                                                                        std::vector<doubleType>& F,
                                                                                        size_t& Findex,
                                                                                        std::vector<double>& G,
                                                                                        const bool& needG)
            {
                //Step 1: virtual electric propellant
                if (this->hasElectricManeuver)
                {
                    //Step 1.1: extract the virtual propellant variable
                    this->virtual_electric_propellant_used = X[Xindex++];

                    //Step 1.2: apply the virtual propellant constraint
                    F[Findex++] = (this->virtual_electric_propellant_used - this->electric_propellant_used)
                        * this->myUniverse->continuity_constraint_scale_factors(6);

                    //Step 1.3: derivatives
                    if (needG)
                    {
                        size_t Gindex = this->Gindex_dVirtualElectricPropellantConstraint_dVirtualElectricPropellant;
                        size_t Xindex = this->jGvar->operator[](Gindex);

                        G[Gindex] = this->X_scale_factors->operator[](Xindex)
                            * this->myUniverse->continuity_constraint_scale_factors(6);
                    }
                }//end virtual electric propellant

                //Step 2: virtual chemical fuel
                if (this->hasBipropManeuver || this->hasMonopropManeuver)
                {
                    //Step 2.1: extract the virtual propellant variable
                    this->virtual_chemical_fuel_used = X[Xindex++];

                    //Step 2.2: apply the virtual propellant constraint
                    F[Findex++] = (this->virtual_chemical_fuel_used - this->chemical_fuel_used)
                        * this->myUniverse->continuity_constraint_scale_factors(6);

                    //Step 2.3: derivatives
                    if (needG)
                    {
                        size_t Gindex = this->Gindex_dVirtualChemicalFuelConstraint_dVirtualChemicalFuel;
                        size_t Xindex = this->jGvar->operator[](Gindex);

                        G[Gindex] = this->X_scale_factors->operator[](Xindex)
                            * this->myUniverse->continuity_constraint_scale_factors(6);
                    }
                }//end virtual chemical fuel

                //Step 3: virtual chemical oxidizer
                if (this->hasBipropManeuver)
                {
                    //Step 3.1: extract the virtual propellant variable
                    this->virtual_chemical_oxidizer_used = X[Xindex++];

                    //Step 3.2: apply the virtual propellant constraint
                    F[Findex++] = (this->virtual_chemical_oxidizer_used - this->chemical_oxidizer_used)
                        * this->myUniverse->continuity_constraint_scale_factors(6);

                    //Step 3.3: derivatives
                    if (needG)
                    {
                        size_t Gindex = this->Gindex_dVirtualChemicalOxidizerConstraint_dVirtualChemicalOxidizer;
                        size_t Xindex = this->jGvar->operator[](Gindex);

                        G[Gindex] = this->X_scale_factors->operator[](Xindex)
                            * this->myUniverse->continuity_constraint_scale_factors(6);
                    }
                }//end virtual chemical oxidizer
            }//end process_virtual_propellant_constraints

			void BoundaryEventBase::compute_orbit_elements_after_event(const bool & generate_derivatives)
			{

                // for each frame that we need, we are going to compute the OE in that frame and store them
                for (ReferenceFrame frame : this->required_orbit_element_frames)
                {
                    static EMTG::math::Matrix <doubleType> rotated_rv(6, 1, 0.0);

                    if (frame != ReferenceFrame::ICRF)
                    {
                        // rotate the boundary state into the requested frame
                        static EMTG::math::Matrix <doubleType> temp_matrix_in(3, 1, 0.0);
                        static EMTG::math::Matrix <doubleType> temp_matrix_out(3, 1, 0.0);
                        static EMTG::math::Matrix <doubleType> rotation_matrix(3, 3, 0.0);
                        static EMTG::math::Matrix <doubleType> augmented_rotation_matrix(6, 6, 0.0);
                        
                        // rotate position
                        for (size_t k = 0; k < 3; ++k)
                        {
                            temp_matrix_in(k, 0) = this->state_after_event(k);
                        }
                        this->myUniverse->LocalFrame.rotate_frame_to_frame(EMTG::ReferenceFrame::ICRF,
                            temp_matrix_in,
                            frame,
                            temp_matrix_out,
                            rotation_matrix,
                            this->state_after_event(7)); // current epoch is entry slot 7 in the state vector
                        for (size_t k = 0; k < 3; ++k)
                        {
                            rotated_rv(k) = temp_matrix_out(k, 0);
                        }

                        // rotate velocity
                        for (size_t k = 0; k < 3; ++k)
                        {
                            temp_matrix_in(k, 0) = this->state_after_event(k + 3);
                        }
                        this->myUniverse->LocalFrame.rotate_frame_to_frame(EMTG::ReferenceFrame::ICRF,
                            temp_matrix_in,
                            frame,
                            temp_matrix_out,
                            rotation_matrix,
                            this->state_after_event(7)); // current epoch is entry slot 7 in the state vector
                        for (size_t k = 0; k < 3; ++k)
                        {
                            rotated_rv(k + 3) = temp_matrix_out(k, 0);
                        }

                        // compute the orbital elements from the boundary state
                        EMTG::Astrodynamics::inertial2COE(rotated_rv,
                            this->myUniverse->mu,
                            this->orbit_elements_after_event[frame],
                            generate_derivatives,
                            this->orbit_elements_after_event_partials_wrt_cartesian_state[frame]);

                        // these derivatives come out of inertial2COE in the wrong orientation
                        this->orbit_elements_after_event_partials_wrt_cartesian_state[frame] = this->orbit_elements_after_event_partials_wrt_cartesian_state[frame].transpose();

                        // we now need to rotate the OE Jacobian info into ICRF
                        // make a 6x6 rotation matrix
                        augmented_rotation_matrix.assign_zeros();
                        for (size_t row = 0; row < 3; ++row)
                        {
                            for (size_t col = 0; col < 3; ++col)
                            {
                                augmented_rotation_matrix(row, col) = rotation_matrix(row, col);
                                augmented_rotation_matrix(row + 3, col + 3) = rotation_matrix(row, col);
                            }
                        }
                        this->orbit_elements_after_event_partials_wrt_cartesian_state[frame] *= augmented_rotation_matrix;
                    }
                    else
                    {
                        for (size_t k = 0; k < 6; ++k)
                        {
                            rotated_rv(k) = this->state_after_event(k);
                        }

                        // compute the orbital elements from the boundary state
                        EMTG::Astrodynamics::inertial2COE(rotated_rv,
                            this->myUniverse->mu,
                            this->orbit_elements_after_event[frame],
                            generate_derivatives,
                            this->orbit_elements_after_event_partials_wrt_cartesian_state[frame]);

                        this->orbit_elements_after_event_partials_wrt_cartesian_state[frame] = this->orbit_elements_after_event_partials_wrt_cartesian_state[frame].transpose();
                    }   
                }
			}

            void BoundaryEventBase::process_specialized_constraints(const std::vector<doubleType>& X,
                    size_t& Xindex,
                    std::vector<doubleType>& F,
                    size_t& Findex,
                    std::vector<double>& G,
                    const bool& needG)
            {

				if (this->compute_orbit_elements)
				{
					this->compute_orbit_elements_after_event(needG);
				}

				for (size_t constraintIndex = 0; constraintIndex < this->mySpecializedConstraints.size(); ++constraintIndex)
				{
					this->mySpecializedConstraints[constraintIndex].process_constraint(X, Xindex, F, Findex, G, needG);
				}
            }//end process_specialized_constraints

        void BoundaryEventBase::output_specialized_constraints(std::ofstream& outputfile)
        {
            for (size_t constraintIndex = 0; constraintIndex < this->mySpecializedConstraints.size(); ++constraintIndex)
                this->mySpecializedConstraints[constraintIndex].output(outputfile);
        }

        void BoundaryEventBase::output_ephemeris(std::ofstream& outputfile)
        {
            //start by copying the state before the event into output_state
            math::Matrix<doubleType> output_state = this->state_before_event;

            //6-state may need to be converted to sun-centered
            if (this->myUniverse->central_body.spice_ID != this->myOptions->forward_integrated_ephemeris_central_body_SPICE_ID)
            {
                double LT_dump;
                double bodyStateDouble[6];
                spkez_c(this->myUniverse->central_body.spice_ID, output_state(7)_GETVALUE - (51544.5 * 86400.0), "J2000", "NONE", this->myOptions->forward_integrated_ephemeris_central_body_SPICE_ID, bodyStateDouble, &LT_dump);

                for (size_t stateIndex = 0; stateIndex < 6; ++stateIndex)
                    output_state(stateIndex) += bodyStateDouble[stateIndex];
            }
            
            this->write_ephemeris_line(outputfile,
                output_state,
                math::Matrix<doubleType>(3, 1, 0.0),//control vector
                0.0,
                0.0,
                0.0,
                0,
                0.0,
                "none");
        }//end output_ephemeris()


        void BoundaryEventBase::print_state_before_event(std::ofstream& outputfile)
        {
            outputfile << this->prefix << "state before event:" << std::endl;

            for (size_t stateIndex = 0; stateIndex < 8; ++stateIndex)
            {
                outputfile.width(20);
                outputfile << this->state_before_event(stateIndex);
            }
            outputfile << std::endl;
        }//end print_state_before_event()

        void BoundaryEventBase::print_state_after_event(std::ofstream& outputfile)
        {
            outputfile << this->prefix << "state after event:" << std::endl;

            for (size_t stateIndex = 0; stateIndex < 8; ++stateIndex)
            {
                outputfile.width(20);
                outputfile << this->state_after_event(stateIndex);
            }
            outputfile << std::endl;
        }//end print_state_after_event()
    }//close namespace events
}//close namespace EMTG