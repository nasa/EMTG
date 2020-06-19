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

//! BoundaryEventBase is the abstract base class from which all boundary events are derived. The base class defines the interfaces common to all boundary events and also the common fields, including the state before and after the event.

#pragma once

#include "doubleType.h"

#include <map>
#include <unordered_set>
#include <string>
#include <vector>

#include <EMTG_Matrix.h>
#include "Spacecraft.h"
#include "missionoptions.h"
#include "universe.h"
#include "writey_thing.h"
#include "sparsey_thing.h"
#include "SpecializedBoundaryConstraintBase.h"

#include "boost/ptr_container/ptr_vector.hpp"

namespace EMTG
{
    namespace BoundaryEvents
    {
        //forward declarations
        namespace SpecializedConstraints
        {
            class SpecializedBoundaryConstraintBase;
        }

        class BoundaryEventBase : virtual public writey_thing, virtual public sparsey_thing
        {
        public:
            //!default constructor
            BoundaryEventBase();

            //!constructor that calls initialize
            BoundaryEventBase(const std::string& name,
                                const size_t& journeyIndex,
                                const size_t& phaseIndex,
                                size_t& stageIndex,
                                Astrodynamics::universe* Universe,
                                HardwareModels::Spacecraft* mySpacecraft,
                                missionoptions* myOptions);

            //!method to initialize properties of the base class, intended to be called from EMTG::Journey
            void initialize(const std::string& name,
                            const size_t& journeyIndex,
                            const size_t& phaseIndex,
                            size_t& stageIndex,
                            Astrodynamics::universe* Universe,
                            HardwareModels::Spacecraft* mySpacecraft,
                            missionoptions* myOptions);

            //!pure abstract method to write a line in the .emtg output file
            virtual void output(std::ofstream& outputfile,
                                const double& launchdate,
                                size_t& eventcount) = 0;

            //method to construct boundary constraints - gets implemented by arrival and departure. This is public so that we can fire it from other places
            virtual void construct_boundary_constraints(std::vector<std::string> givenConstraints = {}) = 0;

            //!method to write a line in the .ephemeris file, containing user-defined fields
            virtual void output_ephemeris(std::ofstream& outputfile);

            //!method to output the values of all specialized constraints held by this boundary event. This method just calls the output method of each specialized constraint attached to the boundary.
            virtual void output_specialized_constraints(std::ofstream& outputfile);

            //!method to print state before event in ICRF
            void print_state_before_event(std::ofstream& outputfile);

            //!method to print state after event in ICRF
            void print_state_after_event(std::ofstream& outputfile);

            //get
            inline math::Matrix<doubleType>& get_state_before_event() { return this->state_before_event; }
            inline math::Matrix<doubleType>& get_state_after_event() { return this->state_after_event; }
            inline math::Matrix<doubleType> get_boundary_state() const { return this->boundary_state; }
			inline math::Matrix<doubleType> get_orbit_elements_after_event(const ReferenceFrame & reference_frame) { return this->orbit_elements_after_event[reference_frame]; }
			inline math::Matrix<doubleType> get_orbit_element_Jacobian_after_event(const ReferenceFrame & reference_frame) { return this->orbit_elements_after_event_partials_wrt_cartesian_state[reference_frame]; }
			inline bool getComputeOrbitElements() const { return this->compute_orbit_elements; }
            inline std::vector< std::tuple<size_t, size_t, double> >& get_Derivatives_of_StateBeforeEvent() { return this->Derivatives_of_StateBeforeEvent; }
            inline std::vector< std::tuple<size_t, size_t, double> >& get_Derivatives_of_StateAfterEvent() { return this->Derivatives_of_StateAfterEvent; }
            inline std::vector< std::tuple<size_t, size_t, double> >& get_Derivatives_of_StateBeforeEvent_wrt_Time() { return this->Derivatives_of_StateBeforeEvent_wrt_Time; }
            inline std::vector< std::tuple<size_t, size_t, double> >& get_Derivatives_of_StateAfterEvent_wrt_Time() { return this->Derivatives_of_StateAfterEvent_wrt_Time; }
            inline size_t getX_index_of_first_decision_variable_in_this_event() const { return this->X_index_of_first_decision_variable_in_this_event; }
            
            inline double get_StatisticalDeltav() const { return this->TCM_magnitude; }
            inline doubleType get_DeterministicDeltav() const { return this->EventDeterministicDeltav; }
            inline bool getLeftBoundaryIsABody() const { return this->LeftBoundaryIsABody; }
            std::string getName() const { return this->name; }
            inline size_t get_stageIndex() const { return this->stageIndex; }
            inline doubleType getC3() const { return this->C3; }
            inline const std::vector<size_t>& get_Xindices_EventLeftEpoch() const { return this->Xindices_EventLeftEpoch; }
            inline const std::vector<size_t>& get_Xindices_EventRightEpoch() const { return this->Xindices_EventRightEpoch; }
            inline const bool& get_hasElectricManeuver() const { return this->hasElectricManeuver; }
            inline const bool& get_hasBipropManeuver() const { return this->hasBipropManeuver; }
            inline const bool& get_hasMonopropManeuver() const { return this->hasMonopropManeuver; }
            inline const bool& get_hasManeuver() const { return this->hasManeuver; }

            inline Astrodynamics::body* getBody() { return this->myBody; }

            //virtual dummy get methods, to be overriden in derived classes that actually use them
            inline size_t getXindex_of_initial_impulse() const { return -1; }
            virtual std::vector<size_t> getdIndex_VbeforeEvent_dVinfinity_in() const { return std::vector<size_t>(); }
            virtual math::Matrix<doubleType> getVinfinityIn() const { return math::Matrix<doubleType>(); }

            //set
            void setName(const std::string& name) { this->name = name; }
			void setComputeOrbitElements(const bool & flag) { this->compute_orbit_elements = flag; }
            void setJourneyOptionsPointer(JourneyOptions* myJourneyOptions) { this->myJourneyOptions = myJourneyOptions; }

            //!method to add an orbit element reference frame to this boundary event.
            /*!
                When the user poses an orbit element constraint to EMTG, the boundary event needs to compute the orbit elements in the frame of the constraint and also their partial derivatives with respect to the decision variables. This method's job is to create a list of frames for which orbit elements are needed.
            */
            void add_orbit_element_reference_frame(const ReferenceFrame & constraint_reference_frame) 
            { 
                this->required_orbit_element_frames.insert(constraint_reference_frame); 

                // scan through both the orbit element set and the orbit element Jacobian set to make sure they are initialized with the appropriate containers
                bool frame_in_element_set = false;
                for (auto entry : this->orbit_elements_after_event)
                {
                    if (entry.first == constraint_reference_frame)
                    {
                        frame_in_element_set = true;
                    }
                }

                if (!frame_in_element_set)
                {
                    this->orbit_elements_after_event.insert(std::pair<ReferenceFrame, math::Matrix<doubleType>>(constraint_reference_frame, math::Matrix<doubleType>(6, 1, 0.0)));
                    this->orbit_elements_after_event_partials_wrt_cartesian_state.insert(std::pair<ReferenceFrame, math::Matrix<doubleType>>(constraint_reference_frame, math::Matrix<doubleType>(6, 6, 0.0)));
                }
            }

            //!method to initialize the various data pointers that will be used by each derived class's calcbounds() methods.
            virtual void setup_calcbounds(
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
                std::vector<double>* A);

            //!abstract prototype for the derived class's calcbounds() method.
            virtual void calcbounds(std::vector<size_t> timeVariables) = 0;

            //!abstract prototype for the derived class's process_event() method
            virtual void process_event(const std::vector<doubleType>& X,
                                       size_t& Xindex,
                                       std::vector<doubleType>& F,
                                       size_t& Findex,
                                       std::vector<double>& G,
                                       const bool& needG) = 0;

            //!method to reset the Event Transition Matrix (ETM). We reset the ETM in each iteration so that we know that the partial derivative information is pristine.
            inline void reset_ETM() { this->ETM = math::Matrix<double>(8, math::identity); }

            //!method to compute orbit elements at the right-hand side of the boundary event. This method is called by specialized orbit element constraints that are children of the boundary event
			void compute_orbit_elements_after_event(const bool & generate_derivatives);

        protected:
            //!abstract prototype of calcbounds_event_left_side()
            virtual void calcbounds_event_left_side(std::vector<size_t> timeVariables) = 0;

            //!method to find all of the time variables that precede this event and put their Xindex into a helper array. That way we can easily compute the epoch at the start of the event.
            void calculate_dependencies_left_epoch(std::vector<size_t> timeVariables);

            //!abstract prototype of calcbounds_event_right_side()
            virtual void calcbounds_event_right_side() = 0;

            //!abstract prototype of calcbounds_event_main()
            virtual void calcbounds_event_main() = 0;

            //!abstract prototype of calcbounds_virtual_propellant_constraints()
            virtual void calcbounds_virtual_propellant_constraints() = 0;

            //!abstract prototype of calcbounds_deltav_contribution()
            virtual void calcbounds_deltav_contribution() {};

            //!abstract prototype of calcbounds_specialized_constraints()
            virtual void calcbounds_specialized_constraints();

            //!This method gets called at the beginning of each process_event() call and ensures that the spacecraft stage is set correctly for this boundary event.
            void process_staging();

            //!abstract prototype for process_event_main()
            virtual void process_event_main(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG) = 0;

            //!abstract prototype for process_event_left_side()
            virtual void process_event_left_side(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG) = 0;

            //!abstract prototype for process_left_epoch()
            void process_left_epoch(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            //!abstract prototype for process_event_right_side()
            virtual void process_event_right_side(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG) = 0;

            //!abstract prototype for process_virtual_propellant_constraints()
            virtual void process_virtual_propellant_constraints(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG) = 0;

            //!virtual prototype for process_deltav_contribution(). This prototype is deliberately not abstract because not all boundary events have a delta-v contribution.
            virtual void process_deltav_contribution(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG) {};

            //!method to evaluate all of the boundary event's specialized constraints. These are stored in a vector of constraint objects at construct-time and then evaluated when process_event() calls process_specialized_constraints()
            virtual void process_specialized_constraints(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);


            //fields
            std::string name;
            size_t journeyIndex;
            size_t phaseIndex;
            size_t stageIndex;
            bool isFirstEventInMission;
            Astrodynamics::body* myBody;
            HardwareModels::Spacecraft* mySpacecraft;
            JourneyOptions* myJourneyOptions;
            size_t X_index_of_first_decision_variable_in_this_event;
            size_t F_index_of_first_constraint_in_this_event;
            bool LeftBoundaryIsABody;
            PropagatorType myPropagatorType;
			bool compute_orbit_elements;

            //times
            doubleType EventTimeWidth;
            doubleType EventLeftEpoch;
            doubleType EventRightEpoch;
            bool EventHasTimeWidth;
            std::vector<size_t> Xindices_EventLeftEpoch;
            std::vector<size_t> Xindices_EventRightEpoch;
            
            //states
            std::vector<std::string> stateNames;
            std::vector<std::string> COENames;
            math::Matrix<doubleType> state_before_event;
            math::Matrix<doubleType> state_after_event;
            math::Matrix<doubleType> state_after_event_propagated;
            math::Matrix<doubleType> boundary_state;
            std::map<ReferenceFrame, math::Matrix<doubleType>> orbit_elements_after_event;
            std::map<ReferenceFrame, math::Matrix<doubleType>> orbit_elements_after_event_partials_wrt_cartesian_state;
            std::unordered_set<ReferenceFrame> required_orbit_element_frames;
            std::vector< std::tuple<size_t, size_t, double> > Derivatives_of_StateBeforeEvent;//Xindex, stateIndex, derivative value
            std::vector< std::tuple<size_t, size_t, double> > Derivatives_of_StateAfterEvent;//Xindex, stateIndex, derivative value
            std::vector< std::tuple<size_t, size_t, double> > Derivatives_of_StateBeforeEvent_wrt_Time;//Xindex, stateIndex, derivative value
            std::vector< std::tuple<size_t, size_t, double> > Derivatives_of_StateAfterEvent_wrt_Time;//Xindex, stateIndex, derivative value
            doubleType C3;

            size_t dIndex_mass_wrt_encodedMass;

            //event transition matrix
            math::Matrix<double> ETM; //8 x 8, last row/column is epoch

            //*********************delta-v and derivatives
            doubleType EventDeterministicDeltav;
            
            //staging
            bool StageBeforeEvent; //if true, then we need to apply a staging mass constraint
      
            //propellant used
            doubleType virtual_electric_propellant_used;
            doubleType virtual_chemical_fuel_used;
            doubleType virtual_chemical_oxidizer_used;
            doubleType electric_propellant_used;
            doubleType chemical_fuel_used;
            doubleType chemical_oxidizer_used;
            bool hasElectricManeuver;
            bool hasMonopropManeuver;
            bool hasBipropManeuver;
            bool hasManeuver;
            PropulsionSystemChoice ChemicalManeuverType;

            size_t Findex_VirtualElectricPropellantConstraint;
            size_t Findex_VirtualChemicalFuelConstraint;
            size_t Findex_VirtualChemicalOxidizerConstraint;

            size_t Gindex_dVirtualElectricPropellantConstraint_dVirtualElectricPropellant;
            size_t Gindex_dVirtualChemicalFuelConstraint_dVirtualChemicalFuel;
            size_t Gindex_dVirtualChemicalOxidizerConstraint_dVirtualChemicalOxidizer;

            //TCM stuff
            bool hasTCM;
            double TCM_magnitude;

            //vector of specialized constraint objects
            boost::ptr_vector< SpecializedConstraints::SpecializedBoundaryConstraintBase > mySpecializedConstraints;
        };
    }//close namespace events
}//close namespace EMTG