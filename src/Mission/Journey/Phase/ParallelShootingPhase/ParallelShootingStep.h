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

//base parallel shooting step for EMTGv9
//Jacob Englander 2-20-2018

#pragma once

#include "doubleType.h"

#include "Spacecraft.h"
#include "EMTG_Matrix.h"
#include "universe.h"
#include "missionoptions.h"
#include "PropagatorBase.h"

#include "ParallelShootingPhase.h"
#include "ParallelShootingStepDistanceConstraint.h"
#include "ParallelShootingStep_maneuver_constraint.h"

#include "StateRepresentation.h"

#include "sparsey_thing.h"
#include "writey_thing.h"

#include "boost/ptr_container/ptr_vector.hpp"

namespace EMTG
{
    namespace Phases
    {
        //forward declare ParallelShootingPhase
        class ParallelShootingPhase;

        //forward declare ParallelShootingStepDistanceConstraint
        class ParallelShootingStepDistanceConstraint;

        //forward declare ParallelShootingStep_maneuver_constraint
        class ParallelShootingStep_maneuver_constraint;

        class ParallelShootingStep : virtual public writey_thing, virtual public sparsey_thing
        {
        public:
            //constructor
            ParallelShootingStep();
            ParallelShootingStep(const std::string& name,
                const size_t& journeyIndex,
                const size_t& phaseIndex,
                const size_t& stepIndex,
                const size_t& stageIndex,
                ParallelShootingStep* previousStep,
                ParallelShootingPhase* myPhase,
                Astrodynamics::universe* myUniverse,
                HardwareModels::Spacecraft* mySpacecraft,
                missionoptions* myOptions);

            virtual void initialize(const std::string& name,
                const size_t& journeyIndex,
                const size_t& phaseIndex,
                const size_t& stepIndex,
                const size_t& stageIndex,
                ParallelShootingStep* previousStep,
                ParallelShootingPhase* myPhase,
                Astrodynamics::universe* myUniverse,
                HardwareModels::Spacecraft* mySpacecraft,
                missionoptions* myOptions);

            //destructor
            virtual ~ParallelShootingStep();

            //clone
            virtual ParallelShootingStep* clone() const = 0;

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

            //output
            virtual void output(std::ofstream& outputfile,
                size_t& eventcount) = 0;

            virtual void output_ephemeris(std::ofstream& outputfile, std::ofstream & acceleration_model_file) = 0;
            
            virtual void output_STM(size_t& STMindex);

            virtual void output_maneuver_and_target_spec(std::ofstream& maneuver_spec_file, std::ofstream& target_spec_file, bool& haveManeuverNeedTarget) = 0;

            virtual void calcbounds_step() = 0;

            virtual void process_step(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG) = 0;

            //get/set
            inline std::string getName() const { return this->name; }
            inline doubleType getstepDeltav() const { return this->stepDeltav; }
            inline double getStepDutyCycle() const { return this->StepDutyCycle; }

            inline math::Matrix<doubleType>& get_StateStepLeftInertial() { return this->StateStepLeftInertial; }
            inline math::Matrix<doubleType>& get_StateStepRightInertial() { return this->StateStepRightInertial; }
            inline std::vector< std::tuple<size_t, size_t, double> >& get_Derivatives_of_StateStepLeftInertial() { return this->Derivatives_of_StateStepLeftInertial; }
            inline std::vector< std::tuple<size_t, size_t, double> >& get_Derivatives_of_StateStepRightInertial() { return this->Derivatives_of_StateStepRightInertial; }
            inline std::vector< std::tuple<size_t, size_t, double> >& get_Derivatives_of_StateStepLeftInertial_wrt_Time() { return this->Derivatives_of_StateStepLeftInertial_wrt_Time; }
            inline std::vector< std::tuple<size_t, size_t, double> >& get_Derivatives_of_StateStepRightInertial_wrt_Time() { return this->Derivatives_of_StateStepRightInertial_wrt_Time; }

            inline std::vector<size_t>& getListOfVariablesAffectingCurrentStepLeftState() { return this->ListOfVariablesAffectingCurrentStepLeftState; }
            inline std::vector<size_t>& getListOfVariablesAffectingCurrentStepLeftPosition() { return this->ListOfVariablesAffectingCurrentStepLeftPosition; }
            inline std::vector<size_t>& getListOfVariablesAffectingCurrentStepLeftVelocity() { return this->ListOfVariablesAffectingCurrentStepLeftVelocity; }
            inline std::vector<size_t>& getListOfTimeVariablesAffectingCurrentStepLeftState() { return this->ListOfTimeVariablesAffectingCurrentStepLeftState; }
            inline std::vector< std::vector< std::tuple<size_t, size_t> > >& getDerivativesOfCurrentStepLeftStateByVariable() { return this->DerivativesOfCurrentStepLeftStateByVariable; }
            inline std::vector< std::vector< std::tuple<size_t, size_t> > >& getDerivativesOfCurrentStepLeftStateByTimeVariable() { return this->DerivativesOfCurrentStepLeftStateByTimeVariable; }

            inline std::vector<size_t>& getXindex_state_elements() { return this->Xindex_state_elements; };
            inline std::vector<size_t>& getXindices_control(const size_t& subStepIndex) { return this->Xindices_control[subStepIndex]; }
            inline math::Matrix<doubleType>& getControlVector(const size_t& subStepIndex) { return this->ControlVector[subStepIndex]; }

        protected:
            //configure propagator
            virtual void configure_propagator() = 0;

            //calcbounds
            void calcbounds_step_left_state();

            void calculate_dependencies_epoch_time();

            void calcbounds_step_control();

            void calcbounds_step_main();

            virtual void calcbounds_step_left_match_point_constraints();

            virtual void calcbounds_step_right_match_point_constraints() {}; //this is blank except for the right-hand step

            void calcbounds_distance_constraints();

            void calcbounds_maneuver_constraints();

            virtual void calcbounds_deltav_contribution() {}; //empty for now

            //process
            virtual void process_epoch_time(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            void process_step_left_state(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            void process_step_control(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            virtual void process_step_main(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG) = 0;

            virtual void process_step_left_match_point_constraints(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            virtual void process_step_right_match_point_constraints(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG) {}; //this is blank except for the right-hand step

            void process_distance_constraints(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            void process_maneuver_constraints(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);
            
            virtual void process_deltav_contribution(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG) {};//empty for now

            void process_derivative_tuples(const std::vector<doubleType>& X,
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
            size_t stepIndex;
            ParallelShootingStep* previousStep;
            HardwareModels::Spacecraft* mySpacecraft;
            JourneyOptions* myJourneyOptions;

            size_t First_X_entry_this_step;

            //pointers to external things
            double dStepTime_dPhaseFlightTime;
            ParallelShootingPhase* myPhase;

            //internal states
            math::Matrix<doubleType> StateStepLeftInertial;
            math::Matrix<doubleType> StateStepLeftEncoded;
            std::vector < math::Matrix<doubleType> > StateAfterSubStepInertial;
            math::Matrix<doubleType> StateStepRightInertial;
            std::vector<std::string> stateNames;

            std::vector< std::tuple<size_t, size_t, double> > Derivatives_of_StateStepLeftInertial; //[dIndex]{Xindex, stateIndex, value}
            std::vector< std::tuple<size_t, size_t, double> > Derivatives_of_StateStepRightInertial; //[dIndex]{Xindex, stateIndex, value}
            std::vector< std::tuple<size_t, size_t, double> > Derivatives_of_StateStepLeftInertial_wrt_Time; //[dIndex]{Xindex, stateIndex, value}
            std::vector< std::tuple<size_t, size_t, double> > Derivatives_of_StateStepRightInertial_wrt_Time; //[dIndex]{Xindex, stateIndex, value}

            //left match point
            size_t numMatchConstraints;
            std::vector<size_t> Findices_left_match_point_constraints;
            std::vector<size_t> ListOfVariablesAffectingCurrentStepLeftState;
            std::vector<size_t> ListOfVariablesAffectingCurrentStepLeftPosition;
            std::vector<size_t> ListOfVariablesAffectingCurrentStepLeftVelocity;
            std::vector< std::vector< std::tuple<size_t, size_t> > > DerivativesOfCurrentStepLeftStateByVariable; //[varIndex][entryIndex]{stateIndex, dIndex}

            std::vector<size_t> ListOfTimeVariablesAffectingCurrentStepLeftState;
            std::vector< std::vector< std::tuple<size_t, size_t> > > DerivativesOfCurrentStepLeftStateByTimeVariable; //[varIndex][entryIndex]{stateIndex, dIndex}

            std::vector<size_t> ListOfVariablesAffectingPreviousStepRightState;
            std::vector< std::vector< std::tuple<size_t, size_t> > > DerivativesOfPreviousStepRightStateByVariable; //[varIndex][entryIndex]{stateIndex, dIndex}            
            std::vector< std::vector<size_t> > Gindices_LeftMatchPoint6state_wrt_PreviousStepRightState; //[varIndex][entryIndex/stateIndex]
            std::vector< std::vector<size_t> > Gindices_LeftMatchPointMassConstraints_wrt_PreviousStepRightState; //[varIndex][entryIndex]

            std::vector<size_t> ListOfTimeVariablesAffectingPreviousStepRightState;
            std::vector< std::vector< std::tuple<size_t, size_t> > > DerivativesOfPreviousStepRightStateByTimeVariable; //[varIndex][entryIndex]{stateIndex, dIndex}
            std::vector< std::vector<size_t> > Gindices_LeftMatchPoint6state_wrt_PreviousStepRightStateTime; //[varIndex][entryIndex/stateIndex]
            std::vector< std::vector<size_t> > Gindices_LeftMatchPointMassConstraints_wrt_PreviousStepRightStateTime; //[varIndex][entryIndex]

            size_t Gindex_StepLeftMatchPoint_wrt_Mass;

            //right state derivatives
            std::vector< std::vector< std::tuple<size_t, size_t> > > DerivativesOfCurrentStepRightStateByVariable; //Xindex, stateIndex, dIndex                                               
            std::vector< std::vector< std::tuple<size_t, size_t> > > DerivativesOfCurrentStepRightStateByTimeVariable; //Xindex, stateIndex, dIndex

            //time
            doubleType PhaseFlightTime;
            doubleType StepFlightTime;
            doubleType StepLeftEpoch;
            doubleType LaunchDate;
            double StepDutyCycle;

            //propagator
            boost::ptr_vector < Astrodynamics::PropagatorBase > myPropagators;
            std::vector< math::Matrix<double> > STM;//substep index
            std::vector< math::Matrix<double> > AugmentedSTM;//substep index
            std::vector< math::Matrix<double> > CumulativeAugmentedSTM;//substep index
            std::vector< math::Matrix<double> > dPropagatedStatedIndependentVariable;//substep index

            //tanks
            doubleType chemical_fuel_used;
            doubleType electric_propellant_used;
            doubleType ACS_fuel_used;
            doubleType virtual_electric_propellant_used;
            doubleType virtual_chemical_fuel_used;

            size_t Findex_VirtualElectricPropellantConstraint;
            size_t Findex_VirtualChemicalFuelConstraint;

            size_t Gindex_StepLeftMatchPoint_wrt_ElectricPropellant;
            size_t Gindex_StepLeftMatchPoint_wrt_ChemicalFuel;

            //control
            bool isVSI;
            size_t num_controls;
            std::vector< std::vector<size_t> > G_indices_control_magnitude_constraint_wrt_Control;
            std::vector< math::Matrix<doubleType> > ControlVector;
            std::vector<doubleType> throttle;
            std::vector< std::vector< std::vector<size_t> > > dIndex_right_state_wrt_control;//substep, state index, control index
            size_t num_interior_control_points;

            //ephemeris resolution
            double EphemerisOutputResolution;

            //output state
            math::Matrix<doubleType> output_state;
            math::Matrix<doubleType> temp_state;

            //delta-v
            doubleType stepDeltav;

            //distance constraints
            std::vector<ParallelShootingStepDistanceConstraint> myDistanceConstraints;

            //maneuver constraints
            boost::ptr_vector<ParallelShootingStep_maneuver_constraint> myManeuverConstraints;

            //state representation
            Astrodynamics::StateRepresentationBase* myStateRepresentation;
            std::vector<bool> encodedStateIsAngle;
            std::vector< std::tuple<std::string, double, double, double, bool> > statesToRepresent; //state name, lower bound, upper bound, scale factor, isAngle

            //scaling
            math::Matrix<double> continuity_constraint_scale_factors;

            //indices
            std::vector<size_t> Xindices_EventLeftEpoch;
            std::vector< std::vector<size_t> > Xindices_control; //substep, control index

            size_t Xindex_mass;
            size_t Xindex_chemical_fuel;
            size_t Xindex_electric_propellant;
            size_t Xindex_PhaseFlightTime;

            std::vector<size_t> Xindex_state_elements; //for example {x, y, z, vx, vy, vz} or {r, RA, DEC, v, vRA, vDEC}

            std::vector<size_t> Gindices_StepLeftMatchPoint_wrt_StepLeftEncodedState;

            size_t dIndex_mass_wrt_mass;
            size_t dIndex_chemicalFuel_wrt_chemicalFuel;
            size_t dIndex_electricPropellant_wrt_electricPropellant;

            std::vector< std::vector<size_t> > dIndex_StateStepLeftInertial_wrt_StateElements;

            std::vector< std::vector<size_t> > dIndex_StateStepRightInertial_wrt_LeftStateVariables;//right stateIndex, variable
            std::vector< std::vector<size_t> > dIndex_StateStepRightInertial_wrt_PreviousTimeVariables;//right stateIndex, variable
            std::vector<size_t> dIndex_StateStepRightInertial_wrt_PhaseFlightTime;
            
            std::vector< std::vector<bool> > TruthTable_StateStepLeftCartesian_wrt_StateStepLeftEncoded; //stateIndex, varIndex
            std::vector< std::vector<bool> > TruthTable_StateStepRightInertial_wrt_StateStepLeftInertial;
            std::vector<bool> TruthTable_StateStepRightInertial_wrt_Control;
            std::vector<bool> TruthTable_StateStepRightInertial_wrt_PhaseFlightTime;
        };

        inline ParallelShootingStep * new_clone(ParallelShootingStep const & other)
        {
            return other.clone();
        }
    }//close namespace Phases
}//close namespace EMTG