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

            //destructor
            virtual ~ParallelShootingStep() {};

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

            inline std::vector<size_t>& getXindices_control(const size_t& subStepIndex) { return this->Xindices_control[subStepIndex]; }
            inline math::Matrix<doubleType>& getControlVector(const size_t& subStepIndex) { return this->ControlVector[subStepIndex]; }
            inline size_t& getXindex_rMag() { return this->Xindex_rMag; }
            inline size_t& getXindex_RA() { return this->Xindex_RA; }
            inline size_t& getXindex_DEC() { return this->Xindex_DEC; }

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

            virtual void calcbounds_waypoint_tracking();

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
            
            virtual void process_waypoint_tracking(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

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
            std::vector< std::vector<size_t> > Gindices_LeftMatchPoint_wrt_PreviousStepRightState; //[varIndex][entryIndex]
            std::vector<size_t> ListOfTimeVariablesAffectingPreviousStepRightState;
            std::vector< std::vector< std::tuple<size_t, size_t> > > DerivativesOfPreviousStepRightStateByTimeVariable; //[varIndex][entryIndex]{stateIndex, dIndex}
            std::vector< std::vector<size_t> > Gindices_LeftMatchPoint_wrt_PreviousStepRightStateTime; //[varIndex][entryIndex]

            std::vector<size_t> Gindices_StepLeftMatchPoint_wrt_r;
            std::vector<size_t> Gindices_StepLeftMatchPoint_wrt_RA;
            std::vector<size_t> Gindices_StepLeftMatchPoint_wrt_DEC;
            std::vector<size_t> Gindices_StepLeftMatchPoint_wrt_v;
            std::vector<size_t> Gindices_StepLeftMatchPoint_wrt_vRA;
            std::vector<size_t> Gindices_StepLeftMatchPoint_wrt_vDEC;
            std::vector<size_t> Gindices_StepLeftMatchPoint_wrt_AZ;
            std::vector<size_t> Gindices_StepLeftMatchPoint_wrt_FPA;
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

            //indices
            std::vector<size_t> Xindices_EventLeftEpoch;
            std::vector< std::vector<size_t> > Xindices_control; //substep, control index

            size_t Xindex_rMag;
            size_t Xindex_RA;
            size_t Xindex_DEC;
            size_t Xindex_vMag;
            size_t Xindex_vRA;
            size_t Xindex_vDEC;
            size_t Xindex_AZ;
            size_t Xindex_FPA;
            size_t Xindex_mass;
            size_t Xindex_chemical_fuel;
            size_t Xindex_electric_propellant;
            size_t Xindex_PhaseFlightTime;

            size_t dIndex_mass_wrt_mass;
            size_t dIndex_chemicalFuel_wrt_chemicalFuel;
            size_t dIndex_electricPropellant_wrt_electricPropellant;
            std::vector<size_t> dIndex_StateStepLeftInertial_wrt_r;//x, y, z
            std::vector<size_t> dIndex_StateStepLeftInertial_wrt_RA;//x, y
            std::vector<size_t> dIndex_StateStepLeftInertial_wrt_DEC;//x, y, z
            std::vector<size_t> dIndex_StateStepLeftInertial_wrt_v;//xdot, ydot, zdot
            std::vector<size_t> dIndex_StateStepLeftInertial_wrt_vRA;//xdot, ydot 
            std::vector<size_t> dIndex_StateStepLeftInertial_wrt_vDEC;//xdot, ydot, zdot
            std::vector<size_t> dIndex_StateStepLeftInertial_wrt_AZ;//xdot, ydot, zdot 
            std::vector<size_t> dIndex_StateStepLeftInertial_wrt_FPA;//xdot, ydot, zdot 
            std::vector< std::vector<size_t> > dIndex_StateStepRightInertial_wrt_LeftStateVariables;//right stateIndex, variable
            std::vector< std::vector<size_t> > dIndex_StateStepRightInertial_wrt_PreviousTimeVariables;//right stateIndex, variable
            std::vector<size_t> dIndex_StateStepRightInertial_wrt_PhaseFlightTime;


            std::vector< std::vector<size_t> > TruthTable_StateStepRightInertial_wrt_StateStepLeftInertial;
            std::vector<size_t> TruthTable_StateStepRightInertial_wrt_Control;
            std::vector<size_t> TruthTable_StateStepRightInertial_wrt_PhaseFlightTime;
            
            //waypoint tracking
            math::Matrix<doubleType> WaypointError;
            math::Matrix<doubleType> CovarianceInverse;
            math::Matrix<double> CovarianceInverseDouble;
            math::Matrix<double> dCovarianceInverse_depoch;
            doubleType mahalanobis_distance;
            doubleType virtual_mahalanobis_distance;
            size_t Gindex_virtual_mahalanobis_distance_wrt_virtual_mahalanobis_distance;
            std::vector<size_t> Gindex_virtual_mahalanobis_distance_wrt_Time;
            size_t Gindex_virtual_mahalanobis_distance_wrt_rMag;
            size_t Gindex_virtual_mahalanobis_distance_wrt_RA;
            size_t Gindex_virtual_mahalanobis_distance_wrt_DEC;
            size_t Gindex_virtual_mahalanobis_distance_wrt_vMag;
            size_t Gindex_virtual_mahalanobis_distance_wrt_vRA;
            size_t Gindex_virtual_mahalanobis_distance_wrt_vDEC;
            size_t Gindex_virtual_mahalanobis_distance_wrt_AZ;
            size_t Gindex_virtual_mahalanobis_distance_wrt_FPA;
            size_t Gindex_virtual_mahalanobis_distance_wrt_mass;
        };

        inline ParallelShootingStep * new_clone(ParallelShootingStep const & other)
        {
            return other.clone();
        }
    }//close namespace Phases
}//close namespace EMTG