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

//ParallelShootingStep body-probe-thrust constraint
//Jacob Englander 4-13-2018

#pragma once

#include "ParallelShootingStep_maneuver_constraint.h"

namespace EMTG
{
    namespace Phases
    {

        class ParallelShootingStep_BPT_angle_constraint : public ParallelShootingStep_maneuver_constraint
        {
        public:
            //constructor
            ParallelShootingStep_BPT_angle_constraint() {};
            ParallelShootingStep_BPT_angle_constraint(const std::string& name,
                const size_t& journeyIndex,
                const size_t& phaseIndex,
                const size_t& stepIndex,
                const size_t& subStepIndex,
                const size_t& stageIndex,
                ParallelShootingStep* myStep,
                Astrodynamics::universe* myUniverse,
                HardwareModels::Spacecraft* mySpacecraft,
                missionoptions* myOptions,
                const std::string& ConstraintDefinition);

            //clone
            virtual ParallelShootingStep_BPT_angle_constraint* clone() const { return new ParallelShootingStep_BPT_angle_constraint(*this); }

            //calcbounds goes in the specialized phase
            void calcbounds();

            //process goes in the specialized phase
            void process_constraint(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

        protected:
            //fields
            Astrodynamics::body* myBody;

            bool isDefinedRelativeToCentralBody;
            math::Matrix<doubleType> spacecraft_to_reference;
            math::Matrix<double> dspacecraft_to_reference_depoch;
            double MinimumAngle_a, MinimumAngle_b, MinimumAngle_c, MinimumAngle_d;
            double MaximumAngle_a, MaximumAngle_b, MaximumAngle_c, MaximumAngle_d;
            bool ConstantMinimumAngle;
            bool ConstantMaximumAngle;

            std::vector<size_t> Gindices_BPT_constraintMinimumAngle_wrt_StepLeftPosition;
            std::vector<size_t> Gindices_BPT_constraintMinimumAngle_wrt_StepLeftTime;
            std::vector<size_t> Gindices_BPT_constraintMinimumAngle_wrt_control;

            std::vector<size_t> Gindices_BPT_constraintMaximumAngle_wrt_StepLeftPosition;
            std::vector<size_t> Gindices_BPT_constraintMaximumAngle_wrt_StepLeftTime;
            std::vector<size_t> Gindices_BPT_constraintMaximumAngle_wrt_control;

            std::vector< std::vector<size_t> > dIndex_MinimumAngle_constraint_wrt_StepLeftPosition; //varIndex, dIndex
            std::vector< std::vector<size_t> > dIndex_MaximumAngle_constraint_wrt_StepLeftPosition; //varIndex, dIndex

            doubleType cosAngleThrustVectorToSun;
        };
    }//close namespace Phases
}//close namespace EMTG