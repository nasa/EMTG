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

//header for Edelbaum spiral
//Jacob Englander 10-24-2017
#pragma once

#include "doubleType.h"

#include "Spacecraft.h"
#include "EMTG_Matrix.h"
#include "universe.h"
#include "missionoptions.h"
#include "EdelbaumSpiralSegment.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        //forward declare EdelbaumSpiralSegment
        class EdelbaumSpiralSegment;

        class EdelbaumSpiral
        {
        public:
            EdelbaumSpiral() {};
            EdelbaumSpiral(const std::string& name,
                const size_t& journeyIndex,
                const size_t& phaseIndex,
                size_t& stageIndex,
                Astrodynamics::universe* Universe,
                Astrodynamics::body* myBody,
                HardwareModels::Spacecraft* mySpacecraft,
                missionoptions* myOptions);
            
            //setup calcbounds
            void setup_calcbounds(
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
            void output(std::ofstream& outputfile,
                const double& launchdate,
                size_t& eventcount);

            //segments
            void initialize_spiral_segments();
            
            void calcbounds(std::vector<size_t> Xindices_leftEpoch);

            void process_spiral(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);



            inline void setSpiralStartRadius(const double& SpiralStartRadius) { this->SpiralStartRadius = SpiralStartRadius; }
            inline void setSpiralEndRadius(const double& SpiralEndRadius) { this->SpiralEndRadius = SpiralEndRadius; }
            inline void setDirection(const double& Direction) { this->Direction = Direction; }
            inline void setInitialMass(const doubleType& InitialMass) { this->InitialMass = InitialMass; }
            inline void setInitialEpoch(const doubleType& InitialEpoch) { this->InitialEpoch = InitialEpoch; }
            inline doubleType getInitialMass() const { return this->InitialMass; }
            inline double getTotalDeltav() const { return this->TotalDeltav; }
            inline double getDirection() const { return this->Direction; }
            inline math::Matrix<doubleType> getStateAfterSpiral() const { return this->state_after_spiral; }
            std::vector<size_t> get_Xindices_SpiralEndEpoch() const { return this->Xindices_SpiralEndEpoch; }

            inline std::vector< std::tuple<size_t, size_t, double> >& get_Derivatives_of_StateAfterSpiral() { return this->Derivatives_of_StateAfterSpiral; }
            inline std::vector< std::tuple<size_t, size_t, double> >& get_Derivatives_of_StateAfterSpiral_wrt_Time() { return this->Derivatives_of_StateAfterSpiral_wrt_Time; }

        protected:
            virtual void process_spiral_end_epoch(const std::vector<doubleType>& X,
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
            Astrodynamics::universe* myUniverse;
            Astrodynamics::body* myBody;
            HardwareModels::Spacecraft* mySpacecraft;
            missionoptions* myOptions;
            JourneyOptions* myJourneyOptions;
            size_t X_index_of_first_decision_variable_in_this_event;

            //calcbounds fields
            std::string prefix;
            std::vector<double>* Xupperbounds;
            std::vector<double>* Xlowerbounds;
            std::vector<double>* X_scale_factors;
            std::vector<double>* Fupperbounds;
            std::vector<double>* Flowerbounds;
            std::vector<double>* F_scale_factors;
            std::vector<std::string>* Xdescriptions;
            std::vector<std::string>* Fdescriptions;
            std::vector<size_t>* iGfun;
            std::vector<size_t>* jGvar;
            std::vector<std::string>* Gdescriptions;
            std::vector<size_t>* iAfun;
            std::vector<size_t>* jAvar;
            std::vector<std::string>* Adescriptions;
            std::vector<double>* A;

            //spiral things
            double SpiralStartRadius;
            double SpiralEndRadius;
            double Direction; //either 1.0 or -1.0, depending on departure vs arrival
            double TotalDeltav;
            doubleType InitialMass;
            doubleType InitialEpoch;
            std::vector<EdelbaumSpiralSegment> SpiralSegments;
            std::vector<double> spiralRadius;            
            math::Matrix<doubleType> state_after_spiral;
            std::vector< std::tuple<size_t, size_t, double> > Derivatives_of_StateAfterSpiral;//Xindex, stateIndex, derivative value
            std::vector< std::tuple<size_t, size_t, double> > Derivatives_of_StateAfterSpiral_wrt_Time;//Xindex, stateIndex, derivative value

            //other things
            std::vector< std::size_t > Xindices_SpiralEndEpoch;
            doubleType LaunchDate;
        };//end class EdelbaumSpiral
    }//end namespace BoundaryEvents
}//end namespace EMTG
