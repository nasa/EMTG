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

//options class for propulsion systems
//Jacob Englander 12-30-2016

#ifndef EMTG_PROPULSIONSYSTEMOPTIONS
#define EMTG_PROPULSIONSYSTEMOPTIONS

#include <string>
#include <vector>

#include "EMTG_enums.h"
#include "file_utilities.h"

namespace EMTG
{
    namespace HardwareModels
    {
        class PropulsionSystemOptions
        {
        public:
            //constructor
            PropulsionSystemOptions();
            PropulsionSystemOptions(const std::string& linestring);

            //destructor
            ~PropulsionSystemOptions() {};

            //methods
            void parse_input_line(const std::string& linestring);
            void write_output_line(std::ofstream& outputfile);
            
            //set functions

            //get
            inline std::string getName() const { return this->name; }
            inline SpacecraftThrusterMode getThrusterMode() const { return this->ThrusterMode; }
            inline double getConstantThrust() const { return this->ConstantThrust; }
            inline double getConstantIsp() const { return this->ConstantIsp; }
            inline double getMinimumOrMonopropIsp() const { return this->MinimumOrMonopropIsp; }
            inline double getFixedEfficiency() const { return this->FixedEfficiency; }
            inline double getPmin() const { return this->Pmin; }
            inline double getPmax() const { return this->Pmax; }
            inline std::string getThrottleTableFile() const {return this->ThrottleTableFile; }
            inline std::vector<double> getThrustCoefficients() const { return this->ThrustCoefficients; }
            inline std::vector<double> getMassFlowCoefficients() const { return this->MassFlowCoefficients; }
            inline double getMassPerString() const { return this->MassPerString; }
            inline size_t getNumberOfStrings() const { return this->NumberOfStrings; }
            inline double getMixtureRatio() const { return this->MixtureRatio; }
            inline double getg0() const { return this->g0; }
            inline double getThrustScaleFactor() const {return this->ThrustScaleFactor; }
			inline double getSharpness() const {return this->Sharpness; }

            inline double getThrustCoefficient(const size_t& p) const { return this->ThrustCoefficients[p]; }
            inline double getMassFlowCoefficient(const size_t& p) const { return this->MassFlowCoefficients[p]; }

            void setName(const std::string& input);
            void setNumberOfStrings(const size_t& input);
            void setThrustScaleFactor(const double& input);
            void setThrusterMode(const SpacecraftThrusterMode& input);
            void setConstantThrust(const double& input);
            void setConstantIsp(const double& input);
            void setMinimumOrMonopropIsp(const double& input);
            void setFixedEfficiency(const double& input);
            void setPmin(const double& input);
            void setPmax(const double& input);
            void setThrottleTableFile(const std::string& input);
            void setThrustCoefficients(const std::vector<double>& input);
            void setMassFlowCoefficients(const std::vector<double>& input);
            void setMassPerString(const double& input);
            void setMixtureRatio(const double& input);
            void setg0(const double& input);
			void setSharpness(const double& input);

            void setThrustCoefficient(const size_t& p, const double& input);
            void setMassFlowCoefficient(const size_t& p, const double& input);

            //fields
        protected:
            std::string name;
            SpacecraftThrusterMode ThrusterMode;
            double Pmin;
            double Pmax;
            std::vector<double> ThrustCoefficients;
            std::vector<double> MassFlowCoefficients;
            std::string ThrottleTableFile;
            double ConstantThrust;
            double ConstantIsp;
            double MinimumOrMonopropIsp;
            double FixedEfficiency;
            double MassPerString;
            size_t NumberOfStrings;
            double MixtureRatio;
            double ThrustScaleFactor;
            double g0;
			double Sharpness;
        };

        class PropulsionSystemOptionsLibrary
        {
        public:
            //constructor
            PropulsionSystemOptionsLibrary();
            PropulsionSystemOptionsLibrary(const std::string& filename);
            PropulsionSystemOptionsLibrary(std::ifstream& inputfile);
            PropulsionSystemOptionsLibrary(const std::vector<PropulsionSystemOptions>& PropulsionSystemOptionsVector);

            //destructor
            ~PropulsionSystemOptionsLibrary() {};//destructor does not need to do anything

            //methods
            void parse_input_file(const std::string& filename);
            void parse_input_file(std::ifstream& inputfile);
            void write_output_file(const std::string& filename, const bool& newFile = true);


            //get functions
            PropulsionSystemOptions getPropulsionSystem(const std::string& name);
            SpacecraftThrusterMode getThrusterMode(const std::string& name);
            double getConstantThrust(const std::string& name);
            double getConstantIsp(const std::string& name);
            double getMinimumOrMonopropIsp(const std::string& name);
            double getFixedEfficiency(const std::string& name);
            double getPmin(const std::string& name);
            double getPmax(const std::string& name);
            std::string getThrottleTableFile(const std::string& name);
            std::vector<double> getThrustCoefficients(const std::string& name);
            std::vector<double> getMassFlowCoefficients(const std::string& name);
            double getMassPerString(const std::string& name);
            size_t getNumberOfStrings(const std::string& name);
            double getMixtureRatio(const std::string& name);
            double getThrustScaleFactor(const std::string& name);


        private:
            //private methods
            size_t getIndex(const std::string& name);

            //fields
            std::vector<PropulsionSystemOptions> PropulsionSystems;
            std::string currentKey;
            size_t currentIndex;
        };
    }//end namespace HardwareModels
}//end namespace EMTG

#endif //EMTG_PROPULSIONSYSTEMOPTIONS