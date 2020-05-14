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

//throttle table class for EMTG
//overloaded for use with algorithmic differentiation
//Jacob Englander 5-14-2016

#include "ThrottleTable.h"

namespace EMTG
{
    namespace HardwareModels
    {
        //comparator
        struct less_than
        {
            inline bool operator() (const ThrottleSetting& LHS, const ThrottleSetting& RHS)
            {
                double A = LHS.getThrusterPower();
                double B = RHS.getThrusterPower();
                return A < B ? true : false;
            }
        };

        //constructor for use with incoming data
        ThrottleTable::ThrottleTable(const std::string& inputfilename, const double& Sharpness) :
            control_type_2d(false)
        {
            this->PPUefficiency = 1.0;
            this->PPUminpower = 0.5;
            this->PPUmaxpower = 10.0;

            this->Sharpness = Sharpness;

            this->ParseThrottleTableFile(inputfilename);

            this->CreateThrottleSets();

            this->FindTransitionValues();

            this->printThrottleSets(inputfilename + "OUTPUT");
            
        }

        //function to parse a throttle table file
        void ThrottleTable::ParseThrottleTableFile(const std::string& inputfilename)
        {

            std::ifstream inputfile(inputfilename.c_str());
            std::string peek;
            std::string line;
            std::vector<std::string> linecell;
            char dump_buffer[1024];

            if (!inputfile.is_open())
            {
                throw std::invalid_argument("Cannot find throttle table file: " + inputfilename + ". Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }

            while (!inputfile.eof())
            {
                peek = inputfile.peek();
                if (peek == "#" || peek == "\r" || peek == "\n")
                {
                    //comment or blank line, do not parse
                    inputfile.getline(dump_buffer, 1024);
                }
                else
                {
                    getline(inputfile, line, '\n');
                    boost::split(linecell, line, boost::is_any_of(","));

					if (linecell[0] == "Throttle level")
					{
						continue;
					}
                    else if (linecell[0] == "PPU efficiency")
                    {
                        this->PPUefficiency = std::stod(linecell[1]);
                    }
                    else if (linecell[0] == "PPU min power (kW)")
                    {
                        this->PPUminpower = std::stod(linecell[1]);
                    }
                    else if (linecell[0] == "PPU max power (kW)")
                    {
                        this->PPUmaxpower = std::stod(linecell[1]);
                    }
                    else if (linecell[0] == "high_thrust_Thrust")
                    {
                        for (size_t k = 0; k < 5; ++k)
                            this->high_thrust_Thrust_coeff[k] = std::stod(linecell[1 + k]);
                    }
                    else if (linecell[0] == "high_thrust_Mdot")
                    {
                        for (size_t k = 0; k < 5; ++k)
                            this->high_thrust_Mdot_coeff[k] = std::stod(linecell[1 + k]);
                    }
                    else if (linecell[0] == "high_Isp_Thrust")
                    {
                        for (size_t k = 0; k < 5; ++k)
                            this->high_Isp_Thrust_coeff[k] = std::stod(linecell[1 + k]);
                    }
                    else if (linecell[0] == "high_Isp_Mdot")
                    {
                        for (size_t k = 0; k < 5; ++k)
                            this->high_Isp_Mdot_coeff[k] = std::stod(linecell[1 + k]);
                    }
                    else if (linecell[0] == "2Dpolyrow1")
                    {
                        for (size_t k = 0; k < 5; ++k)
                            this->poly2D_coefficients[0][k] = std::stod(linecell[1 + k]);
                    }
                    else if (linecell[0] == "2Dpolyrow2")
                    {
                        for (size_t k = 0; k < 5; ++k)
                            this->poly2D_coefficients[1][k] = std::stod(linecell[1 + k]);
                    }
                    else if (linecell[0] == "2Dpolyrow3")
                    {
                        for (size_t k = 0; k < 5; ++k)
                            this->poly2D_coefficients[2][k] = std::stod(linecell[1 + k]);
                    }
                    else if (linecell[0] == "2Dpolyrow4")
                    {
                        for (size_t k = 0; k < 5; ++k)
                            this->poly2D_coefficients[3][k] = std::stod(linecell[1 + k]);
                    }
                    else if (linecell[0] == "2Dpolyrow5")
                    {
                        for (size_t k = 0; k < 5; ++k)
                            this->poly2D_coefficients[4][k] = std::stod(linecell[1 + k]);
                    }
                    else
                    {
                        this->ThrottleSettings.push_back(ThrottleSetting(linecell[0],
                            std::stod(linecell[1]),
                            std::stod(linecell[2]),
                            std::stod(linecell[3]),
                            std::stod(linecell[4]),
                            std::stod(linecell[5]),
                            std::stod(linecell[6]),
                            std::stod(linecell[7]),
                            this->Sharpness));
                    }
                }
            }
        }//end throttle table parse function

        //function to create max-thrust, max-Isp, low-Mdot, and high-Mdot sets
        void ThrottleTable::CreateThrottleSets()
        {
            //create the high-Isp set first
            for (size_t settingIndex = 0; settingIndex < this->ThrottleSettings.size(); ++settingIndex)
            {
                bool isDominated = false;
                for (size_t comparisonSettingIndex = 0; comparisonSettingIndex < this->ThrottleSettings.size(); ++comparisonSettingIndex)
                {
                    if (this->ThrottleSettings[settingIndex].getIsp() < this->ThrottleSettings[comparisonSettingIndex].getIsp()
                        && this->ThrottleSettings[settingIndex].getThrusterPower() >= this->ThrottleSettings[comparisonSettingIndex].getThrusterPower())
                    {
                        isDominated = true;
                        break;
                    }
                }

                if (isDominated == false)
                    this->HighIspSet.push_back(this->ThrottleSettings[settingIndex]);
            }

            //sort the high-Isp set in order of power from least to greatest
            std::sort(this->HighIspSet.begin(), this->HighIspSet.end(), less_than());

            //set the Heaviside step increments for each entry in the high-Isp set
            //first set the first entry
            this->HighIspSet[0].setHeavisideThrustHeight(this->HighIspSet[0].getThrust());
            this->HighIspSet[0].setHeavisideMdotHeight(this->HighIspSet[0].getMdot());
            this->HighIspSet[0].setHeavisideIspHeight(this->HighIspSet[0].getIsp());
            this->HighIspSet[0].setHeavisideThrottleLevelHeight(this->HighIspSet[0].getThrottleLevel());
            this->HighIspSet[0].setHeavisideActivePowerHeight(this->HighIspSet[0].getThrusterPower());
            this->HighIspSet[0].setHeavisideBeamCurrentHeight(this->HighIspSet[0].getBeamCurrent());
            this->HighIspSet[0].setHeavisideBeamVoltageHeight(this->HighIspSet[0].getBeamVoltage());
            //set the remaining entries
            for (size_t settingIndex = 1; settingIndex < this->HighIspSet.size(); ++settingIndex)
            {
                this->HighIspSet[settingIndex].setHeavisideThrustHeight(this->HighIspSet[settingIndex].getThrust() - this->HighIspSet[settingIndex - 1].getThrust());
                this->HighIspSet[settingIndex].setHeavisideMdotHeight(this->HighIspSet[settingIndex].getMdot() - this->HighIspSet[settingIndex - 1].getMdot());
                this->HighIspSet[settingIndex].setHeavisideIspHeight(this->HighIspSet[settingIndex].getIsp() - this->HighIspSet[settingIndex - 1].getIsp());
                this->HighIspSet[settingIndex].setHeavisideThrottleLevelHeight(this->HighIspSet[settingIndex].getThrottleLevel() - this->HighIspSet[settingIndex - 1].getThrottleLevel());
                this->HighIspSet[settingIndex].setHeavisideActivePowerHeight(this->HighIspSet[settingIndex].getThrusterPower() - this->HighIspSet[settingIndex - 1].getThrusterPower());
                this->HighIspSet[settingIndex].setHeavisideBeamCurrentHeight(this->HighIspSet[settingIndex].getBeamCurrent() - this->HighIspSet[settingIndex - 1].getBeamCurrent());
                this->HighIspSet[settingIndex].setHeavisideBeamVoltageHeight(this->HighIspSet[settingIndex].getBeamVoltage() - this->HighIspSet[settingIndex - 1].getBeamVoltage());
            }
            
            //create the high-Thrust set
            for (size_t settingIndex = 0; settingIndex < this->ThrottleSettings.size(); ++settingIndex)
            {
                bool isDominated = false;
                for (size_t comparisonSettingIndex = 0; comparisonSettingIndex < this->ThrottleSettings.size(); ++comparisonSettingIndex)
                {
                    if (this->ThrottleSettings[settingIndex].getThrust() < this->ThrottleSettings[comparisonSettingIndex].getThrust()
                        && this->ThrottleSettings[settingIndex].getThrusterPower() >= this->ThrottleSettings[comparisonSettingIndex].getThrusterPower())
                    {
                        isDominated = true;
                        break;
                    }
                }

                if (isDominated == false)
                {
                    this->HighThrustSet.push_back(this->ThrottleSettings[settingIndex]);
                }
            }

            //sort the high-Thrust set in order of power from least to greatest
            std::sort(this->HighThrustSet.begin(), this->HighThrustSet.end(), less_than());

            //set the Heaviside step increments for each entry in the high-Thrust set
            //first set the first entry
            this->HighThrustSet[0].setHeavisideThrustHeight(this->HighThrustSet[0].getThrust());
            this->HighThrustSet[0].setHeavisideMdotHeight(this->HighThrustSet[0].getMdot());
            this->HighThrustSet[0].setHeavisideIspHeight(this->HighThrustSet[0].getIsp());
            this->HighThrustSet[0].setHeavisideThrottleLevelHeight(this->HighThrustSet[0].getThrottleLevel());
            this->HighThrustSet[0].setHeavisideActivePowerHeight(this->HighThrustSet[0].getThrusterPower());
            this->HighThrustSet[0].setHeavisideBeamCurrentHeight(this->HighThrustSet[0].getBeamCurrent());
            this->HighThrustSet[0].setHeavisideBeamVoltageHeight(this->HighThrustSet[0].getBeamVoltage());
            //set the remaining entries
            for (size_t settingIndex = 1; settingIndex < this->HighThrustSet.size(); ++settingIndex)
            {
                this->HighThrustSet[settingIndex].setHeavisideThrustHeight(this->HighThrustSet[settingIndex].getThrust() - this->HighThrustSet[settingIndex - 1].getThrust());
                this->HighThrustSet[settingIndex].setHeavisideMdotHeight(this->HighThrustSet[settingIndex].getMdot() - this->HighThrustSet[settingIndex - 1].getMdot());
                this->HighThrustSet[settingIndex].setHeavisideIspHeight(this->HighThrustSet[settingIndex].getIsp() - this->HighThrustSet[settingIndex - 1].getIsp());
                this->HighThrustSet[settingIndex].setHeavisideThrottleLevelHeight(this->HighThrustSet[settingIndex].getThrottleLevel() - this->HighThrustSet[settingIndex - 1].getThrottleLevel());
                this->HighThrustSet[settingIndex].setHeavisideActivePowerHeight(this->HighThrustSet[settingIndex].getThrusterPower() - this->HighThrustSet[settingIndex - 1].getThrusterPower());
                this->HighThrustSet[settingIndex].setHeavisideBeamCurrentHeight(this->HighThrustSet[settingIndex].getBeamCurrent() - this->HighThrustSet[settingIndex - 1].getBeamCurrent());
                this->HighThrustSet[settingIndex].setHeavisideBeamVoltageHeight(this->HighThrustSet[settingIndex].getBeamVoltage() - this->HighThrustSet[settingIndex - 1].getBeamVoltage());
            }

            //create the high-Mdot set
            for (size_t settingIndex = 0; settingIndex < this->ThrottleSettings.size(); ++settingIndex)
            {
                bool isDominated = false;
                for (size_t comparisonSettingIndex = 0; comparisonSettingIndex < this->ThrottleSettings.size(); ++comparisonSettingIndex)
                {
                    if (this->ThrottleSettings[settingIndex].getMdot() < this->ThrottleSettings[comparisonSettingIndex].getMdot()
                        && this->ThrottleSettings[settingIndex].getThrusterPower() >= this->ThrottleSettings[comparisonSettingIndex].getThrusterPower())
                    {
                        isDominated = true;
                        break;
                    }
                }

                if (isDominated == false)
                {
                    this->HighMdotSet.push_back(this->ThrottleSettings[settingIndex]);
                }
            }

            //sort the high-Thrust set in order of power from least to greatest
            std::sort(this->HighMdotSet.begin(), this->HighMdotSet.end(), less_than());

            //set the Heaviside step increments for each entry in the high-Thrust set
            //first set the first entry
            this->HighMdotSet[0].setHeavisideThrustHeight(this->HighMdotSet[0].getThrust());
            this->HighMdotSet[0].setHeavisideMdotHeight(this->HighMdotSet[0].getMdot());
            this->HighMdotSet[0].setHeavisideIspHeight(this->HighMdotSet[0].getIsp());
            this->HighMdotSet[0].setHeavisideThrottleLevelHeight(this->HighMdotSet[0].getThrottleLevel());
            this->HighMdotSet[0].setHeavisideActivePowerHeight(this->HighMdotSet[0].getThrusterPower());
            this->HighMdotSet[0].setHeavisideBeamCurrentHeight(this->HighMdotSet[0].getBeamCurrent());
            this->HighMdotSet[0].setHeavisideBeamVoltageHeight(this->HighMdotSet[0].getBeamVoltage());
            //set the remaining entries
            for (size_t settingIndex = 1; settingIndex < this->HighMdotSet.size(); ++settingIndex)
            {
                this->HighMdotSet[settingIndex].setHeavisideThrustHeight(this->HighMdotSet[settingIndex].getThrust() - this->HighMdotSet[settingIndex - 1].getThrust());
                this->HighMdotSet[settingIndex].setHeavisideMdotHeight(this->HighMdotSet[settingIndex].getMdot() - this->HighMdotSet[settingIndex - 1].getMdot());
                this->HighMdotSet[settingIndex].setHeavisideIspHeight(this->HighMdotSet[settingIndex].getIsp() - this->HighMdotSet[settingIndex - 1].getIsp());
                this->HighMdotSet[settingIndex].setHeavisideThrottleLevelHeight(this->HighMdotSet[settingIndex].getThrottleLevel() - this->HighMdotSet[settingIndex - 1].getThrottleLevel());
                this->HighMdotSet[settingIndex].setHeavisideActivePowerHeight(this->HighMdotSet[settingIndex].getThrusterPower() - this->HighMdotSet[settingIndex - 1].getThrusterPower());
                this->HighMdotSet[settingIndex].setHeavisideBeamCurrentHeight(this->HighMdotSet[settingIndex].getBeamCurrent() - this->HighMdotSet[settingIndex - 1].getBeamCurrent());
                this->HighMdotSet[settingIndex].setHeavisideBeamVoltageHeight(this->HighMdotSet[settingIndex].getBeamVoltage() - this->HighMdotSet[settingIndex - 1].getBeamVoltage());
            }

            //create the low-Mdot set
            for (size_t settingIndex = 0; settingIndex < this->ThrottleSettings.size(); ++settingIndex)
            {
                bool isDominated = false;
                for (size_t comparisonSettingIndex = 0; comparisonSettingIndex < this->ThrottleSettings.size(); ++comparisonSettingIndex)
                {
                    if (this->ThrottleSettings[settingIndex].getMdot() > this->ThrottleSettings[comparisonSettingIndex].getMdot()
                        && this->ThrottleSettings[settingIndex].getThrusterPower() <= this->ThrottleSettings[comparisonSettingIndex].getThrusterPower())
                    {
                        isDominated = true;
                        break;
                    }
                }

                if (isDominated == false)
                {
                    this->LowMdotSet.push_back(this->ThrottleSettings[settingIndex]);
                }
            }

            //sort the high-Thrust set in order of power from least to greatest
            std::sort(this->LowMdotSet.begin(), this->LowMdotSet.end(), less_than());

            //set the Heaviside step increments for each entry in the high-Thrust set
            //first set the first entry
            this->LowMdotSet[0].setHeavisideThrustHeight(this->LowMdotSet[0].getThrust());
            this->LowMdotSet[0].setHeavisideMdotHeight(this->LowMdotSet[0].getMdot());
            this->LowMdotSet[0].setHeavisideIspHeight(this->LowMdotSet[0].getIsp());
            this->LowMdotSet[0].setHeavisideThrottleLevelHeight(this->LowMdotSet[0].getThrottleLevel());
            this->LowMdotSet[0].setHeavisideActivePowerHeight(this->LowMdotSet[0].getThrusterPower());
            this->LowMdotSet[0].setHeavisideBeamCurrentHeight(this->LowMdotSet[0].getBeamCurrent());
            this->LowMdotSet[0].setHeavisideBeamVoltageHeight(this->LowMdotSet[0].getBeamVoltage());
            //set the remaining entries
            for (size_t settingIndex = 1; settingIndex < this->LowMdotSet.size(); ++settingIndex)
            {
                this->LowMdotSet[settingIndex].setHeavisideThrustHeight(this->LowMdotSet[settingIndex].getThrust() - this->LowMdotSet[settingIndex - 1].getThrust());
                this->LowMdotSet[settingIndex].setHeavisideMdotHeight(this->LowMdotSet[settingIndex].getMdot() - this->LowMdotSet[settingIndex - 1].getMdot());
                this->LowMdotSet[settingIndex].setHeavisideIspHeight(this->LowMdotSet[settingIndex].getIsp() - this->LowMdotSet[settingIndex - 1].getIsp());
                this->LowMdotSet[settingIndex].setHeavisideThrottleLevelHeight(this->LowMdotSet[settingIndex].getThrottleLevel() - this->LowMdotSet[settingIndex - 1].getThrottleLevel());
                this->LowMdotSet[settingIndex].setHeavisideActivePowerHeight(this->LowMdotSet[settingIndex].getThrusterPower() - this->LowMdotSet[settingIndex - 1].getThrusterPower());
                this->LowMdotSet[settingIndex].setHeavisideBeamCurrentHeight(this->LowMdotSet[settingIndex].getBeamCurrent() - this->LowMdotSet[settingIndex - 1].getBeamCurrent());
                this->LowMdotSet[settingIndex].setHeavisideBeamVoltageHeight(this->LowMdotSet[settingIndex].getBeamVoltage() - this->LowMdotSet[settingIndex - 1].getBeamVoltage());
            }

            //create the high-Efficiency set
            for (size_t settingIndex = 0; settingIndex < this->ThrottleSettings.size(); ++settingIndex)
            {
                bool isDominated = false;
                for (size_t comparisonSettingIndex = 0; comparisonSettingIndex < this->ThrottleSettings.size(); ++comparisonSettingIndex)
                {
                    if (this->ThrottleSettings[settingIndex].getEfficiency() < this->ThrottleSettings[comparisonSettingIndex].getEfficiency()
                        && this->ThrottleSettings[settingIndex].getThrusterPower() >= this->ThrottleSettings[comparisonSettingIndex].getThrusterPower())
                    {
                        isDominated = true;
                        break;
                    }
                }

                if (isDominated == false)
                {
                    this->HighEfficiencySet.push_back(this->ThrottleSettings[settingIndex]);
                }
            }

            //sort the high-Thrust set in order of power from least to greatest
            std::sort(this->HighEfficiencySet.begin(), this->HighEfficiencySet.end(), less_than());

            //set the Heaviside step increments for each entry in the high-Thrust set
            //first set the first entry
            this->HighEfficiencySet[0].setHeavisideThrustHeight(this->HighEfficiencySet[0].getThrust());
            this->HighEfficiencySet[0].setHeavisideMdotHeight(this->HighEfficiencySet[0].getMdot());
            this->HighEfficiencySet[0].setHeavisideIspHeight(this->HighEfficiencySet[0].getIsp());
            this->HighEfficiencySet[0].setHeavisideThrottleLevelHeight(this->HighEfficiencySet[0].getThrottleLevel());
            this->HighEfficiencySet[0].setHeavisideActivePowerHeight(this->HighEfficiencySet[0].getThrusterPower());
            this->HighEfficiencySet[0].setHeavisideBeamCurrentHeight(this->HighEfficiencySet[0].getBeamCurrent());
            this->HighEfficiencySet[0].setHeavisideBeamVoltageHeight(this->HighEfficiencySet[0].getBeamVoltage());
            //set the remaining entries
            for (size_t settingIndex = 1; settingIndex < this->HighEfficiencySet.size(); ++settingIndex)
            {
                this->HighEfficiencySet[settingIndex].setHeavisideThrustHeight(this->HighEfficiencySet[settingIndex].getThrust() - this->HighEfficiencySet[settingIndex - 1].getThrust());
                this->HighEfficiencySet[settingIndex].setHeavisideMdotHeight(this->HighEfficiencySet[settingIndex].getMdot() - this->HighEfficiencySet[settingIndex - 1].getMdot());
                this->HighEfficiencySet[settingIndex].setHeavisideIspHeight(this->HighEfficiencySet[settingIndex].getIsp() - this->HighEfficiencySet[settingIndex - 1].getIsp());
                this->HighEfficiencySet[settingIndex].setHeavisideThrottleLevelHeight(this->HighEfficiencySet[settingIndex].getThrottleLevel() - this->HighEfficiencySet[settingIndex - 1].getThrottleLevel());
                this->HighEfficiencySet[settingIndex].setHeavisideActivePowerHeight(this->HighEfficiencySet[settingIndex].getThrusterPower() - this->HighEfficiencySet[settingIndex - 1].getThrusterPower());
                this->HighEfficiencySet[settingIndex].setHeavisideBeamCurrentHeight(this->HighEfficiencySet[settingIndex].getBeamCurrent() - this->HighEfficiencySet[settingIndex - 1].getBeamCurrent());
                this->HighEfficiencySet[settingIndex].setHeavisideBeamVoltageHeight(this->HighEfficiencySet[settingIndex].getBeamVoltage() - this->HighEfficiencySet[settingIndex - 1].getBeamVoltage());
            }

            //create the high-Efficiency set
            this->FullSet = this->ThrottleSettings;

            //sort the high-Thrust set in order of power from least to greatest
            std::sort(this->FullSet.begin(), this->FullSet.end(), less_than());

            //set the Heaviside step increments for each entry in the high-Thrust set
            //first set the first entry
            this->FullSet[0].setHeavisideThrustHeight(this->FullSet[0].getThrust());
            this->FullSet[0].setHeavisideMdotHeight(this->FullSet[0].getMdot());
            this->FullSet[0].setHeavisideIspHeight(this->FullSet[0].getIsp());
            this->FullSet[0].setHeavisideThrottleLevelHeight(this->FullSet[0].getThrottleLevel());
            this->FullSet[0].setHeavisideActivePowerHeight(this->FullSet[0].getThrusterPower());
            this->FullSet[0].setHeavisideBeamCurrentHeight(this->FullSet[0].getBeamCurrent());
            this->FullSet[0].setHeavisideBeamVoltageHeight(this->FullSet[0].getBeamVoltage());
            //set the remaining entries
            for (size_t settingIndex = 1; settingIndex < this->FullSet.size(); ++settingIndex)
            {
                this->FullSet[settingIndex].setHeavisideThrustHeight(this->FullSet[settingIndex].getThrust() - this->FullSet[settingIndex - 1].getThrust());
                this->FullSet[settingIndex].setHeavisideMdotHeight(this->FullSet[settingIndex].getMdot() - this->FullSet[settingIndex - 1].getMdot());
                this->FullSet[settingIndex].setHeavisideIspHeight(this->FullSet[settingIndex].getIsp() - this->FullSet[settingIndex - 1].getIsp());
                this->FullSet[settingIndex].setHeavisideThrottleLevelHeight(this->FullSet[settingIndex].getThrottleLevel() - this->FullSet[settingIndex - 1].getThrottleLevel());
                this->FullSet[settingIndex].setHeavisideActivePowerHeight(this->FullSet[settingIndex].getThrusterPower() - this->FullSet[settingIndex - 1].getThrusterPower());
                this->FullSet[settingIndex].setHeavisideBeamCurrentHeight(this->FullSet[settingIndex].getBeamCurrent() - this->FullSet[settingIndex - 1].getBeamCurrent());
                this->FullSet[settingIndex].setHeavisideBeamVoltageHeight(this->FullSet[settingIndex].getBeamVoltage() - this->FullSet[settingIndex - 1].getBeamVoltage());
            }
        }//end function to create throttle sets



         //function to find the transition voltage and mdot values for each throttle setting
        void ThrottleTable::FindTransitionValues()
        {             
            //first thing we need to do is build up a list of voltage and Mdot values in the throttle set
            this->Mdot_choices.push_back(this->ThrottleSettings.front().getMdot());
            this->Voltage_choices.push_back(this->ThrottleSettings.front().getBeamVoltage());

            this->Power_max = this->ThrottleSettings[0].getThrusterPower();
            
            for (size_t settingIndex = 1; settingIndex < this->ThrottleSettings.size(); ++settingIndex)
            {
                this->Power_max = this->ThrottleSettings[settingIndex].getThrusterPower() > this->Power_max ? this->ThrottleSettings[settingIndex].getThrusterPower() : this->Power_max;
                
                //is this setting's Mdot in the Mdot_choices vector? If not, add it.
                bool Mdot_already_present = false;
                for (size_t Mdot_index = 0; Mdot_index < this->Mdot_choices.size(); ++Mdot_index)
                {
                    if (fabs(this->ThrottleSettings[settingIndex].getMdot() - this->Mdot_choices[Mdot_index]) < math::SMALL)
                    {
                        Mdot_already_present = true;
                        break;
                    }
                }
                if (!Mdot_already_present)
                    this->Mdot_choices.push_back(this->ThrottleSettings[settingIndex].getMdot());

                //is this setting's BeamVoltage in the BeamVoltage_choices vector? If not, add it.
                bool BeamVoltage_already_present = false;
                for (size_t BeamVoltage_index = 0; BeamVoltage_index < this->Voltage_choices.size(); ++BeamVoltage_index)
                {
                    if (fabs(this->ThrottleSettings[settingIndex].getBeamVoltage() - this->Voltage_choices[BeamVoltage_index]) < math::SMALL)
                    {
                        BeamVoltage_already_present = true;
                        break;
                    }
                }
                if (!BeamVoltage_already_present)
                    this->Voltage_choices.push_back(this->ThrottleSettings[settingIndex].getBeamVoltage());
            }
            std::sort(this->Mdot_choices.begin(), this->Mdot_choices.end());
            std::sort(this->Voltage_choices.begin(), this->Voltage_choices.end());
            
            this->Mdot_max = this->Mdot_choices.back();
            this->Voltage_max = this->Voltage_choices.back();
            
            this->Mdot_choices.push_back(2.0 * this->Mdot_choices.back());
            this->Voltage_choices.push_back(2.0 * this->Voltage_choices.back());

            std::vector< ThrottleSetting > FakeSettings;

            //now set the on and off settings for each throttle setting
            for (size_t settingIndex = 0; settingIndex < this->ThrottleSettings.size(); ++settingIndex)
            {
                double control_setting;
                if (this->control_type_2d)
                {                    
                    control_setting = this->ThrottleSettings[settingIndex].getMdot();
                    this->ThrottleSettings[settingIndex].setMdotOn(control_setting);
                }
                else 
                {
                    control_setting = this->ThrottleSettings[settingIndex].getBeamVoltage();
                    this->ThrottleSettings[settingIndex].setBeamVoltageOn(control_setting);
                }
                
                double currentP = this->ThrottleSettings[settingIndex].getThrusterPower();
                this->ThrottleSettings[settingIndex].setPowerOn(currentP);
      
                double power_delta = 1e100;
                double control_delta = 1e100;
                double minP = 1e100;
                double nextControl = 1e100;
                double nextP = 1.0e+100;
                
                for (size_t settingIndex2 = 0; settingIndex2 < this->ThrottleSettings.size(); ++settingIndex2)
                {
                    double control_setting2;
                    if (control_type_2d)
                    {
                        control_setting2 = this->ThrottleSettings[settingIndex2].getMdot();
                    }
                    else
                    {
                        control_setting2 = this->ThrottleSettings[settingIndex2].getBeamVoltage();
                    }
                    double power_setting2   = this->ThrottleSettings[settingIndex2].getThrusterPower();                    
                    
                    if (control_setting2 > control_setting && control_setting2-control_setting < control_delta)
                    {
                        control_delta = control_setting2-control_setting;
                        nextControl = control_setting2;                        
                    }                    
                    
                    if (control_setting2 == control_setting && power_setting2 > currentP && power_setting2-currentP < power_delta)
                    {
                        power_delta = power_setting2 - currentP;
                        nextP = power_setting2;                      
                    }
                    
                    if (control_setting2 == nextControl && power_setting2 < minP)
                    {
                        minP = power_setting2;
                    }
                }                
                
                if (currentP < minP)
                {
                    if (control_type_2d)
                    {
                        this->ThrottleSettings[settingIndex].setMdotOff(this->Mdot_max * 2);
                    }
                    else
                    {
                        this->ThrottleSettings[settingIndex].setBeamVoltageOff(this->Voltage_max * 2);
                    }
                    
                    if (nextP > minP)
                    {
                        FakeSettings.push_back(ThrottleSetting(
                            this->ThrottleSettings[settingIndex].getThrottleLevelString(),
                            this->ThrottleSettings[settingIndex].getMdot(),
                            this->ThrottleSettings[settingIndex].getBeamCurrent(),
                            this->ThrottleSettings[settingIndex].getBeamVoltage(),
                            this->ThrottleSettings[settingIndex].getThrust(),
                            this->ThrottleSettings[settingIndex].getIsp(),
                            this->ThrottleSettings[settingIndex].getEfficiency(),
                            minP,
                            this->ThrottleSettings[settingIndex].getSharpness()                            
                            ));

                        if (control_type_2d)
                        {
                            FakeSettings.back().setMdotOff(nextControl);
                            FakeSettings.back().setMdotOn(control_setting);
                        }
                        else
                        {
                            FakeSettings.back().setBeamVoltageOff(nextControl); 
                            FakeSettings.back().setBeamVoltageOn(control_setting); 
                        }
                        
                        FakeSettings.back().setPowerOff(nextP);
                        FakeSettings.back().setPowerOn(minP);
                        this->ThrottleSettings[settingIndex].setPowerOff(minP);
                    }  
                    else
                    {
                        this->ThrottleSettings[settingIndex].setPowerOff(nextP); 
                    }                        
                }
                else
                {                
                    this->ThrottleSettings[settingIndex].setPowerOff(nextP); 
                    
                    if (control_type_2d)
                    {
                        this->ThrottleSettings[settingIndex].setMdotOff(nextControl);
                    }
                    else
                    {
                        this->ThrottleSettings[settingIndex].setBeamVoltageOff(nextControl);
                    }
                }
            }
            
            for (size_t settingIndex = 0; settingIndex < FakeSettings.size(); ++settingIndex)
            {
                this->ThrottleSettings.push_back(FakeSettings[settingIndex]);
            }

        }//end function to find on and off switches for all throttle settings

        //function to calculate thruster performance at a given power level
        void ThrottleTable::CalculateThrusterPerformance1D(const doubleType& inputPower, const ThrottleSetDef& PreferredThrottleSet)
        {

            this->HThrust = 0.0;
            this->Hmdot = 0.0;
            this->HIsp = 0.0;
            this->dHthrustdP = 0.0;
            this->dHmdotdP = 0.0;
            this->dHIspdP = 0.0;
            this->dHBeamVoltagedP = 0.0;
            this->HactivePower = 0.0;
            this->HThrottleLevel = 0.0;
            this->HBeamCurrent = 0.0;
            this->HBeamVoltage = 0.0;

            this->ThrusterInputPower = inputPower * this->PPUefficiency;
            std::vector< ThrottleSetting >* ThrottleTableofInterest = nullptr;
            switch (PreferredThrottleSet)
            {
                case ThrottleSetDef::HighIsp:
                    ThrottleTableofInterest = &this->HighIspSet;
                    break;
                case ThrottleSetDef::HighThrust:
                    ThrottleTableofInterest = &this->HighThrustSet;
                    break;
                case ThrottleSetDef::LowMdot:
                    ThrottleTableofInterest = &this->LowMdotSet;
                    break;
                case ThrottleSetDef::HighMdot:
                    ThrottleTableofInterest = &this->HighMdotSet;
                    break;
                case ThrottleSetDef::HighEfficiency:
                    ThrottleTableofInterest = &this->HighEfficiencySet;
                    break;
                case ThrottleSetDef::FullSet:
                    ThrottleTableofInterest = &this->FullSet;
                    break;
            }
            for (ThrottleSetting& thisSetting : (*ThrottleTableofInterest))//size_t settingIndex = 0; settingIndex < ThrottleTableofInterest->size(); ++settingIndex)
            {
                thisSetting.calculateHeavisideValues1D(this->ThrusterInputPower);
                this->HThrust         += thisSetting.getHeavisideThrust();
                this->Hmdot           += thisSetting.getHeavisideMdot();
                this->HIsp            += thisSetting.getHeavisideIsp();
                this->HactivePower    += thisSetting.getHeavisideActivePower();
                this->HThrottleLevel  += thisSetting.getHeavisideThrottleLevel();
                this->HBeamCurrent    += thisSetting.getHeavisideBeamCurrent();
                this->HBeamVoltage    += thisSetting.getHeavisideBeamVoltage();
                this->dHthrustdP      += thisSetting.getHeavisideThrustDerivative() * this->PPUefficiency;
                this->dHmdotdP        += thisSetting.getHeavisideMdotDerivative() * this->PPUefficiency;
                this->dHIspdP         += thisSetting.getHeavisideIspDerivative() * this->PPUefficiency;
                this->dHBeamVoltagedP += thisSetting.getHeavisideBeamVoltage()_GETVALUE * this->PPUefficiency;

                if (thisSetting.getH() > 0.5)
                    this->ThrottleLevelString = thisSetting.getThrottleLevelString();
            }
            // cout<<this->HThrottleLevel<<","<<this->HThrust<<","<<this->Hmdot<<","<<this->HactivePower<<","<<this->HIsp<<endl;
        }//end function to calculate thruster performance at a given power level


        void ThrottleTable::CalculateThrusterPerformance2D(const doubleType& inputPower, const doubleType& u_command)
        {
            doubleType command;

            this->ThrusterInputPower = inputPower * this->PPUefficiency;    
            
            if (this->control_type_2d)
            {
                command = u_command * this->Mdot_max * 1.1;
            }
            else
            {
                command = u_command * this->Voltage_max * 1.1;
            }
            
            this->HThrust = 0.0;
            this->Hmdot = 0.0;
            this->HIsp = 0.0;
            this->dHthrustdP = 0.0;
            this->dHmdotdP = 0.0;
            this->dHIspdP = 0.0;
            this->dHBeamVoltagedP = 0.0;
            this->HactivePower = 0.0;
            this->HThrottleLevel = 0.0;
            this->HBeamCurrent = 0.0;
            this->HBeamVoltage = 0.0;
            this->dHthrust_dcommand = 0.0;
            this->dHmdot_dcommand = 0.0;
            this->dHIsp_dcommand = 0.0;
            this->dHthrust_du_command = 0.0;
            this->dHmdot_du_command = 0.0;
            this->dHIsp_du_command = 0.0;
            
            for (size_t settingIndex = 0; settingIndex < this->ThrottleSettings.size(); ++settingIndex)
            {
                if (this->control_type_2d)
                {
                    this->ThrottleSettings[settingIndex].calculateHeavisideValues2D_mdot(command, this->ThrusterInputPower, this->Mdot_max, this->Power_max);
                    this->dHthrust_du_command += this->ThrottleSettings[settingIndex].getHeavisideThrustCommandedMdotDerivative() ;//* dcommanded_voltage_du_Voltage;
                    this->dHmdot_du_command += this->ThrottleSettings[settingIndex].getHeavisideMdotCommandedMdotDerivative() ;//* dcommanded_voltage_du_Voltage;
                    this->dHIsp_du_command += this->ThrottleSettings[settingIndex].getHeavisideIspCommandedMdotDerivative() ;//* dcommanded_voltage_du_Voltage;
                }
                else
                {
                    this->ThrottleSettings[settingIndex].calculateHeavisideValues2D_voltage(command, this->ThrusterInputPower, this->Voltage_max, this->Power_max);
                    this->dHthrust_du_command += this->ThrottleSettings[settingIndex].getHeavisideThrustCommandedVoltageDerivative() ;//* dcommanded_voltage_du_Voltage;
                    this->dHmdot_du_command += this->ThrottleSettings[settingIndex].getHeavisideMdotCommandedVoltageDerivative() ;//* dcommanded_voltage_du_Voltage;
                    this->dHIsp_du_command += this->ThrottleSettings[settingIndex].getHeavisideIspCommandedVoltageDerivative() ;//* dcommanded_voltage_du_Voltage;
                }

                this->HThrust += this->ThrottleSettings[settingIndex].getHeavisideThrust();
                this->Hmdot += this->ThrottleSettings[settingIndex].getHeavisideMdot();
                this->HIsp += this->ThrottleSettings[settingIndex].getHeavisideIsp();
                this->HactivePower += this->ThrottleSettings[settingIndex].getHeavisideActivePower();
                this->HThrottleLevel += this->ThrottleSettings[settingIndex].getHeavisideThrottleLevel();
                this->HBeamCurrent += this->ThrottleSettings[settingIndex].getHeavisideBeamCurrent();
                this->HBeamVoltage += this->ThrottleSettings[settingIndex].getHeavisideBeamVoltage();
                this->dHthrustdP += this->ThrottleSettings[settingIndex].getHeavisideThrustInputPowerDerivative();
                this->dHmdotdP   += this->ThrottleSettings[settingIndex].getHeavisideMdotInputPowerDerivative();
                this->dHIspdP    += this->ThrottleSettings[settingIndex].getHeavisideIspInputPowerDerivative();
            }
        }

        //get functions
        double ThrottleTable::getPPUminpower() const
        {
            return this->PPUminpower;
        }

        double ThrottleTable::getPPUmaxpower() const
        {
            return this->PPUmaxpower;
        }
        double ThrottleTable::getThrottleSetminpower(const ThrottleSetDef& PreferredThrottleSet)
        {
            std::vector< ThrottleSetting >* ThrottleTableofInterest = nullptr;
            switch (PreferredThrottleSet)
            {
            case ThrottleSetDef::HighIsp:
                ThrottleTableofInterest = &this->HighIspSet;
                break;
            case ThrottleSetDef::HighThrust:
                ThrottleTableofInterest = &this->HighThrustSet;
                break;
            case ThrottleSetDef::LowMdot:
                ThrottleTableofInterest = &this->LowMdotSet;
                break;
            case ThrottleSetDef::HighMdot:
                ThrottleTableofInterest = &this->HighMdotSet;
                break;
            case ThrottleSetDef::HighEfficiency:
                ThrottleTableofInterest = &this->HighEfficiencySet;
                break;
            case ThrottleSetDef::FullSet:
                ThrottleTableofInterest = &this->FullSet;
                break;
            }

            return ThrottleTableofInterest->begin()->getThrusterPower();
        }

        double ThrottleTable::getThrottleSetmaxpower(const ThrottleSetDef& PreferredThrottleSet)
        {
            std::vector< ThrottleSetting >* ThrottleTableofInterest = nullptr;
            switch (PreferredThrottleSet)
            {
            case ThrottleSetDef::HighIsp:
                ThrottleTableofInterest = &this->HighIspSet;
                break;
            case ThrottleSetDef::HighThrust:
                ThrottleTableofInterest = &this->HighThrustSet;
                break;
            case ThrottleSetDef::LowMdot:
                ThrottleTableofInterest = &this->LowMdotSet;
                break;
            case ThrottleSetDef::HighMdot:
                ThrottleTableofInterest = &this->HighMdotSet;
                break;
            case ThrottleSetDef::HighEfficiency:
                ThrottleTableofInterest = &this->HighEfficiencySet;
                break;
            case ThrottleSetDef::FullSet:
                ThrottleTableofInterest = &this->FullSet;
                break;
            }

            return ThrottleTableofInterest->end()->getThrusterPower();
        }
        double ThrottleTable::getPPUefficiency() const
        {
            return this->PPUefficiency;
        }

        doubleType ThrottleTable::getHeavisideThrust() const
        {
            return this->HThrust;
        }

        doubleType ThrottleTable::getHeavisideMdot() const
        {
            return this->Hmdot;
        }

        doubleType ThrottleTable::getHeavisideIsp() const
        {
            return this->HIsp;
        }

        doubleType ThrottleTable::getHeavisideActivePower() const
        {
            return this->HactivePower / this->PPUefficiency;
        }

        double ThrottleTable::getHeavisideThrustDerivative() const
        {
            return this->dHthrustdP;
        }

        double ThrottleTable::getHeavisideMdotDerivative() const
        {
            return this->dHmdotdP;
        }

        double ThrottleTable::getHeavisideIspDerivative() const
        {
            return this->dHIspdP;
        }

        int ThrottleTable::getThrottleLevel() const
        {
            return (int)round(this->HThrottleLevel _GETVALUE);
        }

        double ThrottleTable::getdHthrust_dcommand() const
        {
            return this->dHthrust_du_command;
        }

        double ThrottleTable::getdHmdot_dcommand() const
        {
            return this->dHmdot_du_command;
        }

        double ThrottleTable::getdHIsp_dcommand() const
        {
            return this->dHIsp_du_command;
        }

        //various functions for polynomial throttle modeling
        void ThrottleTable::CalculateThrusterPerformance1D_Polynomial(const doubleType& inputPower, const ThrottleSetDef& PreferredThrottleSet)
        {
            doubleType P2 = inputPower * inputPower;
            doubleType P3 = P2 * inputPower;
            doubleType P4 = P3 * inputPower;


            if (PreferredThrottleSet == HighThrust)
            {
                this->polyThrust_1D = (this->high_thrust_Thrust_coeff[0]
                    + this->high_thrust_Thrust_coeff[1] * inputPower
                    + this->high_thrust_Thrust_coeff[2] * P2
                    + this->high_thrust_Thrust_coeff[3] * P3
                    + this->high_thrust_Thrust_coeff[4] * P4) * 1.0e+3;

                this->poly_dThrust_dP_1D = (this->high_thrust_Thrust_coeff[1]
                    + 2.0 * this->high_thrust_Thrust_coeff[2] * inputPower _GETVALUE
                    + 3.0 * this->high_thrust_Thrust_coeff[3] * P2 _GETVALUE
                    + 4.0 * this->high_thrust_Thrust_coeff[4] * P3 _GETVALUE) * 1.0e+3;

                this->polyMdot_1D = (this->high_thrust_Mdot_coeff[0]
                    + this->high_thrust_Mdot_coeff[1] * inputPower
                    + this->high_thrust_Mdot_coeff[2] * P2
                    + this->high_thrust_Mdot_coeff[3] * P3
                    + this->high_thrust_Mdot_coeff[4] * P4) * 1.0e+6;

                this->poly_dMdotdP_1D = (this->high_thrust_Mdot_coeff[1]
                    + 2.0 * this->high_thrust_Mdot_coeff[2] * inputPower _GETVALUE
                    + 3.0 * this->high_thrust_Mdot_coeff[3] * P2 _GETVALUE
                    + 4.0 * this->high_thrust_Mdot_coeff[4] * P3 _GETVALUE) * 1.0e+6;

                this->d2HighThrust_1D_Mdot_dP2 = (this->high_thrust_Mdot_coeff[2]
                    + 2.0 * this->high_thrust_Mdot_coeff[3] * inputPower _GETVALUE
                    + 6.0 * this->high_thrust_Mdot_coeff[4] * P2 _GETVALUE) * 1.0e+6;
            }
            else
            {
                this->polyThrust_1D = (this->high_Isp_Thrust_coeff[0]
                    + this->high_Isp_Thrust_coeff[1] * inputPower
                    + this->high_Isp_Thrust_coeff[2] * P2
                    + this->high_Isp_Thrust_coeff[3] * P3
                    + this->high_Isp_Thrust_coeff[4] * P4) * 1.0e+3;

                this->poly_dThrust_dP_1D = (this->high_Isp_Thrust_coeff[1]
                    + 2.0 * this->high_Isp_Thrust_coeff[2] * inputPower _GETVALUE
                    + 3.0 * this->high_Isp_Thrust_coeff[3] * P2 _GETVALUE
                    + 4.0 * this->high_Isp_Thrust_coeff[4] * P3 _GETVALUE) * 1.0e+3;

                this->polyMdot_1D = (this->high_Isp_Mdot_coeff[0]
                    + this->high_Isp_Mdot_coeff[1] * inputPower
                    + this->high_Isp_Mdot_coeff[2] * P2
                    + this->high_Isp_Mdot_coeff[3] * P3
                    + this->high_Isp_Mdot_coeff[4] * P4) * 1.0e+6;

                this->poly_dMdotdP_1D = (this->high_Isp_Mdot_coeff[1]
                    + 2.0 * this->high_Isp_Mdot_coeff[2] * inputPower _GETVALUE
                    + 3.0 * this->high_Isp_Mdot_coeff[3] * P2 _GETVALUE
                    + 4.0 * this->high_Isp_Mdot_coeff[4] * P3 _GETVALUE) * 1.0e+6;

                this->d2HighIsp_1D_Mdot_dP2 = (this->high_thrust_Mdot_coeff[2]
                    + 2.0 * this->high_thrust_Mdot_coeff[3] * inputPower _GETVALUE
                    + 6.0 * this->high_thrust_Mdot_coeff[4] * P2 _GETVALUE) * 1.0e+6;
            }
        }

        void ThrottleTable::CalculateThrusterPerformance2D_Polynomial(const doubleType& inputPower, const doubleType& u_Mdot)
        {
            //Step 1: evaluate high-thrust and high-Isp 1D polynomials for Mdot to get bounds on the available Mdot
            //this->CalculateThrusterPerformance1D_Polynomial(inputPower, HighThrust);
            //doubleType HighThrust_1D_Mdot = this->get1D_poly_Mdot();
            //double dHighThrust_1D_Mdot_dP = this->get1D_poly_dMdot_dP();
            this->CalculateThrusterPerformance1D(inputPower, HighMdot);
            doubleType HighThrust_1D_Mdot = this->getHeavisideMdot();
            double dHighThrust_1D_Mdot_dP = this->getHeavisideMdotDerivative();

            //this->CalculateThrusterPerformance1D_Polynomial(inputPower, HighIsp);
            //doubleType HighIsp_1D_Mdot = this->get1D_poly_Mdot();
            //double dHighIsp_1D_Mdot_dP = this->get1D_poly_dMdot_dP();
            this->CalculateThrusterPerformance1D(inputPower, LowMdot);
            doubleType HighIsp_1D_Mdot = this->getHeavisideMdot();
            double dHighIsp_1D_Mdot_dP = this->getHeavisideMdotDerivative();

            //Step 2: compute commanded Mdot
            this->polyMdot_2D = u_Mdot * (HighThrust_1D_Mdot - HighIsp_1D_Mdot) + HighIsp_1D_Mdot;
            this->poly_dMdot_du_Mdot_2D = (HighThrust_1D_Mdot - HighIsp_1D_Mdot)_GETVALUE;
            this->poly_dMdot_dP_2D = u_Mdot _GETVALUE * (dHighThrust_1D_Mdot_dP - dHighIsp_1D_Mdot_dP) + dHighIsp_1D_Mdot_dP;
            double d2poly_dMdot_dP2_2D = u_Mdot _GETVALUE * (this->getd2HighThrust_1D_Mdot_dP2() - this->getd2HighIsp_1D_Mdot_dP2()) + this->getd2HighIsp_1D_Mdot_dP2();

            //Step 3: evaluate the 2D polynomial
            this->polyThrust_2D = 0.0;
            this->poly_dThrust_dP_2D = 0.0;
            this->poly_dThrust_du_Mdot_2D = 0.0;

            doubleType powers_of_P[6];
            doubleType powers_of_commanded_Mdot[6];
            powers_of_P[0] = 0;
            powers_of_commanded_Mdot[0] = 0;
            powers_of_P[1] = 1;
            powers_of_commanded_Mdot[1] = 1;

            for (size_t i = 2; i < 6; ++i)
            {
                powers_of_P[i] = powers_of_P[i - 1] * inputPower;
                powers_of_commanded_Mdot[i] = powers_of_commanded_Mdot[i - 1] * this->polyMdot_2D;
            }

            for (size_t i = 0; i < 5; ++i)
            {
                for (size_t j = 0; j < 5; ++j)
                {
                    this->polyThrust_2D += this->poly2D_coefficients[i][j] * powers_of_P[i + 1] * powers_of_commanded_Mdot[j + 1];
                }
            }

            //Step 4: evaluate the derivatives of the 2D polynomial
            for (size_t i = 0; i < 5; ++i)
            {
                for (size_t j = 0; j < 5; ++j)
                {
                    this->poly_dThrust_dP_2D += this->poly2D_coefficients[i][j] * (powers_of_P[i] * i * powers_of_commanded_Mdot[j + 1] 
                        + powers_of_P[i + 1] * j * powers_of_commanded_Mdot[j] * this->poly_dMdot_dP_2D)_GETVALUE;

                    this->poly_dThrust_du_Mdot_2D += (this->poly2D_coefficients[i][j] * powers_of_P[i + 1] * j * powers_of_commanded_Mdot[j] * this->poly_dMdot_du_Mdot_2D)_GETVALUE;
                }
            }
        }

        doubleType ThrottleTable::get1D_poly_Thrust() const
        {
            return this->polyThrust_1D;
        }

        doubleType ThrottleTable::get1D_poly_Mdot() const
        {
            return this->polyMdot_1D;
        }

        doubleType ThrottleTable::get2D_poly_Thrust() const
        {
            return this->polyThrust_2D;
        }

        doubleType ThrottleTable::get2D_poly_Mdot() const
        {
            return this->polyMdot_2D;
        }

        double ThrottleTable::get2D_poly_dThrust_dP() const
        {
            return this->poly_dThrust_dP_2D;
        }

        double ThrottleTable::get2D_poly_dThrust_duMdot() const
        {
            return this->poly_dThrust_du_Mdot_2D;
        }

        double ThrottleTable::get2D_poly_dMdot_dP() const
        {
            return this->poly_dMdot_dP_2D;
        }

        double ThrottleTable::get2D_poly_dMdot_duMdot() const
        {
            return this->poly_dMdot_du_Mdot_2D;
        }

        double ThrottleTable::get1D_poly_dThrust_dP() const
        {
            return this->poly_dThrust_dP_1D;
        }

        double ThrottleTable::get1D_poly_dMdot_dP() const
        {
            return this->poly_dMdotdP_1D;
        }

        double ThrottleTable::getd2HighThrust_1D_Mdot_dP2() const
        {
            return d2HighThrust_1D_Mdot_dP2;
        }

        double ThrottleTable::getd2HighIsp_1D_Mdot_dP2() const
        {
            return d2HighIsp_1D_Mdot_dP2;
        }

        void ThrottleTable::printThrottleSets(const std::string& ThrottleOutputFileName) const
        {
            std::ofstream ThrottleOutputFile(ThrottleOutputFileName.c_str(), std::ios::trunc);
            ThrottleOutputFile << "PPU efficiency, " << this->PPUefficiency << std::endl;
            ThrottleOutputFile << "PPU min power(kW), " << this->PPUminpower << std::endl;
            ThrottleOutputFile << "PPU max power(kW), " << this->PPUmaxpower << std::endl;
            ThrottleOutputFile << std::endl;
            ThrottleOutputFile << std::endl;
            ThrottleOutputFile << "#high-Isp set" << std::endl;
            ThrottleOutputFile << "#Throttle level,Mass flow rate (mg/s),Beam Current (A),Beam Voltage (V),Thrust (mN),Isp (s),Efficiency,Thruster input power (kW)" << std::endl;
            for (ThrottleSetting setting : this->HighIspSet)
                setting.printThrottleSetting(ThrottleOutputFile, false);

            ThrottleOutputFile << std::endl;
            ThrottleOutputFile << std::endl;
            ThrottleOutputFile << "#high-Thrust set" << std::endl;
            ThrottleOutputFile << "#Throttle level,Mass flow rate (mg/s),Beam Current (A),Beam Voltage (V),Thrust (mN),Isp (s),Efficiency,Thruster input power (kW)" << std::endl;
            for (ThrottleSetting setting : this->HighThrustSet)
                setting.printThrottleSetting(ThrottleOutputFile, false);

            ThrottleOutputFile << std::endl;
            ThrottleOutputFile << std::endl;
            ThrottleOutputFile << "#high-Efficiency set" << std::endl;
            ThrottleOutputFile << "#Throttle level,Mass flow rate (mg/s),Beam Current (A),Beam Voltage (V),Thrust (mN),Isp (s),Efficiency,Thruster input power (kW)" << std::endl;
            for (ThrottleSetting setting : this->HighEfficiencySet)
                setting.printThrottleSetting(ThrottleOutputFile, false);

            ThrottleOutputFile << std::endl;
            ThrottleOutputFile << std::endl;
            ThrottleOutputFile << "#low-Mdot set" << std::endl;
            ThrottleOutputFile << "#Throttle level,Mass flow rate (mg/s),Beam Current (A),Beam Voltage (V),Thrust (mN),Isp (s),Efficiency,Thruster input power (kW)" << std::endl;
            for (ThrottleSetting setting : this->LowMdotSet)
                setting.printThrottleSetting(ThrottleOutputFile, false);

            ThrottleOutputFile << std::endl;
            ThrottleOutputFile << std::endl;
            ThrottleOutputFile << "#Full set" << std::endl;
            ThrottleOutputFile << "#Throttle level,Mass flow rate (mg/s),Beam Current (A),Beam Voltage (V),Thrust (mN),Isp (s),Efficiency,Thruster input power (kW), Mdot_on (mg/s), Mdot_off (mg/s), V_on (V), V_off (V)" << std::endl;
            for (ThrottleSetting setting : this->FullSet)
                setting.printThrottleSetting(ThrottleOutputFile, false);

            ThrottleOutputFile.close();
        }
        
        void ThrottleTable::set2DControlType(const bool & control_type)
        {
            this->control_type_2d = control_type;
        }
    }//end namespace PropulsionSystem
}//end namespace EMTG
