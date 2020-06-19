#EMTG: Evolutionary Mission Trajectory Generator
#An open-source global optimization tool for preliminary mission design
#Provided by NASA Goddard Space Flight Center
#
#Copyright (c) 2014 - 2018 United States Government as represented by the
#Administrator of the National Aeronautics and Space Administration.
#All Other Rights Reserved.
#
#Licensed under the NASA Open Source License (the "License"); 
#You may not use this file except in compliance with the License. 
#You may obtain a copy of the License at:
#https://opensource.org/licenses/NASA-1.3
#Unless required by applicable law or agreed to in writing, software
#distributed under the License is distributed on an "AS IS" BASIS,
#WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either 
#express or implied.   See the License for the specific language
#governing permissions and limitations under the License.

#base class for high-fidelity trajectory driver
#Jacob Englander 1-23-2018

import Convert_TwoPointShootingLowThrust_to_PSFB

#filePath = 'C:/Projects/CAESAR/PhaseA/EMTGv9/HighFidelity/CAESAR_HighFidelityTest_FBLT_perturbed_392018_111229'
filePath = 'C:/emtg/testatron/tests/physics_options' #'C:/Projects/CAESAR/PhaseA/David_Alternate_Inbound/Mars/HighFidelityV2/intermediate/'

#originalMissionFile = 'CAESAR_HighFidelityTest_FBLT_perturbed.emtg'
#originalOptionsFile = 'CAESAR_HighFidelityTest_FBLT_perturbed.emtgopt'
originalMissionFile = 'physicsoptions_PSFBSphericalRADEC.emtg'    #'TradeStudy_Case52_seeded_by_External_TradeStudy_Case52_HighFidelity.emtg'
originalOptionsFile = 'physicsoptions_PSFBSphericalRADEC.emtgopt' #'TradeStudy_Case52_seeded_by_External_TradeStudy_Case52_HighFidelity.emtgopt'
output_path = 'C:/emtg/testatron/tests/physics_options/Test' #'C:/Projects/CAESAR/PhaseA/David_Alternate_Inbound/Mars/HighFidelityV2/'

myPSFB_Converter = Convert_TwoPointShootingLowThrust_to_PSFB.MissionConverter_TPSLT_to_PSFB(filePath + '/' + originalMissionFile,
                                                                                            filePath + '/' + originalOptionsFile,
                                                                                            ['MGALT', 2])

myPSFB_Converter.CreateMission()

myPSFB_Converter.CreateInitialGuess()

newOptions = myPSFB_Converter.getMissionOptions()

newOptions.trialX = myPSFB_Converter.CreateInitialGuess()

newOptions.short_output_file_names = 1

newOptions.mission_type = 9 #variable phase type

newOptions.mission_name = originalMissionFile.replace('.emtg','').replace('MGALT','PSFB')
newOptions.write_options_file(output_path + '/' + newOptions.mission_name + '.emtgopt')
