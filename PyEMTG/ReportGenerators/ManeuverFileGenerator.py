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

#Maneuver file generator
#ingests .emtg files and generates maneuver specs
#intended for use in validating EMTG and for Monte-Carlo analysis
#only works for chemical events for now
#Jacob Englander 4-3-2018

import sys
sys.path.append('/Users/jknittel/EMTG/PyEMTG/')
import Mission
import MissionOptions
import spiceypy

from math import exp

class ManeuverLineGenerator:

    def __init__(self, 
                 myEvent, \
                 SampleNumber = -1, \
                 CentralBody = 'Fuzzy', \
                 ThrustMagnitude = 0.0):

        self.SampleNumber = -1
        self.CentralBody = 'Fuzzy'
        self.Frame = 'EME2000'
        self.EpochETJulian = 0.0
        self.EpochETGregorian = 'tomorrow'
        self.ThrustX = 0.0
        self.ThrustY = 0.0
        self.ThrustZ = 0.0
        self.ThrustMagnitude = 0.0
        self.StartMass = 0.0
        self.Mdot = 0.0
        self.DutyCycle = 1.0
        self.FinalMass = 0.0
        self.DVmagnitude = 0.0
        self.Duration = 0.0

        self.parse_EMTG_event(SampleNumber, CentralBody, myEvent, ThrustMagnitude)

    def parse_EMTG_event(self, SampleNumber, CentralBody, myEvent, ThrustMagnitude):
        self.SampleNumber = SampleNumber
        self.CentralBody = CentralBody
        self.Frame = 'EMO2000'
        self.EpochETJulian = myEvent.JulianDate
        self.EpochETGregorian = spiceypy.timout(spiceypy.unitim(myEvent.JulianDate, 'JDTDB', 'TDB'), 'YYYY Mon DD ::TDB HR:MN:SC')
                
        self.FinalMass = myEvent.Mass
        self.DVmagnitude = myEvent.DVmagorThrottle
        self.StartMass = self.FinalMass * exp(self.DVmagnitude * 1000.0 / myEvent.Isp / 9.80665)

        self.Mdot = ThrustMagnitude / 9.80665 / myEvent.Isp
        self.Duration = (self.StartMass - self.FinalMass) / self.Mdot

        
        self.ThrustX = 0.0
        self.ThrustY = 0.0
        self.ThrustZ = 0.0

        if self.DVmagnitude > 1.0e-25:
            self.ThrustX = myEvent.DeltaVorThrustVectorControl[0] / self.DVmagnitude
            self.ThrustY = myEvent.DeltaVorThrustVectorControl[1] / self.DVmagnitude
            self.ThrustZ = myEvent.DeltaVorThrustVectorControl[2] / self.DVmagnitude

        self.ThrustMagnitude = ThrustMagnitude

    def write_record(self, outputfile):
        outputfile.write(str(self.SampleNumber) + ',')
        outputfile.write("1,")
        outputfile.write(self.Frame + ',')
        outputfile.write(self.EpochETGregorian + ',')
        outputfile.write('%15.15f' % self.ThrustX + ',')
        outputfile.write('%15.15f' % self.ThrustY + ',')
        outputfile.write('%15.15f' % self.ThrustZ + ',')
        outputfile.write('%15.15f' % self.ThrustMagnitude + ',')
        outputfile.write('%15.15f' % self.StartMass + ',')
        outputfile.write('%15.15f' % self.Mdot + ',')
        outputfile.write('%15.15f' % self.DutyCycle + ',')
        outputfile.write('%15.15f' % self.FinalMass + ',')
        outputfile.write('%15.15f' % self.DVmagnitude + ',')
        outputfile.write('%15.15f' % self.Duration + '\n')


#basic plan
#loop over .emtg files in working directory
#for each file loop over maneuvers (chem_burn, departure with an 'impulse' entry under thrust)
#keep count of the number of maneuvers found. If you have one more maneuver that you have entries in ManeuverFileObjectList, then make another one
#every time you hit a maneuver, write it out to the appropriate file
class ManeuverFileGenerator:
    def __init__(self, working_directory = '.'):
        self.ManeuverFileObject = None
        self.working_directory = working_directory
        
    def harvest_and_print(self):
        from os import walk, path

        self.ManeuverFileObject = open(self.working_directory + '/Lucy_HighFidelityPrimary_maxmass.mission_maneuver_spec', 'w')
        self.ManeuverFileObject.write('<SAMPLE>,<NUMBER_OF_MANEUVERS>,<FRAME>,<EPOCH(ET)>,<THRX>,<THRY>,<THRZ>,<THRMAG[N]>,<SMASS[kg]>,<MDOT[kg/s]>,<DUTY>,<FMASS[kg]>,<DV[km/s]>,<DUR[s]>\n')

        MissionBucket = []

        for dirpath, dirnames, filenames in walk(self.working_directory):
            for file in filenames:
                sourcefile = path.join(dirpath, file)
        
                if sourcefile.endswith('.emtg') and 'Hop' not in sourcefile and 'FAILURE' not in sourcefile:
                    myMission = Mission.Mission(sourcefile)
                    myOptions = MissionOptions.MissionOptions(sourcefile.replace('.emtg','.emtgopt'))

                    # myMission.SampleIndex = myOptions.user_data['MonteCarloSampleIndex']
                    myMission.SampleIndex = -1

                    MissionBucket.append(myMission)

        from operator import attrgetter
        MissionBucket = sorted(MissionBucket, key=attrgetter('SampleIndex'))

        for myMission in MissionBucket:
            maneuverIndex = 0

            #print('Writing maneuver specs for Sample #' + str(sampleIndex) + ' \n')

            for journey in myMission.Journeys:
                for event in journey.missionevents:
                    if event.EventType == 'chem_burn' or (event.EventType == 'departure' and event.DVmagorThrottle > 0.0):
                        #increment the maneuver index
                        # maneuverIndex += 1

                        #did we already have a file for this maneuver? If not, make one
                        # if maneuverIndex > len(self.ManeuverFileObjectList):
                            # self.ManeuverFileObjectList.append(open(self.working_directory + '/Maneuver' + str(maneuverIndex) + '.maneuverspec', 'w'))
                            # self.ManeuverFileObjectList[maneuverIndex - 1].write('<SAMPLE>,<FRAME>,<EPOCH(ET)>,<THRX>,<THRY>,<THRZ>,<THRMAG[N]>,<SMASS[kg]>,<MDOT[kg/s]>,<DUTY>,<FMASS[kg]>,<DV[km/s]>,<DUR[s]>\n')


                        #THIS IS TEMPORARY, GET REAL VALUES FROM BRIAN
                        ThrustMagnitude = 220.0 #newtons
                        myManeuverLineGenerator = ManeuverLineGenerator(event, myMission.SampleIndex, 'Sun', ThrustMagnitude)

                        myManeuverLineGenerator.write_record(self.ManeuverFileObject)

        # for fileObject in self.ManeuverFileObjectList:
        self.ManeuverFileObject.close()