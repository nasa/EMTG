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

#Target file generator
#ingests .emtg files and generates target specs
#intended for use in validating EMTG and for Monte-Carlo analysis
#Jacob Englander 4-16-2018

import sys
sys.path.append('c:/emtg/PyEMTG/')
import Mission
import MissionOptions
import spiceypy
import numpy

class TargetLineGenerator:
    def __init__(self, 
                 myEvent, \
                 SampleNumber = -1, \
                 CentralBody = 'Fuzzy'):
        
        self.SampleNumber = -1
        self.CentralBody = 'Fuzzy'
        self.Frame = 'EME2000'
        self.EpochETJulian = 0.0
        self.EpochETGregorian = 'tomorrow'
        self.X = 0.0
        self.Y = 0.0
        self.Z = 0.0
        self.VX = 0.0
        self.VY = 0.0
        self.VZ = 0.0
        self.BR = 0.0
        self.BT = 0.0
        self.Mass = 0.0

        self.parse_EMTG_event(SampleNumber, CentralBody, myEvent)

    def parse_EMTG_event(self, SampleNumber, CentralBody, myEvent):        
        self.SampleNumber = SampleNumber
        self.CentralBody = CentralBody
        self.Frame = 'EMO2000'
        self.EpochETJulian = myEvent.JulianDate
        self.EpochETGregorian = spiceypy.timout(spiceypy.unitim(myEvent.JulianDate, 'JDTDB', 'TDB'), 'YYYY Mon DD ::TDB HR:MN:SC')
        self.Mass = myEvent.Mass
        
        if 'Sun' in self.CentralBody:
            self.X  = myEvent.SpacecraftState[0]
            self.Y  = myEvent.SpacecraftState[1]
            self.Z  = myEvent.SpacecraftState[2]
            self.VX = myEvent.SpacecraftState[3]
            self.VY = myEvent.SpacecraftState[4]
            self.VZ = myEvent.SpacecraftState[5]
            
        elif 'Earth' in self.CentralBody:
            self.X  = myEvent.SpacecraftState[0]
            self.Y  = myEvent.SpacecraftState[1]
            self.Z  = myEvent.SpacecraftState[2]
            self.VX = myEvent.SpacecraftState[3]
            self.VY = myEvent.SpacecraftState[4]
            self.VZ = myEvent.SpacecraftState[5]
            
            rvec = numpy.array(myEvent.SpacecraftState[:3])
            vvec = numpy.array(myEvent.SpacecraftState[3:])
            
            r = numpy.linalg.norm(rvec)
            v = numpy.linalg.norm(vvec)
                        
            h = numpy.cross(rvec,vvec)

            hhat = h / numpy.linalg.norm(h)
            
            r_dot_v = numpy.dot(rvec,vvec)
            
            mu = 398600.435436096
            
            evec = 1.0 / mu * ((v*v - mu / r)*rvec - r_dot_v * vvec)
            
            e = numpy.linalg.norm(evec)
            
            Beta = numpy.arccos(1.0/e)
            
            h_hat_cross_e = numpy.cross(hhat,evec)
            
            #compute b-plane coordinates
            Shat = numpy.cos(Beta)*evec / e + numpy.sin(Beta) * h_hat_cross_e / numpy.linalg.norm(h_hat_cross_e)
            
            khat = numpy.array([0.0, 0.0, -1.0])
            
            T = numpy.cross(Shat,khat)
            
            That = T / numpy.linalg.norm(T)

            Rhat = numpy.cross(Shat, That)

            Bhat = numpy.cross(Shat, hhat)
            
            a = 1.0 / (2.0 / r - v*v / mu)
            
            b = a * numpy.sqrt(e*e - 1.0)
            
            B = Bhat * b
            
            self.BR = numpy.dot(B, Rhat)
            self.BT = numpy.dot(B, That)
            
        else:
            #if the central body is not the Sun then the we are whacking a target
            self.X  = 0.0
            self.Y  = 0.0
            self.Z  = 0.0

            #v-infinity is opposite delta-v
            self.VX = -myEvent.DeltaVorThrustVectorControl[0]
            self.VY = -myEvent.DeltaVorThrustVectorControl[1]
            self.VZ = -myEvent.DeltaVorThrustVectorControl[2]

            #compute the encounter state in the target frame
            encounter_radius = 1000.0
            if ('Polymele' in self.CentralBody):
                encounter_radius = 399.0

            Spacecraft_State_Target_Frame = numpy.array([0.0, 0.0, encounter_radius])

            #compute the basis vectors of the target frame
            Vinfinity = numpy.array([self.VX, self.VY, self.VZ])
            VectorToSun = numpy.array([-myEvent.SpacecraftState[0], -myEvent.SpacecraftState[1], -myEvent.SpacecraftState[2]])

            Xhat = Vinfinity / numpy.linalg.norm(Vinfinity)
            Yhat = numpy.cross(VectorToSun, Xhat) / numpy.linalg.norm(numpy.cross(VectorToSun, Xhat))
            Zhat = numpy.cross(Xhat, Yhat)

            R_from_Ecliptic_to_Target = numpy.matrix([Xhat.T, Yhat.T, Zhat.T])

            R_from_Target_to_Ecliptic = R_from_Ecliptic_to_Target.T

            #rotate the target state into the ecliptic frame
            Spacecraft_State_Ecliptic_Frame = R_from_Target_to_Ecliptic * numpy.matrix(Spacecraft_State_Target_Frame).T

            #set the output states
            self.X = Spacecraft_State_Ecliptic_Frame[0]
            self.Y = Spacecraft_State_Ecliptic_Frame[1]
            self.Z = Spacecraft_State_Ecliptic_Frame[2]

            #compute b-plane coordinates
            Shat = Vinfinity / numpy.linalg.norm(Vinfinity)
            hhat = numpy.cross(numpy.array(Spacecraft_State_Ecliptic_Frame.T), Shat) / encounter_radius

            Bhat = numpy.cross(Shat, hhat)

            khat = numpy.array([0.0, 0.0, 1.0])

            That = numpy.cross(Shat, khat)

            Rhat = numpy.cross(Shat, That)

            Bmagnitude = numpy.linalg.norm(numpy.cross(Shat, numpy.cross(numpy.array(Spacecraft_State_Ecliptic_Frame.T), Shat)))
            self.BR = Bmagnitude * numpy.dot(Bhat, Rhat)
            self.BT = Bmagnitude * numpy.dot(Bhat, That)

    def write_record(self, outputfile):
        outputfile.write(str(self.SampleNumber) + ',')
        outputfile.write(self.CentralBody + ',')
        outputfile.write(self.Frame + ',')
        outputfile.write(self.EpochETGregorian + ',')
        outputfile.write('%15.15f' % self.X + ',')
        outputfile.write('%15.15f' % self.Y + ',')
        outputfile.write('%15.15f' % self.Z + ',')
        outputfile.write('%15.15f' % self.VX + ',')
        outputfile.write('%15.15f' % self.VY + ',')
        outputfile.write('%15.15f' % self.VZ + ',')
        outputfile.write('%15.15f' % self.Mass + ',')
        outputfile.write('%15.15f' % self.BR + ',')
        outputfile.write('%15.15f' % self.BT + '\n')

#basic plan
#loop over .emtg files in working directory
#for each file loop over maneuvers (chem_burn, departure with an 'impulse' entry under thrust)
#keep count of the number of maneuvers found. If you have one more maneuver that you have entries in TargetFileObjectList, then make another one
#every time you hit a maneuver, write roll ahead to the next thing that isn't a coast or a match point and print that
class TargetFileGenerator:
    def __init__(self, working_directory = '.'):
        self.TargetFileObject = None
        self.working_directory = working_directory
        
    def harvest_and_print(self):
        from os import walk, path


        self.TargetFileObject = open(self.working_directory + '/Lucy_HighFidelityPrimary_maxmass.mission_target_spec', 'w')
        self.TargetFileObject.write('<SAMPLE>,<CBODY>,<FRAME>,<EPOCH(ET)>,<X[km]>,<Y[km]>,<Z[km]>,<VX[km/s]>,<VY[km/s]>,<VZ[km/s]>,<MASS[kg]>,<B.R[km]>,<B.T [km]>\n')
        
        
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
            targetIndex = 0

            #print('Writing target specs for Sample #' + str(sampleIndex) + ' \n')

            for journey in myMission.Journeys:
                for event in journey.missionevents:
                    if event.EventType == 'chem_burn' or (event.EventType == 'departure' and event.DVmagorThrottle > 0.0):
                        #increment the target index
                        targetIndex += 1

                        #fast-forward to the next non-coast event
                        found_target = 0
                        for target_journey in myMission.Journeys:
                            for target_event in target_journey.missionevents:
                                if target_event.EventType not in ['coast','match_point','departure'] and target_event.JulianDate > event.JulianDate:

                                    #did we already have a file for this maneuver? If not, make one
                                    # if targetIndex > len(self.TargetFileObjectList):
     #                                    self.TargetFileObjectList.append(open(self.working_directory + '/Target' + str(targetIndex) + '.targetspec', 'w'))
     #                                    self.TargetFileObjectList[targetIndex - 1].write('<SAMPLE>,<CBODY>,<FRAME>,<EPOCH(ET)>,<X[km]>,<Y[km]>,<Z[km]>,<VX[km/s]>,<VY[km/s]>,<VZ[km/s]>,<B.R[km]>,<B.T [km]>\n')

                                    CentralBody = ''

                                    if target_event.Location in ['deep-space', 'free point'] and 'Sun' in target_journey.central_body:
                                        CentralBody = 'Sun'
                                    elif target_event.Location in ['deep-space', 'free point'] and 'Earth' in target_journey.central_body:
                                        CentralBody = 'Earth'
                                    elif target_event.Location == 'Earth_BE':
                                        continue
                                    else:
                                        CentralBody = target_event.Location

                                    myTargetLineGenerator = TargetLineGenerator(target_event, myMission.SampleIndex, CentralBody)

                                    myTargetLineGenerator.write_record(self.TargetFileObject)
                                    
                                    found_target = 1
                                    break
                            if found_target:
                                break

        self.TargetFileObject.close()
            # fileObject.close()