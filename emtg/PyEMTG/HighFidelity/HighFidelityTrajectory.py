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

#base class for high-fidelity trajectory maker
#Jacob Englander 1-23-2018

import HighFidelityJourney
import Mission
import MissionOptions

import copy

class HighFidelityTrajectory(object):
    def __init__(self, originalMissionPath, originalOptionsPath, useSundman=False):
        self.originalMission = Mission.Mission(originalMissionPath)
        self.originalOptions = MissionOptions.MissionOptions(originalOptionsPath)

        self.newOptions = []
        self.HighFidelityJourneys = []

        Transcriptions = [['MGALTS', 0],
                          ['FBLTS', 1],
                          ['MGALT', 2],
                          ['FBLT', 3],
                          ['PSBI', 4],
                          ['PSFB', 5],
                          ['MGAnDSMs', 6],
                          ['CoastPhase', 7],
                          ['SundmanCoastPhase', 8]]

        # if the user already has a Variable phase type mission, account for that
        if self.originalOptions.mission_type == 9:
            # scan through the journey phase types and find the "main" phase type that isn't coast phase
            # at this time we are assuming that only one phase type is allowed not counting the coast phases that may be present
            for Journey in self.originalOptions.Journeys:
                if Journey.phase_type != 7:
                    self.Transcription = Transcriptions[Journey.phase_type]
                    # go with the first non-coast phase type that is found
                    break
        else:
            self.Transcription = Transcriptions[self.originalOptions.mission_type]

        self.useSundman = useSundman

        #set parking orbit altitude over KSC
        self.originalOptions.Journeys[0].PeriapseDeparture_altitude_bounds[0] = 185.0

    def getMissionOptions(self):
        return self.newOptions

    def CreateMission(self):
        self.newOptions = copy.deepcopy(self.originalOptions)

        self.newOptions.Journeys = []

        newJourneyIndex = 0

        for journeyIndex in range(0, self.originalOptions.number_of_journeys):
            if journeyIndex < self.originalOptions.number_of_journeys - 1: #not the last journey
                self.HighFidelityJourneys.append(HighFidelityJourney.HighFidelityJourney(parent = self, 
                                                                                         UniversePath = self.originalOptions.universe_folder,
                                                                                         originalMissionOptions = self.originalOptions,
                                                                                         originalJourney = self.originalMission.Journeys[journeyIndex],
                                                                                         originalJourneyOptions = self.originalOptions.Journeys[journeyIndex],
                                                                                         nextJourney = self.originalMission.Journeys[journeyIndex + 1],
                                                                                         nextJourneyOptions = self.originalOptions.Journeys[journeyIndex + 1],
                                                                                         Transcription = self.Transcription,
                                                                                         useSundman = self.useSundman))
            else: #last journey, different call syntax                
                self.HighFidelityJourneys.append(HighFidelityJourney.HighFidelityJourney(parent = self, 
                                                                                         UniversePath = self.originalOptions.universe_folder,
                                                                                         originalJourney = self.originalMission.Journeys[journeyIndex],
                                                                                         originalMissionOptions = self.originalOptions,
                                                                                         originalJourneyOptions = self.originalOptions.Journeys[journeyIndex],
                                                                                         Transcription = self.Transcription,
                                                                                         useSundman = self.useSundman))
            
            #set the new journey indices    
            newJourneyIndex = self.HighFidelityJourneys[-1].setJourneyIndex(newJourneyIndex)
                
            for entry in self.HighFidelityJourneys[-1].CreateJourney():
                self.newOptions.Journeys.append(entry)

        #refactor mission constraints
        self.newOptions.BoundaryConstraintDefinitions = []
        self.newOptions.ManeuverConstraintDefinitions = []
        self.newOptions.PhaseDistanceConstraintDefinitions = []
        for journey in self.HighFidelityJourneys:
            journeyBoundaryConstraints = journey.createConstraints(self.originalOptions.BoundaryConstraintDefinitions)

            for entry in journeyBoundaryConstraints:
                self.newOptions.BoundaryConstraintDefinitions.append(entry)

            journeyManeuverConstraints = journey.createConstraints(self.originalOptions.ManeuverConstraintDefinitions)

            for entry in journeyManeuverConstraints:
                self.newOptions.ManeuverConstraintDefinitions.append(entry)

            journeyPhaseDistanceConstraints = journey.createConstraints(self.originalOptions.PhaseDistanceConstraintDefinitions)

            for entry in journeyPhaseDistanceConstraints:
                self.newOptions.PhaseDistanceConstraintDefinitions.append(entry)

        self.newOptions.number_of_journeys = len(self.newOptions.Journeys)
        self.newOptions.max_phases_per_journey = 1 #one phase per journey, since 3D flybys are their own journeys

        #set integrated mode appropriately
        if self.Transcription[0] in ['FBLT', 'PSFB']:
            self.newOptions.propagatorType = 1
            self.newOptions.integratorType = 1 #fixed-step

    def CreateInitialGuess(self):
        #construct the old initial guess
        originalDecisionVector = []
        for Xindex in range(0, len(self.originalMission.DecisionVector)):
            originalDecisionVector.append([self.originalMission.Xdescriptions[Xindex].replace('\n','').replace('\r',''), self.originalMission.DecisionVector[Xindex]])

        #construct the new initial guess
        self.newOptions.trialX = []
        for Journey in self.HighFidelityJourneys:
            InitialGuessSet = Journey.CreateInitialGuess(originalDecisionVector)

            for guess in InitialGuessSet:
                for entry in guess:
                    self.newOptions.trialX.append(entry)