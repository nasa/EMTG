#simple SEP Monte Carlo

import EphemerisFileReader
import ConOpsPeriod
import copy
import os

import missionoptions

class SimpleSEPMonteCarloStep(object):
    """single step of a Monte-Carlo - we keep doing these until we're done"""

    def __init__(self, propulatorPath = 'c:/emtg/bin/Propulator.exe', missionoptionsFileName = "default.emtgopt", ephemerisFileName = "default.ephemeris", workingDirectory = "./", stepIndex = 0, sampleIndex = 0, journeyIndex = 0):
        self.propulatorPath = propulatorPath
        self.missionoptionsFileName = missionoptionsFileName
        self.ephemerisFileName = ephemerisFileName
        self.workingDirectory = workingDirectory
        self.stepIndex = stepIndex
        self.sampleIndex = sampleIndex
        self.journeyIndex = journeyIndex

        self.FirstConOpsPeriod = None
        self.SecondConOpsPeriod = None

        self.missionoptions = missionoptions.missionoptions(self.missionoptionsFileName)

        self.myEphemerisFileReader = EphemerisFileReader.EphemerisFileReader(self.ephemerisFileName)

    def execute(self):
        #Step 1: Read the first conops period from the ephemeris file. If it doesn't exist then we've reached the end of the mission, so stop!
        self.FirstConOpsPeriod = self.myEphemerisFileReader.getFirstConopsPeriod()

        if self.FirstConOpsPeriod == None:
            #we're done!
            return "nothing to perturb"

        #Step 2: Apply maneuver execution error
        #Step 2.1: copy the first conops period into a "perturbed" version
        self.perturbedFirstConOpsPeriod = copy.deepcopy(self.FirstConOpsPeriod)
        
        #Step 2.2: here there be perturbing, from Noble
        self.perturbControl(self.perturbedFirstConOpsPeriod)

        #Step 2.3: propulate the first conops period
        self.propulate(self.perturbedFirstConOpsPeriod)

        #Step 3: Read the second conops period from the ephemeris file. If it doesn't exist then we don't have to re-optimize, just apply execution error and see where we are at the end
        self.SecondConOpsPeriod = self.myEphemerisFileReader.getSecondConopsPeriod()

        if self.SecondConOpsPeriod == None:
            #we're done!
            return "nothing to optimize"

        #Step 4: create a new EMTG problem and optimize it
        output_path = "the place where we put the outputs"

        #Step 5: return the path to the EMTG outputs
        return output_path

    def perturbControl(self, ConOpsPeriod):
        #do stuff with the ConOpsPeriod
        print("waiting for stuff from Noble")

    def propulate(self, ConOpsPeriod):
        #Step 1: write a propulator input file for the thrust period
        DutyCycle = 1.0
        if self.missionoptions.Journeys[self.journeyIndex].override_duty_cycle:
            DutyCycle = self.missionoptions.Journeys[self.journeyIndex].duty_cycle
        else:
            DutyCycle = self.missionoptions.engine_duty_cycle

        stepSize = 86400.0
        if self.missionoptions.Journeys[self.journeyIndex].override_integration_step_size:
            stepSize = self.missionoptions.Journeys[self.journeyIndex].integration_step_size
        else:
            stepSize = self.missionoptions.integration_time_step_size

        inputfilepath = self.workingDirectory + '/' + 'PropulatorInput_Sample' + str(self.sampleIndex) + '_Step' + str(self.stepIndex) + '.csv'
        inputfile = open(inputfilepath, 'w')
        inputfile.write('journeyIndex, ' + str(self.journeyIndex) + '\n')
        inputfile.write('DutyCycle, ' + str(DutyCycle) + '\n')
        inputfile.write('stepSize, ' + str(stepSize) + '\n')
        inputfile.write('\n')
        inputfile.write('#JD, x [km/s], y [km/s], z [km/s], xdot [km/s], ydot [km/s], zdot [km/s], mass [kg]\n')
        inputfile.write(str(ConOpsPeriod.StateBeforeThrusting[7]))
        for entry in ConOpsPeriod.StateBeforeThrusting[0:7]:
            inputfile.write(', ' + str(entry))
        for entry in ConOpsPeriod.ControlVector:
            inputfile.write(', ' + str(entry))
        inputfile.write(str(ConOpsPeriod.ThrustTime))
        inputfile.close()

        #Step 2: propulate with control to the end of the thrust period
        outputfilepath = self.workingDirectory + '/' + 'PropulatorOutput_Sample' + str(self.sampleIndex) + '_Step' + str(self.stepIndex) + '.csv'
        os.system(self.propulatorPath + ' ' + self.missionoptionsFilePath + ' ' + inputfilepath + ' ' + outputfilepath)

        #Step 3: read the propulator output file and pull the state at the end of the thrust period
        outputfile = open(outputfilepath, 'r')
        #assume that the interesting stuff is on the sixth line
        lineIndex = 0
        for line in outputfile:
            if lineIndex < 5:
                lineIndex += 1
            else:
                linecell = line.split(',')
                ConOpsPeriod.StateAfterThrusting[0] = float(linecell[7])
                for i in range(0,7):
                    ConOpsPeriod.StateAfterThrusting[i + 1] = float(linecell[i])
        outputfile.close()

        #Step 4: write a propulator input file for the coast period
        inputfile = open(inputfilepath, 'w')
        inputfile.write('journeyIndex,' + str(self.journeyIndex) + '\n')
        inputfile.write('DutyCycle,' + str(DutyCycle) + '\n')
        inputfile.write('stepSize, ' + str(stepSize) + '\n')
        inputfile.write('\n')
        inputfile.write('#JD, x [km/s], y [km/s], z [km/s], xdot [km/s], ydot [km/s], zdot [km/s], mass [kg]\n')
        inputfile.write(str(ConOpsPeriod.StateAfterThrusting[7]))
        for entry in ConOpsPeriod.StateAfterThrusting[0:7]:
            inputfile.write(', ' + str(entry))
        inputfile.write(str(ConOpsPeriod.CoastTime))
        inputfile.close()

        #Step 5: propulate without control to the end of the coast period
        os.system(self.propulatorPath + ' ' + self.missionoptionsFilePath + ' ' + inputfilepath + ' ' + outputfilepath)

        #Step 6: read the propulator output file and pull the state at the end of the coast period
        outputfile = open(outputfilepath, 'r')
        #assume that the interesting stuff is on the sixth line
        lineIndex = 0
        for line in outputfile:
            if lineIndex < 5:
                lineIndex += 1
            else:
                linecell = line.split(',')
                ConOpsPeriod.StateAfterCoasting[0] = float(linecell[7])
                for i in range(0,7):
                    ConOpsPeriod.StateAfterCoasting[i + 1] = float(linecell[i])
        outputfile.close()