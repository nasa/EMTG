#throttle report generator
#ingests .ephemeris files that include thrust, mass flow, and throttle level columns
#creates throttle level reports for propulsion system analysis
#Jacob Englander 10-10-2018


class ThrottleTableReport(object):
    def __init__(self, ephemerisfilename = None, optionsfilename = None, number_of_available_thrusters = 1):
        self.__multiThrusterThroughputPerThrottleLevel__ = {}
        self.__multiThrusterOperatingHoursPerThrottleLevel__ = {}
        self.__multiThrusterTotalPropellantUsed__ = 0.0
        self.__multiThrusterTotalOperatingHours__ = 0.0
        self.__totalThroughtputPerThrottleLevel__ = {}
        self.__totalOperatingHoursPerThrottleLevel__ = {}
        self.__worst_case_singleThrusterThroughput__ = {}
        self.__worst_case_operatingHoursPerThruster__ = {}
        self.__best_case_singleThrusterThroughput__ = {}
        self.__best_case_operatingHoursPerThruster__ = {}
        self.__totalPropellantUsed__ = 0.0
        self.__totalOperatingHours__ = 0.0
        self.__single_thruster_worst_case_throughput__ = 0.0
        self.__single_thruster_worst_case_operating_hours__ = 0.0
        self.__single_thruster_best_case_throughput__ = 0.0
        self.__single_thruster_best_case_operating_hours__ = 0.0
        self.__total_deltav__ = 0.0
        self.__number_of_available_thrusters__ = int(number_of_available_thrusters)

        if not ephemerisfilename == None:
            self.parseEphemerisFile(ephemerisfilename, optionsfilename)


    def parseEphemerisFile(self, ephemerisfilename, optionsfilename):
        import EMTG_ephemeris_reader
        import sys
        sys.path.append('c:/emtg/PyEMTG/')
        import MissionOptions
        from math import log

        myOptions = MissionOptions.MissionOptions(optionsfilename)

        leapsecondspath = myOptions.universe_folder + '/ephemeris_files/' + myOptions.SPICE_leap_seconds_kernel

        myReader = EMTG_ephemeris_reader.EMTG_ephemeris_reader(ephemerisfilename, leapsecondspath)

        for recordIndex in range(0, myReader.getLength()):
            
            myRecord = myReader.getRecord(recordIndex)

            if myRecord['ThrustMagnitude(N)'] != None:
                if myRecord['ThrottleLevel'] != 'none':

                    massExpended = float(myRecord['MassFlowRate(kg/s)']) * myRecord['timeWidth']

                    Isp = float(myRecord['Isp(s)'])
                    startMass = float(myRecord['mass(kg)'])

                    nThrusters = 1
                    throttleLevel = myRecord['ThrottleLevel']
                    if 'x' in  myRecord['ThrottleLevel']:
                        stringcell =  myRecord['ThrottleLevel'].split('x')
                        throttleLevel = stringcell[1]
                        nThrusters = int(stringcell[0])

                    #print('startMass = ', startMass, ', massExpended = ', massExpended, ', ratio = ', startMass / (startMass - massExpended))
                    self.__total_deltav__ += Isp * 9.80665 / 1000.0 * log(startMass / (startMass - massExpended))
            
                    self.__totalPropellantUsed__ += massExpended
                    self.__totalOperatingHours__ += myRecord['timeWidth'] * nThrusters / 3600.0

                    #update total-thruster data fields
                    existingEntry = throttleLevel in self.__totalThroughtputPerThrottleLevel__
                    
                    if existingEntry:
                        self.__totalThroughtputPerThrottleLevel__[throttleLevel] += massExpended
                        self.__totalOperatingHoursPerThrottleLevel__[throttleLevel] += myRecord['timeWidth'] * nThrusters / 3600.0
                    else:
                        self.__totalThroughtputPerThrottleLevel__[throttleLevel] = massExpended
                        self.__totalOperatingHoursPerThrottleLevel__[throttleLevel] = myRecord['timeWidth'] * nThrusters / 3600.0

                    #update multi-thruster
                    existingEntry =  myRecord['ThrottleLevel'] in self.__multiThrusterThroughputPerThrottleLevel__
                    
                    self.__multiThrusterTotalPropellantUsed__ += massExpended
                    self.__multiThrusterTotalOperatingHours__ += myRecord['timeWidth'] / 3600.0
                    if existingEntry:
                        self.__multiThrusterThroughputPerThrottleLevel__[myRecord['ThrottleLevel']] += massExpended
                        self.__multiThrusterOperatingHoursPerThrottleLevel__[myRecord['ThrottleLevel']] += myRecord['timeWidth'] / 3600.0
                    else:
                        self.__multiThrusterThroughputPerThrottleLevel__[myRecord['ThrottleLevel']] = massExpended
                        self.__multiThrusterOperatingHoursPerThrottleLevel__[myRecord['ThrottleLevel']] = myRecord['timeWidth'] / 3600.0

                    ##update single-thruster data fields
                    #if 'x' in myRecord['ThrottleLevel']:
                    #    #make a new throttle level without nx and figure out what n is
                    #    singleThrusterThrottleLevel = myRecord['ThrottleLevel']
                
                    #    if 'x' in singleThrusterThrottleLevel:
                    #        stringcell = singleThrusterThrottleLevel.split('x')
                    #        singleThrusterThrottleLevel = stringcell[1]                        

                    #    #update total thruster table                
                    #    newEntryTotalThruster = singleThrusterThrottleLevel in self.__totalThroughtputPerThrottleLevel__

                    #    if newEntryTotalThruster:
                    #        self.__totalThroughtputPerThrottleLevel__[singleThrusterThrottleLevel] += massExpended
                    #        self.__totalOperatingHoursPerThrottleLevel__[singleThrusterThrottleLevel] += myRecord['timeWidth'] / 3600.0
                    #    else:
                    #        self.__totalThroughtputPerThrottleLevel__[singleThrusterThrottleLevel] = massExpended
                    #        self.__totalOperatingHoursPerThrottleLevel__[singleThrusterThrottleLevel] = myRecord['timeWidth'] / 3600.0


                    #    newEntrySingleThruster = singleThrusterThrottleLevel in self.__worst_case_singleThrusterThroughput__
                
                    #    if newEntrySingleThruster:
                    #        self.__worst_case_singleThrusterThroughput__[singleThrusterThrottleLevel] += massExpended / nThrusters
                    #        self.__worst_case_operatingHoursPerThruster__[singleThrusterThrottleLevel] += myRecord['timeWidth'] * nThrusters / 3600.0 / nThrusters
                    #    else:
                    #        self.__worst_case_singleThrusterThroughput__[singleThrusterThrottleLevel] = massExpended / nThrusters
                    #        self.__worst_case_operatingHoursPerThruster__[singleThrusterThrottleLevel] = myRecord['timeWidth'] * nThrusters / 3600.0 / nThrusters

                        
                    #    newEntrySingleThrusterBestCase = singleThrusterThrottleLevel in self.__best_case_singleThrusterThroughput__
                    #    if newEntrySingleThrusterBestCase:
                    #        self.__best_case_singleThrusterThroughput__[singleThrusterThrottleLevel] += massExpended / self.__number_of_available_thrusters__
                    #        self.__best_case_operatingHoursPerThruster__[singleThrusterThrottleLevel] += myRecord['timeWidth'] * nThrusters / 3600.0 / self.__number_of_available_thrusters__
                    #    else:
                    #        self.__best_case_singleThrusterThroughput__[singleThrusterThrottleLevel] = massExpended / self.__number_of_available_thrusters__
                    #        self.__best_case_operatingHoursPerThruster__[singleThrusterThrottleLevel] = myRecord['timeWidth'] * nThrusters / 3600.0 / self.__number_of_available_thrusters__

                    #    self.__single_thruster_worst_case_throughput__ += massExpended / nThrusters
                    #    self.__single_thruster_worst_case_operating_hours__ += myRecord['timeWidth'] * nThrusters / 3600.0 / nThrusters
                    #    self.__single_thruster_best_case_throughput__ += massExpended / self.__number_of_available_thrusters__
                    #    self.__single_thruster_best_case_operating_hours__ += myRecord['timeWidth'] * nThrusters / 3600.0 / self.__number_of_available_thrusters__
                    #else:
                    
                    #    #update total thruster table                
                    #    newEntryTotalThruster = myRecord['ThrottleLevel'] in self.__totalThroughtputPerThrottleLevel__

                    #    if newEntryTotalThruster:
                    #        self.__totalThroughtputPerThrottleLevel__[myRecord['ThrottleLevel']] += massExpended
                    #        self.__totalOperatingHoursPerThrottleLevel__[myRecord['ThrottleLevel']] += myRecord['timeWidth'] / 3600.0
                    #    else:
                    #        self.__totalThroughtputPerThrottleLevel__[myRecord['ThrottleLevel']] = massExpended
                    #        self.__totalOperatingHoursPerThrottleLevel__[myRecord['ThrottleLevel']] = myRecord['timeWidth'] / 3600.0


                    #    newEntrySingleThruster = myRecord['ThrottleLevel'] in self.__worst_case_singleThrusterThroughput__

                
                    #    if newEntrySingleThruster:
                    #        self.__worst_case_singleThrusterThroughput__[myRecord['ThrottleLevel']] += massExpended / (self.__number_of_available_thrusters__ - 1)
                    #        self.__worst_case_operatingHoursPerThruster__[myRecord['ThrottleLevel']] += myRecord['timeWidth'] * nThrusters / 3600.0 / (self.__number_of_available_thrusters__ - 1)
                    #    else:
                    #        self.__worst_case_singleThrusterThroughput__[myRecord['ThrottleLevel']] = massExpended / (self.__number_of_available_thrusters__ - 1)
                    #        self.__worst_case_operatingHoursPerThruster__[myRecord['ThrottleLevel']] = myRecord['timeWidth'] * nThrusters / 3600.0 / (self.__number_of_available_thrusters__ - 1)

                    #    newEntrySingleThrusterBestCase = myRecord['ThrottleLevel'] in self.__best_case_singleThrusterThroughput__
                
                    #    if newEntrySingleThrusterBestCase:
                    #        self.__best_case_singleThrusterThroughput__[myRecord['ThrottleLevel']] += massExpended / self.__number_of_available_thrusters__
                    #        self.__best_case_operatingHoursPerThruster__[myRecord['ThrottleLevel']] += myRecord['timeWidth'] * nThrusters / 3600.0 / self.__number_of_available_thrusters__
                    #    else:
                    #        self.__best_case_singleThrusterThroughput__[myRecord['ThrottleLevel']] = massExpended / self.__number_of_available_thrusters__
                    #        self.__best_case_operatingHoursPerThruster__[myRecord['ThrottleLevel']] = myRecord['timeWidth'] * nThrusters / 3600.0 / self.__number_of_available_thrusters__

                    #    self.__single_thruster_worst_case_throughput__ += massExpended / (self.__number_of_available_thrusters__ - 1)
                    #    self.__single_thruster_worst_case_operating_hours__ += myRecord['timeWidth'] * nThrusters / 3600.0 / (self.__number_of_available_thrusters__ - 1)
                    #    self.__single_thruster_best_case_throughput__ += massExpended / self.__number_of_available_thrusters__
                    #    self.__single_thruster_best_case_operating_hours__ += myRecord['timeWidth'] * nThrusters / 3600.0 / self.__number_of_available_thrusters__
        

    def printReport(self, reportfilename = 'default_throttle_table_report.csv'):

        with open(reportfilename, 'w') as file:
            #totals
            file.write('all thrusting\n')
            file.write('Throttle setting, throughput (kg), operating hours\n')
            for key in self.__totalThroughtputPerThrottleLevel__.keys():
                if key != "none":
                    file.write(key + ',')
                    file.write(str(self.__totalThroughtputPerThrottleLevel__[key]) + ',')
                    file.write(str(self.__totalOperatingHoursPerThrottleLevel__[key]) + '\n')

            file.write('Total propellant used (kg), ' + str(self.__totalPropellantUsed__) + '\n')
            file.write('Total operating hours, ' + str(self.__totalOperatingHours__) + '\n')
            file.write('\n')

            #multi-thruster
            file.write('multi-thruster\n')
            file.write('Throttle setting, throughput (kg), operating hours\n')
            for key in self.__multiThrusterThroughputPerThrottleLevel__.keys():
                if key != "none":
                    file.write(key + ',')
                    file.write(str(self.__multiThrusterThroughputPerThrottleLevel__[key]) + ',')
                    file.write(str(self.__multiThrusterOperatingHoursPerThrottleLevel__[key]) + '\n')

            file.write('Total propellant used (kg), ' + str(self.__multiThrusterTotalPropellantUsed__) + '\n')
            file.write('Total operating hours, ' + str(self.__multiThrusterTotalOperatingHours__) + '\n')
            file.write('\n')

            #single thruster operations
            file.write('worst-case on a single thruster\n')
            file.write('Throttle setting, throughput on a single thruster (kg), operating hours on a single thruster\n')
            
            #for key in self.__worst_case_singleThrusterThroughput__.keys():
            #    if key != "none":
            #        file.write(key + ',')
            #        file.write(str(self.__worst_case_singleThrusterThroughput__[key]) + ',')
            #        file.write(str(self.__worst_case_operatingHoursPerThruster__[key]) + '\n')
            #file.write('\n')
            #file.write('Worst case throughput on a single thruster (kg), ' + str(self.__single_thruster_worst_case_throughput__) + '\n')
            #file.write('Worst case operating hours on a single thruster, ' + str(self.__single_thruster_worst_case_operating_hours__) + '\n')
            #file.write('\n')

            ##single thruster operations
            #file.write('best-case on a single thruster\n')
            #file.write('Throttle setting, throughput on a single thruster (kg), operating hours on a single thruster\n')
            
            #for key in self.__best_case_singleThrusterThroughput__.keys():
            #    if key != "none":
            #        file.write(key + ',')
            #        file.write(str(self.__best_case_singleThrusterThroughput__[key]) + ',')
            #        file.write(str(self.__best_case_operatingHoursPerThruster__[key]) + '\n')
            #file.write('\n')
            #file.write('Best case throughput on a single thruster (kg), ' + str(self.__single_thruster_best_case_throughput__) + '\n')
            #file.write('Best case operating hours on a single thruster, ' + str(self.__single_thruster_best_case_operating_hours__) + '\n')
            #file.write('\n')

            file.write('\n')
            file.write('Total delta-v (km/s), ' + str(self.__total_deltav__) +'\n')

if __name__ == '__main__':
    from sys import argv
    
    if len(argv) != 5:
        raise Exception('Syntax is "ThrottleTableReport ephemerisFileName optionsFileName reportFileName number_of_thrusters"')

    myThrottleTableReport = ThrottleTableReport(argv[1], argv[2], argv[4])

    myThrottleTableReport.printReport(argv[3])