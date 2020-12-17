#storage class for information on a conops period

import numpy as np
import MissionOptions as MO
import SpiceyPy_Utilities as SpiceyUtil

try:
    import spiceypy as spice
except:
    print("spiceypy not available")

class ConOpsPeriod(object):
    def __init__(self):
        self.StateBeforeThrusting = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]) # 8-state including mass and epoch
        self.ControlVector = np.array([0.0, 0.0, 0.0])
        self.ThrustTime = 0.0
        self.StateAfterThrusting = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]) # 8-state including mass and epoch
        self.CoastTime = 0.0
        self.StateAfterCoasting = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]) # 8-state including mass and epoch

class ConOpsDetector(object):
    def __init__(self):
        self.FirstConOpsPeriod = None
        self.SecondConOpsPeriod = None

    def extractFirstTwoConOps(self, ephemeris_file_rows, maneuver_spec_file_rows, target_spec_file_rows):
        
        # The first ConOps period is simply the first time we detect thrusting in an EMTG .ephemeris file
        # first, find the row in the .ephemeris file that corresponds to the first thrusting period
        ephem_file_index = 0
        self.FirstConOpsPeriod = ConOpsPeriod()
        ephem_file_index = self.extractConOpsPeriodEphemFile(ephemeris_file_rows, 
                                                             maneuver_spec_file_rows, 
                                                             target_spec_file_rows, 
                                                             self.FirstConOpsPeriod, 
                                                             ephem_file_index)

        # Now we need a bit of logic to detect whether or not there is a second ConOps period 
        # If there is a second ConOps cycle, then there will be an entry in the maneuver spec. file that lines up with the Julian date of the end-of-coast in the
        # first ConOps cycle. If the next thrusting entry in the maneuver spec. file is later in time, then we must be at a boundary event, there is no second ConOps
        second_conops_detected = False
        for row in maneuver_spec_file_rows:
            if np.abs(ephemeris_file_rows[ephem_file_index].julian_date - row.maneuvers[0].julian_date) < 0.01:
                second_conops_detected = True           
                self.SecondConOpsPeriod = ConOpsPeriod()
                ephem_file_index = self.extractConOpsPeriodEphemFile(ephemeris_file_rows, 
                                                                     maneuver_spec_file_rows, 
                                                                     target_spec_file_rows, 
                                                                     self.SecondConOpsPeriod, 
                                                                     ephem_file_index)
                break
        
        if not second_conops_detected:
            self.SecondConOpsPeriod = None        

    def extractConOpsPeriodEphemFile(self, ephemeris_file_rows, maneuver_spec_file_rows, target_spec_file_rows, conops_period, ephem_file_index):
        """Locates ConOps period in an EMTG .ephemeris file and extracts its data, returns the row index where the search left off"""

        # we have five different data rows that we need to identify
        ephem_row_begin_thrust = None
        ephem_row_end_thrust = None 
        ephem_row_end_coast = None
        man_spec_row_begin_thrust = None
        tar_spec_row_end_thrust = None

        # The first ConOps period is simply the first time we detect thrusting in an EMTG .ephemeris file
        # first, find the row in the .ephemeris file that corresponds to the first thrusting period
        ephem_row_index = ephem_file_index
        thrusting_detected = False
        for index in range(ephem_file_index, len(ephemeris_file_rows)):
            u_norm = np.sqrt(ephemeris_file_rows[index].control_x**2 + ephemeris_file_rows[index].control_y**2 + ephemeris_file_rows[index].control_z**2)
            if u_norm > 0.0:
                thrusting_detected = True
                # we have our first thrusting event
                ephem_row_begin_thrust = ephemeris_file_rows[index]
                # extract the spacecraft state prior to thrusting
                conops_period.StateBeforeThrusting = np.array([ephem_row_begin_thrust.spacecraft_position_x, \
                                                               ephem_row_begin_thrust.spacecraft_position_y, \
                                                               ephem_row_begin_thrust.spacecraft_position_z, \
                                                               ephem_row_begin_thrust.spacecraft_velocity_x, \
                                                               ephem_row_begin_thrust.spacecraft_velocity_y, \
                                                               ephem_row_begin_thrust.spacecraft_velocity_z, \
                                                               ephem_row_begin_thrust.spacecraft_mass, \
                                                               ephem_row_begin_thrust.julian_date])
                break
            ephem_row_index = ephem_row_index + 1

        if not thrusting_detected:
            print("No thrusting was detected in this .ephemeris file.")        
            conops_period = None
            return

        # second, find the row in the maneuver spec. file that corresponds to the beginning of the first thrusting period
        found_maneuver_spec_begin_thrust = False
        for spec_row in maneuver_spec_file_rows:
            if np.abs(spec_row.maneuvers[0].julian_date - ephem_row_begin_thrust.julian_date) < 0.01:
                found_maneuver_spec_begin_thrust = True
                # we found the corresponding data row in the maneuver spec. file
                man_spec_row_begin_thrust = spec_row
                # for now, make the ControlVector equal to the first maneuver event's control 
                conops_period.ControlVector = np.array([man_spec_row_begin_thrust.maneuvers[0].control_x, man_spec_row_begin_thrust.maneuvers[0].control_y, man_spec_row_begin_thrust.maneuvers[0].control_z])
                # thrust time is just the total of the duration of each maneuver event from the control segment
                conops_period.ThrustTime = 0.0
                for maneuver in man_spec_row_begin_thrust.maneuvers:
                    conops_period.ThrustTime = conops_period.ThrustTime + maneuver.duration
                break
        
        if not found_maneuver_spec_begin_thrust:
            raise Exception("Could not identify a row in the maneuver spec. file that corresponds with the beginning of the thrusting period.")

        # third, find the row in the .ephemeris file that encodes the end of the thrusting period
        # take the starting thrust epoch and add the thrust duration to it, then find the row in the .ephemeris file that matches that new epoch
        thrust_end_epoch = conops_period.StateBeforeThrusting[-1] + conops_period.ThrustTime / 86400.0
        found_ephem_end_thrust = False
        for row_index in range(ephem_row_index, len(ephemeris_file_rows)):
            if np.abs(thrust_end_epoch - ephemeris_file_rows[row_index].julian_date) < 0.01:
                found_ephem_end_thrust = True
                ephem_row_end_thrust = ephemeris_file_rows[row_index]
                ephem_row_index = row_index
                conops_period.StateAfterThrusting = np.array([ephem_row_end_thrust.spacecraft_position_x, \
                                                              ephem_row_end_thrust.spacecraft_position_y, \
                                                              ephem_row_end_thrust.spacecraft_position_z, \
                                                              ephem_row_end_thrust.spacecraft_velocity_x, \
                                                              ephem_row_end_thrust.spacecraft_velocity_y, \
                                                              ephem_row_end_thrust.spacecraft_velocity_z, \
                                                              ephem_row_end_thrust.spacecraft_mass, \
                                                              ephem_row_end_thrust.julian_date])
                break

        if not found_ephem_end_thrust:
            raise Exception("Could not identify a row in the .ephemeris file that corresponds with the end of the thrusting period.")

        # finally, we need to find out how long the nav-coast is and what our state is at the end of that coast
        # the target spec file encodes when the nav-coast ends
        # loop through the target spec. file rows until the event name matches the event name in the maneuver spec. file
        found_maneuver_target_match = False
        for row in target_spec_file_data:
            if row.name == man_spec_row_begin_thrust.name:
                found_maneuver_target_match = True
                tar_spec_row_end_thrust = row
                break

        if not found_maneuver_target_match:
            raise Exception("Could not identify a row in the .mission_target_spec file that corresponds with a row in the .mission_maneuver_spec file.")
        
        # Compute the length of the nav-coast, and determine the epoch when it ends
        conops_period.CoastTime = (tar_spec_row_end_thrust.julian_date - conops_period.StateAfterThrusting[-1]) * 86400.0
        nav_coast_end_epoch = conops_period.StateAfterThrusting[-1] + conops_period.CoastTime / 86400.0
        
        # continue looping through the ephemeris file until we find the row that corresponds with the end of the nav-coast
        found_ephem_end_nav_coast = False
        for row_index in range(ephem_row_index, len(ephemeris_file_rows)):
            if np.abs(nav_coast_end_epoch - ephemeris_file_rows[row_index].julian_date) < 0.01:
                found_ephem_end_nav_coast = True
                ephem_row_end_coast = ephemeris_file_rows[row_index]
                ephem_row_index = row_index
                conops_period.StateAfterCoasting = np.array([ephem_row_end_coast.spacecraft_position_x, \
                                                             ephem_row_end_coast.spacecraft_position_y, \
                                                             ephem_row_end_coast.spacecraft_position_z, \
                                                             ephem_row_end_coast.spacecraft_velocity_x, \
                                                             ephem_row_end_coast.spacecraft_velocity_y, \
                                                             ephem_row_end_coast.spacecraft_velocity_z, \
                                                             ephem_row_end_coast.spacecraft_mass, \
                                                             ephem_row_end_coast.julian_date])
                break

        if not found_ephem_end_nav_coast:
            raise Exception("Could not identify a row in the .ephemeris file that corresponds with the end of the nav-coast period.")
        return ephem_row_index



class EphemerisFileRow(object):
    """Contains the data from one row of an EMTG .ephemeris file"""        
    def __init__(self):
        self.gregorian_date = ""
        self.julian_date = 0.0
        self.spacecraft_position_x = 0.0
        self.spacecraft_position_y = 0.0
        self.spacecraft_position_z = 0.0
        self.spacecraft_velocity_x = 0.0
        self.spacecraft_velocity_y = 0.0
        self.spacecraft_velocity_z = 0.0
        self.spacecraft_mass = 0.0
        self.control_x = 0.0
        self.control_y = 0.0
        self.control_z = 0.0
        self.thrust_magnitude = 0.0
        self.mass_flow_rate = 0.0
        self.specific_impulse = 0.0
        self.number_of_active_thrusters = 0
        self.active_power = 0.0
        self.throttle_level = "none"


class EphemerisFileReader(object):
    def __init__(self):
        self.ephemeris_file = ""
        self.ephemeris_file_rows = []

    def generateCleanEphemerisFileForBSP(self,threshold = 120.0):
        import os
        # comb through each ephemeris file line
        # if the difference in Julian Date between the current line and the previous one is fewer than 60 seconds, don't add the line
        clean_ephem_file_rows = []
        clean_ephem_file_rows.append(self.ephemeris_file_rows[0])
        for index in range(1, len(self.ephemeris_file_rows)-1):
            if (self.ephemeris_file_rows[index].julian_date - clean_ephem_file_rows[-1].julian_date)*86400.0 >= threshold and (self.ephemeris_file_rows[-1].julian_date - self.ephemeris_file_rows[index].julian_date)*86400.0 >= threshold:
                clean_ephem_file_rows.append(self.ephemeris_file_rows[index])
            else:
                print("Row " + str(index) + " removed from .ephemeris file \n")
        clean_ephem_file_rows.append(self.ephemeris_file_rows[-1])
        
        
        # we dont' want to overfit the .bsp, so go through the data again, and if any of the rows are too close together, don't include them
        #if (clean_ephem_file_rows[1].julian_date - clean_ephem_file_rows[0].julian_date < 60.0):
        #    clean_ephem_file_rows = clean_ephem_file_rows[0::2]
        

        # now write out a new, clean .ephemeris file suitable for .bsp conversion by SPICE
        splitty = os.path.split(os.path.abspath(self.ephemeris_file))
        with open(splitty[0]+'/'+splitty[1].split('.')[0]+'_clean'+'.ephemeris', 'w') as clean_ephem_file:
            for row in clean_ephem_file_rows:
                clean_ephem_file.write(row.gregorian_date.split('TDB')[0] + ',' + "{:.15E}".format(row.spacecraft_position_x) + ','
                                               + "{:.15E}".format(row.spacecraft_position_y) + ','
                                               + "{:.15E}".format(row.spacecraft_position_z) + ','
                                               + "{:.15E}".format(row.spacecraft_velocity_x) + ','
                                               + "{:.15E}".format(row.spacecraft_velocity_y) + ','
                                               + "{:.15E}".format(row.spacecraft_velocity_z) + '\n')


            
        
    def parseEMTGephemerisFile(self, ephemeris_file):
        """Parses an EMTG .ephemeris file and returns a List whose entries are EphemerisFileRow objects"""        
        import csv
        import re      

        # parse the EMTG .ephemeris file
        self.ephemeris_file = ephemeris_file
        ephemeris_rows = []
        with open(self.ephemeris_file) as ephem_file:
            csv_reader = csv.reader(ephem_file, delimiter = ',')
            headers = next(csv_reader)
            headers[0] = "gregorian_date"
            headers[:] = [re.sub(r'\(.*\)', '', header) for header in headers]
            if len(headers) > 17:
                raise Exception("EMTG .ephemeris file format has changed, please modify the EphemerisFileReader accordingly!!")

            for row in csv_reader:
                ephem_row = EphemerisFileRow()
                ephem_row.gregorian_date = row[0]

                # if no time system is specified, assume it is TDB
                if "TDB" not in ephem_row.gregorian_date or "TDT" not in ephem_row.gregorian_date or "UTC" not in ephem_row.gregorian_date:
                    ephem_row.gregorian_date = ephem_row.gregorian_date + " TDB"

                ephem_row.julian_date = spice.str2et(ephem_row.gregorian_date) / 86400.0 + 2451545.0
                ephem_row.spacecraft_position_x = float(row[1])
                ephem_row.spacecraft_position_y = float(row[2])
                ephem_row.spacecraft_position_z = float(row[3])
                ephem_row.spacecraft_velocity_x = float(row[4])
                ephem_row.spacecraft_velocity_y = float(row[5])
                ephem_row.spacecraft_velocity_z = float(row[6])
                try:
                    ephem_row.spacecraft_mass = float(row[7])
                    ephem_row.control_x = float(row[8])
                    ephem_row.control_y = float(row[9])
                    ephem_row.control_z = float(row[10])
                    ephem_row.thrust_magnitude = float(row[11])
                    ephem_row.mass_flow_rate = float(row[12])
                    ephem_row.specific_impulse = float(row[13])
                    ephem_row.number_of_active_thrusters = int(row[14])
                    ephem_row.active_power = float(row[15])
                    ephem_row.throttle_level = str(row[16])
                except:
                    pass
                ephemeris_rows.append(ephem_row)
        
        self.ephemeris_file_rows = ephemeris_rows
        return ephemeris_rows

    
class ManeuverData(object):
    """Contains the data for a single maneuver in a single row of an EMTG .mission_maneuver_spec file"""
    def __init__(self):
        self.frame = ""
        self.gregorian_date = ""
        self.julian_date = 0.0
        self.control_x = 0.0
        self.control_y = 0.0
        self.control_z = 0.0
        self.control_magnitude = 0.0
        self.starting_spacecraft_mass = 0.0
        self.mass_flow_rate = 0.0
        self.duty_cycle = 0.0
        self.final_spacecraft_mass = 0.0
        self.delta_v = 0.0
        self.duration = 0.0
        self.epochET_seconds = 0.0
        
class ControlSegment(object):
    """Contains the data from one row of an EMTG .mission_maneuver_spec file"""
    def __init__(self, num_maneuvers_in_segment):
        self.name = ""
        self.num_maneuvers = num_maneuvers_in_segment
        self.maneuvers = [ManeuverData() for _ in range(self.num_maneuvers)]

    def populateManeuverData(self, control_segment_data):
        """Parses the data from a single row in a .mission_meneuver_spec file into multiple maneuver event objects"""
        data_index = 0
        for maneuver in self.maneuvers:
            maneuver.frame = control_segment_data[data_index].strip()
            maneuver.gregorian_date = control_segment_data[data_index + 1].strip()
                
            # if no time system is specified, assume it is TDB
            if "TDB" not in maneuver.gregorian_date or "TDT" not in maneuver.gregorian_date or "UTC" not in maneuver.gregorian_date:
                maneuver.gregorian_date = maneuver.gregorian_date + " TDB"

            maneuver.julian_date = spice.str2et(maneuver.gregorian_date) / 86400.0 + 2451545.0
            maneuver.GMAT_julian_date = SpiceyUtil.greg2GMATJulian(maneuver.gregorian_date)
            maneuver.control_x = float(control_segment_data[data_index + 2])
            maneuver.control_y = float(control_segment_data[data_index + 3])
            maneuver.control_z = float(control_segment_data[data_index + 4])
            maneuver.control_magnitude = float(control_segment_data[data_index + 5])
            maneuver.starting_spacecraft_mass = float(control_segment_data[data_index + 6])
            maneuver.mass_flow_rate = float(control_segment_data[data_index + 7])
            maneuver.duty_cycle = float(control_segment_data[data_index + 8])
            maneuver.final_spacecraft_mass = float(control_segment_data[data_index + 9])
            maneuver.delta_v = float(control_segment_data[data_index + 10])
            maneuver.duration = float(control_segment_data[data_index + 11])
            maneuver.epochET_seconds = float(control_segment_data[data_index + 12])
            data_index = data_index + 13

class ManeuverSpecFileReader(object):

    def __init__(self):
            self.maneuver_spec_file = ""

    def parseEMTGmaneuverSpecFile(self, man_spec_file):
        """Parses an EMTG .mission_maneuver_spec file"""
        import csv
        import re

        self.maneuver_spec_file = man_spec_file
        
        control_segment_data = []

        with open(self.maneuver_spec_file) as man_spec_file:
            csv_reader = csv.reader(man_spec_file, delimiter = ',')
            headers = next(csv_reader)
            headers[:] = [re.sub(r'\(.*\)', '', header) for header in headers]
            headers[:] = [re.sub(r'\[.*\]', '', header) for header in headers]
            if len(headers) != 16:
                raise Exception("EMTG .mission_maneuver_spec file format has changed, please modify the ManeuverSpecFileReader accordingly!!")
            headers = headers[:-1]
            for row in csv_reader:
                if row[0][0] == '#':
                    continue
                num_maneuvers_in_segment = int(row[1])
                control_segment_data.append(ControlSegment(num_maneuvers_in_segment))
                control_segment_data[-1].name = row[0]
                control_segment_data[-1].populateManeuverData(row[2:])

        return control_segment_data


class ManeuverTargetData(object):
    """Contains the data for maneuver target (i.e. a single row of an EMTG .mission_target_spec file)."""
    def __init__(self):
        self.name = ""
        self.central_body = ""
        self.gregorian_date = ""
        self.julian_date = 0.0
        self.frame = ""
        self.target_spacecraft_position_x = 0.0
        self.target_spacecraft_position_y = 0.0
        self.target_spacecraft_position_z = 0.0
        self.target_spacecraft_velocity_x = 0.0
        self.target_spacecraft_velocity_y = 0.0
        self.target_spacecraft_velocity_z = 0.0
        self.target_spacecraft_mass = 0.0
        self.target_Bplane_R = 0.0
        self.target_Bplane_T = 0.0
        self.epochET_seconds = 0.0

    def populateTargetData(self, target_data):
        self.name = target_data[0].lstrip(' ')
        self.sample_number = target_data[1]
        self.central_body = target_data[2].lstrip(' ')
        self.frame = target_data[3].lstrip(' ')
        self.gregorian_date = target_data[4].lstrip(' ')

        # if no time system is specified, assume it is TDB
        if "TDB" not in self.gregorian_date or "TDT" not in self.gregorian_date or "UTC" not in self.gregorian_date:
            self.gregorian_date = self.gregorian_date + " TDB"

        self.julian_date = SpiceyUtil.greg2Julian(self.gregorian_date)
        self.GMAT_julian_date = SpiceyUtil.greg2GMATJulian(self.gregorian_date)
        self.target_spacecraft_position_x = float(target_data[5])
        self.target_spacecraft_position_y = float(target_data[6])
        self.target_spacecraft_position_z = float(target_data[7])
        self.target_spacecraft_velocity_x = float(target_data[8])
        self.target_spacecraft_velocity_y = float(target_data[9])
        self.target_spacecraft_velocity_z = float(target_data[10])
        self.target_spacecraft_mass = float(target_data[11])
        self.target_Bplane_R = float(target_data[12])
        self.target_Bplane_T = float(target_data[13])
        self.epochET_seconds = float(target_data[14])

class TargetSpecFileReader(object):
    
    def __init__(self):
            self.target_spec_file = ""

    def parseEMTGtargetSpecFile(self, target_spec_file):
        """Parses an EMTG .mission_maneuver_spec file"""
        import csv
        import re

        self.target_spec_file = target_spec_file
        
        target_event_data = []

        with open(self.target_spec_file) as tar_spec_file:
            csv_reader = csv.reader(tar_spec_file, delimiter = ',')
            headers = next(csv_reader)
            headers[:] = [re.sub(r'\(.*\)', '', header) for header in headers]
            headers[:] = [re.sub(r'\[.*\]', '', header) for header in headers]
            if len(headers) != 15:
                raise Exception("EMTG .mission_target_spec file format has changed, please modify the TargetSpecFileReader accordingly!!")
            for row in csv_reader:
                if row[0][0] == '#':
                    continue
                target_event_data.append(ManeuverTargetData())
                target_event_data[-1].populateTargetData(row)

        return target_event_data

if __name__ == '__main__':
    #you need to modify this to make it do anything
    my_mission_options = MO.MissionOptions("stuff.emtgopt")
    SPICE_ephem_directory = "C:/emtg/testatron/universe/ephemeris_files/"
    spice_handler = SpiceyUtil.SpiceHandler(SPICE_ephem_directory)
    spice_handler.loadSpiceFiles()
    
    dir = "somewhere"
    mission_name = "myMission"   
    
    my_maneuver_spec_reader = ManeuverSpecFileReader()
    maneuver_spec_file_data = my_maneuver_spec_reader.parseEMTGmaneuverSpecFile(dir + mission_name + ".mission_maneuver_spec")        
    
    my_target_spec_reader = TargetSpecFileReader()
    target_spec_file_data = my_target_spec_reader.parseEMTGtargetSpecFile(dir + mission_name + ".mission_target_spec")        
    
    
    for index in range(0, len(maneuver_spec_file_data)):
        segment = maneuver_spec_file_data[index]
        print("CONTROL SEGMENT " + str(index))
        for man_num in range(0, len(segment.maneuvers)):
            print ("Maneuver " + str(man_num) 
                 + " Greg. date: " + segment.maneuvers[man_num].gregorian_date
                 + " duration: " + str(segment.maneuvers[man_num].duration) 
                 + ", mass flow: " + str(segment.maneuvers[man_num].mass_flow_rate)
                 + ", thrust: " + str(segment.maneuvers[man_num].control_magnitude)
                 + ", ux: " + str(segment.maneuvers[man_num].control_x)
                 + ", uy: " + str(segment.maneuvers[man_num].control_y)
                 + ", uz: " + str(segment.maneuvers[man_num].control_z)
            )


    spice_handler.unloadSpiceFiles()