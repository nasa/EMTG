import spiceypy
try:
    import spiceypy as spice
    import spiceypy.utils.support_types as stypes
except:
    print("spiceypy not available")

class SpiceHandler(object):
    def __init__(self, ephem_directory_in, only_tls=False):
        self.SPICE_ephem_directory = ephem_directory_in
        self.loadedSpiceFiles = []
        self.only_tls = only_tls

    def loadSpiceFiles(self):
        from os import listdir
        from os.path import isfile, join

        # load SPICE ephemeris files
        spiceBinaryFiles = [f for f in listdir(self.SPICE_ephem_directory) if isfile(join(self.SPICE_ephem_directory, f))]
        for spiceFile in spiceBinaryFiles:
            try:
                if self.only_tls:
                    if spiceFile.endswith('.tls'):
                        spice.furnsh(self.SPICE_ephem_directory + spiceFile)
                        self.loadedSpiceFiles.append(spiceFile)
                        print("loaded " + spiceFile)
                else:
                    if spiceFile.endswith('.bsp') or spiceFile.endswith('.tls') or spiceFile.endswith('.tpc'):
                        spice.furnsh(self.SPICE_ephem_directory + spiceFile)
                        self.loadedSpiceFiles.append(spiceFile)
                        print("loaded " + spiceFile)
            except:
                print("Could not furnish ", spiceFile)

    def unloadSpiceFiles(self):
        # unload spice files
        for spiceFile in self.loadedSpiceFiles:
            spice.unload(spiceFile)
        spice.kclear()

def greg2Greg(gregorian_date, time_system_string):
    return spice.timout(spice.str2et(gregorian_date), "YYYY Mon DD ::" + time_system_string + " HR:MN:SC.######")

def gregDateOffsetCalculator(gregorian_date, delta_days, time_system_string):

    # returns the time_system_string Gregorian date, offset from the input date string by delta_days
    # gregorian_date: SPICE compatible Gregorian date string
    # delta_days: number of days that you want to advance the Gregorian date string by
    # time_system_string: string representation of the output requested for the gregorian_date's time system (e.g. "TDB", "UTC" etc.)

    return spice.timout(spice.str2et(str(spice.str2et(gregorian_date) / 86400.0 + 2451545.0 + delta_days) + " JD TDB"), "YYYY Mon DD ::" + time_system_string + " HR:MN:SC.######")

def greg2Julian(gregorian_date):

    # returns the TDB Julian date corresponding to the input Gregorian date
    # gregorian_date: SPICE compatible Gregorian date string

    return spice.str2et(gregorian_date) / 86400.0 + 2451545.0

def greg2GMATJulian(gregorian_date):

    # returns the TDB Julian date corresponding to the input Gregorian date
    # gregorian_date: SPICE compatible Gregorian date string

    return greg2Julian(gregorian_date) - 2430000.0

def julian2Greg(julian_date, timestring='TDB'):
    # returns the TDB Gregorian date
    # julian_date: TDB julian epoch
    # delta_days: number of days that you want to advance the Gregorian date string by
    # time_system_string: string representation of the input gregorian_date's time system (e.g. "TDB", "UTC" etc.)

    return spice.timout((julian_date - 2451545.0) * 86400.0, "YYYY Mon DD ::" + timestring + " HR:MN:SC.###### " + timestring)

def secSinceJ20002Greg(secSinceJ2000):
    # returns the TDB Gregorian date
    # secSinceJ2000: TDB julian seconds since J2000
    # delta_days: number of days that you want to advance the Gregorian date string by
    # time_system_string: string representation of the input gregorian_date's time system (e.g. "TDB", "UTC" etc.)

    return spice.timout(secSinceJ2000, "YYYY Mon DD ::TDB HR:MN:SC.###### TDB")

class OccultationEvent(object):
    def __init__(self):
        self.start = 0
        self.stop = 0
        self.start_JD = 0.0
        self.stop_JD = 0.0
        self.duration = 0.0
        self.type = "FULL"
        self.observer = ''
        self.occulting_body = ''
        self.occulting_body_shape = ''
        self.target = ''
        self.target_body_shape = ''

class OccultingBody(object):
    def __init__(self):
        self.name = '301'
        self.shape = 'ELLIPSOID'
        self.frame = 'IAU_MOON'
        self.occultation_types = ['FULL']
        self.search_start = ''
        self.search_end = ''
        self.search_start_ET_seconds = 0.0
        self.search_end_ET_seconds = 0.0  
        
class OccultationDetector(object):
    def __init__(self, ephem_directory_in):
        #load SPICE  
        try:
            import spiceypy as spice
            import spiceypy.utils.support_types as stypes
        except:
            print("spiceypy not available")      
        self.mySPICEdriver = SpiceHandler(ephem_directory_in)
        self.mySPICEdriver.loadSpiceFiles()
        self.TDBFMT = 'DD-MON-YYYY HR:MN:SC.###### TDB ::TDB'
        self.start = 0
        self.stop = 0
        self.search_step_size = 120.0
        self.SPICE_search_window = stypes.SPICEDOUBLE_CELL(1000000)
        self.results_window = stypes.SPICEDOUBLE_CELL(1000000)
        self.observer = 'EARTH'
        self.target = 'SUN'
        self.target_body_shape = 'ELLIPSOID'
        self.target_frame = 'IAU_SUN'
        self.occulting_bodies = []
        self.aberration_correction = 'NONE'
        self.cumulative_results = {}

    def __del__(self):
        self.mySPICEdriver.unloadSpiceFiles()

    def computeOccultations(self, observer_in, occultation_bodies_in, target_in, target_shape_in, target_frame_in, step_size_in, aberration_correction_in, greg_format_string_in):

        self.greg_format_string = greg_format_string_in
        split_string = self.greg_format_string.split(' ')
        self.time_system_string = [i for i in split_string if '::' in i][0][2:]

        self.observer = str(observer_in)
        self.occultation_bodies = occultation_bodies_in
        self.target = target_in
        self.target_shape = target_shape_in
        self.target_frame = target_frame_in
        self.search_step_size = step_size_in
        self.aberration_correction = aberration_correction_in
        self.cumulative_results = {}
        for body in self.occultation_bodies:
            # add a new key to the results dictionary for this body
            self.cumulative_results[body.name] = []
            for occultation_type in body.occultation_types:
                body.search_start_ET_seconds = spice.str2et( body.search_start )
                body.search_end_ET_seconds = spice.str2et( body.search_end )
                spice.wninsd( body.search_start_ET_seconds, body.search_end_ET_seconds, self.SPICE_search_window )
                spice.gfoclt( occultation_type,
                              body.name, 
                              body.shape, 
                              body.frame,
                              self.target, 
                              self.target_shape, 
                              self.target_frame,   
                              self.aberration_correction,
                              self.observer, 
                              self.search_step_size, 
                              self.SPICE_search_window, 
                              self.results_window )
                winsiz = spice.wncard( self.results_window )
                for body_index in range(winsiz):
                    [intbeg, intend] = spice.wnfetd( self.results_window, body_index )
                    btmstr = spice.timout( intbeg, self.greg_format_string )
                    etmstr = spice.timout( intend, self.greg_format_string )
                    occultation_event = OccultationEvent()
                    occultation_event.start = btmstr
                    occultation_event.stop = etmstr
                    occultation_event.start_JD = greg2Julian(btmstr)
                    occultation_event.stop_JD = greg2Julian(etmstr)
                    occultation_event.duration = intend - intbeg
                    occultation_event.type = occultation_type
                    occultation_event.observer = self.observer
                    occultation_event.occulting_body = body.name
                    occultation_event.occulting_body_shape = body.shape
                    occultation_event.target = self.target
                    occultation_event.target_body_shape = self.target_shape
                    self.cumulative_results[body.name].append(occultation_event)
                

        
        #print( 'Umbra: {:s} : {:s}'.format( btmstr, etmstr ) + '  duration:' + str(intend-intbeg) + ' seconds' )

    def generateOccultationReport(self, report_file_path):
        
        import csv

        with open(report_file_path, mode='w', newline='') as report_file:

            report_writer = csv.writer(report_file, delimiter=',')
            report_writer.writerow(["Occultation Type", "Gregorian start ["+ self.time_system_string + "]", "Gregorian stop ["+ self.time_system_string + "]", "Julian start [TDB]", "Julian stop [TDB]", "Duration [s]"])           
            for name,events in self.cumulative_results.items(): 

                # first sort the eclipse events for this body
                events.sort(key=lambda x: x.start_JD, reverse=False)
                
                report_writer.writerow([name + ' occultation of ' + self.target])

                if len(events) == 0:
                    print("No occultation events found for " + name + ".")
                else:
                    print( str(len(events)) + " occultation events found for " + name + ".")
                    for event in events:
                        if event.type == 'FULL' or event.type == 'Full' or event.type == 'full':
                            if self.target == 'SUN' or self.target == 'sun' or self.target == '10':
                                type = 'Umbra'
                            else:
                                type = 'Full     '
                            report_writer.writerow([type, '{:s}'.format(event.start), '{:s}'.format(event.stop ), str(event.start_JD), str(event.stop_JD), str(event.duration)])
                        elif event.type == 'PARTIAL' or event.type == 'Partial' or event.type == 'partial':
                            if self.target == 'SUN' or self.target == 'sun' or self.target == '10':
                                type = 'Penumbra'
                            else:
                                type = 'Partial'
                            report_writer.writerow([type, '{:s}'.format(event.start), '{:s}'.format(event.stop ), str(event.start_JD), str(event.stop_JD), str(event.duration)])
                        elif event.type == 'ANNULAR' or event.type == 'Annular' or event.type == 'annular':
                            report_writer.writerow([type, '{:s}'.format(event.start), '{:s}'.format(event.stop ), str(event.start_JD), str(event.stop_JD), str(event.duration)])
                report_writer.writerow('\n \n')

class OccultationData(object):
    def __init__(self, source_body, occulting_body, target):
        self.source_body_name = source_body
        self.occulting_body_name = occulting_body
        self.target_name = target
        self.data = []

class OccultationBody(object):
    def __init__(self, name, radius):
        self.name = name
        self.radius = radius

class TrajectoryDataFrame(object):
    def __init__(self):
        self.start_epoch = "01-JAN-2000 00:00:00.000 TDB"
        self.end_epoch = "02-JAN-2000 00:00:00.000 TDB"
        self.greg_time_system_string = 'TDB'
        self.julian_time_system_string = 'TDB'
        self.body_pairs = []
        self.frame = 'J2000' 
        self.time_step = 86400.0
        self.data = {}
        self.data['SEP'] = []
        self.data['SPE'] = []
        self.body_dict = {}
        self.occultation_dicts = None
        self.data['occultations'] = []
        self.spacecraft_SPICE_ID = -123

    def setBodyPairs(self, body_pairs):
        self.body_pairs = body_pairs
        for pair in self.body_pairs:
            if pair[0] in self.data.keys():
                self.data[pair[0]][pair[1]] = [] 
            else:
                self.data[pair[0]] = {pair[1] : []}

    def setOccultations(self, occultation_dicts):
        self.occultation_dicts = occultation_dicts
        for dic in self.occultation_dicts:
            self.data['occultations'].append(OccultationData(dic['source body'], dic['occulting body'], dic['target']))

    def setOccultationBodyInformation(self, body_dict):
        self.body_dict[body_dict['name']] = OccultationBody(body_dict['name'], body_dict['radius'])


    def generateReport(self, report_file_path):

        import csv
        with open(report_file_path, mode='w', newline='') as report_file:
            report_writer = csv.writer(report_file, delimiter=',')

            header_array = ["Gregorian Date ["+ self.greg_time_system_string + "]", "Julian Date ["+ self.julian_time_system_string + "]"]
            # for each body_pair
            for pair in self.body_pairs:
                header_array.append(pair[0] + "->" + pair[1] + ' [' + self.frame + ']'+ ' X [km]')
                header_array.append(pair[0] + "->" + pair[1] + ' [' + self.frame + ']'+ ' Y [km]')
                header_array.append(pair[0] + "->" + pair[1] + ' [' + self.frame + ']'+ ' Z [km]')
                header_array.append(pair[0] + "->" + pair[1] + ' [' + self.frame + ']'+ ' VX [km/s]')
                header_array.append(pair[0] + "->" + pair[1] + ' [' + self.frame + ']'+ ' VY [km/s]')
                header_array.append(pair[0] + "->" + pair[1] + ' [' + self.frame + ']'+ ' VZ [km/s]')
                header_array.append(pair[0] + "->" + pair[1] + ' distance [km]')
                header_array.append(pair[0] + "->" + pair[1] + ' distance [AU]')
            
            # angle header entries
            header_array.append("Sun-Earth-Spacecraft Angle [deg]")
            header_array.append("Sun-Spacecraft-Earth Angle [deg]")

            # eclipse header entries
            for occlt in self.data['occultations']:
                header_array.append(occlt.occulting_body_name + " occultation of " + occlt.source_body_name + ' [%]')
                header_array.append("Occultation type")

            # write the header row to file
            report_writer.writerow(header_array)

            # for each data epoch
            for index in range(0, len(self.data['SEP'])):
                row_array = ['{:s}'.format(self.data['SEP'][index][0]), '{:s}'.format(self.data['SEP'][index][1])]                
                for pair in self.body_pairs:
                    data_entry = self.data[pair[0]][pair[1]][index]
                    epoch_greg = data_entry[0]
                    epoch_JD = data_entry[1]
                    state = data_entry[2]
                    distance_km = data_entry[3]
                    distance_AU = data_entry[4]
                    row_array.extend(['{:.6f}'.format(state[0]),
                                      '{:.6f}'.format(state[1]),
                                      '{:.6f}'.format(state[2]),
                                      '{:.6f}'.format(state[3]),
                                      '{:.6f}'.format(state[4]),
                                      '{:.6f}'.format(state[5]),
                                      '{:.6f}'.format(distance_km),
                                      '{:.9f}'.format(distance_AU)
                                     ])
                row_array.append(self.data['SEP'][index][2])
                row_array.append(self.data['SPE'][index][2])

                for occlt in self.data['occultations']:
                    row_array.append(occlt.data[index][2])
                    row_array.append(occlt.data[index][3])
                
                report_writer.writerow(row_array)
                


class GeometryComputer(object):
    import numpy as np
    def __init__(self, ephem_directory_in):
        
        #load SPICE  
        try:
            import spiceypy as spice
            import spiceypy.utils.support_types as stypes
        except:
            print("spiceypy not available")      
        self.mySPICEdriver = SpiceHandler(ephem_directory_in)
        self.mySPICEdriver.loadSpiceFiles()
        self.AU = 149597870.7 # km
   
        
    def __del__(self):
        self.mySPICEdriver.unloadSpiceFiles()

    def setSPICEnameIDpair(self, name, ID):
        spice.boddef(name, ID)

    def computeDataFrameGeometries(self, data_frames, greg_format_string, julian_format_string):

        for data_frame in data_frames:
            self.greg_format_string = greg_format_string
            self.julian_format_string = julian_format_string
            greg_split_string = self.greg_format_string.split(' ')
            data_frame.greg_time_system_string = [i for i in greg_split_string if '::' in i][0][2:]
            julian_split_string = self.julian_format_string.split(' ')
            data_frame.julian_time_system_string = [i for i in julian_split_string if '::' in i][0][2:]
            current_epoch = data_frame.start_epoch
            current_epoch = spice.str2et(current_epoch)
            stop_epoch = spice.str2et(data_frame.stop_epoch)
            time_remaining = stop_epoch - current_epoch

            while time_remaining >= 0.0:
                greg_date_string = spice.timout(current_epoch, self.greg_format_string)
                JD_date = spice.timout(current_epoch, self.julian_format_string)
                for body_pair in data_frame.body_pairs:
                    # get target body state
                    body_state, light_times = spice.spkez(spice.bodn2c(body_pair[1]), current_epoch, data_frame.frame, 'NONE', spice.bodn2c(body_pair[0]))
                    distance_km = self.np.linalg.norm(body_state[0:3])
                    distance_AU = distance_km / self.AU
                    data_frame.data[body_pair[0]][body_pair[1]].append([greg_date_string, JD_date, body_state, distance_km, distance_AU])

                # compute solar geometries
                SEP_angle, SPE_angle = self.computeSEPandSPEangles(current_epoch, data_frame.spacecraft_SPICE_ID)
                data_frame.data['SEP'].append([greg_date_string, JD_date, SEP_angle])
                data_frame.data['SPE'].append([greg_date_string, JD_date, SPE_angle])

                # compute any eclipses
                self.computeCircleCircleOccultations(current_epoch, greg_date_string, JD_date, data_frame)

                if time_remaining == 0.0:
                    break
                else:
                    if time_remaining >= data_frame.time_step:
                        current_epoch += data_frame.time_step
                        time_remaining -= data_frame.time_step    
                    else:            
                        time_remaining -= stop_epoch - current_epoch
                        current_epoch += stop_epoch - current_epoch
  
    def computeSEPandSPEangles(self, current_epoch, spacecraft_SPICE_ID):
        
        # get position of object w.r.t. Earth
        spacecraft_state_wrt_Earth, light_times = spice.spkez(spacecraft_SPICE_ID, current_epoch, "J2000", 'NONE', 399)

        # get position of Sun w.r.t. Earth
        sun_state_wrt_Earth, light_times = spice.spkez(10, current_epoch, "J2000", 'NONE', 399)

        Earth_state_wrt_spacecraft = [ -x for x in spacecraft_state_wrt_Earth]
        sun_state_wrt_spacecraft = sun_state_wrt_Earth - spacecraft_state_wrt_Earth

        cosAngle_SEP = self.np.dot(sun_state_wrt_Earth[0:3], spacecraft_state_wrt_Earth[0:3])/self.np.linalg.norm(sun_state_wrt_Earth[0:3])/self.np.linalg.norm(spacecraft_state_wrt_Earth[0:3])
        cosAngle_SPE = self.np.dot(sun_state_wrt_spacecraft[0:3], Earth_state_wrt_spacecraft[0:3])/self.np.linalg.norm(sun_state_wrt_spacecraft[0:3])/self.np.linalg.norm(Earth_state_wrt_spacecraft[0:3])

        SEP_angle = self.np.arccos(cosAngle_SEP) * 180.0/self.np.pi
        SPE_angle = self.np.arccos(cosAngle_SPE) * 180.0/self.np.pi

        return SEP_angle, SPE_angle

    def computeCircleCircleOccultations(self, current_epoch, greg_date_string, JD_date, data_frame):
        
        for occlt in data_frame.data['occultations']:
            
            # source body position w.r.t. target
            spacecraft_state_wrt_source, light_times = spice.spkez(spice.bodn2c(occlt.target_name), current_epoch, data_frame.frame, 'NONE', spice.bodn2c(occlt.source_body_name))
            source_range = self.np.linalg.norm(spacecraft_state_wrt_source[0:3])

            # occultation body positon w.r.t. target
            spacecraft_state_wrt_occulter, light_times = spice.spkez(spice.bodn2c(occlt.target_name), current_epoch, data_frame.frame, 'NONE', spice.bodn2c(occlt.occulting_body_name))
            occulter_range = self.np.linalg.norm(spacecraft_state_wrt_occulter[0:3])

            # compute the SPO (source-probe-occulter) angle
            cosAngle_SPO = self.np.dot(spacecraft_state_wrt_source[0:3], spacecraft_state_wrt_occulter[0:3])/source_range/occulter_range
            SPO_angle = self.np.arccos(cosAngle_SPO)

            # compute half-angles
            source_half_angle = self.np.arcsin(data_frame.body_dict[occlt.source_body_name].radius / source_range)
            occulter_half_angle = self.np.arcsin(data_frame.body_dict[occlt.occulting_body_name].radius / occulter_range)
        
            # body areas
            source_area = self.np.pi * source_half_angle**2
            occulter_area = self.np.pi * occulter_half_angle**2

            # compute eclipse %
            eclipse_percentage = 0.0
            eclipse_type = "none"
            # first test if we have an eclipse condition at all
            if source_half_angle + occulter_half_angle < SPO_angle:
                eclipse_percentage = 0.0
            else:
                if SPO_angle < occulter_half_angle - source_half_angle:
                    # we are at full eclipse
                    eclipse_percentage = 100.0
                    eclipse_type = "full"
                else:
                    if SPO_angle < source_half_angle - occulter_half_angle:
                        # annular eclipse
                        eclipse_percentage = 100.0 * (occulter_area / source_area)
                        eclipse_type = "annular"
                    else:
                        # partial eclipse
                        x1 = 1.0 / 2.0 * (SPO_angle**2 - source_half_angle**2 + occulter_half_angle**2) / SPO_angle
                        a1 = occulter_half_angle**2 * self.np.arccos(x1 / occulter_half_angle) - x1 * self.np.sqrt(occulter_half_angle**2 - x1**2)
                        a2 = source_half_angle**2 * self.np.arccos((SPO_angle - x1) / source_half_angle) - (SPO_angle - x1) * self.np.sqrt(source_half_angle**2 - (SPO_angle - x1)**2)
                        eclipse_percentage = (a1 + a2) / source_area
                        eclipse_percentage *= 100.0
                        eclipse_type = "partial"

            occlt.data.append([greg_date_string, JD_date, eclipse_percentage, eclipse_type])