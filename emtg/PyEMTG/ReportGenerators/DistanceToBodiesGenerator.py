import numpy as np
import ConOpsPeriod as ConOps
try:
    import spiceypy as spice
except:
    print("spiceypy not available")

class DistanceToBodiesReport(object):
    def __init__(self):
        pass

class CloseApproach(object):
    def __init__(self):
        self.target_body = ""
        self.target_body_SPICE_ID = -1
        self.Gregorian_date = ""
        self.Julian_date = 0.0
        self.spacecraft_state = np.zeros((6,1))
        self.spacecraft_altitude = 0.0
               

class DistanceReportGenerator(object):
    def __init__(self):
        self.ephemeris_file = ""
        self.report_file = ""
        self.ephemeris_file_reader = ConOps.EphemerisFileReader()
        self.close_approaches = []

    def generateDistanceReport(self, output_file, ephemeris_file, body_SPICE_IDs = []):
        header_row = "Gregorian_date_TDB, Julian_date_TDB, "
        for bodyID in body_SPICE_IDs:
            pass            

        ephem_file_rows = self.ephemeris_file_reader.parseEMTGephemerisFile(ephemeris_file)
        for ephem_row in ephem_file_rows:
            # compute distance from the Sun
            spacecraft_distance_from_sun = np.sqrt(ephem_row.spacecraft_position_x**2 + ephem_row.spacecraft_position_y**2 + ephem_row.spacecraft_position_z**2)

            for bodyID in body_SPICE_IDs:
                # get the distance from the Sun to the body
                current_seconds_past_J2000 = spice.str2et(ephem_row.gregorian_date)
                body_state_ICRF, sun_body_light_times = spice.spkez(10, current_seconds_past_J2000, 'J2000', 'NONE', bodyID)
                spacecraft_distance_from_body = np.sqrt((ephem_row.spacecraft_position_x - body_state_ICRF[0])**2 
                                                      + (ephem_row.spacecraft_position_y - body_state_ICRF[1])**2
                                                      + (ephem_row.spacecraft_position_z - body_state_ICRF[2])**2)


    def detectCloseApproaches(self, output_file, spacecraft_SPICE_ID, body_SPICE_ID, body_radius, epoch_windows):

        for window in epoch_windows:
            # create a close approach object
            self.close_approaches.append(CloseApproach())
            self.close_approaches[-1].target_body_SPICE_ID = body_SPICE_ID
            self.close_approaches[-1].target_body = spice.bodc2n(body_SPICE_ID)

            left_epoch = spice.str2et(window[0])
            right_epoch = spice.str2et(window[1])
            middle_epoch = (right_epoch - left_epoch) / 2.0 + left_epoch
            left_state, left_light_times = spice.spkez(spacecraft_SPICE_ID, left_epoch, 'J2000', 'NONE', body_SPICE_ID)
            right_state, right_light_times = spice.spkez(spacecraft_SPICE_ID, right_epoch, 'J2000', 'NONE', body_SPICE_ID)
            middle_state, middle_light_times = spice.spkez(spacecraft_SPICE_ID, middle_epoch, 'J2000', 'NONE', body_SPICE_ID)
            left_distance = np.sqrt(left_state[0]**2 + left_state[1]**2 + left_state[2]**2)
            right_distance = np.sqrt(right_state[0]**2 + right_state[1]**2 + right_state[2]**2)
            middle_distance = np.sqrt(middle_state[0]**2 + middle_state[1]**2 + middle_state[2]**2)

            # handle the case where left epoch is post periapse and we must be on an ellipse
            if middle_distance > left_distance and middle_distance > right_distance:
                # either left or right must be the closest approach
                if left_distance > right_distance:
                    # right is closest approach
                    self.close_approaches[-1].spacecraft_state = right_state
                else:
                    # left is closest approach
                    self.close_approaches[-1].spacecraft_state = left_state
                self.close_approaches[-1].spacecraft_altitude = left_distance - body_radius
                break

            # middle must be closer than left and right
            # we need to determine if our middle is before or after periapse
            # look a little bit into the future to see if we are getting closer to the body or farther away
            peek_epoch = middle_epoch + 1.0
            peek_state, middle_light_times = spice.spkez(spacecraft_SPICE_ID, peek_epoch, 'J2000', 'NONE', body_SPICE_ID)
            peek_distance = np.sqrt(peek_state[0]**2 + peek_state[1]**2 + peek_state[2]**2)
            if peek_distance > middle_distance:
                # we have passed periapse, therefore we need to bracket with left_epoch
                right_epoch = middle_epoch
            else:
                left_epoch = middle_epoch

            iterations = 0
            tol = 0.001
            while iterations < 100:
                middle_epoch = (right_epoch - left_epoch) / 2.0 + left_epoch
                
                left_state, left_light_times = spice.spkez(spacecraft_SPICE_ID, left_epoch, 'J2000', 'NONE', body_SPICE_ID)
                right_state, right_light_times = spice.spkez(spacecraft_SPICE_ID, right_epoch, 'J2000', 'NONE', body_SPICE_ID)
                middle_state, middle_light_times = spice.spkez(spacecraft_SPICE_ID, middle_epoch, 'J2000', 'NONE', body_SPICE_ID)
                left_distance = np.sqrt(left_state[0]**2 + left_state[1]**2 + left_state[2]**2)
                right_distance = np.sqrt(right_state[0]**2 + right_state[1]**2 + right_state[2]**2)
                middle_distance = np.sqrt(middle_state[0]**2 + middle_state[1]**2 + middle_state[2]**2)


                if np.abs(left_distance - right_distance) < tol:
                    closest_approach_distance = middle_distance
                    break

                # determine where the left, middle and right points are located w.r.t. periapse
                time_step = 0.001
                middle_peek_state, middle_light_times = spice.spkez(spacecraft_SPICE_ID, middle_epoch + time_step, 'J2000', 'NONE', body_SPICE_ID)
                middle_peek_distance = np.sqrt(middle_peek_state[0]**2 + middle_peek_state[1]**2 + middle_peek_state[2]**2)
                
                if middle_peek_distance - middle_distance < 0.0:
                    left_epoch = middle_epoch
                elif middle_peek_distance - middle_distance > 0.0:
                    right_epoch = middle_epoch

            self.close_approaches[-1].Julian_date = middle_epoch / 86400.0 + 2451545.0
            self.close_approaches[-1].spacecraft_state = middle_state
            self.close_approaches[-1].spacecraft_altitude = closest_approach_distance - body_radius


if __name__ == '__main__':
    CAESAR_open_ephem_file = "C:/Users/delliso2/NASA/New_Frontiers/CAESAR/CAESAR_OMC_MIRAGE_renders_10192018/CAESAR_launch_open.ephemeris"
    SPICE_ephem_directory = "C:/Users/delliso2/Utilities/Universes/CAESAR/ephemeris_files/"
    spice_handler = ConOps.SpiceHandler(SPICE_ephem_directory)
    spice_handler.loadSpiceFiles()

    distance_report_generator = DistanceReportGenerator()

    ephemeris_file_data = distance_report_generator.ephemeris_file_reader.parseEMTGephemerisFile(CAESAR_open_ephem_file)
  
    close_approach_output_file = "C:/Users/delliso2/NASA/New_Frontiers/CAESAR/CAESAR_OMC_MIRAGE_renders_10192018/distance.report"
    spacecraft_SPICE_ID = -123
    target_body_SPICE_ID = 399
    target_body_radius = 6378.136
    distance_report_generator.detectCloseApproaches(close_approach_output_file, spacecraft_SPICE_ID, target_body_SPICE_ID, target_body_radius, [("2460999.01083915 JD TDB", "2461005.56488085 JD TDB"),("2465382.85183328 JD TDB", "2465385.55049710 JD TDB")])  

    spice_handler.unloadSpiceFiles()

    print("NOTHING")