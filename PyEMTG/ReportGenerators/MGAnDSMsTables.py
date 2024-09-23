import Mission
import MissionEvent
import numpy as np
import copy
import ConOpsPeriod as ConOps
import SpiceyPy_Utilities as SpiceyUtil
try:
    import spiceypy as spice
except:
    print("spiceypy not available") 


def generateLaunchParameterTable(table_file_name, table_caption, table_label, open_launch_event_data, middle_launch_event_data, close_launch_event_data):

    TableFile = open(table_file_name, mode = 'w')

    open_date_parts = open_launch_event_data[1].split(' ')
    middle_date_parts = middle_launch_event_data[1].split(' ')
    close_date_parts = close_launch_event_data[1].split(' ')


    launch_body = open_launch_event_data[0].parent.central_body 

    state_frame = open_launch_event_data[0].parent.state_frame

    # get the SPICE ID of the launch body
    launch_body_SPICE_ID = spice.bodn2c(launch_body.lower())
        
    # body is a small body or planet 
    if len(str(launch_body_SPICE_ID)) > 3 or (len(str(launch_body_SPICE_ID)) == 3 and str(launch_body_SPICE_ID)[1:] == '99'):
        central_body_ref_body = 'Sun'
        central_body_ref_body_ID = 10        
        
    # body is a moon
    elif len(str(launch_body_SPICE_ID)) == 3:
        central_body_ref_body = spice.bodyc2n(int(str(launch_body_SPICE_ID)[0]) * 100 + 99)
        central_body_ref_body_ID = int(str(launch_body_SPICE_ID)[0]) * 100 + 99
        


    open_launch_state, light_times = spice.spkez(launch_body_SPICE_ID, spice.str2et(str(open_launch_event_data[0].JulianDate) + " JD TDB"), "J2000", 'NONE', central_body_ref_body_ID)
    middle_launch_state, light_times = spice.spkez(launch_body_SPICE_ID, spice.str2et(str(middle_launch_event_data[0].JulianDate) + " JD TDB"), "J2000", 'NONE', central_body_ref_body_ID)
    close_launch_state, light_times = spice.spkez(launch_body_SPICE_ID, spice.str2et(str(close_launch_event_data[0].JulianDate) + " JD TDB"), "J2000", 'NONE', central_body_ref_body_ID)

    open_flyby_body_state_string =   '\multirowcell{6}{[' + str(np.round(open_launch_state[0],6)) + ',\\\\' + str(np.round(open_launch_state[1],6)) + ',\\\\'   + str(np.round(open_launch_state[2],6)) + ',\\\\'   + str(np.round(open_launch_state[3],6)) + ',\\\\'   + str(np.round(open_launch_state[4],6)) + ',\\\\'   + str(np.round(open_launch_state[5],6)) + ']}'
    middle_flyby_body_state_string = '\multirowcell{6}{[' + str(np.round(middle_launch_state[0],6)) + ',\\\\' + str(np.round(middle_launch_state[1],6)) + ',\\\\'   + str(np.round(middle_launch_state[2],6)) + ',\\\\'   + str(np.round(middle_launch_state[3],6)) + ',\\\\'   + str(np.round(middle_launch_state[4],6)) + ',\\\\'   + str(np.round(middle_launch_state[5],6)) + ']}'
    close_flyby_body_state_string =  '\multirowcell{6}{[' + str(np.round(close_launch_state[0],6)) + ',\\\\' + str(np.round(close_launch_state[1],6)) + ',\\\\'   + str(np.round(close_launch_state[2],6)) + ',\\\\'   + str(np.round(close_launch_state[3],6)) + ',\\\\'   + str(np.round(close_launch_state[4],6)) + ',\\\\'   + str(np.round(close_launch_state[5],6)) + ']}'

    open_state_string =   '\multirowcell{6}{[' + str(np.round(open_launch_event_data[0].SpacecraftState[0],6))   + ',\\\\' + str(np.round(open_launch_event_data[0].SpacecraftState[1],6)) + ',\\\\'   + str(np.round(open_launch_event_data[0].SpacecraftState[2],6)) + ',\\\\'   + str(np.round(open_launch_event_data[0].SpacecraftState[3],6)) + ',\\\\'   + str(np.round(open_launch_event_data[0].SpacecraftState[4],6)) + ',\\\\'   + str(np.round(open_launch_event_data[0].SpacecraftState[5],6)) + ']}'
    middle_state_string = '\multirowcell{6}{[' + str(np.round(middle_launch_event_data[0].SpacecraftState[0],6)) + ',\\\\' + str(np.round(middle_launch_event_data[0].SpacecraftState[1],6)) + ',\\\\' + str(np.round(middle_launch_event_data[0].SpacecraftState[2],6)) + ',\\\\' + str(np.round(middle_launch_event_data[0].SpacecraftState[3],6)) + ',\\\\' + str(np.round(middle_launch_event_data[0].SpacecraftState[4],6)) + ',\\\\' + str(np.round(middle_launch_event_data[0].SpacecraftState[5],6)) + ']}'
    close_state_string =  '\multirowcell{6}{[' + str(np.round(close_launch_event_data[0].SpacecraftState[0],6))  + ',\\\\' + str(np.round(close_launch_event_data[0].SpacecraftState[1],6)) + ',\\\\'  + str(np.round(close_launch_event_data[0].SpacecraftState[2],6)) + ',\\\\'  + str(np.round(close_launch_event_data[0].SpacecraftState[3],6)) + ',\\\\'  + str(np.round(close_launch_event_data[0].SpacecraftState[4],6)) + ',\\\\'  + str(np.round(close_launch_event_data[0].SpacecraftState[5],6)) + ']}'

    TableFile.write('       \\begin{longtable}{|l|l|l|l|l|l|l|}\n')
    TableFile.write('           \caption{' + table_caption + '} \\\\')
    TableFile.write('\hline \n')
    TableFile.write('           \\rowcolor{headergrey} \hdb{Event} & \hdb{Reference} & \hdb{Parameter} & \hdb{Open} & \hdb{Middle} & \hdb{Close} & \hdb{Comments} \\\\')
    TableFile.write('\n')
    TableFile.write('           \\rowcolor{headergrey}             & \hdb{Body}      &                 &            &              &             &                \\\\ \hline \n')
    TableFile.write('           \endfirsthead\n')
    TableFile.write('           \caption[]{' + table_caption + '} \\\\ \hline\n')
    TableFile.write('           \\rowcolor{headergrey} \hdb{Event} & \hdb{Reference} & \hdb{Parameter} & \hdb{Open} & \hdb{Middle} & \hdb{Close} & \hdb{Comments} \\\\')
    TableFile.write('\n')
    TableFile.write('           \\rowcolor{headergrey}             & \hdb{Body}      &                 &            &              &             &                \\\\ \n')
    TableFile.write('           \hline\n')
    TableFile.write('           \endhead\n')
    TableFile.write('           \label{tab:'+table_label+'}\n')

    TableFile.write('           \multirow{17}{*}{Launch}  & \multirow{17}{*}{' + launch_body + '} & Date                                                                  & ' + open_date_parts[0]+' '+open_date_parts[1]+' '+open_date_parts[2] + ' & ' + middle_date_parts[0]+' '+middle_date_parts[1]+' '+middle_date_parts[2] + ' & ' + close_date_parts[0]+' '+close_date_parts[1]+' '+close_date_parts[2] + ' & \\\\ \hhline{|~|~|-|-|-|-|-|}')                                                                                                                                                                                                                                               
    TableFile.write('                                     &                                       & Injection TDB                                                         & ' + open_date_parts[4]                                               + ' & ' + middle_date_parts[4]                                                   + ' & ' + close_date_parts[4]                                                 + ' & \\\\ \hhline{|~|~|-|-|-|-|-|}')
    TableFile.write('                                     &                                       & \multirowcell{6}{Cartesian \\\\ state [km, km/s]}                     & ' + open_state_string +                                                ' & ' + middle_state_string +                                                    ' & ' + close_state_string +                                                  ' & \\\\ \hhline{|~|~|~|~|~|~|~|}')
    TableFile.write('                                     &                                       &                                                                       & ' +                                                                    ' & ' +                                                                          ' & ' +                                                                       ' & \\\\ \hhline{|~|~|~|~|~|~|~|}')
    TableFile.write('                                     &                                       &                                                                       & ' +                                                                    ' & ' +                                                                          ' & ' +                                                                       ' & ' + frameNameTranslator(state_frame, launch_body) + ' \\\\ \hhline{|~|~|~|~|~|~|~|}')
    TableFile.write('                                     &                                       &                                                                       & ' +                                                                    ' & ' +                                                                          ' & ' +                                                                       ' & \\\\ \hhline{|~|~|~|~|~|~|~|}')
    TableFile.write('                                     &                                       &                                                                       & ' +                                                                    ' & ' +                                                                          ' & ' +                                                                       ' & \\\\ \hhline{|~|~|~|~|~|~|~|}')
    TableFile.write('                                     &                                       &                                                                       & ' +                                                                    ' & ' +                                                                          ' & ' +                                                                       ' & \\\\ \hhline{|~|~|-|-|-|-|-|}')
    TableFile.write('                                     &                                       & \multirowcell{6}{' + launch_body + ' Cartesian \\\\ state [km, km/s]} & ' + open_flyby_body_state_string +                                     ' & ' + middle_flyby_body_state_string +                                         ' & ' + close_flyby_body_state_string +                                       ' & \\\\ \hhline{|~|~|~|~|~|~|~|}')
    TableFile.write('                                     &                                       &                                                                       & ' +                                                                    ' & ' +                                                                          ' & ' +                                                                       ' & \\\\ \hhline{|~|~|~|~|~|~|~|}')
    TableFile.write('                                     &                                       &                                                                       & ' +                                                                    ' & ' +                                                                          ' & ' +                                                                       ' & ' + frameNameTranslator(state_frame, central_body_ref_body) + ' \\\\ \hhline{|~|~|~|~|~|~|~|}')
    TableFile.write('                                     &                                       &                                                                       & ' +                                                                    ' & ' +                                                                          ' & ' +                                                                       ' & \\\\ \hhline{|~|~|~|~|~|~|~|}')
    TableFile.write('                                     &                                       &                                                                       & ' +                                                                    ' & ' +                                                                          ' & ' +                                                                       ' & \\\\ \hhline{|~|~|~|~|~|~|~|}')
    TableFile.write('                                     &                                       &                                                                       & ' +                                                                    ' & ' +                                                                          ' & ' +                                                                       ' & \\\\ \hhline{|~|~|-|-|-|-|-|}')
    TableFile.write('                                     &                                       & C3 $[\\text{km}^2/\\text{s}^2]$                                       & ' + str(np.round(open_launch_event_data[2],6)) +                       ' & ' + str(np.round(middle_launch_event_data[2],6)) +                           ' & ' + str(np.round(close_launch_event_data[2],6)) +                         ' & \\\\ \hhline{|~|~|-|-|-|-|-|}')
    TableFile.write('                                     &                                       & DLA [degrees]                                                         & ' + str(np.round(open_launch_event_data[4],6)) +                       ' & ' + str(np.round(middle_launch_event_data[4],6)) +                           ' & ' + str(np.round(close_launch_event_data[4],6)) +                         ' & \multirowcell{2}{'+ frameNameTranslator(state_frame, launch_body) +'} \\\\ \hhline{|~|~|-|-|-|-|~|}')
    TableFile.write('                                     &                                       & RLA [degrees]                                                         & ' + str(np.round(open_launch_event_data[3],6)) +                       ' & ' + str(np.round(middle_launch_event_data[3],6)) +                           ' & ' + str(np.round(close_launch_event_data[3],6)) +                         ' & \\\\ \hline')


    TableFile.write('     \end{longtable}\n')

def generateFlybyComparisonTable(table_file_name, table_caption, table_label, open_periapse_event_data, middle_periapse_event_data, close_periapse_event_data):

    TableFile = open(table_file_name, mode = 'w')

    TableFile.write('       \\begin{longtable}{|l|l|l|l|l|l|l|}\n')
    TableFile.write('           \caption{' + table_caption + '} \\\\')
    TableFile.write('\hline \n')
    TableFile.write('           \\rowcolor{headergrey} \hdb{Event} & \hdb{Reference} & \hdb{Parameter} & \hdb{Open} & \hdb{Middle} & \hdb{Close} & \hdb{Comments} \\\\')
    TableFile.write('\n')
    TableFile.write('           \\rowcolor{headergrey}             & \hdb{Body}      &                 &            &              &             &                \\\\ \hline')
    TableFile.write('\n')
    TableFile.write('           \endfirsthead\n')
    TableFile.write('           \caption[]{' + table_caption + '} \\\\ \hline')
    TableFile.write('           \\rowcolor{headergrey} \hdb{Event} & \hdb{Reference} & \hdb{Parameter} & \hdb{Open} & \hdb{Middle} & \hdb{Close} & \hdb{Comments} \\\\')
    TableFile.write('\n')
    TableFile.write('           \\rowcolor{headergrey}             & \hdb{Body}      &                 &            &              &             &                \\\\')    
    TableFile.write('           \hline\n')
    TableFile.write('           \endhead\n')
    TableFile.write('           \label{tab:'+table_label+'}\n')

    for index in range(0, len(open_periapse_event_data)):

        open_date_parts = open_periapse_event_data[index][1].split(' ')
        middle_date_parts = middle_periapse_event_data[index][1].split(' ')
        close_date_parts = close_periapse_event_data[index][1].split(' ')
        open_state_string =   '\multirowcell{6}{[' + str(np.round(open_periapse_event_data[index][0].SpacecraftState[0],6))   + ',\\\\' + str(np.round(open_periapse_event_data[index][0].SpacecraftState[1],6)) + ',\\\\'   + str(np.round(open_periapse_event_data[index][0].SpacecraftState[2],6)) + ',\\\\'   + str(np.round(open_periapse_event_data[index][0].SpacecraftState[3],6)) + ',\\\\'   + str(np.round(open_periapse_event_data[index][0].SpacecraftState[4],6)) + ',\\\\'   + str(np.round(open_periapse_event_data[index][0].SpacecraftState[5],6)) + ']}'
        middle_state_string = '\multirowcell{6}{[' + str(np.round(middle_periapse_event_data[index][0].SpacecraftState[0],6)) + ',\\\\' + str(np.round(middle_periapse_event_data[index][0].SpacecraftState[1],6)) + ',\\\\' + str(np.round(middle_periapse_event_data[index][0].SpacecraftState[2],6)) + ',\\\\' + str(np.round(middle_periapse_event_data[index][0].SpacecraftState[3],6)) + ',\\\\' + str(np.round(middle_periapse_event_data[index][0].SpacecraftState[4],6)) + ',\\\\' + str(np.round(middle_periapse_event_data[index][0].SpacecraftState[5],6)) + ']}'
        close_state_string =  '\multirowcell{6}{[' + str(np.round(close_periapse_event_data[index][0].SpacecraftState[0],6))  + ',\\\\' + str(np.round(close_periapse_event_data[index][0].SpacecraftState[1],6)) + ',\\\\'  + str(np.round(close_periapse_event_data[index][0].SpacecraftState[2],6)) + ',\\\\'  + str(np.round(close_periapse_event_data[index][0].SpacecraftState[3],6)) + ',\\\\'  + str(np.round(close_periapse_event_data[index][0].SpacecraftState[4],6)) + ',\\\\'  + str(np.round(close_periapse_event_data[index][0].SpacecraftState[5],6)) + ']}'

        BdotR_open, BdotT_open, B_angle_open =  filterFormatBplaneParams(open_periapse_event_data[index])   
        BdotR_middle, BdotT_middle, B_angle_middle =  filterFormatBplaneParams(middle_periapse_event_data[index])
        BdotR_close, BdotT_close, B_angle_close =  filterFormatBplaneParams(close_periapse_event_data[index])  

        flyby_body = open_periapse_event_data[index][0].parent.central_body 
        state_frame = open_periapse_event_data[index][0].parent.state_frame

        # get central body state
        flyby_body_SPICE_ID = spice.bodn2c(flyby_body.lower())
        
        # body is a planet or small body
        if len(str(flyby_body_SPICE_ID)) > 3 or (len(str(flyby_body_SPICE_ID)) == 3 and str(flyby_body_SPICE_ID)[1:] == '99'):
            central_body_ref_body = 'Sun'
            central_body_ref_body_ID = 10
        # body is a moon
        elif len(str(flyby_body_SPICE_ID)) == 3:
            central_body_ref_body = spice.bodyc2n(int(str(flyby_body_SPICE_ID)[0]) * 100 + 99)
            central_body_ref_body_ID = int(str(flyby_body_SPICE_ID)[0]) * 100 + 99
            
        open_flyby_state, light_times = spice.spkez(flyby_body_SPICE_ID, spice.str2et(str(open_periapse_event_data[index][0].JulianDate) + " JD TDB"), "J2000", 'NONE', central_body_ref_body_ID)
        middle_flyby_state, light_times = spice.spkez(flyby_body_SPICE_ID, spice.str2et(str(middle_periapse_event_data[index][0].JulianDate) + " JD TDB"), "J2000", 'NONE', central_body_ref_body_ID)
        close_flyby_state, light_times = spice.spkez(flyby_body_SPICE_ID, spice.str2et(str(close_periapse_event_data[index][0].JulianDate) + " JD TDB"), "J2000", 'NONE', central_body_ref_body_ID)
    
        open_flyby_state_string =   '\multirowcell{6}{[' + str(np.round(open_flyby_state[0],6))   + ',\\\\' + str(np.round(open_flyby_state[1],6)) + ',\\\\'   + str(np.round(open_flyby_state[2],6)) + ',\\\\'   + str(np.round(open_flyby_state[3],6)) + ',\\\\'   + str(np.round(open_flyby_state[4],6)) + ',\\\\'   + str(np.round(open_flyby_state[5],6)) + ']}'
        middle_flyby_state_string = '\multirowcell{6}{[' + str(np.round(middle_flyby_state[0],6))   + ',\\\\' + str(np.round(middle_flyby_state[1],6)) + ',\\\\'   + str(np.round(middle_flyby_state[2],6)) + ',\\\\'   + str(np.round(middle_flyby_state[3],6)) + ',\\\\'   + str(np.round(middle_flyby_state[4],6)) + ',\\\\'   + str(np.round(middle_flyby_state[5],6)) + ']}'
        close_flyby_state_string =  '\multirowcell{6}{[' + str(np.round(close_flyby_state[0],6))   + ',\\\\' + str(np.round(close_flyby_state[1],6)) + ',\\\\'   + str(np.round(close_flyby_state[2],6)) + ',\\\\'   + str(np.round(close_flyby_state[3],6)) + ',\\\\'   + str(np.round(close_flyby_state[4],6)) + ',\\\\'   + str(np.round(close_flyby_state[5],6)) + ']}'

        if index > 0:
            TableFile.write('           \pagebreak \n')
        
        TableFile.write('           \multirow{19}{*}{Flyby}   & \multirow{19}{*}{' + flyby_body + '} & Date                                                                & ' + open_date_parts[0]+' '+open_date_parts[1]+' '+open_date_parts[2] + ' & ' + middle_date_parts[0]+' '+middle_date_parts[1]+' '+middle_date_parts[2] + ' & ' + close_date_parts[0]+' '+close_date_parts[1]+' '+close_date_parts[2] + ' & \\\\ \hhline{|~|~|-|-|-|-|-|}')                                                                                                                                                                                                                                               
        TableFile.write('                                     &                                      & Close approach TDB                                                  & ' + open_date_parts[4]                                               + ' & ' + middle_date_parts[4]                                                   + ' & ' + close_date_parts[4]                                                 + ' & \\\\ \hhline{|~|~|-|-|-|-|-|}')
        TableFile.write('                                     &                                      & \multirowcell{6}{Periapse Cartesian \\\\ state [km, km/s]}          & ' + open_state_string +                                                ' & ' + middle_state_string +                                                    ' & ' + close_state_string +                                                  ' & \\\\ \hhline{|~|~|~|~|~|~|~|}')
        TableFile.write('                                     &                                      &                                                                     & ' +                                                                    ' & ' +                                                                          ' & ' +                                                                       ' & \\\\ \hhline{|~|~|~|~|~|~|~|}')
        TableFile.write('                                     &                                      &                                                                     & ' +                                                                    ' & ' +                                                                          ' & ' +                                                                       ' & ' + frameNameTranslator(state_frame, flyby_body) + ' \\\\ \hhline{|~|~|~|~|~|~|~|}')
        TableFile.write('                                     &                                      &                                                                     & ' +                                                                    ' & ' +                                                                          ' & ' +                                                                       ' & \\\\ \hhline{|~|~|~|~|~|~|~|}')
        TableFile.write('                                     &                                      &                                                                     & ' +                                                                    ' & ' +                                                                          ' & ' +                                                                       ' & \\\\ \hhline{|~|~|~|~|~|~|~|}')
        TableFile.write('                                     &                                      &                                                                     & ' +                                                                    ' & ' +                                                                          ' & ' +                                                                       ' & \\\\ \hhline{|~|~|-|-|-|-|-|}')
        TableFile.write('                                     &                                      & \multirowcell{6}{'+ flyby_body + ' Cartesian \\\\ state [km, km/s]} & ' + open_flyby_state_string +                                          ' & ' + middle_flyby_state_string +                                              ' & ' + close_flyby_state_string +                                            ' & \\\\ \hhline{|~|~|~|~|~|~|~|}')
        TableFile.write('                                     &                                      &                                                                     & ' +                                                                    ' & ' +                                                                          ' & ' +                                                                       ' & \\\\ \hhline{|~|~|~|~|~|~|~|}')
        TableFile.write('                                     &                                      &                                                                     & ' +                                                                    ' & ' +                                                                          ' & ' +                                                                       ' & ' + frameNameTranslator(state_frame, central_body_ref_body) + ' \\\\ \hhline{|~|~|~|~|~|~|~|}')
        TableFile.write('                                     &                                      &                                                                     & ' +                                                                    ' & ' +                                                                          ' & ' +                                                                       ' & \\\\ \hhline{|~|~|~|~|~|~|~|}')
        TableFile.write('                                     &                                      &                                                                     & ' +                                                                    ' & ' +                                                                          ' & ' +                                                                       ' & \\\\ \hhline{|~|~|~|~|~|~|~|}')
        TableFile.write('                                     &                                      &                                                                     & ' +                                                                    ' & ' +                                                                          ' & ' +                                                                       ' & \\\\ \hhline{|~|~|-|-|-|-|-|}')
        TableFile.write('                                     &                                      & Flyby C3 $[\\text{km}/\\text{s}]$                                   & ' + str(np.round(open_periapse_event_data[index][2],6)) +              ' & ' + str(np.round(middle_periapse_event_data[index][2],6)) +                  ' & ' + str(np.round(close_periapse_event_data[index][2],6)) +                ' & \\\\ \hhline{|~|~|-|-|-|-|-|}')
        TableFile.write('                                     &                                      & Flyby altitude $[\\text{km}]$                                       & ' + str(np.round(open_periapse_event_data[index][0].Altitude,3)) +     ' & ' + str(np.round(middle_periapse_event_data[index][0].Altitude,3)) +         ' & ' + str(np.round(close_periapse_event_data[index][0].Altitude,3)) +       ' & \\\\ \hhline{|~|~|-|-|-|-|-|}')
        TableFile.write('                                     &                                      & BR [km]                                                             & ' + str(BdotR_open) +                                                  ' & ' + str(BdotR_middle) +                                                      ' & ' + str(BdotR_close) +                                                    ' & \multirowcell{3}{Reference normal:\\\\'+ frameNameTranslator(state_frame, flyby_body) +'-Z vector} \\\\ \hhline{|~|~|-|-|-|-|~|}')
        TableFile.write('                                     &                                      & BT [km]                                                             & ' + str(BdotT_open) +                                                  ' & ' + str(BdotT_middle) +                                                      ' & ' + str(BdotT_close) +                                                    ' & \\\\ \hhline{|~|~|-|-|-|-|~|}')
        TableFile.write('                                     &                                      & B-plane angle [degrees]                                             & ' + str(B_angle_open) +                                                ' & ' + str(B_angle_middle) +                                                    ' & ' + str(B_angle_close) +                                                  ' & \\\\ \hline \n')



    TableFile.write('     \end{longtable}\n')

def generateFlybyTargetTable(table_file_name, table_caption, table_label, periapse_event_data):
    
    TableFile = open(table_file_name, mode = 'w')

    TableFile.write('       \\begin{longtable}{|c|c|c|c|c|c|c|c|c|}\n')
    TableFile.write('           \caption{' + table_caption + '} \\\\')
    TableFile.write('\hline \n')
    TableFile.write('           \\rowcolor{headergrey} \hdb{Flyby Event} & \hdb{Epoch} & \hdb{Altitude} & \hdb{C3}                           & \hdb{BR}   & \hdb{BT}   & \hdb{B-plane Angle}     & \hdb{SEP}       & \hdb{SPE}       \\\\\n')
    TableFile.write('           \\rowcolor{headergrey}                   & \hdb{(TDB)} & \hdb{[km]}     & \hdb{$[\\text{km}^2/\\text{s}^2]$}  & \hdb{[km]} & \hdb{[km]} & \hdb{[degrees]} & \hdb{[degrees]} & \hdb{[degrees]} \\\\')
    TableFile.write('\n')
    TableFile.write('           \endfirsthead\n')
    TableFile.write('           \caption[]{' + table_caption + '} \\\\ \hline')
    TableFile.write('           \\rowcolor{headergrey} \hdb{Flyby Event} & \hdb{Epoch} & \hdb{Altitude} & \hdb{C3}                           & \hdb{BR}   & \hdb{BT}   & \hdb{B-plane Angle}     & \hdb{SEP}       & \hdb{SPE}       \\\\\n')
    TableFile.write('           \\rowcolor{headergrey}                   & \hdb{(TDB)} & \hdb{[km]}     & \hdb{$[\\text{km}^2/\\text{s}^2$]} & \hdb{[km]} & \hdb{[km]} & \hdb{[degrees]} & \hdb{[degrees]} & \hdb{[degrees]} \\\\')
    TableFile.write('           \hline\n')
    TableFile.write('           \endhead\n')
    TableFile.write('           \label{tab:'+table_label+'}\n')

    for event in periapse_event_data:

        date_parts = event[1].split(' ')
        date_string = '\makecell{'+date_parts[0]+' '+date_parts[1]+' '+date_parts[2]+'\\\\'+date_parts[4]+'}'

        BdotR, BdotT, B_angle = filterFormatBplaneParams(event)

        TableFile.write('               ' + event[0].parent.central_body + ' & '
                                          + date_string + ' & '                                            
                                          + str(np.round(event[0].Altitude,3)) + ' & '  
                                          + str(np.round(event[4],6)) + ' & '
                                          + BdotR + ' & '
                                          + BdotT + ' & '
                                          + B_angle + ' & '
                                          + str(np.round(event[2],6)) + ' & '
                                          + str(np.round(event[3],6)) + '\\\\ \hline\n')

    TableFile.write('     \end{longtable}\n')

def generateFlybySpacecraftStateTable(table_file_name, table_caption, table_label, periapse_event_data):

    TableFile = open(table_file_name, mode = 'w')

    TableFile.write('       \\begin{longtable}{|c|c|c|c|c|c|c|c|}\n')
    TableFile.write('           \caption{' + table_caption + '} \\\\')
    TableFile.write('\hline \n')
    TableFile.write('           \\rowcolor{headergrey} \hdb{Flyby Event} & \hdb{Epoch} & \hdb{x}    & \hdb{y}    & \hdb{z}    & \hdb{$\\text{v}_x$} & \hdb{$\\text{v}_y$} & \hdb{$\\text{v}_z$} \\\\\n')
    TableFile.write('           \\rowcolor{headergrey}                   & \hdb{(TDB)} & \hdb{[km]} & \hdb{[km]} & \hdb{[km]} & \hdb{[km/s]}        & \hdb{[km/s]}        & \hdb{[km/s]}        \\\\')
    TableFile.write('\n')
    TableFile.write('           \endfirsthead\n')
    TableFile.write('           \caption[]{' + table_caption + '} \\\\ \hline')
    TableFile.write('           \\rowcolor{headergrey} \hdb{Flyby Event} & \hdb{Epoch} & \hdb{x}    & \hdb{y}    & \hdb{z}    & \hdb{$\\text{v}_x$} & \hdb{$\\text{v}_y$} & \hdb{$\\text{v}_z$} \\\\\n')
    TableFile.write('           \\rowcolor{headergrey}                   & \hdb{(TDB)} & \hdb{[km]} & \hdb{[km]} & \hdb{[km]} & \hdb{[km/s]}        & \hdb{[km/s]}        & \hdb{[km/s]}        \\\\')
    TableFile.write('           \hline\n')
    TableFile.write('           \endhead\n')
    TableFile.write('           \label{tab:'+table_label+'}\n')

    for event in periapse_event_data:

        date_parts = event[1].split(' ')
        date_string = '\makecell{'+date_parts[0]+' '+date_parts[1]+' '+date_parts[2]+'\\\\'+date_parts[4]+'}'

        TableFile.write('               ' + event[0].parent.central_body + ' & '
                                          + date_string + ' & '                                            
                                          + str(np.round(event[0].SpacecraftState[0],6)) + ' & '  
                                          + str(np.round(event[0].SpacecraftState[1],6)) + ' & ' 
                                          + str(np.round(event[0].SpacecraftState[2],6)) + ' & '
                                          + str(np.round(event[0].SpacecraftState[3],6)) + ' & '  
                                          + str(np.round(event[0].SpacecraftState[4],6)) + ' & ' 
                                          + str(np.round(event[0].SpacecraftState[5],6)) + '\\\\ \hline\n')

    TableFile.write('     \end{longtable}\n')

def generateFlybyBodyStateTable(table_file_name, table_caption, table_label, periapse_event_data):

    TableFile = open(table_file_name, mode = 'w')

    TableFile.write('       \\begin{longtable}{|c|c|c|c|c|c|c|c|}\n')
    TableFile.write('           \caption{' + table_caption + '} \\\\')
    TableFile.write('\hline \n')
    TableFile.write('           \\rowcolor{headergrey} \hdb{Flyby Event} & \hdb{Epoch} & \hdb{x}    & \hdb{y}    & \hdb{z}    & \hdb{$\\text{v}_x$} & \hdb{$\\text{v}_y$} & \hdb{$\\text{v}_z$} \\\\\n')
    TableFile.write('           \\rowcolor{headergrey}                   & \hdb{(TDB)} & \hdb{[km]} & \hdb{[km]} & \hdb{[km]} & \hdb{[km/s]}        & \hdb{[km/s]}        & \hdb{[km/s]}        \\\\')
    TableFile.write('\hline \n')
    TableFile.write('           \endfirsthead\n')
    TableFile.write('           \caption[]{' + table_caption + '} \\\\ \hline')
    TableFile.write('           \\rowcolor{headergrey} \hdb{Flyby Event} & \hdb{Epoch} & \hdb{x}    & \hdb{y}    & \hdb{z}    & \hdb{$\\text{v}_x$} & \hdb{$\\text{v}_y$} & \hdb{$\\text{v}_z$} \\\\\n')
    TableFile.write('           \\rowcolor{headergrey}                   & \hdb{(TDB)} & \hdb{[km]} & \hdb{[km]} & \hdb{[km]} & \hdb{[km/s]}        & \hdb{[km/s]}        & \hdb{[km/s]}        \\\\')
    TableFile.write('           \hline\n')
    TableFile.write('           \endhead\n')
    TableFile.write('           \label{tab:'+table_label+'}\n')

    for event in periapse_event_data:

        central_body = event[0].parent.central_body 

        # get central body state
        flyby_body_SPICE_ID = spice.bodn2c(central_body.lower())
        
        state = []
        # body is a planet or small body
        if len(str(flyby_body_SPICE_ID)) > 3 or (len(str(flyby_body_SPICE_ID)) == 3 and str(flyby_body_SPICE_ID)[1:] == '99'):
            central_body_ref_body_ID = 10
        # body is a moon
        elif len(str(flyby_body_SPICE_ID)) == 3:
            central_body_ref_body_ID = int(str(flyby_body_SPICE_ID)[0]) * 100 + 99

        flyby_body_state, light_times = spice.spkez(flyby_body_SPICE_ID, spice.str2et(str(event[0].JulianDate) + " JD TDB"), "J2000", 'NONE', central_body_ref_body_ID)

        date_parts = event[1].split(' ')
        date_string = '\makecell{'+date_parts[0]+' '+date_parts[1]+' '+date_parts[2]+'\\\\'+date_parts[4]+'}'

        TableFile.write('               ' + event[0].parent.central_body + ' & '
                                          + date_string + ' & '                                            
                                          + str(np.round(flyby_body_state[0],6)) + ' & '  
                                          + str(np.round(flyby_body_state[1],6)) + ' & ' 
                                          + str(np.round(flyby_body_state[2],6)) + ' & '
                                          + str(np.round(flyby_body_state[3],6)) + ' & '  
                                          + str(np.round(flyby_body_state[4],6)) + ' & ' 
                                          + str(np.round(flyby_body_state[5],6)) + '\\\\ \hline\n')

    TableFile.write('     \end{longtable}\n')


def generateManeuverComparisonTable(table_file_name, table_caption, table_label, open_maneuver_event_data, middle_maneuver_event_data, close_maneuver_event_data):

    TableFile = open(table_file_name, mode = 'w')

    TableFile.write('       \\begin{longtable}{|l|l|l|l|l|l|l|}\n')
    TableFile.write('           \caption{' + table_caption + '} \\\\')
    TableFile.write('\hline \n')
    TableFile.write('           \\rowcolor{headergrey} \hdb{Event} & \hdb{Reference} & \hdb{Parameter} & \hdb{Open} & \hdb{Middle} & \hdb{Close} & \hdb{Comments} \\\\')
    TableFile.write('\n')
    TableFile.write('           \\rowcolor{headergrey}             & \hdb{Body}      &                 &            &              &             &                \\\\ \hline')
    TableFile.write('\n')
    TableFile.write('           \endfirsthead\n')
    TableFile.write('           \caption[]{' + table_caption + '} \\\\ \hline')
    TableFile.write('           \\rowcolor{headergrey} \hdb{Event} & \hdb{Reference} & \hdb{Parameter} & \hdb{Open} & \hdb{Middle} & \hdb{Close} & \hdb{Comments} \\\\')
    TableFile.write('\n')
    TableFile.write('           \\rowcolor{headergrey}             & \hdb{Body}      &                 &            &              &             &                \\\\')
    TableFile.write('\n')
    TableFile.write('           \hline\n')
    TableFile.write('           \endhead\n')
    TableFile.write('           \label{tab:'+table_label+'}\n')

    for index in range(0, len(open_maneuver_event_data)):

        if open_maneuver_event_data[index][1] != '-':
            open_date_parts = open_maneuver_event_data[index][1].split(' ')
            open_date = open_date_parts[0]+' '+open_date_parts[1]+' '+open_date_parts[2]
            open_time = open_date_parts[4]
            open_state_string =   '\multirowcell{6}{[' + str(np.round(open_maneuver_event_data[index][0].SpacecraftState[0],6))   + ',\\\\' + str(np.round(open_maneuver_event_data[index][0].SpacecraftState[1],6)) + ',\\\\'   + str(np.round(open_maneuver_event_data[index][0].SpacecraftState[2],6)) + ',\\\\'   + str(np.round(open_maneuver_event_data[index][0].SpacecraftState[3] - open_maneuver_event_data[index][0].DeltaVorThrustVectorControl[0],6)) + ',\\\\'   + str(np.round(open_maneuver_event_data[index][0].SpacecraftState[4] - open_maneuver_event_data[index][0].DeltaVorThrustVectorControl[1],6)) + ',\\\\'   + str(np.round(open_maneuver_event_data[index][0].SpacecraftState[5] - open_maneuver_event_data[index][0].DeltaVorThrustVectorControl[2],6)) + ']}'
            open_maneuver_string =   '\multirowcell{3}{[' + str(np.round(1000*open_maneuver_event_data[index][0].DeltaVorThrustVectorControl[0],6))   + ',\\\\' + str(np.round(1000*open_maneuver_event_data[index][0].DeltaVorThrustVectorControl[1],6)) + ',\\\\'   + str(np.round(1000*open_maneuver_event_data[index][0].DeltaVorThrustVectorControl[2],6))   +']}'
            open_SEP = str(np.round(open_maneuver_event_data[index][2],6))
            open_SPE = str(np.round(open_maneuver_event_data[index][3],6))
            open_DVmag = str(np.round(1000*open_maneuver_event_data[index][0].DVmagorThrottle,6))
            central_body_string = open_maneuver_event_data[index][0].parent.central_body
            state_frame = open_maneuver_event_data[index][0].parent.state_frame
            state_string = frameNameTranslator(state_frame, central_body_string)
        else:
            open_date_parts = ['','','','','']
            open_date = '-'
            open_time = '-'
            open_state_string = '\multirowcell{6}{-}'
            open_maneuver_string = '\multirowcell{3}{-}'
            open_SEP = '-'
            open_SPE = '-'
            open_DVmag = '-'

        if middle_maneuver_event_data[index][1] != '-':
            middle_date_parts = middle_maneuver_event_data[index][1].split(' ')
            middle_date = middle_date_parts[0]+' '+middle_date_parts[1]+' '+middle_date_parts[2]
            middle_time = middle_date_parts[4]
            middle_state_string = '\multirowcell{6}{[' + str(np.round(middle_maneuver_event_data[index][0].SpacecraftState[0],6)) + ',\\\\' + str(np.round(middle_maneuver_event_data[index][0].SpacecraftState[1],6)) + ',\\\\' + str(np.round(middle_maneuver_event_data[index][0].SpacecraftState[2],6)) + ',\\\\' + str(np.round(middle_maneuver_event_data[index][0].SpacecraftState[3] - middle_maneuver_event_data[index][0].DeltaVorThrustVectorControl[0],6)) + ',\\\\' + str(np.round(middle_maneuver_event_data[index][0].SpacecraftState[4] - middle_maneuver_event_data[index][0].DeltaVorThrustVectorControl[1],6)) + ',\\\\' + str(np.round(middle_maneuver_event_data[index][0].SpacecraftState[5] - middle_maneuver_event_data[index][0].DeltaVorThrustVectorControl[2],6)) + ']}'        
            middle_maneuver_string = '\multirowcell{3}{[' + str(np.round(1000*middle_maneuver_event_data[index][0].DeltaVorThrustVectorControl[0],6)) + ',\\\\' + str(np.round(1000*middle_maneuver_event_data[index][0].DeltaVorThrustVectorControl[1],6)) + ',\\\\' + str(np.round(1000*middle_maneuver_event_data[index][0].DeltaVorThrustVectorControl[2],6)) + ']}'
            middle_SEP = str(np.round(middle_maneuver_event_data[index][2],6))
            middle_SPE = str(np.round(middle_maneuver_event_data[index][3],6))
            middle_DVmag = str(np.round(1000*middle_maneuver_event_data[index][0].DVmagorThrottle,6))
            central_body_string = middle_maneuver_event_data[index][0].parent.central_body
            state_frame = middle_maneuver_event_data[index][0].parent.state_frame
            state_string = frameNameTranslator(state_frame, central_body_string)
        else:
            middle_date_parts = ['','','','','']
            middle_date = '-'
            middle_time = '-'
            middle_state_string = '\multirowcell{6}{-}'
            middle_maneuver_string = '\multirowcell{3}{-}'
            middle_SEP = '-'
            middle_SPE = '-'
            middle_DVmag = '-'

        if close_maneuver_event_data[index][1] != '-':
            close_date_parts = close_maneuver_event_data[index][1].split(' ')
            close_date = close_date_parts[0]+' '+close_date_parts[1]+' '+close_date_parts[2]
            close_time = close_date_parts[4]
            close_state_string =  '\multirowcell{6}{[' + str(np.round(close_maneuver_event_data[index][0].SpacecraftState[0],6))  + ',\\\\' + str(np.round(close_maneuver_event_data[index][0].SpacecraftState[1],6)) + ',\\\\'  + str(np.round(close_maneuver_event_data[index][0].SpacecraftState[2],6)) + ',\\\\'  + str(np.round(close_maneuver_event_data[index][0].SpacecraftState[3] - close_maneuver_event_data[index][0].DeltaVorThrustVectorControl[0],6)) + ',\\\\'  + str(np.round(close_maneuver_event_data[index][0].SpacecraftState[4] - close_maneuver_event_data[index][0].DeltaVorThrustVectorControl[1],6)) + ',\\\\'  + str(np.round(close_maneuver_event_data[index][0].SpacecraftState[5] - close_maneuver_event_data[index][0].DeltaVorThrustVectorControl[2],6)) + ']}'   
            close_maneuver_string =  '\multirowcell{3}{[' + str(np.round(1000*close_maneuver_event_data[index][0].DeltaVorThrustVectorControl[0],6))  + ',\\\\' + str(np.round(1000*close_maneuver_event_data[index][0].DeltaVorThrustVectorControl[1],6)) + ',\\\\'  + str(np.round(1000*close_maneuver_event_data[index][0].DeltaVorThrustVectorControl[2],6))  + ']}'
            close_SEP = str(np.round(close_maneuver_event_data[index][2],6))
            close_SPE = str(np.round(close_maneuver_event_data[index][3],6))
            close_DVmag = str(np.round(1000*close_maneuver_event_data[index][0].DVmagorThrottle,6))
            central_body_string = close_maneuver_event_data[index][0].parent.central_body
            state_frame = close_maneuver_event_data[index][0].parent.state_frame
            state_string = frameNameTranslator(state_frame, central_body_string)
        else:
            close_date_parts = ['','','','','']
            close_date = '-'
            close_time = '-'
            close_state_string = '\multirowcell{6}{-}'
            close_maneuver_string = '\multirowcell{3}{-}'
            close_SEP = '-'
            close_SPE = '-'
            close_DVmag = '-'
            
        if index > 0:
            TableFile.write('           \pagebreak \n')

        TableFile.write('        \multirow{14}{*}{Maneuver} & \multirow{14}{*}{' + central_body_string + '} & Date                                                               & ' + open_date +            ' & ' + middle_date +            ' & ' + close_date +            ' & \\\\ \hhline{|~|~|-|-|-|-|-|}')                                                                                                                                                                                                                                               
        TableFile.write('                                   &                                               & Maneuver TDB                                                       & ' + open_time +            ' & ' + middle_time +            ' & ' + close_time +            ' & \\\\ \hhline{|~|~|-|-|-|-|-|}')
        TableFile.write('                                   &                                               & \multirowcell{6}{Pre-maneuver Cartesian \\\\ state [km, km/s]}     & ' + open_state_string +    ' & ' + middle_state_string +    ' & ' + close_state_string +    ' & \\\\ \hhline{|~|~|~|~|~|~|~|}')
        TableFile.write('                                   &                                               &                                                                    &                              &                                &                               & \\\\ \hhline{|~|~|~|~|~|~|~|}')
        TableFile.write('                                   &                                               &                                                                    &                              &                                &                               & ' + state_string + ' \\\\ \hhline{|~|~|~|~|~|~|~|}')
        TableFile.write('                                   &                                               &                                                                    &                              &                                &                               & \\\\ \hhline{|~|~|~|~|~|~|~|}')
        TableFile.write('                                   &                                               &                                                                    &                              &                                &                               & \\\\ \hhline{|~|~|~|~|~|~|~|}')
        TableFile.write('                                   &                                               &                                                                    &                              &                                &                               & \\\\ \hhline{|~|~|-|-|-|-|-|}')
        TableFile.write('                                   &                                               & \multirowcell{3}{$\Delta \\text{v vector } [\\text{m}/\\text{s}]$} & ' + open_maneuver_string + ' & ' + middle_maneuver_string + ' & ' + close_maneuver_string + ' & \\\\ \hhline{|~|~|~|~|~|~|~|}')
        TableFile.write('                                   &                                               &                                                                    &                              &                                &                               & ' + state_string + ' \\\\ \hhline{|~|~|~|~|~|~|~|}')
        TableFile.write('                                   &                                               &                                                                    &                              &                                &                               & \\\\ \hhline{|~|~|-|-|-|-|-|}')
        TableFile.write('                                   &                                               & $\Delta \\text{v magnitude } [\\text{m}/\\text{s}]$                & ' + open_DVmag +           ' & ' + middle_DVmag +           ' & ' + close_DVmag +           ' & \\\\ \hhline{|~|~|-|-|-|-|-|}')
        TableFile.write('                                   &                                               & SEP angle [degrees]                                                & ' + open_SEP +             ' & ' + middle_SEP +             ' & ' + close_SEP +             ' & \\\\ \hhline{|~|~|-|-|-|-|-|}')
        TableFile.write('                                   &                                               & SPE angle [degrees]                                                & ' + open_SPE +             ' & ' + middle_SPE +             ' & ' + close_SPE +             ' & \\\\ \hline \n')
        



    TableFile.write('     \end{longtable}\n')

def generateDeltaVTable(table_file_name, table_caption, table_label, maneuver_event_data):
    
    TableFile = open(table_file_name, mode = 'w')

    TableFile.write('       \\begin{longtable}{|c|c|c|c|c|c|c|c|}\n')
    TableFile.write('           \caption{' + table_caption + '} \\\\')
    TableFile.write('\hline \n')
    TableFile.write('           \\rowcolor{headergrey} \hdb{Maneuver Event} & \hdb{Epoch} & \hdb{Magnitude} & \hdb{$\\Delta \\text{v}_x$} & \hdb{$\\Delta \\text{v}_y$} & \hdb{$\\Delta \\text{v}_z$} & \hdb{SEP}       & \hdb{SPE}       \\\\\n')
    TableFile.write('           \\rowcolor{headergrey}                      & \hdb{(TDB)} & \hdb{[m/s]}     & \hdb{[m/s]}                 & \hdb{[m/s]}                 & \hdb{[m/s]}                 & \hdb{[degrees]} & \hdb{[degrees]} \\\\')
    TableFile.write('\n')
    TableFile.write('           \endfirsthead\n')
    TableFile.write('           \caption[]{' + table_caption + '} \\\\ \hline')
    TableFile.write('           \\rowcolor{headergrey} \hdb{Maneuver Event} & \hdb{Epoch} & \hdb{Magnitude} & \hdb{$\\Delta \\text{v}_x$} & \hdb{$\\Delta \\text{v}_y$} & \hdb{$\\Delta \\text{v}_z$} & \hdb{SEP}       & \hdb{SPE}       \\\\\n')
    TableFile.write('           \\rowcolor{headergrey}                      & \hdb{(TDB)} & \hdb{[m/s]}     & \hdb{[m/s]}                 & \hdb{[m/s]}                 & \hdb{[m/s]}                 & \hdb{[degrees]} & \hdb{[degrees]} \\\\')
    TableFile.write('           \hline\n')
    TableFile.write('           \endhead\n')
    TableFile.write('           \label{tab:'+table_label+'}\n')

    maneuver_counter = 1
    for event in maneuver_event_data:

        date_parts = event[1].split(' ')
        date_string = '\makecell{'+date_parts[0]+' '+date_parts[1]+' '+date_parts[2]+'\\\\'+date_parts[4]+'}'

        TableFile.write('               ' + 'DSM' + str(maneuver_counter) + ' & '
                                          + date_string + ' & '                                            
                                          + str(np.round(event[0].DVmagorThrottle * 1000.0,6)) + ' & '  
                                          + str(np.round(event[0].DeltaVorThrustVectorControl[0] * 1000.0,6)) + ' & ' 
                                          + str(np.round(event[0].DeltaVorThrustVectorControl[1] * 1000.0,6)) + ' & ' 
                                          + str(np.round(event[0].DeltaVorThrustVectorControl[2] * 1000.0,6)) + ' & ' 
                                          + str(np.round(event[2],6)) + ' & '
                                          + str(np.round(event[3],6)) + '\\\\ \hline\n')
        maneuver_counter += 1

    TableFile.write('     \end{longtable}\n')

def generateManeuverPositionTable(table_file_name, table_caption, table_label, maneuver_event_data):

    TableFile = open(table_file_name, mode = 'w')

    TableFile.write('       \\begin{longtable}{|c|c|c|c|c|}\n')
    TableFile.write('           \caption{' + table_caption + '} \\\\')
    TableFile.write('\hline \n')
    TableFile.write('           \\rowcolor{headergrey} \hdb{Maneuver Event} & \hdb{Epoch} & \hdb{x}    & \hdb{y}    & \hdb{z}    \\\\\n')
    TableFile.write('           \\rowcolor{headergrey}                      & \hdb{(TDB)} & \hdb{[km]} & \hdb{[km]} & \hdb{[km]} \\\\')
    TableFile.write('\n')                                                                                                     
    TableFile.write('           \endfirsthead\n')                                                                             
    TableFile.write('           \caption[]{' + table_caption + '} \\\\ \hline')                                                      
    TableFile.write('           \\rowcolor{headergrey} \hdb{Maneuver Event} & \hdb{Epoch} & \hdb{x}    & \hdb{y}    & \hdb{z}    \\\\\n')
    TableFile.write('           \\rowcolor{headergrey}                      & \hdb{(TDB)} & \hdb{[km]} & \hdb{[km]} & \hdb{[km]} \\\\')
    TableFile.write('           \hline\n')
    TableFile.write('           \endhead\n')
    TableFile.write('           \label{tab:'+table_label+'}\n')

    for event in maneuver_event_data:

        date_parts = event[1].split(' ')
        date_string = '\makecell{'+date_parts[0]+' '+date_parts[1]+' '+date_parts[2]+'\\\\'+date_parts[4]+'}'

        TableFile.write('               ' + event[0].parent.central_body + ' & '
                                          + date_string + ' & '                                            
                                          + str(np.round(event[0].SpacecraftState[0],6)) + ' & '  
                                          + str(np.round(event[0].SpacecraftState[1],6)) + ' & ' 
                                          + str(np.round(event[0].SpacecraftState[2],6)) + '\\\\ \hline\n')

    TableFile.write('     \end{longtable}\n')

def generateManeuverVelocityTable(table_file_name, table_caption, table_label, maneuver_event_data):

    TableFile = open(table_file_name, mode = 'w')

    TableFile.write('       \\begin{longtable}{|c|c|c|c|c|c|c|c|}\n')
    TableFile.write('           \caption{' + table_caption + '} \\\\')
    TableFile.write('\hline \n')
    TableFile.write('           \\rowcolor{headergrey}                      &             &              \multicolumn{3}{c}{Pre-maneuver velocity}            &              \multicolumn{3}{c}{Post-maneuver velocity}  \\\\\n')
    TableFile.write('           \\rowcolor{headergrey} \hdb{Maneuver Event} & \hdb{Epoch} & \hdb{$\\text{v}_x$} & \hdb{$\\text{v}_y$}   & \hdb{$\\text{v}_z$} & \hdb{$\\text{v}_x$} & \hdb{$\\text{v}_y$}   & \hdb{$\\text{v}_z$} \\\\\n')
    TableFile.write('           \\rowcolor{headergrey}                      & \hdb{(TDB)} & \hdb{[km]}          & \hdb{[km]}            & \hdb{[km]}          & \hdb{[km/s]}        & \hdb{[km/s]}          & \hdb{[km/s]}        \\\\')
    TableFile.write('\n')
    TableFile.write('           \endfirsthead\n')
    TableFile.write('           \caption[]{' + table_caption + '} \\\\ \hline')
    TableFile.write('           \\rowcolor{headergrey}                      &             &              \multicolumn{3}{c}{Pre-maneuver velocity}            &              \multicolumn{3}{c}{Post-maneuver velocity}  \\\\\n')
    TableFile.write('           \\rowcolor{headergrey} \hdb{Maneuver Event} & \hdb{Epoch} & \hdb{$\\text{v}_x$} & \hdb{$\\text{v}_y$}   & \hdb{$\\text{v}_z$} & \hdb{$\\text{v}_x$} & \hdb{$\\text{v}_y$}   & \hdb{$\\text{v}_z$} \\\\\n')
    TableFile.write('           \\rowcolor{headergrey}                      & \hdb{(TDB)} & \hdb{[km]}          & \hdb{[km]}            & \hdb{[km]}          & \hdb{[km/s]}        & \hdb{[km/s]}          & \hdb{[km/s]}        \\\\')
    TableFile.write('           \hline\n')
    TableFile.write('           \endhead\n')
    TableFile.write('           \label{tab:'+table_label+'}\n')

    for event in maneuver_event_data:

        date_parts = event[1].split(' ')
        date_string = '\makecell{'+date_parts[0]+' '+date_parts[1]+' '+date_parts[2]+'\\\\'+date_parts[4]+'}'

        TableFile.write('               ' + event[0].parent.central_body + ' & '
                                          + date_string + ' & '                                            
                                          + str(np.round(event[0].SpacecraftState[3] - event[0].DeltaVorThrustVectorControl[0],6)) + ' & '  
                                          + str(np.round(event[0].SpacecraftState[4] - event[0].DeltaVorThrustVectorControl[1],6)) + ' & ' 
                                          + str(np.round(event[0].SpacecraftState[5] - event[0].DeltaVorThrustVectorControl[2],6)) + ' & '
                                          + str(np.round(event[0].SpacecraftState[3],6)) + ' & '  
                                          + str(np.round(event[0].SpacecraftState[4],6)) + ' & ' 
                                          + str(np.round(event[0].SpacecraftState[5],6)) + '\\\\ \hline\n')

    TableFile.write('     \end{longtable}\n')

def frameNameTranslator(frame_string, body_name):
    
    state_frame_string = 'ERROR'
    body_name = body_name.upper()
    if frame_string == 'ICRF':
        state_frame_string = body_name + '\_ICRF'
    elif frame_string == 'J2000_BCI':
        state_frame_string = body_name + 'IAU'
    elif frame_string == 'J2000_BCF':
        state_frame_string = 'IAU_Fixed\_' + body_name
    elif frame_string == 'TrueOfDate_BCI':
        state_frame_string = body_name + 'TOD'
    elif frame_string == 'TrueOfDate_BCF':
        state_frame_string = 'TOD' + '\_Fixed\_' + body_name
    elif frame_string == 'Topocentric' or frame_string == 'PrincipleAxes' or frame_string == 'Polar':
        state_frame_string = frame_string
            
    return state_frame_string
    

def filterFormatBplaneParams(periapse_event):

    if periapse_event[0].BdotR != '-':
        BdotR = str(np.round(periapse_event[0].BdotR,6))
    else:
        BdotR = periapse_event[0].BdotR

    if periapse_event[0].BdotT != '-':
        BdotT = str(np.round(periapse_event[0].BdotT,6))
    else:
        BdotT = periapse_event[0].BdotT

    if BdotR == '-' or BdotT == '-':
        B_angle = '-'
    else:
        B_angle = str(np.round(np.arctan2(periapse_event[0].BdotR, periapse_event[0].BdotT) * 180.0 / np.pi,6))

    return BdotR, BdotT, B_angle

def extractPeriapseEvents(EMTG_mission):

    periapse_events = []
    for journey in EMTG_mission.Journeys:
        for event in journey.missionevents:
            if event.EventType == 'periapse':
                periapse_events.append([event, spice.timout(spice.str2et(str(event.JulianDate) + " JD TDB"), "YYYY Mon DD ::TDB HR:MN:SC.###")])

    return periapse_events

def extractManeuverEvents(EMTG_mission):

    maneuver_events = []
    for journey in EMTG_mission.Journeys:
        for event in journey.missionevents:
            if event.EventType == 'chem_burn':
                maneuver_events.append([event, spice.timout(spice.str2et(str(event.JulianDate) + " JD TDB"), "YYYY Mon DD ::TDB HR:MN:SC.###")])

    return maneuver_events

def extractLaunchEvent(EMTG_mission):
    
    for journey in EMTG_mission.Journeys:
        for event in journey.missionevents:
            if event.EventType == 'launch':
                launch_event_data = [event, spice.timout(spice.str2et(str(event.JulianDate) + " JD TDB"), "YYYY Mon DD ::TDB HR:MN:SC.###")]

    return launch_event_data

def computeC3RLADLA(mission_event):

    mu = mission_event[0].parent.mu
    elements = spice.oscelt(mission_event[0].SpacecraftState, spice.str2et(str(mission_event[0].JulianDate) + " JD TDB"), mu)
    SMA = elements[0] / (1.0 - elements[1])
    C3 = -mission_event[0].parent.mu / SMA

    mission_event.append(C3)

    rvec = mission_event[0].SpacecraftState[0:3]
    vvec = mission_event[0].SpacecraftState[3:6]
    r = np.linalg.norm(rvec)
    v = np.linalg.norm(vvec)
    ECCvec = np.multiply((np.multiply((v*v - mu / r), rvec) - np.multiply(np.dot(rvec, vvec), vvec)), 1.0 / mu)
    hvec = np.cross(rvec, vvec)
    h = np.linalg.norm(hvec)

    shat = (1.0 / (1.0 + C3 * (h / mu)**2)) * ((np.sqrt(C3) / mu) * np.cross(hvec, ECCvec) - ECCvec)
    
    RLA = np.arctan2(shat[1], shat[0]) * 180.0 / np.pi
    if RLA < 0.0:
        RLA += 360.0
    DLA = np.arcsin(shat[2]) * 180.0 / np.pi

    mission_event.append(RLA)
    mission_event.append(DLA)

    return mission_event

def computeC3(mission_events):
    
    for event in mission_events:
        mu = event[0].parent.mu
        elements = spice.oscelt(event[0].SpacecraftState, spice.str2et(str(event[0].JulianDate) + " JD TDB"), mu)
        SMA = elements[0] / (1.0 - elements[1])
        C3 = -event[0].parent.mu / SMA
        event.append(C3)

    return mission_events

def computeEventSEPandSPEangles(mission_events):

    mission_event_data = []
    for event in mission_events:
        temp_event = copy.deepcopy(event)
        # get the central body SPICE ID
        central_body_SPICE_ID = spice.bodn2c(temp_event[0].parent.central_body.lower())
        if central_body_SPICE_ID == 399:
            spacecraft_state_wrt_Earth = temp_event[0].SpacecraftState
        else:
            # get position of central body w.r.t. Earth
            state, light_times = spice.spkez(central_body_SPICE_ID, spice.str2et(str(temp_event[0].JulianDate) + " JD TDB"), "J2000", 'NONE', 399)
            central_body_state_wrt_Earth = state
            spacecraft_state_wrt_Earth = temp_event[0].SpacecraftState + central_body_state_wrt_Earth

        Earth_state_wrt_spacecraft = [ -x for x in spacecraft_state_wrt_Earth]
        sun_state_wrt_Earth, light_times = spice.spkez(10, spice.str2et(str(temp_event[0].JulianDate) + " JD TDB"), "J2000", 'NONE', 399)
        sun_state_wrt_spacecraft = sun_state_wrt_Earth - spacecraft_state_wrt_Earth

        cosAngle_SEP = np.dot(sun_state_wrt_Earth[0:3], spacecraft_state_wrt_Earth[0:3])/np.linalg.norm(sun_state_wrt_Earth[0:3])/np.linalg.norm(spacecraft_state_wrt_Earth[0:3])
        cosAngle_SPE = np.dot(sun_state_wrt_spacecraft[0:3], Earth_state_wrt_spacecraft[0:3])/np.linalg.norm(sun_state_wrt_spacecraft[0:3])/np.linalg.norm(Earth_state_wrt_spacecraft[0:3])
        temp_event.extend((np.arccos(cosAngle_SEP) * 180.0/np.pi, np.arccos(cosAngle_SPE) * 180.0/np.pi))
        mission_event_data.append(temp_event)

    return mission_event_data

def filterManeuvers(open_maneuver_event_data, middle_maneuver_event_data, close_maneuver_event_data):
    
    open_maneuvers = []
    middle_maneuvers = []
    close_maneuvers = []
    # this assumes that Open/Middle/Close have the same number of maneuvers in the EMTG problem structure
    for index in range(0, len(open_maneuver_event_data)):
        # filter out maneuvers that are small
        if open_maneuver_event_data[index][0].DVmagorThrottle > 0.00001:
            open_maneuvers.append(open_maneuver_event_data[index])
        else:
            open_maneuvers.append(['-', '-', '-', '-'])

        if middle_maneuver_event_data[index][0].DVmagorThrottle > 0.00001:
            middle_maneuvers.append(middle_maneuver_event_data[index])
        else:
            middle_maneuvers.append(['-', '-', '-', '-'])

        if close_maneuver_event_data[index][0].DVmagorThrottle > 0.00001:
            close_maneuvers.append(close_maneuver_event_data[index])
        else:
            close_maneuvers.append(['-', '-', '-', '-'])

        if open_maneuvers[-1] == ['-', '-', '-', '-'] and middle_maneuvers[-1] == ['-', '-', '-', '-'] and close_maneuvers[-1] == ['-', '-', '-', '-']:
            open_maneuvers.pop()
            middle_maneuvers.pop()
            close_maneuvers.pop()

    return open_maneuvers, middle_maneuvers, close_maneuvers
        

def generateMissionData(EMTG_mission):

    launch_event_data = extractLaunchEvent(EMTG_mission)
    periapse_event_data = extractPeriapseEvents(EMTG_mission)
    maneuver_event_data = extractManeuverEvents(EMTG_mission)

    # filter out maneuvers that are small
    #maneuver_event_data = [(i) for i in maneuver_event_data if i[0].DVmagorThrottle > 0.00001]

    periapse_event_data = computeEventSEPandSPEangles(periapse_event_data)
    maneuver_event_data = computeEventSEPandSPEangles(maneuver_event_data)

    periapse_event_data = computeC3(periapse_event_data)
    launch_event_data = computeC3RLADLA(launch_event_data)
    
    for index in range(0, len(periapse_event_data)):
        if periapse_event_data[index][4] <= 0.0:
            periapse_event_data[index][0].BdotR = '-'
            periapse_event_data[index][0].BdotT = '-'

    # launch_event_data is like    [MissionEvent, TDBG, SEP_anlge, SPE_angle, C3, RLA, DLA]
    # periapse_event_data is like [MissionEvent, TDBG, SEP_angle, SPE_angle, C3]
    # maneuver_event_data is like [MissionEvent, TDBG, SEP_angle, SPE_angle]
    return launch_event_data, periapse_event_data, maneuver_event_data

def replaceEMTGwSPICE(mission_data, spacecraft_SPICE_ID, bsp_file):

    # furnish the spacecraft BSP
    spice.furnsh(bsp_file)
    for journey_index in range(0, len(mission_data.Journeys)):
        event_central_body = mission_data.Journeys[journey_index].central_body 
        state_frame = mission_data.Journeys[journey_index].state_frame
        for event_index in range(0, len(mission_data.Journeys[journey_index].missionevents)):

            # get the SPICE ID of the event central body
            event_central_body_SPICE_ID = spice.bodn2c(event_central_body.lower())

            try:
                state, light_times = spice.spkez(spacecraft_SPICE_ID, spice.str2et(str(mission_data.Journeys[journey_index].missionevents[event_index].JulianDate) + " JD TDB"), "J2000", 'NONE', event_central_body_SPICE_ID)
            except:
                print("spacecraft state not available in " + bsp_file + " at epoch " + SpiceyUtil.julian2Greg(mission_data.Journeys[journey_index].missionevents[event_index].JulianDate))
                continue

            mission_data.Journeys[journey_index].missionevents[event_index].SpacecraftState = state

    spice.unload(bsp_file)