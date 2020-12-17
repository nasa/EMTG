#function to convert old-style periapse launch (spherical coordinates) to new-style outgoing B-plane
#Jacob Englander 7-27-2020

def convert_launch(myJourneyOptions, mu=1.0):
    import StateConverter

    #step 1: convert all instances of "PeriapseLaunchOrImpulsiveDeparture" to "PeriapseLaunch"
    Xindex = 0
    has_launch = False
    X = myJourneyOptions.trialX

    while Xindex < len(X):
        if 'PeriapseLaunchOrImpulsiveDeparture' in X[Xindex][0]:
            X[Xindex][0] = X[Xindex][0].replace('PeriapseLaunchOrImpulsiveDeparture', 'PeriapseLaunch')
            has_launch = True

        Xindex += 1

    if has_launch:

        #step 2: step through the decision vector, looking for the spherical state entries and renaming as we go
        stateSpherical = [0.0]*6
        stateSphericalType = 'have not decided yet'
        Xindices_6state = []
        Xindex = 0

        while Xindex < len(X):
            description = X[Xindex][0]
            value = X[Xindex][1]
        
            if 'PeriapseLaunch' in description:
                if 'event left state r ' in description:
                    stateSpherical[0] = float(value)
                    Xindices_6state.append(Xindex)
                    X[Xindex][0] = X[Xindex][0].replace('event left state r ', 'event left state VINFout')
                elif 'event left state RA' in description:
                    stateSpherical[1] = float(value)
                    Xindices_6state.append(Xindex)
                    X[Xindex][0] = X[Xindex][0].replace('event left state RA', 'event left state RHAout')
                elif 'event left state DEC' in description:
                    stateSpherical[2] = float(value)
                    Xindices_6state.append(Xindex)
                    X[Xindex][0] = X[Xindex][0].replace('event left state DEC', 'event left state DHAout')
                elif 'event left state v ' in description:
                    stateSpherical[3] = float(value)
                    Xindices_6state.append(Xindex)
                    X[Xindex][0] = X[Xindex][0].replace('event left state v ', 'event left state BRADIUSout')
                elif 'event left state vRA' in description:
                    stateSphericalType = 'SphericalRADEC'
                    stateSpherical[4] = float(value)
                    Xindices_6state.append(Xindex)
                    X[Xindex][0] = X[Xindex][0].replace('event left state vRA', 'event left state BTHETAout')
                elif 'event left state vDEC' in description:
                    stateSpherical[5] = float(value)
                    Xindices_6state.append(Xindex)
                    X[Xindex][0] = X[Xindex][0].replace('event left state vDEC', 'event left state TAout')
                elif 'event left state AZ' in description:
                    stateSphericalType = 'SphericalAZFPA'
                    stateSpherical[4] = float(value)
                    Xindices_6state.append(Xindex)
                    X[Xindex][0] = X[Xindex][0].replace('event left state AZ', 'event left state BTHETAout')
                elif 'event left state FPA' in description:
                    stateSpherical[5] = float(value)
                    Xindices_6state.append(Xindex)
                    X[Xindex][0] = X[Xindex][0].replace('event left state FPA', 'event left state TAout')
                elif 'magnitude of periapse maneuver' in description:
                    stateSpherical[3] += float(value)

            Xindex += 1

        #step 3: convert to outgoingBplane
        stateOutgoingBplane = []
        myStateConverter = StateConverter.StateConverter()

        if stateSphericalType == 'SphericalRADEC':
            stateOutgoingBplane = myStateConverter.SphericalRADECtoOutgoingBplane(stateSpherical, mu)
        else:
            stateOutgoingBplane = myStateConverter.SphericalAZFPAtoOutgoingBplane(stateSpherical, mu)

        #step 4: put the new entries in the decision vector
        for stateIndex in range(0, 6):
            X[Xindices_6state[stateIndex]][1] = stateOutgoingBplane[stateIndex]

        #step 5: convert the upper and lower bounds on the journey initial velocity bounds
        from copy import deepcopy
        oldBounds = deepcopy(myJourneyOptions.initial_impulse_bounds)

        r_planet = 6378.1363
        try:                                                                                                                                                                                                  
            import Universe                                                                                                                                                                                   
            myUniverse = Universe.Universe(myJourneyOptions.universe_folder + "/" + myJourneyOptions.journey_central_body + ".emtg_universe")                                                                                         
            r_planet = myUniverse.central_body_radius                                                                                                                                                                                
        except:                                                                                                                                                                                               
            print("Failed to find " + myJourneyOptions.universe_folder + "/" + myJourneyOptions.journey_central_body + ".emtg_universe" + "  Cannot find appropriate mu for decision vector conversion. Using 1.0. Good luck.")       

        periapse_altitude_min = r_planet + myJourneyOptions.PeriapseDeparture_altitude_bounds[0]
        periapse_altitude_max = r_planet + myJourneyOptions.PeriapseDeparture_altitude_bounds[1]

        parking_velocity_min_altitude = (mu / periapse_altitude_min) ** 0.5
        parking_velocity_max_altitude = (mu / periapse_altitude_max) ** 0.5

        if (oldBounds[0] + parking_velocity_min_altitude)**2 - 2*mu/periapse_altitude_min > 1.0e-8:
            myJourneyOptions.initial_impulse_bounds[0] = ((oldBounds[0] + parking_velocity_min_altitude)**2 - 2*mu/periapse_altitude_min) ** 0.5
        else:
            myJourneyOptions.initial_impulse_bounds[0] = 1.0e-8

            
        if (oldBounds[1] + parking_velocity_max_altitude)**2 - 2*mu/periapse_altitude_max > 1.0e-8:
            myJourneyOptions.initial_impulse_bounds[1] = ((oldBounds[1] + parking_velocity_max_altitude)**2 - 2*mu/periapse_altitude_max) ** 0.5
        else:
            myJourneyOptions.initial_impulse_bounds[1] = 1.0e-8

    return myJourneyOptions