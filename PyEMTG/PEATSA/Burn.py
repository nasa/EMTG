"""
Burn.py

Holds burn class, which is a PEATSA class for holding a maneuver.

TODO: rename class to maneuver

"""
import Util

class burn(object):
    """
    PEATSA class for holding a maneuver.
    
    Parameters
    ----------
    None.

    Returns
    -------
    None.
    
    """
    # Initialize the maneuver data
    def __init__(self):
        self.start_epoch = 0 # JD
        """At what epoch does the maneuver begin? (JD TDB)"""
        self.duration = 0.0 # second
        """How long does the maneuver last? (seconds)"""
        self.mass_flow = 0.0 # kg / s
        """What is the engine mass flow during the maneuver? (kg/s)"""
        self.ra = 0.0 # deg
        """What is the right ascension of the burn direction? (deg)"""
        self.dec = 0.0 # deg
        """What is the declination of the burn direction? (deg)"""
        self.thrust = 0.0 # N
        """What is the thrust magnitude of the maneuver? (N)"""
        
        # The following variables are used to help setup optimization problems

        
        self.assumed_start = 0
        """What is the initial guess for start time of this maneuver? (JD TDB) The start_epoch will be optimized relative to this value."""
        self.assumed_duration = 0
        """What is the assumed duration of the maneuver? (seconds) The duration will be optimized relative to this value."""
        self.linked = 0
        """Is this maneuver linked to the previous? (0 for false, 1 for true) As in, does this burn likely need to start within a few hours of the previous burn? If so, this start epoch will be optimized relative to the end of the previous maneuver rather than relative to assumed_start."""
        
        # Everything below here is just used as a reference, and is not used in propagation, or updated after optimization
        
        
        self.name = ""
        """String name of maneuver. Used for reference."""
        self.thrust_dir = [0,0,0]
        """Cartesian components of thrust direction. (3-element list)"""
        self.start_mass = 0
        """Mass of spacecraft at start of maneuver (kg)"""
        self.final_mass = 0
        """Mass of spacecraft at end of maneuver (kg)"""
        self.deltaV = 0.0
        """Total Delta v of maneuver (km/s)"""
        self.frame = ""
        """String name of frame that ra, dec, and thrust_dir are specified in."""
        
    def parse_maneuver_spec_line(self,linestring,which_burn = 0):
        """
        Initializes the maneuver data from a .mission_maneuver_spec file.
        


        Parameters
        ----------
        linestring : String
            A line of a maneuver spec file, organized as
            
            <EVENTNAME>,<NUMBER_OF_MANEUVERS>,<FRAME>,<EPOCH(ET)>,<THRX>,<THRY>,<THRZ>,<THRMAG[N]>,<SMASS[kg]>,<MDOT[kg/s]>,<DUTY>,<FMASS[kg]>,<DV[km/s]>,<DUR[s]>, repeat...  
        which_burn : Integer, options
            Description. The default is 0.

        Returns
        -------
        nBurns : Integer
            Description. Only returned if which_burn == 0.
        """
        linesplit = linestring.split(",")
        
        self.name = linesplit[0]
        
        if which_burn == 0:
            nBurns = int(linesplit[1])
        
        start_idx = 12 * which_burn
        
        self.frame = linesplit[start_idx + 2]
        self.start_epoch = Util.convert_time_string_to_JD( linesplit[start_idx +3] )
        self.assumed_start = self.start_epoch
        
        self.thrust_dir[0] = float(linesplit[start_idx +4])
        self.thrust_dir[1] = float(linesplit[start_idx +5])
        self.thrust_dir[2] = float(linesplit[start_idx +6])
        
        # Make sure that the thrust_dir is a unit vector
        thrust_norm = self.thrust_dir[0]*self.thrust_dir[0]+self.thrust_dir[1]*self.thrust_dir[1]+self.thrust_dir[2]*self.thrust_dir[2]
        if abs(thrust_norm - 1.0) > 1e-3:
            raise Exception("Thrust is not unit")
        
        self.thrust = float(linesplit[start_idx +7])
        self.start_mass = float(linesplit[start_idx +8])
        self.mass_flow = float(linesplit[start_idx +9])
        self.final_mass = float(linesplit[start_idx +11])
        self.deltaV = float(linesplit[start_idx +12])
        self.duration = float(linesplit[start_idx +13]) * float(linesplit[start_idx +10])
        self.assumed_duration = self.duration        
                
        if which_burn == 0:
            return nBurns
            