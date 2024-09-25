"""
State.py
============

Contains state class, which holds data related to spacecraft state.
"""

import Util

# Create a class to hold a state object
class state(object):
    """
    Holds data related to spacecraft state. Constructor sets class variables for states equal to indata.
    
    Parameters
    ----------
    indata : 8-element list of real numbers, optional
        The state of the spacecraft. Organized as
            
        [epoch (JD TDB), x (km), y (km), z (km), vx (km/s), vy (km/s), vz (km/s), mass (kg)].
        
        The default is [0,0,0,0,0,0,0,0].

    Returns
    -------
    None.
    """
    # Initialization method
    def __init__(self,indata = [0,0,0,0,0,0,0,0]):
        # Initialize the state variables
        self.epoch = indata[0] # Julian Days
        """Current epoch (JD TDB)"""
        self.x = indata[1] # km
        """Cartesian x state (km)"""
        self.y = indata[2] # km
        """Cartesian y state (km)"""
        self.z = indata[3] # km
        """Cartesian z state (km)"""
        self.vx = indata[4] # km/s
        """Cartesian vx state (km/s)"""
        self.vy = indata[5] # km/s
        """Cartesian vy state (km/s)"""
        self.vz = indata[6] # km/s
        """Cartesian vz state (km/s)"""
        self.mass = indata[7] # kg
        """Mass (kg)"""
        # In certain other cases we might want more details about the state, like the name, frame, central body, etc.
        self.name = ""
        """Name of state (string)"""
        self.frame = ""
        """Frame in which state is specified (string)"""
        self.central_body = ""
        """Central body relative to which the state is specified (string)"""
        self.b_dot_r = 0 # Not currently used for anything, but probably km
        """B plane B dot R parameter (km). NOT CURRENTLY SET OR USED"""
        self.b_dot_t = 0 # Not currently used for anything, but probably km
        """B plane B dot T parameter (km). NOT CURRENTLY SET OR USED"""
    
    def updateStates(self,indata):
        """
        Set the state class variables to be the contents of indata.

        Parameters
        ----------
        indata : 8-element list of real numbers
            The state of the spacecraft. Organized as 
        
            [epoch (JD TDB), x (km), y (km), z (km), vx (km/s), vy (km/s), vz (km/s), mass (kg)]

        Returns
        -------
        None.

        """
        self.epoch = indata[0] # Julian Days
        self.x = indata[1] # km
        self.y = indata[2] # km
        self.z = indata[3] # km
        self.vx = indata[4] # km/s
        self.vy = indata[5] # km/s
        self.vz = indata[6] # km/s
        self.mass = indata[7] # kg
        
    def updateStatesFromState(self,right):
        """
        Set the state class variables to be the contents of right.

        Parameters
        ----------
        right : state object.
            The state of the spacecraft. Organized as
                
            right.epoch (JD TDB)
            
            right.x (km)
            
            right.y (km)
            
            right.z (km)
            
            right.vx (km/s)
            
            right.vz (km/s)
            
            right.mass (kg)

        Returns
        -------
        None.

        """
        self.epoch = right.epoch
        self.x = right.x
        self.y = right.y
        self.z = right.z
        self.vx = right.vx
        self.vy = right.vy
        self.vz = right.vz
        self.mass = right.mass
    
    def getList(self):
        """
        Return the position and velocity state as a list rather than as an object.
        
        Parameters
        ----------
        None.
        
        Returns
        -------
        List : List
            [x (km), y (km), z (km), vx (km/s), vy (km/s), vz (km/s)]
            
        """
        return [self.x,self.y,self.z,self.vx,self.vy,self.vz]
      
    def printState(self):
        """
        Print the 8 state elements to standard output. Time is written as JD and date string.
        
        Parameters
        ----------
        None.
        
        Returns
        -------
        None.
        
        """
        print("t = " + str(self.epoch) + "(" + Util.jd2datestr(self.epoch) +")")
        print("X = " + str(self.x))
        print("Y = " + str(self.y))
        print("Z = " + str(self.z))
        print("VX = " + str(self.vx))
        print("VY = " + str(self.vy))
        print("VZ = " + str(self.vz))
        print("Mass = " + str(self.mass))
    
    def getState(self,idx):
        """
        Return a single state element.

        Parameters
        ----------
        idx : Integer
            The index of the state element to return [0-6]. CANNOT GET EPOCH, so element 0 is x.

        Returns
        -------
        Real
            The single returned state element.

        """
        if idx == 0:
            return self.x 
        elif idx == 1:
            return self.y 
        elif idx == 2:
            return self.z
        elif idx == 3:
            return self.vx 
        elif idx == 4:
            return self.vy
        elif idx == 5:
            return self.vz
        elif idx == 6:
            return self.mass
            
    def getStateName(self,idx):
        """
        Return string name of a single state element.

        Parameters
        ----------
        idx : Integer
            The index of the state whose name is to be returned [0-6]. EXCLUDES EPOCH, so element 0 is x.

        Returns
        -------
        str
            Single string name of the state to be returned. 

        """
        if idx == 0:
            return "x" 
        elif idx == 1:
            return "y"
        elif idx == 2:
            return "z"
        elif idx == 3:
            return "vx"
        elif idx == 4:
            return "vy"
        elif idx == 5:
            return "vz"
        elif idx == 6:
            return "mass"
            
    def parse_ephem_line(self,linestring):
        """
        Initializes the state variables from a line in a .ephemeris file. It is assumed that the line has a certain structure of 
        
        Year Mon Day Hour:minute:second.fractional_seconds, x,y,z,vx,vy,vz,mass

        Parameters
        ----------
        linestring : String
            Line from an ephemeris file organized as
            
            Year Mon Day Hour:minute:second.fractional_seconds, x,y,z,vx,vy,vz,mass

        Returns
        -------
        None.

        """
        linesplit = linestring.split(",")
        
        self.central_body = "Sun"
        self.frame = "EME2000"
        self.epoch = Util.convert_time_string_to_JD(linesplit[0])
        self.x = float(linesplit[1])
        self.y = float(linesplit[2])
        self.z = float(linesplit[3])
        self.vx = float(linesplit[4])
        self.vy = float(linesplit[5])
        self.vz = float(linesplit[6])
        self.mass = float(linesplit[7])
        
    # A function to 
    def parse_target_spec_line(self,linestring):
        """
        Initialize the state data from a .mission_target_spec line. It is assumed that the line has a certain structure of
        
        <EVENTNAME>,<SAMPLE>,<CBODY>,<FRAME>,<EPOCH(ET)>,<X[km]>,<Y[km]>,<Z[km]>,<VX[km/s]>,<VY[km/s]>,<VZ[km/s]>,<MASS[kg]>,<B.R[km]>,<B.T[km]>

        Parameters
        ----------
        linestring : String
            Line from a .mission_target_spec file, organized as
            
            <EVENTNAME>,<SAMPLE>,<CBODY>,<FRAME>,<EPOCH(ET)>,<X[km]>,<Y[km]>,<Z[km]>,<VX[km/s]>,<VY[km/s]>,<VZ[km/s]>,<MASS[kg]>,<B.R[km]>,<B.T[km]>


        Returns
        -------
        None.

        """
        linesplit = linestring.split(",")
    
        self.name = linesplit[0]
        self.central_body = linesplit[1]
        self.frame = linesplit[2]
        self.epoch = Util.convert_time_string_to_JD(linesplit[3])
        self.x = float(linesplit[4])
        self.y = float(linesplit[5])
        self.z = float(linesplit[6])
        self.vx = float(linesplit[7])
        self.vy = float(linesplit[8])
        self.vz = float(linesplit[9])
        self.mass = float(linesplit[10])
        self.b_dot_r = float(linesplit[11])
        self.b_dot_t = float(linesplit[12])