"""
Util.py
========

File contains standalone PEATSA utility functions.

"""

# Helper function to convert a Year,Month,day,time of day list into a julian date
def date_list_to_JD(date_list):
    from numpy import sign, floor
    
    # Extract the year
    K = date_list[0]
    # Extract the month
    M = date_list[1]
    # Extract the day
    I = date_list[2]
    # Extract the time of day, in hours
    UT = date_list[3]
    
    # Calculate the Julian Date. 
    # From: http://aa.usno.navy.mil/faq/docs/JD_Formula.php
    return 367*K - floor( ( 7 * ( K + floor((M+9.0)/12.0) ) ) / 4.0 ) + floor( (275*M) / 9.0 ) + I + 1721013.5 + UT/24 - 0.5*sign(100 * K + M - 190002.5) + 0.5
    
# Helper function to convert a month string to an integer month
def convert_month_string_to_int(input):
    
    # I dont think its necessary to comment this function. I think a fifth grader can understand this
    if input == "jan" or input == "Jan" or input == "JAN" or input == "JANUARY" or input == "january" or input == "January" or input == "1" or input == "01":
        return 1
    if input == "feb" or input == "Feb" or input == "FEB" or input == "FEBRUARY" or input == "february" or input == "February" or input == "2" or input == "02":
        return 2
    if input == "mar" or input == "Mar" or input == "MAR" or input == "MARCH" or input == "march" or input == "March" or input == "3" or input == "03":
        return 3
    if input == "apr" or input == "Apr" or input == "APR" or input == "APRIL" or input == "april" or input == "April" or input == "4" or input == "04":
        return 4
    if input == "may" or input == "May" or input == "MAY" or input == "5" or input == "05":
        return 5
    if input == "jun" or input == "Jun" or input == "JUN" or input == "JUNE" or input == "june" or input == "June" or input == "6" or input == "06":
        return 6
    if input == "jul" or input == "Jul" or input == "JUL" or input == "JULY" or input == "july" or input == "July" or input == "7" or input == "07":
        return 7
    if input == "aug" or input == "Aug" or input == "AUG" or input == "AUGUST" or input == "august" or input == "August" or input == "8" or input == "08":
        return 8
    if input == "sep" or input == "Sep" or input == "SEP" or input == "SEPTEMBER" or input == "september" or input == "September" or input == "9" or input == "09":
        return 9
    if input == "oct" or input == "Oct" or input == "OCT" or input == "OCTOBER" or input == "october" or input == "October" or input == "10" or input == "10":
        return 10
    if input == "nov" or input == "Nov" or input == "NOV" or input == "NOVEMBER" or input == "november" or input == "November" or input == "11" or input == "11":
        return 11
    if input == "dec" or input == "Dec" or input == "DEC" or input == "DECEMBER" or input == "december" or input == "December" or input == "12" or input == "12":
        return 12

def convert_time_string_to_JD(input):
    
    return date_list_to_JD(convert_time_string_to_date_list(input))
        
def convert_time_string_to_date_list(input):
    
    date_split1 = input.lstrip(" ").rstrip(" ").replace("  "," ").split(" ")
    date_split2 = input.lstrip(" ").rstrip(" ").replace("  "," ").split("-")
    
    if len(date_split1) == 4:
        date_split = date_split1
    elif len(date_split2) == 4:
        date_split = date_split2
    elif len(date_split1) == 3:
        date_split = date_split1
    elif len(date_split2) == 3:
        last_split = date_split2[-1].split(" ")
        if len(last_split) > 1:
            date_split2[-1] = last_split[0]
            for entry in last_split[1:]:
                date_split2.append(entry)
        date_split = date_split2
                
    else:
        raise Exception("Unsure of time string format: " + input)
    
    if len(date_split) == 4:
        hours = convert_time_string_to_hours_float(date_split[3])
    else:
        hours = 12.0
    
    if int(date_split[0]) < 32 and int(date_split[2]) > 32:
        days = int(date_split[0])
        years = int(date_split[2])
    elif int(date_split[0]) > 32 and int(date_split[2]) < 32: 
        days = int(date_split[2])
        years = int(date_split[0])
        
    return [years, # year
            convert_month_string_to_int(date_split[1]), # month
            days, # days
            hours]
    
# A method to convert a hour:minute:second.fractional_second string into a float of hours past midnight
def convert_time_string_to_hours_float(input):
    # Split the line by colon, assuming it is hours:minutes:seconds.seconds_fraction
    input_split = input.split(":")
    
    # Calculate the fractional hours past midnight
    return float(input_split[0]) + float(input_split[1]) / 60.0 + float(input_split[2]) / 3600.0
            
# Helper function to convert a julian date (or modified julian date) to a Year Mon Day Hour:minute:second string
def jd2datestr(dateval,ifMJD = 0):
    import math
    
    JD = dateval
    if ifMJD:
        JD += 2400000.5
    Q = JD+0.5
    Z = math.floor(Q)
    W = math.floor((Z - 1867216.25)/36524.25)
    X = math.floor(W/4.0)
    A = Z+1+W-X
    B = A+1524
    C = math.floor((B-122.1)/365.25)
    D = math.floor(365.25*C)
    E = math.floor((B-D)/30.6001)
    F = math.floor(30.6001*E)
    day  = B-D-F+(Q-Z)
    month = E-1
    if month >= 13:
        month -=12
    if month == 1 or month == 2:
        year = C-4715
    else:
        year = C-4716
    
    month_list = ["JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC"]
    
    from math import floor
    
    hr = (Q - floor(Q)) * 24.0 
    
    minut = (hr - floor(hr)) * 60.0
    
    sec = (minut - floor(minut)) * 60.0
    
    return str(int(floor(day))) + "-" + month_list[int(floor(month))-1] + "-" + str(int(year)) + " " + str(int(floor(hr))) + ":" + str(int(floor(minut))) + ":" + str(sec)[:10]

def cart2kep(r,v,mu):
    from math  import pi, acos
    from numpy import linalg, cross, dot
    
    r_mag = linalg.norm(r)
    v_mag = linalg.norm(v)
    h = cross(r,v)    
    h_mag = linalg.norm(h)    
    h_hat = [x/h_mag for x in h]
    n = cross([0,0,1],h_hat)    
    n_mag = linalg.norm(n)
    crossVH = cross(v,h)
    e_vec = [crossVH[0]/mu-r[0]/r_mag,crossVH[1]/mu-r[1]/r_mag,crossVH[2]/mu-r[2]/r_mag]    
    E = v_mag*v_mag/2 - mu/r_mag
    
    a = -mu/(2*E)  
    e = linalg.norm(e_vec) 
    i = 0;
    OM = 0;
    omega = 0;
    nu = 0;
        
    delta = 1e-7
    if (e > delta and n_mag>delta):
        i = acos(h[2]/h_mag)
        OM = acos(n[0]/n_mag)
        if (n[1] < 0):
            OM = 2*pi - OM
            
        omega = acos(dot(n,e_vec)/(n_mag*e))
        if (e_vec[2] < 0):
            omega = 2*pi - omega
        
        insides = dot(e_vec,r)/(e*r_mag)
        if abs(insides) > 1:
            if abs(insides) - 1 < 1e-6:
                insides = round(insides)
        nu = acos(insides)
        if (dot(r,v) < 0):
            nu = 2*pi - nu
            
    elif (e < delta and n_mag > delta):
        i = acos(h[2]/h_mag)
        omega = 0
    
        OM = acos(n[0]/n_mag)
        if (n[1] < 0):
            OM = 2*pi - OM   
            
        insides = dot(n,r)/(n_mag*r_mag)
        if abs(insides) > 1:
            if abs(insides) - 1 < 1e-6:
                insides = round(insides)
        nu = acos(insides)
        if (r[2] < 0):
            nu = 2*pi - nu

    elif (e < delta and n_mag < delta):
        i = 0
        omega = 0
        OM = 0
        nu = acos(r[2]/r_mag)
        if (r[2] < 0):
            nu = 2*pi - nu
            
    elif (1 > e and e > delta and n_mag < delta):
        omega = acos(e_vec[0]/e)
        if (e_vec[1] < 0):
            omega = 2*pi - omega
            
        i = 0
        OM = 0
    
        insides = dot(e_vec,r)/(e*r_mag)
        if abs(insides) > 1:
            if abs(insides) - 1 < 1e-6:
                insides = round(insides)
        nu = acos(insides)
        if (dot(r,v) < 0):
            nu = 2*pi - nu

    elif (e > 1 and n_mag < delta):
        omega = acos(e_vec[0]/e)
        if (e_vec[1] < 0):    
            omega = 2*pi - omega
        
        i = 0
        OM = 0
        
        insides = dot(e_vec,r)/(e*r_mag)
        if abs(insides) > 1:
            if abs(insides) - 1 < 1e-6:
                insides = round(insides)
        nu = acos(insides)
        if (dot(r,v) < 0):
            nu = 2*pi - nu
    
    oe = [a,e,i,OM,omega,nu]
    
    return oe

def kep2cart(oe, mu):
    #  import statements
    from math  import cos, sin, sqrt
    from numpy import matrix, array, zeros
    
    #  if 'a' is set to zero, then body is at rest in the
    #+ current frame of reference.
    #+ return a zero vector for r and v
    if oe[0] == 0.0: return zeros(3,1), zeros(3,1)
    
    a  = oe[0]
    e  = oe[1]
    i  = oe[2]
    Om = oe[3]
    om = oe[4]
    f  = oe[5]
    
    p  = a*(1 - e*e)
    r  = p/(1 + e*cos(f))
    rv = matrix([r*cos(f), r*sin(f),   0])
    vv = matrix([-sin(f),  e + cos(f), 0])
    vv = sqrt(mu/p)*vv
    
    c0 = cos(Om); s0 = sin(Om)
    co = cos(om); so = sin(om)
    ci = cos(i);  si = sin(i)
    
    R  = matrix([[c0*co - s0*so*ci, -c0*so - s0*co*ci,  s0*si],
                 [s0*co + c0*so*ci, -s0*so + c0*co*ci, -c0*si],
                 [so*si,             co*si,             ci]])
                 
    ri = array(R*rv.T); ri = ri.reshape(3)
    vi = array(R*vv.T); vi = vi.reshape(3)
    
    return ri, vi
    