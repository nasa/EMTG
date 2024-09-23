#EMTG Cmake pre-build configuration file.
#This file contains some pre-defined addresses that you may need to customize for your own system.  Setting the 
#below correctly can assist with building EMTG quickly and easily.

#-----CSPICE Directory-----

#This is the location of your CSPICE file.  It must be uncommented, and you must change it to where you put cspice.  If you follow the readme instructions and unzip cspice into the EMTG/cspice folder, the below can remain unchanged.

set(CSPICE_DIR /archive/Utilities/CSPICE)
#set(CSPICE_DIR /archive/Utilities/CSPICE)


#---------SNOPT hints------------
#Change this next line to point towards your snopt7 root directory.  This is not your cppsrc or your lib directory, but up one level from that.
#If SNOPT has been installed on the system path, this hint is probably unnecessary. 
#You can alternatively add and set the SNOPT_INCLUDEDIR and SNOPT_LIBDIR as direct paths to your cppsrc and appropriate lib directory.

#Update this line and change it to be your appropriate path (works for Windows or Unix-based systems)
	set(SNOPT_ROOT_DIR /archive/Libraries/snopt76)
    
#If you are on Windows and have installed Visual Fortran, you will need to point to your Visual Fortran library directory
#    set(INTEL_FORTRAN_DIR "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries_2016.1.146/windows/compiler/lib/intel64")

#-------BOOST HINTS------------------
#cmake usually has an easy time finding a properly installed boost.  If it cannot find your copy of boost, uncomment the next two lines and 
#appropriately modify them to point at your boost distribution.

	set(BOOST_ROOT /archive/Utilities/boost_1_69_0) #set a hint to find boost, but it'll search in the obvious spots too
	set(BOOST_INCLUDE_DIR ${BOOST_ROOT}/boost) #frequently this is where the headers are
	set(BOOST_LIBRARY_DIRS ${BOOST_ROOT}/stage/lib)   # this is sometimes {BOOST_ROOT}/stage/lib	
	
	#This option will tell boost NOT to look at the system path.  If you want to force boost to use some local version, this is how
	#set(Boost_NO_SYSTEM_PATHS ON)
	set (Boost_NO_BOOST_CMAKE ON)
	#set (Boost_USE_MULTITHREADED OFF)
    
#----------GSL path
#If you are using SplineEphem, which requires GSL, supply your path to GSL
set(GSL_PATH /archive/Utilities/gsl-2.2.1/build)

#----------GSAD path
set(GSAD_PATH F:/research/GSAD)

#----------Python versions
#If you plan to use boost::python, you need to specify which version of Python you plan to use. These need to be the same as whatever you linked
#your boost::python against
set(PYTHON_VERSIONS_REQUIRED 3.7)


#