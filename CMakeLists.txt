cmake_minimum_required(VERSION 3.8)

#pull in the user-defined configurations
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR})
include(EMTG-Config.cmake)

#identify the project
project(EMTG)

link_directories(${PROJECT_SOURCE_DIR}/Release)
link_directories(${PROJECT_SOURCE_DIR}/Debug)


# enable multicore build on Windows
if (MSVC)
    set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} " /MP")
    set(CMAKE_C_FLAGS ${CMAKE_C_FLAGS} " /MP")
endif()

# Boost spews auto_ptr deprecation warnings
add_definitions(-D_HAS_AUTO_PTR_ETC)
if (NOT MSVC)
    set(CMAKE_C_FLAGS ${CMAKE_C_FLAGS} " -Wno-deprecated-declarations")
    set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} " -Wno-deprecated-declarations")
endif()

# Boost doesn't recognize newer compilers
add_definitions(-DBOOST_CONFIG_SUPPRESS_OUTDATED_MESSAGE)

# Boost doesn't like std::iterator either
add_definitions(-D_SILENCE_CXX17_ITERATOR_BASE_CLASS_DEPRECATION_WARNING)

# specify where the executable will be installed
if(MSVC)
	STRING(REGEX REPLACE "/" "\\\\" OUTPUTDESTINATION "${PROJECT_SOURCE_DIR}\\bin")
else(MSVC)
	set(OUTPUTDESTINATION ${PROJECT_SOURCE_DIR}/bin)
endif(MSVC)
if(NOT IS_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
	file(MAKE_DIRECTORY ${OUTPUTDESTINATION})
endif(NOT IS_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)


#compiler-time Define-level options  for EMTG
option(QUIET_SNOPT "Makes SNOPT not print things." OFF)
option(SAFE_SNOPT "Trap things that could blow up SNOPT" ON)
option(SNOPT72 "Set SNOPT 7.2" OFF)
option(SNOPT75 "Use SNOPT 7.5" OFF)
option(SNOPT76 "Use SNOPT 7.6" ON)
option(SNOPT_MINGW_DLL "Use SNOPT built with MINGW in Windows" ON)
option(DEBUG_CHECK_VARIABLES_AND_CONSTRAINTS_IN_USER_FUNCTION "Check decision and objective/constraint vectors for NaN and inf" OFF)
option(FAST_EMTG_MATRIX "Enable fast EMTG matrix class (removes bounds checking)" ON)
option(SPLINE_EPHEM "Enable SplineEphem" ON)
option(BUILD_EMTG_TESTBED "Build the EMTG testing framework" OFF)
option(BUILD_PROPULATOR "Build the Propulator" OFF)
option(BUILD_PYHARDWARE "Build PyHardware" OFF)

# We want Background mode to default to ON for NavLab runs
if (MSVC)
    option(BACKGROUND_MODE "This suppresses 'press any key to continue' at the end of an EMTG run" OFF)
else()
    option(BACKGROUND_MODE "This suppresses 'press any key to continue' at the end of an EMTG run" ON)
endif()

#switches for files that open source users don't have access to
option(HAS_PROBEENTRYPHASE "Build with ProbeEntryPhase?" OFF)
option(HAS_BUILT_IN_THRUSTERS "Build with built-in thruster library?" OFF)
mark_as_advanced(FORCE HAS_PROBEENTRYPHASE HAS_BUILT_IN_THRUSTERS BUILD_EMTG_TESTBED BUILD_PROPULATOR BUILD_PYHARDWARE SPLINE_EPHEM DEBUG_CHECK_VARIABLES_AND_CONSTRAINTS_IN_USER_FUNCTION FAST_EMTG_MATRIX)

set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build in Debug or Release mode")
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release)

message(STATUS "Creating build files for " ${CMAKE_BUILD_TYPE} " mode... \n")

if(CMAKE_BUILD_TYPE MATCHES Release)
    if (${CMAKE_CXX_COMPILER_ID} MATCHES "MSVC")
        set(CMAKE_CXX_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE} " /Ot /Ox /GL")
    elseif (${CMAKE_CXX_COMPILER_ID} MATCHES "GNU")
        set(CMAKE_CXX_FLAGS_RELEASE " -g -O2")
    endif()
	message(STATUS "The following compiler optimization flags have been set: " ${CMAKE_CXX_FLAGS_RELEASE} "\n")
elseif(CMAKE_BUILD_TYPE MATCHES Debug)
	message(STATUS "The following compiler optimization flags have been set: " ${CMAKE_CXX_FLAGS_DEBUG} "\n")
endif()

if(QUIET_SNOPT)
	add_definitions(-DQUIET_SNOPT)
endif(QUIET_SNOPT)

if(BACKGROUND_MODE)
    message(STATUS "Background mode activated")
	add_definitions(-DBACKGROUND_MODE)
endif(BACKGROUND_MODE)

if(SAFE_SNOPT)
	add_definitions(-DSAFE_SNOPT)
endif(SAFE_SNOPT)

if(DEBUG_CHECK_VARIABLES_AND_CONSTRAINTS_IN_USER_FUNCTION)
	add_definitions(-DDEBUG_CHECK_VARIABLES_AND_CONSTRAINTS_IN_USER_FUNCTION)
endif(DEBUG_CHECK_VARIABLES_AND_CONSTRAINTS_IN_USER_FUNCTION)

if(SNOPT75)
    add_definitions(-DSNOPT75)
endif(SNOPT75)

if(SNOPT76)
    add_definitions(-DSNOPT76)
endif(SNOPT76)

if(FAST_EMTG_MATRIX)
    add_definitions(-DFAST_EMTG_MATRIX)
endif(FAST_EMTG_MATRIX)


if(HAS_PROBEENTRYPHASE)
    add_definitions(-DHAS_PROBEENTRYPHASE)
endif(HAS_PROBEENTRYPHASE)

if(HAS_BUILT_IN_THRUSTERS)
    add_definitions(-DHAS_BUILT_IN_THRUSTERS)
endif(HAS_BUILT_IN_THRUSTERS)

if(SPLINE_EPHEM)
    add_definitions(-DSPLINE_EPHEM)
	
	# GNU Scientific Library
    set(GSL_INCLUDE_DIRS ${GSL_PATH})    
    include_directories(${GSL_INCLUDE_DIRS})
    
	add_library(gsl STATIC IMPORTED)
	add_library(gslcblas STATIC IMPORTED)
    if(MSVC)
		set_property(TARGET gsl PROPERTY IMPORTED_LOCATION ${GSL_PATH}/${CMAKE_CFG_INTDIR}/gsl.lib)
		set_property(TARGET gslcblas PROPERTY IMPORTED_LOCATION ${GSL_PATH}/${CMAKE_CFG_INTDIR}/gslcblas.lib)        
    else()
        set_property(TARGET gsl PROPERTY IMPORTED_LOCATION ${GSL_PATH}/${CMAKE_CFG_INTDIR}/libgsl.a)
		set_property(TARGET gslcblas PROPERTY IMPORTED_LOCATION ${GSL_PATH}/${CMAKE_CFG_INTDIR}/libgslcblas.a)
    endif()

	set(GSL_LIBRARIES gsl gslcblas)
	set(EMTG_LIBRARIES ${EMTG_LIBRARIES} ${GSL_LIBRARIES})
endif(SPLINE_EPHEM)

#----------------CSpice Work------------------------------------------------------
message(STATUS "Now checking for CSpice...")
if (MSVC) # aka 'windows'
	set(CSPICE_FILE_TO_BE_CORRECTED "${CSPICE_DIR}/src/cspice/mkprodct.bat")
	set(CSPICE_MAKEFILE "makeall.bat")
else (MSVC) #aka Linux or Mac, as both of those use the mkprodct.csh system
	set(CSPICE_FILE_TO_BE_CORRECTED "${CSPICE_DIR}/src/cspice/mkprodct.csh")
	set(CSPICE_MAKEFILE "makeall.csh")
endif(MSVC)


if (EXISTS "${CSPICE_DIR}/${CSPICE_MAKEFILE}")
	message(STATUS "Cspice located and unzipped correctly.")
else()
	message(FATAL_ERROR "CSPICE is not unzipped into the emtg/cspice folder correctly.")
endif()

if (NOT MSVC) #cspice doesn't label the file correctly when it is done; easier to just rename it
	if (NOT EXISTS "${CSPICE_DIR}/lib/libcspice.a")
		message(STATUS "Renaming cspice.a to libcspice.a")
		file(RENAME "${CSPICE_DIR}/lib/cspice.a" "${CSPICE_DIR}/lib/libcspice.a")
		if (NOT EXISTS "${CSPICE_DIR}/lib/libcspice.a")
			message(FATAL_ERROR "CSpice not built properly for a Unix System.")
		endif()
	endif()
else()
	if(NOT EXISTS "${CSPICE_DIR}/lib/cspice.lib")	
		message(FATAL_ERROR "CSpice not built properly for a Windows System.  Please run the cspice/makeall.bat from within the Visual Studio Command Line tools and then reconfigure cmake.")
	endif()	
endif()

set(CSPICE_LIBRARY cspice)
set(EMTG_LIBRARIES ${EMTG_LIBRARIES} ${CSPICE_LIBRARY})
include_directories("${CSPICE_DIR}/include")
link_directories("${CSPICE_DIR}/lib")

#----------------------------SNOPT Stuff-------------------------------------
message(STATUS "Now checking on SNOPT...")

#assume newest version of SNOPT
if (IS_DIRECTORY ${SNOPT_ROOT_DIR}/interfaces)
	include_directories("${SNOPT_ROOT_DIR}/interfaces/include")
	
	set(SNOPT_INTERFACE_HEADERS ${SNOPT_ROOT_DIR}/interfaces/include/snoptProblem.hpp ${SNOPT_ROOT_DIR}/interfaces/include/snopt.h)
	set(SNOPT_INTERFACE_SOURCE ${SNOPT_ROOT_DIR}/interfaces/src/snoptProblem.cpp)
    
    if (MSVC AND SNOPT_MINGW_DLL)
    
        set(SNOPT_LIBRARIES ${SNOPT_ROOT_DIR}/build/${CMAKE_CFG_INTDIR}/libsnopt.lib)
        
        add_library(snopt_interface STATIC IMPORTED)
        set_property(TARGET snopt_interface PROPERTY IMPORTED_LOCATION "${SNOPT_ROOT_DIR}/build/${CMAKE_CFG_INTDIR}/snopt_interface.lib")
        link_directories(${SNOPT_ROOT_DIR}/build/${CMAKE_CFG_INTDIR})
        set(SNOPT_LIBRARIES ${SNOPT_LIBRARIES} snopt_interface)
    
    else()
	
        #PULL IN SNOPT.  THIS ASSUMES SNOPT HAS ALREADY BEEN COMPILED AND IS IN A VERY SPECIFIC SPOT.
        if (NOT MSVC)
            add_library(snopt SHARED IMPORTED)
            add_library(snopt_cpp SHARED IMPORTED)
            set_property(TARGET snopt PROPERTY IMPORTED_LOCATION "${SNOPT_ROOT_DIR}/lib/.libs/libsnopt7.so")
            set_property(TARGET snopt_cpp PROPERTY IMPORTED_LOCATION "${SNOPT_ROOT_DIR}/lib/.libs/libsnopt7_cpp.so")
        else(NOT MSVC)
            add_library(snopt STATIC IMPORTED)
            set_property(TARGET snopt PROPERTY IMPORTED_LOCATION "${SNOPT_ROOT_DIR}/lib/snopt.lib")
            set(EMTG_LIBRARIES ${EMTG_LIBRARIES} ifmodintr.lib)
            link_directories("${INTEL_FORTRAN_DIR}")
            
        endif(NOT MSVC)
        

        if (MSVC)
            link_directories("${PROJECT_SOURCE_DIR}/lib")
            link_directories(${SNOPT_ROOT_DIR}/lib)
        else()
            link_directories("${PROJECT_SOURCE_DIR}/lib")
            link_directories(${SNOPT_ROOT_DIR}/lib/.libs)
        endif()
        add_library(snopt_interface ${SNOPT_INTERFACE_HEADERS} ${SNOPT_INTERFACE_SOURCE})
        if (MSVC)
            set(SNOPT_LIBRARIES snopt_interface snopt)
        else()
            target_link_libraries(snopt_interface gfortran)
            set(SNOPT_LIBRARIES snopt_interface snopt_cpp snopt)
        endif()
        
        set_property(TARGET snopt_interface PROPERTY POSITION_INDEPENDENT_CODE ON)
        
    endif()

	
#okay, do we have an older version of SNOPT?

elseif(IS_DIRECTORY ${SNOPT_ROOT_DIR}/cppsrc)
	include_directories("${SNOPT_ROOT_DIR}/cppsrc")
	include_directories("${SNOPT_ROOT_DIR}/win32")
	set(SNOPT_INTERFACE_HEADERS ${SNOPT_ROOT_DIR}/win32/f2c.h ${SNOPT_ROOT_DIR}/cppsrc/snoptProblem.hh ${SNOPT_ROOT_DIR}/cppsrc/snopt.hh ${SNOPT_ROOT_DIR}/cppsrc/snfilewrapper.hh)
	set(SNOPT_INTERFACE_SOURCE ${SNOPT_ROOT_DIR}/cppsrc/snfilewrapper.c ${SNOPT_ROOT_DIR}/cppsrc/snoptProblem.cc)
	
	#PULL IN SNOPT.  THIS ASSUMES SNOPT HAS ALREADY BEEN COMPILED AND IS IN A VERY SPECIFIC SPOT.
	add_library(f2c STATIC IMPORTED)
	add_library(snopt_cpp STATIC IMPORTED)
	add_library(snopt_c STATIC IMPORTED)
	add_library(snoptsnprint STATIC IMPORTED)
	
	if (NOT MSVC)
		set_property(TARGET f2c PROPERTY IMPORTED_LOCATION "${SNOPT_ROOT_DIR}/lib/libf2c.a")
		set_property(TARGET snopt_cpp PROPERTY IMPORTED_LOCATION "${SNOPT_ROOT_DIR}/lib/libsnopt_cpp.a")
		set_property(TARGET snoptsnprint PROPERTY IMPORTED_LOCATION "${SNOPT_ROOT_DIR}/lib/libsnprint_c.a")
		set_property(TARGET snopt_c PROPERTY IMPORTED_LOCATION "${SNOPT_ROOT_DIR}/lib/libsnopt_c.a")
		set(SNOPT_LIBRARIES snopt_c snoptsnprint snopt_cpp f2c)
	else(NOT MSVC)
		set_property(TARGET f2c PROPERTY IMPORTED_LOCATION "${SNOPT_ROOT_DIR}/win32/snopt_csrc/Release/vcf2c.lib")
		set_property(TARGET snopt_c PROPERTY IMPORTED_LOCATION "${SNOPT_ROOT_DIR}/win32/snopt_csrc/Release/snopt.lib")
		set(SNOPT_LIBRARIES snopt_c f2c)
	endif(NOT MSVC)
	
	link_directories("${PROJECT_SOURCE_DIR}/lib")
	
	add_definitions(-DSNOPT72 "Use SNOPT 7.2-12 or older" ON)

else(IS_DIRECTORY ${SNOPT_ROOT_DIR}/interfaces)
	message(FATAL_ERROR "Could not find snopt")
endif(IS_DIRECTORY ${SNOPT_ROOT_DIR}/interfaces)

include_directories(${SNOPT_INCLUDE_DIR})

set(EMTG_LIBRARIES ${EMTG_LIBRARIES} ${SNOPT_LIBRARIES})

#-------------------------Boost related work------------------------------
message(STATUS "Now looking for BOOST..")
include(FindBoost)
if(MSVC)
	Set(Boost_NAMESPACE libboost)
endif(MSVC)

find_package(Boost 1.60 COMPONENTS filesystem system serialization REQUIRED)

include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})

set(EMTG_LIBRARIES ${EMTG_LIBRARIES} ${Boost_LIBRARIES})

#pull in EMTG itself
add_subdirectory(src)


