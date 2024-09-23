cmake_minimum_required(VERSION 3.0)

set (CMAKE_CXX_STANDARD 11)

#name the project, trigger the compiler search

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(EMTG_TESTBED)

option(RUN_ACCELERATION_MODEL_TESTBED "Do you want to run the acceleration model testbed?" OFF)
option(RUN_MISSION_TESTBED "Do you want to run the mission testbed?" ON)
option(RUN_MISSIONOPTIONS_TESTBED "Do you want to run the missionoptions testbed?" ON)
option(USE_AD_INSTRUMENTATION "Do you want to perform derivative accuracy checks with AD?" ON)

if (USE_AD_INSTRUMENTATION AND (RUN_MISSION_TESTBED OR RUN_ACCELERATION_MODEL_TESTBED))
    set(TESTBED_LIBRARY ON)

    add_definitions(-DAD_INSTRUMENTATION)
    # GSAD location
    include_directories(${GSAD_PATH})
    
endif()

if (TESTBED_LIBRARY)
    # build up the source files
    add_subdirectory(../Astrodynamics ${CMAKE_BINARY_DIR}/TESTBED/Astrodynamics)
    add_subdirectory(../Core ${CMAKE_BINARY_DIR}/TESTBED/Core)
    add_subdirectory(../Executable ${CMAKE_BINARY_DIR}/TESTBED/Executable)
    add_subdirectory(../HardwareModels ${CMAKE_BINARY_DIR}/TESTBED/HardwareModels)
    add_subdirectory(../Integration ${CMAKE_BINARY_DIR}/TESTBED/Integration)
    add_subdirectory(../Math ${CMAKE_BINARY_DIR}/TESTBED/Math)
    add_subdirectory(../Mission ${CMAKE_BINARY_DIR}/TESTBED/Mission)
    add_subdirectory(../Propagation ${CMAKE_BINARY_DIR}/TESTBED/Propagation)
    add_subdirectory(../Utilities ${CMAKE_BINARY_DIR}/TESTBED/Utilities)
    add_subdirectory(../Scalatron ${CMAKE_BINARY_DIR}/TESTBED/Scalatron)
    if(SPLINE_EPHEM)
        add_subdirectory(../SplineEphem ${CMAKE_BINARY_DIR}/TESTBED/SplineEphem)
        include_directories(../SplineEphem)
        include_directories(${GSL_INCLUDE_DIRS})
    endif(SPLINE_EPHEM)
    include_directories(Astrodynamics 
                        Core 
                        Executable
                        HardwareModels 
                        InnerLoop 
                        Integration 
                        Math
                        Propagation 
                        Utilities
                        Scalatron)

    include_directories(Astrodynamics/AccelerationModel
                        Astrodynamics/EquationsOfMotion)
    include_directories(Mission 
                        Mission/Journey 
                        Mission/Journey/Phase 
                        Mission/Journey/Phase/TwoPointShootingPhase 
                        Mission/Journey/Phase/TwoPointShootingPhase/TwoPointShootingLowThrustPhase
                        Mission/Journey/Phase/TwoPointShootingPhase/TwoPointShootingLowThrustPhase/FBLT 
                        Mission/Journey/Phase/TwoPointShootingPhase/TwoPointShootingLowThrustPhase/FBLT/FBLTS
                        Mission/Journey/Phase/TwoPointShootingPhase/TwoPointShootingLowThrustPhase/MGALT 
                        Mission/Journey/Phase/TwoPointShootingPhase/TwoPointShootingLowThrustPhase/MGALT/MGALTS
                        Mission/Journey/Phase/TwoPointShootingPhase/MGAnDSMs
                        Mission/Journey/Phase/TwoPointShootingPhase/MGAnDSMs/ManeuverConstraints
                        Mission/Journey/Phase/TwoPointShootingPhase/CoastPhase
                        Mission/Journey/Phase/TwoPointShootingPhase/CoastPhase/SundmanCoastPhase
                        Mission/Journey/Phase/TwoPointShootingPhase/CoastPhase/ControlLawThrustPhase
                        Mission/Journey/Phase/ParallelShootingPhase
                        Mission/Journey/Phase/ParallelShootingPhase/PSFBphase
                        Mission/Journey/Phase/ParallelShootingPhase/PSBIphase
                        Mission/Journey/Phase/ParallelShootingPhase/ParallelShootingConstraints
                        Mission/Journey/Phase/BoundaryEvents 
                        Mission/Journey/Phase/BoundaryEvents/Arrival
                        Mission/Journey/Phase/BoundaryEvents/Arrival/EphemerisPeggedArrival
                        Mission/Journey/Phase/BoundaryEvents/Arrival/EphemerisReferencedArrival
                        Mission/Journey/Phase/BoundaryEvents/Arrival/FreePointArrival
                        Mission/Journey/Phase/BoundaryEvents/Arrival/PeriapseArrival
                        Mission/Journey/Phase/BoundaryEvents/Departure
                        Mission/Journey/Phase/BoundaryEvents/Departure/EphemerisPeggedDeparture
                        Mission/Journey/Phase/BoundaryEvents/Departure/EphemerisReferencedDeparture
                        Mission/Journey/Phase/BoundaryEvents/Departure/FreePointDeparture
                        Mission/Journey/Phase/BoundaryEvents/Departure/PeriapseDeparture
                        Mission/Journey/Phase/BoundaryEvents/SpecializedBoundaryConstraints
                        Mission/Journey/Phase/BoundaryEvents/SpecializedBoundaryConstraints/OrbitElementsConstraints
                        Mission/Journey/Phase/BoundaryEvents/SpecializedBoundaryConstraints/BoundaryAngleConstraints
                        Mission/ObjectiveFunctions)
    
    if(HAS_PROBEENTRYPHASE)
        include_directories(Mission/Journey/Phase/TwoPointShootingPhase/MGAnDSMs/ProbeEntryPhase)
    endif()

    set(TESTBED_FILES ${TESTBED_FILES} ${ASTRODYNAMICS_FILES}
                                       ${CORE_FILES}
                                       ${HARDWARE_MODEL_FILES}
                                       ${INTEGRATION_FILES}
                                       ${KEPLER_FILES}
                                       ${MATH_FILES}
                                       ${MISSION_FILES}
                                       ${PROPAGATION_FILES}
                                       ${UTILITIES_FILES}
                                       ${SCALATRON_FILES})

    if(SPLINE_EPHEM)
        set(TESTBED_FILES ${TESTBED_FILES} ${SPLINE_EPHEM_FILES})
    endif()

    # Make source groups for IDEs
    set(root_EMTG_src_path ${CMAKE_CURRENT_SOURCE_DIR})
    string(REPLACE "/" "\\" windows_src_path_fixed "${root_EMTG_src_path}")
    foreach(source IN LISTS TESTBED_FILES)
        get_filename_component(source_path "${source}" PATH)
        string(REPLACE "/" "\\" source_path_msvc "${source_path}")
        string(REPLACE "${windows_src_path_fixed}\\" "" source_path_msvc_relative ${source_path_msvc})
        source_group("${source_path_msvc_relative}" FILES "${source}")
    endforeach()

    set(TESTBED_DEPENDENCIES ${TESTBED_DEPENDENCIES} 
                             ${Boost_LIBRARIES}
                             ${GSL_LIBRARIES}
                             ${CSPICE_LIBRARY}) #${Boost_LIBRARIES} is set by findBoost

    # create the EMTG library
    
    add_library(testbed_library ${TESTBED_FILES})
    set_property(TARGET testbed_library PROPERTY CXX_STANDARD 17)
    target_link_libraries(testbed_library ${TESTBED_DEPENDENCIES})
    #target_compile_definitions(testbed_library PRIVATE AD_INSTRUMENTATION)
endif()

if (RUN_ACCELERATION_MODEL_TESTBED)
    add_subdirectory(acceleration_model_tests)
endif()

if (RUN_MISSION_TESTBED)
    add_subdirectory(mission_tests)
endif()

if (RUN_MISSIONOPTIONS_TESTBED)
    add_subdirectory(missionoptions_tests)
endif()