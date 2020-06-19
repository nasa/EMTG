set(mission_test_executable_files   ${mission_test_executable_files}
                                    ${CMAKE_CURRENT_SOURCE_DIR}/testbed_driver.cpp
                                    ${CMAKE_CURRENT_SOURCE_DIR}/mission_testbed.h
                                    ${CMAKE_CURRENT_SOURCE_DIR}/mission_testbed.cpp
                                    ${CMAKE_CURRENT_SOURCE_DIR}/MTMtestbed.h
                                    ${CMAKE_CURRENT_SOURCE_DIR}/MTMtestbed.cpp
                                    ${CMAKE_CURRENT_SOURCE_DIR}/propagatorTestbed.h
                                    ${CMAKE_CURRENT_SOURCE_DIR}/propagatorTestbed.cpp
                                    ${CMAKE_CURRENT_SOURCE_DIR}/eventTestbed.h
                                    ${CMAKE_CURRENT_SOURCE_DIR}/eventTestbed.cpp
                                    ${CMAKE_CURRENT_SOURCE_DIR}/mathTestbed.h
                                    ${CMAKE_CURRENT_SOURCE_DIR}/mathTestbed.cpp
									${CMAKE_CURRENT_SOURCE_DIR}/MIRAGE_verification.h
                                    ${CMAKE_CURRENT_SOURCE_DIR}/MIRAGE_verification.cpp
                                    ${CMAKE_CURRENT_SOURCE_DIR}/frame_testbed.h
                                    ${CMAKE_CURRENT_SOURCE_DIR}/frame_testbed.cpp
                                    ${CMAKE_CURRENT_SOURCE_DIR}/StateRepresentation_testbed.h
                                    ${CMAKE_CURRENT_SOURCE_DIR}/StateRepresentation_testbed.cpp)

source_group("Executable" FILES ${mission_test_executable_files})


include_directories(Astrodynamics 
                    Core 
                    Executable
                    HardwareModels 
                    Integration 
                    Kepler 
                    Math 
                    Propagation 
                    Utilities)

include_directories(Astrodynamics/AccelerationModel
                    Astrodynamics/EquationsOfMotion
                    Astrodynamics/StateRepresentation)
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
                    Mission/Journey/Phase/BoundaryEvents 
                    Mission/Journey/Phase/BoundaryEvents/Arrival
                    Mission/Journey/Phase/BoundaryEvents/Arrival/EphemerisPeggedArrival
                    Mission/Journey/Phase/BoundaryEvents/Arrival/EphemerisReferencedArrival
                    Mission/Journey/Phase/BoundaryEvents/Arrival/FreePointArrival
                    Mission/Journey/Phase/BoundaryEvents/Departure
                    Mission/Journey/Phase/BoundaryEvents/Departure/EphemerisPeggedDeparture
                    Mission/Journey/Phase/BoundaryEvents/Departure/EphemerisReferencedDeparture
                    Mission/Journey/Phase/BoundaryEvents/Departure/FreePointDeparture
                    Mission/Journey/Phase/BoundaryEvents/SpecializedBoundaryConstraints
                    Mission/ObjectiveFunctions)
                    
if(SPLINE_EPHEM)
	include_directories(../../SplineEphem)
endif(SPLINE_EPHEM)

if(GSAD_INSTRUMENTATION)
	include_directories(${GSAD_PATH})
endif(GSAD_INSTRUMENTATION)

link_directories(${PROJECT_BINARY_DIR}/src/Testing/${CMAKE_CFG_INTDIR})

# create the executable
add_executable(mission_testbed  ${mission_test_executable_files})
target_link_libraries(mission_testbed testbed_library)