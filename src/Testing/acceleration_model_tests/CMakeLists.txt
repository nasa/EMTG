set(acceleration_model_tests_executable_files   ${acceleration_model_tests_executable_files}
                                                ${CMAKE_CURRENT_SOURCE_DIR}/STM_test.cpp)

source_group("Executable" FILES ${acceleration_model_tests_executable_files})


include_directories(Astrodynamics 
                    Core 
                    Executable
                    HardwareModels
                    InnerLoop
                    Integration
                    Kepler
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
                    
if(SPLINE_EPHEM)
	include_directories(../../SplineEphem)
endif(SPLINE_EPHEM)

link_directories(${PROJECT_BINARY_DIR}/src/Testing/${CMAKE_CFG_INTDIR})

# create the executable
add_executable(acceleration_model_tests  ${acceleration_model_tests_executable_files})
if (USE_AD_INSTRUMENTATION)
    target_link_libraries(acceleration_model_tests testbed_library)
else()
    target_link_libraries(acceleration_model_tests emtg)
endif()