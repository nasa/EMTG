add_subdirectory(FBLT)
add_subdirectory(MGALT)

set(TWO_POINT_SHOOTING_LOW_THRUST_PHASE_HEADERS ${TWO_POINT_SHOOTING_LOW_THRUST_PHASE_HEADERS} 
                                                ${CMAKE_CURRENT_SOURCE_DIR}/TwoPointShootingLowThrustPhase.h
												${CMAKE_CURRENT_SOURCE_DIR}/PhaseDistanceConstraintFactory.h)
																	  
set(TWO_POINT_SHOOTING_LOW_THRUST_PHASE_SOURCE ${TWO_POINT_SHOOTING_LOW_THRUST_PHASE_SOURCE} 
											   ${CMAKE_CURRENT_SOURCE_DIR}/TwoPointShootingLowThrustPhase.cpp
                                               ${CMAKE_CURRENT_SOURCE_DIR}/PhaseDistanceConstraintFactory.cpp)

set(TWO_POINT_SHOOTING_LOW_THRUST_PHASE_FILES ${TWO_POINT_SHOOTING_LOW_THRUST_PHASE_FILES} 
                                              ${TWO_POINT_SHOOTING_LOW_THRUST_PHASE_HEADERS} 
                                              ${TWO_POINT_SHOOTING_LOW_THRUST_PHASE_SOURCE} 
                                              ${FBLT_FILES}
                                              ${MGALT_FILES} PARENT_SCOPE)
