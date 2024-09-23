add_subdirectory(MGAnDSMs)
add_subdirectory(TwoPointShootingLowThrustPhase)
add_subdirectory(CoastPhase)

set(TWO_POINT_SHOOTING_PHASE_HEADERS ${TWO_POINT_SHOOTING_PHASE_HEADERS} 
                                     ${CMAKE_CURRENT_SOURCE_DIR}/TwoPointShootingPhase.h)
																	  
set(TWO_POINT_SHOOTING_PHASE_SOURCE ${TWO_POINT_SHOOTING_PHASE_SOURCE} 
                                    ${CMAKE_CURRENT_SOURCE_DIR}/TwoPointShootingPhase.cpp)

set(TWO_POINT_SHOOTING_PHASE_FILES ${TWO_POINT_SHOOTING_PHASE_FILES} 
                                   ${TWO_POINT_SHOOTING_PHASE_HEADERS} 
                                   ${TWO_POINT_SHOOTING_PHASE_SOURCE} 
                                   ${MGANDSMS_FILES}
                                   ${COASTPHASE_FILES}
                                   ${TWO_POINT_SHOOTING_LOW_THRUST_PHASE_FILES} PARENT_SCOPE)
