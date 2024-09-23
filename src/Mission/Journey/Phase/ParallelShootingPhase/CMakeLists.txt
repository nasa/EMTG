add_subdirectory(PSFBphase)
add_subdirectory(PSBIphase)
add_subdirectory(ParallelShootingConstraints)

set(PARALLEL_SHOOTING_PHASE_HEADERS ${PARALLEL_SHOOTING_PHASE_HEADERS} 
                                    ${CMAKE_CURRENT_SOURCE_DIR}/ParallelShootingPhase.h
                                    ${CMAKE_CURRENT_SOURCE_DIR}/ParallelShootingStep.h
                                    ${CMAKE_CURRENT_SOURCE_DIR}/ParallelShootingFirstStep.h
                                    ${CMAKE_CURRENT_SOURCE_DIR}/ParallelShootingLastStep.h
                                    ${CMAKE_CURRENT_SOURCE_DIR}/ParallelShootingOneStepToRuleThemAll.h)
																	  
set(PARALLEL_SHOOTING_PHASE_SOURCE ${PARALLEL_SHOOTING_PHASE_SOURCE} 
                                   ${CMAKE_CURRENT_SOURCE_DIR}/ParallelShootingPhase.cpp
                                   ${CMAKE_CURRENT_SOURCE_DIR}/ParallelShootingStep.cpp
                                   ${CMAKE_CURRENT_SOURCE_DIR}/ParallelShootingFirstStep.cpp
                                   ${CMAKE_CURRENT_SOURCE_DIR}/ParallelShootingLastStep.cpp
                                   ${CMAKE_CURRENT_SOURCE_DIR}/ParallelShootingOneStepToRuleThemAll.cpp)

set(PARALLEL_SHOOTING_PHASE_FILES ${PARALLEL_SHOOTING_PHASE_FILES} 
                                  ${PARALLEL_SHOOTING_PHASE_HEADERS} 
                                  ${PARALLEL_SHOOTING_PHASE_SOURCE} 
                                  ${PSFB_PHASE_FILES}
                                  ${PSBI_PHASE_FILES} PARENT_SCOPE)
