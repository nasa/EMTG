set(PARALLEL_SHOOTING_PHASE_HEADERS ${PARALLEL_SHOOTING_PHASE_HEADERS} 
                                    ${CMAKE_CURRENT_SOURCE_DIR}/ParallelShootingStep_maneuver_constraint.h
                                    ${CMAKE_CURRENT_SOURCE_DIR}/ParallelShootingStep_BPT_angle_constraint.h
                                    ${CMAKE_CURRENT_SOURCE_DIR}/ParallelShootingStepDistanceConstraint.h
                                    ${CMAKE_CURRENT_SOURCE_DIR}/ParallelShootingStep_maneuver_constraint_factory.h)
                     
set(PARALLEL_SHOOTING_PHASE_SOURCE ${PARALLEL_SHOOTING_PHASE_SOURCE} 
                                   ${CMAKE_CURRENT_SOURCE_DIR}/ParallelShootingStep_maneuver_constraint.cpp
                                   ${CMAKE_CURRENT_SOURCE_DIR}/ParallelShootingStep_BPT_angle_constraint.cpp
                                   ${CMAKE_CURRENT_SOURCE_DIR}/ParallelShootingStepDistanceConstraint.cpp
                                   ${CMAKE_CURRENT_SOURCE_DIR}/ParallelShootingStep_maneuver_constraint_factory.cpp)
                     
set(PARALLEL_SHOOTING_PHASE_FILES ${PARALLEL_SHOOTING_PHASE_FILES} 
                                  ${PARALLEL_SHOOTING_PHASE_HEADERS} 
                                  ${PARALLEL_SHOOTING_PHASE_SOURCE} PARENT_SCOPE)
