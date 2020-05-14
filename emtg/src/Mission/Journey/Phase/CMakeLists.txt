add_subdirectory(BoundaryEvents)
add_subdirectory(TwoPointShootingPhase)
add_subdirectory(ParallelShootingPhase)

set(PHASE_HEADERS ${PHASE_HEADERS} ${CMAKE_CURRENT_SOURCE_DIR}/phase.h)
																	  
set(PHASE_SOURCE ${PHASE_SOURCE} ${CMAKE_CURRENT_SOURCE_DIR}/phase.cpp)

set(PHASE_FILES ${PHASE_FILES} ${PHASE_HEADERS} 
                               ${PHASE_SOURCE} 
                               ${BOUNDARY_EVENTS_FILES}
                               ${TWO_POINT_SHOOTING_PHASE_FILES}
                               ${PARALLEL_SHOOTING_PHASE_FILES} PARENT_SCOPE)
