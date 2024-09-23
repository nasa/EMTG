add_subdirectory(Phase)

set(JOURNEY_HEADERS  ${JOURNEY_HEADERS} ${CMAKE_CURRENT_SOURCE_DIR}/journey.h)
																	  
set(JOURNEY_SOURCE  ${JOURNEY_SOURCE} ${CMAKE_CURRENT_SOURCE_DIR}/journey.cpp)

set(JOURNEY_FILES ${MISSION_FILES} ${JOURNEY_HEADERS} 
                                   ${JOURNEY_SOURCE} 
                                   ${PHASE_FILES} PARENT_SCOPE)
