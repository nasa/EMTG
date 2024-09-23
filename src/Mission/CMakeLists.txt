add_subdirectory(Journey)
add_subdirectory(ObjectiveFunctions)

set(MISSION_HEADERS  ${MISSION_HEADERS} ${CMAKE_CURRENT_SOURCE_DIR}/mission.h)
																	  
set(MISSION_SOURCE  ${MISSION_SOURCE} ${CMAKE_CURRENT_SOURCE_DIR}/mission.cpp)

set(MISSION_FILES ${MISSION_FILES} ${MISSION_HEADERS} 
                                   ${MISSION_SOURCE} 
                                   ${JOURNEY_FILES} 
                                   ${OBJECTIVE_FUNCTIONS_FILES} PARENT_SCOPE)
