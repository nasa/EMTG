add_subdirectory(MGALTS)

set(MGALT_HEADERS ${MGALT_HEADERS} ${CMAKE_CURRENT_SOURCE_DIR}/BackwardBoundedImpulseManeuver.h
                                   ${CMAKE_CURRENT_SOURCE_DIR}/BoundedImpulseManeuver.h
                                   ${CMAKE_CURRENT_SOURCE_DIR}/ForwardBoundedImpulseManeuver.h
                                   ${CMAKE_CURRENT_SOURCE_DIR}/MGALTphase.h)
																	  
set(MGALT_SOURCE ${MGALT_SOURCE} ${CMAKE_CURRENT_SOURCE_DIR}/BackwardBoundedImpulseManeuver.cpp
                                 ${CMAKE_CURRENT_SOURCE_DIR}/BoundedImpulseManeuver.cpp
                                 ${CMAKE_CURRENT_SOURCE_DIR}/ForwardBoundedImpulseManeuver.cpp
                                 ${CMAKE_CURRENT_SOURCE_DIR}/MGALTphase.cpp)

set(MGALT_FILES ${MGALT_FILES} ${MGALT_HEADERS} ${MGALT_SOURCE} ${MGALTS_FILES} PARENT_SCOPE)
