add_subdirectory(Arrival)
add_subdirectory(Departure)
add_subdirectory(SpecializedBoundaryConstraints)

set(BOUNDARY_EVENTS_HEADERS ${BOUNDARY_EVENTS_HEADERS} ${CMAKE_CURRENT_SOURCE_DIR}/BoundaryEventBase.h
                                                       ${CMAKE_CURRENT_SOURCE_DIR}/EphemerisPeggedBoundary.h
                                                       ${CMAKE_CURRENT_SOURCE_DIR}/FreePointBoundary.h
                                                       ${CMAKE_CURRENT_SOURCE_DIR}/EphemerisReferencedBoundary.h
													   ${CMAKE_CURRENT_SOURCE_DIR}/PeriapseBoundary.h
                                                       ${CMAKE_CURRENT_SOURCE_DIR}/EdelbaumSpiral.h
                                                       ${CMAKE_CURRENT_SOURCE_DIR}/EdelbaumSpiralSegment.h)
																	  
set(BOUNDARY_EVENTS_SOURCE ${BOUNDARY_EVENTS_SOURCE} ${CMAKE_CURRENT_SOURCE_DIR}/BoundaryEventBase.cpp
                                                     ${CMAKE_CURRENT_SOURCE_DIR}/EphemerisPeggedBoundary.cpp
                                                     ${CMAKE_CURRENT_SOURCE_DIR}/FreePointBoundary.cpp
                                                     ${CMAKE_CURRENT_SOURCE_DIR}/EphemerisReferencedBoundary.cpp
													 ${CMAKE_CURRENT_SOURCE_DIR}/PeriapseBoundary.cpp
                                                     ${CMAKE_CURRENT_SOURCE_DIR}/EdelbaumSpiral.cpp
                                                     ${CMAKE_CURRENT_SOURCE_DIR}/EdelbaumSpiralSegment.cpp)

set(BOUNDARY_EVENTS_FILES ${BOUNDARY_EVENTS_FILES} ${BOUNDARY_EVENTS_HEADERS} 
                                                   ${BOUNDARY_EVENTS_SOURCE} 
                                                   ${ARRIVAL_FILES} 
                                                   ${DEPARTURE_FILES}
                                                   ${SPECIALIZED_BOUNDARY_CONSTRAINTS_FILES} PARENT_SCOPE)
