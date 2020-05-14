add_subdirectory(Exterior)
add_subdirectory(Interior)

set(EPHEMERIS_REFERENCED_DEPARTURE_HEADERS ${EPHEMERIS_REFERENCED_DEPARTURE_HEADERS} 
                               ${CMAKE_CURRENT_SOURCE_DIR}/EphemerisReferencedDeparture.h)

set(EPHEMERIS_REFERENCED_DEPARTURE_SOURCE ${FREE_POINT_DEPARTURE_SOURCE} 
                                        ${CMAKE_CURRENT_SOURCE_DIR}/EphemerisReferencedDeparture.cpp)

set(EPHEMERIS_REFERENCED_DEPARTURE_FILES ${EPHEMERIS_REFERENCED_DEPARTURE_FILES} 
                                       ${EPHEMERIS_REFERENCED_DEPARTURE_HEADERS} 
                                       ${EPHEMERIS_REFERENCED_DEPARTURE_SOURCE}
                                       ${EPHEMERIS_REFERENCED_DEPARTURE_EXTERIOR_FILES}
                                       ${EPHEMERIS_REFERENCED_DEPARTURE_INTERIOR_FILES}
                                       PARENT_SCOPE)
