add_subdirectory(Exterior)
add_subdirectory(Interior)

set(EPHEMERIS_REFERENCED_ARRIVAL_HEADERS ${EPHEMERIS_REFERENCED_ARRIVAL_HEADERS} 
                               ${CMAKE_CURRENT_SOURCE_DIR}/EphemerisReferencedArrival.h)

set(EPHEMERIS_REFERENCED_ARRIVAL_SOURCE ${FREE_POINT_ARRIVAL_SOURCE} 
                                        ${CMAKE_CURRENT_SOURCE_DIR}/EphemerisReferencedArrival.cpp)

set(EPHEMERIS_REFERENCED_ARRIVAL_FILES ${EPHEMERIS_REFERENCED_ARRIVAL_FILES} 
                                       ${EPHEMERIS_REFERENCED_ARRIVAL_HEADERS} 
                                       ${EPHEMERIS_REFERENCED_ARRIVAL_SOURCE}
                                       ${EPHEMERIS_REFERENCED_ARRIVAL_EXTERIOR_FILES}
                                       ${EPHEMERIS_REFERENCED_ARRIVAL_INTERIOR_FILES}
                                       PARENT_SCOPE)
