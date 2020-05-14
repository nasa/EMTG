add_subdirectory(EphemerisPeggedDeparture)
add_subdirectory(EphemerisReferencedDeparture)
add_subdirectory(FreePointDeparture)
add_subdirectory(PeriapseDeparture)

set(DEPARTURE_HEADERS ${DEPARTURE_HEADERS} ${CMAKE_CURRENT_SOURCE_DIR}/departure.h)

set(DEPARTURE_SOURCE ${DEPARTURE_SOURCE} ${CMAKE_CURRENT_SOURCE_DIR}/departure.cpp)

set(DEPARTURE_FILES ${DEPARTURE_FILES} ${DEPARTURE_HEADERS} 
                                       ${DEPARTURE_SOURCE} 
                                       ${EPHEMERIS_PEGGED_DEPARTURE_FILES} 
                                       ${EPHEMERIS_REFERENCED_DEPARTURE_FILES}
                                       ${PERIAPSE_DEPARTURE_FILES}
                                       ${FREE_POINT_DEPARTURE_FILES} PARENT_SCOPE)
