
set(EPHEMERIS_REFERENCED_ARRIVAL_INTERIOR_HEADERS ${EPHEMERIS_REFERENCED_ARRIVAL_INTERIOR_HEADERS} 
                               ${CMAKE_CURRENT_SOURCE_DIR}/EphemerisReferencedArrivalInterior.h
                               ${CMAKE_CURRENT_SOURCE_DIR}/EphemerisReferencedLTRendezvousInterior.h
                               ${CMAKE_CURRENT_SOURCE_DIR}/EphemerisReferencedArrivalWithVinfinityInterior.h
                               #${CMAKE_CURRENT_SOURCE_DIR}/EphemerisReferencedChemRendezvousInterior.h
                               ${CMAKE_CURRENT_SOURCE_DIR}/EphemerisReferencedInterceptInterior.h)

set(EPHEMERIS_REFERENCED_ARRIVAL_INTERIOR_SOURCE ${EPHEMERIS_REFERENCED_ARRIVAL_INTERIOR_SOURCE} 
                              ${CMAKE_CURRENT_SOURCE_DIR}/EphemerisReferencedArrivalInterior.cpp
                              ${CMAKE_CURRENT_SOURCE_DIR}/EphemerisReferencedLTRendezvousInterior.cpp
                              ${CMAKE_CURRENT_SOURCE_DIR}/EphemerisReferencedArrivalWithVinfinityInterior.cpp
                              #${CMAKE_CURRENT_SOURCE_DIR}/EphemerisReferencedChemRendezvousInterior.cpp
                              ${CMAKE_CURRENT_SOURCE_DIR}/EphemerisReferencedInterceptInterior.cpp)

set(EPHEMERIS_REFERENCED_ARRIVAL_INTERIOR_FILES ${EPHEMERIS_REFERENCED_ARRIVAL_INTERIOR_FILES} 
                             ${EPHEMERIS_REFERENCED_ARRIVAL_INTERIOR_HEADERS} 
                             ${EPHEMERIS_REFERENCED_ARRIVAL_INTERIOR_SOURCE} PARENT_SCOPE)
