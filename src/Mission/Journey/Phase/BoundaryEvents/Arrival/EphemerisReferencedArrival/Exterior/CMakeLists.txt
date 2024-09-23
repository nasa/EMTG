
set(EPHEMERIS_REFERENCED_ARRIVAL_EXTERIOR_HEADERS ${EPHEMERIS_REFERENCED_ARRIVAL_EXTERIOR_HEADERS} 
                               ${CMAKE_CURRENT_SOURCE_DIR}/EphemerisReferencedArrivalExterior.h
                               ${CMAKE_CURRENT_SOURCE_DIR}/EphemerisReferencedLTRendezvousExterior.h
                               ${CMAKE_CURRENT_SOURCE_DIR}/EphemerisReferencedArrivalWithVinfinityExterior.h
                               #${CMAKE_CURRENT_SOURCE_DIR}/EphemerisReferencedChemRendezvousExterior.h
                               ${CMAKE_CURRENT_SOURCE_DIR}/EphemerisReferencedInterceptExterior.h)

set(EPHEMERIS_REFERENCED_ARRIVAL_EXTERIOR_SOURCE ${EPHEMERIS_REFERENCED_ARRIVAL_EXTERIOR_SOURCE} 
                              ${CMAKE_CURRENT_SOURCE_DIR}/EphemerisReferencedArrivalExterior.cpp
                              ${CMAKE_CURRENT_SOURCE_DIR}/EphemerisReferencedLTRendezvousExterior.cpp
                              ${CMAKE_CURRENT_SOURCE_DIR}/EphemerisReferencedArrivalWithVinfinityExterior.cpp
                              #${CMAKE_CURRENT_SOURCE_DIR}/EphemerisReferencedChemRendezvousExterior.cpp
                              ${CMAKE_CURRENT_SOURCE_DIR}/EphemerisReferencedInterceptExterior.cpp)

set(EPHEMERIS_REFERENCED_ARRIVAL_EXTERIOR_FILES ${EPHEMERIS_REFERENCED_ARRIVAL_EXTERIOR_FILES} 
                             ${EPHEMERIS_REFERENCED_ARRIVAL_EXTERIOR_HEADERS} 
                             ${EPHEMERIS_REFERENCED_ARRIVAL_EXTERIOR_SOURCE} PARENT_SCOPE)
