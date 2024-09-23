add_subdirectory(EphemerisPeggedArrival)
add_subdirectory(EphemerisReferencedArrival)
add_subdirectory(FreePointArrival)
add_subdirectory(PeriapseArrival)

set(ARRIVAL_HEADERS ${ARRIVAL_HEADERS} ${CMAKE_CURRENT_SOURCE_DIR}/arrival.h)
																	  
set(ARRIVAL_SOURCE ${ARRIVAL_SOURCE} ${CMAKE_CURRENT_SOURCE_DIR}/arrival.cpp)

set(ARRIVAL_FILES ${ARRIVAL_FILES} ${ARRIVAL_HEADERS} 
                                   ${ARRIVAL_SOURCE} 
                                   ${EPHEMERIS_PEGGED_ARRIVAL_FILES} 
                                   ${EPHEMERIS_REFERENCED_ARRIVAL_FILES}
								   ${PERIAPSE_ARRIVAL_FILES}
                                   ${FREE_POINT_ARRIVAL_FILES} PARENT_SCOPE)
