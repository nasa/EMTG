
set(EPHEMERIS_PEGGED_DEPARTURE_HEADERS ${EPHEMERIS_PEGGED_DEPARTURE_HEADERS} 
                                       ${CMAKE_CURRENT_SOURCE_DIR}/EphemerisPeggedDeparture.h
                                       ${CMAKE_CURRENT_SOURCE_DIR}/EphemerisPeggedFlybyOut.h
                                       ${CMAKE_CURRENT_SOURCE_DIR}/EphemerisPeggedFreeDirectDeparture.h
                                       ${CMAKE_CURRENT_SOURCE_DIR}/EphemerisPeggedLaunchDirectInsertion.h
                                       ${CMAKE_CURRENT_SOURCE_DIR}/EphemerisPeggedPoweredFlyby.h
                                       ${CMAKE_CURRENT_SOURCE_DIR}/EphemerisPeggedUnpoweredFlyby.h
								       ${CMAKE_CURRENT_SOURCE_DIR}/EphemerisPeggedZeroTurnFlyby.h
                                       ${CMAKE_CURRENT_SOURCE_DIR}/EphemerisPeggedSpiralDeparture.h)
																	  
set(EPHEMERIS_PEGGED_DEPARTURE_SOURCE ${EPHEMERIS_PEGGED_DEPARTURE_SOURCE} 
                                      ${CMAKE_CURRENT_SOURCE_DIR}/EphemerisPeggedDeparture.cpp
                                      ${CMAKE_CURRENT_SOURCE_DIR}/EphemerisPeggedFlybyOut.cpp
                                      ${CMAKE_CURRENT_SOURCE_DIR}/EphemerisPeggedFreeDirectDeparture.cpp
                                      ${CMAKE_CURRENT_SOURCE_DIR}/EphemerisPeggedLaunchDirectInsertion.cpp
                                      ${CMAKE_CURRENT_SOURCE_DIR}/EphemerisPeggedPoweredFlyby.cpp
                                      ${CMAKE_CURRENT_SOURCE_DIR}/EphemerisPeggedUnpoweredFlyby.cpp
								      ${CMAKE_CURRENT_SOURCE_DIR}/EphemerisPeggedZeroTurnFlyby.cpp
                                      ${CMAKE_CURRENT_SOURCE_DIR}/EphemerisPeggedSpiralDeparture.cpp)

set(EPHEMERIS_PEGGED_DEPARTURE_FILES ${EPHEMERIS_PEGGED_DEPARTURE_FILES} 
                                     ${EPHEMERIS_PEGGED_DEPARTURE_HEADERS} 
                                     ${EPHEMERIS_PEGGED_DEPARTURE_SOURCE} PARENT_SCOPE)
