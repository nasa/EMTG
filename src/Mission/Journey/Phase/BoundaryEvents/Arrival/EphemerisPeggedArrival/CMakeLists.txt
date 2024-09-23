
set(EPHEMERIS_PEGGED_ARRIVAL_HEADERS ${EPHEMERIS_PEGGED_ARRIVAL_HEADERS} 
                                     ${CMAKE_CURRENT_SOURCE_DIR}/EphemerisPeggedArrival.h
                                     ${CMAKE_CURRENT_SOURCE_DIR}/EphemerisPeggedArrivalWithVinfinity.h
                                     ${CMAKE_CURRENT_SOURCE_DIR}/EphemerisPeggedChemMatchVinf.h
                                     ${CMAKE_CURRENT_SOURCE_DIR}/EphemerisPeggedChemRendezvous.h
                                     ${CMAKE_CURRENT_SOURCE_DIR}/EphemerisPeggedFlybyIn.h
                                     ${CMAKE_CURRENT_SOURCE_DIR}/EphemerisPeggedIntercept.h
									 ${CMAKE_CURRENT_SOURCE_DIR}/EphemerisPeggedMomentumTransfer.h
                                     ${CMAKE_CURRENT_SOURCE_DIR}/EphemerisPeggedLTMatchVinf.h
                                     ${CMAKE_CURRENT_SOURCE_DIR}/EphemerisPeggedLTRendezvous.h
                                     ${CMAKE_CURRENT_SOURCE_DIR}/EphemerisPeggedMatchVinf.h
                                     ${CMAKE_CURRENT_SOURCE_DIR}/EphemerisPeggedOrbitInsertion.h
                                     ${CMAKE_CURRENT_SOURCE_DIR}/EphemerisPeggedSpiralArrival.h)
																	  
set(EPHEMERIS_PEGGED_ARRIVAL_SOURCE ${EPHEMERIS_PEGGED_ARRIVAL_SOURCE} 
                                    ${CMAKE_CURRENT_SOURCE_DIR}/EphemerisPeggedArrival.cpp
                                    ${CMAKE_CURRENT_SOURCE_DIR}/EphemerisPeggedArrivalWithVinfinity.cpp
                                    ${CMAKE_CURRENT_SOURCE_DIR}/EphemerisPeggedChemRendezvous.cpp
                                    ${CMAKE_CURRENT_SOURCE_DIR}/EphemerisPeggedFlybyIn.cpp
                                    ${CMAKE_CURRENT_SOURCE_DIR}/EphemerisPeggedIntercept.cpp
									${CMAKE_CURRENT_SOURCE_DIR}/EphemerisPeggedMomentumTransfer.cpp
                                    ${CMAKE_CURRENT_SOURCE_DIR}/EphemerisPeggedLTRendezvous.cpp
                                    ${CMAKE_CURRENT_SOURCE_DIR}/EphemerisPeggedOrbitInsertion.cpp
                                    ${CMAKE_CURRENT_SOURCE_DIR}/EphemerisPeggedSpiralArrival.cpp)

set(EPHEMERIS_PEGGED_ARRIVAL_FILES ${EPHEMERIS_PEGGED_ARRIVAL_FILES} 
                                   ${EPHEMERIS_PEGGED_ARRIVAL_HEADERS} 
                                   ${EPHEMERIS_PEGGED_ARRIVAL_SOURCE} PARENT_SCOPE)
