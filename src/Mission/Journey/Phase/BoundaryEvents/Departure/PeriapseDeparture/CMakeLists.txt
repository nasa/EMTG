
set(PERIAPSE_DEPARTURE_HEADERS ${PERIAPSE_DEPARTURE_HEADERS} 
                                 ${CMAKE_CURRENT_SOURCE_DIR}/PeriapseDeparture.h
                                 ${CMAKE_CURRENT_SOURCE_DIR}/PeriapseLaunch.h)
                                 #${CMAKE_CURRENT_SOURCE_DIR}/PeriapseFreeDirectDeparture.h)
																	  
set(PERIAPSE_DEPARTURE_SOURCE ${PERIAPSE_DEPARTURE_SOURCE} 
                                ${CMAKE_CURRENT_SOURCE_DIR}/PeriapseDeparture.cpp
                                ${CMAKE_CURRENT_SOURCE_DIR}/PeriapseLaunch.cpp)
                                #${CMAKE_CURRENT_SOURCE_DIR}/PeriapseFreeDirectDeparture.cpp)

set(PERIAPSE_DEPARTURE_FILES ${PERIAPSE_DEPARTURE_FILES} 
                               ${PERIAPSE_DEPARTURE_HEADERS} 
                               ${PERIAPSE_DEPARTURE_SOURCE} PARENT_SCOPE)
