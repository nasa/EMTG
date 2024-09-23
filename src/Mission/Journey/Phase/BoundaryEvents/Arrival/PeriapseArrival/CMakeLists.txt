
set(PERIAPSE_ARRIVAL_HEADERS ${PERIAPSE_ARRIVAL_HEADERS} 
                             ${CMAKE_CURRENT_SOURCE_DIR}/PeriapseArrival.h
                             ${CMAKE_CURRENT_SOURCE_DIR}/PeriapseFlybyIn.h)
																	  
set(PERIAPSE_ARRIVAL_SOURCE ${PERIAPSE_ARRIVAL_SOURCE} 
                            ${CMAKE_CURRENT_SOURCE_DIR}/PeriapseArrival.cpp
                            ${CMAKE_CURRENT_SOURCE_DIR}/PeriapseFlybyIn.cpp)

set(PERIAPSE_ARRIVAL_FILES ${PERIAPSE_ARRIVAL_FILES} 
                           ${PERIAPSE_ARRIVAL_HEADERS} 
                           ${PERIAPSE_ARRIVAL_SOURCE} PARENT_SCOPE)
