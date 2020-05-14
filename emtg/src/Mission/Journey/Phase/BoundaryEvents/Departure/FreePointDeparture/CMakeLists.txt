
set(FREE_POINT_DEPARTURE_HEADERS ${FREE_POINT_DEPARTURE_HEADERS} 
                                 ${CMAKE_CURRENT_SOURCE_DIR}/FreePointDeparture.h
                                 ${CMAKE_CURRENT_SOURCE_DIR}/FreePointFreeDirectDeparture.h
                                 ${CMAKE_CURRENT_SOURCE_DIR}/FreePointDirectInsertion.h)
																	  
set(FREE_POINT_DEPARTURE_SOURCE ${FREE_POINT_DEPARTURE_SOURCE} 
                                ${CMAKE_CURRENT_SOURCE_DIR}/FreePointDeparture.cpp
                                ${CMAKE_CURRENT_SOURCE_DIR}/FreePointFreeDirectDeparture.cpp
                                ${CMAKE_CURRENT_SOURCE_DIR}/FreePointDirectInsertion.cpp)

set(FREE_POINT_DEPARTURE_FILES ${FREE_POINT_DEPARTURE_FILES} 
                               ${FREE_POINT_DEPARTURE_HEADERS} 
                               ${FREE_POINT_DEPARTURE_SOURCE} PARENT_SCOPE)
