
set(FREE_POINT_ARRIVAL_HEADERS ${FREE_POINT_ARRIVAL_HEADERS} 
                               ${CMAKE_CURRENT_SOURCE_DIR}/FreePointArrival.h
                               ${CMAKE_CURRENT_SOURCE_DIR}/FreePointArrivalWithVinfinity.h
                               ${CMAKE_CURRENT_SOURCE_DIR}/FreePointChemRendezvous.h
                               ${CMAKE_CURRENT_SOURCE_DIR}/FreePointIntercept.h
                               ${CMAKE_CURRENT_SOURCE_DIR}/FreePointLTRendezvous.h)
																	  
set(FREE_POINT_ARRIVAL_SOURCE ${FREE_POINT_ARRIVAL_SOURCE} 
                              ${CMAKE_CURRENT_SOURCE_DIR}/FreePointArrival.cpp
                              ${CMAKE_CURRENT_SOURCE_DIR}/FreePointArrivalWithVinfinity.cpp
                              ${CMAKE_CURRENT_SOURCE_DIR}/FreePointChemRendezvous.cpp
                              ${CMAKE_CURRENT_SOURCE_DIR}/FreePointIntercept.cpp
                              ${CMAKE_CURRENT_SOURCE_DIR}/FreePointLTRendezvous.cpp)

set(FREE_POINT_ARRIVAL_FILES ${FREE_POINT_ARRIVAL_FILES} 
                             ${FREE_POINT_ARRIVAL_HEADERS} 
                             ${FREE_POINT_ARRIVAL_SOURCE} PARENT_SCOPE)
