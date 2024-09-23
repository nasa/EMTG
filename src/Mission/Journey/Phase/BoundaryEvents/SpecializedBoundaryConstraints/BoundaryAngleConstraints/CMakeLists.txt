
set(BOUNDARY_ANGLE_CONSTRAINTS_HEADERS ${BOUNDARY_ANGLE_CONSTRAINTS_HEADERS} 
                                       ${CMAKE_CURRENT_SOURCE_DIR}/RBP.h
                                       ${CMAKE_CURRENT_SOURCE_DIR}/RPB.h
                                       ${CMAKE_CURRENT_SOURCE_DIR}/RRP.h
                                       ${CMAKE_CURRENT_SOURCE_DIR}/RPR.h
                                       ${CMAKE_CURRENT_SOURCE_DIR}/angularMomentumReferenceAngle.h)
																	  
set(BOUNDARY_ANGLE_CONSTRAINTS_SOURCE ${BOUNDARY_ANGLE_CONSTRAINTS_SOURCE} 
                                       ${CMAKE_CURRENT_SOURCE_DIR}/RBP.cpp
                                       ${CMAKE_CURRENT_SOURCE_DIR}/RPB.cpp
                                       ${CMAKE_CURRENT_SOURCE_DIR}/RRP.cpp
                                       ${CMAKE_CURRENT_SOURCE_DIR}/RPR.cpp
                                       ${CMAKE_CURRENT_SOURCE_DIR}/angularMomentumReferenceAngle.cpp)

set(BOUNDARY_ANGLE_CONSTRAINTS_FILES ${BOUNDARY_ANGLE_CONSTRAINTS_FILES} 
                                     ${BOUNDARY_ANGLE_CONSTRAINTS_HEADERS} 
                                     ${BOUNDARY_ANGLE_CONSTRAINTS_SOURCE} PARENT_SCOPE)
