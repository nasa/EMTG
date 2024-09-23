
set(ORBIT_ELEMENTS_CONSTRAINTS_HEADERS ${ORBIT_ELEMENTS_CONSTRAINTS_HEADERS} 
                                       ${CMAKE_CURRENT_SOURCE_DIR}/AOPConstraint.h
                                       ${CMAKE_CURRENT_SOURCE_DIR}/ApoapseDistanceConstraint.h
                                       ${CMAKE_CURRENT_SOURCE_DIR}/EccentricityConstraint.h
                                       ${CMAKE_CURRENT_SOURCE_DIR}/InclinationConstraint.h
                                       ${CMAKE_CURRENT_SOURCE_DIR}/OrbitElementConstraintBase.h
                                       ${CMAKE_CURRENT_SOURCE_DIR}/OrbitPeriodConstraint.h
                                       ${CMAKE_CURRENT_SOURCE_DIR}/PeriapseDistanceConstraint.h
                                       ${CMAKE_CURRENT_SOURCE_DIR}/RAANConstraint.h
                                       ${CMAKE_CURRENT_SOURCE_DIR}/SemimajorAxisConstraint.h
                                       ${CMAKE_CURRENT_SOURCE_DIR}/TrueAnomalyConstraint.h)
																	  
set(ORBIT_ELEMENTS_CONSTRAINTS_SOURCE ${ORBIT_ELEMENTS_CONSTRAINTS_SOURCE} 
                                      ${CMAKE_CURRENT_SOURCE_DIR}/AOPConstraint.cpp
                                      ${CMAKE_CURRENT_SOURCE_DIR}/ApoapseDistanceConstraint.cpp
                                      ${CMAKE_CURRENT_SOURCE_DIR}/EccentricityConstraint.cpp
                                      ${CMAKE_CURRENT_SOURCE_DIR}/InclinationConstraint.cpp
                                      ${CMAKE_CURRENT_SOURCE_DIR}/OrbitElementConstraintBase.cpp
                                      ${CMAKE_CURRENT_SOURCE_DIR}/OrbitPeriodConstraint.cpp
                                      ${CMAKE_CURRENT_SOURCE_DIR}/PeriapseDistanceConstraint.cpp
                                      ${CMAKE_CURRENT_SOURCE_DIR}/RAANConstraint.cpp
                                      ${CMAKE_CURRENT_SOURCE_DIR}/SemimajorAxisConstraint.cpp
                                      ${CMAKE_CURRENT_SOURCE_DIR}/TrueAnomalyConstraint.cpp)

set(ORBIT_ELEMENTS_CONSTRAINTS_FILES ${ORBIT_ELEMENTS_CONSTRAINTS_FILES} 
                                     ${ORBIT_ELEMENTS_CONSTRAINTS_HEADERS} 
                                     ${ORBIT_ELEMENTS_CONSTRAINTS_SOURCE} PARENT_SCOPE)
