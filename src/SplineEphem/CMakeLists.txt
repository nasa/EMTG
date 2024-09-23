
set(SPLINE_EPHEM_HEADERS ${SPLINE_EPHEM_HEADERS} ${CMAKE_CURRENT_SOURCE_DIR}/SplineEphem_body.h
                                                 ${CMAKE_CURRENT_SOURCE_DIR}/SplineEphem_universe.h)
                                                 
set(SPLINE_EPHEM_SOURCE ${SPLINE_EPHEM_SOURCE} ${CMAKE_CURRENT_SOURCE_DIR}/SplineEphem_body.cpp
                                               ${CMAKE_CURRENT_SOURCE_DIR}/SplineEphem_universe.cpp) 
set(SPLINE_EPHEM_FILES ${SPLINE_EPHEM_FILES} ${SPLINE_EPHEM_HEADERS} ${SPLINE_EPHEM_SOURCE} PARENT_SCOPE)

# FBLT forward half-phase fixed-step files
if (RUN_FBLT_TESTBED)
	set(FBLT_TESTBED_SPLINE_EPHEM_HEADERS ${FBLT_TESTBED_SPLINE_EPHEM_HEADERS} 
										  ${CMAKE_CURRENT_SOURCE_DIR}/SplineEphem_body.h
										  ${CMAKE_CURRENT_SOURCE_DIR}/SplineEphem_universe.h)
													 
	set(FBLT_TESTBED_SPLINE_EPHEM_SOURCE ${FBLT_TESTBED_SPLINE_EPHEM_SOURCE} 
										 ${CMAKE_CURRENT_SOURCE_DIR}/SplineEphem_body.cpp
										 ${CMAKE_CURRENT_SOURCE_DIR}/SplineEphem_universe.cpp) 
	set(FBLT_TESTBED_SPLINE_EPHEM_FILES ${FBLT_TESTBED_SPLINE_EPHEM_FILES} 
						 			    ${FBLT_TESTBED_SPLINE_EPHEM_HEADERS} 
						         		${FBLT_TESTBED_SPLINE_EPHEM_SOURCE} PARENT_SCOPE)
endif()
   
# acceleration model unit test files
if (RUN_ACCELERATION_MODEL_TESTBED)
    set(ACCELERATION_MODEL_TESTBED_SPLINE_EPHEM_HEADERS ${ACCELERATION_MODEL_TESTBED_SPLINE_EPHEM_HEADERS} 
									                    ${CMAKE_CURRENT_SOURCE_DIR}/SplineEphem_body.h
										                ${CMAKE_CURRENT_SOURCE_DIR}/SplineEphem_universe.h)
												  
	set(ACCELERATION_MODEL_TESTBED_SPLINE_EPHEM_SOURCE ${ACCELERATION_MODEL_TESTBED_SPLINE_EPHEM_SOURCE} 
								                       ${CMAKE_CURRENT_SOURCE_DIR}/SplineEphem_body.cpp
										               ${CMAKE_CURRENT_SOURCE_DIR}/SplineEphem_universe.cpp) 
                                                        
    set(ACCELERATION_MODEL_TESTBED_SPLINE_EPHEM_FILES ${ACCELERATION_MODEL_TESTBED_SPLINE_EPHEM_FILES} 
                                                      ${ACCELERATION_MODEL_TESTBED_SPLINE_EPHEM_HEADERS}
								                      ${ACCELERATION_MODEL_TESTBED_SPLINE_EPHEM_SOURCE} PARENT_SCOPE)
endif()