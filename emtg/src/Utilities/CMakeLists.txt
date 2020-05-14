# main EMTG library files
set(UTILITIES_HEADERS ${UTILITIES_HEADERS} ${CMAKE_CURRENT_SOURCE_DIR}/EMTG_solver_utilities.h
                                           ${CMAKE_CURRENT_SOURCE_DIR}/EMTG_string_utilities.h
                                           ${CMAKE_CURRENT_SOURCE_DIR}/EMTG_time_utilities.h
                                           ${CMAKE_CURRENT_SOURCE_DIR}/file_utilities.h
                                           ${CMAKE_CURRENT_SOURCE_DIR}/mjd_to_mdyhms.h                                           
                                           ${CMAKE_CURRENT_SOURCE_DIR}/sparsey_thing.h                                           
                                           ${CMAKE_CURRENT_SOURCE_DIR}/writey_thing.h
                                           ${CMAKE_CURRENT_SOURCE_DIR}/ephemeris_reader.h
                                           ${CMAKE_CURRENT_SOURCE_DIR}/covariance_reader.h
                                           ${CMAKE_CURRENT_SOURCE_DIR}/target_spec_line.h
                                           ${CMAKE_CURRENT_SOURCE_DIR}/maneuver_spec_line.h
                                           ${CMAKE_CURRENT_SOURCE_DIR}/maneuver_spec_item.h)
                    
set(UTILITIES_SOURCE ${UTILITIES_SOURCE} ${CMAKE_CURRENT_SOURCE_DIR}/EMTG_solver_utilities.cpp
                                         ${CMAKE_CURRENT_SOURCE_DIR}/EMTG_string_utilities.cpp
                                         ${CMAKE_CURRENT_SOURCE_DIR}/EMTG_time_utilities.cpp
                                         ${CMAKE_CURRENT_SOURCE_DIR}/file_utilities.cpp
                                         ${CMAKE_CURRENT_SOURCE_DIR}/mjd_to_mdyhms.cpp                                            
                                         ${CMAKE_CURRENT_SOURCE_DIR}/sparsey_thing.cpp                                          
                                         ${CMAKE_CURRENT_SOURCE_DIR}/writey_thing.cpp
                                         ${CMAKE_CURRENT_SOURCE_DIR}/ephemeris_reader.cpp
                                         ${CMAKE_CURRENT_SOURCE_DIR}/covariance_reader.cpp
                                         ${CMAKE_CURRENT_SOURCE_DIR}/target_spec_line.cpp
                                         ${CMAKE_CURRENT_SOURCE_DIR}/maneuver_spec_line.cpp
                                         ${CMAKE_CURRENT_SOURCE_DIR}/maneuver_spec_item.cpp)
                                         
set(UTILITIES_FILES ${UTILITIES_FILES} ${UTILITIES_HEADERS} ${UTILITIES_SOURCE} PARENT_SCOPE)


# FBLT forward half-phase fixed-step files
if (RUN_FBLT_TESTBED)
    set(FBLT_TESTBED_UTILITIES_HEADERS ${FBLT_TESTBED_UTILITIES_HEADERS} 
                                       ${CMAKE_CURRENT_SOURCE_DIR}/EMTG_time_utilities.h
                                       ${CMAKE_CURRENT_SOURCE_DIR}/file_utilities.h
                                       ${CMAKE_CURRENT_SOURCE_DIR}/ephemeris_reader.h)
                             
    set(FBLT_TESTBED_UTILITIES_SOURCE ${FBLT_TESTBED_UTILITIES_SOURCE} 
									  ${CMAKE_CURRENT_SOURCE_DIR}/EMTG_time_utilities.cpp
                                      ${CMAKE_CURRENT_SOURCE_DIR}/file_utilities.cpp
                                      ${CMAKE_CURRENT_SOURCE_DIR}/ephemeris_reader.cpp)
                                                        
    set(FBLT_TESTBED_UTILITIES_FILES ${FBLT_TESTBED_UTILITIES_FILES} 
                                     ${FBLT_TESTBED_UTILITIES_HEADERS} 
                                     ${FBLT_TESTBED_UTILITIES_SOURCE} PARENT_SCOPE)
endif()

# acceleration model unit test files
if (RUN_ACCELERATION_MODEL_TESTBED)
    set(ACCELERATION_MODEL_TESTBED_UTILITIES_HEADERS ${ACCELERATION_MODEL_TESTBED_UTILITIES_HEADERS} 
									                 ${CMAKE_CURRENT_SOURCE_DIR}/EMTG_time_utilities.h
                                                     ${CMAKE_CURRENT_SOURCE_DIR}/file_utilities.h
                                                     ${CMAKE_CURRENT_SOURCE_DIR}/ephemeris_reader.h)
												  
	set(ACCELERATION_MODEL_TESTBED_UTILITIES_SOURCE ${ACCELERATION_MODEL_TESTBED_UTILITIES_SOURCE} 
								                    ${CMAKE_CURRENT_SOURCE_DIR}/EMTG_time_utilities.cpp
                                                    ${CMAKE_CURRENT_SOURCE_DIR}/file_utilities.cpp
                                                    ${CMAKE_CURRENT_SOURCE_DIR}/ephemeris_reader.cpp) 
                                                        
    set(ACCELERATION_MODEL_TESTBED_UTILITIES_FILES ${ACCELERATION_MODEL_TESTBED_UTILITIES_FILES} 
                                                   ${ACCELERATION_MODEL_TESTBED_UTILITIES_HEADERS}
								                   ${ACCELERATION_MODEL_TESTBED_UTILITIES_SOURCE} PARENT_SCOPE)
endif()