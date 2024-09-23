# main EMTG library files
set(CORE_HEADERS ${CORE_HEADERS} ${CMAKE_CURRENT_SOURCE_DIR}/chinchilla.h
                                 ${CMAKE_CURRENT_SOURCE_DIR}/doubleType.h
								 ${CMAKE_CURRENT_SOURCE_DIR}/EMTG_enums.h
								 ${CMAKE_CURRENT_SOURCE_DIR}/journeyoptions.h
								 ${CMAKE_CURRENT_SOURCE_DIR}/missionoptions.h
								 ${CMAKE_CURRENT_SOURCE_DIR}/problem.h)
						
set(CORE_SOURCE ${CORE_SOURCE} ${CMAKE_CURRENT_SOURCE_DIR}/journeyoptions.cpp
							   ${CMAKE_CURRENT_SOURCE_DIR}/missionoptions.cpp
							   ${CMAKE_CURRENT_SOURCE_DIR}/problem.cpp)
					
set(CORE_FILES ${CORE_FILES} ${CORE_HEADERS} ${CORE_SOURCE} PARENT_SCOPE)

# FBLT forward half-phase fixed-step files
if (RUN_FBLT_TESTBED)
    set(FBLT_TESTBED_CORE_HEADERS ${FBLT_TESTBED_CORE_HEADERS}
                                  ${CMAKE_CURRENT_SOURCE_DIR}/doubleType.h
                                  ${CMAKE_CURRENT_SOURCE_DIR}/EMTG_enums.h
								  ${CMAKE_CURRENT_SOURCE_DIR}/journeyoptions.h
								  ${CMAKE_CURRENT_SOURCE_DIR}/missionoptions.h)
                                  
    set(FBLT_TESTBED_CORE_SOURCE ${FBLT_TESTBED_CORE_SOURCE} 
                                 ${CMAKE_CURRENT_SOURCE_DIR}/journeyoptions.cpp
                                 ${CMAKE_CURRENT_SOURCE_DIR}/missionoptions.cpp)
                                 
    set(FBLT_TESTBED_CORE_FILES ${FBLT_TESTBED_CORE_FILES} 
                                ${FBLT_TESTBED_CORE_HEADERS} 
                                ${FBLT_TESTBED_CORE_SOURCE} PARENT_SCOPE)
endif()

# acceleration model unit test files
if (RUN_ACCELERATION_MODEL_TESTBED)
    set(ACCELERATION_MODEL_TESTBED_CORE_HEADERS ${ACCELERATION_MODEL_TESTBED_CORE_HEADERS}
                                                ${CMAKE_CURRENT_SOURCE_DIR}/doubleType.h
                                                ${CMAKE_CURRENT_SOURCE_DIR}/EMTG_enums.h
								                ${CMAKE_CURRENT_SOURCE_DIR}/journeyoptions.h
								                ${CMAKE_CURRENT_SOURCE_DIR}/missionoptions.h)
                                  
    set(ACCELERATION_MODEL_TESTBED_CORE_SOURCE ${ACCELERATION_MODEL_TESTBED_CORE_SOURCE} 
                                               ${CMAKE_CURRENT_SOURCE_DIR}/journeyoptions.cpp
                                               ${CMAKE_CURRENT_SOURCE_DIR}/missionoptions.cpp)
                                 
    set(ACCELERATION_MODEL_TESTBED_CORE_FILES ${ACCELERATION_MODEL_TESTBED_CORE_FILES} 
                                              ${ACCELERATION_MODEL_TESTBED_CORE_HEADERS} 
                                              ${ACCELERATION_MODEL_TESTBED_CORE_SOURCE} PARENT_SCOPE)
endif()