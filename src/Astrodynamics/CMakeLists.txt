add_subdirectory(AccelerationModel)
add_subdirectory(EquationsOfMotion)
add_subdirectory(StateRepresentation)

set(ASTRODYNAMICS_HEADERS  ${ASTRODYNAMICS_HEADERS} ${CMAKE_CURRENT_SOURCE_DIR}/body.h
                                                    ${CMAKE_CURRENT_SOURCE_DIR}/CentralBody.h
							                        ${CMAKE_CURRENT_SOURCE_DIR}/bplane.h
                                                    ${CMAKE_CURRENT_SOURCE_DIR}/frame.h
                                                    ${CMAKE_CURRENT_SOURCE_DIR}/orbit_element_conversions.h
                                                    ${CMAKE_CURRENT_SOURCE_DIR}/universe.h
                                                    ${CMAKE_CURRENT_SOURCE_DIR}/atmosphere.h
                                                    ${CMAKE_CURRENT_SOURCE_DIR}/ExponentialAtmosphere.h
                                                    ${CMAKE_CURRENT_SOURCE_DIR}/BodydeticConversions.h
                                                    )
																	  
set(ASTRODYNAMICS_SOURCE  ${ASTRODYNAMICS_SOURCE} ${CMAKE_CURRENT_SOURCE_DIR}/universe.cpp
                                                  ${CMAKE_CURRENT_SOURCE_DIR}/body.cpp
                                                  ${CMAKE_CURRENT_SOURCE_DIR}/bplane.cpp
                                                  ${CMAKE_CURRENT_SOURCE_DIR}/CentralBody.cpp
                                                  ${CMAKE_CURRENT_SOURCE_DIR}/frame.cpp
                                                  ${CMAKE_CURRENT_SOURCE_DIR}/atmosphere.cpp
                                                  ${CMAKE_CURRENT_SOURCE_DIR}/ExponentialAtmosphere.cpp
                                                  ${CMAKE_CURRENT_SOURCE_DIR}/BodydeticConversions.cpp
                                                  )
                                                  
# add the proprietary engine models if you have them
if(HAS_BUILT_IN_THRUSTERS)
    set(ASTRODYNAMICS_HEADERS  ${ASTRODYNAMICS_HEADERS} ${CMAKE_CURRENT_SOURCE_DIR}/engine_model.h)
endif()
    
set(ASTRODYNAMICS_FILES ${ASTRODYNAMICS_FILES} ${ASTRODYNAMICS_HEADERS} 
                                               ${ASTRODYNAMICS_SOURCE} 
                                               ${ACCELERATION_MODEL_FILES} 
                                               ${STATEREPRESENTATION_FILES}
                                               ${EOM_FILES} PARENT_SCOPE)

# FBLT forward half-phase fixed-step files
if (RUN_FBLT_TESTBED)
    set(FBLT_TESTBED_ASTRODYNAMICS_HEADERS ${FBLT_TESTBED_ASTRODYNAMICS_HEADERS}
                                           ${CMAKE_CURRENT_SOURCE_DIR}/body.h
                                           ${CMAKE_CURRENT_SOURCE_DIR}/CentralBody.h
							               ${CMAKE_CURRENT_SOURCE_DIR}/bplane.h
                                           ${CMAKE_CURRENT_SOURCE_DIR}/FBLT_acceleration_model.h
                                           ${CMAKE_CURRENT_SOURCE_DIR}/FBLT_EOM.h
                                           ${CMAKE_CURRENT_SOURCE_DIR}/frame.h
                                           ${CMAKE_CURRENT_SOURCE_DIR}/orbit_element_conversions.h
                                           ${CMAKE_CURRENT_SOURCE_DIR}/universe.h)
                                           
    set(FBLT_TESTBED_ASTRODYNAMICS_SOURCE ${FBLT_TESTBED_ASTRODYNAMICS_SOURCE}
                                          ${CMAKE_CURRENT_SOURCE_DIR}/body.cpp
                                          ${CMAKE_CURRENT_SOURCE_DIR}/bplane.cpp
                                          ${CMAKE_CURRENT_SOURCE_DIR}/CentralBody.cpp
                                          ${CMAKE_CURRENT_SOURCE_DIR}/FBLT_acceleration_model.cpp
                                          ${CMAKE_CURRENT_SOURCE_DIR}/FBLT_EOM.cpp
                                          ${CMAKE_CURRENT_SOURCE_DIR}/frame.cpp
                                          ${CMAKE_CURRENT_SOURCE_DIR}/universe.cpp)
                                           
    set(FBLT_TESTBED_ASTRODYNAMICS_FILES ${FBLT_TESTBED_ASTRODYNAMICS_FILES} 
                                         ${FBLT_TESTBED_ASTRODYNAMICS_HEADERS}
                                         ${FBLT_TESTBED_ASTRODYNAMICS_SOURCE} PARENT_SCOPE)
endif()

# acceleration model unit test files
if (RUN_ACCELERATION_MODEL_TESTBED)
    set(ACCELERATION_MODEL_TESTBED_ASTRODYNAMICS_HEADERS ${ACCELERATION_MODEL_TESTBED_ASTRODYNAMICS_HEADERS}
                                                         ${CMAKE_CURRENT_SOURCE_DIR}/body.h
                                                         ${CMAKE_CURRENT_SOURCE_DIR}/CentralBody.h
							                             ${CMAKE_CURRENT_SOURCE_DIR}/bplane.h
                                                         ${CMAKE_CURRENT_SOURCE_DIR}/frame.h
                                                         ${CMAKE_CURRENT_SOURCE_DIR}/orbit_element_conversions.h
                                                         ${CMAKE_CURRENT_SOURCE_DIR}/universe.h)
                                           
    set(ACCELERATION_MODEL_TESTBED_ASTRODYNAMICS_SOURCE ${ACCELERATION_MODEL_TESTBED_ASTRODYNAMICS_SOURCE}
                                                        ${CMAKE_CURRENT_SOURCE_DIR}/body.cpp
                                                        ${CMAKE_CURRENT_SOURCE_DIR}/bplane.cpp
                                                        ${CMAKE_CURRENT_SOURCE_DIR}/CentralBody.cpp
                                                        ${CMAKE_CURRENT_SOURCE_DIR}/frame.cpp
                                                        ${CMAKE_CURRENT_SOURCE_DIR}/universe.cpp)
                                           
    set(ACCELERATION_MODEL_TESTBED_ASTRODYNAMICS_FILES ${ACCELERATION_MODEL_TESTBED_ASTRODYNAMICS_FILES} 
                                                       ${ACCELERATION_MODEL_TESTBED_ASTRODYNAMICS_HEADERS}
                                                       ${ACCELERATION_MODEL_TESTBED_ASTRODYNAMICS_SOURCE}
                                                       ${ACCELERATION_MODEL_TESTBED_ACCELERATION_MODEL_FILES}
                                                       ${ACCELERATION_MODEL_TESTBED_EOM_FILES} PARENT_SCOPE)
endif()