
set(ACCELERATION_MODEL_HEADERS  ${ACCELERATION_MODEL_HEADERS} ${CMAKE_CURRENT_SOURCE_DIR}/AccelerationModel.h
                                                              ${CMAKE_CURRENT_SOURCE_DIR}/CentralBodyGravityTerm.h
															  ${CMAKE_CURRENT_SOURCE_DIR}/GravityTerm.h
															  ${CMAKE_CURRENT_SOURCE_DIR}/SolarRadiationPressureTerm.h
                                                              ${CMAKE_CURRENT_SOURCE_DIR}/SpacecraftAccelerationModel.h
                                                              ${CMAKE_CURRENT_SOURCE_DIR}/SpacecraftAccelerationModelTerm.h
                                                              ${CMAKE_CURRENT_SOURCE_DIR}/SphericalHarmonicTerm.h
															  ${CMAKE_CURRENT_SOURCE_DIR}/ThrustTerm.h
                                                              ${CMAKE_CURRENT_SOURCE_DIR}/AerodynamicDragTerm.h)
																	  
set(ACCELERATION_MODEL_SOURCE  ${ACCELERATION_MODEL_SOURCE} ${CMAKE_CURRENT_SOURCE_DIR}/AccelerationModel.cpp
                                                            ${CMAKE_CURRENT_SOURCE_DIR}/CentralBodyGravityTerm.cpp
															${CMAKE_CURRENT_SOURCE_DIR}/GravityTerm.cpp
															${CMAKE_CURRENT_SOURCE_DIR}/SolarRadiationPressureTerm.cpp
                                                            ${CMAKE_CURRENT_SOURCE_DIR}/SpacecraftAccelerationModel.cpp
                                                            ${CMAKE_CURRENT_SOURCE_DIR}/SpacecraftAccelerationModelTerm.cpp
                                                            ${CMAKE_CURRENT_SOURCE_DIR}/SphericalHarmonicTerm.cpp
															${CMAKE_CURRENT_SOURCE_DIR}/ThrustTerm.cpp
                                                            ${CMAKE_CURRENT_SOURCE_DIR}/AerodynamicDragTerm.cpp)
    
set(ACCELERATION_MODEL_FILES ${ACCELERATION_MODEL_FILES} ${ACCELERATION_MODEL_HEADERS} ${ACCELERATION_MODEL_SOURCE} PARENT_SCOPE)


#if (RUNNING_FBLT_TESTBED)
#    set(FBLT_TESTBED_ASTRODYNAMICS_HEADERS ${FBLT_TESTBED_ASTRODYNAMICS_HEADERS}
#                                           ${CMAKE_CURRENT_SOURCE_DIR}/body.h
#                                           ${CMAKE_CURRENT_SOURCE_DIR}/engine_model.h
#                                           ${CMAKE_CURRENT_SOURCE_DIR}/FBLT_acceleration_model.h
#                                           ${CMAKE_CURRENT_SOURCE_DIR}/FBLT_EOM.h
#                                           ${CMAKE_CURRENT_SOURCE_DIR}/frame.h
#                                           ${CMAKE_CURRENT_SOURCE_DIR}/orbit_element_conversions.h
#                                           ${CMAKE_CURRENT_SOURCE_DIR}/universe.h)
#                                           
#    set(FBLT_TESTBED_ASTRODYNAMICS_SOURCE ${FBLT_TESTBED_ASTRODYNAMICS_SOURCE}
#                                          ${CMAKE_CURRENT_SOURCE_DIR}/body.cpp
#                                          ${CMAKE_CURRENT_SOURCE_DIR}/FBLT_acceleration_model.cpp
#                                          ${CMAKE_CURRENT_SOURCE_DIR}/FBLT_EOM.cpp
#                                          ${CMAKE_CURRENT_SOURCE_DIR}/frame.cpp
#                                          ${CMAKE_CURRENT_SOURCE_DIR}/universe.cpp)
#                                           
#    set(FBLT_TESTBED_ASTRODYNAMICS_FILES ${FBLT_TESTBED_ASTRODYNAMICS_FILES} 
#                                         ${FBLT_TESTBED_ASTRODYNAMICS_HEADERS}
#                                         ${FBLT_TESTBED_ASTRODYNAMICS_SOURCE} PARENT_SCOPE)
#endif()

# acceleration model unit test files
if (RUN_ACCELERATION_MODEL_TESTBED)
    set(ACCELERATION_MODEL_TESTBED_ACCELERATION_MODEL_HEADERS ${ACCELERATION_MODEL_TESTBED_ACCELERATION_MODEL_HEADERS} 
                                                              ${CMAKE_CURRENT_SOURCE_DIR}/AccelerationModel.h
                                                              ${CMAKE_CURRENT_SOURCE_DIR}/CentralBodyGravityTerm.h
                                                              ${CMAKE_CURRENT_SOURCE_DIR}/GravityTerm.h
                                                              ${CMAKE_CURRENT_SOURCE_DIR}/SolarRadiationPressureTerm.h
                                                              ${CMAKE_CURRENT_SOURCE_DIR}/SpacecraftAccelerationModel.h
                                                              ${CMAKE_CURRENT_SOURCE_DIR}/SpacecraftAccelerationModelTerm.h
                                                              ${CMAKE_CURRENT_SOURCE_DIR}/SphericalHarmonicTerm.h
                                                              ${CMAKE_CURRENT_SOURCE_DIR}/ThrustTerm.h
                                                              ${CMAKE_CURRENT_SOURCE_DIR}/AerodynamicDragTerm.h)
                                                                          
    set(ACCELERATION_MODEL_TESTBED_ACCELERATION_MODEL_SOURCE ${ACCELERATION_MODEL_TESTBED_ACCELERATION_MODEL_SOURCE} 
                                                             ${CMAKE_CURRENT_SOURCE_DIR}/AccelerationModel.cpp
                                                             ${CMAKE_CURRENT_SOURCE_DIR}/CentralBodyGravityTerm.cpp
                                                             ${CMAKE_CURRENT_SOURCE_DIR}/GravityTerm.cpp
                                                             ${CMAKE_CURRENT_SOURCE_DIR}/SolarRadiationPressureTerm.cpp
                                                             ${CMAKE_CURRENT_SOURCE_DIR}/SpacecraftAccelerationModel.cpp
                                                             ${CMAKE_CURRENT_SOURCE_DIR}/SpacecraftAccelerationModelTerm.cpp
                                                             ${CMAKE_CURRENT_SOURCE_DIR}/SphericalHarmonicTerm.cpp
                                                             ${CMAKE_CURRENT_SOURCE_DIR}/ThrustTerm.cpp
                                                             ${CMAKE_CURRENT_SOURCE_DIR}/AerodynamicDragTerm.cpp)
        
    set(ACCELERATION_MODEL_TESTBED_ACCELERATION_MODEL_FILES ${ACCELERATION_MODEL_TESTBED_ACCELERATION_MODEL_FILES} 
                                                            ${ACCELERATION_MODEL_TESTBED_ACCELERATION_MODEL_HEADERS} 
                                                            ${ACCELERATION_MODEL_TESTBED_ACCELERATION_MODEL_SOURCE} PARENT_SCOPE)
endif()