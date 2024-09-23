# main EMTG library files
    
set(HARDWARE_MODEL_HEADERS ${HARDWARE_MODEL_HEADERS} ${CMAKE_CURRENT_SOURCE_DIR}/ChemicalPropulsionSystem.h
													 ${CMAKE_CURRENT_SOURCE_DIR}/ElectricPropulsionSystem.h
													 ${CMAKE_CURRENT_SOURCE_DIR}/HardwareBase.h
													 ${CMAKE_CURRENT_SOURCE_DIR}/LaunchVehicle.h
													 ${CMAKE_CURRENT_SOURCE_DIR}/LaunchVehicleOptions.h
													 ${CMAKE_CURRENT_SOURCE_DIR}/PowerSystem.h
													 ${CMAKE_CURRENT_SOURCE_DIR}/PowerSystemOptions.h
													 ${CMAKE_CURRENT_SOURCE_DIR}/PropulsionSystem.h
													 ${CMAKE_CURRENT_SOURCE_DIR}/PropulsionSystemOptions.h
													 ${CMAKE_CURRENT_SOURCE_DIR}/Spacecraft.h
													 ${CMAKE_CURRENT_SOURCE_DIR}/SpacecraftOptions.h
													 ${CMAKE_CURRENT_SOURCE_DIR}/Stage.h
													 ${CMAKE_CURRENT_SOURCE_DIR}/StageOptions.h
													 ${CMAKE_CURRENT_SOURCE_DIR}/ThrottleSetting.h
													 ${CMAKE_CURRENT_SOURCE_DIR}/ThrottleTable.h
                                                     ${CMAKE_CURRENT_SOURCE_DIR}/SpacecraftOptionsFactory.h
                                                     ${CMAKE_CURRENT_SOURCE_DIR}/LaunchVehicleOptionsFactory.h)
                             
set(HARDWARE_MODEL_SOURCE ${HARDWARE_MODEL_SOURCE} ${CMAKE_CURRENT_SOURCE_DIR}/ChemicalPropulsionSystem.cpp
												   ${CMAKE_CURRENT_SOURCE_DIR}/ElectricPropulsionSystem.cpp
												   ${CMAKE_CURRENT_SOURCE_DIR}/LaunchVehicle.cpp
												   ${CMAKE_CURRENT_SOURCE_DIR}/LaunchVehicleOptions.cpp
												   ${CMAKE_CURRENT_SOURCE_DIR}/PowerSystem.cpp
												   ${CMAKE_CURRENT_SOURCE_DIR}/PowerSystemOptions.cpp
												   ${CMAKE_CURRENT_SOURCE_DIR}/PropulsionSystemOptions.cpp
												   ${CMAKE_CURRENT_SOURCE_DIR}/Spacecraft.cpp
												   ${CMAKE_CURRENT_SOURCE_DIR}/SpacecraftOptions.cpp
												   ${CMAKE_CURRENT_SOURCE_DIR}/Stage.cpp
												   ${CMAKE_CURRENT_SOURCE_DIR}/StageOptions.cpp
												   ${CMAKE_CURRENT_SOURCE_DIR}/ThrottleSetting.cpp
												   ${CMAKE_CURRENT_SOURCE_DIR}/ThrottleTable.cpp
                                                   ${CMAKE_CURRENT_SOURCE_DIR}/SpacecraftOptionsFactory.cpp
                                                   ${CMAKE_CURRENT_SOURCE_DIR}/LaunchVehicleOptionsFactory.cpp)
                                                    
set(HARDWARE_MODEL_FILES ${HARDWARE_MODEL_FILES} ${HARDWARE_MODEL_HEADERS} ${HARDWARE_MODEL_SOURCE} PARENT_SCOPE)
    
# FBLT forward half-phase fixed-step files
if (RUN_FBLT_TESTBED)
    set(FBLT_TESTBED_HARDWARE_MODEL_HEADERS ${FBLT_TESTBED_HARDWARE_MODEL_HEADERS} 
                                            ${CMAKE_CURRENT_SOURCE_DIR}/ChemicalPropulsionSystem.h
										    ${CMAKE_CURRENT_SOURCE_DIR}/ElectricPropulsionSystem.h
										    ${CMAKE_CURRENT_SOURCE_DIR}/HardwareBase.h
										    ${CMAKE_CURRENT_SOURCE_DIR}/LaunchVehicle.h
										    ${CMAKE_CURRENT_SOURCE_DIR}/LaunchVehicleOptions.h
										    ${CMAKE_CURRENT_SOURCE_DIR}/PowerSystem.h
										    ${CMAKE_CURRENT_SOURCE_DIR}/PowerSystemOptions.h
										    ${CMAKE_CURRENT_SOURCE_DIR}/PropulsionSystem.h
										    ${CMAKE_CURRENT_SOURCE_DIR}/PropulsionSystemOptions.h
										    ${CMAKE_CURRENT_SOURCE_DIR}/Spacecraft.h
										    ${CMAKE_CURRENT_SOURCE_DIR}/SpacecraftOptions.h
										    ${CMAKE_CURRENT_SOURCE_DIR}/Stage.h
										    ${CMAKE_CURRENT_SOURCE_DIR}/StageOptions.h
										    ${CMAKE_CURRENT_SOURCE_DIR}/ThrottleSetting.h
										    ${CMAKE_CURRENT_SOURCE_DIR}/ThrottleTable.h
                                            ${CMAKE_CURRENT_SOURCE_DIR}/SpacecraftOptionsFactory.h
                                            ${CMAKE_CURRENT_SOURCE_DIR}/LaunchVehicleOptionsFactory.h)
                             
    set(FBLT_TESTBED_HARDWARE_MODEL_SOURCE ${FBLT_TESTBED_HARDWARE_MODEL_SOURCE} 
										   ${CMAKE_CURRENT_SOURCE_DIR}/ChemicalPropulsionSystem.cpp
                                           ${CMAKE_CURRENT_SOURCE_DIR}/ElectricPropulsionSystem.cpp
                                           ${CMAKE_CURRENT_SOURCE_DIR}/LaunchVehicle.cpp
                                           ${CMAKE_CURRENT_SOURCE_DIR}/LaunchVehicleOptions.cpp
                                           ${CMAKE_CURRENT_SOURCE_DIR}/PowerSystem.cpp
                                           ${CMAKE_CURRENT_SOURCE_DIR}/PowerSystemOptions.cpp
                                           ${CMAKE_CURRENT_SOURCE_DIR}/PropulsionSystemOptions.cpp
                                           ${CMAKE_CURRENT_SOURCE_DIR}/Spacecraft.cpp
                                           ${CMAKE_CURRENT_SOURCE_DIR}/SpacecraftOptions.cpp
                                           ${CMAKE_CURRENT_SOURCE_DIR}/Stage.cpp
                                           ${CMAKE_CURRENT_SOURCE_DIR}/StageOptions.cpp
                                           ${CMAKE_CURRENT_SOURCE_DIR}/ThrottleSetting.cpp
                                           ${CMAKE_CURRENT_SOURCE_DIR}/ThrottleTable.cpp
                                           ${CMAKE_CURRENT_SOURCE_DIR}/SpacecraftOptionsFactory.cpp
                                           ${CMAKE_CURRENT_SOURCE_DIR}/LaunchVehicleOptionsFactory.cpp)
                                                        
    set(FBLT_TESTBED_HARDWARE_MODEL_FILES ${FBLT_TESTBED_HARDWARE_MODEL_FILES} 
                                          ${FBLT_TESTBED_HARDWARE_MODEL_HEADERS} 
                                          ${FBLT_TESTBED_HARDWARE_MODEL_SOURCE} PARENT_SCOPE)
endif()

# acceleration model unit test files
if (RUN_ACCELERATION_MODEL_TESTBED)
    set(ACCELERATION_MODEL_TESTBED_HARDWARE_MODEL_HEADERS ${ACCELERATION_MODEL_TESTBED_HARDWARE_MODEL_HEADERS} 
                                                          ${CMAKE_CURRENT_SOURCE_DIR}/ChemicalPropulsionSystem.h
										                  ${CMAKE_CURRENT_SOURCE_DIR}/ElectricPropulsionSystem.h
										                  ${CMAKE_CURRENT_SOURCE_DIR}/HardwareBase.h
										                  ${CMAKE_CURRENT_SOURCE_DIR}/LaunchVehicle.h
										                  ${CMAKE_CURRENT_SOURCE_DIR}/LaunchVehicleOptions.h
										                  ${CMAKE_CURRENT_SOURCE_DIR}/PowerSystem.h
										                  ${CMAKE_CURRENT_SOURCE_DIR}/PowerSystemOptions.h
										                  ${CMAKE_CURRENT_SOURCE_DIR}/PropulsionSystem.h
										                  ${CMAKE_CURRENT_SOURCE_DIR}/PropulsionSystemOptions.h
										                  ${CMAKE_CURRENT_SOURCE_DIR}/Spacecraft.h
										                  ${CMAKE_CURRENT_SOURCE_DIR}/SpacecraftOptions.h
										                  ${CMAKE_CURRENT_SOURCE_DIR}/Stage.h
										                  ${CMAKE_CURRENT_SOURCE_DIR}/StageOptions.h
										                  ${CMAKE_CURRENT_SOURCE_DIR}/ThrottleSetting.h
										                  ${CMAKE_CURRENT_SOURCE_DIR}/ThrottleTable.h
                                                          ${CMAKE_CURRENT_SOURCE_DIR}/SpacecraftOptionsFactory.h
                                                          ${CMAKE_CURRENT_SOURCE_DIR}/LaunchVehicleOptionsFactory.h)
                             
    set(ACCELERATION_MODEL_TESTBED_HARDWARE_MODEL_SOURCE ${ACCELERATION_MODEL_TESTBED_HARDWARE_MODEL_SOURCE} 
										                 ${CMAKE_CURRENT_SOURCE_DIR}/ChemicalPropulsionSystem.cpp
                                                         ${CMAKE_CURRENT_SOURCE_DIR}/ElectricPropulsionSystem.cpp
                                                         ${CMAKE_CURRENT_SOURCE_DIR}/LaunchVehicle.cpp
                                                         ${CMAKE_CURRENT_SOURCE_DIR}/LaunchVehicleOptions.cpp
                                                         ${CMAKE_CURRENT_SOURCE_DIR}/PowerSystem.cpp
                                                         ${CMAKE_CURRENT_SOURCE_DIR}/PowerSystemOptions.cpp
                                                         ${CMAKE_CURRENT_SOURCE_DIR}/PropulsionSystemOptions.cpp
                                                         ${CMAKE_CURRENT_SOURCE_DIR}/Spacecraft.cpp
                                                         ${CMAKE_CURRENT_SOURCE_DIR}/SpacecraftOptions.cpp
                                                         ${CMAKE_CURRENT_SOURCE_DIR}/Stage.cpp
                                                         ${CMAKE_CURRENT_SOURCE_DIR}/StageOptions.cpp
                                                         ${CMAKE_CURRENT_SOURCE_DIR}/ThrottleSetting.cpp
                                                         ${CMAKE_CURRENT_SOURCE_DIR}/ThrottleTable.cpp
                                                         ${CMAKE_CURRENT_SOURCE_DIR}/SpacecraftOptionsFactory.cpp
                                                         ${CMAKE_CURRENT_SOURCE_DIR}/LaunchVehicleOptionsFactory.cpp)
                                                        
    set(ACCELERATION_MODEL_TESTBED_HARDWARE_MODEL_FILES ${ACCELERATION_MODEL_TESTBED_HARDWARE_MODEL_FILES} 
                                                        ${ACCELERATION_MODEL_TESTBED_HARDWARE_MODEL_HEADERS} 
                                                        ${ACCELERATION_MODEL_TESTBED_HARDWARE_MODEL_SOURCE} PARENT_SCOPE)
endif()