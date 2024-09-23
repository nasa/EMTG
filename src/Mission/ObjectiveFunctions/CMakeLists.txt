
set(OBJECTIVE_FUNCTIONS_HEADERS  ${OBJECTIVE_FUNCTIONS_HEADERS}
                                 ${CMAKE_CURRENT_SOURCE_DIR}/ObjectiveFunctionFactory.h
                                 ${CMAKE_CURRENT_SOURCE_DIR}/MinimizeDeltavObjective.h
                                 ${CMAKE_CURRENT_SOURCE_DIR}/MaximizeMassObjective.h
                                 ${CMAKE_CURRENT_SOURCE_DIR}/MaximizeLogeMassObjective.h
                                 ${CMAKE_CURRENT_SOURCE_DIR}/MaximizeLog10MassObjective.h
                                 ${CMAKE_CURRENT_SOURCE_DIR}/MinimizeTimeObjective.h
                                 ${CMAKE_CURRENT_SOURCE_DIR}/MaximizeInitialMassObjective.h
                                 ${CMAKE_CURRENT_SOURCE_DIR}/ObjectiveFunctionBase.h
                                 ${CMAKE_CURRENT_SOURCE_DIR}/ArriveAsEarlyAsPossibleObjective.h
                                 ${CMAKE_CURRENT_SOURCE_DIR}/ArriveAsLateAsPossibleObjective.h
                                 ${CMAKE_CURRENT_SOURCE_DIR}/DepartAsEarlyAsPossibleObjective.h
                                 ${CMAKE_CURRENT_SOURCE_DIR}/DepartAsLateAsPossibleObjective.h
                                 ${CMAKE_CURRENT_SOURCE_DIR}/MaximizeDryMassObjective.h
                                 ${CMAKE_CURRENT_SOURCE_DIR}/MaximizeLog10DryMassObjective.h
                                 ${CMAKE_CURRENT_SOURCE_DIR}/MaximizeLogeDryMassObjective.h
                                 ${CMAKE_CURRENT_SOURCE_DIR}/MinimizeChemicalFuelObjective.h
                                 ${CMAKE_CURRENT_SOURCE_DIR}/MinimizeChemicalOxidizerObjective.h
                                 ${CMAKE_CURRENT_SOURCE_DIR}/MinimizeElectricPropellantObjective.h
                                 ${CMAKE_CURRENT_SOURCE_DIR}/MinimizeTotalPropellantObjective.h
                                 ${CMAKE_CURRENT_SOURCE_DIR}/MinimizeInitialImpulseObjective.h
                                 ${CMAKE_CURRENT_SOURCE_DIR}/MinimizeWaypointTrackingErrorObjective.h
                                 ${CMAKE_CURRENT_SOURCE_DIR}/MaximizeDistanceFromCentralBodyObjective.h)
																	  
set(OBJECTIVE_FUNCTIONS_SOURCE  ${OBJECTIVE_FUNCTIONS_SOURCE}
                                ${CMAKE_CURRENT_SOURCE_DIR}/ObjectiveFunctionFactory.cpp
                                ${CMAKE_CURRENT_SOURCE_DIR}/MinimizeDeltavObjective.cpp
                                ${CMAKE_CURRENT_SOURCE_DIR}/MaximizeMassObjective.cpp
                                ${CMAKE_CURRENT_SOURCE_DIR}/MaximizeLogeMassObjective.cpp
                                ${CMAKE_CURRENT_SOURCE_DIR}/MaximizeLog10MassObjective.cpp
                                ${CMAKE_CURRENT_SOURCE_DIR}/MinimizeTimeObjective.cpp
                                ${CMAKE_CURRENT_SOURCE_DIR}/MaximizeInitialMassObjective.cpp
                                ${CMAKE_CURRENT_SOURCE_DIR}/ObjectiveFunctionBase.cpp
                                ${CMAKE_CURRENT_SOURCE_DIR}/ArriveAsEarlyAsPossibleObjective.cpp
                                ${CMAKE_CURRENT_SOURCE_DIR}/ArriveAsLateAsPossibleObjective.cpp
                                ${CMAKE_CURRENT_SOURCE_DIR}/DepartAsEarlyAsPossibleObjective.cpp
                                ${CMAKE_CURRENT_SOURCE_DIR}/DepartAsLateAsPossibleObjective.cpp
                                ${CMAKE_CURRENT_SOURCE_DIR}/MaximizeDryMassObjective.cpp
                                ${CMAKE_CURRENT_SOURCE_DIR}/MaximizeLog10DryMassObjective.cpp
                                ${CMAKE_CURRENT_SOURCE_DIR}/MaximizeLogeDryMassObjective.cpp
                                ${CMAKE_CURRENT_SOURCE_DIR}/MinimizeChemicalFuelObjective.cpp
                                ${CMAKE_CURRENT_SOURCE_DIR}/MinimizeChemicalOxidizerObjective.cpp
                                ${CMAKE_CURRENT_SOURCE_DIR}/MinimizeElectricPropellantObjective.cpp
                                ${CMAKE_CURRENT_SOURCE_DIR}/MinimizeTotalPropellantObjective.cpp
                                ${CMAKE_CURRENT_SOURCE_DIR}/MinimizeInitialImpulseObjective.cpp
                                ${CMAKE_CURRENT_SOURCE_DIR}/MinimizeWaypointTrackingErrorObjective.cpp
                                ${CMAKE_CURRENT_SOURCE_DIR}/MaximizeDistanceFromCentralBodyObjective.cpp)

set(OBJECTIVE_FUNCTIONS_FILES ${OBJECTIVE_FUNCTIONS_FILES} 
                              ${OBJECTIVE_FUNCTIONS_HEADERS} 
                              ${OBJECTIVE_FUNCTIONS_SOURCE} PARENT_SCOPE)
