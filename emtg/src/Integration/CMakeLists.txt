# main EMTG library files
set(INTEGRATION_HEADERS ${INTEGRATION_HEADERS} ${CMAKE_CURRENT_SOURCE_DIR}/Integrand.h
                                               ${CMAKE_CURRENT_SOURCE_DIR}/IntegrationScheme.h
                                               ${CMAKE_CURRENT_SOURCE_DIR}/IntegrationSchemeFactory.h
						                       ${CMAKE_CURRENT_SOURCE_DIR}/ExplicitRungeKutta.h
						                       ${CMAKE_CURRENT_SOURCE_DIR}/IntegrationCoefficients.h
						                       ${CMAKE_CURRENT_SOURCE_DIR}/RungeKutta4Tableau.h
						                       ${CMAKE_CURRENT_SOURCE_DIR}/RungeKuttaDP87Tableau.h
						                       ${CMAKE_CURRENT_SOURCE_DIR}/RungeKuttaTableau.h)
                                               
set(INTEGRATION_SOURCE ${INTEGRATION_SOURCE} ${CMAKE_CURRENT_SOURCE_DIR}/Integrand.cpp 
                                             ${CMAKE_CURRENT_SOURCE_DIR}/IntegrationScheme.cpp
                                             ${CMAKE_CURRENT_SOURCE_DIR}/IntegrationSchemeFactory.cpp
						                     ${CMAKE_CURRENT_SOURCE_DIR}/ExplicitRungeKutta.cpp
						                     ${CMAKE_CURRENT_SOURCE_DIR}/IntegrationCoefficients.cpp
						                     ${CMAKE_CURRENT_SOURCE_DIR}/RungeKutta4Tableau.cpp
					                         ${CMAKE_CURRENT_SOURCE_DIR}/RungeKuttaDP87Tableau.cpp
						                     ${CMAKE_CURRENT_SOURCE_DIR}/RungeKuttaTableau.cpp)
                                               
set(INTEGRATION_FILES ${INTEGRATION_FILES} ${INTEGRATION_HEADERS} 
                                           ${INTEGRATION_SOURCE} PARENT_SCOPE)

# acceleration model unit test files
if (RUN_ACCELERATION_MODEL_TESTBED)
    set(ACCELERATION_MODEL_TESTBED_INTEGRATION_HEADERS ${ACCELERATION_MODEL_TESTBED_INTEGRATION_HEADERS}
                                                       ${CMAKE_CURRENT_SOURCE_DIR}/Integrand.h
                                                       ${CMAKE_CURRENT_SOURCE_DIR}/IntegrationScheme.h
							                           ${CMAKE_CURRENT_SOURCE_DIR}/ExplicitRungeKutta.h
							                           ${CMAKE_CURRENT_SOURCE_DIR}/IntegrationCoefficients.h
							                           ${CMAKE_CURRENT_SOURCE_DIR}/RungeKutta4Tableau.h
							                           ${CMAKE_CURRENT_SOURCE_DIR}/RungeKuttaDP87Tableau.h
							                           ${CMAKE_CURRENT_SOURCE_DIR}/RungeKuttaTableau.h)
                                                       
    set(ACCELERATION_MODEL_TESTBED_INTEGRATION_SOURCE ${ACCELERATION_MODEL_TESTBED_INTEGRATION_SOURCE}
                                                      ${CMAKE_CURRENT_SOURCE_DIR}/Integrand.cpp
                                                      ${CMAKE_CURRENT_SOURCE_DIR}/IntegrationScheme.cpp
							                          ${CMAKE_CURRENT_SOURCE_DIR}/ExplicitRungeKutta.cpp
							                          ${CMAKE_CURRENT_SOURCE_DIR}/IntegrationCoefficients.cpp
							                          ${CMAKE_CURRENT_SOURCE_DIR}/RungeKutta4Tableau.cpp
							                          ${CMAKE_CURRENT_SOURCE_DIR}/RungeKuttaDP87Tableau.cpp
							                          ${CMAKE_CURRENT_SOURCE_DIR}/RungeKuttaTableau.cpp)
                                                        
    set(ACCELERATION_MODEL_TESTBED_INTEGRATION_FILES ${ACCELERATION_MODEL_TESTBED_INTEGRATION_FILES} 
                                                     ${ACCELERATION_MODEL_TESTBED_INTEGRATION_HEADERS} 
                                                     ${ACCELERATION_MODEL_TESTBED_INTEGRATION_SOURCE} PARENT_SCOPE)
endif()