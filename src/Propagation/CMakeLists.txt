
set(PROPAGATION_HEADERS ${PROPAGATION_HEADERS} ${CMAKE_CURRENT_SOURCE_DIR}/IntegratedAdaptiveStepPropagator.h
                                               ${CMAKE_CURRENT_SOURCE_DIR}/IntegratedFixedStepPropagator.h
                                               ${CMAKE_CURRENT_SOURCE_DIR}/IntegratedPropagator.h
                                               ${CMAKE_CURRENT_SOURCE_DIR}/KeplerPropagator.h
                                               ${CMAKE_CURRENT_SOURCE_DIR}/KeplerPropagatorTimeDomain.h
                                               ${CMAKE_CURRENT_SOURCE_DIR}/PropagatorBase.h
                                               ${CMAKE_CURRENT_SOURCE_DIR}/PropagatorFactory.h)
												  
set(PROPAGATION_SOURCE ${PROPAGATION_SOURCE} ${CMAKE_CURRENT_SOURCE_DIR}/IntegratedAdaptiveStepPropagator.cpp
                                             ${CMAKE_CURRENT_SOURCE_DIR}/IntegratedPropagator.cpp
                                             ${CMAKE_CURRENT_SOURCE_DIR}/IntegratedFixedStepPropagator.cpp
                                             ${CMAKE_CURRENT_SOURCE_DIR}/KeplerPropagator.cpp
                                             ${CMAKE_CURRENT_SOURCE_DIR}/KeplerPropagatorTimeDomain.cpp
                                             ${CMAKE_CURRENT_SOURCE_DIR}/PropagatorBase.cpp
                                             ${CMAKE_CURRENT_SOURCE_DIR}/PropagatorFactory.cpp)
                                           
set(PROPAGATION_FILES ${PROPAGATION_FILES} ${PROPAGATION_HEADERS} ${PROPAGATION_SOURCE} PARENT_SCOPE)
    

if (RUN_ACCELERATION_MODEL_TESTBED)
    
    set(ACCELERATION_MODEL_TESTBED_PROPAGATION_HEADERS ${ACCELERATION_MODEL_TESTBED_PROPAGATION_HEADERS}
                                                       ${CMAKE_CURRENT_SOURCE_DIR}/IntegratedAdaptiveStepPropagator.h
                                                       ${CMAKE_CURRENT_SOURCE_DIR}/IntegratedFixedStepPropagator.h
                                                       ${CMAKE_CURRENT_SOURCE_DIR}/IntegratedPropagator.h
                                                       ${CMAKE_CURRENT_SOURCE_DIR}/PropagatorBase.h
                                                       ${CMAKE_CURRENT_SOURCE_DIR}/PropagatorFactory.h)
                                                       
    set(ACCELERATION_MODEL_TESTBED_PROPAGATION_SOURCE ${ACCELERATION_MODEL_TESTBED_PROPAGATION_SOURCE}
                                                      ${CMAKE_CURRENT_SOURCE_DIR}/IntegratedAdaptiveStepPropagator.cpp
                                                      ${CMAKE_CURRENT_SOURCE_DIR}/IntegratedFixedStepPropagator.cpp
                                                      ${CMAKE_CURRENT_SOURCE_DIR}/IntegratedPropagator.cpp
                                                      ${CMAKE_CURRENT_SOURCE_DIR}/PropagatorBase.cpp
                                                      ${CMAKE_CURRENT_SOURCE_DIR}/PropagatorFactory.cpp)
                                                      
    set(ACCELERATION_MODEL_TESTBED_PROPAGATION_FILES ${ACCELERATION_MODEL_TESTBED_PROPAGATION_FILES}
                                                     ${ACCELERATION_MODEL_TESTBED_PROPAGATION_HEADERS}
                                                     ${ACCELERATION_MODEL_TESTBED_PROPAGATION_SOURCE} PARENT_SCOPE)
    
endif()
    
