
set(EOM_HEADERS  ${EOM_HEADERS} ${CMAKE_CURRENT_SOURCE_DIR}/SundmanSpacecraftEOM.h
                                ${CMAKE_CURRENT_SOURCE_DIR}/TimeDomainSpacecraftEOM.h
                                ${CMAKE_CURRENT_SOURCE_DIR}/EOMFactory.h)

set(EOM_SOURCE  ${EOM_SOURCE} ${CMAKE_CURRENT_SOURCE_DIR}/SundmanSpacecraftEOM.cpp
                              ${CMAKE_CURRENT_SOURCE_DIR}/TimeDomainSpacecraftEOM.cpp
                              ${CMAKE_CURRENT_SOURCE_DIR}/EOMFactory.cpp)
    
set(EOM_FILES ${EOM_FILES} ${EOM_HEADERS} ${EOM_SOURCE} PARENT_SCOPE)

if (RUN_ACCELERATION_MODEL_TESTBED)
    
    set(ACCELERATION_MODEL_TESTBED_EOM_HEADERS ${ACCELERATION_MODEL_TESTBED_EOM_HEADERS}
                                               ${CMAKE_CURRENT_SOURCE_DIR}/SundmanSpacecraftEOM.h
                                               ${CMAKE_CURRENT_SOURCE_DIR}/TimeDomainSpacecraftEOM.h
                                               ${CMAKE_CURRENT_SOURCE_DIR}/EOMFactory.h)
                                               
                                               
                                                       
    set(ACCELERATION_MODEL_TESTBED_EOM_SOURCE ${ACCELERATION_MODEL_TESTBED_EOM_SOURCE}
                                              ${CMAKE_CURRENT_SOURCE_DIR}/SundmanSpacecraftEOM.cpp
                                              ${CMAKE_CURRENT_SOURCE_DIR}/TimeDomainSpacecraftEOM.cpp
                                              ${CMAKE_CURRENT_SOURCE_DIR}/EOMFactory.cpp)
                                              
                                              
                                                      
    set(ACCELERATION_MODEL_TESTBED_EOM_FILES ${ACCELERATION_MODEL_TESTBED_EOM_FILES}
                                             ${ACCELERATION_MODEL_TESTBED_EOM_HEADERS}
                                             ${ACCELERATION_MODEL_TESTBED_EOM_SOURCE} PARENT_SCOPE)
    
endif()


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