set(PSBI_PHASE_HEADERS ${PSBI_PHASE_HEADERS} 
                       ${CMAKE_CURRENT_SOURCE_DIR}/PSBIphase.h
                       ${CMAKE_CURRENT_SOURCE_DIR}/PSBIstep.h
                       ${CMAKE_CURRENT_SOURCE_DIR}/PSBIfirststep.h
                       ${CMAKE_CURRENT_SOURCE_DIR}/PSBIlaststep.h
                       ${CMAKE_CURRENT_SOURCE_DIR}/PSBIOneStepToRuleThemAll.h
                       ${CMAKE_CURRENT_SOURCE_DIR}/PSBIstep_factory.h)
																	  
set(PSBI_PHASE_SOURCE ${PSBI_PHASE_SOURCE} 
                      ${CMAKE_CURRENT_SOURCE_DIR}/PSBIphase.cpp
                      ${CMAKE_CURRENT_SOURCE_DIR}/PSBIstep.cpp
                      ${CMAKE_CURRENT_SOURCE_DIR}/PSBIfirststep.cpp
                      ${CMAKE_CURRENT_SOURCE_DIR}/PSBIlaststep.cpp
                      ${CMAKE_CURRENT_SOURCE_DIR}/PSBIOneStepToRuleThemAll.cpp
                      ${CMAKE_CURRENT_SOURCE_DIR}/PSBIstep_factory.cpp)

set(PSBI_PHASE_FILES ${PSBI_PHASE_FILES} 
                     ${PSBI_PHASE_HEADERS} 
                     ${PSBI_PHASE_SOURCE} PARENT_SCOPE)
