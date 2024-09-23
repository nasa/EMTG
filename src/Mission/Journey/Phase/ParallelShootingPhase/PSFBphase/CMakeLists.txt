set(PSFB_PHASE_HEADERS ${PSFB_PHASE_HEADERS} 
                       ${CMAKE_CURRENT_SOURCE_DIR}/PSFBphase.h
                       ${CMAKE_CURRENT_SOURCE_DIR}/PSFBstep.h
                       ${CMAKE_CURRENT_SOURCE_DIR}/PSFBfirststep.h
                       ${CMAKE_CURRENT_SOURCE_DIR}/PSFBlaststep.h
                       ${CMAKE_CURRENT_SOURCE_DIR}/PSFBOneStepToRuleThemAll.h
                       ${CMAKE_CURRENT_SOURCE_DIR}/PSFB_HifiDuty_step.h
                       ${CMAKE_CURRENT_SOURCE_DIR}/PSFB_HifiDuty_firststep.h
                       ${CMAKE_CURRENT_SOURCE_DIR}/PSFB_HifiDuty_laststep.h
                       ${CMAKE_CURRENT_SOURCE_DIR}/PSFB_HifiDuty_OneStepToRuleThemAll.h
                       ${CMAKE_CURRENT_SOURCE_DIR}/PSFBstep_factory.h)
																	  
set(PSFB_PHASE_SOURCE ${PSFB_PHASE_SOURCE} 
                      ${CMAKE_CURRENT_SOURCE_DIR}/PSFBphase.cpp
                      ${CMAKE_CURRENT_SOURCE_DIR}/PSFBstep.cpp
                      ${CMAKE_CURRENT_SOURCE_DIR}/PSFBfirststep.cpp
                      ${CMAKE_CURRENT_SOURCE_DIR}/PSFBlaststep.cpp
                      ${CMAKE_CURRENT_SOURCE_DIR}/PSFBOneStepToRuleThemAll.cpp
                      ${CMAKE_CURRENT_SOURCE_DIR}/PSFB_HifiDuty_step.cpp
                      ${CMAKE_CURRENT_SOURCE_DIR}/PSFB_HifiDuty_firststep.cpp
                      ${CMAKE_CURRENT_SOURCE_DIR}/PSFB_HifiDuty_laststep.cpp
                      ${CMAKE_CURRENT_SOURCE_DIR}/PSFB_HifiDuty_OneStepToRuleThemAll.cpp
                      ${CMAKE_CURRENT_SOURCE_DIR}/PSFBstep_factory.cpp)

set(PSFB_PHASE_FILES ${PSFB_PHASE_FILES} 
                     ${PSFB_PHASE_HEADERS} 
                     ${PSFB_PHASE_SOURCE} PARENT_SCOPE)
