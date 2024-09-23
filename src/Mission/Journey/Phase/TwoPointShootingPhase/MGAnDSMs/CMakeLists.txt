add_subdirectory(ManeuverConstraints)

if(HAS_PROBEENTRYPHASE)
    add_subdirectory(ProbeEntryPhase)
endif()

set(MGANDSMS_HEADERS ${MGANDSMS_HEADERS} 
                     ${CMAKE_CURRENT_SOURCE_DIR}/MGAnDSMs_phase.h
                     ${CMAKE_CURRENT_SOURCE_DIR}/MGAnDSMs_subphase.h
                     ${CMAKE_CURRENT_SOURCE_DIR}/Forward_MGAnDSMs_subphase.h
                     ${CMAKE_CURRENT_SOURCE_DIR}/Backward_MGAnDSMs_subphase.h)
																	  
set(MGANDSMS_SOURCE ${MGANDSMS_SOURCE} 
                    ${CMAKE_CURRENT_SOURCE_DIR}/MGAnDSMs_phase.cpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/MGAnDSMs_subphase.cpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/Forward_MGAnDSMs_subphase.cpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/Backward_MGAnDSMs_subphase.cpp)

set(MGANDSMS_FILES ${MGANDSMS_FILES} 
                   ${MGANDSMS_HEADERS} 
                   ${MGANDSMS_SOURCE}
                   ${PROBE_ENTRY_PHASE_FILES} PARENT_SCOPE)
