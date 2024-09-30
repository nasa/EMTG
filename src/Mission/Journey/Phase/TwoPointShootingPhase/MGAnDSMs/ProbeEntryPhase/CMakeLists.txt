set(PROBE_ENTRY_PHASE_HEADERS ${PROBE_ENTRY_PHASE_HEADERS} 
                              ${CMAKE_CURRENT_SOURCE_DIR}/ProbeEntryPhase.h)
																	  
set(PROBE_ENTRY_PHASE_SOURCE ${PROBE_ENTRY_PHASE_SOURCE} 
                             ${CMAKE_CURRENT_SOURCE_DIR}/ProbeEntryPhase.cpp)

set(PROBE_ENTRY_PHASE_FILES ${PROBE_ENTRY_PHASE_FILES} 
                            ${PROBE_ENTRY_PHASE_HEADERS} 
                            ${PROBE_ENTRY_PHASE_SOURCE} PARENT_SCOPE)
