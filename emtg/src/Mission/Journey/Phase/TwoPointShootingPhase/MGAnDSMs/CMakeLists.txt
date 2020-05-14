add_subdirectory(ManeuverConstraints)

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
                   ${MGANDSMS_SOURCE} PARENT_SCOPE)
