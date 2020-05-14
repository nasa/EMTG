set(MGANDSMS_HEADERS ${MGANDSMS_HEADERS} 
                     ${CMAKE_CURRENT_SOURCE_DIR}/MGAnDSMs_maneuver_constraint.h
                     ${CMAKE_CURRENT_SOURCE_DIR}/MGAnDSMs_maneuver_epoch_constraint.h
                     ${CMAKE_CURRENT_SOURCE_DIR}/MGAnDSMs_maneuver_magnitude_constraint.h
                     ${CMAKE_CURRENT_SOURCE_DIR}/MGAnDSMs_maneuver_constraint_factory.h
                     ${CMAKE_CURRENT_SOURCE_DIR}/MGAnDSMs_backward_maneuver_epoch_constraint_absolute.h
                     ${CMAKE_CURRENT_SOURCE_DIR}/MGAnDSMs_backward_maneuver_epoch_constraint_wrt_left_boundary.h
					 ${CMAKE_CURRENT_SOURCE_DIR}/MGAnDSMs_backward_maneuver_epoch_constraint_wrt_next_event.h
                     ${CMAKE_CURRENT_SOURCE_DIR}/MGAnDSMs_backward_maneuver_epoch_constraint_wrt_previous_event.h
                     ${CMAKE_CURRENT_SOURCE_DIR}/MGAnDSMs_backward_maneuver_epoch_constraint_wrt_right_boundary.h
                     ${CMAKE_CURRENT_SOURCE_DIR}/MGAnDSMs_forward_maneuver_epoch_constraint_absolute.h
                     ${CMAKE_CURRENT_SOURCE_DIR}/MGAnDSMs_forward_maneuver_epoch_constraint_wrt_left_boundary.h
					 ${CMAKE_CURRENT_SOURCE_DIR}/MGAnDSMs_forward_maneuver_epoch_constraint_wrt_next_event.h
                     ${CMAKE_CURRENT_SOURCE_DIR}/MGAnDSMs_forward_maneuver_epoch_constraint_wrt_previous_event.h
                     ${CMAKE_CURRENT_SOURCE_DIR}/MGAnDSMs_forward_maneuver_epoch_constraint_wrt_right_boundary.h
                     ${CMAKE_CURRENT_SOURCE_DIR}/MGAnDSMs_maneuver_cone_angle_exclusion_constraint.h)
                     
set(MGANDSMS_SOURCE ${MGANDSMS_SOURCE} 
                     ${CMAKE_CURRENT_SOURCE_DIR}/MGAnDSMs_maneuver_constraint.cpp
                     ${CMAKE_CURRENT_SOURCE_DIR}/MGAnDSMs_maneuver_epoch_constraint.cpp
                     ${CMAKE_CURRENT_SOURCE_DIR}/MGAnDSMs_maneuver_magnitude_constraint.cpp
                     ${CMAKE_CURRENT_SOURCE_DIR}/MGAnDSMs_maneuver_constraint_factory.cpp
                     ${CMAKE_CURRENT_SOURCE_DIR}/MGAnDSMs_backward_maneuver_epoch_constraint_absolute.cpp
                     ${CMAKE_CURRENT_SOURCE_DIR}/MGAnDSMs_backward_maneuver_epoch_constraint_wrt_left_boundary.cpp
					 ${CMAKE_CURRENT_SOURCE_DIR}/MGAnDSMs_backward_maneuver_epoch_constraint_wrt_next_event.cpp
                     ${CMAKE_CURRENT_SOURCE_DIR}/MGAnDSMs_backward_maneuver_epoch_constraint_wrt_previous_event.cpp
                     ${CMAKE_CURRENT_SOURCE_DIR}/MGAnDSMs_backward_maneuver_epoch_constraint_wrt_right_boundary.cpp
                     ${CMAKE_CURRENT_SOURCE_DIR}/MGAnDSMs_forward_maneuver_epoch_constraint_absolute.cpp
                     ${CMAKE_CURRENT_SOURCE_DIR}/MGAnDSMs_forward_maneuver_epoch_constraint_wrt_left_boundary.cpp
					 ${CMAKE_CURRENT_SOURCE_DIR}/MGAnDSMs_forward_maneuver_epoch_constraint_wrt_next_event.cpp
                     ${CMAKE_CURRENT_SOURCE_DIR}/MGAnDSMs_forward_maneuver_epoch_constraint_wrt_previous_event.cpp
                     ${CMAKE_CURRENT_SOURCE_DIR}/MGAnDSMs_forward_maneuver_epoch_constraint_wrt_right_boundary.cpp
                     ${CMAKE_CURRENT_SOURCE_DIR}/MGAnDSMs_maneuver_cone_angle_exclusion_constraint.cpp)
                     
set(MGANDSMS_FILES ${MGANDSMS_FILES} 
                   ${MGANDSMS_HEADERS} 
                   ${MGANDSMS_SOURCE} PARENT_SCOPE)
