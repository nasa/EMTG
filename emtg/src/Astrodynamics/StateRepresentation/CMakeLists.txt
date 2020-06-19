set(STATEREPRESENTATION_HEADERS  ${STATEREPRESENTATION_HEADERS} ${CMAKE_CURRENT_SOURCE_DIR}/StateRepresentation.h
                                                                ${CMAKE_CURRENT_SOURCE_DIR}/CartesianStateRepresentation.h
                                                                ${CMAKE_CURRENT_SOURCE_DIR}/SphericalRADECStateRepresentation.h
                                                                ${CMAKE_CURRENT_SOURCE_DIR}/SphericalAZFPAStateRepresentation.h
                                                                ${CMAKE_CURRENT_SOURCE_DIR}/COEStateRepresentation.h
                                                                ${CMAKE_CURRENT_SOURCE_DIR}/MEEStateRepresentation.h
                                                                ${CMAKE_CURRENT_SOURCE_DIR}/IncomingBplaneStateRepresentation.h
                                                                ${CMAKE_CURRENT_SOURCE_DIR}/OutgoingBplaneStateRepresentation.h
                                                                ${CMAKE_CURRENT_SOURCE_DIR}/StateRepresentationFactory.h
                                                                )
																	  
set(STATEREPRESENTATION_SOURCE  ${STATEREPRESENTATION_SOURCE} ${CMAKE_CURRENT_SOURCE_DIR}/StateRepresentation.cpp
                                                              ${CMAKE_CURRENT_SOURCE_DIR}/CartesianStateRepresentation.cpp
                                                              ${CMAKE_CURRENT_SOURCE_DIR}/SphericalRADECStateRepresentation.cpp
                                                              ${CMAKE_CURRENT_SOURCE_DIR}/SphericalAZFPAStateRepresentation.cpp
                                                              ${CMAKE_CURRENT_SOURCE_DIR}/COEStateRepresentation.cpp
                                                              ${CMAKE_CURRENT_SOURCE_DIR}/MEEStateRepresentation.cpp
                                                              ${CMAKE_CURRENT_SOURCE_DIR}/IncomingBplaneStateRepresentation.cpp
                                                              ${CMAKE_CURRENT_SOURCE_DIR}/OutgoingBplaneStateRepresentation.cpp
                                                              ${CMAKE_CURRENT_SOURCE_DIR}/StateRepresentationFactory.cpp
                                                              )
    
set(STATEREPRESENTATION_FILES ${STATEREPRESENTATION_FILES} ${STATEREPRESENTATION_HEADERS} ${STATEREPRESENTATION_SOURCE} PARENT_SCOPE)