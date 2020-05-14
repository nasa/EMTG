
set(SCALATRON_HEADERS ${SCALATRON_HEADERS} ${CMAKE_CURRENT_SOURCE_DIR}/ScalatronBase.h
                                           ${CMAKE_CURRENT_SOURCE_DIR}/ScalatronPJRN.h
                                           ${CMAKE_CURRENT_SOURCE_DIR}/ScalatronFactory.h)

set(SCALATRON_SOURCE ${SCALATRON_SOURCE} ${CMAKE_CURRENT_SOURCE_DIR}/ScalatronBase.cpp
                                         ${CMAKE_CURRENT_SOURCE_DIR}/ScalatronPJRN.cpp
                                         ${CMAKE_CURRENT_SOURCE_DIR}/ScalatronFactory.cpp) 

set(SCALATRON_FILES ${SCALATRON_FILES} ${SCALATRON_HEADERS} ${SCALATRON_SOURCE} PARENT_SCOPE)