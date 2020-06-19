set(missionoptions_test_executable_files   ${missionoptions_test_executable_files}
                                    ${CMAKE_CURRENT_SOURCE_DIR}/missionoptions_testbed_driver.cpp)

set(missionoptions_test_source      ${missionoptions_test_source}
                                    ${CMAKE_CURRENT_SOURCE_DIR}/../../Utilities/file_utilities.h
                                    ${CMAKE_CURRENT_SOURCE_DIR}/../../Utilities/file_utilities.cpp
                                    ${CMAKE_CURRENT_SOURCE_DIR}/../../core/missionoptions.h
                                    ${CMAKE_CURRENT_SOURCE_DIR}/../../core/missionoptions.cpp
                                    ${CMAKE_CURRENT_SOURCE_DIR}/../../core/journeyoptions.h
                                    ${CMAKE_CURRENT_SOURCE_DIR}/../../core/journeyoptions.cpp)

source_group("Executable" FILES ${missionoptions_test_executable_files} ${missionoptions_test_source})


link_directories(${PROJECT_BINARY_DIR}/src/Testing/${CMAKE_CFG_INTDIR})

# create the executable
if(RUN_MISSIONOPTIONS_TESTBED)
    add_executable(missionoptions_testbed  ${missionoptions_test_executable_files} ${missionoptions_test_source})
endif()

#create the boost::python module
if(BUILD_PyMISSIONOPTIONS)
    source_group("PyMissionOptions" FILES ${missionoptions_test_source})
    
    SET(Python_ADDITIONAL_VERSIONS ${PYTHON_VERSIONS_REQUIRED}$ )
        
    FIND_PACKAGE(PythonInterp)

    if (PYTHONINTERP_FOUND)
      if (UNIX AND NOT APPLE)
        if (PYTHON_VERSION_MAJOR EQUAL 3)
            FIND_PACKAGE(Boost COMPONENTS python${PYTHON_VERSION_SUFFIX})
            FIND_PACKAGE(PythonInterp 3)
            FIND_PACKAGE(PythonLibs 3 REQUIRED)
        else()
            FIND_PACKAGE(Boost COMPONENTS python)
            FIND_PACKAGE(PythonInterp)
            FIND_PACKAGE(PythonLibs REQUIRED)
        endif()
      else()	
        if (PYTHON_VERSION_MAJOR EQUAL 3)
            FIND_PACKAGE(Boost COMPONENTS python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR})
            FIND_PACKAGE(PythonInterp 3)
            FIND_PACKAGE(PythonLibs 3 REQUIRED)
        else()
            FIND_PACKAGE(Boost COMPONENTS python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR})
            FIND_PACKAGE(PythonInterp)
            FIND_PACKAGE(PythonLibs REQUIRED)
        endif()
      endif()
    else()
        message("Python not found")
    endif()

    message(STATUS "PYTHON_LIBRARIES = ${PYTHON_LIBRARIES}")
    message(STATUS "PYTHON_EXECUTABLE = ${PYTHON_EXECUTABLE}")
    message(STATUS "PYTHON_INCLUDE_DIRS = ${PYTHON_INCLUDE_DIRS}")
    message(STATUS "Boost_LIBRARIES = ${Boost_LIBRARIES}")

    add_library(PyMissionOptions SHARED ${missionoptions_test_source})
    target_include_directories(PyMissionOptions PRIVATE ${PYTHON_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})
    #turn on the flag to compile the python interface, and the flag to disable boost auto-link
    target_compile_definitions(PyMissionOptions PRIVATE EMTG_OPTIONS_PYTHON_INTERFACE BOOST_ALL_NO_LIB BOOST_PYTHON_STATIC_LIB)
    target_link_libraries(PyMissionOptions ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})

    add_custom_command(TARGET PyMissionOptions POST_BUILD
                               COMMAND ${CMAKE_COMMAND} -E copy
                                   "${CMAKE_CFG_INTDIR}\\PyMissionOptions.dll"
                                   "${OUTPUTDESTINATION}\\PyMissionOptions.pyd")
endif(BUILD_PyMISSIONOPTIONS)