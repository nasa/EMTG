
set(INNER_LOOP_HEADERS ${INNER_LOOP_HEADERS} ${CMAKE_CURRENT_SOURCE_DIR}/EMTG_innerloop_solution.h
                                             ${CMAKE_CURRENT_SOURCE_DIR}/EMTG_WORHP_interface.h
                                             ${CMAKE_CURRENT_SOURCE_DIR}/FilamentWalker.h
                                             ${CMAKE_CURRENT_SOURCE_DIR}/monotonic_basin_hopping.h
                                             ${CMAKE_CURRENT_SOURCE_DIR}/NLPoptions.h
                                             ${CMAKE_CURRENT_SOURCE_DIR}/NLP_interface.h
                                             ${CMAKE_CURRENT_SOURCE_DIR}/SNOPT_interface.h
                                             ${CMAKE_CURRENT_SOURCE_DIR}/snoptProblemExtension.h)
												  
set(INNER_LOOP_SOURCE ${INNER_LOOP_SOURCE} ${CMAKE_CURRENT_SOURCE_DIR}/EMTG_innerloop_solution.cpp
								           ${CMAKE_CURRENT_SOURCE_DIR}/EMTG_WORHP_interface.cpp
                                           ${CMAKE_CURRENT_SOURCE_DIR}/FilamentWalker.cpp
                                           ${CMAKE_CURRENT_SOURCE_DIR}/monotonic_basin_hopping.cpp
                                           ${CMAKE_CURRENT_SOURCE_DIR}/NLPoptions.cpp
                                           ${CMAKE_CURRENT_SOURCE_DIR}/NLP_interface.cpp
                                           ${CMAKE_CURRENT_SOURCE_DIR}/SNOPT_interface.cpp)
                                           
set(INNER_LOOP_FILES ${INNER_LOOP_FILES} ${INNER_LOOP_HEADERS} ${INNER_LOOP_SOURCE} PARENT_SCOPE)
    

    
