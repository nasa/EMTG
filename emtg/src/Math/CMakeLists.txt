# main EMTG library files
set(MATH_HEADERS ${MATH_HEADERS} ${CMAKE_CURRENT_SOURCE_DIR}/EMTG_math.h
                                 ${CMAKE_CURRENT_SOURCE_DIR}/EMTG_Matrix.h
                                 ${CMAKE_CURRENT_SOURCE_DIR}/EMTG_Tensor.h
                                 ${CMAKE_CURRENT_SOURCE_DIR}/randutils.hpp
								 ${CMAKE_CURRENT_SOURCE_DIR}/interpolator.h)
								 
set(MATH_SOURCE ${MATH_SOURCE} ${CMAKE_CURRENT_SOURCE_DIR}/interpolator.cpp
                               ${CMAKE_CURRENT_SOURCE_DIR}/EMTG_Matrix.cpp
                               ${CMAKE_CURRENT_SOURCE_DIR}/EMTG_Tensor.cpp)
								 
set(MATH_FILES ${MATH_FILES} ${MATH_HEADERS} ${MATH_SOURCE} PARENT_SCOPE)

# FBLT forward half-phase fixed-step files
if (RUN_FBLT_TESTBED)
    set(FBLT_TESTBED_MATH_HEADERS ${FBLT_TESTBED_MATH_HEADERS} 
								  ${CMAKE_CURRENT_SOURCE_DIR}/EMTG_math.h
                                  ${CMAKE_CURRENT_SOURCE_DIR}/EMTG_Matrix.h
                                  ${CMAKE_CURRENT_SOURCE_DIR}/EMTG_Tensor.h)
                                  
    set(FBLT_TESTBED_MATH_SOURCE  ${FBLT_TESTBED_MATH_SOURCE} 
                                  ${CMAKE_CURRENT_SOURCE_DIR}/EMTG_Matrix.cpp
                                  ${CMAKE_CURRENT_SOURCE_DIR}/EMTG_Tensor.cpp)
                                                        
    set(FBLT_TESTBED_MATH_FILES ${FBLT_TESTBED_MATH_FILES} 
                                ${FBLT_TESTBED_MATH_SOURCE} 
                                ${FBLT_TESTBED_MATH_HEADERS} PARENT_SCOPE)
endif()
    
# acceleration model unit test files
if (RUN_ACCELERATION_MODEL_TESTBED)
    set(ACCELERATION_MODEL_TESTBED_MATH_HEADERS ${ACCELERATION_MODEL_TESTBED_MATH_HEADERS} 
									            ${CMAKE_CURRENT_SOURCE_DIR}/EMTG_math.h
                                                ${CMAKE_CURRENT_SOURCE_DIR}/EMTG_Matrix.h
                                                ${CMAKE_CURRENT_SOURCE_DIR}/EMTG_Tensor.h)
												  
	set(ACCELERATION_MODEL_TESTBED_MATH_SOURCE ${ACCELERATION_MODEL_TESTBED_MATH_SOURCE} 
								               ${CMAKE_CURRENT_SOURCE_DIR}/EMTG_Matrix.cpp
                                               ${CMAKE_CURRENT_SOURCE_DIR}/EMTG_Tensor.cpp)
                                                        
    set(ACCELERATION_MODEL_TESTBED_MATH_FILES ${ACCELERATION_MODEL_TESTBED_MATH_FILES} 
                                              ${ACCELERATION_MODEL_TESTBED_MATH_HEADERS}
								              ${ACCELERATION_MODEL_TESTBED_MATH_SOURCE} PARENT_SCOPE)
endif()