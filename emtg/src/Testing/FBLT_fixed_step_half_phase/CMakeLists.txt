set(FBLT_fixed_step_half_phase_executable_files ${FBLT_fixed_step_half_phase_executable_files}
                                                ${CMAKE_CURRENT_SOURCE_DIR}/FBLT_fixed_step_half_phase.cpp)

source_group("Executable" FILES ${FBLT_fixed_step_half_phase_executable_files})

include_directories(../../Astrodynamics 
					../../Core 
					../../HardwareModels 
					../../Integration 
					../../Kepler 
					../../Math 
					../../Utilities)
					
if(SPLINE_EPHEM)
	include_directories(../../SplineEphem)
endif(SPLINE_EPHEM)

if(GSAD_INSTRUMENTATION)
	include_directories(${GSAD_PATH})
endif(GSAD_INSTRUMENTATION)

link_directories(${PROJECT_BINARY_DIR}/src/Testing/${CMAKE_CFG_INTDIR})

# create the executable
add_executable(FBLT_fixed_step_half_phase_test  ${FBLT_fixed_step_half_phase_executable_files})
target_link_libraries(FBLT_fixed_step_half_phase_test emtg_acceleration_model_library)