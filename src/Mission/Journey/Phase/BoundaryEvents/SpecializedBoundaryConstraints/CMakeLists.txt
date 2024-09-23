add_subdirectory(OrbitElementsConstraints)
add_subdirectory(BoundaryAngleConstraints)

set(SPECIALIZED_BOUNDARY_CONSTRAINTS_HEADERS ${SPECIALIZED_BOUNDARY_CONSTRAINTS_HEADERS} 
                                             ${CMAKE_CURRENT_SOURCE_DIR}/SpecializedBoundaryConstraintBase.h
                                             ${CMAKE_CURRENT_SOURCE_DIR}/BoundaryDistanceConstraint.h
                                             ${CMAKE_CURRENT_SOURCE_DIR}/BoundaryIsApseConstraint.h
                                             ${CMAKE_CURRENT_SOURCE_DIR}/BoundaryMassConstraint.h
                                             ${CMAKE_CURRENT_SOURCE_DIR}/BoundaryLongitudeConstraint.h
                                             ${CMAKE_CURRENT_SOURCE_DIR}/BoundaryBCFlatitudeConstraint.h
                                             ${CMAKE_CURRENT_SOURCE_DIR}/BoundaryDeticLatitudeConstraint.h
                                             ${CMAKE_CURRENT_SOURCE_DIR}/BoundaryDeticAltitudeConstraint.h
                                             ${CMAKE_CURRENT_SOURCE_DIR}/BoundaryOrbitalEnergyConstraint.h
                                             ${CMAKE_CURRENT_SOURCE_DIR}/BoundaryVelocityDeclinationConstraint.h
                                             ${CMAKE_CURRENT_SOURCE_DIR}/BoundaryVelocityMagnitudeConstraint.h
                                             ${CMAKE_CURRENT_SOURCE_DIR}/BoundaryVelocitySphericalAzimuthConstraint.h
                                             ${CMAKE_CURRENT_SOURCE_DIR}/BoundaryRelativeVelocityMagnitudeConstraint.h
                                             ${CMAKE_CURRENT_SOURCE_DIR}/BoundaryRelativeVelocityAzimuthConstraint.h
                                             ${CMAKE_CURRENT_SOURCE_DIR}/BoundaryRelativeVelocityHFPAConstraint.h
                                             ${CMAKE_CURRENT_SOURCE_DIR}/BoundaryTargetBodyDeticElevationConstraint.h
                                             ${CMAKE_CURRENT_SOURCE_DIR}/BoundaryVelocityVFPAconstraint.h
                                             ${CMAKE_CURRENT_SOURCE_DIR}/SpecializedBoundaryConstraintFactory.h)
																	  
set(SPECIALIZED_BOUNDARY_CONSTRAINTS_SOURCE ${SPECIALIZED_BOUNDARY_CONSTRAINTS_SOURCE} 
                                            ${CMAKE_CURRENT_SOURCE_DIR}/SpecializedBoundaryConstraintBase.cpp
                                            ${CMAKE_CURRENT_SOURCE_DIR}/BoundaryDistanceConstraint.cpp
                                            ${CMAKE_CURRENT_SOURCE_DIR}/BoundaryIsApseConstraint.cpp
                                            ${CMAKE_CURRENT_SOURCE_DIR}/BoundaryMassConstraint.cpp
                                            ${CMAKE_CURRENT_SOURCE_DIR}/BoundaryLongitudeConstraint.cpp
                                            ${CMAKE_CURRENT_SOURCE_DIR}/BoundaryBCFlatitudeConstraint.cpp
                                            ${CMAKE_CURRENT_SOURCE_DIR}/BoundaryDeticLatitudeConstraint.cpp
                                            ${CMAKE_CURRENT_SOURCE_DIR}/BoundaryDeticAltitudeConstraint.cpp
                                            ${CMAKE_CURRENT_SOURCE_DIR}/BoundaryOrbitalEnergyConstraint.cpp
                                            ${CMAKE_CURRENT_SOURCE_DIR}/BoundaryVelocityDeclinationConstraint.cpp
                                            ${CMAKE_CURRENT_SOURCE_DIR}/BoundaryVelocityMagnitudeConstraint.cpp
                                            ${CMAKE_CURRENT_SOURCE_DIR}/BoundaryVelocitySphericalAzimuthConstraint.cpp
                                            ${CMAKE_CURRENT_SOURCE_DIR}/BoundaryRelativeVelocityMagnitudeConstraint.cpp
                                            ${CMAKE_CURRENT_SOURCE_DIR}/BoundaryRelativeVelocityAzimuthConstraint.cpp
                                            ${CMAKE_CURRENT_SOURCE_DIR}/BoundaryRelativeVelocityHFPAConstraint.cpp
                                            ${CMAKE_CURRENT_SOURCE_DIR}/BoundaryTargetBodyDeticElevationConstraint.cpp
                                            ${CMAKE_CURRENT_SOURCE_DIR}/BoundaryVelocityVFPAconstraint.cpp
                                            ${CMAKE_CURRENT_SOURCE_DIR}/SpecializedBoundaryConstraintFactory.cpp)

set(SPECIALIZED_BOUNDARY_CONSTRAINTS_FILES ${SPECIALIZED_BOUNDARY_CONSTRAINTS_FILES} 
                                           ${SPECIALIZED_BOUNDARY_CONSTRAINTS_HEADERS} 
                                           ${SPECIALIZED_BOUNDARY_CONSTRAINTS_SOURCE} 
                                           ${SPECIALIZED_EPHEMERIS_PEGGED_BOUNDARY_CONSTRAINTS_FILES} 
                                           ${ORBIT_ELEMENTS_CONSTRAINTS_FILES} 
                                           ${BOUNDARY_ANGLE_CONSTRAINTS_FILES} PARENT_SCOPE)
