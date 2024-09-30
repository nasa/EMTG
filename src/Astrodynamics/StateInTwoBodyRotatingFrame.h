// EMTG: Evolutionary Mission Trajectory Generator
// An open-source global optimization tool for preliminary mission design
// Provided by NASA Goddard Space Flight Center
//
// Copyright (c) 2013 - 2024 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// All Other Rights Reserved.

// Licensed under the NASA Open Source License (the "License");
// You may not use this file except in compliance with the License.
// You may obtain a copy of the License at:
// https://opensource.org/licenses/NASA-1.3
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
// express or implied.   See the License for the specific language
// governing permissions and limitations under the License.

//header file for StateInTwoBodyRotatingFrame

#pragma once

#include "EMTG_math.h"
#include "EMTG_Matrix.h"
#include "EMTG_Tensor.h"

#include<vector>
#include<iomanip>

//#include "doubleType.h"

//#include <string>
//#include <vector>
//#include <fstream>
//#include <iostream>
//
//#include "missionoptions.h"
//#include "CentralBody.h"
//#include "body.h"
//#include "frame.h"
//#include "EMTG_math.h"
//#include "atmosphere.h"
//
//#include "SpiceUsr.h"
//
//#include "boost/ptr_container/ptr_vector.hpp"
//#include "boost/algorithm/string.hpp"
//
//#ifdef SPLINE_EPHEM
//#include "SplineEphem_universe.h"
//#endif
//
//#ifdef AD_INSTRUMENTATION
//#include <set>
//#endif

//#define DEBUG_ROTFRAME_B2WrtB1AngularVelocityVector
//#define DEBUG_ROTFRAME_TwoBodyRotatingFrameUnitVectors
//#define DEBUG_ROTFRAME_CalculateStateInRotatingFrameWrtB2

namespace EMTG
{
    namespace Astrodynamics
    {
		namespace TwoBodyRotatingFrame
		{
			//------------------------------------------------------------------------------
			// B2WrtB1AngularVelocityVector(const math::Matrix<T>& stateB1WrtB2Inertial, const math::Matrix<T>& dStateB2WrtB1InertialDt, 
			// math::Matrix<T>& omega, const bool& generateDerivatives, math::Matrix<T>& dOmegaDt)
			//------------------------------------------------------------------------------
			/**
			 * Calculate instantaneous angular velocity of body B2 with respect to
			 * another body B1, given the inertial position, velocity, and acceleration
			 * vectors of B2 with respect to B1. Optionally return time derivative of the
			 * angular velocity vector. Assumes the state of B2 with respect to B1 --
			 * and thus the angular velocity vector -- is a function of time only.
			 *
			 * @param stateB2WrtB1Inertial 6x1 EMTG matrix. Contains position and velocity vectors of B2 with respect to B1 in an inertial frame.
			 * @param dstateB2WrtB1Inertial 6x1 EMTG matrix. Contains time derivatives of stateB2WrtB1Inertial. Only used if generateDerivatives it true.
			 * @return omega 3x1 EMTG matrix. Instantaneous angular velocity vector of B2 about B1.
			 * @param generateDerivatives Boolean. If true, generate (time) derivatives of angular velocity vector. If false, don't.
			 * @return dOmegaDt 3x1 EMTG matrix. Time derivative of angular velocity vector: d [omega] / d [t]
			 */
			 //------------------------------------------------------------------------------
			template<class T> void B2WrtB1AngularVelocityVector(const math::Matrix<T>& stateB2WrtB1Inertial, const math::Matrix<T>& dStateB2WrtB1InertialDt,
				math::Matrix<T>& omega, const bool& generateDerivatives, math::Matrix<T>& dOmegaDt)
			{
				// unpack state vector
				math::Matrix<T> r(3, 1), v(3, 1), drdt(3, 1), dvdt(3, 1);
				r = stateB2WrtB1Inertial.getSubMatrix1D(0, 2);
				T rMag = r.norm();
				T oneByRSquared = 1.0 / (rMag * rMag);
				v = stateB2WrtB1Inertial.getSubMatrix1D(3, 5);
				drdt = dStateB2WrtB1InertialDt.getSubMatrix1D(0, 2);
				dvdt = dStateB2WrtB1InertialDt.getSubMatrix1D(3, 5);

				// specific angular momentum vector
				math::Matrix<T> h(3, 1);
				h = r.cross(v);

				// angular velocity vector itself
				omega = h * oneByRSquared;

				if (generateDerivatives)
				{
					// build up with chain rule
					T dOneByRSquaredDt = (-2.0 * oneByRSquared * oneByRSquared) * r.dot(drdt);

					math::Matrix<T> dhdt(3, 1), dhdr(3, 3), dhdv(3, 3), dhdx(6, 3);
					dhdx = r.crossDerivative(v); // contains dh/dr and dh/dv
					for (size_t i = 0; i < 3; ++i)
					{
						for (size_t j = 0; j < 3; ++j)
						{
							dhdr(j, i) = dhdx(i, j); // note index flipping! emtg math organizes things in a way that may be counterintuitive
							dhdv(j, i) = dhdx(i + 3, j); // note index flipping! emtg math organizes things in a way that may be counterintuitive
						}
					}
					dhdt = dhdr * drdt + dhdv * dvdt;
					dOmegaDt = h * dOneByRSquaredDt + dhdt * oneByRSquared;

#ifdef DEBUG_ROTFRAME_B2WrtB1AngularVelocityVector
					std::cout << std::setprecision(16);
					/*std::cout << "dh/dr analytical:\n\n";
					for (size_t i = 0; i < 3; ++i)
					{
						for (size_t j = 0; j < 3; ++j)
						{
							std::cout << "(" << i << ", " << j << "): " << dhdr(i, j) << "\n";
						}
					}
					std::cout << "dh/dv analytical:\n\n";
					for (size_t i = 0; i < 3; ++i)
					{
						for (size_t j = 0; j < 3; ++j)
						{
							std::cout << "(" << i << ", " << j << "): " << dhdv(i, j) << "\n";
						}
					}
					std::cout << "dh/dr AD:\n\n";
					for (size_t i = 0; i < 3; ++i)
					{
						for (size_t j = 0; j < 3; ++j)
						{
							std::cout << "(" << i << ", " << j << "): " << h(i).getDerivative(j + 1) << "\n";
						}
					}
					std::cout << "dh/dv AD:\n\n";
					for (size_t i = 0; i < 3; ++i)
					{
						for (size_t j = 0; j < 3; ++j)
						{
							std::cout << "(" << i << ", " << j << "): " << h(i).getDerivative(j + 4) << "\n";
						}
					}*/
					std::cout << "\n\nB2WrtB1AngularVelocityVector derivatives checking:\n\n";
					std::cout << "dh/dr analytical - AD:\n";
					for (size_t i = 0; i < 3; ++i)
					{
						for (size_t j = 0; j < 3; ++j)
						{
							std::cout << "(" << i << ", " << j << "): " << dhdr(i, j) - h(i).getDerivative(j + 1) << "\n";
						}
					}
					std::cout << "\ndh/dv analytical - AD:\n";
					for (size_t i = 0; i < 3; ++i)
					{
						for (size_t j = 0; j < 3; ++j)
						{
							std::cout << "(" << i << ", " << j << "): " << dhdv(i, j) - h(i).getDerivative(j + 4) << "\n";
						}
					}
					std::cout << "\ndh/dt analytical - AD:\n";
					for (size_t i = 0; i < 3; ++i)
					{
						std::cout << "(" << i << "): " << dhdt(i) - h(i).getDerivative(0) << "\n";
					}
					std::cout << "\nd (1/r^2)/dt analytical - AD:\n";
					std::cout << dOneByRSquaredDt - oneByRSquared.getDerivative(0) << "\n";
#endif
				}

				return;
			}

	

			//------------------------------------------------------------------------------
			// TwoBodyRotatingFrameUnitVectors(const math::Matrix<T>& stateB2WrtB1Inertial, const math::Matrix<T>& dStateB2WrtB1InertialDt,
			// math::Matrix<T>& unitVectors, const bool& generateDerivatives, math::Matrix<T>& dUnitVectorsDt)
			//------------------------------------------------------------------------------
			/**
			 * Calculate the units vectors of a reference frame defined by the
			 * instantaneous angular velocity of body B2 with respect to
			 * another body B1, given the inertial position, velocity, and acceleration
			 * vectors of B2 with respect to B1. Optionally return time derivatives of the
			 * unit vectors. Assumes the state of B2 with respect to B1 --
			 * and thus the unit vectors of the reference frame -- are functions of time only.
			 *
			 * The frame is defined such that the kR axis is in the direction of the angular momentum
			 * vector of B2 with respect to B1. The iR axis is directed along the position vector from
			 * B1 to B2. The jR axis completes the right-handed set (kR \times iR).
			 *
			 * @param stateB2WrtB1Inertial 6x1 EMTG matrix. Contains position and velocity vectors of B2 with respect to B1 in an inertial frame.
			 * @param dStateB2WrtB1InertialDt 6x1 EMTG matrix. Contains time derivative of stateB2WrtB1Inertial. Only used if generateDerivatives it true.
			 * @return unitVectors 3x3 EMTG matrix. Instantaneous unit vectors of the reference frame: [iR, jR, kR]. Note that this is equivalent to the
			 * transpose of the transformation matrix used to relate the inertial frame to the rotating frame:
			 * vector_{rotating frame} = unitVectors^T * vector_{inertial frame}
			 * @param generateDerivatives Boolean. If true, generate (time) derivatives of unit vectors. If false, don't.
			 * @return dUnitVectorsDt 3x3 EMTG matrix. Time derivative of unit vectors: d [unitVectors] / d [t]
			 */
			 //------------------------------------------------------------------------------
			template<class T> void TwoBodyRotatingFrameUnitVectors(const math::Matrix<T>& stateB2WrtB1Inertial, const math::Matrix<T>& dStateB2WrtB1InertialDt,
				math::Matrix<T>& unitVectors, const bool& generateDerivatives, math::Matrix<T>& dUnitVectorsDt)
			{
				// get omega
				math::Matrix<T> omega(3, 1, 0.0);

				math::Matrix<T> dOmegaDt(3, 1);
				B2WrtB1AngularVelocityVector(stateB2WrtB1Inertial, dStateB2WrtB1InertialDt,
					omega, generateDerivatives, dOmegaDt);

			#ifdef DEBUG_ROTFRAME_TwoBodyRotatingFrameUnitVectors
				size_t adIndex = 7;
				omega(0).setDerivative(adIndex++, 1.0);
				omega(1).setDerivative(adIndex++, 1.0);
				omega(2).setDerivative(adIndex++, 1.0);
			#endif

				// construct the unit vectors
				math::Matrix<T> iR(3, 1), jR(3, 1), kR(3, 1), r(3, 1), v(3, 1), drdt(3, 1), dvdt(3, 1);
				r = stateB2WrtB1Inertial.getSubMatrix1D(0, 2);
				v = stateB2WrtB1Inertial.getSubMatrix1D(3, 5);
				
				iR = r.unitize();
				kR = omega.unitize();

			#ifdef DEBUG_ROTFRAME_TwoBodyRotatingFrameUnitVectors
				iR(0).setDerivative(adIndex++, 1.0);
				iR(1).setDerivative(adIndex++, 1.0);
				iR(2).setDerivative(adIndex++, 1.0);
			#endif

			#ifdef DEBUG_ROTFRAME_TwoBodyRotatingFrameUnitVectors
				kR(0).setDerivative(adIndex++, 1.0);
				kR(1).setDerivative(adIndex++, 1.0);
				kR(2).setDerivative(adIndex++, 1.0);
			#endif
				jR = kR.cross(iR);

				// put unit vectors in matrix
				for (size_t i = 0; i < 3; ++i)
				{
					unitVectors(i, 0) = iR(i);
					unitVectors(i, 1) = jR(i);
					unitVectors(i, 2) = kR(i);
				}
				if (generateDerivatives)
				{
					drdt = dStateB2WrtB1InertialDt.getSubMatrix1D(0, 2);
					dvdt = dStateB2WrtB1InertialDt.getSubMatrix1D(3, 5);

					math::Matrix<T> diRdR(3, 3), dkRdOmega(3, 3), djRdiR(3, 3), djRdkR(3, 3), djRdx(6, 3);
					math::Matrix<T> diRdt(3, 1), djRdt(3, 1), dkRdt(3, 1);
					diRdR = iR.unitDerivative(r); // operates on the unit vector and takes as input the un-unitized vector
					diRdt = diRdR * drdt;

					dkRdOmega = kR.unitDerivative(omega); // operates on the unit vector and takes as input the un-unitized vector
					dkRdt = dkRdOmega * dOmegaDt;

					djRdx = kR.crossDerivative(iR); // contains djR/dkR and djR/diR

					for (size_t i = 0; i < 3; ++i)
					{
						for (size_t j = 0; j < 3; ++j)
						{
							djRdkR(j, i) = djRdx(i, j); // note index flipping! emtg math organizes things in a way that may be counterintuitive
							djRdiR(j, i) = djRdx(i + 3, j); // note index flipping! emtg math organizes things in a way that may be counterintuitive
						}
					}

					djRdt = djRdkR * dkRdt + djRdiR * diRdt;

					for (size_t i = 0; i < 3; ++i)
					{
						dUnitVectorsDt(i, 0) = diRdt(i);
						dUnitVectorsDt(i, 1) = djRdt(i);
						dUnitVectorsDt(i, 2) = dkRdt(i);
					}

			#ifdef DEBUG_ROTFRAME_TwoBodyRotatingFrameUnitVectors
					std::cout << std::setprecision(16);
					std::cout << "\n\nTwoBodyRotatingFrameUnitVectors derivatives checking:\n\n";
					std::cout << "diR/dr analytical - AD:\n";
					for (size_t i = 0; i < 3; ++i)
					{
						for (size_t j = 0; j < 3; ++j)
						{
							std::cout << "(" << i << ", " << j << "): " << diRdR(i, j) - iR(i).getDerivative(j + 1) << "\n";
						}
					}
					std::cout << "\ndiR/dv analytical - AD:\n";
					for (size_t i = 0; i < 3; ++i)
					{
						for (size_t j = 0; j < 3; ++j)
						{
							std::cout << "(" << i << ", " << j << "): " << 0.0 - iR(i).getDerivative(j + 4) << "\n"; // manually put in 0 because, analytically, it should be 0
						}
					}
					std::cout << "\ndiR/dt analytical - AD:\n";
					for (size_t i = 0; i < 3; ++i)
					{
						std::cout << "(" << i << "): " << diRdt(i) - iR(i).getDerivative(0) << "\n";
					}
					std::cout << "\ndkR/dOmega analytical - AD:\n";
					for (size_t i = 0; i < 3; ++i)
					{
						for (size_t j = 0; j < 3; ++j)
						{
							std::cout << "(" << i << ", " << j << "): " << dkRdOmega(i, j) - kR(i).getDerivative(j + 7) << "\n";
						}
					}
					std::cout << "\ndkR/dt analytical - AD:\n";
					for (size_t i = 0; i < 3; ++i)
					{
						std::cout << "(" << i << "): " << dkRdt(i) - kR(i).getDerivative(0) << "\n";
					}
					std::cout << "\ndjR/diR analytical - AD:\n";
					for (size_t i = 0; i < 3; ++i)
					{
						for (size_t j = 0; j < 3; ++j)
						{
							std::cout << "(" << i << ", " << j << "): " << djRdiR(i, j) - jR(i).getDerivative(j + 10) << "\n";
						}
					}
					std::cout << "\ndjR/dkR analytical - AD:\n";
					for (size_t i = 0; i < 3; ++i)
					{
						for (size_t j = 0; j < 3; ++j)
						{
							std::cout << "(" << i << ", " << j << "): " << djRdkR(i, j) - jR(i).getDerivative(j + 13) << "\n";
						}
					}
					std::cout << "\ndjR/dt analytical - AD:\n";
					for (size_t i = 0; i < 3; ++i)
					{
						std::cout << "(" << i << "): " << djRdt(i) - jR(i).getDerivative(0) << "\n";
					}
			#endif
				}
				return;
			}


			//------------------------------------------------------------------------------
			// CalculateStateInRotatingFrameWrtB2(const math::Matrix<T>& stateB1WrtCBInertial, const math::Matrix<T>& dStateB1WrtCBInertialDt,
			// const math::Matrix<T>& stateB2WrtCBInertial, const math::Matrix<T>& dStateB2WrtCBInertialDt,
			//	const math::Matrix<T>& stateSCWrtCBInertial,
			//	math::Matrix<T>& stateSCWrtB2Rot,
			//	const bool& generateDerivatives, math::Matrix<T>& dStateSCWrtB2RotDStateInertial, math::Matrix<T>& dStateSCWrtB2RotDt)
			//------------------------------------------------------------------------------
			/**
			 * Calculate the position and velocity vector of the spacecraft with respect to
			 * body 2 in the frame rotating with the angular velocity of body 2 about body 1.
			 * Optionally return derivatives of position and velocity vectors with respect
			 * to the inertial position and velocity vectors of the spacecraft w/r/t/
			 * the central body and time.
			 *
			 * @param stateB1WrtCBInertial 6x1 EMTG matrix. Contains position and velocity vectors of B1 with respect to the central body in an inertial frame.
			 * @param dStateB1WrtCBInertialDt 6x1 EMTG matrix. Contains time derivative of stateB1WrtCBInertial. Only used if generateDerivatives it true.
			 * @param stateB2WrtCBInertial 6x1 EMTG matrix. Contains position and velocity vectors of B2 with respect to the central body in an inertial frame.
			 * @param dStateB2WrtCBInertialDt 6x1 EMTG matrix. Contains time derivative of stateB2WrtCBInertial. Only used if generateDerivatives it true.
			 * @param stateSCWrtCBInertial 6x1 EMTG matrix. Contains position and velocity vectors of the spacecraft with respect to the central body in an inertial frame.
			 * @return stateSCWrtB2Rot 6x1 EMTG matrix. Contains position and velocity vectors of the spacecraft with respect to B2 in the rotating frame.
			 * @param generateDerivatives Boolean. If true, generate derivatives of the returned values. If false, don't.
			 * @return dStateSCWrtB2RotDStateInertial 6x6 EMTG matrix. Derivatives of stateSCWrtB2Rot w/r/t/ stateSCWrtCBInertial.
			 * @return dUnitVectorsDt 6x1 EMTG matrix. Derivatives of stateSCWrtB2Rot w/r/t/ time.
			 */
			 //------------------------------------------------------------------------------
			template<class T> void CalculateStateInRotatingFrameWrtB2(const math::Matrix<T>& stateB1WrtCBInertial, const math::Matrix<T>& dStateB1WrtCBInertialDt,
				const math::Matrix<T>& stateB2WrtCBInertial, const math::Matrix<T>& dStateB2WrtCBInertialDt,
				const math::Matrix<T>& stateSCWrtCBInertial,
				math::Matrix<T>& stateSCWrtB2Rot,
				const bool& generateDerivatives, math::Matrix<T>& dStateSCWrtB2RotDStateInertial, math::Matrix<T>& dStateSCWrtB2RotDt)
			{
				// unpack state vectors into positions and velocities
				math::Matrix<T> rB1WrtCBInertial(3, 1), vB1WrtCBInertial(3, 1), drB1WrtCBInertialDt(3, 1), dvB1WrtCBInertialDt(3, 1);
				math::Matrix<T> rB2WrtCBInertial(3, 1), vB2WrtCBInertial(3, 1), drB2WrtCBInertialDt(3, 1), dvB2WrtCBInertialDt(3, 1);
				math::Matrix<T> rSCWrtCBInertial(3, 1), vSCWrtCBInertial(3, 1);

				rB1WrtCBInertial = stateB1WrtCBInertial.getSubMatrix1D(0, 2);
				vB1WrtCBInertial = stateB1WrtCBInertial.getSubMatrix1D(3, 5);
				drB1WrtCBInertialDt = dStateB1WrtCBInertialDt.getSubMatrix1D(0, 2);
				dvB1WrtCBInertialDt = dStateB1WrtCBInertialDt.getSubMatrix1D(3, 5);
				rB2WrtCBInertial = stateB2WrtCBInertial.getSubMatrix1D(0, 2);
				vB2WrtCBInertial = stateB2WrtCBInertial.getSubMatrix1D(3, 5);
				drB2WrtCBInertialDt = dStateB2WrtCBInertialDt.getSubMatrix1D(0, 2);
				dvB2WrtCBInertialDt = dStateB2WrtCBInertialDt.getSubMatrix1D(3, 5);
				rSCWrtCBInertial = stateSCWrtCBInertial.getSubMatrix1D(0, 2);
				vSCWrtCBInertial = stateSCWrtCBInertial.getSubMatrix1D(3, 5);

				// basic intermediate quantities

				// position, velocity, dr/dt, dv/dt of B2 w/r/t/ B1
				math::Matrix<T> stateB2WrtB1Inertial(6, 1), dStateB2WrtB1InertialDt(6, 1), rB2WrtB1Inertial(3, 1), vB2WrtB1Inertial(3, 1), drB2WrtB1InertialDt(3, 1), dvB2WrtB1InertialDt(3, 1);
				stateB2WrtB1Inertial = stateB2WrtCBInertial - stateB1WrtCBInertial;
				rB2WrtB1Inertial = stateB2WrtB1Inertial.getSubMatrix1D(0, 2);
				vB2WrtB1Inertial = stateB2WrtB1Inertial.getSubMatrix1D(3, 5);
				dStateB2WrtB1InertialDt = dStateB2WrtCBInertialDt - dStateB1WrtCBInertialDt;
				drB2WrtB1InertialDt = dStateB2WrtB1InertialDt.getSubMatrix1D(0, 2);
				dvB2WrtB1InertialDt = dStateB2WrtB1InertialDt.getSubMatrix1D(3, 5);

				// inertial position and velocity of s/c w/r/t/ B1
				math::Matrix<T> stateSCWrtB1Inertial(6, 1), rSCWrtB1Inertial(3, 1), vSCWrtB1Inertial(3, 1);
				stateSCWrtB1Inertial = stateSCWrtCBInertial - stateB1WrtCBInertial;
				rSCWrtB1Inertial = stateSCWrtB1Inertial.getSubMatrix1D(0, 2);
				vSCWrtB1Inertial = stateSCWrtB1Inertial.getSubMatrix1D(3, 5);

				// inertial position and velocity of s/c w/r/t/ B2
				math::Matrix<T> stateSCWrtB2Inertial(6, 1), rSCWrtB2Inertial(3, 1), vSCWrtB2Inertial(3, 1);
				stateSCWrtB2Inertial = stateSCWrtCBInertial - stateB2WrtCBInertial;
				rSCWrtB2Inertial = stateSCWrtB2Inertial.getSubMatrix1D(0, 2);
				vSCWrtB2Inertial = stateSCWrtB2Inertial.getSubMatrix1D(3, 5);

				// the angular velocity vector
				math::Matrix<T> omega(3, 1), dOmegaDt(3, 1), omegaCross(3, 3, 0.0);
				B2WrtB1AngularVelocityVector(stateB2WrtB1Inertial, dStateB2WrtB1InertialDt,
					omega, generateDerivatives, dOmegaDt);
				omegaCross(0, 1) = -omega(2);
				omegaCross(0, 2) = omega(1);
				omegaCross(1, 0) = omega(2);
				omegaCross(1, 2) = -omega(0);
				omegaCross(2, 0) = -omega(1);
				omegaCross(2, 1) = omega(0);

				// the unit vectors of the rotating frame
				math::Matrix<T> unitVectors(3, 3), dUnitVectorsDt(3, 3);
				TwoBodyRotatingFrameUnitVectors(stateB2WrtB1Inertial, dStateB2WrtB1InertialDt,
					unitVectors, generateDerivatives, dUnitVectorsDt);

				// transformation matrix from inertial to rotating
				math::Matrix<T> RinertialToRotating(3, 3);
				RinertialToRotating = unitVectors.transpose();

				// position of s/c w/r/t/ B2 in rotating frame
				math::Matrix<T> rSCWrtB2Rotating(3, 1);
				rSCWrtB2Rotating = RinertialToRotating * rSCWrtB2Inertial;

				// velocity of s/c w/r/t/ B2 in rotating frame, expressed in rotating frame
				math::Matrix<T> vSCWrtB2Rotating(3, 1);
				vSCWrtB2Rotating = RinertialToRotating * (vSCWrtB2Inertial - omega.cross(rSCWrtB2Inertial));

				
				// pack up rotating state vector
				for (size_t i = 0; i < 3; ++i)
				{
					stateSCWrtB2Rot(i) = rSCWrtB2Rotating(i);
					stateSCWrtB2Rot(i + 3) = vSCWrtB2Rotating(i);
				}

				// derivatives, if requested
				if (generateDerivatives)
				{
					for (size_t i = 0; i < 6; ++i)
					{
						dStateSCWrtB2RotDt(i) = 0.0;
						for (size_t j = 0; j < 6; ++j)
						{
							dStateSCWrtB2RotDStateInertial(i, j) = 0.0;
						}
					}

					// time derivative of transformation matrix
					math::Matrix<T> dRinertialToRotatingDt(3, 3);
					dRinertialToRotatingDt = dUnitVectorsDt.transpose();

					// time derivative of position of s/c w/r/t/ B2 in rotating frame
					math::Matrix<T> drSCWrtB2RotDt(3, 1), drSCWrtB2DtInertial(3, 1);
					drSCWrtB2DtInertial = -drB2WrtCBInertialDt;
					drSCWrtB2RotDt = dRinertialToRotatingDt * rSCWrtCBInertial - (dRinertialToRotatingDt * rB2WrtCBInertial + RinertialToRotating * drB2WrtCBInertialDt);


					// state derivatives of position of s/c w/r/t/ B2 in rotating frame
					// position w/r/t/ position is the transformation matrix; others are 0
					for (size_t i = 0; i < 3; ++i)
					{
						for (size_t j = 0; j < 3; ++j)
						{
							dStateSCWrtB2RotDStateInertial(i, j) = RinertialToRotating(i, j);
						}
					}

					// derivatives of velocity of s/c w/r/t/ B2 in rotating frame, expressed in rotating frame
					// there are several intermediate terms we need to build up, then apply the chain rule for time and state

					// derivative of inertial velocity of s/c w/r/t/ central body
					math::Matrix<T> dvSCWrtCBInertialDstate(3, 6, 0.0); // 0 w/r/t position, identity w/r/t/ velocity
					dvSCWrtCBInertialDstate(0, 3) = 1.0;
					dvSCWrtCBInertialDstate(1, 4) = 1.0;
					dvSCWrtCBInertialDstate(2, 5) = 1.0;
					math::Matrix<T> dvSCWrtCBInertialDt(3, 1, 0.0); // comes from state vector, so independent of time

					// derivative of inertial velocity of B2 w/r/t/ central body
					//math::Matrix<T> dvB1WrtCBInertialDt(3, 1, 0.0);
					//dvB1WrtCBInertialDt = dvB1WrtCBInertialDt;
					math::Matrix<T> dvB2WrtCBInertialDstate(3, 6, 0.0); // ephemeris lookup, so independent of state

					// derivative of inertial velocity of s/c w/r/t/ B1
					math::Matrix<T> dvSCWrtB2InertialDt(3, 1, 0.0), dvSCWrtB2InertialDstate(3, 6);
					dvSCWrtB2InertialDt = dvSCWrtCBInertialDt - dvB2WrtCBInertialDt;
					dvSCWrtB2InertialDstate = dvSCWrtCBInertialDstate - dvB2WrtCBInertialDstate;

					// derivative of omega \times position vector of s/c w/r/t/ B2

					// derivative of omega \times position vector of s/c w/r/t/ B2 w/r/t/ its arguments omega and position vector of s/c w/r/t/ B2
					math::Matrix<T> dOmegaCrossRSCWrtB2Dx(6, 3), dOmegaCrossRSCWrtB2Domega(3, 3), dOmegaCrossRSCWrtB2DrSCWrtB2Inertial(3, 3);
					dOmegaCrossRSCWrtB2Dx = omega.crossDerivative(rSCWrtB2Inertial);

					for (size_t i = 0; i < 3; ++i)
					{
						for (size_t j = 0; j < 3; ++j)
						{
							dOmegaCrossRSCWrtB2Domega(j, i) = dOmegaCrossRSCWrtB2Dx(i, j); // note index flipping! emtg math organizes things in a way that may be counterintuitive
							dOmegaCrossRSCWrtB2DrSCWrtB2Inertial(j, i) = dOmegaCrossRSCWrtB2Dx(i + 3, j); // note index flipping! emtg math organizes things in a way that may be counterintuitive
						}
					}

					// already have derivative of omega w/r/t/ time (dOmegaDt)

					// derivative of omega w/r/t/ state is 0
					math::Matrix<T> dOmegaDstate(3, 6, 0.0);

					// derivative of position vector of s/c w/r/t/ B2
					math::Matrix<T> drSCWrtB2InertialDt(3, 1), drSCWrtB2InertialDstate(3, 6, 0.0);
					drSCWrtB2InertialDt = -drB2WrtCBInertialDt; // position of s/c w/r/t/ CB is independent of time

					// position of B2 w/r/t/ CB is independent of state. position of s/c w/r/t/ CB is identity for position submatrix and zero for velocity submatrix
					drSCWrtB2InertialDstate(0, 0) = 1.0;
					drSCWrtB2InertialDstate(1, 1) = 1.0;
					drSCWrtB2InertialDstate(2, 2) = 1.0;

					// chain rule
					// derivative of cross product w/r/t/ time
					math::Matrix<T> dOmegaCrossRSCWrtB2Dt(3, 1);
					dOmegaCrossRSCWrtB2Dt = dOmegaCrossRSCWrtB2Domega * dOmegaDt + dOmegaCrossRSCWrtB2DrSCWrtB2Inertial * drSCWrtB2InertialDt;

					// derivative of cross product w/r/t/ state
					math::Matrix<T> dOmegaCrossRSCWrtB2Dstate(3, 6);
					dOmegaCrossRSCWrtB2Dstate = dOmegaCrossRSCWrtB2Domega * dOmegaDstate + dOmegaCrossRSCWrtB2DrSCWrtB2Inertial * drSCWrtB2InertialDstate;

					// full derivative of velocity in rotating frame w/r/t time
					//*******************************
					math::Matrix<T> dvSCWrtB2RotDt(3, 1);
					//dvSCWrtB2RotDt = dvSCWrtB1InertialDt - dOmegaCrossRSCWrtB1Dt;
					dvSCWrtB2RotDt = dRinertialToRotatingDt * (vSCWrtB2Inertial - omega.cross(rSCWrtB2Inertial)) + RinertialToRotating * (dvSCWrtB2InertialDt - dOmegaCrossRSCWrtB2Dt);

					// full derivative of velocity in rotating frame w/r/t state
					//**********************************
					math::Matrix<T> dvSCWrtB2RotDStateInertial(3, 6);
					//dvSCWrtB2RotDStateInertial = dvSCWrtB1InertialDstate - dOmegaCrossRSCWrtB1Dstate;
					dvSCWrtB2RotDStateInertial = RinertialToRotating * (dvSCWrtB2InertialDstate - dOmegaCrossRSCWrtB2Dstate); // no state dependence of transformation matrix

					// pack up time derivatives
					for (size_t i = 0; i < 3; ++i)
					{
						dStateSCWrtB2RotDt(i) = drSCWrtB2RotDt(i);
						dStateSCWrtB2RotDt(i + 3) = dvSCWrtB2RotDt(i);
					}

					// pack up state derivatives (only need for derivatives of velocity; derivatives of position were already packed up)
					for (size_t i = 0; i < 3; ++i)
					{
						for (size_t j = 0; j < 6; ++j)
						{
							dStateSCWrtB2RotDStateInertial(i + 3, j) = dvSCWrtB2RotDStateInertial(i, j);
						}
					}

			#ifdef DEBUG_ROTFRAME_CalculateStateInRotatingFrameWrtB2
					std::cout << std::setprecision(16);
					std::cout << "drSCWrtB2InertialDt analytical - AD:\n";
					for (size_t i = 0; i < 3; ++i)
					{
						std::cout << "(" << i << "): " << drSCWrtB2DtInertial(i) - rSCWrtB2Inertial(i).getDerivative(0) << "\n";
					}
					std::cout << "\ndRinertialToRotating/dt analytical - AD:\n";
					for (size_t i = 0; i < 3; ++i)
					{
						for (size_t j = 0; j < 3; ++j)
						{
							std::cout << "(" << i << ", " << j << "): " << dRinertialToRotatingDt(i, j) - RinertialToRotating(i, j).getDerivative(0) << "\n";
						}
					}

					std::cout << "\n\nstateSCWrtB2Rot derivatives checking:\n\n";
					std::cout << "dStateSCWrtB2RotDt analytical - AD:\n";
					for (size_t i = 0; i < 6; ++i)
					{
						std::cout << "(" << i << "): " << dStateSCWrtB2RotDt(i) - stateSCWrtB2Rot(i).getDerivative(0) << "\n";
					}
					std::cout << "\ndStateSCWrtB2RotDStateInertial analytical - AD:\n";
					for (size_t i = 0; i < 6; ++i)
					{
						for (size_t j = 0; j < 6; ++j)
						{
							std::cout << "(" << i << ", " << j << "): " << dStateSCWrtB2RotDStateInertial(i, j) - stateSCWrtB2Rot(i).getDerivative(j + 1) << "\n";
						}
					}
			#endif
				}
				return;
			}
		} // close namespace TwoBodyRotatingFrame
    }//close namespace Astrodynamics
}//close namespace EMTG