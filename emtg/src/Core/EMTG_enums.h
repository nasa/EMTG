// EMTG: Evolutionary Mission Trajectory Generator
// An open-source global optimization tool for preliminary mission design
// Provided by NASA Goddard Space Flight Center
//
// Copyright (c) 2013 - 2020 United States Government as represented by the
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

//header file containing EMTG enums
//Jacob Englander 7-6-2017
//this file exists so that we do not have to remember what all of the switches mean!

#pragma once

namespace EMTG
{
    //***********************************************
    //solver options enums

    //***********************************************
    //global options enums
    enum PhaseType {MGALTS, FBLTS, MGALT, FBLT, PSBI, PSFB, MGAnDSMs, CoastPhase, SundmanCoastPhase, VARIABLE_PHASE_TYPE};

    enum StateRepresentation{Cartesian, SphericalRADEC, SphericalAZFPA, COE};

    enum ControlMagnitudeType {UpToUnitMagnitude, UnitMagnitude, ZeroMagnitude};

    enum ObjectiveFunctionType {MINIMIZE_DELTAV,
                                MINIMIZE_TIME,
                                MAXIMIZE_MASS,
                                MAXIMIZE_INITIAL_MASS,
                                DEPART_AS_LATE_AS_POSSIBLE,
                                DEPART_AS_EARLY_AS_POSSIBLE,
                                MAXIMIZE_ORBIT_ENERGY,
                                MINIMIZE_LAUNCH_MASS,
                                ARRIVE_AS_EARLY_AS_POSSIBLE,
                                ARRIVE_AS_LATE_AS_POSSIBLE,
                                MINIMIZE_PROPELLANT,
                                MAXIMIZE_WET_DRY_RATIO,
                                MAXIMIZE_ARRIVAL_KINETIC_ENERGY,
                                MINIMIZE_BOL_POWER,
                                MAXIMIZE_LOG10_MASS,
                                MAXIMIZE_LOGE_MASS,
                                MAXIMIZE_DRY_MASS_MARGIN,
                                MAXIMIZE_DRY_MASS,
                                MAXIMIZE_LOG10_DRY_MASS,
                                MAXIMIZE_LOGE_DRY_MASS,
                                MINIMIZE_CHEMICAL_FUEL,
                                MINIMIZE_CHEMICAL_OXIDIZER,
                                MINIMIZE_ELECTRIC_PROPELLANT,
                                MINIMIZE_TOTAL_PROPELLANT,
                                MINIMIZE_WAYPOINT_TRACKING_ERROR,
                                MINIMIZE_INITIAL_IMPULSE,
                                MAXIMIZE_FINAL_PERIAPSIS};

    //***********************************************
    //hardware options enums
    enum SpacecraftModelInputType {AssembleFromLibraries, ReadSpacecraftFile, AssembleFromMissionOptions};
    enum SpacecraftBusPowerType {TypeA_Quadratic, TypeB_Conditional};
    enum SpacecraftPowerSupplyType {ConstantPower, Solar};
    enum SpacecraftPowerSupplyCurveType {Sauer, Polynomial};

    enum SpacecraftThrusterMode {ConstantThrustIsp, 
                                 FixedEfficiencyCSI, 
                                 FixedEfficiencyVSI,
                                 Poly1D, 
                                 SteppedHThrust1D, 
                                 SteppedLMdot1D,
                                 SteppedHisp1D,
                                 SteppedHEfficiency1D,
                                 SteppedFullSet1D,
                                 Stepped2D,
                                 Poly2D};

    enum ThrottleLogic{MaxThrusters, MinThrusters};
    enum PropulsionSystemChoice{ Monoprop, Biprop, EP};

    //***********************************************
    //journey options enums

    enum DepartureType{LAUNCH_OR_DIRECT_INSERTION, 
                       DEPART_PARKING_ORBIT,
                       FREE_DIRECT_DEPARTURE,
                       FLYBY,
                       FLYBY_FIXED_VINF,
                       SPIRAL_ESCAPE, 
                       ZERO_TURN_FLYBY};

    enum ArrivalType{INSERTION_INTO_PARKING, 
                     CHEM_RENDEZVOUS, 
                     INTERCEPT, 
                     LT_RENDEZVOUS,
                     CHEM_MATCH_VINF,
                     LT_MATCH_VINF, 
                     SPIRAL_CAPTURE,
					 MOMENTUM_EXCHANGE};

    enum BoundaryClass{EphemerisPegged, 
                       FreePoint,
                       EphemerisReferenced,
                       Periapse};

	enum MomentumExchangeType{Kinetic,
							  ExplosiveDevice};

    //***********************************************

    //physics enums
    enum PropagatorType {KeplerPropagator, IntegratedPropagator };
    enum IntegratorType {rk7813m_adaptive, rk8_fixed};
	enum PropagationDomain {Time, Sundman};
    enum DutyCycleType {Averaged, Realistic};
    enum AtmosphereModel { Exponential };

    //frame rotation enums
    enum ReferenceFrame {ICRF, J2000_BCI, J2000_BCF, TrueOfDate_BCI, TrueOfDate_BCF, PrincipleAxes, Topocentric, Polar, SAM};

    //***********************************************
    //solver enums
    enum InnerLoopSolverType {RUN_TRIALX, MBH, ACDE, NLP, FilamentWalker};
    enum NLPMode {FeasiblePoint, Optimize, FilamentFinder};
}