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

#include <vector>
#include <string>

namespace EMTG
{
    //***********************************************
    //solver options enums

    //***********************************************
    //global options enums
    enum PhaseType {MGALTS, FBLTS, MGALT, FBLT, PSBI, PSFB, MGAnDSMs, CoastPhase, SundmanCoastPhase, VARIABLE_PHASE_TYPE, ProbeEntryPhase, ControlLawThrustPhase};
    const std::vector<std::string> PhaseTypeStrings({ "MGALTS", "FBLTS", "MGALT", "FBLT", "PSBI", "PSFB", "MGAnDSMs", "CoastPhase", "SundmanCoastPhase", "VARIABLE_PHASE_TYPE", "ProbeEntryPhase", "ControlLawThrustPhase" });

    enum StateRepresentation {Cartesian, SphericalRADEC, SphericalAZFPA, COE, MEE, IncomingBplane, OutgoingBplane };
    const std::vector<std::string> StateRepresentationStrings({ "Cartesian", "SphericalRADEC", "SphericalAZFPA", "COE", "MEE", "IncomingBplane", "OutgoingBplane" });

    enum ConstraintStateRepresentation {CartesianConstraint, Native}; //"native" means "use the same state representation as the encoded states"
    const std::vector<std::string> ConstraintStateRepresentationStrings({ "CartesianConstraint", "NativeConstraint" });

    enum ControlMagnitudeType {UpToUnitMagnitude, UnitMagnitude, ZeroMagnitude};
    const std::vector<std::string> ControlMagnitudeTypeString({ "UpToUnitMagnitude", "UnitMagnitude", "ZeroMagnitude" });

    enum ThrustControlLaw {CartesianUnit, Velocity, AntiVelocity};
    const std::vector<std::string> ThrustControlLawStrings({ "CartesianUnit", "Velocity", "AntiVelocity" });

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
                                PLACEHOLDER_OBJECTIVE_FUNCTION,
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

    const std::vector<std::string> ObjectiveFunctionTypeStrings ({
                                "MINIMIZE_DELTAV",
                                "MINIMIZE_TIME",
                                "MAXIMIZE_MASS",
                                "MAXIMIZE_INITIAL_MASS",
                                "DEPART_AS_LATE_AS_POSSIBLE",
                                "DEPART_AS_EARLY_AS_POSSIBLE",
                                "MAXIMIZE_ORBIT_ENERGY",
                                "MINIMIZE_LAUNCH_MASS",
                                "ARRIVE_AS_EARLY_AS_POSSIBLE",
                                "ARRIVE_AS_LATE_AS_POSSIBLE",
                                "PLACEHOLDER_OBJECTIVE_FUNCTION",
                                "MAXIMIZE_WET_DRY_RATIO",
                                "MAXIMIZE_ARRIVAL_KINETIC_ENERGY",
                                "MINIMIZE_BOL_POWER",
                                "MAXIMIZE_LOG10_MASS",
                                "MAXIMIZE_LOGE_MASS",
                                "MAXIMIZE_DRY_MASS_MARGIN",
                                "MAXIMIZE_DRY_MASS",
                                "MAXIMIZE_LOG10_DRY_MASS",
                                "MAXIMIZE_LOGE_DRY_MASS",
                                "MINIMIZE_CHEMICAL_FUEL",
                                "MINIMIZE_CHEMICAL_OXIDIZER",
                                "MINIMIZE_ELECTRIC_PROPELLANT",
                                "MINIMIZE_TOTAL_PROPELLANT",
                                "MINIMIZE_WAYPOINT_TRACKING_ERROR",
                                "MINIMIZE_INITIAL_IMPULSE",
                                "MAXIMIZE_FINAL_PERIAPSIS"
                            });

    //***********************************************
    //hardware options enums
    enum SpacecraftModelInputType {AssembleFromLibraries, ReadSpacecraftFile, AssembleFromMissionOptions};
    const std::vector<std::string> SpacecraftModelInputTypeStrings ({ "AssembleFromLibraries", "ReadSpacecraftFile", "AssembleFromMissionOptions" });

    enum SpacecraftBusPowerType {TypeA_Quadratic, TypeB_Conditional};
    const std::vector<std::string> SpacecraftBusPowerTypeStrings ({ "TypeA_Quadratic", "TypeB_Conditional" });

    enum SpacecraftPowerSupplyType {ConstantPower, Solar};
    const std::vector<std::string> SpacecraftPowerSupplyTypeStrings ({ "ConstantPower", "Solar" });

    enum SpacecraftPowerSupplyCurveType {Sauer, Polynomial};
    const std::vector<std::string> SpacecraftPowerSupplyCurveTypeStrings ({ "Sauer", "Polynomial" });

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
    const std::vector<std::string> SpacecraftThrusterModeStrings ({
                                "ConstantThrustIsp",
                                "FixedEfficiencyCSI",
                                "FixedEfficiencyVSI",
                                "Poly1D",
                                "SteppedHThrust1D",
                                "SteppedLMdot1D",
                                "SteppedHisp1D",
                                "SteppedHEfficiency1D",
                                "SteppedFullSet1D",
                                "Stepped2D",
                                "Poly2D"
                            });

    enum ThrottleLogic{MaxThrusters, MinThrusters};
    const std::vector<std::string> ThrottleLogicStrings ({ "MaxThrusters", "MinThrusters" });

    enum PropulsionSystemChoice{ Monoprop, Biprop, EP};
    const std::vector<std::string> PropulsionSystemChoiceStrings ({ "Monoprop", "Biprop", "EP" });

    //***********************************************
    //journey options enums

    enum DepartureType{LAUNCH_OR_DIRECT_INSERTION, 
                       DEPART_PARKING_ORBIT,
                       FREE_DIRECT_DEPARTURE,
                       FLYBY,
                       FLYBY_FIXED_VINF,
                       SPIRAL_ESCAPE, 
                       ZERO_TURN_FLYBY};
    const std::vector<std::string> DepartureTypeStrings ({
                       "LAUNCH_OR_DIRECT_INSERTION",
                       "DEPART_PARKING_ORBIT",
                       "FREE_DIRECT_DEPARTURE",
                       "FLYBY",
                       "FLYBY_FIXED_VINF",
                       "SPIRAL_ESCAPE",
                       "ZERO_TURN_FLYBY"
                    });

    enum ArrivalType{INSERTION_INTO_PARKING, 
                     CHEM_RENDEZVOUS, 
                     INTERCEPT, 
                     LT_RENDEZVOUS,
                     CHEM_MATCH_VINF,
                     LT_MATCH_VINF, 
                     SPIRAL_CAPTURE,
					 MOMENTUM_EXCHANGE};
    const std::vector<std::string> ArrivalTypeStrings ({
                    "INSERTION_INTO_PARKING",
                    "CHEM_RENDEZVOUS",
                    "INTERCEPT",
                    "LT_RENDEZVOUS",
                    "CHEM_MATCH_VINF",
                    "LT_MATCH_VINF",
                    "SPIRAL_CAPTURE",
                    "MOMENTUM_EXCHANGE"
                });

    enum BoundaryClass{EphemerisPegged, 
                       FreePoint,
                       EphemerisReferenced,
                       Periapse};
    const std::vector<std::string> BoundaryClassStrings ({
                       "EphemerisPegged",
                       "FreePoint",
                       "EphemerisReferenced",
                       "Periapse"
                 });

	enum MomentumExchangeType{Kinetic,
							  ExplosiveDevice};
    const std::vector<std::string> MomentumExchangeTypeStrings ({
                             "Kinetic",
                             "ExplosiveDevice"
                        });

    //***********************************************

    //physics enums
    enum PropagatorType {KeplerPropagator, IntegratedPropagator };
    const std::vector<std::string> PropagatorTypeStrings ({ "KeplerPropagator", "IntegratedPropagator" });

    enum IntegratorType {rk7813m_adaptive, rk8_fixed};
    const std::vector<std::string> IntegratorTypeStrings ({ "rk7813m_adaptive", "rk8_fixed" });

    enum IntegrationCoefficientsType {rk4, rkdp87};
    const std::vector<std::string> IntegrationCoefficientsTypeStrings({ "rk4", "rkdp87" });

	enum PropagationDomain {Time, Sundman};
    const std::vector<std::string> PropagationDomainStrings ({ "Time", "Sundman" });

    enum DutyCycleType {Averaged, Realistic};
    const std::vector<std::string> DutyCycleTypeStrings ({ "Averaged", "Realistic" });
    
    enum AtmosphereModel { Exponential };
    const std::vector<std::string> AtmosphereModelStrings({ "Exponential" });

    //frame rotation enums
    enum ReferenceFrame {ICRF, J2000_BCI, J2000_BCF, TrueOfDate_BCI, TrueOfDate_BCF, PrincipleAxes, Topocentric, Polar, SAM, ObjectReferenced};
    const std::vector<std::string> ReferenceFrameStrings ({ "ICRF", "J2000_BCI", "J2000_BCF", "TrueOfDate_BCI", "TrueOfDate_BCF", "PrincipleAxes", "Topocentric", "Polar", "SAM", "ObjectReferenced" });

    //***********************************************
    //solver enums
    enum InnerLoopSolverType {RUN_TRIALX, MBH, ACDE, NLP, FilamentWalker};
    const std::vector<std::string> InnerLoopSolverTypeStrings ({ "RUN_TRIALX", "MBH", "ACDE", "NLP", "FilamentWalker" });

    enum NLPMode {FeasiblePoint, Optimize, FilamentFinder};
    const std::vector<std::string> NLPModeStrings ({ "FeasiblePoint", "Optimize", "FilamentFinder" });
}