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

//boundary constraint factory for ephemeris pegged boundary events
//11-1-2017

#include "SpecializedBoundaryConstraintFactory.h"

#include "BoundaryDistanceConstraint.h"
#include "BoundaryIsApseConstraint.h"
#include "BoundaryMassConstraint.h"
#include "BoundaryLongitudeConstraint.h"
#include "BoundaryBCFlatitudeConstraint.h"
#include "BoundaryDeticLatitudeConstraint.h"
#include "BoundaryDeticAltitudeConstraint.h"
#include "BoundaryOrbitalEnergyConstraint.h"
#include "BoundaryVelocityDeclinationConstraint.h"
#include "BoundaryVelocityMagnitudeConstraint.h"
#include "BoundaryVelocitySphericalAzimuthConstraint.h"
#include "BoundaryRelativeVelocityMagnitudeConstraint.h"
#include "BoundaryRelativeVelocityAzimuthConstraint.h"
#include "BoundaryRelativeVelocityHFPAConstraint.h"
#include "BoundaryTargetBodyDeticElevationConstraint.h"
#include "BoundaryVelocityVFPAconstraint.h"
#include "SemimajorAxisConstraint.h"
#include "InclinationConstraint.h"
#include "EccentricityConstraint.h"
#include "RAANConstraint.h"
#include "AOPConstraint.h"
#include "TrueAnomalyConstraint.h"
#include "OrbitPeriodConstraint.h"
#include "PeriapseDistanceConstraint.h"
#include "ApoapseDistanceConstraint.h"
#include "RBP.h"
#include "RPB.h"
#include "RRP.h"
#include "RPR.h"
#include "angularMomentumReferenceAngle.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        namespace SpecializedConstraints
        {
            SpecializedBoundaryConstraintBase* create_boundary_event_constraint(const std::string& name,
                const size_t& journeyIndex,
                const size_t& phaseIndex,
                const size_t& stageIndex,
                BoundaryEventBase* myBoundaryEvent,
                Astrodynamics::universe* Universe,
                HardwareModels::Spacecraft* mySpacecraft,
                missionoptions* myOptions,
                const std::string& ConstraintDefinition,
                const std::string& EventDefinition)
            {
                std::vector<std::string> ConstraintDefinitionCell;

                std::string splitty_thing = boost::algorithm::to_lower_copy(ConstraintDefinition);
                boost::split(ConstraintDefinitionCell, splitty_thing, boost::is_any_of("_"), boost::token_compress_on);

                if (EventDefinition == "departure")
                {
                    if (ConstraintDefinitionCell[2].find("distanceconstraint") < 1024)
                    {
                        return new BoundaryDistanceConstraint(name,
                            journeyIndex,
                            phaseIndex,
                            stageIndex,
                            Universe,
                            mySpacecraft,
                            myOptions,
                            myBoundaryEvent,
                            ConstraintDefinition);
                    }
                    else if (ConstraintDefinitionCell[2].find("isapse") < 1024)
                    {
                        return new BoundaryIsApseConstraint(name,
                            journeyIndex,
                            phaseIndex,
                            stageIndex,
                            Universe,
                            mySpacecraft,
                            myOptions,
                            myBoundaryEvent,
                            ConstraintDefinition);
                    }
                    else if (ConstraintDefinitionCell[2].find("massconstraint") < 1024)
                    {
                        return new BoundaryMassConstraint(name,
                            journeyIndex,
                            phaseIndex,
                            stageIndex,
                            Universe,
                            mySpacecraft,
                            myOptions,
                            myBoundaryEvent,
                            ConstraintDefinition);
                    }
                    else if (ConstraintDefinitionCell[2].find("sma") < 1024)
                    {
                        return new SemimajorAxisConstraint(name,
                            journeyIndex,
                            phaseIndex,
                            stageIndex,
                            Universe,
                            mySpacecraft,
                            myOptions,
                            myBoundaryEvent,
                            ConstraintDefinition);
                    }
                    else if (ConstraintDefinitionCell[2].find("inc") < 1024)
                    {
                        // if a frame is specified, then inform the constraint what the user wants
                        // otherwise terminate program
                        if (ConstraintDefinitionCell.size() > 5)
                        {
                            ReferenceFrame user_specified_frame;
                            if (ConstraintDefinitionCell[5].find("icrf") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::ICRF;
                            }
                            else if (ConstraintDefinitionCell[5].find("j2000bci") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::J2000_BCI;
                            }
                            else
                            {
                                throw std::invalid_argument("Currently only ICRF and J2000BCI frames are supported: inclination boundary constraint " + name);
                            }

                            return new InclinationConstraint(name,
                                journeyIndex,
                                phaseIndex,
                                stageIndex,
                                Universe,
                                mySpacecraft,
                                myOptions,
                                user_specified_frame,
                                myBoundaryEvent,
                                ConstraintDefinition);
                        }
                        else
                        {
                            throw std::invalid_argument("Reference frame must be specified for inclination boundary constraint " + name);
                        }
                    }
                    else if (ConstraintDefinitionCell[2].find("ecc") < 1024)
                    {
                        return new EccentricityConstraint(name,
                            journeyIndex,
                            phaseIndex,
                            stageIndex,
                            Universe,
                            mySpacecraft,
                            myOptions,
                            myBoundaryEvent,
                            ConstraintDefinition);
                    }
                    else if (ConstraintDefinitionCell[2].find("raan") < 1024)
                    {
                        // if a frame is specified, then inform the constraint what the user wants
                        // otherwise terminate program
                        if (ConstraintDefinitionCell.size() > 5)
                        {
                            ReferenceFrame user_specified_frame;
                            if (ConstraintDefinitionCell[5].find("icrf") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::ICRF;
                            }
                            else if (ConstraintDefinitionCell[5].find("j2000bci") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::J2000_BCI;
                            }
                            else
                            {
                                throw std::invalid_argument("Currently only ICRF and J2000BCI frames are supported: RAAN boundary constraint " + name);
                            }

                            return new RAANConstraint(name,
                                journeyIndex,
                                phaseIndex,
                                stageIndex,
                                Universe,
                                mySpacecraft,
                                myOptions,
                                user_specified_frame,
                                myBoundaryEvent,
                                ConstraintDefinition);
                        }
                        else
                        {
                            throw std::invalid_argument("Reference frame must be specified for RAAN boundary constraint " + name);
                        }
                    }
                    else if (ConstraintDefinitionCell[2].find("aop") < 1024)
                    {
                        // if a frame is specified, then inform the constraint what the user wants
                        // otherwise terminate program
                        if (ConstraintDefinitionCell.size() > 5)
                        {
                            ReferenceFrame user_specified_frame;
                            if (ConstraintDefinitionCell[5].find("icrf") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::ICRF;
                            }
                            else if (ConstraintDefinitionCell[5].find("j2000bci") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::J2000_BCI;
                            }
                            else
                            {
                                throw std::invalid_argument("Currently only ICRF and J2000BCI frames are supported: AOP boundary constraint " + name);
                            }

                            return new AOPConstraint(name,
                                journeyIndex,
                                phaseIndex,
                                stageIndex,
                                Universe,
                                mySpacecraft,
                                myOptions,
                                user_specified_frame,
                                myBoundaryEvent,
                                ConstraintDefinition);
                        }
                        else
                        {
                            throw std::invalid_argument("Reference frame must be specified for AOP boundary constraint " + name);
                        }
                    }
                    else if (ConstraintDefinitionCell[2].find("trueanomaly") < 1024)
                    {
                        return new TrueAnomalyConstraint(name,
                            journeyIndex,
                            phaseIndex,
                            stageIndex,
                            Universe,
                            mySpacecraft,
                            myOptions,
                            myBoundaryEvent,
                            ConstraintDefinition);
                    }
                    else if (ConstraintDefinitionCell[2].find("orbitperiod") < 1024)
                    {
                        return new OrbitPeriodConstraint(name,
                            journeyIndex,
                            phaseIndex,
                            stageIndex,
                            Universe,
                            mySpacecraft,
                            myOptions,
                            myBoundaryEvent,
                            ConstraintDefinition);
                    }
                    else if (ConstraintDefinitionCell[2].find("orbitalenergy") < 1024)
                    {
                        return new BoundaryOrbitalEnergyConstraint(name,
                            journeyIndex,
                            phaseIndex,
                            stageIndex,
                            Universe,
                            mySpacecraft,
                            myOptions,
                            myBoundaryEvent,
                            ConstraintDefinition);
                    }
                    else if (ConstraintDefinitionCell[2].find("periapsedistance") < 1024)
                    {
                        return new PeriapseDistanceConstraint(name,
                            journeyIndex,
                            phaseIndex,
                            stageIndex,
                            Universe,
                            mySpacecraft,
                            myOptions,
                            myBoundaryEvent,
                            ConstraintDefinition);
                    }
                    else if (ConstraintDefinitionCell[2].find("apoapsedistance") < 1024)
                    {
                        return new ApoapseDistanceConstraint(name,
                            journeyIndex,
                            phaseIndex,
                            stageIndex,
                            Universe,
                            mySpacecraft,
                            myOptions,
                            myBoundaryEvent,
                            ConstraintDefinition);
                    }
                    else if (ConstraintDefinitionCell[2].find("bcflatitude") < 1024)
                    {
                        return new BoundaryBCFlatitudeConstraint(name,
                            journeyIndex,
                            phaseIndex,
                            stageIndex,
                            Universe,
                            mySpacecraft,
                            myOptions,
                            myBoundaryEvent,
                            ConstraintDefinition);
                    }
                    else if (ConstraintDefinitionCell[2].find("longitude") < 1024)
                    {
                        // if a frame is specified, then inform the constraint what the user wants
                        // otherwise terminate program
                        if (ConstraintDefinitionCell.size() > 5)
                        {
                            ReferenceFrame user_specified_frame;
                            if (ConstraintDefinitionCell[5].find("j2000bcf") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::J2000_BCF;
                            }
                            else if (ConstraintDefinitionCell[5].find("trueofdatebcf") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::TrueOfDate_BCF;
                            }
                            else if (ConstraintDefinitionCell[5].find("icrf") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::ICRF;
                            }
                            else
                            {
                                throw std::invalid_argument("Currently only ICFR, J2000BCF and TrueOfDateBCF frames are supported: longitude boundary constraint " + name);
                            }

                            return new BoundaryLongitudeConstraint(name,
                                journeyIndex,
                                phaseIndex,
                                stageIndex,
                                Universe,
                                mySpacecraft,
                                myOptions,
                                user_specified_frame,
                                myBoundaryEvent,
                                ConstraintDefinition);
                        }
                        else
                        {
                            throw std::invalid_argument("Reference frame must be specified for longitude boundary constraint " + name);
                        }
                    }
                    else if (ConstraintDefinitionCell[2].find("deticlatitude") < 1024)
                    {
                        // if a frame is specified, then inform the constraint what the user wants
                        // otherwise terminate program
                        if (ConstraintDefinitionCell.size() > 5)
                        {
                            ReferenceFrame user_specified_frame;
                            if (ConstraintDefinitionCell[5].find("j2000bcf") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::J2000_BCF;
                            }
                            else if (ConstraintDefinitionCell[5].find("trueofdatebcf") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::TrueOfDate_BCF;
                            }
                            else
                            {
                                throw std::invalid_argument("Currently only J2000BCF and TrueOfDateBCF frames are supported: detic latitude boundary constraint " + name);
                            }

                            return new BoundaryDeticLatitudeConstraint(name,
                                journeyIndex,
                                phaseIndex,
                                stageIndex,
                                Universe,
                                mySpacecraft,
                                myOptions,
                                user_specified_frame,
                                myBoundaryEvent,
                                ConstraintDefinition);
                        }
                        else
                        {
                            throw std::invalid_argument("Reference frame must be specified for detic latitude boundary constraint " + name);
                        }
                    }
                    else if (ConstraintDefinitionCell[2].find("deticaltitude") < 1024)
                    {
                        // if a frame is specified, then inform the constraint what the user wants
                        // otherwise terminate program
                        if (ConstraintDefinitionCell.size() > 5)
                        {
                            ReferenceFrame user_specified_frame;
                            if (ConstraintDefinitionCell[5].find("j2000bcf") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::J2000_BCF;
                            }
                            else if (ConstraintDefinitionCell[5].find("trueofdatebcf") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::TrueOfDate_BCF;
                            }
                            else
                            {
                                throw std::invalid_argument("Currently only J2000BCF and TrueOfDateBCF frames are supported: detic altitude boundary constraint " + name);
                            }

                            return new BoundaryDeticAltitudeConstraint(name,
                                journeyIndex,
                                phaseIndex,
                                stageIndex,
                                Universe,
                                mySpacecraft,
                                myOptions,
                                user_specified_frame,
                                myBoundaryEvent,
                                ConstraintDefinition);
                        }
                        else
                        {
                            throw std::invalid_argument("Reference frame must be specified for detic altitude boundary constraint " + name);
                        }
                    }
                    else if (ConstraintDefinitionCell[2].find("deticelevation") < 1024)
                    {
                        // if a frame is specified, then inform the constraint what the user wants
                        // otherwise terminate program
                        if (ConstraintDefinitionCell.size() > 6)
                        {
                            ReferenceFrame user_specified_frame;
                            if (ConstraintDefinitionCell[6].find("trueofdatebcf") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::TrueOfDate_BCF;
                            }
                            else
                            {
                                throw std::invalid_argument("Currently only TrueOfDateBCF frames are supported: detic elevation boundary constraint " + name);
                            }

                            return new BoundaryTargetBodyDeticElevationConstraint(name,
                                journeyIndex,
                                phaseIndex,
                                stageIndex,
                                Universe,
                                mySpacecraft,
                                myOptions,
                                user_specified_frame,
                                myBoundaryEvent,
                                ConstraintDefinition);
                        }
                        else
                        {
                            throw std::invalid_argument("Reference frame must be specified for detic elevation boundary constraint " + name);
                        }
                    }
                    else if (ConstraintDefinitionCell[2].find("velocitydeclination") < 1024)
                    {
                        // if a frame is specified, then inform the constraint what the user wants
                        // otherwise terminate program
                        if (ConstraintDefinitionCell.size() > 5)
                        {
                            ReferenceFrame user_specified_frame;
                            if (ConstraintDefinitionCell[5].find("icrf") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::ICRF;
                            }
                            else if (ConstraintDefinitionCell[5].find("j2000bci") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::J2000_BCI;
                            }
                            else if (ConstraintDefinitionCell[5].find("j2000bcf") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::J2000_BCF;
                            }
                            else if (ConstraintDefinitionCell[5].find("trueofdatebci") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::TrueOfDate_BCI;
                            }
                            else if (ConstraintDefinitionCell[5].find("trueofdatebcf") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::TrueOfDate_BCF;
                            }
                            else
                            {
                                throw std::invalid_argument("Currently only ICRF, J2000BCI, J2000BCF, TrueOfDateBCI, and TrueOfDateBCF frames are supported: velocity declination boundary constraint " + name);
                            }

                            return new BoundaryVelocityDeclinationConstraint(name,
                                journeyIndex,
                                phaseIndex,
                                stageIndex,
                                Universe,
                                mySpacecraft,
                                myOptions,
                                user_specified_frame,
                                myBoundaryEvent,
                                ConstraintDefinition);
                        }
                    }
                    else if (ConstraintDefinitionCell[2].find("velocitymagnitude") < 1024)
                    {
                        // if a frame is specified, then inform the constraint what the user wants
                        // otherwise terminate program
                        if (ConstraintDefinitionCell.size() > 5)
                        {
                            ReferenceFrame user_specified_frame;
                            if (ConstraintDefinitionCell[5].find("icrf") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::ICRF;
                            }
                            else if (ConstraintDefinitionCell[5].find("j2000bci") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::J2000_BCI;
                            }
                            else if (ConstraintDefinitionCell[5].find("j2000bcf") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::J2000_BCF;
                            }
                            else if (ConstraintDefinitionCell[5].find("trueofdatebci") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::TrueOfDate_BCI;
                            }
                            else if (ConstraintDefinitionCell[5].find("trueofdatebcf") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::TrueOfDate_BCF;
                            }
                            else
                            {
                                throw std::invalid_argument("Currently only ICRF, J2000BCI, J2000BCF, TrueOfDateBCI, and TrueOfDateBCF frames are supported: velocity declination boundary constraint " + name);
                            }

                            return new BoundaryVelocityMagnitudeConstraint(name,
                                journeyIndex,
                                phaseIndex,
                                stageIndex,
                                Universe,
                                mySpacecraft,
                                myOptions,
                                user_specified_frame,
                                myBoundaryEvent,
                                ConstraintDefinition);
                        }
                        else
                        {
                            throw std::invalid_argument("Reference frame must be specified for velocity magnitude boundary constraint " + name);
                        }
                    }
                    else if (ConstraintDefinitionCell[2].find("sphericalazimuth") < 1024)
                    {
                        // if a frame is specified, then inform the constraint what the user wants
                        // otherwise terminate program
                        if (ConstraintDefinitionCell.size() > 5)
                        {
                            ReferenceFrame user_specified_frame;
                            if (ConstraintDefinitionCell[5].find("icrf") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::ICRF;
                            }
                            else if (ConstraintDefinitionCell[5].find("j2000bci") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::J2000_BCI;
                            }
                            else if (ConstraintDefinitionCell[5].find("j2000bcf") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::J2000_BCF;
                            }
                            else if (ConstraintDefinitionCell[5].find("trueofdatebci") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::TrueOfDate_BCI;
                            }
                            else if (ConstraintDefinitionCell[5].find("trueofdatebcf") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::TrueOfDate_BCF;
                            }
                            else
                            {
                                throw std::invalid_argument("Currently only ICRF, J2000BCI, J2000BCF, TrueOfDateBCI, and TrueOfDateBCF frames are supported: velocity spherical Azimuth boundary constraint " + name);
                            }

                            return new BoundaryVelocitySphericalAzimuthConstraint(name,
                                journeyIndex,
                                phaseIndex,
                                stageIndex,
                                Universe,
                                mySpacecraft,
                                myOptions,
                                user_specified_frame,
                                myBoundaryEvent,
                                ConstraintDefinition);
                        }
                        else
                        {
                            throw std::invalid_argument("Reference frame must be specified for velocity spherical Azimuth boundary constraint " + name);
                        }
                    }
                    else if (ConstraintDefinitionCell[2].find("relativevmagnitude") < 1024)
                    {
                        // if a frame is specified, then inform the constraint what the user wants
                        // otherwise terminate program
                        if (ConstraintDefinitionCell.size() > 5)
                        {
                            ReferenceFrame user_specified_frame;
                            if (ConstraintDefinitionCell[5].find("j2000bcf") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::J2000_BCF;
                            }
                            else if (ConstraintDefinitionCell[5].find("trueofdatebcf") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::TrueOfDate_BCF;
                            }
                            else
                            {
                                throw std::invalid_argument("Currently only J2000BCF and TrueOfDateBCF frames are supported: relative velocity magnitude " + name);
                            }

                            return new BoundaryRelativeVelocityMagnitudeConstraint(name,
                                journeyIndex,
                                phaseIndex,
                                stageIndex,
                                Universe,
                                mySpacecraft,
                                myOptions,
                                user_specified_frame,
                                myBoundaryEvent,
                                ConstraintDefinition);
                        }
                        else
                        {
                            throw std::invalid_argument("Reference frame must be specified for relative velocity magnitude boundary constraint " + name);
                        }
                    }
                    else if (ConstraintDefinitionCell[2].find("relativevazimuth") < 1024)
                    {
                        // if a frame is specified, then inform the constraint what the user wants
                        // otherwise terminate program
                        if (ConstraintDefinitionCell.size() > 5)
                        {
                            ReferenceFrame user_specified_frame;
                            if (ConstraintDefinitionCell[5].find("j2000bcf") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::J2000_BCF;
                            }
                            else if (ConstraintDefinitionCell[5].find("trueofdatebcf") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::TrueOfDate_BCF;
                            }
                            else
                            {
                                throw std::invalid_argument("Currently only J2000BCF and TrueOfDateBCF frames are supported: relative velocity Azimuth " + name);
                            }

                            return new BoundaryRelativeVelocityAzimuthConstraint(name,
                                journeyIndex,
                                phaseIndex,
                                stageIndex,
                                Universe,
                                mySpacecraft,
                                myOptions,
                                user_specified_frame,
                                myBoundaryEvent,
                                ConstraintDefinition);
                        }
                        else
                        {
                            throw std::invalid_argument("Reference frame must be specified for relative velocity Azimuth boundary constraint " + name);
                        }
                    }
                    else if (ConstraintDefinitionCell[2].find("relativevhfpa") < 1024)
                    {
                        // if a frame is specified, then inform the constraint what the user wants
                        // otherwise terminate program
                        if (ConstraintDefinitionCell.size() > 5)
                        {
                            ReferenceFrame user_specified_frame;
                            if (ConstraintDefinitionCell[5].find("j2000bcf") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::J2000_BCF;
                            }
                            else if (ConstraintDefinitionCell[5].find("trueofdatebcf") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::TrueOfDate_BCF;
                            }
                            else
                            {
                                throw std::invalid_argument("Currently only J2000BCF and TrueOfDateBCF frames are supported: relative velocity HFPA " + name);
                            }

                            return new BoundaryRelativeVelocityHFPAConstraint(name,
                                journeyIndex,
                                phaseIndex,
                                stageIndex,
                                Universe,
                                mySpacecraft,
                                myOptions,
                                user_specified_frame,
                                myBoundaryEvent,
                                ConstraintDefinition);
                        }
                        else
                        {
                            throw std::invalid_argument("Reference frame must be specified for relative velocity HFPA boundary constraint " + name);
                        }
                    }
                    else if (ConstraintDefinitionCell[2].find("vfpa") < 1024)
                    {
                        // if a frame is specified, then inform the constraint what the user wants
                        // otherwise terminate program
                        if (ConstraintDefinitionCell.size() > 5)
                        {
                            ReferenceFrame user_specified_frame;
                            if (ConstraintDefinitionCell[5].find("icrf") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::ICRF;
                            }
                            else if (ConstraintDefinitionCell[5].find("j2000bci") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::J2000_BCI;
                            }
                            else if (ConstraintDefinitionCell[5].find("j2000bcf") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::J2000_BCF;
                            }
                            else if (ConstraintDefinitionCell[5].find("trueofdatebci") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::TrueOfDate_BCI;
                            }
                            else if (ConstraintDefinitionCell[5].find("trueofdatebcf") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::TrueOfDate_BCF;
                            }
                            else
                            {
                                throw std::invalid_argument("Currently only ICRF, J2000BCI, J2000BCF, TrueOfDateBCI, and TrueOfDateBCF frames are supported: velocity declination boundary constraint " + name);
                            }

                            return new BoundaryVelocityVFPAconstraint(name,
                                journeyIndex,
                                phaseIndex,
                                stageIndex,
                                Universe,
                                mySpacecraft,
                                myOptions,
                                user_specified_frame,
                                myBoundaryEvent,
                                ConstraintDefinition);
                        }
                        else
                        {
                            throw std::invalid_argument("Reference frame must be specified for vertical flight path angle boundary constraint " + name);
                        }
                    }
                    else if (ConstraintDefinitionCell[2].find("angularmomentumreferenceangle") < 1024)
                    {
                        return new angularMomentumReferenceAngle(name,
                            journeyIndex,
                            phaseIndex,
                            stageIndex,
                            Universe,
                            mySpacecraft,
                            myOptions,
                            myBoundaryEvent,
                            ConstraintDefinition);
                    }
                    else
                    {
                        throw std::invalid_argument("Invalid boundary constraint, " + ConstraintDefinition);
                        return NULL;
                    }
                }//end departure
                else //arrival events
                {
                    if (ConstraintDefinitionCell[2].find("distanceconstraint") < 1024)
                    {
                        return new BoundaryDistanceConstraint(name,
                            journeyIndex,
                            phaseIndex,
                            stageIndex,
                            Universe,
                            mySpacecraft,
                            myOptions,
                            myBoundaryEvent,
                            ConstraintDefinition);
                    }
                    else if (ConstraintDefinitionCell[2].find("isapse") < 1024)
                    {
                        return new BoundaryIsApseConstraint(name,
                            journeyIndex,
                            phaseIndex,
                            stageIndex,
                            Universe,
                            mySpacecraft,
                            myOptions,
                            myBoundaryEvent,
                            ConstraintDefinition);
                    }
                    else if (ConstraintDefinitionCell[2].find("massconstraint") < 1024)
                    {
                        return new BoundaryMassConstraint(name,
                            journeyIndex,
                            phaseIndex,
                            stageIndex,
                            Universe,
                            mySpacecraft,
                            myOptions,
                            myBoundaryEvent,
                            ConstraintDefinition);
                    }
                    else if (ConstraintDefinitionCell[2].find("sma") < 1024)
                    {
                        return new SemimajorAxisConstraint(name,
                            journeyIndex,
                            phaseIndex,
                            stageIndex,
                            Universe,
                            mySpacecraft,
                            myOptions,
                            myBoundaryEvent,
                            ConstraintDefinition);
                    }
                    else if (ConstraintDefinitionCell[2].find("inc") < 1024)
                    {
                        // if a frame is specified, then inform the constraint what the user wants
                        // otherwise terminate program
                        if (ConstraintDefinitionCell.size() > 5)
                        {
                            ReferenceFrame user_specified_frame;
                            if (ConstraintDefinitionCell[5].find("icrf") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::ICRF;
                            }
                            else if (ConstraintDefinitionCell[5].find("j2000bci") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::J2000_BCI;
                            }
                            else
                            {
                                throw std::invalid_argument("Currently only ICRF and J2000BCI frames are supported: inclination boundary constraint " + name);
                            }

                            return new InclinationConstraint(name,
                                journeyIndex,
                                phaseIndex,
                                stageIndex,
                                Universe,
                                mySpacecraft,
                                myOptions,
                                user_specified_frame,
                                myBoundaryEvent,
                                ConstraintDefinition);
                        }
                        else
                        {
                            throw std::invalid_argument("Reference frame must be specified for inclination boundary constraint " + name);
                        }
                    }
                    else if (ConstraintDefinitionCell[2].find("ecc") < 1024)
                    {
                        return new EccentricityConstraint(name,
                            journeyIndex,
                            phaseIndex,
                            stageIndex,
                            Universe,
                            mySpacecraft,
                            myOptions,
                            myBoundaryEvent,
                            ConstraintDefinition);
                    }
                    else if (ConstraintDefinitionCell[2].find("raan") < 1024)
                    {
                        // if a frame is specified, then inform the constraint what the user wants
                        // otherwise terminate program
                        if (ConstraintDefinitionCell.size() > 5)
                        {
                            ReferenceFrame user_specified_frame;
                            if (ConstraintDefinitionCell[5].find("icrf") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::ICRF;
                            }
                            else if (ConstraintDefinitionCell[5].find("j2000bci") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::J2000_BCI;
                            }
                            else
                            {
                                throw std::invalid_argument("Currently only ICRF and J2000BCI frames are supported: RAAN boundary constraint " + name);
                            }

                            return new RAANConstraint(name,
                                journeyIndex,
                                phaseIndex,
                                stageIndex,
                                Universe,
                                mySpacecraft,
                                myOptions,
                                user_specified_frame,
                                myBoundaryEvent,
                                ConstraintDefinition);
                        }
                        else
                        {
                            throw std::invalid_argument("Reference frame must be specified for RAAN boundary constraint " + name);
                        }
                    }
                    else if (ConstraintDefinitionCell[2].find("aop") < 1024)
                    {
                        // if a frame is specified, then inform the constraint what the user wants
                        // otherwise terminate program
                        if (ConstraintDefinitionCell.size() > 5)
                        {
                            ReferenceFrame user_specified_frame;
                            if (ConstraintDefinitionCell[5].find("icrf") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::ICRF;
                            }
                            else if (ConstraintDefinitionCell[5].find("j2000bci") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::J2000_BCI;
                            }
                            else
                            {
                                throw std::invalid_argument("Currently only ICRF and J2000BCI frames are supported: AOP boundary constraint " + name);
                            }

                            return new AOPConstraint(name,
                                journeyIndex,
                                phaseIndex,
                                stageIndex,
                                Universe,
                                mySpacecraft,
                                myOptions,
                                user_specified_frame,
                                myBoundaryEvent,
                                ConstraintDefinition);
                        }
                        else
                        {
                            throw std::invalid_argument("Reference frame must be specified for AOP boundary constraint " + name);
                        }
                    }
                    else if (ConstraintDefinitionCell[2].find("trueanomaly") < 1024)
                    {
                        return new TrueAnomalyConstraint(name,
                            journeyIndex,
                            phaseIndex,
                            stageIndex,
                            Universe,
                            mySpacecraft,
                            myOptions,
                            myBoundaryEvent,
                            ConstraintDefinition);
                    }
                    else if (ConstraintDefinitionCell[2].find("orbitperiod") < 1024)
                    {
                        return new OrbitPeriodConstraint(name,
                            journeyIndex,
                            phaseIndex,
                            stageIndex,
                            Universe,
                            mySpacecraft,
                            myOptions,
                            myBoundaryEvent,
                            ConstraintDefinition);
                    }
                    else if (ConstraintDefinitionCell[2].find("orbitalenergy") < 1024)
                    {
                        return new BoundaryOrbitalEnergyConstraint(name,
                            journeyIndex,
                            phaseIndex,
                            stageIndex,
                            Universe,
                            mySpacecraft,
                            myOptions,
                            myBoundaryEvent,
                            ConstraintDefinition);
                    }
                    else if (ConstraintDefinitionCell[2].find("periapsedistance") < 1024)
                    {
                        return new PeriapseDistanceConstraint(name,
                            journeyIndex,
                            phaseIndex,
                            stageIndex,
                            Universe,
                            mySpacecraft,
                            myOptions,
                            myBoundaryEvent,
                            ConstraintDefinition);
                    }
                    else if (ConstraintDefinitionCell[2].find("apoapsedistance") < 1024)
                    {
                        return new ApoapseDistanceConstraint(name,
                            journeyIndex,
                            phaseIndex,
                            stageIndex,
                            Universe,
                            mySpacecraft,
                            myOptions,
                            myBoundaryEvent,
                            ConstraintDefinition);
                    }
                    else if (ConstraintDefinitionCell[2].find("bcflatitude") < 1024)
                    {
                        return new BoundaryBCFlatitudeConstraint(name,
                            journeyIndex,
                            phaseIndex,
                            stageIndex,
                            Universe,
                            mySpacecraft,
                            myOptions,
                            myBoundaryEvent,
                            ConstraintDefinition);
                    }
                    else if (ConstraintDefinitionCell[2].find("longitude") < 1024)
                    {
                        // if a frame is specified, then inform the constraint what the user wants
                        // otherwise terminate program
                        if (ConstraintDefinitionCell.size() > 5)
                        {
                            ReferenceFrame user_specified_frame;
                            if (ConstraintDefinitionCell[5].find("j2000bcf") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::J2000_BCF;
                            }
                            else if (ConstraintDefinitionCell[5].find("trueofdatebcf") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::TrueOfDate_BCF;
                            }
                            else if (ConstraintDefinitionCell[5].find("icrf") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::ICRF;
                            }
                            else
                            {
                                throw std::invalid_argument("Currently only ICRF, J2000BCF and TrueOfDateBCF frames are supported: longitude boundary constraint " + name);
                            }

                            return new BoundaryLongitudeConstraint(name,
                                journeyIndex,
                                phaseIndex,
                                stageIndex,
                                Universe,
                                mySpacecraft,
                                myOptions,
                                user_specified_frame,
                                myBoundaryEvent,
                                ConstraintDefinition);
                        }
                        else
                        {
                            throw std::invalid_argument("Reference frame must be specified for longitude boundary constraint " + name);
                        }
                    }
                    else if (ConstraintDefinitionCell[2].find("deticlatitude") < 1024)
                    {
                        // if a frame is specified, then inform the constraint what the user wants
                        // otherwise terminate program
                        if (ConstraintDefinitionCell.size() > 5)
                        {
                            ReferenceFrame user_specified_frame;
                            if (ConstraintDefinitionCell[5].find("j2000bcf") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::J2000_BCF;
                            }
                            else if (ConstraintDefinitionCell[5].find("trueofdatebcf") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::TrueOfDate_BCF;
                            }
                            else
                            {
                                throw std::invalid_argument("Currently only J2000BCF and TrueOfDateBCF frames are supported: detic latitude boundary constraint " + name);
                            }

                            return new BoundaryDeticLatitudeConstraint(name,
                                journeyIndex,
                                phaseIndex,
                                stageIndex,
                                Universe,
                                mySpacecraft,
                                myOptions,
                                user_specified_frame,
                                myBoundaryEvent,
                                ConstraintDefinition);
                        }
                        else
                        {
                            throw std::invalid_argument("Reference frame must be specified for detic latitude boundary constraint " + name);
                        }
                    }
                    else if (ConstraintDefinitionCell[2].find("deticaltitude") < 1024)
                    {
                        // if a frame is specified, then inform the constraint what the user wants
                        // otherwise terminate program
                        if (ConstraintDefinitionCell.size() > 5)
                        {
                            ReferenceFrame user_specified_frame;
                            if (ConstraintDefinitionCell[5].find("j2000bcf") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::J2000_BCF;
                            }
                            else if (ConstraintDefinitionCell[5].find("trueofdatebcf") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::TrueOfDate_BCF;
                            }
                            else
                            {
                                throw std::invalid_argument("Currently only J2000BCF and TrueOfDateBCF frames are supported: detic altitude boundary constraint " + name);
                            }

                            return new BoundaryDeticAltitudeConstraint(name,
                                journeyIndex,
                                phaseIndex,
                                stageIndex,
                                Universe,
                                mySpacecraft,
                                myOptions,
                                user_specified_frame,
                                myBoundaryEvent,
                                ConstraintDefinition);
                        }
                        else
                        {
                            throw std::invalid_argument("Reference frame must be specified for detic altitude boundary constraint " + name);
                        }
                    }
                    else if (ConstraintDefinitionCell[2].find("deticelevation") < 1024)
                    {
                        // if a frame is specified, then inform the constraint what the user wants
                        // otherwise terminate program
                        if (ConstraintDefinitionCell.size() > 6)
                        {
                            ReferenceFrame user_specified_frame;
                            if (ConstraintDefinitionCell[6].find("trueofdatebcf") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::TrueOfDate_BCF;
                            }
                            else
                            {
                                throw std::invalid_argument("Currently only TrueOfDateBCF frames are supported: detic elevation boundary constraint " + name);
                            }

                            return new BoundaryTargetBodyDeticElevationConstraint(name,
                                journeyIndex,
                                phaseIndex,
                                stageIndex,
                                Universe,
                                mySpacecraft,
                                myOptions,
                                user_specified_frame,
                                myBoundaryEvent,
                                ConstraintDefinition);
                        }
                        else
                        {
                            throw std::invalid_argument("Reference frame must be specified for detic elevation boundary constraint " + name);
                        }
                    }
                    else if (ConstraintDefinitionCell[2].find("velocitydeclination") < 1024)
                    {
                        // if a frame is specified, then inform the constraint what the user wants
                        // otherwise terminate program
                        if (ConstraintDefinitionCell.size() > 5)
                        {
                            ReferenceFrame user_specified_frame;
                            if (ConstraintDefinitionCell[5].find("icrf") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::ICRF;
                            }
                            else if (ConstraintDefinitionCell[5].find("j2000bci") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::J2000_BCI;
                            }
                            else if (ConstraintDefinitionCell[5].find("j2000bcf") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::J2000_BCF;
                            }
                            else if (ConstraintDefinitionCell[5].find("trueofdatebci") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::TrueOfDate_BCI;
                            }
                            else if (ConstraintDefinitionCell[5].find("trueofdatebcf") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::TrueOfDate_BCF;
                            }
                            else
                            {
                                throw std::invalid_argument("Currently only ICRF, J2000BCI, J2000BCF, TrueOfDateBCI, and TrueOfDateBCF frames are supported: velocity declination boundary constraint " + name);
                            }

                            return new BoundaryVelocityDeclinationConstraint(name,
                                journeyIndex,
                                phaseIndex,
                                stageIndex,
                                Universe,
                                mySpacecraft,
                                myOptions,
                                user_specified_frame,
                                myBoundaryEvent,
                                ConstraintDefinition);
                        }
                    }
                    else if (ConstraintDefinitionCell[2].find("velocitymagnitude") < 1024)
                    {
                        // if a frame is specified, then inform the constraint what the user wants
                        // otherwise terminate program
                        if (ConstraintDefinitionCell.size() > 5)
                        {
                            ReferenceFrame user_specified_frame;
                            if (ConstraintDefinitionCell[5].find("icrf") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::ICRF;
                            }
                            else if (ConstraintDefinitionCell[5].find("j2000bci") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::J2000_BCI;
                            }
                            else if (ConstraintDefinitionCell[5].find("j2000bcf") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::J2000_BCF;
                            }
                            else if (ConstraintDefinitionCell[5].find("trueofdatebci") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::TrueOfDate_BCI;
                            }
                            else if (ConstraintDefinitionCell[5].find("trueofdatebcf") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::TrueOfDate_BCF;
                            }
                            else
                            {
                                throw std::invalid_argument("Currently only ICRF, J2000BCI, J2000BCF, TrueOfDateBCI, and TrueOfDateBCF frames are supported: velocity declination boundary constraint " + name);
                            }

                            return new BoundaryVelocityMagnitudeConstraint(name,
                                journeyIndex,
                                phaseIndex,
                                stageIndex,
                                Universe,
                                mySpacecraft,
                                myOptions,
                                user_specified_frame,
                                myBoundaryEvent,
                                ConstraintDefinition);
                        }
                        else
                        {
                            throw std::invalid_argument("Reference frame must be specified for velocity magnitude boundary constraint " + name);
                        }
                    }
                    else if (ConstraintDefinitionCell[2].find("sphericalazimuth") < 1024)
                    {
                        // if a frame is specified, then inform the constraint what the user wants
                        // otherwise terminate program
                        if (ConstraintDefinitionCell.size() > 5)
                        {
                            ReferenceFrame user_specified_frame;
                            if (ConstraintDefinitionCell[5].find("icrf") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::ICRF;
                            }
                            else if (ConstraintDefinitionCell[5].find("j2000bci") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::J2000_BCI;
                            }
                            else if (ConstraintDefinitionCell[5].find("j2000bcf") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::J2000_BCF;
                            }
                            else if (ConstraintDefinitionCell[5].find("trueofdatebci") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::TrueOfDate_BCI;
                            }
                            else if (ConstraintDefinitionCell[5].find("trueofdatebcf") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::TrueOfDate_BCF;
                            }
                            else
                            {
                                throw std::invalid_argument("Currently only ICRF, J2000BCI, J2000BCF, TrueOfDateBCI, and TrueOfDateBCF frames are supported: velocity spherical Azimuth boundary constraint " + name);
                            }

                            return new BoundaryVelocitySphericalAzimuthConstraint(name,
                                journeyIndex,
                                phaseIndex,
                                stageIndex,
                                Universe,
                                mySpacecraft,
                                myOptions,
                                user_specified_frame,
                                myBoundaryEvent,
                                ConstraintDefinition);
                        }
                        else
                        {
                            throw std::invalid_argument("Reference frame must be specified for velocity spherical Azimuth boundary constraint " + name);
                        }
                    }
                    else if (ConstraintDefinitionCell[2].find("relativevmagnitude") < 1024)
                    {
                        // if a frame is specified, then inform the constraint what the user wants
                        // otherwise terminate program
                        if (ConstraintDefinitionCell.size() > 5)
                        {
                            ReferenceFrame user_specified_frame;
                            if (ConstraintDefinitionCell[5].find("j2000bcf") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::J2000_BCF;
                            }
                            else if (ConstraintDefinitionCell[5].find("trueofdatebcf") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::TrueOfDate_BCF;
                            }
                            else
                            {
                                throw std::invalid_argument("Currently only J2000BCF and TrueOfDateBCF frames are supported: relative velocity magnitude boundary constraint " + name);
                            }

                            return new BoundaryRelativeVelocityMagnitudeConstraint(name,
                                journeyIndex,
                                phaseIndex,
                                stageIndex,
                                Universe,
                                mySpacecraft,
                                myOptions,
                                user_specified_frame,
                                myBoundaryEvent,
                                ConstraintDefinition);
                        }
                        else
                        {
                            throw std::invalid_argument("Reference frame must be specified for relative velocity magnitude boundary constraint " + name);
                        }
                    }
                    else if (ConstraintDefinitionCell[2].find("relativevazimuth") < 1024)
                    {
                        // if a frame is specified, then inform the constraint what the user wants
                        // otherwise terminate program
                        if (ConstraintDefinitionCell.size() > 5)
                        {
                            ReferenceFrame user_specified_frame;
                            if (ConstraintDefinitionCell[5].find("j2000bcf") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::J2000_BCF;
                            }
                            else if (ConstraintDefinitionCell[5].find("trueofdatebcf") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::TrueOfDate_BCF;
                            }
                            else
                            {
                                throw std::invalid_argument("Currently only J2000BCF and TrueOfDateBCF frames are supported: relative velocity Azimuth " + name);
                            }

                            return new BoundaryRelativeVelocityAzimuthConstraint(name,
                                journeyIndex,
                                phaseIndex,
                                stageIndex,
                                Universe,
                                mySpacecraft,
                                myOptions,
                                user_specified_frame,
                                myBoundaryEvent,
                                ConstraintDefinition);
                        }
                        else
                        {
                            throw std::invalid_argument("Reference frame must be specified for relative velocity Azimuth boundary constraint " + name);
                        }
                    }
                    else if (ConstraintDefinitionCell[2].find("relativevhfpa") < 1024)
                    {
                        // if a frame is specified, then inform the constraint what the user wants
                        // otherwise terminate program
                        if (ConstraintDefinitionCell.size() > 5)
                        {
                            ReferenceFrame user_specified_frame;
                            if (ConstraintDefinitionCell[5].find("j2000bcf") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::J2000_BCF;
                            }
                            else if (ConstraintDefinitionCell[5].find("trueofdatebcf") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::TrueOfDate_BCF;
                            }
                            else
                            {
                                throw std::invalid_argument("Currently only J2000BCF and TrueOfDateBCF frames are supported: relative velocity HFPA " + name);
                            }

                            return new BoundaryRelativeVelocityHFPAConstraint(name,
                                journeyIndex,
                                phaseIndex,
                                stageIndex,
                                Universe,
                                mySpacecraft,
                                myOptions,
                                user_specified_frame,
                                myBoundaryEvent,
                                ConstraintDefinition);
                        }
                        else
                        {
                            throw std::invalid_argument("Reference frame must be specified for relative velocity HFPA boundary constraint " + name);
                        }
                    }
                    else if (ConstraintDefinitionCell[2].find("vfpa") < 1024)
                    {
                        // if a frame is specified, then inform the constraint what the user wants
                        // otherwise terminate program
                        if (ConstraintDefinitionCell.size() > 5)
                        {
                            ReferenceFrame user_specified_frame;
                            if (ConstraintDefinitionCell[5].find("icrf") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::ICRF;
                            }
                            else if (ConstraintDefinitionCell[5].find("j2000bci") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::J2000_BCI;
                            }
                            else if (ConstraintDefinitionCell[5].find("j2000bcf") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::J2000_BCF;
                            }
                            else if (ConstraintDefinitionCell[5].find("trueofdatebci") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::TrueOfDate_BCI;
                            }
                            else if (ConstraintDefinitionCell[5].find("trueofdatebcf") < 1024)
                            {
                                user_specified_frame = ReferenceFrame::TrueOfDate_BCF;
                            }
                            else
                            {
                                throw std::invalid_argument("Currently only ICRF, J2000BCI, J2000BCF, TrueOfDateBCI, and TrueOfDateBCF frames are supported: velocity declination boundary constraint " + name);
                            }

                            return new BoundaryVelocityVFPAconstraint(name,
                                journeyIndex,
                                phaseIndex,
                                stageIndex,
                                Universe,
                                mySpacecraft,
                                myOptions,
                                user_specified_frame,
                                myBoundaryEvent,
                                ConstraintDefinition);
                        }
                        else
                        {
                            throw std::invalid_argument("Reference frame must be specified for vertical flight path angle boundary constraint " + name);
                        }
                    }
                    else if (ConstraintDefinitionCell[2].find("rbp") < 1024)
                    {
                        if (myOptions->Journeys[journeyIndex].arrival_class == BoundaryClass::EphemerisPegged
                            && (myOptions->Journeys[journeyIndex].arrival_type == ArrivalType::INTERCEPT
                                || myOptions->Journeys[journeyIndex].arrival_type == ArrivalType::CHEM_RENDEZVOUS
                                || myOptions->Journeys[journeyIndex].arrival_type == ArrivalType::MOMENTUM_EXCHANGE
                                || myOptions->Journeys[journeyIndex].arrival_type == ArrivalType::INSERTION_INTO_PARKING))
                        {
                            return new RBPconstraint(name,
                                journeyIndex,
                                phaseIndex,
                                stageIndex,
                                Universe,
                                mySpacecraft,
                                myOptions,
                                myBoundaryEvent,
                                ConstraintDefinition);
                        }
                        else
                        {
                            throw std::invalid_argument("RBP constraint may only be applied to Ephemeris Pegged Intercept, Chemical Rendezvous, Momentum Exchange, or Insertion into Parking Orbit arrival events. Invalid constraint, " + ConstraintDefinition);
                        }
                    }
                    else if (ConstraintDefinitionCell[2].find("rpb") < 1024)
                    {
                        if (myOptions->Journeys[journeyIndex].arrival_class == BoundaryClass::EphemerisPegged
                            && (myOptions->Journeys[journeyIndex].arrival_type == ArrivalType::INTERCEPT
                                || myOptions->Journeys[journeyIndex].arrival_type == ArrivalType::CHEM_RENDEZVOUS
                                || myOptions->Journeys[journeyIndex].arrival_type == ArrivalType::MOMENTUM_EXCHANGE
                                || myOptions->Journeys[journeyIndex].arrival_type == ArrivalType::INSERTION_INTO_PARKING))
                        {
                            return new RPBconstraint(name,
                                journeyIndex,
                                phaseIndex,
                                stageIndex,
                                Universe,
                                mySpacecraft,
                                myOptions,
                                myBoundaryEvent,
                                ConstraintDefinition);
                        }
                        else
                        {
                            throw std::invalid_argument("RPB constraint may only be applied to Ephemeris Pegged Intercept, Chemical Rendezvous, Momentum Exchange, or Insertion into Parking Orbit arrival events. Invalid constraint, " + ConstraintDefinition);
                        }
                    }
                    else if (ConstraintDefinitionCell[2].find("rrp") < 1024)
                    {

                        return new RRPconstraint(name,
                            journeyIndex,
                            phaseIndex,
                            stageIndex,
                            Universe,
                            mySpacecraft,
                            myOptions,
                            myBoundaryEvent,
                            ConstraintDefinition);
                    }
                    else if (ConstraintDefinitionCell[2].find("rpr") < 1024)
                    {

                        return new RPRconstraint(name,
                            journeyIndex,
                            phaseIndex,
                            stageIndex,
                            Universe,
                            mySpacecraft,
                            myOptions,
                            myBoundaryEvent,
                            ConstraintDefinition);
                    }
                    else if (ConstraintDefinitionCell[2].find("angularmomentumreferenceangle") < 1024)
                    {
                        return new angularMomentumReferenceAngle(name,
                            journeyIndex,
                            phaseIndex,
                            stageIndex,
                            Universe,
                            mySpacecraft,
                            myOptions,
                            myBoundaryEvent,
                            ConstraintDefinition);
                    }
                    else
                    {
                        throw std::invalid_argument("Invalid boundary constraint, " + ConstraintDefinition);
                        return NULL;
                    }
                }//end arrival
                throw std::invalid_argument("You have reached the end of the specialized boundary constraint factory without either making a constraint or throwing an error for an invalid constraint. This should be impossible. Constraint is:" + ConstraintDefinition);
                return NULL;
            }

        }//end namespace SpecializedConstraints
    }//end namespace BoundaryEvents
}//end namespace EMTG