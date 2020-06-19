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

//event testbed class

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <random>
#include <vector>

#include "EphemerisPeggedDeparture/EphemerisPeggedLaunchDirectInsertion.h"
#include "EphemerisPeggedDeparture/EphemerisPeggedFreeDirectDeparture.h"
#include "EphemerisPeggedDeparture/EphemerisPeggedSpiralDeparture.h"

#include "EphemerisPeggedArrival/EphemerisPeggedIntercept.h"
#include "EphemerisPeggedArrival/EphemerisPeggedFlybyIn.h"
#include "EphemerisPeggedArrival/EphemerisPeggedLTRendezvous.h"
#include "EphemerisPeggedArrival/EphemerisPeggedChemRendezvous.h"

#include "FreePointDeparture/FreePointFreeDirectDeparture.h"
#include "FreePointDeparture/FreePointDirectInsertion.h"

#include "FreePointArrival/FreePointLTRendezvous.h"
#include "FreePointArrival/FreePointChemRendezvous.h"
#include "FreePointArrival/FreePointIntercept.h"

#include "EphemerisReferencedArrival/Exterior/EphemerisReferencedLTRendezvousExterior.h"
#include "EphemerisReferencedArrival/Exterior/EphemerisReferencedInterceptExterior.h"

#include "EphemerisReferencedArrival/Interior/EphemerisReferencedLTRendezvousInterior.h"
#include "EphemerisReferencedArrival/Interior/EphemerisReferencedInterceptInterior.h"

#include "PeriapseDeparture/PeriapseLaunch.h"

#include "PeriapseArrival/PeriapseFlybyIn.h"

#include "eventTestbed.h"

void eventTestbed(EMTG::missionoptions& options,
    std::vector< EMTG::Astrodynamics::universe > TheUniverse,
    EMTG::HardwareModels::Spacecraft& mySpacecraft,
    EMTG::HardwareModels::LaunchVehicle& myLaunchVehicle,
    std::mt19937 RNG,
    std::uniform_real_distribution<> UniformDouble)
{
    //Create the boundary/derivative vectors. They are empty until populated by the event objects.
    std::vector<double> Xupperbounds;
    std::vector<double> Xlowerbounds;
    std::vector<double> X_scale_factors;
    std::vector<double> Fupperbounds;
    std::vector<double> Flowerbounds;
	std::vector<double> F_scale_factors;
    std::vector<std::string> Xdescriptions;
    std::vector<std::string> Fdescriptions;
    std::vector<size_t> iGfun;
    std::vector<size_t> jGvar;
    std::vector<std::string> Gdescriptions;
    std::vector<size_t> iAfun;
    std::vector<size_t> jAvar;
    std::vector<std::string> Adescriptions;
    std::vector<double> A;

    //dummy phase pointer
    std::string phaseName("j0p0");
    size_t stageIndex = 0;
    std::string boundaryname;

    //********************************************************************************event setup
    //----------------------------------Ephemeris Pegged Spiral Departure event
    //name the event
    boundaryname = phaseName + "EphemerisPeggedSpiralDeparture";

    //create the event
    EMTG::BoundaryEvents::EphemerisPeggedSpiralDeparture myEphemerisPeggedSpiralDeparture(boundaryname,//name
        0,//journeyIndex
        0,//phaseIndex
        stageIndex,//stageIndex
        &TheUniverse.front(),//Universe
        &mySpacecraft,//spacecraft
        &options,//options
        NULL);//previous phase arrival event

    //setup bounds pointers
    myEphemerisPeggedSpiralDeparture.setup_calcbounds(&Xupperbounds,
        &Xlowerbounds,
        &X_scale_factors,
        &Fupperbounds,
        &Flowerbounds,
        &F_scale_factors,
        &Xdescriptions,
        &Fdescriptions,
        &iGfun,
        &jGvar,
        &Gdescriptions,
        &iAfun,
        &jAvar,
        &Adescriptions,
        &A);


    //call calcbounds
    myEphemerisPeggedSpiralDeparture.calcbounds(std::vector<size_t>({}));

    ////----------------------------------Periapse Launch event
    ////name the event
    //boundaryname = phaseName + "PeriapseLaunch";

    ////create the event
    //EMTG::BoundaryEvents::PeriapseLaunch myPeriapseLaunch(boundaryname,//name
    //    0,//journeyIndex
    //    0,//phaseIndex
    //    stageIndex,//stageIndex
    //    &TheUniverse.front(),//Universe
    //    &mySpacecraft,//spacecraft
    //    &myLaunchVehicle,//launch vehicle
    //    &options,//options
    //    NULL);//previous phase arrival event

    //          //setup bounds pointers
    //myPeriapseLaunch.setup_calcbounds(&Xupperbounds,
    //    &Xlowerbounds,
    //    &X_scale_factors,
    //    &Fupperbounds,
    //    &Flowerbounds,
    //    &Xdescriptions,
    //    &Fdescriptions,
    //    &iGfun,
    //    &jGvar,
    //    &Gdescriptions,
    //    &iAfun,
    //    &jAvar,
    //    &Adescriptions,
    //    &A);


    ////call calcbounds
    //myPeriapseLaunch.calcbounds();

  //  //----------------------------------Periapse flyby in event
  //  //name the event
  //  boundaryname = phaseName + "PeriapseFlybyIn";

  //  //create the event
  //  EMTG::BoundaryEvents::PeriapseFlybyIn myPeriapseFlybyIn(boundaryname,//name
  //      0,//journeyIndex
  //      0,//phaseIndex
  //      stageIndex,//stageIndex
  //      &TheUniverse.front(),//Universe
  //      &mySpacecraft,//spacecraft
  //      &options);

  //  //setup bounds pointers
  //  myPeriapseFlybyIn.setup_calcbounds(&Xupperbounds,
  //      &Xlowerbounds,
  //      &X_scale_factors,
  //      &Fupperbounds,
  //      &Flowerbounds,
		//&F_scale_factors,
  //      &Xdescriptions,
  //      &Fdescriptions,
  //      &iGfun,
  //      &jGvar,
  //      &Gdescriptions,
  //      &iAfun,
  //      &jAvar,
  //      &Adescriptions,
  //      &A);


  //  //call calcbounds
  //  myPeriapseFlybyIn.calcbounds();




  //  //----------------------------------Ephemeris referenced intercept interior event
  //  //name the event
  //  boundaryname = phaseName + "EphemerisReferencedInterceptInterior";

  //  //create the event
  //  EMTG::BoundaryEvents::EphemerisReferencedInterceptInterior myEphemerisReferencedInterceptInterior(boundaryname,//name
  //      0,//journeyIndex
  //      0,//phaseIndex
  //      stageIndex,//stageIndex
  //      &TheUniverse.front(),//Universe
  //      &mySpacecraft,//spacecraft
  //      &options);

  //  //setup bounds pointers
  //  myEphemerisReferencedInterceptInterior.setup_calcbounds(&Xupperbounds,
  //      &Xlowerbounds,
  //      &X_scale_factors,
  //      &Fupperbounds,
  //      &Flowerbounds,
		//&F_scale_factors,
  //      &Xdescriptions,
  //      &Fdescriptions,
  //      &iGfun,
  //      &jGvar,
  //      &Gdescriptions,
  //      &iAfun,
  //      &jAvar,
  //      &Adescriptions,
  //      &A);


  //  //call calcbounds
  //  myEphemerisReferencedInterceptInterior.calcbounds();

  //  //----------------------------------Free point direct insertion event
  //  //name the event
  //  boundaryname = phaseName + "FreePointDirectInsertion";

  //  //create the event
  //  EMTG::BoundaryEvents::FreePointDirectInsertion myFreePointDirectInsertion(boundaryname,//name
  //      0,//journeyIndex
  //      0,//phaseIndex
  //      stageIndex,//stageIndex
  //      &TheUniverse.front(),//Universe
  //      &mySpacecraft,//spacecraft
  //      &options,//options
  //      NULL);//previous phase arrival event

  //  //setup bounds pointers
  //  myFreePointDirectInsertion.setup_calcbounds(&Xupperbounds,
  //      &Xlowerbounds,
  //      &X_scale_factors,
  //      &Fupperbounds,
  //      &Flowerbounds,
		//&F_scale_factors,
  //      &Xdescriptions,
  //      &Fdescriptions,
  //      &iGfun,
  //      &jGvar,
  //      &Gdescriptions,
  //      &iAfun,
  //      &jAvar,
  //      &Adescriptions,
  //      &A);


  //  //call calcbounds
  //  myFreePointDirectInsertion.calcbounds();

    ////----------------------------------Ephemeris Pegged Free Direct Departure event
    ////name the event
    //boundaryname = phaseName + "EphemerisPeggedFreeDirectDeparture";

    ////create the event
    //EMTG::BoundaryEvents::EphemerisPeggedFreeDirectDeparture<GSAD::adouble> myEphemerisPeggedFreeDirectDeparture(boundaryname,//name
    //    0,//journeyIndex
    //    0,//phaseIndex
    //    stageIndex,//stageIndex
    //    &TheUniverse.front(),//Universe
    //    &mySpacecraft,//spacecraft
    //    &options,//options
    //    NULL);//previous phase arrival event

    ////setup bounds pointers
    //myEphemerisPeggedFreeDirectDeparture.setup_calcbounds(&Xupperbounds,
    //    &Xlowerbounds,
    //    &X_scale_factors,
    //    &Fupperbounds,
    //    &Flowerbounds,
    //    &Xdescriptions,
    //    &Fdescriptions,
    //    &iGfun,
    //    &jGvar,
    //    &Gdescriptions,
    //    &iAfun,
    //    &jAvar,
    //    &Adescriptions,
    //    &A);


    ////call calcbounds
    //myEphemerisPeggedFreeDirectDeparture.calcbounds();

  //  //----------------------------------Ephemeris Pegged LT Rendezvous event
  //  //name the event
  //  boundaryname = phaseName + "FreePointLTRendezvous";

  //  //create the event
  //  EMTG::BoundaryEvents::FreePointLTRendezvous myFreePointLTRendezvous(boundaryname,//name
  //      0,//journeyIndex
  //      0,//phaseIndex
  //      stageIndex,//stageIndex
  //      &TheUniverse.front(),//Universe
  //      &mySpacecraft,//spacecraft
  //      &options);//options

  //  //setup bounds pointers
  //  myFreePointLTRendezvous.setup_calcbounds(&Xupperbounds,
  //      &Xlowerbounds,
  //      &X_scale_factors,
  //      &Fupperbounds,
  //      &Flowerbounds,
		//&F_scale_factors,
  //      &Xdescriptions,
  //      &Fdescriptions,
  //      &iGfun,
  //      &jGvar,
  //      &Gdescriptions,
  //      &iAfun,
  //      &jAvar,
  //      &Adescriptions,
  //      &A);

  //  //call calcbounds
  //  myFreePointLTRendezvous.calcbounds();


    ////----------------------------------Ephemeris Pegged Intercept event
    ////name the event
    //boundaryname = phaseName + "EphemerisPeggedIntercept";


    //EMTG::BoundaryEvents::EphemerisPeggedIntercept<GSAD::adouble> myEphemerisPeggedIntercept(boundaryname,//name
    //    0,//journeyIndex
    //    0,//phaseIndex
    //    stageIndex,//stageIndex
    //    &TheUniverse.front(),//Universe
    //    &mySpacecraft,//spacecraft
    //    &options);//options

    //              //setup bounds pointers
    //myEphemerisPeggedIntercept.setup_calcbounds(&Xupperbounds,
    //    &Xlowerbounds,
    //    &X_scale_factors,
    //    &Fupperbounds,
    //    &Flowerbounds,
    //    &Xdescriptions,
    //    &Fdescriptions,
    //    &iGfun,
    //    &jGvar,
    //    &Gdescriptions,
    //    &iAfun,
    //    &jAvar,
    //    &Adescriptions,
    //    &A);


    ////call calcbounds
    //myEphemerisPeggedIntercept.calcbounds();

    ////----------------------------------Ephemeris Pegged FlybyIn event
    ////name the event
    //boundaryname = phaseName + "EphemerisPeggedFlybyIn";


    //EMTG::BoundaryEvents::EphemerisPeggedFlybyIn<GSAD::adouble> myEphemerisPeggedFlybyIn(boundaryname,//name
    //    0,//journeyIndex
    //    0,//phaseIndex
    //    stageIndex,//stageIndex
    //    &TheUniverse.front(),//Universe
    //    &mySpacecraft,//spacecraft
    //    &options);//options

    //              //setup bounds pointers
    //myEphemerisPeggedFlybyIn.setup_calcbounds(&Xupperbounds,
    //    &Xlowerbounds,
    //    &X_scale_factors,
    //    &Fupperbounds,
    //    &Flowerbounds,
    //    &Xdescriptions,
    //    &Fdescriptions,
    //    &iGfun,
    //    &jGvar,
    //    &Gdescriptions,
    //    &iAfun,
    //    &jAvar,
    //    &Adescriptions,
    //    &A);


    ////call calcbounds
    //myEphemerisPeggedFlybyIn.calcbounds();

    ////----------------------------------Ephemeris Pegged ChemRendezvous event
    ////name the event
    //boundaryname = phaseName + "EphemerisPeggedChemRendezvous";


    //EMTG::BoundaryEvents::EphemerisPeggedChemRendezvous myEphemerisPeggedChemRendezvous(boundaryname,//name
    //    0,//journeyIndex
    //    0,//phaseIndex
    //    stageIndex,//stageIndex
    //    &TheUniverse.front(),//Universe
    //    &mySpacecraft,//spacecraft
    //    &options);//options

    ////setup bounds pointers
    //myEphemerisPeggedChemRendezvous.setup_calcbounds(&Xupperbounds,
    //    &Xlowerbounds,
    //    &X_scale_factors,
    //    &Fupperbounds,
    //    &Flowerbounds,
    //    &Xdescriptions,
    //    &Fdescriptions,
    //    &iGfun,
    //    &jGvar,
    //    &Gdescriptions,
    //    &iAfun,
    //    &jAvar,
    //    &Adescriptions,
    //    &A);


    ////call calcbounds
    //myEphemerisPeggedChemRendezvous.calcbounds();


    //pass scale factors to the Universes
    for (EMTG::Astrodynamics::universe& thisUniverse : TheUniverse)
        thisUniverse.set_X_scale_factors(&X_scale_factors);

    //********************************************************************************create decision and constraint vectors

    std::vector<GSAD::adouble> F(Fdescriptions.size(), 0.0);
    std::vector<double> G(jGvar.size(), 1.0e+100);
    std::vector<GSAD::adouble> X;
    std::vector< EMTG::math::Matrix<GSAD::adouble> > StateBeforeEvent;
    std::vector< EMTG::math::Matrix<GSAD::adouble> > StateAfterEvent;
    std::vector< std::tuple<size_t, size_t, double> > derivative_container1;
    std::vector< std::tuple<size_t, size_t, double> > derivative_container2;
    std::vector< std::vector< std::tuple<size_t, size_t, double> > > dStateBeforeEvent_dX;
    std::vector< std::vector< std::tuple<size_t, size_t, double> > > dStateAfterEvent_dX;
    std::vector<std::string> eventNames;

    
    for (size_t Xindex = 0; Xindex < Xdescriptions.size(); ++Xindex)
    {
        std::uniform_real_distribution<double> myDistribution(Xlowerbounds[Xindex] / X_scale_factors[Xindex],
            Xupperbounds[Xindex] / X_scale_factors[Xindex]);
        X.push_back(myDistribution(RNG));
        X[Xindex].setDerivative(Xindex, 1.0);
        X[Xindex] = X[Xindex] * X_scale_factors[Xindex];

        if (X.back() < Xlowerbounds[Xindex])
            X.back() = Xlowerbounds[Xindex] + 1.0e-8;
        else if (X.back() > Xupperbounds[Xindex])
            X.back() = Xupperbounds[Xindex] - 1.0e-8;
    }

    //evaluate!
    size_t Xindex = 0;
    size_t Findex = 0;
    myEphemerisPeggedSpiralDeparture.process_event(X, Xindex, F, Findex, G, true);
    eventNames.push_back(myEphemerisPeggedSpiralDeparture.getName());
    StateBeforeEvent.push_back(myEphemerisPeggedSpiralDeparture.get_state_before_event());
    derivative_container1 = myEphemerisPeggedSpiralDeparture.get_Derivatives_of_StateBeforeEvent();
    derivative_container2 = myEphemerisPeggedSpiralDeparture.get_Derivatives_of_StateBeforeEvent_wrt_Time();
    derivative_container1.insert(derivative_container1.end(), derivative_container2.begin(), derivative_container2.end());
    dStateBeforeEvent_dX.push_back(derivative_container1);
    StateAfterEvent.push_back(myEphemerisPeggedSpiralDeparture.get_state_after_event());
    derivative_container1 = myEphemerisPeggedSpiralDeparture.get_Derivatives_of_StateAfterEvent();
    derivative_container2 = myEphemerisPeggedSpiralDeparture.get_Derivatives_of_StateAfterEvent_wrt_Time();
    derivative_container1.insert(derivative_container1.end(), derivative_container2.begin(), derivative_container2.end());
    dStateAfterEvent_dX.push_back(derivative_container1);

    /*myPeriapseLaunch.process_event(X, Xindex, F, Findex, G, true);
    eventNames.push_back(myPeriapseLaunch.getName());
    StateBeforeEvent.push_back(myPeriapseLaunch.get_state_before_event());
    derivative_container1 = myPeriapseLaunch.get_Derivatives_of_StateBeforeEvent();
    derivative_container2 = myPeriapseLaunch.get_Derivatives_of_StateBeforeEvent_wrt_Time();
    derivative_container1.insert(derivative_container1.end(), derivative_container2.begin(), derivative_container2.end());
    dStateBeforeEvent_dX.push_back(derivative_container1);
    StateAfterEvent.push_back(myPeriapseLaunch.get_state_after_event());
    derivative_container1 = myPeriapseLaunch.get_Derivatives_of_StateAfterEvent();
    derivative_container2 = myPeriapseLaunch.get_Derivatives_of_StateAfterEvent_wrt_Time();
    derivative_container1.insert(derivative_container1.end(), derivative_container2.begin(), derivative_container2.end());
    dStateAfterEvent_dX.push_back(derivative_container1);*/

    /*myPeriapseFlybyIn.process_event(X, Xindex, F, Findex, G, true);
    eventNames.push_back(myPeriapseFlybyIn.getName());
    StateBeforeEvent.push_back(myPeriapseFlybyIn.get_state_before_event());
    derivative_container1 = myPeriapseFlybyIn.get_Derivatives_of_StateBeforeEvent();
    derivative_container2 = myPeriapseFlybyIn.get_Derivatives_of_StateBeforeEvent_wrt_Time();
    derivative_container1.insert(derivative_container1.end(), derivative_container2.begin(), derivative_container2.end());
    dStateBeforeEvent_dX.push_back(derivative_container1);
    StateAfterEvent.push_back(myPeriapseFlybyIn.get_state_after_event());
    derivative_container1 = myPeriapseFlybyIn.get_Derivatives_of_StateAfterEvent();
    derivative_container2 = myPeriapseFlybyIn.get_Derivatives_of_StateAfterEvent_wrt_Time();
    derivative_container1.insert(derivative_container1.end(), derivative_container2.begin(), derivative_container2.end());
    dStateAfterEvent_dX.push_back(derivative_container1);

    myEphemerisReferencedInterceptInterior.process_event(X, Xindex, F, Findex, G, true);
    eventNames.push_back(myEphemerisReferencedInterceptInterior.getName());
    StateBeforeEvent.push_back(myEphemerisReferencedInterceptInterior.get_state_before_event());
    derivative_container1 = myEphemerisReferencedInterceptInterior.get_Derivatives_of_StateBeforeEvent();
    derivative_container2 = myEphemerisReferencedInterceptInterior.get_Derivatives_of_StateBeforeEvent_wrt_Time();
    derivative_container1.insert(derivative_container1.end(), derivative_container2.begin(), derivative_container2.end());
    dStateBeforeEvent_dX.push_back(derivative_container1);
    StateAfterEvent.push_back(myEphemerisReferencedInterceptInterior.get_state_after_event());
    derivative_container1 = myEphemerisReferencedInterceptInterior.get_Derivatives_of_StateAfterEvent();
    derivative_container2 = myEphemerisReferencedInterceptInterior.get_Derivatives_of_StateAfterEvent_wrt_Time();
    derivative_container1.insert(derivative_container1.end(), derivative_container2.begin(), derivative_container2.end());
    dStateAfterEvent_dX.push_back(derivative_container1);*/

 /*   myFreePointDirectInsertion.process_event(X, Xindex, F, Findex, G, true);
    eventNames.push_back(myFreePointDirectInsertion.getName());
    StateBeforeEvent.push_back(myFreePointDirectInsertion.get_state_before_event());
    derivative_container1 = myFreePointDirectInsertion.get_Derivatives_of_StateBeforeEvent();
    derivative_container2 = myFreePointDirectInsertion.get_Derivatives_of_StateBeforeEvent_wrt_Time();
    derivative_container1.insert(derivative_container1.end(), derivative_container2.begin(), derivative_container2.end());
    dStateBeforeEvent_dX.push_back(derivative_container1);
    StateAfterEvent.push_back(myFreePointDirectInsertion.get_state_after_event());
    derivative_container1 = myFreePointDirectInsertion.get_Derivatives_of_StateAfterEvent();
    derivative_container2 = myFreePointDirectInsertion.get_Derivatives_of_StateAfterEvent_wrt_Time();
    derivative_container1.insert(derivative_container1.end(), derivative_container2.begin(), derivative_container2.end());
    dStateAfterEvent_dX.push_back(derivative_container1);
*/
    /*
    myEphemerisPeggedFreeDirectDeparture.process_event(X, Xindex, F, Findex, G, true);
    eventNames.push_back(myEphemerisPeggedFreeDirectDeparture.getName());
    StateBeforeEvent.push_back(myEphemerisPeggedFreeDirectDeparture.get_state_before_event());
    dStateBeforeEvent_dX.push_back(myEphemerisPeggedFreeDirectDeparture.get_Derivatives_of_StateBeforeEvent());
    StateAfterEvent.push_back(myEphemerisPeggedFreeDirectDeparture.get_state_after_event());
    dStateAfterEvent_dX.push_back(myEphemerisPeggedFreeDirectDeparture.get_Derivatives_of_StateAfterEvent());*/

    //myFreePointLTRendezvous.process_event(X, Xindex, F, Findex, G, true);
    //eventNames.push_back(myFreePointLTRendezvous.getName());
    //StateBeforeEvent.push_back(myFreePointLTRendezvous.get_state_before_event());
    //derivative_container1 = myFreePointLTRendezvous.get_Derivatives_of_StateBeforeEvent();
    //derivative_container2 = myFreePointLTRendezvous.get_Derivatives_of_StateBeforeEvent_wrt_Time();
    //derivative_container1.insert(derivative_container1.end(), derivative_container2.begin(), derivative_container2.end());
    //dStateBeforeEvent_dX.push_back(derivative_container1);
    //StateAfterEvent.push_back(myFreePointLTRendezvous.get_state_after_event());
    //derivative_container1 = myFreePointLTRendezvous.get_Derivatives_of_StateAfterEvent();
    //derivative_container2 = myFreePointLTRendezvous.get_Derivatives_of_StateAfterEvent_wrt_Time();
    //derivative_container1.insert(derivative_container1.end(), derivative_container2.begin(), derivative_container2.end());
    //dStateAfterEvent_dX.push_back(derivative_container1);

    /*myEphemerisPeggedIntercept.process_event(X, Xindex, F, Findex, G, true);
    eventNames.push_back(myEphemerisPeggedIntercept.getName());
    StateBeforeEvent.push_back(myEphemerisPeggedIntercept.get_state_before_event());
    dStateBeforeEvent_dX.push_back(myEphemerisPeggedIntercept.get_Derivatives_of_StateBeforeEvent());
    StateAfterEvent.push_back(myEphemerisPeggedIntercept.get_state_after_event());
    dStateAfterEvent_dX.push_back(myEphemerisPeggedIntercept.get_Derivatives_of_StateAfterEvent());

    myEphemerisPeggedFlybyIn.process_event(X, Xindex, F, Findex, G, true);
    eventNames.push_back(myEphemerisPeggedFlybyIn.getName());
    StateBeforeEvent.push_back(myEphemerisPeggedFlybyIn.get_state_before_event());
    dStateBeforeEvent_dX.push_back(myEphemerisPeggedFlybyIn.get_Derivatives_of_StateBeforeEvent());
    StateAfterEvent.push_back(myEphemerisPeggedFlybyIn.get_state_after_event());
    dStateAfterEvent_dX.push_back(myEphemerisPeggedFlybyIn.get_Derivatives_of_StateAfterEvent());*/
/*
    myEphemerisPeggedChemRendezvous.process_event(X, Xindex, F, Findex, G, true);
    eventNames.push_back(myEphemerisPeggedChemRendezvous.getName());
    StateBeforeEvent.push_back(myEphemerisPeggedChemRendezvous.get_state_before_event());
    dStateBeforeEvent_dX.push_back(myEphemerisPeggedChemRendezvous.get_Derivatives_of_StateBeforeEvent());
    StateAfterEvent.push_back(myEphemerisPeggedChemRendezvous.get_state_after_event());
    dStateAfterEvent_dX.push_back(myEphemerisPeggedChemRendezvous.get_Derivatives_of_StateAfterEvent());*/


    //print XF
    std::ofstream XFout("tests/EventXFout.csv", std::ios::trunc);
    XFout << "index, Name, Lowerbound, Upperbound, Scale, Value" << std::endl;
    for (size_t Xindex = 0; Xindex < Xdescriptions.size(); ++Xindex)
        XFout << Xindex << "," << Xdescriptions[Xindex] << "," << Xlowerbounds[Xindex] << "," << Xupperbounds[Xindex] << "," << X_scale_factors[Xindex] << "," << X[Xindex].getValue() << std::endl;
    XFout << "index, Name, Lowerbound, Upperbound, Value" << std::endl;
    for (size_t Findex = 0; Findex < Fdescriptions.size(); ++Findex)
        XFout << Findex << "," << Fdescriptions[Findex] << "," << Flowerbounds[Findex] << "," << Fupperbounds[Findex] << "," << F[Findex].getValue() << std::endl;
    XFout.close();
    
    //print (sparse) Jacobian
    std::ofstream Gout("tests/EventGout.csv", std::ios::trunc);
    Gout << "index, iGfun, jGvar, Name, Analytical, Algorithmic, Absolute error, Relative error, an/al, al/an" << std::endl;
    for (size_t Gindex = 0; Gindex < Gdescriptions.size(); ++Gindex)
    {
        size_t Findex = iGfun[Gindex];
        size_t Xindex = jGvar[Gindex];
        Gout << Gindex << "," << Findex << "," << Xindex << "," << Gdescriptions[Gindex];
        double Galgorithmic = F[Findex].getDerivative(Xindex);
        double Ganalytical = G[Gindex];
        double abserror = Ganalytical - Galgorithmic;
        double relerror = abserror / Galgorithmic;
        double an_al = Ganalytical / Galgorithmic;
        double al_an = Galgorithmic / Ganalytical;
        Gout << "," << Ganalytical << "," << Galgorithmic << "," << abserror << "," << relerror << "," << an_al << "," << al_an << std::endl;
    }
    Gout.close();


    //search for missing Jacobian entries
    std::vector< std::tuple<size_t, size_t, double> > missingDerivatives; //Findex, Xindex, magnitude
    for (size_t Findex = 0; Findex < Fdescriptions.size(); ++Findex)
    {
        std::vector<size_t> derivativeIndices = F[Findex].getDerivativeIndicies();

        for (size_t Xindex : derivativeIndices)
        {
            bool found = false;
            for (size_t Gindex = 0; Gindex < Gdescriptions.size(); ++Gindex)
            {
                if (jGvar[Gindex] == Xindex && iGfun[Gindex] == Findex)
                {
                    found = true;
                    break;
                }
            }

            if (!found)
            {
                double magnitude = F[Findex].getDerivative(Xindex);
                missingDerivatives.push_back({ Findex, Xindex, magnitude });
            }
        }
    }

    std::ofstream Mout("tests/Event_MissingEntries.csv", std::ios::trunc);
    Mout << "Findex, Xindex, description, magnitude" << std::endl;

    for (std::tuple<size_t, size_t, double> missingDerivative : missingDerivatives)
    {
        size_t Findex = std::get<0>(missingDerivative);
        size_t Xindex = std::get<1>(missingDerivative);
        Mout << Findex << "," << Xindex << ",";
        Mout << "Derivative of " << Fdescriptions[Findex] << " F[" << Findex << "] with respect to X[" << Xindex << "]: " << Xdescriptions[Xindex];
        Mout << "," << std::get<2>(missingDerivative);
        Mout << std::endl;
    }
    Mout.close();

    //print derivatives of state before event
    std::vector<std::string> stateNames(std::vector<std::string>({ "x", "y", "z", "xdot", "ydot", "zdot", "mass", "epoch" }));
    std::ofstream StateBeforeEventOut("tests/EventStateBeforeEvent.csv", std::ios::trunc);
    StateBeforeEventOut << "Name, Analytical, Algorithmic, Absolute error, Relative error, an/al, al/an" << std::endl;
    StateBeforeEventOut.precision(15);
    for (size_t eventIndex = 0; eventIndex < dStateBeforeEvent_dX.size(); ++eventIndex)
    {
        for (size_t dIndex = 0; dIndex < dStateBeforeEvent_dX[eventIndex].size(); ++dIndex)
        {
            size_t Xindex = std::get<0>(dStateBeforeEvent_dX[eventIndex][dIndex]);
            size_t stateIndex = std::get<1>(dStateBeforeEvent_dX[eventIndex][dIndex]);
            std::string Name = "Derivative of " + eventNames[eventIndex] + " state_before_event(" + std::to_string(stateIndex) + "): " + stateNames[stateIndex] + " with respect to " + Xdescriptions[Xindex];
            double dAnalytical = std::get<2>(dStateBeforeEvent_dX[eventIndex][dIndex]) * X_scale_factors[Xindex];
            double dAlgorithmic = StateBeforeEvent[eventIndex](stateIndex).getDerivative(Xindex);
            double abserror = dAnalytical - dAlgorithmic;
            double relerror = abserror / dAlgorithmic;
            double an_al = dAnalytical / dAlgorithmic;
            double al_an = dAlgorithmic / dAnalytical;
            StateBeforeEventOut << Name << "," << dAnalytical << "," << dAlgorithmic << "," << abserror << "," << relerror << "," << an_al << "," << al_an << std::endl;
        }
    }
    //print derivatives of state after event
    std::ofstream StateAfterEventOut("tests/EventStateAfterEvent.csv", std::ios::trunc);
    StateAfterEventOut << "Name, Analytical, Algorithmic, Absolute error, Relative error, an/al, al/an" << std::endl;
    StateAfterEventOut.precision(15);
    for (size_t eventIndex = 0; eventIndex < dStateAfterEvent_dX.size(); ++eventIndex)
    {
        for (size_t dIndex = 0; dIndex < dStateAfterEvent_dX[eventIndex].size(); ++dIndex)
        {
            size_t Xindex = std::get<0>(dStateAfterEvent_dX[eventIndex][dIndex]);
            size_t stateIndex = std::get<1>(dStateAfterEvent_dX[eventIndex][dIndex]);
            std::string Name = "Derivative of " + eventNames[eventIndex] + " state_after_event(" + std::to_string(stateIndex) + "): " + stateNames[stateIndex] + " with respect to " + Xdescriptions[Xindex];
            double dAnalytical = std::get<2>(dStateAfterEvent_dX[eventIndex][dIndex]) * X_scale_factors[Xindex];
            double dAlgorithmic = StateAfterEvent[eventIndex](stateIndex).getDerivative(Xindex);
            double abserror = dAnalytical - dAlgorithmic;
            double relerror = abserror / dAlgorithmic;
            double an_al = dAnalytical / dAlgorithmic;
            double al_an = dAlgorithmic / dAnalytical;
            StateAfterEventOut << Name << "," << dAnalytical << "," << dAlgorithmic << "," << abserror << "," << relerror << "," << an_al << "," << al_an << std::endl;
        }
    }

    //print the event - we have to at least see if printing works
    std::ofstream outputcheck("tests/Event_outputcheck.emtg", std::ios::trunc);
    size_t eventcount = 1;

    //next, column headers
    {

        //column headers line 1
        outputcheck.width(5); outputcheck << "#";
        outputcheck.width(3); outputcheck << " | ";
        outputcheck.width(16); outputcheck << "JulianDate";
        outputcheck.width(3); outputcheck << " | ";
        outputcheck.width(11); outputcheck << "MM/DD/YYYY";
        outputcheck.width(3); outputcheck << " | ";
        outputcheck.width(12); outputcheck << "event type";
        outputcheck.width(3); outputcheck << " | ";
        outputcheck.width(25); outputcheck << "location";
        outputcheck.width(3); outputcheck << " | ";
        outputcheck.width(15); outputcheck << "step size";
        outputcheck.width(3); outputcheck << " | ";
        outputcheck.width(19); outputcheck << "altitude";
        outputcheck.width(3); outputcheck << " | ";
        outputcheck.width(19); outputcheck << "BdotR";
        outputcheck.width(3); outputcheck << " | ";
        outputcheck.width(19); outputcheck << "BdotT";
        outputcheck.width(3); outputcheck << " | ";
        outputcheck.width(8); outputcheck << "RA";
        outputcheck.width(3); outputcheck << " | ";
        outputcheck.width(8); outputcheck << "DEC";
        outputcheck.width(3); outputcheck << " | ";
        outputcheck.width(14); outputcheck << "C3";
        outputcheck.width(3); outputcheck << " | ";
        outputcheck.width(19); outputcheck << " x";
        outputcheck.width(3); outputcheck << " | ";
        outputcheck.width(19); outputcheck << " y";
        outputcheck.width(3); outputcheck << " | ";
        outputcheck.width(19); outputcheck << " z";
        outputcheck.width(3); outputcheck << " | ";
        outputcheck.width(19); outputcheck << " xdot";
        outputcheck.width(3); outputcheck << " | ";
        outputcheck.width(19); outputcheck << " ydot";
        outputcheck.width(3); outputcheck << " | ";
        outputcheck.width(19); outputcheck << " zdot";
        outputcheck.width(3); outputcheck << " | ";
        outputcheck.width(19); outputcheck << " dV_x";
        outputcheck.width(3); outputcheck << " | ";
        outputcheck.width(19); outputcheck << " dV_y";
        outputcheck.width(3); outputcheck << " | ";
        outputcheck.width(19); outputcheck << " dV_z";
        outputcheck.width(3); outputcheck << " | ";
        outputcheck.width(19); outputcheck << " T_x";
        outputcheck.width(3); outputcheck << " | ";
        outputcheck.width(19); outputcheck << " T_y";
        outputcheck.width(3); outputcheck << " | ";
        outputcheck.width(19); outputcheck << " T_z";
        outputcheck.width(3); outputcheck << " | ";
        outputcheck.width(17); outputcheck << "|dV| (km/s)";
        outputcheck.width(3); outputcheck << " | ";
        outputcheck.width(14); outputcheck << "Avail. Thrust";
        outputcheck.width(3); outputcheck << " | ";
        outputcheck.width(14); outputcheck << "Isp";
        outputcheck.width(3); outputcheck << " | ";
        outputcheck.width(14); outputcheck << "Avail. Power";
        outputcheck.width(3); outputcheck << " | ";
        outputcheck.width(19); outputcheck << "Mass Flow";
        outputcheck.width(3); outputcheck << " | ";
        outputcheck.width(14); outputcheck << "mass";
        outputcheck.width(3); outputcheck << " | ";
        outputcheck.width(14); outputcheck << "number of";
        outputcheck.width(3); outputcheck << " | ";
        outputcheck.width(14); outputcheck << "active power";
        outputcheck << std::endl;

        //column headers line 2
        outputcheck.width(5); outputcheck << "";
        outputcheck.width(3); outputcheck << " | ";
        outputcheck.width(16); outputcheck << " (ET)";
        outputcheck.width(3); outputcheck << " | ";
        outputcheck.width(11); outputcheck << "";
        outputcheck.width(3); outputcheck << " | ";
        outputcheck.width(12); outputcheck << "";
        outputcheck.width(3); outputcheck << " | ";
        outputcheck.width(25); outputcheck << "";
        outputcheck.width(3); outputcheck << " | ";
        outputcheck.width(15); outputcheck << "(days)";
        outputcheck.width(3); outputcheck << " | ";
        outputcheck.width(19); outputcheck << "(km)";
        outputcheck.width(3); outputcheck << " | ";
        outputcheck.width(19); outputcheck << "(km)";
        outputcheck.width(3); outputcheck << " | ";
        outputcheck.width(19); outputcheck << "(km)";
        outputcheck.width(3); outputcheck << " | ";
        outputcheck.width(8); outputcheck << "degrees";
        outputcheck.width(3); outputcheck << " | ";
        outputcheck.width(8); outputcheck << "degrees";
        outputcheck.width(3); outputcheck << " | ";
        outputcheck.width(14); outputcheck << "(km^2/s^2)";
        outputcheck.width(3); outputcheck << " | ";
        outputcheck.width(19); outputcheck << "(km)";
        outputcheck.width(3); outputcheck << " | ";
        outputcheck.width(19); outputcheck << "(km)";
        outputcheck.width(3); outputcheck << " | ";
        outputcheck.width(19); outputcheck << "(km)";
        outputcheck.width(3); outputcheck << " | ";
        outputcheck.width(19); outputcheck << "(km/s)";
        outputcheck.width(3); outputcheck << " | ";
        outputcheck.width(19); outputcheck << "(km/s)";
        outputcheck.width(3); outputcheck << " | ";
        outputcheck.width(19); outputcheck << "(km/s)";
        outputcheck.width(3); outputcheck << " | ";
        outputcheck.width(19); outputcheck << "(km/s)";
        outputcheck.width(3); outputcheck << " | ";
        outputcheck.width(19); outputcheck << "(km/s)";
        outputcheck.width(3); outputcheck << " | ";
        outputcheck.width(19); outputcheck << "(km/s)";
        outputcheck.width(3); outputcheck << " | ";
        outputcheck.width(19); outputcheck << "(N)";
        outputcheck.width(3); outputcheck << " | ";
        outputcheck.width(19); outputcheck << "(N)";
        outputcheck.width(3); outputcheck << " | ";
        outputcheck.width(19); outputcheck << "(N)";
        outputcheck.width(3); outputcheck << " | ";
        outputcheck.width(17); outputcheck << "throttle (0-1)";
        outputcheck.width(3); outputcheck << " | ";
        outputcheck.width(14); outputcheck << "(N)";
        outputcheck.width(3); outputcheck << " | ";
        outputcheck.width(14); outputcheck << "(s)";
        outputcheck.width(3); outputcheck << " | ";
        outputcheck.width(14); outputcheck << "(kW)";
        outputcheck.width(3); outputcheck << " | ";
        outputcheck.width(19); outputcheck << "rate (kg/s)";
        outputcheck.width(3); outputcheck << " | ";
        outputcheck.width(14); outputcheck << "(kg)";
        outputcheck.width(3); outputcheck << " | ";
        outputcheck.width(14); outputcheck << "active engines";
        outputcheck.width(3); outputcheck << " | ";
        outputcheck.width(14); outputcheck << "(kW)";
        outputcheck << std::endl;


        for (size_t k = 0; k < 615; ++k)
            outputcheck << "-";
        outputcheck << std::endl;
    }

    outputcheck << "Testing " << myEphemerisPeggedSpiralDeparture.getName() << std::endl;
    outputcheck << std::endl;
    myEphemerisPeggedSpiralDeparture.output(outputcheck, options.launch_window_open_date, eventcount);
    outputcheck << std::endl;
/*
    outputcheck << "Testing " << myPeriapseLaunch.getName() << std::endl;
    outputcheck << std::endl;
    myPeriapseLaunch.output(outputcheck, options.launch_window_open_date, eventcount);
    outputcheck << std::endl;*/

    /*outputcheck << "Testing " << myPeriapseFlybyIn.getName() << std::endl;
    outputcheck << std::endl;
    myPeriapseFlybyIn.output(outputcheck, options.launch_window_open_date, eventcount);
    outputcheck << std::endl;

    outputcheck << "Testing " << myEphemerisReferencedInterceptInterior.getName() << std::endl;
    outputcheck << std::endl;
    myEphemerisReferencedInterceptInterior.output(outputcheck, options.launch_window_open_date, eventcount);
    outputcheck << std::endl;*/

    /*outputcheck << "Testing " << myFreePointDirectInsertion.getName() << std::endl;
    outputcheck << std::endl;
    myFreePointDirectInsertion.output(outputcheck, options.launch_window_open_date, eventcount);
    outputcheck << std::endl;*/

    /*
    outputcheck << "Testing " << myEphemerisPeggedFreeDirectDeparture.getName() << std::endl;
    outputcheck << std::endl;
    myEphemerisPeggedFreeDirectDeparture.output(outputcheck, options.launch_window_open_date, eventcount);
    outputcheck << std::endl;*/

   /* outputcheck << "Testing " << myFreePointLTRendezvous.getName() << std::endl;
    outputcheck << std::endl;
    myFreePointLTRendezvous.output(outputcheck, options.launch_window_open_date, eventcount);
    outputcheck << std::endl;*/

    /*outputcheck << "Testing " << myEphemerisPeggedIntercept.getName() << std::endl;
    outputcheck << std::endl;
    myEphemerisPeggedIntercept.output(outputcheck, options.launch_window_open_date, eventcount);
    outputcheck << std::endl;

    outputcheck << "Testing " << myEphemerisPeggedFlybyIn.getName() << std::endl;
    outputcheck << std::endl;
    myEphemerisPeggedFlybyIn.output(outputcheck, options.launch_window_open_date, eventcount);
    outputcheck << std::endl;*/
/*
    outputcheck << "Testing " << myEphemerisPeggedChemRendezvous.getName() << std::endl;
    outputcheck << std::endl;
    myEphemerisPeggedChemRendezvous.output(outputcheck, options.launch_window_open_date, eventcount);
    outputcheck << std::endl;
*/

    outputcheck.close();

}//end main