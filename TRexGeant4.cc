/*
 * TRexGeant4.cc
 *
 *  Created on: Jun 15, 2014
 *      Author: sklupp
 * 
 * modified to use custom QGSP_BIC physics list, which includes
 * nuclear recoils on GenericIon class.
 * dhymers 2017/06/12 and than by Leyla 2018/06/12
 */

//#ifdef G4MULTITHREADED
//#include "G4MTRunManager.hh"
//#else
#include "G4RunManager.hh"
//#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#endif


#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#include "G4OpenGLStoredX.hh"
#include "G4OpenGLStoredQt.hh"
#endif

#include "TRexDetectorConstruction.hh"
//#include "TRexPhysicsList.hh"
#include "PhysicsList.hh"

#include "my_QGSP_BIC.hh"

#include "TRexPrimaryGeneratorAction.hh"
#include "TRexRunAction.hh"
#include "TRexEventAction.hh"
#include "TRexData.hh"
#include "TRexSettings.hh"

#include "MiniBallHistoManager.hh"
#include "MiniBallHistoGenerator.hh"
#include "MiniBallRootGenerator.hh"

#include "TRexActionInitialization.hh"

#include "Randomize.hh"

#include "G4RadioactiveDecayPhysics.hh"

static __thread int foo;

int main(int argc,char** argv) {
  // class to organize command line
  TRexData data(argc, argv);

  // class which contains the settings stored in the settings file
  TRexSettings::Get()->ReadSettingsFile(data.GetSettingsFileName());

  // choose the Random engine
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
  G4long seed;
  seed = time(0);
  CLHEP::HepRandom::setTheSeed(seed);

  // Construct the default run manager  //
  //#ifdef G4MULTITHREADED
  //  G4cout << "RUNNING MULTITHREADED" << G4endl;
  //  G4int nThreads = 1;
  //  G4MTRunManager * runManager = new G4MTRunManager;
  //	runManager->SetNumberOfThreads(nThreads);
  //#else
  //Miniball structure restricts us to single threaded only
  G4RunManager * runManager = new G4RunManager;
  //#endif

  // initalize Miniball
  MiniBallRootGenerator* miniballHistoGen = new MiniBallRootGenerator(data.GetOutputFileName(), 5);
  MiniBallHistoManager* miniballHistoManager = new MiniBallHistoManager(miniballHistoGen);
  miniballHistoManager->SetMBTupleType(EVENT);
  miniballHistoManager->SetEnergyUnit(CLHEP::keV);

  // set mandatory initialization classes  //
  G4VUserDetectorConstruction* detector = new TRexDetectorConstruction(miniballHistoManager);
  runManager->SetUserInitialization(detector);

  //G4VUserPhysicsList* physics = new TRexPhysicsList;
  //using custom physics lists
  G4VUserPhysicsList* physics = new myQGSP_BIC();
  runManager->SetUserInitialization(physics);

  TRexActionInitialization* action = new TRexActionInitialization((TRexDetectorConstruction*)detector, &data, miniballHistoManager);
  runManager->SetUserInitialization(action);

  runManager->Initialize();

  // visualization
#ifdef G4VIS_USE
  G4VisManager* visManager = new G4VisExecutive;
  visManager->RegisterGraphicsSystem(new G4OpenGLStoredX);
  visManager->RegisterGraphicsSystem(new G4OpenGLStoredQt);
  visManager->Initialize();
#endif

  // define UI terminal
  G4UIsession* session = new G4UIterminal(new G4UItcsh);

  // run command
  data.run(runManager);

  // Job termination
  //
  // Free the store: user actions, physics_list and detector_description are
  //                 owned and deleted by the run manager, so they should not
  //                 be deleted in the main() program !
  //

#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;

  delete session;

  return 0;
}





