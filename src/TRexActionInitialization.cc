/*
 * ActionInitialization class used as part of migration to Geant 4.10.03
 * Used to initialize all optional user action classes. Currently set up
 * only for sequential mode.
 * 
 * dhymers 2017/06/12
 */

#include "TRexActionInitialization.hh"

#include "MiniBallHistoManager.hh"
#include "MiniBallHistoGenerator.hh"
#include "MiniBallRootGenerator.hh"

#include "TRexPrimaryGeneratorAction.hh"
#include "TRexRunAction.hh"
#include "TRexEventAction.hh"
#include "TRexStackingAction.hh"
#include "HadronicStackingAction.hh"
#include "MySteppingAction.hh"

TRexActionInitialization::TRexActionInitialization(TRexDetectorConstruction* detConst, TRexData* data, MiniBallHistoManager* histoMan) : G4VUserActionInitialization(), fDetConst(detConst), fData(data), fHistoMan(histoMan)
{

}

TRexActionInitialization::~TRexActionInitialization(){

}

// worker threads and single threaded
void TRexActionInitialization::Build() const {
    TRexPrimaryGeneratorAction* gen_action = new TRexPrimaryGeneratorAction();
    SetUserAction(gen_action);

    TRexEventAction* event_action = new TRexEventAction(fDetConst, fHistoMan);
    SetUserAction(event_action);

    TRexRunAction* run_action = new TRexRunAction(fData, gen_action, event_action);
    SetUserAction(run_action);
    
    //TRexStackingAction* stackAction = new TRexStackingAction();
    //SetUserAction(stackAction);

    HadronicStackingAction* myAction = new HadronicStackingAction();
    SetUserAction(myAction);
    
    MySteppingAction* mystep_action = new MySteppingAction();
    SetUserAction(mystep_action);
}

//master thread in multithreaded run
void TRexActionInitialization::BuildForMaster() const{
    //TRexPrimaryGeneratorAction* gen_action = new TRexPrimaryGeneratorAction;

    //TRexEventAction* event_action = new TRexEventAction(fDetConst, fHistoMan);

    //TRexRunAction* run_action = new TRexRunAction(fData, gen_action, event_action);
    //SetUserAction(run_action);
}
