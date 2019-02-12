#include "SteppingAction.hh"
#include "TRexEventAction.hh"
#include "TRexDetectorConstruction.hh"

#include "Analysis.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4RunManager.hh"

SteppingAction::SteppingAction( const TRexDetectorConstruction* detectorConstruction,TRexEventAction* eventAction)
	: G4UserSteppingAction(),
      fDetConstruction(detectorConstruction),
	  fEventAction(eventAction){
}

SteppingAction::~SteppingAction(){
}

void SteppingAction::UserSteppingAction(const G4Step* step){

	if (step->GetTrack()->GetParticleDefinition()->GetParticleName() == "gamma"){//Check if the track belongs to a gamma

		if (step->GetPostStepPoint()->GetStepStatus() == fWorldBoundary){//Check for the case where gamma leaves world, failure to filter leads to seg fault
				//std::cout << "world boundary" << std::endl;
		}
		else if (step->IsLastStepInVolume() && step->GetPostStepPoint()->GetMaterial()->GetName() == "HPGermanium"){//If we enter into a germanium detector take the full energy
			G4double ekin = step->GetPreStepPoint()->GetKineticEnergy();
			//std::cout<<ekin<<std::endl;
			fEventAction->setGamTotal(ekin);
			//std::cout<<fEventAction->getGamTotal()<<std::endl;
	       	//G4String secondVolume = step->GetPostStepPoint()->GetMaterial()->GetName();
			//std::cout<<secondVolume<<std::endl;
			step->GetTrack()->SetTrackStatus(fStopAndKill);
			
		}
	}

}
