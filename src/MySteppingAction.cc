
#include "MySteppingAction.hh"
#include "TRexRunAction.hh"
#include "TRexDetectorConstruction.hh"
#include "TRexPrimaryGeneratorAction.hh"

#include "G4SteppingManager.hh"
#include "G4VTouchable.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// RunAction, DetectorConstruction, PrimaryGeneratorAction objets in argument
//MySteppingAction::MySteppingAction(TRexRunAction* run,TRexDetectorConstruction* det,TRexPrimaryGeneratorAction* pri) ***
//:Run(run),Detector(det),Primary(pri)***

MySteppingAction::MySteppingAction()
{ 
	
//myFile=fopen("beamIn_132Sn_792MeV_gas_traprob1_1000mb.dat","a");
myFile=fopen("beamIn_132Sn_792MeV_1000mb_target20cm_foilBe8u_ver2.txt","a");	
	
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

MySteppingAction::~MySteppingAction()
{
	
fclose (myFile);	
	
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//void MySteppingAction::UserSteppingAction(const G4Step* s)
void MySteppingAction::UserSteppingAction(const G4Step* s)
  
{ 

// Dose incrementation

 /*if (s->GetPreStepPoint()->GetPhysicalVolume()->GetName()  == "Target")
	{ 
   	 G4double dose = (e_SI*(s->GetTotalEnergyDeposit()/eV))/(Run->GetMassTarget());
   	 Run->AddDoseTarget(dose);
	}*/

// Beam spread after target

 if ( (s->GetTrack()->GetParentID() == 0)
   //(s->GetTrack()->GetDynamicParticle()->GetDefinition() ->GetParticleName() == "proton")
   //&& (s->GetPreStepPoint()->GetPhysicalVolume()->GetName()  == "Target")
   //&& (s->GetPostStepPoint()->GetPhysicalVolume()->GetName() == "World") 
    )
	{ 
 		
	fprintf 
	 /*( myFile,"%e%s%e%s%e%s%e%s",
	  (s->GetTrack()->GetPosition().x()) / CLHEP::micrometer,
	  "\t",
	  (s->GetTrack()->GetPosition().y()) / CLHEP::micrometer,
	  "\t",
	  (s->GetTrack()->GetPosition().z()) / CLHEP::micrometer,
	  "\t",
	  (s->GetTrack()->GetKineticEnergy()) / CLHEP::MeV,	  
	  "\n"
	 );*/
	 
	 ( myFile,"%e%s%e%s%e%s%e%s",
	  (s->GetPostStepPoint()->GetPosition().x()) / CLHEP::micrometer,
	  "\t",
	  (s->GetPostStepPoint()->GetPosition().y()) / CLHEP::micrometer,
	  "\t",
	  (s->GetPostStepPoint()->GetPosition().z()) / CLHEP::micrometer,
	  "\t",
	  (s->GetPostStepPoint()->GetKineticEnergy()) / CLHEP::MeV,	  
	  "\n"
	 );
	
	}
	
	
}
