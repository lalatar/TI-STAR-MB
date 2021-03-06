/*
 * TRexDetectorConstruction.cc
 *
 *  Created on: Jun 15, 2014
 *      Author: sklupp
 */

#include "TRexDetectorConstruction.hh"
#include "TRexMaterials.hh"
#include "TRexSettings.hh"
#include "TRexColour.hh"
#include "TRexVacuumChamberCylinder.hh"
#include "TRexVacuumChamberSphere.hh"
#include "TRexVacuumChamberBox.hh"
#include "TRexVacuumChamberGasCylinder.hh"
#include "TRexVacuumChamberGasSphere.hh"
#include "TRexVacuumChamberGasBox.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4UnitsTable.hh"
#include "globals.hh"

#include "TMath.h"

#include <cstdlib>
#include <iomanip>

TRexDetectorConstruction::TRexDetectorConstruction() {
	fMiniballHistoManager = NULL;
}

TRexDetectorConstruction::TRexDetectorConstruction(MiniBallHistoManager* miniballHistoManager) :
	fMiniballHistoManager(miniballHistoManager) {
	fSDMan = G4SDManager::GetSDMpointer();
	DefineDetectors();
}

TRexDetectorConstruction::~TRexDetectorConstruction() {
	// TODO Auto-generated destructor stub
}

G4VPhysicalVolume* TRexDetectorConstruction::Construct() {
	ConstructExperimentalHall();

	// Miniball for gammas
	if(TRexSettings::Get()->SimulateGammas()) {
		DefineMiniball();
	}

	// vacuum chamber around silicon detectors and vacuum chamber gas into which all detectors and the target aree placed
	ConstructVacuumChamber();

	// Silicon particle detectors
	for(unsigned int i = 0; i < fDetectors.size(); i++) {
		std::cout << "Constructing detector " << fDetectors[i]->GetName() << " ....\n" << std::endl;
		fDetectors[i]->Construct(fVacuumChamberGas_log, fSDMan);
	}

	ConstructTarget();

	return fExperimentalHall_phys; 
}

void TRexDetectorConstruction::ConstructTarget() {
		
	// experimental hall material
	G4Material* targetMaterial = TRexMaterials::Get()->GetMaterial("target");

	G4double innerRadius = 0.0*CLHEP::mm;
	G4double outerRadius = TRexSettings::Get()->GetTargetDiameter()/2.;
	G4double startAngle = 0.*CLHEP::deg;
	G4double spanningAngle = 360.*CLHEP::deg;
	G4double halfThickness = TRexSettings::Get()->GetTargetPhysicalLength()/2;
	
	// ################### Target Entrance Foil added by Leila begin ###################
	
	G4Material* targetEntranceFoilMaterial = TRexMaterials::Get()->GetMaterial("beryllium");
	
	G4double targetEntranceFoilHalfThickness = 4.0*CLHEP::um;// total thickness is 8 um
	 
	G4Tubs* targetEntranceFoil_tube = new G4Tubs("targetEntranceFoil",
			innerRadius, outerRadius,
			targetEntranceFoilHalfThickness,
			startAngle,	spanningAngle);
			
	G4LogicalVolume* targetEntranceFoil_log = new G4LogicalVolume(targetEntranceFoil_tube,targetEntranceFoilMaterial,"targetEntranceFoil_log", 0, 0, 0);
	
	//G4UnionSolid* beFoil_target = new G4UnionSolid("beFoil_target", targetEntranceFoil_tube,target_tube);
	
	new G4PVPlacement(0, G4ThreeVector(0,0,-(TRexSettings::Get()->GetTargetPhysicalLength()/2 + targetEntranceFoilHalfThickness)*CLHEP::mm), targetEntranceFoil_log, "target_beFoil",fVacuumChamberGas_log,false,0);
	new G4PVPlacement(0, G4ThreeVector(0,0,(TRexSettings::Get()->GetTargetPhysicalLength()/2 + targetEntranceFoilHalfThickness)*CLHEP::mm), targetEntranceFoil_log, "target_beFoil",fVacuumChamberGas_log,false,0);
	
	if(TRexSettings::Get()->Colours()) {targetEntranceFoil_log->SetVisAttributes(TRexColour::Get()->cyan);}
				
	// ################### Target Entrance Foil added by Leila  end  ###################
	
	// ################### Target Mantel (Foil) added by Leila start ###################
	
	G4Material* targetMantelFoilMaterial = TRexMaterials::Get()->GetMaterial("mylar");
	
	G4double targetMantelFoilHalfThickness = TRexSettings::Get()->GetTargetPhysicalLength()/2;
	G4double targetMantelinnerRadius = TRexSettings::Get()->GetTargetDiameter()/2.;
	G4double targetMantelouterRadius = TRexSettings::Get()->GetTargetDiameter()/2. + TRexSettings::Get()->GetFBarrelDeltaESingleFoilThickness()/2.;
	G4double targetMantelstartAngle = 0.*CLHEP::deg;
	G4double targetMantelspanningAngle = 360.*CLHEP::deg;			
	
	G4Tubs* targetMantelFoil_tube = new G4Tubs("targetMantelFoil",
			targetMantelinnerRadius, targetMantelouterRadius,
			targetMantelFoilHalfThickness,
			targetMantelstartAngle,	targetMantelspanningAngle);
			
	G4LogicalVolume* targetMantelFoil_log = new G4LogicalVolume(targetMantelFoil_tube,targetMantelFoilMaterial,"targetMantelFoil_log", 0, 0, 0);
	
	//new G4PVPlacement(0, G4ThreeVector(0,0,-80.0*CLHEP::mm), targetMantelFoil_log, "target_MantelFoil",fVacuumChamberGas_log,false,0);
	new G4PVPlacement(0, G4ThreeVector(), targetMantelFoil_log, "target_MantelFoil",fVacuumChamberGas_log,false,0);
	
	if(TRexSettings::Get()->Colours()) {targetMantelFoil_log->SetVisAttributes(TRexColour::Get()->magenta);}
	
	// ################### Target Mantel (Foil) added by Leila  end  ###################
	
	// ################### Al shielding added by Leila  start  ###################
	
	/***G4Material* shieldMaterial = TRexMaterials::Get()->GetMaterial("aluminium");
	
	G4double shieldHalfLengthX = 7.5* CLHEP::mm;
	G4double shieldHalfLengthY = 1.0* CLHEP::mm;
	G4double shieldHalfLengthZ = 4.0* CLHEP::mm;
	
	G4Box* shiled = new G4Box("Shielding", shieldHalfLengthX, shieldHalfLengthY, shieldHalfLengthZ);
	
	G4LogicalVolume* shield_log = new G4LogicalVolume(shiled,shieldMaterial,"shield_log", 0, 0, 0);
	
	new G4PVPlacement(0, G4ThreeVector(0.0*CLHEP::mm,7.5*CLHEP::mm,-TRexSettings::Get()->GetTargetPhysicalLength()/2*CLHEP::mm + shieldHalfLengthZ), shield_log, "Shielding",fVacuumChamberGas_log,false,0);
	new G4PVPlacement(0, G4ThreeVector(0.0*CLHEP::mm,7.5*CLHEP::mm,TRexSettings::Get()->GetTargetPhysicalLength()/2*CLHEP::mm - shieldHalfLengthZ), shield_log, "Shielding",fVacuumChamberGas_log,false,0);
	new G4PVPlacement(0, G4ThreeVector(0.0*CLHEP::mm,-7.5*CLHEP::mm,-TRexSettings::Get()->GetTargetPhysicalLength()/2*CLHEP::mm + shieldHalfLengthZ), shield_log, "Shielding",fVacuumChamberGas_log,false,0);
	new G4PVPlacement(0, G4ThreeVector(0.0*CLHEP::mm,-7.5*CLHEP::mm,TRexSettings::Get()->GetTargetPhysicalLength()/2*CLHEP::mm - shieldHalfLengthZ), shield_log, "Shielding",fVacuumChamberGas_log,false,0);
	
	if(TRexSettings::Get()->Colours()) {shield_log->SetVisAttributes(TRexColour::Get()->green);}***/
		
	// ################### Al shielding added by Leila  end  ###################

	G4Tubs* target_tube = new G4Tubs("detector",
			innerRadius, outerRadius,
			halfThickness,
			startAngle,	spanningAngle);

	G4LogicalVolume* target_log = new G4LogicalVolume(target_tube,
			targetMaterial, "target_log", 0, 0, 0);

	//G4VPhysicalVolume* target_phys =
	new G4PVPlacement(0, G4ThreeVector(), target_log, "target",fVacuumChamberGas_log,false,0);

	if(TRexSettings::Get()->Colours()) {
		target_log->SetVisAttributes(TRexColour::Get()->lightblue);
	}

	if(targetMaterial->GetState() == kStateGas) {
		TRexSettings::Get()->SetTargetThickness(target_log->GetMass()/(TMath::Pi()*TMath::Power(outerRadius,2)));
		TRexSettings::Get()->SetTargetMaterialDensity(target_log->GetMass()/(TMath::Pi()*TMath::Power(outerRadius,2)*2.*halfThickness));
		std::cout<<"Built gas target with pressure "<<targetMaterial->GetPressure()/CLHEP::bar*1000.<<" mbar, and "<<2.*halfThickness/CLHEP::cm<<" cm length => area density = "<<target_log->GetMass()/(TMath::Pi()*TMath::Power(outerRadius,2))/(CLHEP::mg/CLHEP::cm2)<<" mg/cm2"<<" target mass from file: "<<target_log->GetMass()<<std::endl;
		
	} else {
		std::cout<<"Built solid target with thickness "<<2.*halfThickness/um<<" um => area density = "<<target_log->GetMass()/(TMath::Pi()*TMath::Power(outerRadius,2))/(CLHEP::mg/CLHEP::cm2)<<" mg/cm2 and density "<<TRexSettings::Get()->GetTargetMaterialDensity()/(CLHEP::g/CLHEP::cm3)<<" g/cm3"<<std::endl;
	}		
	
}

void TRexDetectorConstruction::ConstructVacuumChamber() {
	if(TRexSettings::Get()->IncludeVacuumChamber()) {
		if(TRexSettings::Get()->GetVacuumChamberType() == "cylinder") {
			std::cout << "Constructing vacuum chamber with type " << TRexSettings::Get()->GetVacuumChamberType() << " ....\n" << std::endl;
			fVacuumChamber = new TRexVacuumChamberCylinder();
		} else if(TRexSettings::Get()->GetVacuumChamberType() == "sphere") {
			std::cout << "Constructing vacuum chamber with type " << TRexSettings::Get()->GetVacuumChamberType() << " ....\n" << std::endl;
			fVacuumChamber = new TRexVacuumChamberSphere();
		} else {
			std::cerr << "Vacuum chamber with type " << TRexSettings::Get()->GetVacuumChamberType() << " not implemented yet !" << std::endl;
			exit(1);
		}
		fVacuumChamber->ConstructChamber(fExperimentalHall_log);
	}
	if(TRexSettings::Get()->GetVacuumChamberType() == "cylinder") {
		std::cout << "Constructing vacuum chamber gas with type " << TRexSettings::Get()->GetVacuumChamberType() << " ....\n" << std::endl;
		fVacuumChamberGas = new TRexVacuumChamberGasCylinder();
	} else if(TRexSettings::Get()->GetVacuumChamberType() == "sphere") {
		std::cout << "Constructing vacuum chamber gas with type " << TRexSettings::Get()->GetVacuumChamberType() << " ....\n" << std::endl;
		fVacuumChamberGas = new TRexVacuumChamberGasSphere();
	} else if(TRexSettings::Get()->GetVacuumChamberType() == "box") {
		std::cout << "Constructing vacuum chamber gas with type " << TRexSettings::Get()->GetVacuumChamberType() << " ....\n" << std::endl;
		fVacuumChamberGas = new TRexVacuumChamberGasBox();
	} else {
		std::cerr << "Vacuum chamber with type " << TRexSettings::Get()->GetVacuumChamberType() << " not implemented yet !" << std::endl;
		exit(1);
	}
	fVacuumChamberGas_log = fVacuumChamberGas->ConstructChamberGas(fExperimentalHall_log);
}

void TRexDetectorConstruction::ConstructExperimentalHall() {
	// experimental hall material
	//G4Material* experimentalHallMaterial = TRexMaterials::Get()->vacuum;
	G4Material* experimentalHallMaterial = MiniBallMaterial::Get()->GetMaterialVacuum();


	/////////////////////////////////////////////////////////////////////////////
	// experimental hall
	/////////////////////////////////////////////////////////////////////////////
	G4double expHall_x = 1.0*CLHEP::m;
	G4double expHall_y = 1.0*CLHEP::m;
	G4double expHall_z = 1.0*CLHEP::m;

	fExperimentalHall_box = new G4Box("expHall_box",
			expHall_x,expHall_y,expHall_z);

	fExperimentalHall_log = new G4LogicalVolume(fExperimentalHall_box,
			experimentalHallMaterial,
			"expHall_log", 0, 0, 0);

	fExperimentalHall_phys = new G4PVPlacement(0, G4ThreeVector(),
			fExperimentalHall_log,
			"expHall",0,false,0);

	// set experimental hall invisible in the visualization
	fExperimentalHall_log->SetVisAttributes(G4VisAttributes::Invisible);
}

void TRexDetectorConstruction::DefineDetectors() {
	// Forward CD deltaE
	if(TRexSettings::Get()->UseFCDdeltaE() && !TRexSettings::Get()->UseFCDErest()) {
		// loop over all FCDErest single detectors
		//for(int i = 0; i < TRexSettings::Get()->GetNbOfFCDdeltaESingles(); i++) {
		//	fDetectors.push_back(new TRexCDdeltaESingle(Form("FCDdeltaESingle%d", i), i);
		//}
		fDetectors.push_back(new TRexCDdeltaE("FCDdeltaE", "forward", TRexSettings::Get()->GetNbOfFCDdeltaESingles()));
	}
	// Forward CD Erest
	if(TRexSettings::Get()->UseFCDErest() && !TRexSettings::Get()->UseFCDdeltaE()) {
		// loop over all FCDErest single detectors
		//for(int i = 0; i < TRexSettings::Get()->GetNbOfFCDErestSingles(); i++) {
		//	fDetectors.push_back(new TRexCDErestSingle(Form("FCDErestSingle%d", i), i);
		//}
		fDetectors.push_back(new TRexCDErest("FCDErest", "forward", TRexSettings::Get()->GetNbOfFCDErestSingles()));
	}
	// Forward CD deltaE and forward CD Erest
	if(TRexSettings::Get()->UseFCDdeltaE() && TRexSettings::Get()->UseFCDErest()) {
		// loop over all FCDErest single detectors
		//for(int i = 0; i < TRexSettings::Get()->GetNbOfFCDdeltaESingles(); i++) {
		//	fDetectors.push_back(new TRexCDdeltaESingle(Form("FCDdeltaESingle%d", i), i);
		//}
		fDetectors.push_back(new TRexCD("FCD", "forward", TRexSettings::Get()->GetNbOfFCDdeltaESingles()));
	}

	// Backward CD deltaE
	if(TRexSettings::Get()->UseBCDdeltaE() && !TRexSettings::Get()->UseBCDErest()) {
		fDetectors.push_back(new TRexCDdeltaE("BCDdeltaE", "backward", TRexSettings::Get()->GetNbOfBCDdeltaESingles()));
	}
	// Backward CD Erest
	if(TRexSettings::Get()->UseBCDErest() && !TRexSettings::Get()->UseBCDdeltaE()) {
		fDetectors.push_back(new TRexCDErest("BCDErest", "backward", TRexSettings::Get()->GetNbOfBCDErestSingles()));
	}
	// Backward CD deltaE and backward CD Erest
	if(TRexSettings::Get()->UseBCDdeltaE() && TRexSettings::Get()->UseBCDErest()) {
		fDetectors.push_back(new TRexCD("BCD", "backward", TRexSettings::Get()->GetNbOfBCDdeltaESingles()));
	}

	// Forward Barrel deltaE
	if(TRexSettings::Get()->UseFBarrelDeltaE()) {
		fDetectors.push_back(new TRexBarrelDeltaE("FBarrelDeltaE", "forward", TRexSettings::Get()->GetNbOfFBarrelDeltaESingles()));
	}
	// Forward Barrel Erest
	if(TRexSettings::Get()->UseFBarrelErest()) {
		fDetectors.push_back(new TRexBarrelErest("FBarrelErest", "forward", TRexSettings::Get()->GetNbOfFBarrelErestSingles()));
	}
	// 2nd Forward Barrel deltaE
	if(TRexSettings::Get()->UseSecondFBarrelDeltaE()) {
		fDetectors.push_back(new TRexBarrelDeltaE("SecondFBarrelDeltaE", "forward", TRexSettings::Get()->GetNbOfSecondFBarrelDeltaESingles()));
	}

	// Middle Barrel deltaE
	if(TRexSettings::Get()->UseMBarrelDeltaE()) {
		fDetectors.push_back(new TRexBarrelDeltaE("MBarrelDeltaE", "middle", TRexSettings::Get()->GetNbOfMBarrelDeltaESingles()));
	}
	// Middle Barrel Erest
	if(TRexSettings::Get()->UseMBarrelErest()) {
		fDetectors.push_back(new TRexBarrelErest("MBarrelErest", "middle", TRexSettings::Get()->GetNbOfMBarrelErestSingles()));
	}

	// Backward Barrel deltaE
	if(TRexSettings::Get()->UseBBarrelDeltaE()) {
		fDetectors.push_back(new TRexBarrelDeltaE("BBarrelDeltaE", "backward", TRexSettings::Get()->GetNbOfBBarrelDeltaESingles()));
	}
	// Backward Barrel Erest
	if(TRexSettings::Get()->UseBBarrelErest()) {
		fDetectors.push_back(new TRexBarrelErest("BBarrelErest", "backward", TRexSettings::Get()->GetNbOfBBarrelErestSingles()));
	}
	// 2nd Forward Barrel deltaE
	if(TRexSettings::Get()->UseSecondBBarrelDeltaE()) {
		fDetectors.push_back(new TRexBarrelDeltaE("SecondBBarrelDeltaE", "backward", TRexSettings::Get()->GetNbOfSecondBBarrelDeltaESingles()));
	}

	// Forward Trapezoid deltaE
	if(TRexSettings::Get()->UseFTrapezoidDeltaE()) {
		fDetectors.push_back(new TRexTrapezoidDeltaE("FTrapezoidDeltaE", "forward", TRexSettings::Get()->GetNbOfFTrapezoidDeltaESingles()));
	}
	// Forward Trapezoid deltaE
	if(TRexSettings::Get()->UseFTrapezoidErest()) {
		fDetectors.push_back(new TRexTrapezoidErest("FTrapezoidErest", "forward", TRexSettings::Get()->GetNbOfFTrapezoidErestSingles()));
	}

	// Backward Trapezoid deltaE
	if(TRexSettings::Get()->UseBTrapezoidDeltaE()) {
		fDetectors.push_back(new TRexTrapezoidDeltaE("BTrapezoidDeltaE", "backward", TRexSettings::Get()->GetNbOfBTrapezoidDeltaESingles()));
	}
	// Backward Trapezoid deltaE
	if(TRexSettings::Get()->UseBTrapezoidErest()) {
		fDetectors.push_back(new TRexTrapezoidErest("BTrapezoidErest", "backward", TRexSettings::Get()->GetNbOfBTrapezoidErestSingles()));
	}
}

void TRexDetectorConstruction::DefineMiniball() {
	std::cout << "Constructing MINIBALL ....\n" << std::endl;

	G4bool MiniballCapsule = true;
	G4bool MiniballCryoWall = true;

	fMiniball = new MiniBallDetectorArray(fMiniballHistoManager, fExperimentalHall_phys, TRexSettings::Get()->GetNbOfMiniballClusters(),
			TRexSettings::Get()->GetMiniballClusterNames(), TRexSettings::Get()->GetMiniballClusterTypes(), TRexSettings::Get()->GetMiniballClusterRho(),
			TRexSettings::Get()->GetMiniballClusterTheta(), TRexSettings::Get()->GetMiniballClusterPhi(), TRexSettings::Get()->GetMiniballClusterSpin(),
			MiniballCapsule, MiniballCryoWall);

	fMiniball->SetAllDeadlayer(0.7*CLHEP::micrometer, DEADLAYER_OUTER_CONTACT);

/*		for(i = 0; i < 8; i++) {
		if(data->settings->MiniballRelativeEfficiency[i][0] > 0.) {
		fMiniball->SetRelativeEfficiencyDetA(i, data->settings->MiniballRelativeEfficiency[i][0]);
		}
		if(data->settings->MiniballRelativeEfficiency[i][1] > 0.) {
		fMiniball->SetRelativeEfficiencyDetB(i, data->settings->MiniballRelativeEfficiency[i][1]);
		}
		if(data->settings->MiniballRelativeEfficiency[i][2] > 0.) {
		fMiniball->SetRelativeEfficiencyDetC(i, data->settings->MiniballRelativeEfficiency[i][2]);
		}
		}*/

	fMiniball->Construct();
}
