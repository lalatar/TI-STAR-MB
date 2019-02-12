/*
 * TRexVacuumChamberCylinder.cc
 *
 *  Created on: Jun 16, 2014
 *      Author: sklupp
 */

#include "TRexVacuumChamberCylinder.hh"
#include "TRexSettings.hh"

TRexVacuumChamberCylinder::TRexVacuumChamberCylinder() {
}

TRexVacuumChamberCylinder::~TRexVacuumChamberCylinder() {
	// TODO Auto-generated destructor stub
}

void TRexVacuumChamberCylinder::ConstructChamber(G4LogicalVolume* experimentalHall_log) {
	G4Material* chamberMaterial = TRexMaterials::Get()->GetMaterial("aluminium");

	G4Tubs *chamber_tubs;
	//G4Tubs *chamber_tubs1;
	//G4Tubs *chamber_tubs2;
	//G4Tubs *chamber_tubs3;
	G4Cons *chamber_cons1;
	G4Cons *chamber_cons2;

	if(TRexSettings::Get()->VisualizationCut()) {
		//chamber_tubs = new G4Tubs("ChamberTubs",70.*CLHEP::mm,72.*CLHEP::mm,62.*CLHEP::mm,270.*CLHEP::degree,360.*CLHEP::degree);
		chamber_tubs = new G4Tubs("ChamberTubs",70.*CLHEP::mm,72.*CLHEP::mm,112.*CLHEP::mm,270.*CLHEP::degree,180.*CLHEP::degree);
		//chamber_tubs1 = new G4Tubs("ChamberTubs1",70.*CLHEP::mm,72.*CLHEP::mm,34.*CLHEP::mm,270.*CLHEP::degree,180.*CLHEP::degree);
		//chamber_tubs2 = new G4Tubs("ChamberTubs2",70.*CLHEP::mm,90.*CLHEP::mm,62.5*CLHEP::mm,270.*CLHEP::degree,180.*CLHEP::degree);
		//chamber_tubs3 = new G4Tubs("ChamberTubs3",51.*CLHEP::mm,54.*CLHEP::mm,41.5*CLHEP::mm,270.*CLHEP::degree,180.*CLHEP::degree);
		chamber_cons1 = new G4Cons("ChamberCons1",59.*CLHEP::mm,61.*CLHEP::mm,70.*CLHEP::mm,72.*CLHEP::mm,10.*CLHEP::mm,270.*CLHEP::degree,180.*CLHEP::degree);
		chamber_cons2 = new G4Cons("ChamberCons2",12.5*CLHEP::mm,14.5*CLHEP::mm,59.*CLHEP::mm,61.*CLHEP::mm,9.5*CLHEP::mm,270.*CLHEP::degree,180.*CLHEP::degree);
	} else {
		//chamber_tubs = new G4Tubs("ChamberTubs",70.*CLHEP::mm,72.*CLHEP::mm,62.*CLHEP::mm,0.*CLHEP::degree,360.*CLHEP::degree);
		//chamber_tubs = new G4Tubs("ChamberTubs",70.*CLHEP::mm,72.*CLHEP::mm,112.*CLHEP::mm,0.*CLHEP::degree,360.*CLHEP::degree);// original
		chamber_tubs = new G4Tubs("ChamberTubs",80.*CLHEP::mm,82.*CLHEP::mm,112.*CLHEP::mm,0.*CLHEP::degree,360.*CLHEP::degree);// changed by Leila
		//chamber_tubs1 = new G4Tubs("ChamberTubs1",70.*CLHEP::mm,72.*CLHEP::mm,34.*CLHEP::mm,0.*CLHEP::degree,360.*CLHEP::degree);
		//chamber_tubs2 = new G4Tubs("ChamberTubs2",70.*CLHEP::mm,90.*CLHEP::mm,62.5*CLHEP::mm,0.*CLHEP::degree,360.*CLHEP::degree);
		//chamber_tubs3 = new G4Tubs("ChamberTubs3",51.*CLHEP::mm,54.*CLHEP::mm,41.5*CLHEP::mm,0.*CLHEP::degree,360.*CLHEP::degree); 
		chamber_cons1 = new G4Cons("ChamberCons1",59.*CLHEP::mm,61.*CLHEP::mm,70.*CLHEP::mm,72.*CLHEP::mm,10.*CLHEP::mm,0.*CLHEP::degree,360.*CLHEP::degree);
		chamber_cons2 = new G4Cons("ChamberCons2",12.5*CLHEP::mm,14.5*CLHEP::mm,59.*CLHEP::mm,61.*CLHEP::mm,9.5*CLHEP::mm,0.*CLHEP::degree,360.*CLHEP::degree);
	}

	G4UnionSolid* chamber_cons = new G4UnionSolid("chamber_cons", chamber_cons1, chamber_cons2, new G4RotationMatrix(),G4ThreeVector(0,0,-19.5*CLHEP::mm));

	//G4UnionSolid* chamber_solid = new G4UnionSolid("chamber_solid", chamber_tubs, chamber_cons, new G4RotationMatrix(),G4ThreeVector(0,0,-72.05*CLHEP::mm));
	G4UnionSolid* chamber_solid = new G4UnionSolid("chamber_solid", chamber_tubs, chamber_cons, new G4RotationMatrix(),G4ThreeVector(0,0,-122.05*CLHEP::mm));

	// logical volume
	//G4LogicalVolume* chamber_log = new G4LogicalVolume(chamber_solid, chamberMaterial, "chamber_log");//original
	G4LogicalVolume* chamber_log = new G4LogicalVolume(chamber_tubs, chamberMaterial, "chamber_log");// changed bei Leila

	// physical volume
	new G4PVPlacement(0, G4ThreeVector(0,0,0), chamber_log, "chamber", experimentalHall_log, false, 0);

	// grey color
	//chamber_log->SetVisAttributes(new G4VisAttributes(true, G4Colour(0.7,0.7,0.7)));
	chamber_log->SetVisAttributes(TRexColour::Get()->gray);
	//chamber_log->SetVisAttributes(TRexColour::Get()->yellow);
}
