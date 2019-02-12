/*
 * TRexVacuumChamberCylinder.cc
 *
 *  Created on: Jun 16, 2014
 *      Author: sklupp
 */

#include "TRexVacuumChamberBox.hh"
#include "TRexSettings.hh"

TRexVacuumChamberBox::TRexVacuumChamberBox() {
}

TRexVacuumChamberBox::~TRexVacuumChamberBox() {
	// TODO Auto-generated destructor stub
}

void TRexVacuumChamberBox::ConstructChamber(G4LogicalVolume* experimentalHall_log) {
	G4Material* chamberMaterial = TRexMaterials::Get()->GetMaterial("aluminium");
	
	G4double boxWidth = 150 * CLHEP::mm; // in x direction
	G4double boxHeight = 150. * CLHEP::mm;// in y direction
	G4double boxLength = 200. * CLHEP::mm;// in z direction
	G4double boxThickness = 2. * CLHEP::mm;

	G4Box* chamber_box1 = new G4Box("chamber_box1", boxWidth / 2., boxHeight / 2., boxLength / 2.);
	G4Box* chamber_box2 = new G4Box("chamber_box2", (boxWidth-boxThickness) / 2., (boxHeight-boxThickness) / 2., (boxLength-boxThickness) / 2.);
	
	G4SubtractionSolid* chamber_box = new G4SubtractionSolid("chamber_box", chamber_box1, chamber_box2, 0, G4ThreeVector(0.,0.,0.));

	// logical volume
	G4LogicalVolume* chamber_log = new G4LogicalVolume(chamber_box, chamberMaterial, "chamber_log");

	// physical volume
	new G4PVPlacement(0, G4ThreeVector(0,0,0), chamber_log, "chamber", experimentalHall_log, false, 0);

	// grey color
	//chamber_log->SetVisAttributes(new G4VisAttributes(true, G4Colour(0.7,0.7,0.7)));
	chamber_log->SetVisAttributes(TRexColour::Get()->gray);
}
