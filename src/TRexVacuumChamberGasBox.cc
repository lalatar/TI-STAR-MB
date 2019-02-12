/*
 * TRexVacuumChamberGasCylinder.cc
 *
 *  Created on: Jan 31, 2017
 *      Author: vbildste
 */

#include "TRexVacuumChamberGasBox.hh"
#include "TRexSettings.hh"

TRexVacuumChamberGasBox::TRexVacuumChamberGasBox() {
}

TRexVacuumChamberGasBox::~TRexVacuumChamberGasBox() {
	// TODO Auto-generated destructor stub
}

G4LogicalVolume* TRexVacuumChamberGasBox::ConstructChamberGas(G4LogicalVolume* experimentalHall_log) {
	G4Material* chamberGasMaterial = TRexMaterials::Get()->GetMaterial("chamberGas");
	
	G4double boxWidth = 150 * CLHEP::mm; // in x direction
	G4double boxHeight = 150. * CLHEP::mm; // 1.6
	G4double boxLength = 200. * CLHEP::mm;// in z direction
	G4double boxThickness = 2. * CLHEP::mm;

	G4Box* chamberGas_box = new G4Box("chamberGas_box", (boxWidth-boxThickness) / 2., (boxHeight-boxThickness) / 2., (boxLength-boxThickness) / 2.);

	// logical volume
	G4LogicalVolume* chamberGas_log = new G4LogicalVolume(chamberGas_box, chamberGasMaterial, "chambergas_log");// 

	// physical volume
	new G4PVPlacement(0, G4ThreeVector(0,0,0), chamberGas_log, "chamber_gas", experimentalHall_log, false, 0);

	return chamberGas_log;
}
