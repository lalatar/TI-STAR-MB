/*
 * TRexVacuumChamberGasCylinder.cc
 *
 *  Created on: Jan 31, 2017
 *      Author: vbildste
 */

#include "TRexVacuumChamberGasCylinder.hh"
#include "TRexSettings.hh"

TRexVacuumChamberGasCylinder::TRexVacuumChamberGasCylinder() {
}

TRexVacuumChamberGasCylinder::~TRexVacuumChamberGasCylinder() {
	// TODO Auto-generated destructor stub
}

G4LogicalVolume* TRexVacuumChamberGasCylinder::ConstructChamberGas(G4LogicalVolume* experimentalHall_log) {
	G4Material* chamberGasMaterial = TRexMaterials::Get()->GetMaterial("chamberGas");

	G4Tubs *chamberGas_tubs;
	G4Cons *chamberGas_cons1;
	G4Cons *chamberGas_cons2;

	if(TRexSettings::Get()->VisualizationCut()) {
		chamberGas_tubs = new G4Tubs("ChamberTubs",0.*CLHEP::mm,70.*CLHEP::mm,112.*CLHEP::mm,270.*CLHEP::degree,180.*CLHEP::degree);//original
	    chamberGas_cons1 = new G4Cons("ChamberCons1",0.*CLHEP::mm,59.*CLHEP::mm,0.*CLHEP::mm,70.*CLHEP::mm,10.*CLHEP::mm,270.*CLHEP::degree,180.*CLHEP::degree);//original
		chamberGas_cons2 = new G4Cons("ChamberCons2",0.*CLHEP::mm,12.5*CLHEP::mm,0.*CLHEP::mm,59.*CLHEP::mm,9.5*CLHEP::mm,270.*CLHEP::degree,180.*CLHEP::degree);//original 
		//chamberGas_tubs = new G4Tubs("ChamberTubs",0.*CLHEP::mm,0.1*CLHEP::mm,0.1*CLHEP::mm,0.*CLHEP::degree,0.1*CLHEP::degree);//changed by Leila no chamber
		//chamberGas_cons1 = new G4Cons("ChamberCons1",0.*CLHEP::mm,0.1*CLHEP::mm,0.*CLHEP::mm,0.1*CLHEP::mm,0.1*CLHEP::mm,0.*CLHEP::degree,0.1*CLHEP::degree);
		//chamberGas_cons2 = new G4Cons("ChamberCons2",0.*CLHEP::mm,0.1*CLHEP::mm,0.*CLHEP::mm,0.1*CLHEP::mm,0.1*CLHEP::mm,0.*CLHEP::degree,0.1*CLHEP::degree);
	} else {
		//chamberGas_tubs = new G4Tubs("ChamberTubs",0.*CLHEP::mm,70.*CLHEP::mm,112.*CLHEP::mm,0.*CLHEP::degree,360.*CLHEP::degree);// original
		chamberGas_tubs = new G4Tubs("ChamberTubs",0.*CLHEP::mm,80.*CLHEP::mm,112.*CLHEP::mm,0.*CLHEP::degree,360.*CLHEP::degree);// changed by  Leila
		chamberGas_cons1 = new G4Cons("ChamberCons1",0.*CLHEP::mm,59.*CLHEP::mm,0.*CLHEP::mm,70.*CLHEP::mm,10.*CLHEP::mm,0.*CLHEP::degree,360.*CLHEP::degree);
		chamberGas_cons2 = new G4Cons("ChamberCons2",0.*CLHEP::mm,12.5*CLHEP::mm,0.*CLHEP::mm,59.*CLHEP::mm,9.5*CLHEP::mm,0.*CLHEP::degree,360.*CLHEP::degree);
	}

	G4UnionSolid* chamberGas_cons = new G4UnionSolid("chambergas_cons", chamberGas_cons1, chamberGas_cons2, new G4RotationMatrix(),G4ThreeVector(0,0,-19.5*CLHEP::mm));

	G4UnionSolid* chamberGas_solid = new G4UnionSolid("chambergas_solid", chamberGas_tubs, chamberGas_cons, new G4RotationMatrix(),G4ThreeVector(0,0,-122.05*CLHEP::mm));

	// logical volume
	//G4LogicalVolume* chamberGas_log = new G4LogicalVolume(chamberGas_solid, chamberGasMaterial, "chambergas_log"); //original
	G4LogicalVolume* chamberGas_log = new G4LogicalVolume(chamberGas_tubs, chamberGasMaterial, "chambergas_log");// 

	// physical volume
	new G4PVPlacement(0, G4ThreeVector(0,0,0), chamberGas_log, "chamber_gas", experimentalHall_log, false, 0);

	return chamberGas_log;
}
