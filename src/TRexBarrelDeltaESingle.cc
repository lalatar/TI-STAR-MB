/*
 * TRexBarrelDeltaESingle.cc
 *
 *  Created on: Jun 16, 2014
 *      Author: sklupp
 */

#include "TRexBarrelDeltaESingle.hh"
#include "TRexBarrelDeltaESingleSensitiveDetector.hh"
#include "TRexSettings.hh"

TRexBarrelDeltaESingle::TRexBarrelDeltaESingle() {
}

TRexBarrelDeltaESingle::TRexBarrelDeltaESingle(std::string name, std::string direction, int id) {
	fName = name;
	fId = id;
	fDirection = direction;

	fBaseName = fName;
	fBaseName.erase(fBaseName.end()-1, fBaseName.end());

	//if(fDirection == "forward") {
	if(fBaseName == "FBarrelDeltaESingle") {
		fDetectorLengthX = TRexSettings::Get()->GetFBarrelDeltaESingleLengthX();
		fDetectorLengthY = TRexSettings::Get()->GetFBarrelDeltaESingleLengthY();
		fThicknessDetector = TRexSettings::Get()->GetFBarrelDeltaESingleThickness()[fId];
		fStartAngleDetector = TRexSettings::Get()->GetFBarrelDeltaESingleStartAngle()[fId];

		fDeadLayer = TRexSettings::Get()->GetFBarrelDeltaESingleDeadLayer();

		fPos = G4ThreeVector(TRexSettings::Get()->GetFBarrelDeltaESingleDistanceToBeam()[fId] * cos(fStartAngleDetector / CLHEP::rad) + (13.0/3.) * CLHEP::mm, // ### target shifting -13.0/3 | -12.24 mm
									TRexSettings::Get()->GetFBarrelDeltaESingleDistanceToBeam()[fId] * sin(fStartAngleDetector / CLHEP::rad),
									TRexSettings::Get()->GetFBarrelDeltaESinglePosZ()[fId]);

		fFoilThickness = TRexSettings::Get()->GetFBarrelDeltaESingleFoilThickness();
	} else if(fBaseName == "SecondFBarrelDeltaESingle") {
		fDetectorLengthX = TRexSettings::Get()->GetSecondFBarrelDeltaESingleLengthX();
		fDetectorLengthY = TRexSettings::Get()->GetSecondFBarrelDeltaESingleLengthY();
		fThicknessDetector = TRexSettings::Get()->GetSecondFBarrelDeltaESingleThickness()[fId];
		fStartAngleDetector = TRexSettings::Get()->GetSecondFBarrelDeltaESingleStartAngle()[fId];

		fDeadLayer = TRexSettings::Get()->GetSecondFBarrelDeltaESingleDeadLayer();

		fPos = G4ThreeVector(TRexSettings::Get()->GetSecondFBarrelDeltaESingleDistanceToBeam()[fId] * cos(fStartAngleDetector / CLHEP::rad) + 13.0 * CLHEP::mm,// shifting -13.0 mm
									TRexSettings::Get()->GetSecondFBarrelDeltaESingleDistanceToBeam()[fId] * sin(fStartAngleDetector / CLHEP::rad),
									TRexSettings::Get()->GetSecondFBarrelDeltaESinglePosZ()[fId]);
	std::cout<<"fPos x "<<fPos.x()<<", y "<<fPos.y()<<", z "<<fPos.z()<<std::endl;
	

		fFoilThickness = TRexSettings::Get()->GetSecondFBarrelDeltaESingleFoilThickness();
	} else if(fBaseName == "MBarrelDeltaESingle") {  
		fDetectorLengthX = TRexSettings::Get()->GetMBarrelDeltaESingleLengthX();
		fDetectorLengthY = TRexSettings::Get()->GetMBarrelDeltaESingleLengthY();
		fThicknessDetector = TRexSettings::Get()->GetMBarrelDeltaESingleThickness()[fId];
		fStartAngleDetector = TRexSettings::Get()->GetMBarrelDeltaESingleStartAngle()[fId];

		fDeadLayer = TRexSettings::Get()->GetMBarrelDeltaESingleDeadLayer();

		fPos = G4ThreeVector(TRexSettings::Get()->GetMBarrelDeltaESingleDistanceToBeam()[fId] * cos(fStartAngleDetector / CLHEP::rad),
									TRexSettings::Get()->GetMBarrelDeltaESingleDistanceToBeam()[fId] * sin(fStartAngleDetector / CLHEP::rad),
									TRexSettings::Get()->GetMBarrelDeltaESinglePosZ()[fId]);

		fFoilThickness = TRexSettings::Get()->GetMBarrelDeltaESingleFoilThickness();
	} else if(fBaseName == "BBarrelDeltaESingle") {
		fDetectorLengthX = TRexSettings::Get()->GetBBarrelDeltaESingleLengthX();
		fDetectorLengthY = TRexSettings::Get()->GetBBarrelDeltaESingleLengthY();
		fThicknessDetector = TRexSettings::Get()->GetBBarrelDeltaESingleThickness()[fId];
		fStartAngleDetector = TRexSettings::Get()->GetBBarrelDeltaESingleStartAngle()[fId];

		fDeadLayer = TRexSettings::Get()->GetBBarrelDeltaESingleDeadLayer();

		fPos = G4ThreeVector(TRexSettings::Get()->GetBBarrelDeltaESingleDistanceToBeam()[fId] * cos(fStartAngleDetector / CLHEP::rad) + (13.0/3.) * CLHEP::mm, // ### target shifting -12.24 mm
									TRexSettings::Get()->GetBBarrelDeltaESingleDistanceToBeam()[fId] * sin(fStartAngleDetector / CLHEP::rad),
									TRexSettings::Get()->GetBBarrelDeltaESinglePosZ()[fId]);

		fFoilThickness = TRexSettings::Get()->GetBBarrelDeltaESingleFoilThickness();
	} else if(fBaseName == "SecondBBarrelDeltaESingle") {
		fDetectorLengthX = TRexSettings::Get()->GetSecondBBarrelDeltaESingleLengthX();
		fDetectorLengthY = TRexSettings::Get()->GetSecondBBarrelDeltaESingleLengthY();
		fThicknessDetector = TRexSettings::Get()->GetSecondBBarrelDeltaESingleThickness()[fId];
		fStartAngleDetector = TRexSettings::Get()->GetSecondBBarrelDeltaESingleStartAngle()[fId];

		fDeadLayer = TRexSettings::Get()->GetSecondBBarrelDeltaESingleDeadLayer();

		fPos = G4ThreeVector(TRexSettings::Get()->GetSecondBBarrelDeltaESingleDistanceToBeam()[fId] * cos(fStartAngleDetector / CLHEP::rad),
									TRexSettings::Get()->GetSecondBBarrelDeltaESingleDistanceToBeam()[fId] * sin(fStartAngleDetector / CLHEP::rad),
									TRexSettings::Get()->GetSecondBBarrelDeltaESinglePosZ()[fId]);

		fFoilThickness = TRexSettings::Get()->GetSecondBBarrelDeltaESingleFoilThickness();
	} else {
		std::cerr << "Direction " << fDirection << " is wrong! Use forward or backward." << std::endl;
	}

	fRotMatrix = new G4RotationMatrix;
	fRotMatrix->rotateZ(fStartAngleDetector);
}

TRexBarrelDeltaESingle::~TRexBarrelDeltaESingle() {
	// TODO Auto-generated destructor stub
}

void TRexBarrelDeltaESingle::Construct(G4LogicalVolume* experimentalHall_log, G4SDManager *SDMan) {
	// active detector
	ConstructSilicon(experimentalHall_log, SDMan);

	// PCB
	if(TRexSettings::Get()->ConstructPCB()) {
		ConstructPCB(experimentalHall_log);
	}

	// include dead layer ?
	if(fDeadLayer > 0.1*CLHEP::um) {
		ConstructDeadLayer(experimentalHall_log);
	}

	// protection foils ?
	if(fFoilThickness > 0.1*CLHEP::um) {
		ConstructFoil(experimentalHall_log);
	}
}

void TRexBarrelDeltaESingle::ConstructSilicon(G4LogicalVolume* experimentalHall_log, G4SDManager *SDMan) {
	G4Material* detectorMaterial = TRexMaterials::Get()->GetMaterial("silicon");

	fSolid = new G4Box("detector", fThicknessDetector / 2. - fDeadLayer, fDetectorLengthX / 2., fDetectorLengthY / 2.);

	fLogicalVolume = new G4LogicalVolume(fSolid, detectorMaterial, fName + "_log", 0,0 ,0);

	// does not work! -> gets local hit position wrong !!!!
	//fPhysicalVolume = new G4PVPlacement(fRotMatrix, fPos,
	//		fLogicalVolume, "BarrelDeltaE", experimentalHall_log, false, 0);
	// does work! -> gets local hit position correct !!!!
	fPhysicalVolume = new G4PVPlacement(G4Transform3D(*fRotMatrix, fPos),
			fLogicalVolume, "BarrelDeltaE", experimentalHall_log, false, 0);

	// register sensitive detector
	SDMan = G4SDManager::GetSDMpointer();
	fBarrelDeltaESingleSensitiveDetector = new TRexBarrelDeltaESingleSensitiveDetector(G4String(fName), fDirection, fId);
	SDMan->AddNewDetector(fBarrelDeltaESingleSensitiveDetector);
	fLogicalVolume->SetSensitiveDetector(fBarrelDeltaESingleSensitiveDetector);

	if(TRexSettings::Get()->Colours()) {
		fLogicalVolume->SetVisAttributes(TRexColour::Get()->yellow);
	}
}

void TRexBarrelDeltaESingle::ConstructPCB(G4LogicalVolume* experimentalHall_log) {
	G4Material* pcb = TRexMaterials::Get()->GetMaterial("pcb");

	G4double pcbWidth = 80 * CLHEP::mm;
	G4double pcbThickness = 0.6 * CLHEP::mm; // 1.6
	G4double pcbLength = 75 * CLHEP::mm;
	G4double pcbRecessWidth = 23 * CLHEP::mm;
	G4double pcbRecessLength = 22 * CLHEP::mm;
	G4double barrelDisplacementX = -12.8 * CLHEP::mm;// target shifted --> 12.8
	G4double barrelDisplacementZ = 10.25 * CLHEP::mm;// target shifted --> -10.25

	G4Box* PCB_dE = new G4Box("PCB_dE", pcbThickness / 2., pcbWidth / 2., pcbLength / 2.);

	G4Box* PCB_hole = new G4Box("PCB_hole", pcbThickness/1.9, fDetectorLengthX / 2., fDetectorLengthY / 2.);
	//G4Box* PCB_corner = new G4Box("PCB_corner", pcbThickness / 1.9, pcbRecessWidth /2. + 0.1, pcbRecessLength /2. + 0.1);commented out by Leila

	// subtract hole for detector
	G4SubtractionSolid* PCBQuadrant_solid = new G4SubtractionSolid("deltaPCB_solid", PCB_dE, PCB_hole, 0,//new G4RotationMatrix(),
			G4ThreeVector(0.0, barrelDisplacementX, barrelDisplacementZ));//G4ThreeVector(0.0, barrelDisplacement, .0));

	// subtract corner
	//PCBQuadrant_solid = new G4SubtractionSolid("deltaPCB_solid", PCBQuadrant_solid, PCB_corner, 0,//new G4RotationMatrix(), G4ThreeVector(0,  -(pcbWidth / 2.- pcbRecessWidth /2.), pcbLength /2. - pcbRecessLength/2.)); commented out by Leila

	G4LogicalVolume* PcbBarrel_log = new G4LogicalVolume(PCBQuadrant_solid, pcb, "PCBForwardBarrel_log");

	fRotMatrixPcb = fRotMatrix;

	G4ThreeVector pcbPos = fPos;
	
	pcbPos.setX(fPos.x() - barrelDisplacementX * sin(fStartAngleDetector / CLHEP::rad));
	pcbPos.setY(fPos.y() - barrelDisplacementX * cos(fStartAngleDetector / CLHEP::rad));
	pcbPos.setZ(fPos.z() + barrelDisplacementZ * CLHEP::mm);// target shifted --> - barrelDisplacementZ * CLHEP::mm

	std::cout<<"fPos x "<<fPos.x()<<", y "<<fPos.y()<<", z "<<fPos.z()<<std::endl;
	std::cout<<"x "<<pcbPos.x()<<", y "<<pcbPos.y()<<", z "<<pcbPos.z()<<std::endl;
	
	//if(fabs(fPos.x()) < 0.1) {
		
		fRotMatrixPcb->rotateX(180.*CLHEP::degree);
		pcbPos.setX(fPos.x() + barrelDisplacementX * sin(fStartAngleDetector / CLHEP::rad));
		pcbPos.setY(fPos.y() + barrelDisplacementX * cos(fStartAngleDetector / CLHEP::rad));
	
	//}

	std::cout<<"x "<<pcbPos.x()<<", y "<<pcbPos.y()<<", z "<<pcbPos.z()<<std::endl;
	
	if(fDirection == "backward") {
		
		pcbPos.setX(fPos.x() - barrelDisplacementX * sin(fStartAngleDetector / CLHEP::rad));// target shifted --> + barrelDisplacementX * CLHEP::mm
		pcbPos.setY(fPos.y() - barrelDisplacementX * cos(fStartAngleDetector / CLHEP::rad));
		pcbPos.setZ(fPos.z() - barrelDisplacementZ * CLHEP::mm);// target shifted --> + barrelDisplacementZ * CLHEP::mm

		if(fabs(fPos.y()) < 0.1) {
			fRotMatrixPcb->rotateY(180.*CLHEP::degree);
			pcbPos.setX(fPos.x() + barrelDisplacementX * sin(fStartAngleDetector / CLHEP::rad));
			pcbPos.setY(fPos.y() + barrelDisplacementX * cos(fStartAngleDetector / CLHEP::rad));
		}

		//fRotMatrixPcb->rotateZ(-180.*CLHEP::degree); // original
		fRotMatrixPcb->rotateX(-180.*CLHEP::degree);
	}
	
	std::cout<<"x "<<pcbPos.x()<<", y "<<pcbPos.y()<<", z "<<pcbPos.z()<<std::endl;
	
	// because the 2. layer is single not double as the first layer
	
	if(fBaseName == "SecondFBarrelDeltaESingle") { // added by Leila start
		
		pcbWidth = 140. * CLHEP::mm;
		pcbThickness = 0.6 * CLHEP::mm; // 1.6
		pcbLength = 140. * CLHEP::mm;
		barrelDisplacementX = -13.0 * CLHEP::mm;// before shifting: -13.3
		barrelDisplacementZ = 18.5 * CLHEP::mm;// before shifting 18.5
		
		std::cout<<" 2. forward layer x : "<<fDetectorLengthX<<" y: "<<fDetectorLengthY<<" fName: "<<fName<<std::endl;
		// 2. forward layer x : 100 y: 100 fName: SecondFBarrelDeltaESingle0
		// 2. forward layer x : 100 y: 100 fName: SecondFBarrelDeltaESingle1

		G4Box* PCB_dE2 = new G4Box("PCB_dE2", pcbThickness / 2., pcbWidth / 2., pcbLength / 2.);
		
		G4Box* PCB_hole2 = new G4Box("PCB_hole2", pcbThickness/1.9, fDetectorLengthX / 2., fDetectorLengthY / 2.);
		
		G4SubtractionSolid* PCBQuadrant_solid2 = new G4SubtractionSolid("deltaPCB_solid2", PCB_dE2, PCB_hole2, 0, G4ThreeVector(0.0, barrelDisplacementX, barrelDisplacementZ));// before shifting  + barrelDisplacementZ
		
		G4LogicalVolume* PcbBarrel_log2= new G4LogicalVolume(PCBQuadrant_solid2, pcb, "PCBForwardBarrel_log2");
		
		pcbPos.setX(fPos.x() - barrelDisplacementX * sin(fStartAngleDetector / CLHEP::rad)-2*13*CLHEP::mm);// before shifting  - barrelDisplacementX *
		pcbPos.setY(fPos.y() - barrelDisplacementX * cos(fStartAngleDetector / CLHEP::rad));
		pcbPos.setZ(fPos.z() + barrelDisplacementZ * CLHEP::mm);// before shifting + barrelDisplacementZ * CLHEP::mm
		
		//fRotMatrixPcb->rotateX(-180.*CLHEP::degree);

		if(fabs(fPos.x()) < 0.1) { // before shifting fPos.x < 0.1
		
		pcbPos.setX(fPos.x() + barrelDisplacementX * sin(fStartAngleDetector / CLHEP::rad));// before shifting  + barrelDisplacementX *
		pcbPos.setY(fPos.y() + barrelDisplacementX * cos(fStartAngleDetector / CLHEP::rad));
	
	}
	
	std::cout<<" after x SecondFBarrelDeltaESingle0 "<<pcbPos.x()<<", y "<<pcbPos.y()<<", z "<<pcbPos.z()<<std::endl;
	
	
	if(fName == "SecondFBarrelDeltaESingle0") new G4PVPlacement(fRotMatrixPcb, pcbPos, PcbBarrel_log2, "PCBForwardBarrel", experimentalHall_log, false, 0);
	if(TRexSettings::Get()->Colours()) PcbBarrel_log2->SetVisAttributes(TRexColour::Get()->darkgreen);
	
		} // added by Leila end
		
		if(fName == "SecondFBarrelDeltaESingle1"){
			
		barrelDisplacementZ = -18.5 * CLHEP::mm;// before shifting -18.5 mm
		
		G4Box* PCB_dE3 = new G4Box("PCB_dE3", pcbThickness / 2., pcbWidth / 2., pcbLength / 2.);
		G4Box* PCB_hole3 = new G4Box("PCB_hole3", pcbThickness/1.9, fDetectorLengthX / 2., fDetectorLengthY / 2.);
		G4SubtractionSolid* PCBQuadrant_solid3 = new G4SubtractionSolid("deltaPCB_solid3", PCB_dE3, PCB_hole3, 0, G4ThreeVector(0.0, barrelDisplacementX, barrelDisplacementZ));
		
		G4LogicalVolume* PcbBarrel_log3= new G4LogicalVolume(PCBQuadrant_solid3, pcb, "PCBForwardBarrel_log3");
		
		std::cout<<" before x SecondFBarrelDeltaESingle1 "<<pcbPos.x()<<", y "<<pcbPos.y()<<", z "<<pcbPos.z()<<std::endl;
			
		fRotMatrixPcb->rotateX(180.*CLHEP::degree);
		
		std::cout<<" after rotation x SecondFBarrelDeltaESingle1 "<<pcbPos.x()<<", y "<<pcbPos.y()<<", z "<<pcbPos.z()<<std::endl;
		
		pcbPos.setX(fPos.x() - barrelDisplacementX * sin(fStartAngleDetector / CLHEP::rad));// before shifting  - barrelDisplacementX *
		pcbPos.setZ(fPos.z() + barrelDisplacementZ * CLHEP::mm);

		if(fabs(fPos.x()) < 0.1) {
		
		pcbPos.setX(fPos.x() - barrelDisplacementX * sin(fStartAngleDetector / CLHEP::rad));// before shifting - barrelDisplacementX *
		pcbPos.setY(fPos.y() + barrelDisplacementX * cos(fStartAngleDetector / CLHEP::rad));
	
	}		
	
	if(fName == "SecondFBarrelDeltaESingle1") new G4PVPlacement(fRotMatrixPcb, pcbPos, PcbBarrel_log3, "PCBForwardBarrel", experimentalHall_log, false, 0);		
		if(TRexSettings::Get()->Colours()) PcbBarrel_log3->SetVisAttributes(TRexColour::Get()->darkgreen);
	
	}
		
	std::cout<<"x "<<pcbPos.x()<<", y "<<pcbPos.y()<<", z "<<pcbPos.z()<<std::endl;
	

	//G4VPhysicalVolume* phys_vol =
	
	//if(fBaseName == "FBarrelDeltaESingle") 
	new G4PVPlacement(fRotMatrixPcb, pcbPos, PcbBarrel_log, "PCBForwardBarrel", experimentalHall_log, false, 0);

	if(TRexSettings::Get()->Colours()) {
		PcbBarrel_log->SetVisAttributes(TRexColour::Get()->darkgreen);
	}
}

void TRexBarrelDeltaESingle::ConstructDeadLayer(G4LogicalVolume* experimentalHall_log) {
	G4Material* detectorMaterial = TRexMaterials::Get()->GetMaterial("silicon");

	// total detector volume
	G4Box* deadLayerWithDetector_solid = new G4Box("detectorAndDeadLayer", fThicknessDetector / 2.,
			fDetectorLengthX / 2., fDetectorLengthY / 2.);

	// subtract the active detector volume to get only the dead layer volume
	G4SubtractionSolid* deadLayer_solid = new G4SubtractionSolid("deadLayer_solid", deadLayerWithDetector_solid, fSolid);

	G4LogicalVolume* deadLayer_log = new G4LogicalVolume(deadLayer_solid, detectorMaterial, fName + "_deadLayer_log");

	//G4VPhysicalVolume* phys_vol =
	new G4PVPlacement(G4Transform3D(*fRotMatrix, fPos), deadLayer_log, "deadLayerBarreldeltaE", experimentalHall_log, false, 0);

	if(TRexSettings::Get()->Colours()) {
		fLogicalVolume->SetVisAttributes(TRexColour::Get()->yellow);
	}
}


void TRexBarrelDeltaESingle::ConstructFoil(G4LogicalVolume* experimentalHall_log) {
	G4Material* foilMaterial = TRexMaterials::Get()->GetMaterial("mylar");

	G4double foilDistance = (1.0 + 2.0) * CLHEP::mm;
	//G4double foilWidth = (TRexSettings::Get()->GetBBarrelDeltaESingleDistanceToBeam()[fId] -  1 * (foilDistance - fFoilThickness / 2.)) * 2;
	G4double foilWidth = (TRexSettings::Get()->GetBBarrelDeltaESingleDistanceToBeam()[fId] - 1 * (foilDistance - fFoilThickness / 2.)) * 2;
	//G4double foilLength = 54 * CLHEP::mm;
	G4double foilLength = 4 * CLHEP::mm;

	if(fDirection == "forward") {
		foilLength += TRexSettings::Get()->GetFBarrelDeltaESingleLengthY();
	} else if(fDirection == "middle") {
		foilLength += TRexSettings::Get()->GetMBarrelDeltaESingleLengthY();
	} else {
		foilLength += TRexSettings::Get()->GetBBarrelDeltaESingleLengthY();
	}

	G4Box* foil_solid = new G4Box("BarrelFoil_dE", fFoilThickness / 2., foilWidth / 2., foilLength / 2.);

	G4LogicalVolume* foil_log = new G4LogicalVolume(foil_solid, foilMaterial, "foilCD_log");


	// position of the foil
	G4ThreeVector position = fPos;
	position.setX(fPos.x() - (foilDistance - fFoilThickness / 2.) * cos(fStartAngleDetector / CLHEP::rad));
	position.setY(fPos.y() - (foilDistance - fFoilThickness / 2.) * sin(fStartAngleDetector / CLHEP::rad));

	//new G4PVPlacement(G4Transform3D(*fRotMatrix, position), foil_log, "foil", experimentalHall_log, false, 0); // commented out by Leila

	if(TRexSettings::Get()->Colours()) {
		foil_log->SetVisAttributes(TRexColour::Get()->darkblue);//silver);
	}
}


std::vector<ParticleMC>* TRexBarrelDeltaESingle::GetParticleMCvector() {
	std::vector<ParticleMC>* particleMCvector = new std::vector<ParticleMC>;

	particleMCvector->push_back(*fBarrelDeltaESingleSensitiveDetector->GetParticleMC());

	return particleMCvector;
}
