/*
 * TRexDetectorConstruction.hh
 *
 *  Created on: Jun 15, 2014
 *      Author: sklupp
 */

#ifndef TREXDETECTORCONSTRUCTION_HH_
#define TREXDETECTORCONSTRUCTION_HH_

#include "TRexBaseDetector.hh"
#include "TRexCD.hh"
#include "TRexCDdeltaE.hh"
//#include "TRexCDdeltaESingle.hh"
#include "TRexCDErest.hh"
//#include "TRexCDErestSingle.hh"
#include "TRexBarrel.hh"
#include "TRexBarrelDeltaE.hh"
#include "TRexBarrelErest.hh"
//#include "TRexTrapezoidDeltaE.hh"
#include "TRexTrapezoid.hh"
#include "TRexVacuumChamber.hh"
#include "TRexVacuumChamberGas.hh"

#include "MiniBallDetectorArray.hh"
#include "MiniBallMaterial.hh"
#include "MiniBallWorld.hh"

#include "G4VUserDetectorConstruction.hh"

#include <vector>

class G4Box;
class G4Tubs;
class G4VPhysicalVolume;
class G4LogicalVolume;

class TRexDetectorConstruction : public G4VUserDetectorConstruction {
	public:
		TRexDetectorConstruction();
		TRexDetectorConstruction(MiniBallHistoManager* miniballHistoManager);
		virtual ~TRexDetectorConstruction();

		G4VPhysicalVolume* Construct();

		//void SetTree(TTree *tree);

		std::vector<TRexBaseDetector*> GetAllDetectors() {return fDetectors; }

	private:
		void ConstructExperimentalHall();
		void ConstructTarget();
		void ConstructVacuumChamber();
		void DefineDetectors();
		void DefineMiniball();

		// Miniball histo manager (needed to make Miniball a sensitive detector)
		MiniBallHistoManager* fMiniballHistoManager;

		// experimental hall
		G4Box* fExperimentalHall_box;
		G4LogicalVolume* fExperimentalHall_log;
		G4VPhysicalVolume* fExperimentalHall_phys;

		// stores all detectors and their position
		std::vector<TRexBaseDetector*> fDetectors;

		// sensitive detector manager
		G4SDManager* fSDMan;

		MiniBallDetectorArray* fMiniball;

		// vacuum chamber
		TRexVacuumChamber* fVacuumChamber;

		// vacuum chamber gas
		TRexVacuumChamberGas* fVacuumChamberGas;
		G4LogicalVolume* fVacuumChamberGas_log;
};

#endif /* TREXDETECTORCONSTRUCTION_HH_ */
