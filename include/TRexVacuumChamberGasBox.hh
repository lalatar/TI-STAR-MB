/*
 * TRexVacuumChamberGasBox.hh
 *
 *  Created on: March 28, 2018
 *      Author: latar
 */

#ifndef TREXVACUUMCHAMBERGASBOX_HH_
#define TREXVACUUMCHAMBERGASBOX_HH_

#include "TRexVacuumChamberGas.hh"
#include "TRexMaterials.hh"
#include "TRexColour.hh"

#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

class TRexVacuumChamberGasBox : public TRexVacuumChamberGas {
	public:
		TRexVacuumChamberGasBox();
		~TRexVacuumChamberGasBox();

		G4LogicalVolume* ConstructChamberGas(G4LogicalVolume* experimentalHall_log);

	private:

};

#endif /* TREXVACUUMCHAMBERGASBOX_HH_ */
