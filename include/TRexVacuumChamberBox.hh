/*
 * TRexVacuumChamberCylinder.hh
 *
 *  Created on: Jun 16, 2014
 *      Author: sklupp
 */

#ifndef TREXVACUUMCHAMBERBOX_HH_
#define TREXVACUUMCHAMBERBOX_HH_

#include "TRexVacuumChamber.hh"
#include "TRexMaterials.hh"
#include "TRexColour.hh"

#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

class TRexVacuumChamberBox : public TRexVacuumChamber {
	public:
		TRexVacuumChamberBox();
		virtual ~TRexVacuumChamberBox();

		void ConstructChamber(G4LogicalVolume* experimentalHall_log);

	private:

};

#endif /* TREXVACUUMCHAMBERBOX_HH_ */
