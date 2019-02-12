/*
 * TRexEventAction.hh
 *
 *  Created on: Jun 17, 2014
 *      Author: sklupp
 */

#ifndef TREXEVENTACTION_HH_
#define TREXEVENTACTION_HH_

#include "TRexDetectorConstruction.hh"

#include "MiniBallEventAction.hh"
#include "MiniBallHistoManager.hh"

#include "Germanium.hh"

#include "G4UserEventAction.hh"

#include <vector>

class G4Event;

class TRexEventAction : public MiniBallEventAction {
	public:
		TRexEventAction(TRexDetectorConstruction* detectorConst, MiniBallHistoManager* miniballHistoManager);
		virtual ~TRexEventAction();

		void BeginOfEventAction(const G4Event*);
		void EndOfEventAction(const G4Event*);

		void SetTree(TTree *tree);
		void CreateBranches();
		
		void setGamTotal(G4double eTotal){
			//std::cout<<"Pre Add energy: "<<gamTotalE<<std::endl;
			gamTotalE = eTotal;
			//std::cout<<"In Func energy: "<<gamTotalE<<std::endl;
			}
			
		G4double getGamTotal (){
				return gamTotalE;
			}
		

	private:
		void MiniballEndOfEventAction(const G4Event*);

		void CollectDataFromDetectors();
		//! clears all data coming from the detectors
		void ClearDetectorData();
		void CollectMiniballData();

		TRexDetectorConstruction* fDetectorConst;
		MiniBallHistoManager* fMbHistoMan;
		
		G4double gamTotalE;

		TTree* fTree;

		//std::vector<ParticleMC*> fDataOfDetectors;
		std::vector<std::vector<ParticleMC>* > fDataOfDetectors;
		std::vector<Germanium>* fMiniball;

		bool fIsHit;
		bool fMiniballHit;
		// fill always the detector tree regardless if a detector was hit or not
		// (useful to correlate the generator and the detector tree)
		bool fWriteAllEvents;
};

#endif /* TREXEVENTACTION_HH_ */
