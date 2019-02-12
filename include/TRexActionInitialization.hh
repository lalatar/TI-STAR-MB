/*
 * ActionInitialization class used as part of migration to Geant 4.10.03
 * Used to initialize all optional user action classes. Currently set up
 * only for sequential mode.
 * 
 * dhymers 2017/05/26
 */


#ifndef TRexActionInitialization_h 
#define TRexActionInitialization_h 1

#include "G4VUserActionInitialization.hh"

#include "TRexDetectorConstruction.hh"
#include "TRexData.hh"
#include "MiniBallHistoManager.hh"

class TRexActionInitialization : public G4VUserActionInitialization
{
	public:
		// constructor accepts arguments which are used for one or more
		// user action classes
		TRexActionInitialization(TRexDetectorConstruction*, TRexData* data, MiniBallHistoManager*);
		~TRexActionInitialization();

		// for master thread
		virtual void BuildForMaster() const;

		// for worker threads, or single thread in sequential mode
		virtual void Build() const;

	private:
		// used to be directly used in user action classes, get passed
		// to constructor so they can still be used
		TRexDetectorConstruction* fDetConst;
		TRexData* fData;
		MiniBallHistoManager* fHistoMan;
};

#endif
