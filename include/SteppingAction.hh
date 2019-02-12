#ifndef SteppingAction_h
#define SteppingAction_h 1

#include "G4UserSteppingAction.hh"
//#include map
class TRexDetectorConstruction;
class TRexEventAction;

/// Stepping action class.
///
/// In UserSteppingAction() there are collected the energy deposit and track
/// lengths of charged particles in Absober and Gap layers and
/// updated in B4aEventAction.

class SteppingAction : public G4UserSteppingAction
{
public:
  SteppingAction(const TRexDetectorConstruction* detectorConstruction,
                    TRexEventAction* eventAction);
  virtual ~SteppingAction();

  virtual void UserSteppingAction(const G4Step* step);

private:
  const TRexDetectorConstruction* fDetConstruction;
  TRexEventAction*  fEventAction;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
