//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4UserSteppingAction.hh 95247 2016-02-02 09:36:27Z gcosmo $
//
//---------------------------------------------------------------
//
//  G4UserSteppingAction.hh
//
// class description:
//    This class represents actions taken place by the user at each
//    end of stepping. 
//
// Contact:
//   Questions and comments to this code should be sent to
//     Katsuya Amako  (e-mail: Katsuya.Amako@kek.jp)
//     Takashi Sasaki (e-mail: Takashi.Sasaki@kek.jp)
//
//---------------------------------------------------------------
 
#ifndef MySteppingAction_h
#define MySteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "G4ParticleTable.hh"

class G4Step;
class G4SteppingManager;            // Forward declaration
 
///////////////////////////
class MySteppingAction : public G4UserSteppingAction

{ 

   /*public:
 
   G4UserSteppingAction(); 
   virtual ~G4UserSteppingAction();*/
    
	public:
	
	MySteppingAction();
	~MySteppingAction();
	    
   
// Member functions
    //virtual void SetSteppingManagerPointer(G4SteppingManager* pValue);
    virtual void UserSteppingAction(const G4Step*);
    
    //void MySteppingAction(const G4Step*);
    private:
    
    FILE* myFile;
    
    //const G4Step* s;
 
 //-----------
     protected:
 //-----------
  
 // Member data
    G4SteppingManager* fpSteppingManager;
 
 };
 
 #endif
 
 
