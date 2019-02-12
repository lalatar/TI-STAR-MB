/*
 * TRexBeamIn.cc
 *
 *  Created on: Jun 16, 2014
 *      Author: sklupp
 * 
 * modified to build projectile based on settings file Z, A
 * dhymers 2017/06/12
 */

#include "TRexBeamIn.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "TRexSettings.hh"
#include "G4ParticleGun.hh"
#include "G4Alpha.hh"
#include "Randomize.hh"

#include "TRexSettings.hh"

using namespace CLHEP;

TRexBeamIn::TRexBeamIn() {
	fParticleGun = new G4ParticleGun(1);
}

TRexBeamIn::~TRexBeamIn() {
	// TODO Auto-generated destructor stub
}

/************************************************************
 *
 * Simulates a quadruple alpha source
 *
 ************************************************************/
void TRexBeamIn::GeneratePrimaries(G4Event *anEvent) {
	
	
	ShootReactionPosition();
	G4double fReactionEnergy = TRexSettings::Get()->GetBeamEnergy();  
	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();//Creating Projectile Particle 
	G4ParticleDefinition* beam = particleTable->GetIonTable()->GetIon(
		TRexSettings::Get()->GetProjectileZ(), //Atomic Number
		TRexSettings::Get()->GetProjectileA(), //Atomic Mass
		0                               //Excitation Energy
	);
	fParticleGun->SetParticleDefinition(beam);
	fParticleGun->SetParticlePosition(G4ThreeVector(fReactionX, fReactionY, fReactionZ));
	fParticleGun->SetParticleEnergy(fReactionEnergy);
	
	// isotropic distribution, if reaction_z < 0 shoot only to negative z's (theta < 90) otherwise only to positive z's (theta > 90)
	CreateIsotropicDistribution();

	G4ThreeVector direction;
	direction.set(0,0,1);
	direction.setMag(1.);
	direction.setTheta(fThetaCM);
	direction.setPhi(fPhi);
	fParticleGun->SetParticleMomentumDirection(direction);

	fParticleGun->GeneratePrimaryVertex(anEvent);
	//std::cout << "***********fBeamEnergy***************"<<fReactionEnergy<<std::endl;
	FillTree();
}

void TRexBeamIn::ShootReactionPosition() {
	G4double BeamDiameter = TRexSettings::Get()->GetBeamWidth() / mm;
	fReactionZ = TRexSettings::Get()->GetGasTargetLength() / cm ;
	fReactionZ = -(fReactionZ/2.) - 0.0008;
	fReactionZ *= cm;
	
	//hard-code start position
	//fReactionZ = -1 * cm; ???

	do {
		fReactionX = CLHEP::RandFlat::shoot(-BeamDiameter / 2., BeamDiameter / 2.);
		fReactionY = CLHEP::RandFlat::shoot(-BeamDiameter / 2., BeamDiameter / 2.);
	} while(sqrt(pow(fReactionX,2) + pow(fReactionY,2)) > BeamDiameter / 2.);

	fReactionX *= mm;
	fReactionY *= mm;
}

void TRexBeamIn::CreateIsotropicDistribution() {
	//fThetaCM = CLHEP::RandFlat::shoot(0., M_PI/180.)*radian;
        fThetaCM =(0.0*M_PI/180.)*radian;
	fPhi = CLHEP::RandFlat::shoot(-M_PI,M_PI)*radian;//M_PI=180 deg
	//fPhi = CLHEP::RandFlat::shoot(0., 2.* M_PI)*radian;
	//fPhi = CLHEP::RandFlat::shoot(-M_PI / 2.,M_PI + M_PI / 2.)*radian;
}

void TRexBeamIn::CreateTreeBranches() {
	if(!fTree) {
		std::cout << "\n\n\nTRexBeamIn: Tree doesn't exist!\n\n" << std::endl;
	}
	fTree->Branch("reactionEnergy", &fReactionEnergy, "reactionEnergy/D");
	fTree->Branch("reactionX", &fReactionX, "reactionX/D");
	fTree->Branch("reactionY", &fReactionY, "reactionY/D");
	fTree->Branch("reactionZ", &fReactionZ, "reactionZ/D");
	fTree->Branch("thetaCM", &fThetaCM, "thetaCM/D");
	fTree->Branch("phi", &fPhi, "phi/D");
}




