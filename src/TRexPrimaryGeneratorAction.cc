/*
 * TRexPrimaryGeneratorAction.cc
 *
 *  Created on: Jun 15, 2014
 *      Author: sklupp
 * 
 * added option for BeamIn to allow use of TRexBeamIn primary generator
 * dhymers 2017/05/26
 */

#include "TRexPrimaryGeneratorAction.hh"
#include "TRexSettings.hh"

TRexPrimaryGeneratorAction::TRexPrimaryGeneratorAction() {
	// set the primary generator
	SetGenerator();
}

TRexPrimaryGeneratorAction::~TRexPrimaryGeneratorAction() {
	// TODO Auto-generated destructor stub
}

void TRexPrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent) {
	fCurrentGenerator->GeneratePrimaries(anEvent);
}

void TRexPrimaryGeneratorAction::SetGenerator() {
	std::string generatorName = TRexSettings::Get()->GetPrimaryGenerator();
	if(generatorName == "TestSource") {
		std::cout<<std::endl<<"Using test source ....\n"<<std::endl;
		fCurrentGenerator = new TRexTestSource;
	} else if(generatorName == "Rutherford") {
		std::cout<<std::endl<<"Using Rutherford scattering ....\n"<<std::endl;
		fCurrentGenerator = new TRexRutherford;
	} else if(generatorName == "AngularDistribution") {
		std::cout<<std::endl<<"Using given angular distribution ....\n"<<std::endl;
		fCurrentGenerator = new TRexAngularDistribution;
	} else if(generatorName == "AlphaSource") {
		std::cout<<std::endl<<"Using alpha source ....\n"<<std::endl;
		fCurrentGenerator = new TRexAlphaSource;
	} else if(generatorName == "BeamIn") {
		std::cout<<std::endl<<"Using beamIn source ....\n"<<std::endl;
		fCurrentGenerator = new TRexBeamIn;
	}
}

void TRexPrimaryGeneratorAction::SetTree(TTree* tree) {
	fCurrentGenerator->SetTree(tree);
	fCurrentGenerator->CreateTreeBranches();
}
