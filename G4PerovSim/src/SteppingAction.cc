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
/// \file PerovSim/src/SteppingAction.cc
/// \brief Implementation of the PervoSim::SteppingAction class

#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "G4ParticleDefinition.hh"
#include "G4OpticalPhoton.hh"
#include "G4Track.hh"
#include "G4TrackVector.hh"
#include "G4ProcessType.hh"

namespace G4PerovSim
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(EventAction* eventAction)
: fEventAction(eventAction)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{
  if (!fScoringVolume) {
    const auto detConstruction = static_cast<const DetectorConstruction*>(
      G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    fScoringVolume = detConstruction->GetScoringVolume();
  }

  // get volume of the current step
  //G4LogicalVolume* volume
  //  = step->GetPreStepPoint()->GetTouchableHandle()
  //    ->GetVolume()->GetLogicalVolume();

  // check if we are in scoring volume
  //if (volume != fScoringVolume) return;

  // want to include muons entering and exiting the crystal
  G4LogicalVolume* preVol = nullptr;
  G4LogicalVolume* postVol = nullptr;

  if (step->GetPreStepPoint()->GetTouchableHandle()->GetVolume())
    preVol = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();

  if (step->GetPostStepPoint()->GetTouchableHandle()->GetVolume())
    postVol = step->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();

  // If the particle is in or entering or exiting the scoring volume, process the step
  if (preVol != fScoringVolume && postVol != fScoringVolume) return;


  G4Track* track = step->GetTrack();
  G4String processName = step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
  //if (track->GetDefinition()->GetParticleName() == "gamma") {
    //G4cout << "Gamma interacted via: " << processName << G4endl;

  // troubleshooting muon transportation energy deposits
  G4StepPoint* preStepPoint = step->GetPreStepPoint();
  const G4VProcess* process = preStepPoint->GetProcessDefinedStep();
  G4StepPoint* postStepPoint = step->GetPostStepPoint();

  //if (process) {
    //G4cout << "Process: " << process->GetProcessName() << G4endl;
  //  G4cout << "Pre-volume: " << preStepPoint->GetTouchableHandle()->GetVolume()->GetName() << G4endl;
  //  G4cout << "Post-volume: " << postStepPoint->GetTouchableHandle()->GetVolume()->GetName() << G4endl;
    //G4cout << "Track status: " << track->GetTrackStatus() << G4endl;
  //}

  //}

    // Get secondaries generated in this step
    auto secondaries = step->GetSecondaryInCurrentStep();
    
    for (const auto& secondary : *secondaries) {
        const G4String& creatorProcess = secondary->GetCreatorProcess() ? secondary->GetCreatorProcess()->GetProcessName() : "Unknown";

        //if (creatorProcess == "Scintillation") {
            G4cout << "✨ Secondary particle ["
                   << secondary->GetParticleDefinition()->GetParticleName()
                   << "] created by " << creatorProcess
                   << " at position " << secondary->GetPosition()
                   << " with energy " << secondary->GetKineticEnergy() / CLHEP::keV << " keV"
                   << G4endl;
        //}
    }


  // collect energy deposited in this step
  G4double edepStep = step->GetTotalEnergyDeposit();
  if (edepStep > 0.) { //Debug step//
    G4Track* track = step->GetTrack();
    G4String processName = step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
    G4cout << "Edep: " << edepStep / CLHEP::keV << " keV"
             << " by " << track->GetDefinition()->GetParticleName()
             << ", Track ID: " << track->GetTrackID()
             << ", Parent ID: " << track->GetParentID()
             << ", interacted via: " << processName
             << G4endl;
    //if (track->GetParentID() == 0) {
      //G4cout << "Edep: " << edepStep
             //<< " by " << track->GetDefinition()->GetParticleName()
             //<< ", Track ID: " << track->GetTrackID()
             //<< ", Parent ID: " << track->GetParentID()
             //<< ", interacted via: " << processName << G4endl;
             //<< G4endl;
  
      // Accumulate only if it's a primary
      //fEventAction->AddEdep(edepStep);
    //}
  } 
  fEventAction->AddEdep(edepStep);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
