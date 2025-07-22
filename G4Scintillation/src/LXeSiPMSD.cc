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
/// \file optical/LXe/src/LXeSiPMSD.cc
/// \brief Implementation of the LXeSiPMSD class
//
//
#include "LXeSiPMSD.hh"
#include "LXeSiPMHit.hh"

#include "LXeDetectorConstruction.hh"
#include "LXeUserTrackInformation.hh"

#include "G4LogicalVolume.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4SDManager.hh"
#include "G4Step.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VProcess.hh"
#include "G4VTouchable.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LXeSiPMSD::LXeSiPMSD(G4String name) : G4VSensitiveDetector(name)
{
  collectionName.insert("sipmCollection");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void LXeSiPMSD::Initialize(G4HCofThisEvent* hitsCE)
{
  fSiPMCollection = new LXeSiPMHitsCollection(SensitiveDetectorName, collectionName[0]);

  if (fHitsCID < 0) {
    fHitsCID = G4SDManager::GetSDMpointer()->GetCollectionID(fSiPMCollection);
  }
  hitsCE->AddHitsCollection(fHitsCID, fSiPMCollection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool LXeSiPMSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
  G4double SiPMedep = aStep->GetTotalEnergyDeposit();
  if (SiPMedep == 0.) return false;  // No edep so don't count as hit

  G4StepPoint* thePrePoint = aStep->GetPreStepPoint();
  auto theTouchable = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
  G4VPhysicalVolume* thePrePV = theTouchable->GetVolume();

  G4StepPoint* thePostPoint = aStep->GetPostStepPoint();

  auto sipmHit = new LXeSiPMHit(thePrePV);

  sipmHit->SetSiPMEdep(SiPMedep);

  fSiPMCollection->insert(sipmHit);

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool LXeSiPMSD::ProcessHits_boundary(const G4Step* aStep, G4TouchableHistory*)
{
  // need to know if this is an optical photon
  if (aStep->GetTrack()->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition())
    return false;

  // User number 1 since we only have one SiPM
  G4int SiPMNumber = 1;
  G4VPhysicalVolume* physVol = aStep->GetPostStepPoint()->GetTouchable()->GetVolume(1);

  // Find the correct hit collection
  size_t n = fSiPMCollection->entries();
  LXeSiPMHit* hit = nullptr;
  for (size_t i = 0; i < n; ++i) {
    if ((*fSiPMCollection)[i]->GetSiPMNumber() == SiPMNumber) {
      hit = (*fSiPMCollection)[i];
      break;
    }
  }

  if (hit == nullptr) {  // the SiPM wasn't previously hit in this event
    hit = new LXeSiPMHit();  // so create new hit
    hit->SetSiPMNumber(SiPMNumber);
    hit->SetSiPMPhysVol(physVol);
    fSiPMCollection->insert(hit);
    //hit->SetSiPMPos((*fSiPMPositionsX)[pmtNumber], (*fPMTPositionsY)[pmtNumber],
    //               (*fPMTPositionsZ)[pmtNumber]);
  }

  hit->IncSiPMPhotonCount();  // increment hit for the selected SiPM

  /*if (!LXeDetectorConstruction::GetSphereOn()) {
    hit->SetDrawit(true);
    // If the sphere is disabled then this hit is automaticaly drawn
  }
  else {  // sphere enabled
    auto trackInfo = (LXeUserTrackInformation*)aStep->GetTrack()->GetUserInformation();
    if (trackInfo->GetTrackStatus() & hitSphere)
      // only draw this hit if the photon has hit the sphere first
      hit->SetDrawit(true);
  }*/

  return true;
}