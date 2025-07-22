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
/// \file optical/LXe/include/LXeSiPMHit.hh
/// \brief Definition of the LXeSiPMHit class
//
//
#ifndef LXeSiPMHit_h
#define LXeSiPMHit_h 1

#include "G4Allocator.hh"
#include "G4THitsCollection.hh"
#include "G4ThreeVector.hh"
#include "G4VHit.hh"
#include "G4VPhysicalVolume.hh"

class LXeSiPMHit : public G4VHit
{
    public:
        LXeSiPMHit() = default;
        //LXeSiPMHit(const LXeSiPMHit& right);
        LXeSiPMHit(G4VPhysicalVolume* pVol);
        ~LXeSiPMHit() override = default;

        LXeSiPMHit(const LXeSiPMHit& right);
        const LXeSiPMHit& operator=(const LXeSiPMHit& right);
        G4bool operator==(const LXeSiPMHit& right) const;

        inline void* operator new(size_t);
        inline void operator delete(void* aHit);

        inline void SetSiPMEdep(G4double SiPMde) { fSiPMEdep = SiPMde; }
        inline void AddSiPMEdep(G4double SiPMde) { fSiPMEdep += SiPMde; }
        inline G4double GetSiPMEdep() { return fSiPMEdep; }

        inline void IncSiPMPhotonCount() { ++fPhotonCount; }
        inline G4int GetSiPMPhotonCount() { return fPhotonCount; }

        inline void SetSiPMNumber(G4int n) { fSiPMNumber = n; }
        inline G4int GetSiPMNumber() { return fSiPMNumber; }

        inline void SetSiPMPhysVol(G4VPhysicalVolume* physVol) { fPhysVol = physVol; }
        inline G4VPhysicalVolume* GetSiPMPhysVol() { return fPhysVol; }

    private:
        G4double fSiPMEdep = 0.;
        G4int fPhotonCount = 0.;
        G4int fSiPMNumber = -1;
        G4VPhysicalVolume* fPhysVol = nullptr;
};

typedef G4THitsCollection<LXeSiPMHit> LXeSiPMHitsCollection;

extern G4ThreadLocal G4Allocator<LXeSiPMHit>* LXeSiPMHitAllocator;

inline void*LXeSiPMHit::operator new(size_t) {
  if (!LXeSiPMHitAllocator) LXeSiPMHitAllocator = new G4Allocator<LXeSiPMHit>;
  return LXeSiPMHitAllocator->MallocSingle();
}

inline void LXeSiPMHit::operator delete(void* hit) {
  LXeSiPMHitAllocator->FreeSingle((LXeSiPMHit*)hit);
}

#endif