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
/// \file optical/LXe/src/LXeDetectorConstruction.cc
/// \brief Implementation of the LXeDetectorConstruction class
//
//
#include "LXeDetectorConstruction.hh"

#include "LXeDetectorMessenger.hh"
#include "LXeMainVolume.hh"
#include "LXePMTSD.hh"
#include "LXeScintSD.hh"
#include "LXeWLSSlab.hh"

#include "G4Box.hh"
#include "G4GeometryManager.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4OpticalSurface.hh"
#include "G4PVPlacement.hh"
#include "G4PhysicalConstants.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4SolidStore.hh"
#include "G4Sphere.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Tubs.hh"
#include "G4UImanager.hh"
#include "G4VisAttributes.hh"
#include "globals.hh"

G4bool LXeDetectorConstruction::fSphereOn = true;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LXeDetectorConstruction::LXeDetectorConstruction()
{
  SetDefaults();
  DefineMaterials();
  fDetectorMessenger = new LXeDetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LXeDetectorConstruction::~LXeDetectorConstruction()
{
  delete fMainVolume;
  delete fLXe_mt;
  delete fDetectorMessenger;
  delete fMPTPStyrene;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::DefineMaterials()
{
  //
  // Perovskite Crystal
  //

  G4double z, a, density;
  G4String name, symbol;
  G4int ncomponents, natoms;

  //MAPbBr3 CH3NH3PbBr3 
  a = 1.01*g/mole;

  G4Element* elH  = new G4Element(name="Hydrogen",symbol="H" , z= 1., a);

  a = 12.01*g/mole;
  G4Element* elC  = new G4Element(name="Carbon"  ,symbol="C" , z= 6., a);

  a = 14.01*g/mole;
  G4Element* elN  = new G4Element(name="Azote"  ,symbol="N" , z= 7., a);

  a = 207.2*g/mole;
  G4Element* elPb  = new G4Element(name="Lead"  ,symbol="Pb" , z= 82., a);

  a = 79.90*g/mole;
  G4Element* elBr  = new G4Element(name="Bromium"  ,symbol="Br" , z= 35., a);


  density = 3.83*g/cm3;
  G4Material* crystal_mat = new G4Material(name="MAPbBr3",density,ncomponents=6);
  crystal_mat->AddElement(elC,  natoms=1);
  crystal_mat->AddElement(elH,  natoms=3);
  crystal_mat->AddElement(elN,  natoms=1);
  crystal_mat->AddElement(elH,  natoms=3);
  crystal_mat->AddElement(elPb, natoms=1);
  crystal_mat->AddElement(elBr, natoms=3);

  //***Material properties tables

  std::vector<G4double> crystal_Energy = {7.0 * eV, 7.07 * eV, 7.14 * eV};

  std::vector<G4double> crystal_SCINT = {0.1, 1.0, 0.1};
  std::vector<G4double> crystal_RIND = {159, 1.57, 1.54};
  std::vector<G4double> crystal_ABSL = {35. * cm, 35. * cm, 35. * cm};
  crystal_prop = new G4MaterialPropertiesTable();
  crystal_prop->AddProperty("SCINTILLATIONCOMPONENT1", crystal_Energy, crystal_SCINT);
  crystal_prop->AddProperty("SCINTILLATIONCOMPONENT2", crystal_Energy, crystal_SCINT);
  crystal_prop->AddProperty("RINDEX", crystal_Energy, crystal_RIND);
  crystal_prop->AddProperty("ABSLENGTH", crystal_Energy, crystal_ABSL);
  crystal_prop->AddConstProperty("SCINTILLATIONYIELD", 12000. / MeV);
  crystal_prop->AddConstProperty("RESOLUTIONSCALE", 1.0);
  crystal_prop->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 20. * ns);
  crystal_prop->AddConstProperty("SCINTILLATIONTIMECONSTANT2", 45. * ns);
  crystal_prop->AddConstProperty("SCINTILLATIONYIELD", 12000. / MeV);
  crystal_prop->AddConstProperty("SCINTILLATIONYIELD1", 1.0);
  crystal_prop->AddConstProperty("SCINTILLATIONYIELD2", 0.0);
  crystal_mat->SetMaterialPropertiesTable(crystal_prop);

  // Set the Birks Constant for the crystal scintillator
  crystal_mat->GetIonisation()->SetBirksConstant(0.126 * mm / MeV);

  //
  // Air
  //

  auto vacuum_mt = new G4MaterialPropertiesTable();
  vacuum_mt->AddProperty("RINDEX", "Air");
  fVacuum->SetMaterialPropertiesTable(vacuum_mt);
  fAir->SetMaterialPropertiesTable(vacuum_mt);  // Give air the same rindex

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* LXeDetectorConstruction::Construct()
{
  // The experimental hall walls are all 1m away from housing walls
  G4double expHall_x = fScint_x + fD_mtl + 1. * m;
  G4double expHall_y = fScint_y + fD_mtl + 1. * m;
  G4double expHall_z = fScint_z + fD_mtl + 1. * m;

  // Create experimental hall
  fExperimentalHall_box = new G4Box("expHall_box", expHall_x, expHall_y, expHall_z);
  fExperimentalHall_log = new G4LogicalVolume(fExperimentalHall_box, fVacuum, "expHall_log");
  fExperimentalHall_phys = new G4PVPlacement(nullptr, G4ThreeVector(), fExperimentalHall_log,
                                             "expHall", nullptr, false, 0);

  fExperimentalHall_log->SetVisAttributes(G4VisAttributes::GetInvisible());

  // Place the main volume
  if (fMainVolumeOn) {
    fMainVolume =
      new LXeMainVolume(nullptr, G4ThreeVector(), fExperimentalHall_log, false, 0, this);
  }

  // Place the WLS slab
  if (fWLSslab) {
    G4VPhysicalVolume* slab =
      new LXeWLSSlab(nullptr, G4ThreeVector(0., 0., -fScint_z / 2. - fSlab_z - 1. * cm),
                     fExperimentalHall_log, false, 0, this);

    // Surface properties for the WLS slab
    auto scintWrap = new G4OpticalSurface("ScintWrap");

    new G4LogicalBorderSurface("ScintWrap", slab, fExperimentalHall_phys, scintWrap);

    scintWrap->SetType(dielectric_metal);
    scintWrap->SetFinish(polished);
    scintWrap->SetModel(glisur);

    std::vector<G4double> pp = {2.0 * eV, 3.5 * eV};
    std::vector<G4double> reflectivity = {1.0, 1.0};
    std::vector<G4double> efficiency = {0.0, 0.0};

    auto scintWrapProperty = new G4MaterialPropertiesTable();

    scintWrapProperty->AddProperty("REFLECTIVITY", pp, reflectivity);
    scintWrapProperty->AddProperty("EFFICIENCY", pp, efficiency);
    scintWrap->SetMaterialPropertiesTable(scintWrapProperty);
  }

  return fExperimentalHall_phys;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::ConstructSDandField()
{
  if (!fMainVolume) return;

  // PMT SD

  LXePMTSD* pmt = fPmt_SD.Get();
  if (!pmt) {
    // Created here so it exists as pmts are being placed
    G4cout << "Construction /LXeDet/pmtSD" << G4endl;
    auto pmt_SD = new LXePMTSD("/LXeDet/pmtSD");
    fPmt_SD.Put(pmt_SD);

    pmt_SD->InitPMTs();
    pmt_SD->SetPmtPositions(fMainVolume->GetPmtPositions());
  }
  else {
    pmt->InitPMTs();
    pmt->SetPmtPositions(fMainVolume->GetPmtPositions());
  }
  G4SDManager::GetSDMpointer()->AddNewDetector(fPmt_SD.Get());
  // sensitive detector is not actually on the photocathode.
  // processHits gets done manually by the stepping action.
  // It is used to detect when photons hit and get absorbed & detected at the
  // boundary to the photocathode (which doesn't get done by attaching it to a
  // logical volume.
  // It does however need to be attached to something or else it doesn't get
  // reset at the begining of events

  SetSensitiveDetector(fMainVolume->GetLogPhotoCath(), fPmt_SD.Get());

  // Scint SD

  if (!fScint_SD.Get()) {
    G4cout << "Construction /LXeDet/scintSD" << G4endl;
    auto scint_SD = new LXeScintSD("/LXeDet/scintSD");
    fScint_SD.Put(scint_SD);
  }
  G4SDManager::GetSDMpointer()->AddNewDetector(fScint_SD.Get());
  SetSensitiveDetector(fMainVolume->GetLogScint(), fScint_SD.Get());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::SetDimensions(G4ThreeVector dims)
{
  fScint_x = dims[0];
  fScint_y = dims[1];
  fScint_z = dims[2];
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::SetHousingThickness(G4double d_mtl)
{
  fD_mtl = d_mtl;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::SetNX(G4int nx)
{
  fNx = nx;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::SetNY(G4int ny)
{
  fNy = ny;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::SetNZ(G4int nz)
{
  fNz = nz;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::SetPMTRadius(G4double outerRadius_pmt)
{
  fOuterRadius_pmt = outerRadius_pmt;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::SetDefaults()
{
  // Resets to default values
  fD_mtl = 0.0635 * cm;

  fScint_x = 17.8 * cm;
  fScint_y = 17.8 * cm;
  fScint_z = 22.6 * cm;

  fNx = 2;
  fNy = 2;
  fNz = 3;

  fOuterRadius_pmt = 2.3 * cm;

  fSphereOn = true;
  fRefl = 1.0;

  fNfibers = 15;
  fWLSslab = false;
  fMainVolumeOn = true;
  fMainVolume = nullptr;
  fSlab_z = 2.5 * mm;

  G4UImanager::GetUIpointer()->ApplyCommand("/LXe/detector/scintYieldFactor 1.");

  if (fLXe_mt) fLXe_mt->AddConstProperty("SCINTILLATIONYIELD", 12000. / MeV);
  if (fMPTPStyrene) fMPTPStyrene->AddConstProperty("SCINTILLATIONYIELD", 10. / keV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::SetSphereOn(G4bool b)
{
  fSphereOn = b;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::SetHousingReflectivity(G4double r)
{
  fRefl = r;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::SetWLSSlabOn(G4bool b)
{
  fWLSslab = b;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::SetMainVolumeOn(G4bool b)
{
  fMainVolumeOn = b;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::SetNFibers(G4int n)
{
  fNfibers = n;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::SetMainScintYield(G4double y)
{
  fLXe_mt->AddConstProperty("SCINTILLATIONYIELD", y / MeV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::SetWLSScintYield(G4double y)
{
  fMPTPStyrene->AddConstProperty("SCINTILLATIONYIELD", y / MeV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::SetSaveThreshold(G4int save)
{
  // Sets the save threshold for the random number seed. If the number of
  // photons generated in an event is lower than this, then save the seed for
  // this event in a file called run###evt###.rndm

  fSaveThreshold = save;
  G4RunManager::GetRunManager()->SetRandomNumberStore(true);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
