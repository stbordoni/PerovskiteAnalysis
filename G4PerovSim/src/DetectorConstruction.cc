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
/// \file PerovSim/src/DetectorConstruction.cc copied B1/src/DetectorConstruction.cc
/// \brief Implementation of the PerovSim::DetectorConstruction class

#include "DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4OpticalSurface.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4VisAttributes.hh"

#include "G4UserLimits.hh"

namespace G4PerovSim
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();

  // Envelope parameters
  //
  G4double env_sizeXY = 10*cm, env_sizeZ = 10*cm; //was 40x40x40
  G4Material* env_mat = nist->FindOrBuildMaterial("G4_Galactic"); // changed from "G4_AIR" to "G4_Galactic"

  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  //
  // World
  //
  
  //G4double world_sizeXY = 1.2*env_sizeXY;
  //G4double world_sizeZ  = 1.2*env_sizeZ;
  G4double world_sizeX  = 50*cm;
  G4double world_sizeY  = 50*cm;
  G4double world_sizeZ  = 50*cm;
  
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_Galactic"); // changed from "G4_AIR" to "G4_Galactic"

  auto solidWorld = new G4Box("World",                           // its name
			      world_sizeX, world_sizeY, world_sizeZ);  // its size

  auto logicWorld = new G4LogicalVolume(solidWorld,  // its solid
    world_mat,                                       // its material
    "World");                                        // its name

  auto worldVisAttr = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0, 0.1)); // White with 10% opacity
  worldVisAttr->SetForceSolid(true); // Optional: show as solid rather than wireframe
  logicWorld->SetVisAttributes(worldVisAttr);

  auto physWorld = new G4PVPlacement(nullptr,  // no rotation
    G4ThreeVector(),                           // at (0,0,0)
    logicWorld,                                // its logical volume
    "World",                                   // its name
    nullptr,                                   // its mother  volume
    false,                                     // no boolean operation
    0,                                         // copy number
    checkOverlaps);                            // overlaps checking

  
  //
  // Envelope
  //

  auto solidEnv = new G4Box("Envelope",                    // its name
    env_sizeXY, env_sizeXY, env_sizeZ);  // its size

  auto logicEnv = new G4LogicalVolume(solidEnv,  // its solid
    env_mat,                                     // its material
    "Envelope");                                 // its name

  auto envVisAttr = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0, 0.1));
  envVisAttr->SetForceSolid(true);
  logicEnv->SetVisAttributes(envVisAttr);

  new G4PVPlacement(nullptr,  // no rotation
    G4ThreeVector(),          // at (0,0,0)
    logicEnv,                 // its logical volume
    "Envelope",               // its name
    logicWorld,               // its mother  volume
    false,                    // no boolean operation
    0,                        // copy number
    checkOverlaps);           // overlaps checking
  
  
  //
  // SiPM
  //

  G4double SiPM_sizeX  = 20*mm;
  G4double SiPM_sizeY  = 2.5*mm;
  G4double SiPM_sizeZ  = 20*mm;

  G4Material* SiPM_mat = nist->FindOrBuildMaterial("G4_Si");

  auto solidSiPM = new G4Box("SiPM",                    // its name
    SiPM_sizeX, SiPM_sizeY, SiPM_sizeZ);  // its size

  auto logicSiPM = new G4LogicalVolume(solidSiPM,  // its solid
    SiPM_mat,                                     // its material
    "SiPM");                                 // its name

  G4VisAttributes* sipmVisAttr = new G4VisAttributes(G4Colour(0.0, 0.4, 0.0, 0.7)); // green, 70% opaque
  sipmVisAttr->SetForceSolid(true);
  logicSiPM->SetVisAttributes(sipmVisAttr);

  auto physSiPM = new G4PVPlacement(nullptr,  // no rotation
    G4ThreeVector(0,-6.25*mm,0),          // at (0,0,0)
    logicSiPM,                 // its logical volume
    "SiPM",               // its name
    logicWorld,               // its mother  volume
    false,                    // no boolean operation
    0,                        // copy number
    checkOverlaps);           // overlaps checking


  //
  // Crystal
  //
  //G4Material* crystal_mat = nist->FindOrBuildMaterial("G4_A-150_TISSUE");
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

  //
  // Scintillator Properties Tables
  //
  
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
  
  //G4ThreeVector pos1 = G4ThreeVector(0, -20*cm, 0); // position the crystal toward the bottom of the world
  G4ThreeVector pos1 = G4ThreeVector(0, 0, 0); // position the crystal at the center o

  // Rectangular prism section shape
  G4double crystal_x =  20*mm, shape1_rmaxa = 2.*cm; //was 6.03mm, changed to 20mm
  G4double crystal_y =  5*mm, shape1_rmaxb = 4.*cm; //was 2.72mm, changed to 5mm
  G4double crystal_z =  20*mm; //was 6.02mm, changed to 20mm for testing

  auto solidcrystal = new G4Box("Crystal",                           // its name
                              crystal_x, crystal_y, crystal_z);  // its size  


  auto logiccrystal = new G4LogicalVolume(solidcrystal,  // its solid
    crystal_mat,                                        // its material
    "Crystal");                                         // its name

  G4VisAttributes* crystalVisAttr = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0, 0.3));
  crystalVisAttr->SetForceSolid(true);
  logiccrystal->SetVisAttributes(crystalVisAttr);

  auto physcrystal = new G4PVPlacement(nullptr,  // no rotation
    pos1,                     // at position
    logiccrystal,              // its logical volume
    "Crystal",                 // its name
    logicEnv,                 // its mother  volume
    false,                    // no boolean operation
    0,                        // copy number
    checkOverlaps);           // overlaps checking


  //
  // Optical Interface
  //

  G4OpticalSurface* sipmSurface = new G4OpticalSurface("SiPMOpticalSurface");
  sipmSurface->SetType(dielectric_metal);
  sipmSurface->SetFinish(polished);
  sipmSurface->SetModel(unified);

  G4MaterialPropertiesTable* sipmMPT = new G4MaterialPropertiesTable();

  // Wavelength range
  const G4int num = 2;
  G4double photonEnergy[num] = {1.5*eV, 3.5*eV};  // ~350–830 nm
  G4double efficiency[num] = {0.2, 0.2};
  G4double reflectivity[num] = {0.0, 0.0};

  sipmMPT->AddProperty("EFFICIENCY", photonEnergy, efficiency, num);
  sipmMPT->AddProperty("REFLECTIVITY", photonEnergy, reflectivity, num);
  sipmSurface->SetMaterialPropertiesTable(sipmMPT);

  new G4LogicalBorderSurface("SiPMBorder",
                           physcrystal,  // Volume where photon comes from
                           physSiPM,          // Volume photon enters
                           sipmSurface);      // The surface to apply

  
  // After defining your logical volume for the crystal:
  G4double maxStep = 0.1 * mm;  // or even smaller
  logiccrystal->SetUserLimits(new G4UserLimits(maxStep));

  // Set Shape1 as scoring volume  // is the envelope necessary? the scoring volume can be the world? 
  //
  fScoringVolume = logiccrystal; //changed logicEnv to logiccrystal *****//

  //
  //always return the physical World
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
