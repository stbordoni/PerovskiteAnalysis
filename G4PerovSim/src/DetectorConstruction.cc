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
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

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
  G4Material* env_mat = nist->FindOrBuildMaterial("G4_AIR");

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
  
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");

  auto solidWorld = new G4Box("World",                           // its name
			      world_sizeX, world_sizeY, world_sizeZ);  // its size

  auto logicWorld = new G4LogicalVolume(solidWorld,  // its solid
    world_mat,                                       // its material
    "World");                                        // its name

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

  new G4PVPlacement(nullptr,  // no rotation
    G4ThreeVector(),          // at (0,0,0)
    logicEnv,                 // its logical volume
    "Envelope",               // its name
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
  
  
  //G4ThreeVector pos1 = G4ThreeVector(0, -20*cm, 0); // position the crystal toward the bottom of the world
  G4ThreeVector pos1 = G4ThreeVector(0, 0, 0); // position the crystal at the center o

  // Conical section shape
  G4double crystal_x =  6.03*mm, shape1_rmaxa = 2.*cm;
  G4double crystal_y =  5*mm, shape1_rmaxb = 4.*cm; //was 2.72mm, now 5mm
  G4double crystal_z =  6.02*mm;

  auto solidcrystal = new G4Box("Crystal",                           // its name
                              crystal_x, crystal_y, crystal_z);  // its size  


  auto logiccrystal = new G4LogicalVolume(solidcrystal,  // its solid
    crystal_mat,                                        // its material
    "Crystal");                                         // its name

  new G4PVPlacement(nullptr,  // no rotation
    pos1,                     // at position
    logiccrystal,              // its logical volume
    "Crystal",                 // its name
    logicEnv,                 // its mother  volume
    false,                    // no boolean operation
    0,                        // copy number
    checkOverlaps);           // overlaps checking

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
