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
/// \file optical/LXe/src/LXeMainVolume.cc
/// \brief Implementation of the LXeMainVolume class
//
//
#include "LXeMainVolume.hh"

#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4OpticalSurface.hh"
#include "G4Sphere.hh"
#include "G4SystemOfUnits.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LXeMainVolume::LXeMainVolume(G4RotationMatrix* pRot, const G4ThreeVector& tlate,
                             G4LogicalVolume* pMotherLogical, G4bool pMany, G4int pCopyNo,
                             LXeDetectorConstruction* c)
  // Pass info to the G4PVPlacement constructor
  : G4PVPlacement(
      pRot, tlate,
      // Temp logical volume must be created here
      new G4LogicalVolume(new G4Box("temp", 1, 1, 1), G4Material::GetMaterial("Vacuum"), "temp"),
      "housing", pMotherLogical, pMany, pCopyNo),
    fConstructor(c)
{
  CopyValues();
  G4LogicalVolume* temp_log = GetLogicalVolume();

  G4double housing_x = 15;//fScint_x + fD_mtl;
  G4double housing_y = 15;//fScint_y + fD_mtl;
  G4double housing_z = 15;//fScint_z + 5. * fD_mtl; // changed from fScint_z + 2. * fD_mtl for all

  //*************************** housing and scintillator
  fScint_box = new G4Box("scint_box", fScint_x / 2., fScint_y / 2., fScint_z / 2.);
  fHousing_box = new G4Box("housing_box", housing_x / 2., housing_y / 2., housing_z / 2.);

  fScint_log = new G4LogicalVolume(fScint_box, G4Material::GetMaterial("LXe"), "scint_log");
  fHousing_log = new G4LogicalVolume(fHousing_box, G4Material::GetMaterial("Vacuum"), "housing_log"); // Changed Al to Vacuum

  auto envVisAttr = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0, 0.6));
  envVisAttr->SetForceSolid(true);
  fScint_log->SetVisAttributes(envVisAttr);
  
  // Define SiPM dimensions (half thickness)
  G4double sipmThickness = 0.68 * mm;
  G4double sipmLength = 6.02 * mm;
  G4double sipmWidth = 6.03 * mm;

  G4double pocketThickness = 0.01 * mm;

  // Coating and Gold foil
  G4double pdms_thickness = 1.0 * mm;
  G4double pdms_x = fScint_x + pdms_thickness;
  G4double pdms_y = fScint_y + pdms_thickness;
  G4double pdms_z = fScint_z + 2 * sipmThickness; //+ 2 * pdms_thickness;

  G4double gold_thickness = 0.05 * mm;
  G4double gold_x = gold_thickness;
  G4double gold_y = fScint_y;
  G4double gold_z = fScint_z;

  fPDMS_box = new G4Box("PDMS_box", pdms_x/2, pdms_y/2, pdms_z/2);
  fAg_Box = new G4Box("Ag_box", gold_x / 2, gold_y / 2, gold_z / 2);

  fPDMS_log = new G4LogicalVolume(fPDMS_box, G4Material::GetMaterial("PDMS"), "PDMS_log");
  fAg_log = new G4LogicalVolume(fAg_Box, G4Material::GetMaterial("Ag"), "Ag_log");

  G4VisAttributes* pdmsVisAttr = new G4VisAttributes(G4Colour(1.0, 0.5, 0.0, 0.75));
  pdmsVisAttr->SetVisibility(true);
  pdmsVisAttr->SetForceSolid(true);
  fPDMS_log->SetVisAttributes(pdmsVisAttr);

  G4VisAttributes* AgVisAttr = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0, 1));
  AgVisAttr->SetVisibility(true);
  AgVisAttr->SetForceSolid(true);
  fAg_log->SetVisAttributes(AgVisAttr);

  new G4PVPlacement(nullptr, G4ThreeVector(), fHousing_log, "housing_phys", temp_log, false, 0);

  new G4PVPlacement(nullptr, G4ThreeVector(), fPDMS_log, "pdms_phys", fHousing_log, false, 0);

  new G4PVPlacement(nullptr, G4ThreeVector(-fScint_x / 2 - gold_thickness / 2, 0, 0), fAg_log, "gold_negX", fPDMS_log, false, 0);
  new G4PVPlacement(nullptr, G4ThreeVector(fScint_x / 2 + gold_thickness / 2, 0, 0), fAg_log, "gold_posX", fPDMS_log, false, 1);

  // Define SiPM solid and logical volume
  fSipm_box = new G4Box("SiPM_solid", sipmLength/2, sipmWidth/2, sipmThickness/2);
  fSipm_log = new G4LogicalVolume(fSipm_box, G4Material::GetMaterial("Si"), "SiPM_log");
  G4VisAttributes* greenVisAttr = new G4VisAttributes(G4Colour(0.0, 0.5, 0.0, 0.9));
  greenVisAttr->SetVisibility(true);
  greenVisAttr->SetForceSolid(true);
  fSipm_log->SetVisAttributes(greenVisAttr);

  // Define Air pocket btwn SiPM and Crystal
  auto fPocket_box = new G4Box("Pocket_solid", sipmLength/2, sipmWidth/2, pocketThickness/2);
  auto fPocket_log = new G4LogicalVolume(fPocket_box, G4Material::GetMaterial("Air"), "Pocket_log");
  G4VisAttributes* airVisAttr = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0, 0.4));
  airVisAttr->SetVisibility(true);
  airVisAttr->SetForceSolid(true);
  fPocket_log->SetVisAttributes(airVisAttr);

  // Place the Crystal
  auto scintPhysical = new G4PVPlacement(nullptr, G4ThreeVector(), fScint_log, "scintillator", fPDMS_log, false, 0);
  
  // Position SiPM right behind crystal along z
  G4double sipmPosZ = 1.36 * mm + sipmThickness / 2;
  auto sipmPhysical = new G4PVPlacement(nullptr, G4ThreeVector(0, 0, sipmPosZ), fSipm_log, "SiPM", fPDMS_log, false, 0);

  // Position Air Pocket right between crystal and SiPM
  G4double pocketPosZ = 1.355 * mm; //1.36 * mm + pocketThickness / 2;
  auto pocketPhysical = new G4PVPlacement(nullptr, G4ThreeVector(0, 0, pocketPosZ), fPocket_log, "Pocket", fScint_log, false, 0);


  //*************** Miscellaneous sphere to demonstrate skin surfaces
  fSphere = new G4Sphere("sphere", 0., 2. * cm, 0. * deg, 360. * deg, 0. * deg, 360. * deg);
  fSphere_log = new G4LogicalVolume(fSphere, G4Material::GetMaterial("Al"), "sphere_log");
  if (fSphereOn)
    new G4PVPlacement(nullptr, G4ThreeVector(5. * cm, 5. * cm, 5. * cm), fSphere_log, "sphere",
                      fScint_log, false, 0);

  //****************** Build PMTs
  G4double innerRadius_pmt = 0.;
  G4double height_pmt = fD_mtl / 2.;
  G4double startAngle_pmt = 0.;
  G4double spanningAngle_pmt = 360. * deg;

  fPmt = new G4Tubs("pmt_tube", innerRadius_pmt, fOuterRadius_pmt, height_pmt, startAngle_pmt,
                    spanningAngle_pmt);

  // the "photocathode" is a metal slab at the back of the glass that
  // is only a very rough approximation of the real thing since it only
  // absorbs or detects the photons based on the efficiency set below
  fPhotocath = fSipm_box; //new G4Tubs("photocath_tube", innerRadius_pmt, fOuterRadius_pmt, height_pmt / 2.,
                          //startAngle_pmt, spanningAngle_pmt); //
                        

  fPmt_log = new G4LogicalVolume(fPmt, G4Material::GetMaterial("Glass"), "pmt_log");
  fPhotocath_log =  fSipm_log; //new G4LogicalVolume(fPhotocath, G4Material::GetMaterial("Si"), "housing_log"); //

  //new G4PVPlacement(nullptr, G4ThreeVector(0., 0., -height_pmt / 2.), fPhotocath_log, "photocath",
  //                  fPmt_log, false, 0);
  //***********Arrange pmts around the outside of housing**********

  G4double dx = fScint_x / fNx;
  G4double dy = fScint_y / fNy;
  G4double dz = fScint_z / fNz;

  G4double x, y, z;
  G4double xmin = -fScint_x / 2. - dx / 2.;
  G4double ymin = -fScint_y / 2. - dy / 2.;
  G4double zmin = -fScint_z / 2. - dz / 2.;
  G4int k = 0;

  z = -fScint_z / 2. - height_pmt;  // front
  PlacePMTs(fPmt_log, nullptr, x, y, dx, dy, xmin, ymin, fNx, fNy, x, y, z, k);

  auto rm_z = new G4RotationMatrix();
  rm_z->rotateY(180. * deg);
  z = fScint_z / 2. + height_pmt;  // back
  PlacePMTs(fPmt_log, rm_z, x, y, dx, dy, xmin, ymin, fNx, fNy, x, y, z, k);

  auto rm_y1 = new G4RotationMatrix();
  rm_y1->rotateY(-90. * deg);
  x = -fScint_x / 2. - height_pmt;  // left
  PlacePMTs(fPmt_log, rm_y1, y, z, dy, dz, ymin, zmin, fNy, fNz, x, y, z, k);

  auto rm_y2 = new G4RotationMatrix();
  rm_y2->rotateY(90. * deg);
  x = fScint_x / 2. + height_pmt;  // right
  PlacePMTs(fPmt_log, rm_y2, y, z, dy, dz, ymin, zmin, fNy, fNz, x, y, z, k);

  auto rm_x1 = new G4RotationMatrix();
  rm_x1->rotateX(90. * deg);
  y = -fScint_y / 2. - height_pmt;  // bottom
  PlacePMTs(fPmt_log, rm_x1, x, z, dx, dz, xmin, zmin, fNx, fNz, x, y, z, k);

  auto rm_x2 = new G4RotationMatrix();
  rm_x2->rotateX(-90. * deg);
  y = fScint_y / 2. + height_pmt;  // top
  PlacePMTs(fPmt_log, rm_x2, x, z, dx, dz, xmin, zmin, fNx, fNz, x, y, z, k);

  VisAttributes();
  SurfaceProperties();
  SetLogicalVolume(fHousing_log);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeMainVolume::CopyValues()
{
  fScint_x = fConstructor->GetScintX();
  fScint_y = fConstructor->GetScintY();
  fScint_z = fConstructor->GetScintZ();
  fD_mtl = fConstructor->GetHousingThickness();
  fNx = fConstructor->GetNX();
  fNy = fConstructor->GetNY();
  fNz = fConstructor->GetNZ();
  fOuterRadius_pmt = fConstructor->GetPMTRadius();
  fSphereOn = fConstructor->GetSphereOn();
  fRefl = fConstructor->GetHousingReflectivity();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeMainVolume::PlacePMTs(G4LogicalVolume* pmt_log, G4RotationMatrix* rot, G4double& a,
                              G4double& b, G4double da, G4double db, G4double amin, G4double bmin,
                              G4int na, G4int nb, G4double& x, G4double& y, G4double& z, G4int& k)
{
  /*  PlacePMTs : a different way to parameterize placement that does not depend
   * on calculating the position from the copy number
   *
   *  pmt_log = logical volume for pmts to be placed
   *  rot = rotation matrix to apply
   *  a,b = coordinates to vary(ie. if varying in the xy plane then pass x,y)
   *  da,db = value to increment a,b by
   *  amin,bmin = start values for a,b
   *  na,nb = number of repitions in a and b
   *  x,y,z = just pass x,y, and z by reference (the same ones passed for a,b)
   *  k = copy number to start with
   *  sd = sensitive detector for pmts
   */
  a = amin;
  for (G4int j = 1; j <= na; ++j) {
    a += da;
    b = bmin;
    for (G4int i = 1; i <= nb; ++i) {
      b += db;
      //new G4PVPlacement(rot, G4ThreeVector(x, y, z), pmt_log, "pmt", fHousing_log, false, k);
      fPmtPositions.push_back(G4ThreeVector(x, y, z));
      ++k;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeMainVolume::VisAttributes()
{
  auto housing_va = new G4VisAttributes(G4Colour(0.8, 0.8, 0.8));
  fHousing_log->SetVisAttributes(housing_va);

  auto sphere_va = new G4VisAttributes();
  sphere_va->SetForceSolid(true);
  fSphere_log->SetVisAttributes(sphere_va);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeMainVolume::SurfaceProperties()
{
  std::vector<G4double> ephoton = {7.0 * eV, 7.14 * eV};

  //**Scintillator housing properties
  std::vector<G4double> reflectivity = {fRefl, fRefl};
  std::vector<G4double> efficiency = {0.0, 0.0};
  auto scintHsngPT = new G4MaterialPropertiesTable();
  scintHsngPT->AddProperty("REFLECTIVITY", ephoton, reflectivity);
  scintHsngPT->AddProperty("EFFICIENCY", ephoton, efficiency);
  auto OpScintHousingSurface =
    new G4OpticalSurface("HousingSurface", unified, polished, dielectric_metal);
  OpScintHousingSurface->SetMaterialPropertiesTable(scintHsngPT);

  //**Sphere surface properties
  std::vector<G4double> sphereReflectivity = {1.0, 1.0};
  std::vector<G4double> sphereEfficiency = {0.0, 0.0};
  auto spherePT = new G4MaterialPropertiesTable();
  spherePT->AddProperty("REFLECTIVITY", ephoton, sphereReflectivity);
  spherePT->AddProperty("EFFICIENCY", ephoton, sphereEfficiency);
  auto OpSphereSurface = new G4OpticalSurface("SphereSurface", unified, polished, dielectric_metal);
  OpSphereSurface->SetMaterialPropertiesTable(spherePT);

  //**Photocathode surface properties
  std::vector<G4double> photocath_EFF = {1., 1.};
  std::vector<G4double> photocath_REF = {0., 0.};
  std::vector<G4double> photocath_ReR = {1.92, 1.92}; //from 1.92, 1.92
  std::vector<G4double> photocath_ImR = {1.69, 1.69}; //from 1.69, 1.69
  auto photocath_mt = new G4MaterialPropertiesTable();
  photocath_mt->AddProperty("EFFICIENCY", ephoton, photocath_EFF);
  photocath_mt->AddProperty("REALRINDEX", ephoton, photocath_ReR);
  photocath_mt->AddProperty("IMAGINARYRINDEX", ephoton, photocath_ImR);
  //photocath_mt->AddProperty("REFLECTIVITY", ephoton, photocath_REF);
  auto photocath_opsurf =
    new G4OpticalSurface("photocath_opsurf", unified, polished, dielectric_metal);
  photocath_opsurf->SetMaterialPropertiesTable(photocath_mt);

  //**Create logical skin surfaces
  //new G4LogicalSkinSurface("photocath_surf", fHousing_log, OpScintHousingSurface);
  new G4LogicalSkinSurface("sphere_surface", fSphere_log, OpSphereSurface);
  new G4LogicalSkinSurface("photocath_surf", fPhotocath_log, photocath_opsurf);
}