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
// $Id: LabDetectorConstruction.hh,v 1.8 2006/06/29 17:47:30 gunter Exp $
// GEANT4 tag $Name: geant4-09-00 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef LabDetectorConstruction_h
#define LabDetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "LabDetectorSD.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"

class G4Box;
class G4Tubs;
class G4Cons;
class G4Sphere;
class G4Ellipsoid;
class G4Trd;
class G4GenericTrap;
class G4Polyhedra;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4VSolid;
class G4UnionSolid;
class G4SubtractionSolid;
class G4IntersectionSolid;
class G4Material;
class LabDetectorMessenger;
class G4AssemblyVolume;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class LabDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
	LabDetectorConstruction();
    ~LabDetectorConstruction();
	void DefineMaterials();
	void UpdateGeometry();

  public:
  
	G4VPhysicalVolume* Construct();
	G4double GetWorldFullLength()   {return WorldLength;}; 
		 
  private:
	
	G4VPhysicalVolume* ConstructLabDetector();
    
    G4Box*              SolidWorld;                 // pointer to the solid world
    G4LogicalVolume*    LVWorld;                    // pointer to the logical world
    G4VPhysicalVolume*  PVWorld;                    // pointer to the physical world
    
    G4Tubs*             SolidSiLi;                  // pointer to the solid SiLi Detector
    G4LogicalVolume*    DetectorLog;                     // pointer to the logical SiLi Detector
    G4VPhysicalVolume*  DetectorPhys;                     // pointer to the physical SiLi Detector

    G4Tubs*             TopHat;
    G4Tubs*             OuterRing;

    G4Box*              ABgroove;
    G4Box*              BCgroove;
    G4Box*              CDgroove;
    G4Box*              DEgroove;
    G4Box*              EFgroove;
    G4Box*              FGgroove;
    G4Box*              GHgroove;

    G4SubtractionSolid* Brim;
    G4UnionSolid*       Astrip;
    G4UnionSolid*       Bstrip;
    G4UnionSolid*       Cstrip;
    G4UnionSolid*       Dstrip;
    G4UnionSolid*       Estrip;
    G4UnionSolid*       Fstrip;
    G4UnionSolid*       grooves;
    G4SubtractionSolid* detector;

    G4Box*              mountPlate;
    G4GenericTrap*      cornerOne;
    G4GenericTrap*      cornerTwo;
    G4GenericTrap*      cornerThr;
    G4GenericTrap*      cornerFou;
    G4Tubs*             bigHole;
    G4Cons*             slope;
    G4Box*              bottCut;
    G4Tubs*             detectHole;
    G4Box*              topCut;
    G4UnionSolid*       topCorners;
    G4UnionSolid*       botCorners;
    G4UnionSolid*       corners;
    G4UnionSolid*       bottom;
    G4UnionSolid*       detecterHold;
    G4SubtractionSolid* cornerCut;
    G4SubtractionSolid* hole;
    G4SubtractionSolid* holeSlope;
    G4SubtractionSolid* holder;

    G4LogicalVolume*    mountLog;
    G4VPhysicalVolume*  mountPhys;

    G4Tubs*             housingCyl;
    G4Tubs*             cavity;
    G4Cons*             cone;
    G4Tubs*             sourceCav;
    G4Tubs*             ceram;
    G4Tubs*             sourceSolid;
    G4SubtractionSolid* housingCav;
    G4SubtractionSolid* housingCon;
    G4SubtractionSolid* housing;
    G4SubtractionSolid* source;

    G4LogicalVolume*    HousingLog;
    G4LogicalVolume*    SourceLog;
    G4LogicalVolume*    CeramicLog;

    G4AssemblyVolume*   CompleteSource;

    G4Box*              faradSolid;
    G4Box*              faradCavit;
    G4SubtractionSolid* faradayCage;

    G4LogicalVolume*    faradayLog;
    G4VPhysicalVolume*  faradayPhys;
    
    G4Material* Air;
    G4Material* Si;
    G4Material* Al;
    G4Material* Cu;
    G4Material* Fe;
    G4Material* Ge;
    G4Material* Pb;
    G4Material* Am;
    G4Material* StainlessSteel;
    
    G4Material* Vacuum;
    G4Material* Ceramic;
	

	G4VSensitiveDetector*	SiLiSD;
    G4VSensitiveDetector*	AlSD;
    G4VSensitiveDetector*   FaradSD;
	LabDetectorMessenger* detectorMessenger;        // pointer to the Messenger
    
    G4Material*     MaterialWorld;                  // pointer to the world material
    G4Material*     MaterialSiLi;                   // pointer to the Ge Detector material
    G4Material*     MaterialAlTop;                  // pointer to the Chamber material
    G4Material*     MaterialPbBrick;                // pointer to the Chamber material
    G4Material*     MaterialCeramicPlate;           // pointer to the Chamber material
    G4Material*     MaterialCuPlate;                // pointer to the Chamber material
    G4Material*     MaterialAlPlate;                // pointer to the Chamber material
    G4Material*     MaterialChamber;                // pointer to the Chamber material
    
    G4double WorldLength;                           // Full length of the world volume
    G4int overlap;
    
    G4double rSiLi;                                 // radius of SiLi Detector
    G4double tSiLi;                                 // thickness of SiLi Detector
  
    // G4int overlap;
    G4int CopyWorld;
    G4int CopySiLi;
    G4int CopyChamber;
    G4int CopyAlTop;
    G4int CopyPbBrick;
    G4int CopyCeramicPlate;
    G4int CopyCuPlate;
    G4int CopyAlPlate;
    
    G4VisAttributes* WhiteVisAtt;
    G4VisAttributes* RedVisAtt;
    G4VisAttributes* GreenVisAtt;
    G4VisAttributes* BlueVisAtt;
    G4VisAttributes* GrayVisAtt;
    G4VisAttributes* CyanVisAtt;
    G4VisAttributes* MagentaVisAtt;
    G4VisAttributes* YellowVisAtt;
	G4VisAttributes* TransparentVisAtt;
    G4VisAttributes* DetectorColor;
    G4VisAttributes* FaradyColor;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
