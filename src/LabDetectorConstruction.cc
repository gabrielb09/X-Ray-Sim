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
// $Id: LabDetectorConstruction.cc,v 1.19 2007/05/11 14:35:01 gcosmo Exp $
// GEANT4 tag $Name: geant4-09-00 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "LabDetectorConstruction.hh"
#include "LabDetectorMessenger.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Ellipsoid.hh"
#include "G4Trd.hh"
#include "G4Polyhedra.hh"
#include "G4Cons.hh"
#include "G4ExtrudedSolid.hh"
#include "G4GenericTrap.hh"
#include "G4LogicalVolume.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4AssemblyVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"
#include "G4RunManager.hh"

#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4ProductionCuts.hh"

#include "G4String.hh"
#include "G4UserLimits.hh"
#include "G4Colour.hh"
#include "G4ios.hh"
#include "global.h"
#include "G4SystemOfUnits.hh"

#include <vector> 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LabDetectorConstruction::LabDetectorConstruction()
:SolidWorld(0),LVWorld(0),PVWorld(0)
{
  DefineMaterials();
	detectorMessenger = new LabDetectorMessenger(this);

    extern global_struct global;
	// sprintf(global.outdir, "tmp");
	// sprintf(global.outfile, "test");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LabDetectorConstruction::~LabDetectorConstruction()
{
  delete detectorMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LabDetectorConstruction::DefineMaterials()
{
    // Material definition with NIST Manager

    G4NistManager* nistManager = G4NistManager::Instance();

    Air = nistManager->FindOrBuildMaterial("G4_AIR");
    Al = nistManager->FindOrBuildMaterial("G4_Al");
    Si = nistManager->FindOrBuildMaterial("G4_Si");
    Ge = nistManager->FindOrBuildMaterial("G4_Ge");
    Fe = nistManager->FindOrBuildMaterial("G4_Fe");
    Cu = nistManager->FindOrBuildMaterial("G4_Cu");
    Pb = nistManager->FindOrBuildMaterial("G4_Pb");
    Am = nistManager->FindOrBuildMaterial("G4_Am");

    G4Element* C = nistManager->FindOrBuildElement("C");
    G4Element* Sil = nistManager->FindOrBuildElement("Si");
    G4Element* Cr = nistManager->FindOrBuildElement("Cr");
    G4Element* Mn = nistManager->FindOrBuildElement("Mn");
    G4Element* Iro = nistManager->FindOrBuildElement("Fe");
    G4Element* Ni = nistManager->FindOrBuildElement("Ni");


    G4double a, z;
    G4double density, temperature, pressure;
    G4int nel;

    // Vacuum
    density= 2.376e-15*g/cm3;
    temperature= 300*kelvin;
    pressure= 1.0e-8*bar;
    Vacuum = new G4Material("Vacuum", density, nel=1, kStateGas,temperature,pressure);
    Vacuum-> AddMaterial(Air, 100*perCent);

    //Ceramic

    G4Element* elO = new G4Element("Oxigen","O", z=8., a = 16.00*g/mole);
    G4Element* elAl = new G4Element("Aluminum", "Al", z=13., a = 26.98*g/mole);
    density = 2.88*g/cm3;
    Ceramic = new G4Material("Ceramic",density,nel=2);
    Ceramic->AddElement(elAl,2);
    Ceramic->AddElement(elO,3);

    //Type 304 stainless steel
    G4double FractionMass;
    G4int nComponents;
    StainlessSteel = new G4Material("StainlessSteel", density = .803*g/cm3, nComponents=6);
    StainlessSteel-> AddElement(C, FractionMass = .0008);
    StainlessSteel-> AddElement(Sil, FractionMass = .0075);
    StainlessSteel-> AddElement(Cr, FractionMass = .19);
    StainlessSteel-> AddElement(Mn, FractionMass = .02);
    StainlessSteel-> AddElement(Ni, FractionMass = .1);
    StainlessSteel-> AddElement(Iro, FractionMass = .6817);


    // Print materials
    G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

G4VPhysicalVolume* LabDetectorConstruction::Construct()
{
	return ConstructLabDetector();
}
G4VPhysicalVolume* LabDetectorConstruction::ConstructLabDetector()
{
	extern global_struct global;
	overlap = global.CheckOverlap;

    //------------------------------------------------
	// Sensitive detectors
	//------------------------------------------------

    G4SDManager* SDman = G4SDManager::GetSDMpointer();
	if(!SiLiSD)
	{
		SiLiSD = new LabDetectorSD("SiLiSD");
        AlSD = new LabDetectorSD("AlSD");
        FaradSD = new LabDetectorSD("FaradSD");
		SDman->AddNewDetector(SiLiSD);
        SDman->AddNewDetector(AlSD);
        SDman->AddNewDetector(FaradSD);
	}

    //------------------------------------------------
    // Copy Number
    //------------------------------------------------

    CopySiLi = 1;
    CopyChamber = 2;
    CopyAlTop = 3;
    CopyPbBrick = 4;
    CopyCeramicPlate = 5;
    CopyCuPlate = 6;
    CopyAlPlate = 7;
    CopyWorld = -1;

    //------------------------------------------------
    // Material
    //------------------------------------------------

    MaterialWorld = Vacuum;
    MaterialSiLi = Si;
    MaterialAlTop = Al;
    MaterialPbBrick = Pb;
    MaterialCeramicPlate = Ceramic;
    MaterialCuPlate = Cu;
    MaterialAlPlate = Al;
    MaterialChamber = Fe;

	//--------- Definitions of Solids, Logical Volumes, Physical Volumes ---------

    //------------------------------------------------
    // World
    //------------------------------------------------

    WorldLength = 5.0*m;
    // G4GeometryManager::GetInstance()->SetWorldMaximumExtent(WorldLength);
    //  G4cout << "Computed tolerance = "
    //         << G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()/mm
    //         << " mm" << G4endl;

    SolidWorld = new G4Box("SolidWorld",WorldLength*0.5,WorldLength*0.5,WorldLength*0.5);
    LVWorld= new G4LogicalVolume(SolidWorld, MaterialWorld, "World", 0, 0, 0);

    //  Must place the World Physical volume unrotated at (0,0,0).
    //
    PVWorld = new G4PVPlacement(0,                  // no rotation
                                G4ThreeVector(),    // at (0,0,0)
                                LVWorld,            // its logical volume
                                "World",            // its name
                                0,                  // its mother  volume
                                false,              // no boolean operations
                                CopyWorld,          // copy number
                                true);

    //------------------------------------------------
    // SiLi
    //------------------------------------------------

    // Initial Wafer
    rSiLi = 5.05*cm;
    tSiLi = 0.25*cm;
    SolidSiLi = new G4Tubs("SolidSiLi",0.0,rSiLi,tSiLi*0.5,0*deg,360*deg);


    // Top Hat
    G4double innerHatRadius = 4.85*cm; //inner radius of the top hat
    G4double outterHatRadus = 5.06*cm; //outer raidus of the top hat
    G4double hatDepth = 0.15*cm;   //depth of the top hat

    TopHat = new G4Tubs("TopHat", innerHatRadius, outterHatRadus, hatDepth*0.5, 0, 360 * degree);
    Brim = new G4SubtractionSolid("Brim", SolidSiLi, TopHat, 0, G4ThreeVector(0,0, .05*cm)); //subtract the top hat from the wafer

    // guard ring
    G4double radiusIn = 4.50*cm; //inner radius of the guard ring
    G4double radiusOut = 4.60*cm; //outer radius of the guard ring
    G4double grooveDepth = .03*0.5*cm; //depth of the groove

    OuterRing = new G4Tubs("OuterRing", radiusIn, radiusOut, grooveDepth, 0, 360 * degree);
    //G4SubtractionSolid* GaurdRing = new G4SubtractionSolid("GaurdRing", Brim, OuterRing, 0, G4ThreeVector(0, 0, .11*cm)); //subtract the ring from the wafer

    // groove width
    G4double grooveWidth = 0.1*.5*cm;

    //--------------------------------------------------------------------------------------------------------
    // grooves all arranged as follows:
    //  define length of groove: sqrt((45.5)^2 - (x)^2) where x is the distance from the center of the wafer
    //  construct a "groove solid", a rectanugalr prism with length width and depth of a groove
    //  subtract the groove solid from the wafer
    //--------------------------------------------------------------------------------------------------------

    // groove AB
    G4double ABgrooveLen = 3.4907*cm;
    ABgroove = new G4Box("ABgroove", grooveWidth, ABgrooveLen, grooveDepth);
    Astrip = new G4UnionSolid("Astrip", OuterRing, ABgroove, 0, G4ThreeVector(-2.89*cm, 0., 0.));

    // groove BC
    G4double BCgrooveLen =  4.1287*cm;
    BCgroove = new G4Box("BCgroove", grooveWidth, BCgrooveLen, grooveDepth);
    Bstrip = new G4UnionSolid("Bstrip", Astrip, BCgroove, 0, G4ThreeVector(-1.84*cm, 0., 0.));
    // groove CD
    G4double CDgrooveLen = 4.4190*cm;
    CDgroove = new G4Box("CDgroove", grooveWidth, CDgrooveLen, grooveDepth);
    Cstrip = new G4UnionSolid("Cstrip", Bstrip, CDgroove, 0, G4ThreeVector(-0.9*cm, 0., 0.));

    // groove DE
    G4double DEgrooveLen = 4.50*cm;
    DEgroove = new G4Box("DEgroove", grooveWidth, DEgrooveLen, grooveDepth);
    Dstrip = new G4UnionSolid("Dstrip", Cstrip, DEgroove, 0, G4ThreeVector(0.0*cm, 0., 0.));


    // groove EF
    G4double EFgrooveLen = 4.4190*cm;
    EFgroove = new G4Box("EFgroove", grooveWidth, EFgrooveLen, grooveDepth);
    Estrip = new G4UnionSolid("Estrip", Dstrip, EFgroove, 0, G4ThreeVector(0.9*cm, 0., 0.));

    // groove FG
    G4double FGgrooveLen =  4.1287*cm;
    FGgroove = new G4Box("FGgroove", grooveWidth, FGgrooveLen, grooveDepth);
    Fstrip = new G4UnionSolid("Fstrip", Estrip, FGgroove, 0, G4ThreeVector(1.84*cm, 0., 0.));

    // groove GH
    G4double GHgrooveLen = 3.4907*cm;
    GHgroove = new G4Box("GHgroove", grooveWidth, GHgrooveLen, grooveDepth);
    grooves = new G4UnionSolid("grooves", Fstrip, GHgroove, 0, G4ThreeVector(2.89*cm, 0., 0.));
    
    detector = new G4SubtractionSolid("detector", Brim, grooves, 0, G4ThreeVector(0*cm, 0, .11*cm));

    // transform the solid into a logical volume
    DetectorLog = new G4LogicalVolume(detector, Si, "DetectorLog");
    // transform the solid into a physcial volume
    DetectorPhys = new G4PVPlacement(0, G4ThreeVector(0, 0, 0),DetectorLog, "DetectorPhys", LVWorld, false, CopySiLi);

    // make the detector a sensitive detector
    DetectorLog->SetSensitiveDetector(SiLiSD);

    
    //------------------------------------------------
    // X-Ray Source
    //------------------------------------------------
    housingCyl = new G4Tubs("housingCyl", 0.0, 0.3*cm, 0.4*cm, 0*deg, 360*deg); //solid cylinder

    cavity = new G4Tubs("cavity", 0.0, 0.125*cm, 0.075*cm, 0*deg, 360*deg);
    housingCav = new G4SubtractionSolid("housingCav", housingCyl, cavity, 0, G4ThreeVector(0, 0, -0.145*cm));// subtract the cylindrical part of the cavity from the hosuing
    
    cone = new G4Cons("cone",0.0,0.125*cm,0.0,0.0,0.035*cm,0*deg,360*deg);
    housingCon = new G4SubtractionSolid("housingCon", housingCav, cone, 0, G4ThreeVector(0, 0, -0.035*cm));// subtract the conical part of the cavity

    sourceCav = new G4Tubs("sourceCav", 0.0, 0.16*cm, 0.115*cm, 0*deg, 360*deg);
    housing = new G4SubtractionSolid("housing", housingCon, sourceCav, 0, G4ThreeVector(0, 0, .185*cm));// subtract the cavity for the source

    G4double radioThicc = .001*cm;// thickenss of the radiactive source
    G4double ceramRadius = 0.16*cm - radioThicc;// radius of ceramic
    G4double ceramHeight = 0.115*cm - radioThicc;// height of cermaic

    ceram  = new G4Tubs("ceram", 0.0, ceramRadius, ceramHeight, 0*deg, 360*deg);// ceramic component of source
    sourceSolid = new G4Tubs("sourceSolid", 0.0, 0.16*cm, 0.115*cm, 0*deg, 360*deg);// full cylinder of americium
    source = new G4SubtractionSolid("source", sourceSolid, ceram, 0, G4ThreeVector(0,0,0));// cut out the space for the ceramic compnenet

    HousingLog = new G4LogicalVolume(housing, StainlessSteel, "HousingLog");
    SourceLog  = new G4LogicalVolume(source, Am, "SourceLog");
    CeramicLog = new G4LogicalVolume(ceram, Ceramic, "CeramicLog");
   
    CompleteSource = new G4AssemblyVolume();//assemble the logical volumes into one complete structure

    G4RotationMatrix Rot;
    G4ThreeVector Tlate;
    G4Transform3D Trans;

    Tlate = G4ThreeVector(0,0,0);
    Trans = G4Transform3D(Rot,Tlate);
    CompleteSource->AddPlacedVolume(HousingLog, Trans);
    Tlate = G4ThreeVector(0,0,.185*cm);
    Trans = G4Transform3D(Rot,Tlate);
    CompleteSource->AddPlacedVolume(CeramicLog, Trans);
    CompleteSource->AddPlacedVolume(SourceLog,  Trans);

    Tlate = G4ThreeVector(0,0,5.4*cm);
    Rot.rotateX(180*deg);
    Trans = G4Transform3D(Rot,Tlate);
    CompleteSource->MakeImprint(LVWorld, Trans);

    //------------------------------------------------
    // Test Mount
    //------------------------------------------------
    mountPlate = new G4Box("mountPlate", 60*mm, 60*mm, 6*mm);

    std::vector<G4TwoVector> triangle(8);
    triangle[0] = G4TwoVector( 7*mm, 7*mm);
    triangle[1] = G4TwoVector( 7*mm,-7*mm);
    triangle[2] = G4TwoVector(-7*mm,-7*mm);
    triangle[3] = G4TwoVector(-7*mm, 7*mm);    
    triangle[4] = G4TwoVector( 7*mm, 7*mm);
    triangle[5] = G4TwoVector( 7*mm,-7*mm);
    triangle[6] = G4TwoVector(-7*mm,-7*mm);
    triangle[7] = G4TwoVector(-7*mm, 7*mm); 

    triangle[2] = G4TwoVector(-7*mm, 7*mm);
    triangle[6] = G4TwoVector(-7*mm, 7*mm);
    cornerOne = new G4GenericTrap("cornerOne", 7*mm, triangle);
    triangle[2] = G4TwoVector(-7*mm,-7*mm);
    triangle[6] = G4TwoVector(-7*mm,-7*mm);

    triangle[1] = G4TwoVector(-7*mm,-7*mm); 
    triangle[5] = G4TwoVector(-7*mm,-7*mm); 
    cornerTwo = new G4GenericTrap("cornerTwo", 7*mm, triangle);
    triangle[1] = G4TwoVector( 7*mm,-7*mm); 
    triangle[5] = G4TwoVector( 7*mm,-7*mm); 

    triangle[0] = G4TwoVector( 7*mm, -7*mm);
    triangle[4] = G4TwoVector( 7*mm, -7*mm);
    cornerThr = new G4GenericTrap("cornerThr", 7*mm, triangle);
    triangle[0] = G4TwoVector( 7*mm, 7*mm);
    triangle[4] = G4TwoVector( 7*mm, 7*mm);

    triangle[3] = G4TwoVector( 7*mm, 7*mm);
    triangle[7] = G4TwoVector( 7*mm, 7*mm);
    cornerFou = new G4GenericTrap("cornerFou", 7*mm, triangle);
    triangle[3] = G4TwoVector(-7*mm, 7*mm);
    triangle[7] = G4TwoVector(-7*mm, 7*mm);

    topCorners = new G4UnionSolid("topCorners", cornerOne, cornerTwo, 0, G4ThreeVector(-108*mm, 0., 0.));
    botCorners = new G4UnionSolid("botCorners", cornerFou, cornerThr, 0, G4ThreeVector(-108*mm, 0., 0.));
    corners = new G4UnionSolid("corners", topCorners, botCorners, 0, G4ThreeVector(0., -108*mm, 0.));

    cornerCut = new G4SubtractionSolid("cornerCut", mountPlate, corners, 0, G4ThreeVector(54*mm, 54*mm, 0.));

    bigHole = new G4Tubs("bigHole", 0, 47.5*mm, 7*mm, 0, 360*deg);
    hole = new G4SubtractionSolid("hole", cornerCut, bigHole, 0, G4ThreeVector(0., 0., 0.));

    slope = new G4Cons("slope", 0, 54*mm, 0, 47.5*mm, 2.5*mm, 0, 360*deg);
    bottCut = new G4Box("bottCut", 35*mm, 35*mm, 2.25*mm);
    bottom = new G4UnionSolid("bottom", slope, bottCut, 0, G4ThreeVector(60*mm, 0, -0.125*mm));
    holeSlope = new G4SubtractionSolid("holeSlope", hole, bottom, 0, G4ThreeVector(0, 0, -3.5*mm));

    detectHole = new G4Tubs("detectHole", 0, 57*mm, 3.50*mm, 0, 360*deg);
    topCut = new G4Box("topCut", 4*mm,4*mm, 3.50*mm);
    detecterHold = new G4UnionSolid("detecterHold", detectHole, topCut, 0, G4ThreeVector(-58.5*mm, 0., 0.));
    holder = new G4SubtractionSolid("holder", holeSlope, detecterHold, 0, G4ThreeVector(0., 0., 2.5*mm));

    //G4SubtractionSolid* topNotch = new G4SubtractionSolid("topNotch", holder, topCut, 0, G4ThreeVector(-58.5, 0, 2.5*mm));

    // transform the solid into a logical volume
    mountLog = new G4LogicalVolume(holder, Al, "mountLog");
    // transform the solid into a physcial volume
    mountPhys = new G4PVPlacement(0, G4ThreeVector(0., 0., -0.25*mm), mountLog, "mountPhys", LVWorld, false, 2);

    //------------------------------------------------
    // Faraday Cage
    //------------------------------------------------
    faradSolid = new G4Box("faradSolid", 6.01*cm, 6.01*cm, 2.8175*cm);
    faradCavit = new G4Box("faradCavit", 6.00*cm, 6.00*cm, 2.8075*cm);
    faradayCage = new G4SubtractionSolid("faradayCage", faradSolid, faradCavit, 0, G4ThreeVector(0., 0., 0.));

    faradayLog = new G4LogicalVolume(faradayCage, Al, "faradayLog");
    faradayPhys = new G4PVPlacement(0, G4ThreeVector(0., 0., 2.1825*cm), faradayLog, "faradayPhys", LVWorld, false, 3);

    faradayLog->SetSensitiveDetector(FaradSD);

    //------------------------------------------------
    // Visualization attributes
    //------------------------------------------------

    WhiteVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0,1.0));
    RedVisAtt= new G4VisAttributes(G4Colour(1.0,0.0,0.0,1.0));
    GreenVisAtt= new G4VisAttributes(G4Colour(0.0,1.0,0.0,1.0));
    BlueVisAtt= new G4VisAttributes(G4Colour(0.0,0.0,1.0,1.0));
    GrayVisAtt= new G4VisAttributes(G4Colour(0.5,0.5,0.5,1.0));
    CyanVisAtt= new G4VisAttributes(G4Colour(0.0,1.0,1.0,1.0));
    MagentaVisAtt= new G4VisAttributes(G4Colour(1.0,0.0,1.0,1.0));
    YellowVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,0.0,1.0));
    TransparentVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,1.0,0.0));
    DetectorColor = new G4VisAttributes(G4Colour(0.8, 0.75, 0.33, 1.0));
    FaradyColor = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0, 0.0));

    LVWorld->SetVisAttributes(TransparentVisAtt);
    //LVChamber->SetVisAttributes(YellowVisAtt);
    DetectorLog->SetVisAttributes(DetectorColor);
    //LVCuPlate->SetVisAttributes(RedVisAtt);
    //LVAlPlate->SetVisAttributes(GreenVisAtt);
    //LVCeramicPlate->SetVisAttributes(GrayVisAtt);
    SourceLog->SetVisAttributes(RedVisAtt);
    CeramicLog->SetVisAttributes(WhiteVisAtt);
    HousingLog->SetVisAttributes(GrayVisAtt);
    mountLog->SetVisAttributes(GrayVisAtt);
    faradayLog->SetVisAttributes(FaradyColor);
	return PVWorld;
}

void LabDetectorConstruction::UpdateGeometry()
{
	// Cleanup old geometry
	G4GeometryManager::GetInstance()->OpenGeometry();
	G4PhysicalVolumeStore::GetInstance()->Clean();
	G4LogicalVolumeStore::GetInstance()->Clean();
	G4SolidStore::GetInstance()->Clean();

	G4RunManager::GetRunManager()->DefineWorldVolume(ConstructLabDetector());
	G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
