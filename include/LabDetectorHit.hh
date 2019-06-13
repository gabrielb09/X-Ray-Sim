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
// $Id: LabDetectorHit.hh,v 1.8 2006/06/29 17:47:53 gunter Exp $
// GEANT4 tag $Name: geant4-09-00 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef LabDetectorHit_h
#define LabDetectorHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4String.hh"
#include "G4VProcess.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class LabDetectorHit : public G4VHit
{
  public:

      	LabDetectorHit();
       ~LabDetectorHit();
      	LabDetectorHit(const LabDetectorHit&);
      	const LabDetectorHit& operator=(const LabDetectorHit&);
      	G4int operator==(const LabDetectorHit&) const;

      	inline void* operator new(size_t);
      	inline void  operator delete(void*);

      	void Draw();
      	void Print();
		void fPrint();

  public:
  		void SetTrackID  		(G4int track)      	{ trackID_ = track; };
		void SetLabID  			(G4int Lab)      	{ LabID_ = Lab; };
		void SetParentID  		(G4int parent)      { parentID_ = parent; };
		void SetTime     		(G4double t)	  	{ time_ = t; };
		void SetParticleID  	(G4int PDGid)      	{ particleID_ = PDGid; };
		void SetMomentum 		(G4ThreeVector xyz) { momentum_ = xyz; };
		void SetEnergy   		(G4double kinE)     { kEnergy_ = kinE; };
		void SetParentProcess	(G4String proc) 	{ parentProcess_ = proc; };
		void SetProcessName		(G4String proc) 	{ processName_ = proc; };
		void SetVertexPos		(G4ThreeVector xyz)	{ vertexPos_ = xyz; };
		void SetCopyNb			(G4int copy)      	{ copyNb_ = copy; };
		void SetCopyNb1			(G4int copy1)      	{ copyNb1_ = copy1; };
		void SetMaterialName	(G4String material) { materialName_ = material; };   
      	void SetEdep			(G4double de)      	{ edep_ = de; };
		void SetStepLength     	(G4double sl)      	{ stepLength_ = sl; };
      	void SetPos      		(G4ThreeVector xyz)	{ pos_ = xyz; };
/*
      G4int GetTrackID()    { return trackID; };
			G4int GetLabID()    { return LabID; };
			G4int GetParentID()    { return parentID; };
			G4double GetTime()    { return time; };
			G4int GetPartID()     { return partID; };
			G4ThreeVector GetMomentum()	{ return Momentum; };
			G4double GetEnergy()  { return kEnergy; };
			G4String GetParentProc(){ return parentProc; };
			G4String GetProcName(){ return procName; };
			G4ThreeVector GetVertex(){ return vertexPos; };
			G4int GetCopyNb()  { return copyNb; };
			G4int GetCopyNb1()  { return copyNb1; };
			G4String GetMaterialName()  { return materialName; };
      G4double GetEdep()    { return edep; };      
      G4ThreeVector GetPos(){ return pos; };
*/      
  private:
  
    	G4int        	trackID_;
		G4int        	LabID_;
		G4double	  	time_;
		G4int         	particleID_;
		G4int         	parentID_;
		G4ThreeVector 	momentum_;
		G4double      	kEnergy_;
		G4String	  	parentProcess_;
		G4String	  	processName_;
		G4ThreeVector 	vertexPos_;
		G4int         	copyNb_;
		G4int         	copyNb1_;
		G4String    	materialName_;
      	G4double      	edep_;
		G4double      	stepLength_;
      	G4ThreeVector 	pos_;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

typedef G4THitsCollection<LabDetectorHit> LabDetectorHitsCollection;

extern G4Allocator<LabDetectorHit> LabDetectorHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* LabDetectorHit::operator new(size_t)
{
  void *aHit;
  aHit = (void *) LabDetectorHitAllocator.MallocSingle();
  return aHit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void LabDetectorHit::operator delete(void *aHit)
{
  LabDetectorHitAllocator.FreeSingle((LabDetectorHit*) aHit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
