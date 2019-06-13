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
// $Id: LabRunAction.cc,v 1.9 2006/06/29 17:48:16 gunter Exp $
// GEANT4 tag $Name: geant4-09-00 $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "LabRunAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
/*
#include "G4ios.hh"
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
*/
#include "global.h"
#include <fstream>
#include <TTree.h>
#include <TFile.h>

using namespace std;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LabRunAction::LabRunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LabRunAction::~LabRunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LabRunAction::BeginOfRunAction(const G4Run* aRun)
{
  extern global_struct global;
  char fname[100];
	
  ((G4Run *)(aRun))->SetRunID(global.runnum);
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;

    if(global.OutputFormat == 1)
    {
        sprintf(fname, "%s/%s.root", global.outdir, global.outfile );
        global.fROOT=new TFile(fname,"RECREATE");
        global.tree = new TTree("tree", "");
        
        global.tree -> Branch("eventID",&global.eventID,"eventID/I");
        global.tree -> Branch("trackID",&global.trackID,"trackID/I");
        global.tree -> Branch("parentID",&global.parentID,"parentID/I");
        global.tree -> Branch("particleID",&global.particleID,"particleID/I");
        global.tree -> Branch("copyNb1",&global.copyNb1,"copyNb1/I");
        global.tree -> Branch("copyNb",&global.copyNb,"copyNb/I");
        global.tree -> Branch("material",global.material,"material[10]/C");
        global.tree -> Branch("processName",global.processName,"processName[30]/C");
        global.tree -> Branch("parentProcess",global.parentProcess,"parentProcess[30]/C");
        global.tree -> Branch("time",&global.time,"time/F");
        global.tree -> Branch("energy",&global.energy,"energy/F");
        global.tree -> Branch("eDep",&global.eDep,"eDep/F");
        global.tree -> Branch("px",&global.px,"px/F");
        global.tree -> Branch("py",&global.py,"py/F");
        global.tree -> Branch("pz",&global.pz,"pz/F");
        global.tree -> Branch("stepLength",&global.stepLength,"stepLength/F");
        global.tree -> Branch("x",&global.x,"x/F");
        global.tree -> Branch("y",&global.y,"y/F");
        global.tree -> Branch("z",&global.z,"z/F");
    }
    if(global.OutputFormat == 0)
    {
        sprintf(fname, "%s/%s.dat", global.outdir, global.outfile );
        G4cout << "Output file: " << fname << G4endl;
        global.output.open (fname);
    }
    if(global.SimulationType >= 2 && global.SimulationType != 5) // simulation for pbar/dbar stop event
    {
        global.inputX.clear();
        global.inputY.clear();
        global.inputZ.clear();
        sprintf(fname, "%s/%s", global.indir, global.infile );
        G4cout << "Input file for primary particles: " << fname << G4endl;
        global.input.open (fname);
        // Read input files for primary particles
        float X,Y,Z;
        int nLine = 0;
        while (global.input.good())
        {
            global.input >> X >> Y >> Z;
            //				cout << X << " " << Y << " " <<  Z << endl;
            global.inputX.push_back(X);
            global.inputY.push_back(Y);
            global.inputZ.push_back(Z);
            nLine++;
        }
        cout << "there are " << nLine << " lines in the input file" << endl;
        global.nLine = nLine;
        global.input.close();
    }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LabRunAction::EndOfRunAction(const G4Run*)
{
	extern global_struct global;
	if(global.OutputFormat == 0) global.output.close();
	if(global.OutputFormat == 1)
	{
		global.fROOT->Print();
		global.fROOT->Write();
		global.fROOT->Close();
	}
	G4cout << "Run end  " << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



