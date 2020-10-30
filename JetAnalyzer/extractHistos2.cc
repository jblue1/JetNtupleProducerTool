#include "TROOT.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TPaveLabel.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"

#include <iostream>

	
void extractHistos2(){
	
	TCanvas * c = new TCanvas();
	c->SetCanvasSize(1200, 600);
	
	TH1::AddDirectory(kFALSE);
	
	TFile * f = TFile::Open("JetNtuple_RunIISummer16_13TeV_MC.root");
	TDirectory * d = f->GetDirectory("AK4jets");
	
	TH1D * genJetPT = nullptr;
	TH1D * recoJetPT = nullptr;
	TH1D * algJetPT = nullptr;
	TH1D *  fullAlgJetPT = nullptr;

	
	d->GetObject("genJetPT", genJetPT);
	d->GetObject("recoJetPT", recoJetPT);
	d->GetObject("algJetPT", algJetPT);
	d->GetObject("fullAlgJetPT", fullAlgJetPT);

	
	genJetPT->Draw("HIST");
	c->Print("genJetPT.png");
	c->Clear();

	recoJetPT->Draw("HIST");
	c->Print("recoJetPT.png");
	c->Clear();
	
	algJetPT->Draw("HIST");
	c->Print("algJetPT.png");
	c->Clear();
	
	fullAlgJetPT->Draw("HIST");
	c->Print("fullAlgJetPT.png");
	c->Clear();
	
	
	

	
	
}
