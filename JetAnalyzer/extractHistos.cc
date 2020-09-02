#include "TROOT.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TPaveLabel.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"

#include <iostream>

	
void extractHistos(){
	
	TCanvas * c = new TCanvas();
	c->SetCanvasSize(1200, 600);
	
	TH1::AddDirectory(kFALSE);
	
	TFile * f = TFile::Open("JetNtuple_RunIISummer16_13TeV_MC.root");
	TDirectory * d = f->GetDirectory("AK4jets");
	
	TH1D * matchDR = nullptr;
	TH1D * matchDPT = nullptr;
	TH2D * matchDRDPT = nullptr;
	TH1D * genDR = nullptr;
	TH1D * genDPhi = nullptr;
	TH1D * genDEta = nullptr;
	
	TH1D * matchPercent = nullptr;
	TH1D * matchNumber = nullptr;
	TH1D * genNumber = nullptr;
	
	d->GetObject("matchDR", matchDR);
	d->GetObject("matchDPT", matchDPT);
	d->GetObject("matchDRDPT", matchDRDPT);
	d->GetObject("genDR", genDR);
	d->GetObject("genDPhi", genDPhi);
	d->GetObject("genDEta", genDEta);
	
	d->GetObject("matchPercent", matchPercent);
	d->GetObject("matchNumber", matchNumber);
	d->GetObject("genNumber", genNumber);
	
	matchDR->Draw("HIST");
	c->Print("matchDR.png");
	c->Clear();

	matchDPT->Draw("HIST");
	c->Print("matchDPT.png");
	c->Clear();
	
	
	matchDRDPT->Draw("COLZ");
	matchDRDPT->SetStats(0);
	c->Print("matchDRDPT.png");
	c->Clear();
	
	matchDRDPT->SetMaximum(1e7);
	matchDRDPT->Draw("COLZ");
	c->Print("matchDRDPTCutoff.png");
	c->Clear();
	
	genDR->Draw("HIST");
	c->Print("genDR.png");
	c->Clear();
	
	genDPhi->Draw("HIST");
	c->Print("genDPhi.png");
	c->Clear();
	
	genDEta->Draw("HIST");
	c->Print("genDEta.png");
	c->Clear();
	
	
	matchPercent->Draw("HIST");
	c->Print("matchPercent.png");
	c->Clear();
	
	matchNumber->Draw("HIST");
	c->Print("matchNumber.png");
	c->Clear();
	
	genNumber->Draw("HIST");
	c->Print("genNumber.png");
	c->Clear();
	
	
}