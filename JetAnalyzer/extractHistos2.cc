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

	TH1D * genJetPT50_100 = nullptr;
	TH1D * recoJetPT50_100 = nullptr;
	TH1D * algJetPT50_100 = nullptr;
	TH1D *  fullAlgJetPT50_100 = nullptr;
	
	TH1D * genJetPT100_150 = nullptr;
	TH1D * recoJetPT100_150 = nullptr;
	TH1D * algJetPT100_150 = nullptr;
	TH1D *  fullAlgJetPT100_150 = nullptr;


	
	d->GetObject("genJetPT50_100", genJetPT50_100);
	d->GetObject("recoJetPT50_100", recoJetPT50_100);
	d->GetObject("algJetPT50_100", algJetPT50_100);
	d->GetObject("fullAlgJetPT50_100", fullAlgJetPT50_100);

	
	genJetPT-50_100>Draw("HIST");
	c->Print("genJetPT50_100.png");
	c->Clear();

	recoJetPT50_100->Draw("HIST");
	c->Print("recoJetPT50_100.png");
	c->Clear();
	
	algJetPT50_100->Draw("HIST");
	c->Print("algJetPT50_100.png");
	c->Clear();
	
	fullAlgJetPT50_100->Draw("HIST");
	c->Print("fullAlgJetPT50_100.png");
	c->Clear();
	
	
	d->GetObject("genJetPT100_150", genJetPT100_150);
	d->GetObject("recoJetPT100_150", recoJetPT100_150);
	d->GetObject("algJetPT100_150", algJetPT100_150);
	d->GetObject("fullAlgJetPT100_150", fullAlgJetPT100_150);

	
	genJetPT100_150 ->Draw("HIST");
	c->Print("genJetPT100_150.png");
	c->Clear();

	recoJetPT100_150 ->Draw("HIST");
	c->Print("recoJetPT100_150.png");
	c->Clear();
	
	algJetPT100_150 ->Draw("HIST");
	c->Print("algJetPT100_150.png");
	c->Clear();
	
	fullAlgJetPT100_150 ->Draw("HIST");
	c->Print("fullAlgJetPT100_150.png");
	c->Clear();
	
	
	

	
	
}
