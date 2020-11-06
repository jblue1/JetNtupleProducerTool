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
	
	


	
	d->GetObject("genJetPT", genJetPT);
	d->GetObject("recoJetPT", recoJetPT);
	d->GetObject("algJetPT", algJetPT);
	d->GetObject("fullAlgJetPT", fullAlgJetPT);

	recoJetPT->SetLineColor(kRed);
	recoJetPT->SetStats(0);
	recoJetPT->Draw("SAME HIST");

	genJetPT->SetLineColor(kBlue);
	genJetPT->SetStats(0);
	genJetPT->Draw("HIST");
	

	
	algJetPT->SetLineColor(kGreen);
	algJetPT->SetStats(0);
	algJetPT->Draw("SAME HIST");
	
	fullAlgJetPT->SetLineColor(kBlack);
	fullAlgJetPT->SetStats(0);
	fullAlgJetPT->Draw("SAME HIST");
	
	TLegend * legendJetPT = new TLegend();
	legendJetPT->AddEntry(genJetPT, "gen", "L");
	legendJetPT->AddEntry(recoJetPT, "reco", "L");
	legendJetPT->AddEntry(algJetPT, "algorithm", "L");
	legendJetPT->AddEntry(fullAlgJetPT, "full algorithm", "L");
	legendJetPT->Draw();
	
	c->Print("JetPT.png");
	c->Clear();
	
	
	TH1D * genJetPT0_50 = nullptr;
	TH1D * recoJetPT0_50 = nullptr;
	TH1D * algJetPT0_50 = nullptr;
	TH1D *  fullAlgJetPT0_50 = nullptr;
	
	recoJetPT0_50->SetLineColor(kRed);
	recoJetPT0_50->SetStats(0);
	recoJetPT0_50->Draw("SAME HIST");
	
	genJetPT0_50->SetLineColor(kBlue);
	genJetPT0_50->SetStats(0);
	genJetPT0_50->Draw("HIST");
	
	
	
	algJetPT0_50->SetLineColor(kGreen);
	algJetPT0_50->SetStats(0);
	algJetPT0_50->Draw("SAME HIST");
	
	fullAlgJetPT0_50->SetLineColor(kBlack);
	fullAlgJetPT0_50->SetStats(0);
	fullAlgJetPT0_50->Draw("SAME HIST");
	
	TLegend * legendJetPT0_50 = new TLegend();
	legendJetPT0_50->AddEntry(genJetPT0_50, "gen", "L");
	legendJetPT0_50->AddEntry(recoJetPT0_50, "reco", "L");
	legendJetPT0_50->AddEntry(algJetPT0_50, "algorithm", "L");
	legendJetPT0_50->AddEntry(fullAlgJetPT0_50, "full algorithm", "L");
	legendJetPT0_50->Draw();
	
	c->Print("JetPT0_50.png");
	c->Clear();
	
	
	d->GetObject("genJetPT50_100", genJetPT50_100);
	d->GetObject("recoJetPT50_100", recoJetPT50_100);
	d->GetObject("algJetPT50_100", algJetPT50_100);
	d->GetObject("fullAlgJetPT50_100", fullAlgJetPT50_100);

	
	recoJetPT50_100->SetLineColor(kRed);
	recoJetPT50_100->SetStats(0);
	recoJetPT50_100->Draw("SAME HIST");	


	genJetPT50_100->SetLineColor(kBlue);
	genJetPT50_100->SetStats(0);
	genJetPT50_100->Draw("HIST");
	

	
	algJetPT50_100->SetLineColor(kGreen);
	algJetPT50_100->SetStats(0);
	algJetPT50_100->Draw("SAME HIST");
	
	fullAlgJetPT50_100->SetLineColor(kBlack);
	fullAlgJetPT50_100->SetStats(0);
	fullAlgJetPT50_100->Draw("SAME HIST");
	
	TLegend * legendJetPT50_100 = new TLegend();
	legendJetPT50_100->AddEntry(genJetPT50_100, "gen", "L");
	legendJetPT50_100->AddEntry(recoJetPT50_100, "reco", "L");
	legendJetPT50_100->AddEntry(algJetPT50_100, "algorithm", "L");
	legendJetPT50_100->AddEntry(fullAlgJetPT50_100, "full algorithm", "L");
	legendJetPT50_100->Draw();
	
	c->Print("JetPT50_100.png");
	c->Clear();
	
	
	TH1D * genJetPT100_150 = nullptr;
	TH1D * recoJetPT100_150 = nullptr;
	TH1D * algJetPT100_150 = nullptr;
	TH1D *  fullAlgJetPT100_150 = nullptr;
	
	d->GetObject("genJetPT100_150", genJetPT100_150);
	d->GetObject("recoJetPT100_150", recoJetPT100_150);
	d->GetObject("algJetPT100_150", algJetPT100_150);
	d->GetObject("fullAlgJetPT100_150", fullAlgJetPT100_150);

	
	recoJetPT100_150->SetLineColor(kRed);
	recoJetPT100_150->SetStats(0);
	recoJetPT100_150->Draw("SAME HIST");
	
	genJetPT100_150->SetLineColor(kBlue);
	genJetPT100_150->SetStats(0);
	genJetPT100_150->Draw("HIST");
	
	
	
	algJetPT100_150->SetLineColor(kGreen);
	algJetPT100_150->SetStats(0);
	algJetPT100_150->Draw("SAME HIST");
	
	fullAlgJetPT100_150->SetLineColor(kBlack);
	fullAlgJetPT100_150->SetStats(0);
	fullAlgJetPT100_150->Draw("SAME HIST");
	
	TLegend * legendJetPT100_150 = new TLegend();
	legendJetPT100_150->AddEntry(genJetPT100_150, "gen", "L");
	legendJetPT100_150->AddEntry(recoJetPT100_150, "reco", "L");
	legendJetPT100_150->AddEntry(algJetPT100_150, "algorithm", "L");
	legendJetPT100_150->AddEntry(fullAlgJetPT100_150, "full algorithm", "L");
	legendJetPT100_150->Draw();
	
	c->Print("JetPT100_150.png");
	c->Clear();
	
	d->GetObject("genJetPT100_150", genJetPT100_150);
	d->GetObject("recoJetPT100_150", recoJetPT100_150);
	d->GetObject("algJetPT100_150", algJetPT100_150);
	d->GetObject("fullAlgJetPT100_150", fullAlgJetPT100_150);

	
	TH1D * genJetPT150_200 = nullptr;
	TH1D * recoJetPT150_200 = nullptr;
	TH1D * algJetPT150_200 = nullptr;
	TH1D *  fullAlgJetPT150_200 = nullptr;
	
	recoJetPT150_200->SetLineColor(kRed);
	recoJetPT150_200->SetStats(0);
	recoJetPT150_200->Draw("SAME HIST");
	
	genJetPT150_200->SetLineColor(kBlue);
	genJetPT150_200->SetStats(0);
	genJetPT150_200->Draw("HIST");
	
	
	
	algJetPT150_200->SetLineColor(kGreen);
	algJetPT150_200->SetStats(0);
	algJetPT150_200->Draw("SAME HIST");
	
	fullAlgJetPT150_200->SetLineColor(kBlack);
	fullAlgJetPT150_200->SetStats(0);
	fullAlgJetPT150_200->Draw("SAME HIST");
	
	TLegend * legendJetPT150_200 = new TLegend();
	legendJetPT150_200->AddEntry(genJetPT150_200, "gen", "L");
	legendJetPT150_200->AddEntry(recoJetPT150_200, "reco", "L");
	legendJetPT150_200->AddEntry(algJetPT150_200, "algorithm", "L");
	legendJetPT150_200->AddEntry(fullAlgJetPT150_200, "full algorithm", "L");
	legendJetPT150_200->Draw();
	
	c->Print("JetPT150_200.png");
	c->Clear();
	
	

	
	
}
