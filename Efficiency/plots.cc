#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TEfficiency.h>
#include <TCanvas.h>

void plots()
{
	TFile* file = TFile::Open("ana.root", "READ");
	TTree* save = (TTree*)file->Get("save");
	int n = save->GetEntries();

	int detectorID, elementID_exp, elementID_closest;

	save->SetBranchAddress("detectorID", &detectorID);
	save->SetBranchAddress("elementID_exp", &elementID_exp);
	save->SetBranchAddress("elementID_closest", &elementID_closest);

	TH1D* hist1 = new TH1D("hist1", "; element id; counts", 18, 0.0, 18.0);
	TH1D* hist2 = new TH1D("hist2", "; element id; counts", 18, 0.0, 18.0);

	TH1D* hist3 = new TH1D("hist3", "; paddel diff; counts", 20, -10.0, 10.0);

	for(int i = 0; i < n; i++)
	{
		save->GetEntry(i);
		if(elementID_closest > 0)
		{
			if(detectorID == 41)
			{
				hist1->Fill(elementID_exp);
				if(abs(elementID_exp - elementID_closest) < 5){hist2->Fill(elementID_exp);}
				hist3->Fill(elementID_exp - elementID_closest);
			}
		}
	}

	TEfficiency* effi = new TEfficiency(*hist2, *hist1);

	TCanvas* can = new TCanvas();
	effi->Draw("APE1");
	can->Draw();

}