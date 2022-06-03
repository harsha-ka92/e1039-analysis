/*
base class for hodoscope debugging
*/
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TEfficiency.h>
#include <TGraphErrors.h>
#include <TString.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <map>
#include <vector>
#include <iostream>

#include "hodo_debug.h"

using namespace std;

hodo_debug::hodo_debug()
{
    file = TFile::Open("ana.root", "READ");
    tree = (TTree*)file->Get("save");
    events = tree->GetEntries();
    
    tree->SetBranchAddress("detectorID", &detectorID);
    tree->SetBranchAddress("elementID_exp", &elementID_exp);
    tree->SetBranchAddress("elementID_closest", &elementID_closest);
    tree->SetBranchAddress("x_exp", &x_exp);
    tree->SetBranchAddress("y_exp", &y_exp);
    tree->SetBranchAddress("nim1", &nim1);
    tree->SetBranchAddress("nim2", &nim2);
    tree->SetBranchAddress("nim4", &nim4);
    tree->SetBranchAddress("fpga5", &fpga5);
    tree->SetBranchAddress("mask_h1x", &mask_h1x);
    tree->SetBranchAddress("mask_h1y", &mask_h1y);
    tree->SetBranchAddress("mask_h2x", &mask_h2x);
    tree->SetBranchAddress("mask_h2y", &mask_h2y);
    tree->SetBranchAddress("mask_h3x", &mask_h3x);
    tree->SetBranchAddress("mask_h4x", &mask_h4x);
    tree->SetBranchAddress("mask_h4y1", &mask_h4y1);
    tree->SetBranchAddress("mask_h4y2", &mask_h4y2);
}

void hodo_debug::reg_det(TString det_name)
{
    det_names.push_back(det_name);
}

void hodo_debug::ana()
{
    TH1D* h1;
    TH1D* h2;
    TH1D* h3;
    TEfficiency* e1;
    
    for(unsigned int i = 0; i < det_names.size(); i++)
    {
        auto det_name = det_names.at(i);
        
        if(det_name.EqualTo("H4Bu"))
        {
            TString aname = Form("hall_%s", det_name.Data());
            TString title = Form("; ele_id [det_name = %s]; counts", det_name.Data());
            h1 = new TH1D(aname.Data(), title.Data(), 16, 0.5, 16.5);
            
            TString pname = Form("hpass_%s", det_name.Data());
            h2 = new TH1D(pname.Data(), title.Data(), 16, 0.5, 16.5);
            
            TString dname = Form("hdiff_%s", det_name.Data());
            h3 = new TH1D(dname.Data(), title.Data(), 20, -10.0, 10.0);
            
            for(int j = 0; j < events; j++)
            {
                tree->GetEntry(j);
                for(unsigned int k = 0; k < detectorID->size(); k++)
                {
                    if(detectorID->at(k) == det_id[det_name.Data()])
                    {
                        if(nim4 > 0 && mask_h4y1 > 0 && y_exp->at(k) > -90.)
                        {
                            h1->Fill(elementID_exp->at(k));
                            if(elementID_closest->at(k) > 0)
                            {
                                h2->Fill(elementID_exp->at(k));
                                h3->Fill(elementID_exp->at(k) - elementID_closest->at(k));
                            }
                        }
                    }
                }
            }// end of tree
            
            e1 = new TEfficiency(*h2, *h1);
            
            TString ename = Form("effi_%s", det_name.Data());
            TString etitle = Form("; elemet_id [det_name = %s]; efficiency", det_name.Data());
            e1->SetName(ename.Data());
            e1->SetTitle(etitle.Data());
        }
        //
        //
        if(det_name.EqualTo("H4Y1Rr"))
        {
            TString aname = Form("hall_%s", det_name.Data());
            TString title = Form("; ele_id [det_name = %s]; counts", det_name.Data());
            h1 = new TH1D(aname.Data(), title.Data(), 16, 0.5, 16.5);
            
            TString pname = Form("hpass_%s", det_name.Data());
            h2 = new TH1D(pname.Data(), title.Data(), 16, 0.5, 16.5);
            
            TString dname = Form("hdiff_%s", det_name.Data());
            h3 = new TH1D(dname.Data(), title.Data(), 20, -10.0, 10.0);
            
            for(int j = 0; j < events; j++)
            {
                tree->GetEntry(j);
                for(unsigned int k = 0; k < detectorID->size(); k++)
                {
                    if(detectorID->at(k) == det_id[det_name.Data()])
                    {
                        if(nim4 > 0 && mask_h4x > 0 && x_exp->at(k) > -65.)
                        {
                            h1->Fill(elementID_exp->at(k));
                            if(elementID_closest->at(k) > 0)
                            {
                                h2->Fill(elementID_exp->at(k));
                                h3->Fill(elementID_exp->at(k) - elementID_closest->at(k));
                            }
                        }
                    }
                }
            }// end of tree
            
            e1 = new TEfficiency(*h2, *h1);
            
            TString ename = Form("effi_%s", det_name.Data());
            TString etitle = Form("; elemet_id [det_name = %s]; efficiency", det_name.Data());
            e1->SetName(ename.Data());
            e1->SetTitle(etitle.Data());
        }
        //
        //
        if(det_name.EqualTo("H4Y1Lr"))
        {
            TString aname = Form("hall_%s", det_name.Data());
            TString title = Form("; ele_id [det_name = %s]; counts", det_name.Data());
            h1 = new TH1D(aname.Data(), title.Data(), 16, 0.5, 16.5);
            
            TString pname = Form("hpass_%s", det_name.Data());
            h2 = new TH1D(pname.Data(), title.Data(), 16, 0.5, 16.5);
            
            TString dname = Form("hdiff_%s", det_name.Data());
            h3 = new TH1D(dname.Data(), title.Data(), 20, -10.0, 10.0);
            
            for(int j = 0; j < events; j++)
            {
                tree->GetEntry(j);
                for(unsigned int k = 0; k < detectorID->size(); k++)
                {
                    if(detectorID->at(k) == det_id[det_name.Data()])
                    {
                        if(nim4 > 0 && mask_h4x > 0 && x_exp->at(k) > 65.)
                        {
                            h1->Fill(elementID_exp->at(k));
                            if(elementID_closest->at(k) > 0)
                            {
                                h2->Fill(elementID_exp->at(k));
                                h3->Fill(elementID_exp->at(k) - elementID_closest->at(k));
                            }
                        }
                    }
                }
            }// end of tree
            
            e1 = new TEfficiency(*h2, *h1);
            
            TString ename = Form("effi_%s", det_name.Data());
            TString etitle = Form("; elemet_id [det_name = %s]; efficiency", det_name.Data());
            e1->SetName(ename.Data());
            e1->SetTitle(etitle.Data());
        }
        //
        //
        if(det_name.EqualTo("H4Y2Rr"))
        {
            TString aname = Form("hall_%s", det_name.Data());
            TString title = Form("; ele_id [det_name = %s]; counts", det_name.Data());
            h1 = new TH1D(aname.Data(), title.Data(), 16, 0.5, 16.5);
            
            TString pname = Form("hpass_%s", det_name.Data());
            h2 = new TH1D(pname.Data(), title.Data(), 16, 0.5, 16.5);
            
            TString dname = Form("hdiff_%s", det_name.Data());
            h3 = new TH1D(dname.Data(), title.Data(), 20, -10.0, 10.0);
            
            for(int j = 0; j < events; j++)
            {
                tree->GetEntry(j);
                for(unsigned int k = 0; k < detectorID->size(); k++)
                {
                    if(detectorID->at(k) == det_id[det_name.Data()])
                    {
                        if(nim4 > 0 && mask_h4x > 0 && x_exp->at(k) > -65.)
                        {
                            h1->Fill(elementID_exp->at(k));
                            if(elementID_closest->at(k) > 0)
                            {
                                h2->Fill(elementID_exp->at(k));
                                h3->Fill(elementID_exp->at(k) - elementID_closest->at(k));
                            }
                        }
                    }
                }
            }// end of tree
            
            e1 = new TEfficiency(*h2, *h1);
            
            TString ename = Form("effi_%s", det_name.Data());
            TString etitle = Form("; elemet_id [det_name = %s]; efficiency", det_name.Data());
            e1->SetName(ename.Data());
            e1->SetTitle(etitle.Data());
        }
        //
        //
        hall.push_back(h1);
        hpass.push_back(h2);
        hdiff.push_back(h3);
        effi.push_back(e1);
    }// end of detectors
}

void hodo_debug::print()
{
    //gROOT->SetStyle("ATLAS");
    TGraphAsymmErrors* G;
    for(unsigned int i = 0; i < det_names.size(); i++)
    {
        can = new TCanvas(effi.at(i)->GetName(), "", 1000, 500);
        can->Divide(2);
        can->cd(1);
        
        effi.at(i)->SetMarkerColor(2);
        effi.at(i)->SetStatisticOption(TEfficiency::kBBayesian);
        effi.at(i)->SetConfidenceLevel(0.68);
        
        effi.at(i)->Draw("APE1");
        can->Update();
        
        G = effi.at(i)->GetPaintedGraph();
        G->SetMinimum(0.0);
        G->SetMaximum(1.05);
        can->Update();
        
        can->cd(2);
        hdiff.at(i)->SetFillColor(kAzure+10);
        hdiff.at(i)->Draw();
        
        can->Draw();
    }
}