/*
base class for hodoscope efficiency
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
#include <string>
#include <fstream>


#include "hodo_effi.h"

using namespace std;

hodo_effi::hodo_effi()
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

void hodo_effi::reg_det(TString det_name)
{
    det_names.push_back(det_name);
}

void hodo_effi::ana()
{
    TH1D* h1;
    TH1D* h2;
    TH1D* h3;
    TEfficiency* e1;
    
    for(unsigned int i = 0; i < det_names.size(); i++)
    {
        auto det_name = det_names.at(i);
        //
        //
        if(det_name.EqualTo("H1B"))
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
                        if((nim2 > 0 || fpga5 > 0) && mask_h1y > 0)
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
        if(det_name.EqualTo("H1T"))
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
                        if((nim2 > 0 || fpga5 > 0) && mask_h1y > 0)
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
        if(det_name.EqualTo("H1L"))
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
                        if((nim2 > 0 || fpga5 > 0) && mask_h1x > 0)
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
        if(det_name.EqualTo("H1R"))
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
                        if((nim2 > 0 || fpga5 > 0) && mask_h1x > 0)
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
        if(det_name.EqualTo("H2T"))
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
                        if((nim2 > 0 || fpga5 > 0) && mask_h2y > 0)
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
        if(det_name.EqualTo("H2B"))
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
                        if((nim2 > 0 || fpga5 > 0) && mask_h2y > 0)
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
        if(det_name.EqualTo("H2L"))
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
                        if((nim2 > 0 || fpga5 > 0) && mask_h2x > 0)
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
        if(det_name.EqualTo("H2R"))
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
                        if((nim2 > 0 || fpga5 > 0) && mask_h2x > 0)
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
        if(det_name.EqualTo("H3T"))
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
                        if(nim4 > 0 && mask_h4y1 > 0)
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
        if(det_name.EqualTo("H3B"))
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
                        if(nim4 > 0 && mask_h4y1 > 0)
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
        if(det_name.EqualTo("H4B"))
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
                        if(nim4 > 0 && mask_h4y1 > 0)
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
        if(det_name.EqualTo("H4T"))
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
                        if(nim4 > 0 && mask_h4y1 > 0)
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
        if(det_name.EqualTo("H4Y1L"))
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
                        if(nim4 > 0 && mask_h4x > 0)
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
        if(det_name.EqualTo("H4Y1R"))
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
                        if(nim4 > 0 && mask_h4x > 0)
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
        if(det_name.EqualTo("H4Y2L"))
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
                        if(nim4 > 0 && mask_h4x > 0)
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
        if(det_name.EqualTo("H4Y2R"))
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
                        if(nim4 > 0 && mask_h4x > 0)
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
        
        hall.push_back(h1);
        hpass.push_back(h2);
        hdiff.push_back(h3);
        effi.push_back(e1);
    }// end of detectors
}

void hodo_effi::print()
{
    //gROOT->SetStyle("ATLAS");
    //Get the Plane ID and PMT number 
    TGraphAsymmErrors* pmt_eff = new TGraphAsymmErrors();
    char x;
    int y;
    char str[80];
    cout << "Enter detector plane (Available IDs: ):\n";
    cout << "Enter the PMT number:\n";
    cin >> y;
    cin >> x;
    
    // get the voltages from txt file and store in an array (are vectors supported in TGraphAsymmErrors?).
    
    double hv_step[5]; 
    short entry=0; 
    double value; 
    ifstream myfile ("RunVoltages.txt"); 
        if (myfile.is_open()) 
        {
            while (! myfile.eof() ) 
            {
                getline (myfile,value); 
                hv_step[entry] = value;
                entry++;
            }
            myfile.close(); 
        }
        else cout << "can't open the file"; 
    
    //Plot efficeincies of the requested paddle
    for (unsigned int runs=0;runs<5; runs++)
    {
        
        for(unsigned int i = 0; i < det_names.size(); i++)
         {  
            if (effi.at(i)->GetName()==x)
            {
             can = new TCanvas("paddle efficiency" , "", 1000, 500);
             pmt_eff->SetPoint(runs, hvStep[runs], eff.at(i)->GetEfficiency(y-1));
			 pmt_eff->SetPointError(runs, 0., 0., eff->GetEfficiencyErrorUp(y-1), eff->GetEfficiencyErrorLow(y-1));
			
			
		    TCanvas* c1 = new TCanvas("c1", "", 1000, 500);		
   		    c1->SetFillColor(42);
   		    c1->SetGrid();
   		    c1->GetFrame()->SetFillColor(21);
   		    c1->GetFrame()->SetBorderSize(12);

		    pmt_eff->SetTitle("Efficiency variation with voltage");
   		    pmt_eff->SetMarkerColor(4);
   		    pmt_eff->SetMarkerStyle(21);
   		    pmt_eff->Draw("SAME");
            }
        }
    }
}
