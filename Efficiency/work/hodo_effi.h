/*
base class for hodoscope efficiency
*/
#ifndef _HODO_EFFI__H_
#define _HODO_EFFI__H_

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

using namespace std;

class hodo_effi
{
    public:
    TFile* file;
    TTree* tree;
    int events;
    
    vector<int>* detectorID = 0;
    vector<int>* elementID_exp = 0;
    vector<int>* elementID_closest = 0;
    vector<double>* x_exp = 0;
    vector<double>* y_exp = 0;
    int nim1;
    int nim2;
    int nim4;
    int fpga5;
    int mask_h1x;
    int mask_h1y;
    int mask_h2x;
    int mask_h2y;
    int mask_h3x;
    int mask_h4x;
    int mask_h4y1;
    int mask_h4y2;
    
    map<string, int> det_id
    {{"H1B", 31},
     {"H1T", 32},
     {"H1L", 33},
     {"H1R", 34},
     {"H2L", 35},
     {"H2R", 36},
     {"H2T", 38},
     {"H2B", 37},
     {"H3B", 39},
     {"H3T", 40},
     {"H4Y1L", 41},
     {"H4Y1R", 42},
     {"H4Y2L", 43},
     {"H4Y2R", 44},
     {"H4B", 45},
     {"H4T", 46},
    };
    
    vector<TString> det_names;
    
    vector<TH1D*> hall;
    vector<TH1D*> hpass;
    vector<TH1D*> hdiff;
    vector<TEfficiency*> effi;
    
    TCanvas* can;
    
    hodo_effi();
    void reg_det(TString det_name);
    void ana();
    void print();
};

#endif /* _HODO_EFFI__H_ */