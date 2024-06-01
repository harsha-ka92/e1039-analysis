#include <iomanip>
#include <fstream>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TFitResult.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <interface_main/SQRun.h>
#include <interface_main/SQEvent.h>
#include <interface_main/SQHitVector.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHIODataNode.h>
#include <phool/getClass.h>
#include <geom_svc/GeomSvc.h>
#include <UtilAna/UtilSQHit.h>
#include "AnaHodoHit.h"
using namespace std;

/// List of detectors that you want to analyze
const vector<string> AnaHodoHit::m_list_det_name = { "H1T", "H1B", "H1L", "H1R", "H2T", "H2B", "H2L", "H2R", "H3T", "H3B", "H4Y1L", "H4Y1R",  "H4Y2L", "H4Y2R", "H4T", "H4B"};

AnaHodoHit::AnaHodoHit(const std::string& name)
  : SubsysReco(name)
  , m_file(0)
  , m_tree(0)
{
  ;
}

int AnaHodoHit::Init(PHCompositeNode* topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaHodoHit::InitRun(PHCompositeNode* topNode)
{
  ///
  /// Input
  ///
  m_evt     = findNode::getClass<SQEvent    >(topNode, "SQEvent");
  m_hit_vec = findNode::getClass<SQHitVector>(topNode, "SQHitVector");
  if (!m_evt || !m_hit_vec) return Fun4AllReturnCodes::ABORTEVENT;

  ///
  /// Output
  ///
  gSystem->mkdir("TDCtime", true);

  m_file = new TFile("result/output.root", "RECREATE");
  m_tree = new TTree("tree", "Created by AnaHodoHit");
  m_tree->Branch("det_name", &b_det_name, "det_name/C");
  m_tree->Branch("det"     , &b_det     ,      "det/I");
  m_tree->Branch("ele"     , &b_ele     ,      "ele/I");
  m_tree->Branch("time"    , &b_time    ,     "time/D");
  
  int nbins; double bin_max; double bin_min; const double DT=12.0/9.0;
  ostringstream oss;
  GeomSvc* geom = GeomSvc::instance();

  for (unsigned int i_det = 0; i_det < m_list_det_name.size(); i_det++) {
    string name = m_list_det_name[i_det];
    int id = geom->getDetectorID(name);
    if (id <= 0) {
      cerr << "!ERROR!  AnaHodoHit::InitRun():  Invalid ID (" << id << ").  Probably the detector name that you specified in 'list_det_name' (" << name << ") is not valid.  Abort." << endl;
      exit(1);
    }
    m_list_det_id.push_back(id);
    int n_ele = geom->getPlaneNElements(id);
    cout << "  " << setw(6) << name << " = " << id << endl;

    std::cout<<"index 2 is"<< name[1] <<std::endl;

    if (name[1]=='1' || name[1]=='2'){bin_min = 616.5*DT; bin_max = 696.5*DT; nbins = 80;}
    /*else if (name[1]=='2'){bin_min = 650.*DT; bin_max = 673.*DT; nbins = 23;}
    else if (name[1]=='3'){bin_min = 706.*DT; bin_max = 729.*DT; nbins = 23;}
    else if (name[1]=='4' && name[2]!='Y'){bin_min = 703.*DT; bin_max = 722.*DT; nbins = 19;}
    else if (name[1]=='4' && name[2]=='Y'){bin_min = 706.*DT; bin_max = 737.*DT; nbins = 31;}*/
    else {bin_min =681.5*DT; bin_max=771.5*DT; nbins = 90;}
    
    oss.str("");
    oss << "TDC_time_of_hits_" << name;
    m_h1_tdc[i_det] = new TH1D(oss.str().c_str(), "", nbins, bin_min, bin_max);
    oss.str("");
    oss << name << ";TDC time(ns);Hit count";
    m_h1_tdc[i_det]->SetTitle(oss.str().c_str());

    oss.str("");
    oss << "h1_nhit_" << name;
    m_h1_nhit[i_det] = new TH1D(oss.str().c_str(), "", 10, -0.5, 9.5);
    oss.str("");
    oss << name << ";N of hits/plane/event;Hit count";
    m_h1_nhit[i_det]->SetTitle(oss.str().c_str());

    oss.str("");
    oss << "Padle_by_padle_tdc_" << name;
    m_h2_tdc[i_det] = new TH2D(oss.str().c_str(), "", n_ele, 0.5, n_ele+0.5, nbins, bin_min, bin_max);
    oss.str("");
    oss << name << ";Element ID;TDC time;Hit count";
    m_h2_tdc[i_det]->SetTitle(oss.str().c_str());

  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaHodoHit::process_event(PHCompositeNode* topNode)
{
  //int spill_id = m_evt->get_spill_id();
  //int event_id = m_evt->get_event_id();

  ///
  ///Event selection: select events from NIM2 for now
  ///
  ///
  ///if (! m_evt->get_trigger(SQEvent::NIM2)) {
  ///  return Fun4AllReturnCodes::EVENT_OK;
  ///}

  ///
  /// Get & fill the hit info
  ///
  for (unsigned int i_det = 0; i_det < m_list_det_name.size(); i_det++) {
    strncpy(b_det_name, m_list_det_name[i_det].c_str(), sizeof(b_det_name));
    b_det = m_list_det_id[i_det];
    shared_ptr<SQHitVector> hv(UtilSQHit::FindHits(m_hit_vec, b_det));
    for (SQHitVector::ConstIter it = hv->begin(); it != hv->end(); it++) {
      b_ele  = (*it)->get_element_id();
      b_time = (*it)->get_tdc_time  ();
      m_tree->Fill();

      m_h1_tdc[i_det]->Fill(b_time);
      m_h2_tdc[i_det]->Fill(b_ele,b_time);
    }
    m_h1_nhit[i_det]->Fill(hv->size());
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaHodoHit::End(PHCompositeNode* topNode)
{
  std::ofstream outfile;
  outfile.open("test.txt", std::ios_base::app);

  //gStyle->SetOptStat(0);
  gStyle->SetOptFit(0111);
  ostringstream oss;
  TCanvas* c1 = new TCanvas("c1", "");
  c1->SetGrid();
  for (unsigned int i_det = 0; i_det < m_list_det_id.size(); i_det++) {
    
    m_h1_tdc[i_det]->Draw();
    
    int max_bin = m_h1_tdc[i_det]->GetMaximumBin();
    double max_bin_x = m_h1_tdc[i_det]->GetXaxis()->GetBinCenter(max_bin);
    
    TFitResultPtr r  = m_h1_tdc[i_det]->Fit("gaus","","", max_bin_x -10, max_bin_x+10);
    c1->Update();
    
    oss.str("");
    oss<< r->Parameter(1)<<","<<r->ParError(1)<<"/"<<r->Parameter(2)<<","<<r->ParError(2) <<"\n";
    r->Write("outfile");

    oss.str("");
    oss << "TDCtime/" << m_h1_tdc[i_det]->GetName() << ".png";
    c1->SaveAs(oss.str().c_str());

    m_h1_nhit[i_det]->Draw();
    oss.str("");
    oss << "TDCtime/" << m_h1_nhit[i_det]->GetName() << ".png";
    c1->SaveAs(oss.str().c_str());

    m_h2_tdc[i_det]->Draw("colz");
    oss.str("");
    oss << "TDCtime/" << m_h2_tdc[i_det]->GetName() << ".png";
    c1->SaveAs(oss.str().c_str());
  }
  delete c1;

  m_file->cd();
  m_file->Write();
  m_file->Close();
  
  return Fun4AllReturnCodes::EVENT_OK;
}
