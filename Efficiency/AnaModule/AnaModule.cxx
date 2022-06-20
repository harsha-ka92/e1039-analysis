#include <iomanip>
#include <TFile.h>
#include <TTree.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>
#include <interface_main/SQHitVector_v1.h>
#include <interface_main/SQTrackVector_v1.h>
#include <interface_main/SQDimuonVector_v1.h>


#include "AnaModule.h"

AnaModule::AnaModule(const std::string& name): SubsysReco(name), p_geomSvc(GeomSvc::instance())
{}

AnaModule::~AnaModule()
{}

int AnaModule::Init(PHCompositeNode* topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaModule::InitRun(PHCompositeNode* topNode)
{
  int ret = GetNodes(topNode);
  if(ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  eventID = 0;
  MakeTree();
  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaModule::process_event(PHCompositeNode* topNode)
{
  int nTracklets = trackletVec->size();
  for(int i = 0; i < nTracklets; ++i)
  {
    Tracklet* tracklet = trackletVec->at(i);
    nHits = tracklet->getNHits();
    chisq = tracklet->getChisq();

    //very loose cuts here
    if(nHits < 9) continue;
    if(chisq > 20.) continue;

    effi_h4(tracklet);

    saveTree->Fill();

    detectorID.clear();
    elementID_exp.clear();
    elementID_closest.clear();
  }

  ++eventID;
  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaModule::End(PHCompositeNode* topNode)
{
  saveFile->cd();
  saveTree->Write();
  saveFile->Close();

  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaModule::GetNodes(PHCompositeNode* topNode)
{
  event = findNode::getClass<SQEvent>(topNode, "SQEvent");
  hitVector   = findNode::getClass<SQHitVector>(topNode, "SQHitVector");
  trackletVec = findNode::getClass<TrackletVector>(topNode, "TrackletVector");
  if(!event || !hitVector || !trackletVec)
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void AnaModule::MakeTree()
{
  saveFile = new TFile(saveName, "RECREATE");

  saveTree = new TTree("save", "Efficiency tree Created by AnaModule");
  saveTree->Branch("eventID", &eventID, "eventID/I");
  saveTree->Branch("detectorID", &detectorID);
  saveTree->Branch("elementID_exp", &elementID_exp);
  saveTree->Branch("elementID_closest", &elementID_closest);
  saveTree->Branch("nHits", &nHits, "nHits/I");
  saveTree->Branch("chisq", &chisq, "chisq/D");
}

SQHit* AnaModule::findHit(int detID, int eleID)
{
  int delta = 999;
  SQHit* hit = nullptr;
  for(SQHitVector::Iter it = hitVector->begin(); it != hitVector->end(); ++it)
  {
    if((*it)->get_detector_id() != detID) continue;
    int delta_temp = abs(eleID - (*it)->get_element_id());
    if(delta > delta_temp)
    {
      delta = delta_temp;
      hit = (*it);
    }
  }

  return hit;
}

int AnaModule::fit_prop(int det_id, Tracklet* tracklet)
{
  std::vector<int> track3 = {19, 21, 22, 51, 52, 53, 54};

  TGraphErrors* gx = new TGraphErrors();
  TGraphErrors* gy = new TGraphErrors();

  int ndet = track3.size();

  for(int i = 0; i < ndet; i++)
  {
    double zz0 = p_geomSvc->getPlanePosition(track3.at(i));
    double xx0 = tracklet->getExpPositionX(zz0);
    double yy0 = tracklet->getExpPositionY(zz0);
    double exx0 = tracklet->getExpPosErrorX(zz0);
    double eyy0 = tracklet->getExpPosErrorY(zz0);

    // set x points
    gx->SetPoint(i, xx0, zz0);
    gx->SetPointError(i, exx0, 0.);

    // set y points
    gy->SetPoint(i, yy0, zz0);
    gy->SetPointError(i, eyy0, 0.);
  }

  // fit functions
  TF1* fx = new TF1("fx", "[0]* x + [1]", 1900., 2400.);
  TF1* fy = new TF1("fy", "[0]* x + [1]", 1900., 2400.);

  gx->Fit("fx");
  gy->Fit("fy");

  double axx = fx->GetParameter(0);
  double cxx = fx->GetParameter(1);

  double ayy = fy->GetParameter(0);
  double cyy = fy->GetParameter(1);

  double zz1 = p_geomSvc->getPlanePosition(det_id);
  double xx1 = axx* zz1 + cxx;
  double yy1 = ayy* zz1 + cyy;

  if(!p_geomSvc->isInPlane(det_id, xx1, yy1)) continue;

  double pos = p_geomSvc->getCostheta(det_id)*xx1 + p_geomSvc->getSintheta(det_id)*yy1;

  return p_geomSvc->getExpElementID(det_id, pos);
}

void AnaModule::effi_h4(Tracklet* tracklet)
{
  // only NIM4 events are considered
  if(!event->get_trigger(SQEvent::NIM4)) continue;

  std::vector<int> hodo3 = {41, 42, 43, 44, 45, 46};

  int nhodo = hodo3.size();

  for(int i = 0; i < nhodo; i++)
  {
    int det_id = hodo3.at(i);
    int exp_id = fit_prop(det_id, tracklet);

    SQHit* hit = findHit(det_id, exp_id);
    int close_id = hit == nullptr ? -1 : hit->get_element_id();

    detectorID.push_back(det_id);
    elementID_exp.push_back(exp_id);
    elementID_closest.push_back(close_id);
  }

}
