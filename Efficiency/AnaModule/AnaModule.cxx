#include <iomanip>
#include <TFile.h>
#include <TTree.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>
#include <interface_main/SQEvent.h>
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
  //
  // trigger info.
  //
  nim1 = event->get_trigger(SQEvent::NIM1);
  nim2 = event->get_trigger(SQEvent::NIM2);
  nim3 = event->get_trigger(SQEvent::NIM3);
  nim4 = event->get_trigger(SQEvent::NIM4);
  fpga5 = event->get_trigger(SQEvent::MATRIX5);

  int nTracklets = trackletVec->size();
  for(int i = 0; i < nTracklets; ++i)
  {
    //if(nTracklets < 0) continue;
    Tracklet* tracklet = trackletVec->at(i);
    nHits = tracklet->getNHits();
    chisq = tracklet->getChisq();

    //very loose cuts here
    if(nHits < 10) continue;
    if(chisq > 10.) continue;

    
    /*std::vector<int> vec_detID;
    std::vector<int> vec_closeID;
    std::vector<int> vec_expID;
    std::vector<double> vec_xexp;
    std::vector<double> vec_yexp;*/
   
    for(auto it = detectorIDs.begin(); it != detectorIDs.end(); ++it)
    {
      detectorID = *it;
      double z_exp = p_geomSvc->getPlanePosition(detectorID);
      x_exp = tracklet->getExpPositionX(z_exp);
      y_exp = tracklet->getExpPositionY(z_exp);
      if(!p_geomSvc->isInPlane(detectorID, x_exp, y_exp)) continue;

      elementID_exp = p_geomSvc->getExpElementID(detectorID, tracklet->getExpPositionW(detectorID));
      if(elementID_exp < 1 || elementID_exp > p_geomSvc->getPlaneNElements(detectorID)) continue;

      SQHit* hit = findHit(detectorID, elementID_exp);
      elementID_closest = hit == nullptr ? -1 : hit->get_element_id();

	  vec_detID.push_back(detectorID);
	  vec_expID.push_back(elementID_exp);
	  vec_closeID.push_back(elementID_closest);
      vec_xexp.push_back(x_exp);
      vec_yexp.push_back(y_exp);

      //saveTree->Fill();
    }

    //
    // masking condition for hodoscopes
    //
	for(int i = 0; i < vec_detID.size(); ++i)
	{
      if((vec_detID.at(i) == 31 && vec_expID.at(i) > 0)||(vec_detID.at(i) == 32 && vec_expID.at(i) > 0)) mask_h1x = 1;
	  if((vec_detID.at(i) == 33 && vec_expID.at(i) > 0)||(vec_detID.at(i) == 34 && vec_expID.at(i) > 0)) mask_h1y = 1;
      if((vec_detID.at(i) == 37 && vec_expID.at(i) > 0)||(vec_detID.at(i) == 38 && vec_expID.at(i) > 0)) mask_h2x = 1;
      if((vec_detID.at(i) == 35 && vec_expID.at(i) > 0)||(vec_detID.at(i) == 36 && vec_expID.at(i) > 0)) mask_h2y = 1;
      if((vec_detID.at(i) == 39 && vec_expID.at(i) > 0)||(vec_detID.at(i) == 40 && vec_expID.at(i) > 0)) mask_h3x = 1;
      if((vec_detID.at(i) == 45 && vec_expID.at(i) > 0)||(vec_detID.at(i) == 46 && vec_expID.at(i) > 0)) mask_h4x = 1;
      if((vec_detID.at(i) == 41 && vec_expID.at(i) > 0)||(vec_detID.at(i) == 42 && vec_expID.at(i) > 0)) mask_h4y1 = 1;
      if((vec_detID.at(i) == 43 && vec_expID.at(i) > 0)||(vec_detID.at(i) == 44 && vec_expID.at(i) > 0)) mask_h4y2 = 1;
	}

  }

  saveTree->Fill();
 
  std::cout << "eventID : "<< eventID << " nTracklets : "<< nTracklets << " N detoctors : " << vec_detID.size() << std::endl;

  
  vec_detID.clear();
  vec_closeID.clear();
  vec_expID.clear();
  vec_xexp.clear();
  vec_yexp.clear();

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
  saveTree->Branch("detectorID", &vec_detID);
  saveTree->Branch("elementID_exp", &vec_expID);
  saveTree->Branch("elementID_closest", &vec_closeID);
  saveTree->Branch("x_exp", &vec_xexp);
  saveTree->Branch("y_exp", &vec_yexp);
  saveTree->Branch("nHits", &nHits, "nHits/I");
  saveTree->Branch("chisq", &chisq, "chisq/D");
  // trigger info
  saveTree->Branch("nim1", &nim1, "nim1/I");
  saveTree->Branch("nim2", &nim2, "nim2/I");
  saveTree->Branch("nim3", &nim3, "nim3/I");// only used in noice study
  saveTree->Branch("nim4", &nim4, "nim4/I");
  saveTree->Branch("fpga5", &fpga5, "fpga5/I");
  // hodo mask
  saveTree->Branch("mask_h1x", &mask_h1x, "mask_h1x/I");
  saveTree->Branch("mask_h1y", &mask_h1y, "mask_h1y/I");
  saveTree->Branch("mask_h2x", &mask_h2x, "mask_h2x/I");
  saveTree->Branch("mask_h2y", &mask_h2y, "mask_h2y/I");
  saveTree->Branch("mask_h3x", &mask_h3x, "mask_h3x/I");
  saveTree->Branch("mask_h4x", &mask_h3x, "mask_h4x/I");
  saveTree->Branch("mask_h4y1", &mask_h4y1, "mask_h4y1/I");
  saveTree->Branch("mask_h4y2", &mask_h4y2, "mask_h4y2/I");
}

void AnaModule::registerDetector(TString name)
{
  detectorIDs.insert(p_geomSvc->getDetectorID(name.Data()));
  std::cout << name.Data() << " " << p_geomSvc->getDetectorID(name.Data()) << std::endl;
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
