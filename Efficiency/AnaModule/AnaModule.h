#ifndef _ANA_Module__H_
#define _ANA_Module__H_

#include <map>
#include <set>
#include <fun4all/SubsysReco.h>
#include <TString.h>
#include <TVector3.h>
#include <interface_main/SQEvent.h>
#include <ktracker/SRecEvent.h>
#include <ktracker/FastTracklet.h>
#include <geom_svc/GeomSvc.h>
#include <interface_main/SQHit_v1.h>


class TFile;
class TTree;
class SQHitVector;
class SQTrackVector;
class SQDimuonVector;

class AnaModule: public SubsysReco 
{
public:
  AnaModule(const std::string& name = "AnaModule");
  virtual ~AnaModule();

  int Init(PHCompositeNode* topNode);
  int InitRun(PHCompositeNode* topNode);
  int process_event(PHCompositeNode* topNode);
  int End(PHCompositeNode* topNode);

  void set_output_filename(const TString& n) { saveName = n; }
  void registerDetector(TString name);

private:
  int GetNodes(PHCompositeNode* topNode);
  void MakeTree();

  SQHit* findHit(int detectorID, int elementID);
  std::set<int> detectorIDs;

  GeomSvc* p_geomSvc;

  // Input
  SQHitVector*    hitVector;
  TrackletVector* trackletVec;
  SQEvent* event;

  // Output
  TString saveName;
  TFile*  saveFile;
  TTree*  saveTree;

  int eventID;
  int detectorID;
  int elementID_exp;
  int elementID_closest;
  double x_exp;
  double y_exp;
  int nHits;
  double chisq;

  // event info.
  std::vector<int> vec_detID;
  std::vector<int> vec_closeID;
  std::vector<int> vec_expID;
  std::vector<double> vec_xexp;
  std::vector<double> vec_yexp;

  // trigger info.
  int nim1;
  int nim2;
  int nim3;
  int nim4;
  int fpga5;
  // hodo masking
  int mask_h1x;
  int mask_h1y;
  int mask_h2x;
  int mask_h2y;
  int mask_h3x;
  int mask_h4x;
  int mask_h4y1;
  int mask_h4y2;
};

#endif
