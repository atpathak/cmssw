// -*- C++ -*-
//
// Package:    SiPixelPhase1bumpbonding
// Class:      SiPixelPhase1bumpbonding
//

// Author: Atanu Pathak

// C++ stuff
#include <iostream>
#include <string>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>

// CMSSW stuff
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

// DQM Stuff
#include "DQMServices/Core/interface/DQMStore.h"

// Input data stuff
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"

// PixelDQM Framework
#include "DQM/SiPixelPhase1Common/interface/SiPixelPhase1Base.h"

using namespace std;
using namespace edm;


namespace {

  class SiPixelPhase1bumpbonding : public DQMEDHarvester {
  public:
    explicit SiPixelPhase1bumpbonding(const edm::ParameterSet& conf);
    ~SiPixelPhase1bumpbonding() override;


  protected:
    // BeginJob
    //  void beginJob(void) override;
    
    // BeginRun 
    void beginRun(edm::Run const& run, edm::EventSetup const& eSetup) override;

    // EndJob
    void dqmEndJob(DQMStore::IBooker& iBooker, DQMStore::IGetter& iGetter) override;
    
  private:
    std::string topFolderName_;
    int minHits_;
    edm::ParameterSet conf_;
    
    std::map<std::string, MonitorElement*> residuals_;
    std::map<std::string, MonitorElement*> DRnR_;
    
    //Book Monitoring Elements
    void bookMEs(DQMStore::IBooker& iBooker);
  
    //Fill Monitoring Elements
    void fillMEs(DQMStore::IBooker& iBooker, DQMStore::IGetter& iGetter);

  }; //class

  //---------------------------
  SiPixelPhase1bumpbonding::SiPixelPhase1bumpbonding(const edm::ParameterSet& iConfig)
    : DQMEDHarvester(iConfig), conf_(iConfig) {
  LogInfo("PixelDQM") << "SiPixelPhase1bumpbonding::SiPixelPhase1bumpbonding: Got DQM BackEnd interface" << endl;
  topFolderName_ = conf_.getParameter<std::string>("TopFolderName");
  minHits_ = conf_.getParameter<int>("MinHits");
  }

  SiPixelPhase1bumpbonding::~SiPixelPhase1bumpbonding() {
    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)
    LogInfo("PixelDQM") << "SiPixelPhase1bumpbonding::~SiPixelPhase1bumpbonding: Destructor" << endl;
  }
  
  void SiPixelPhase1bumpbonding::beginRun(edm::Run const& run, edm::EventSetup const& eSetup) {}
  
  
  void SiPixelPhase1bumpbonding::dqmEndJob(DQMStore::IBooker& iBooker, DQMStore::IGetter& iGetter) {
    bookMEs(iBooker);
    fillMEs(iBooker, iGetter);
  }
  //-------------------------


  //------------------------------------------------------------------
  // Used to book the MEs
  //------------------------------------------------------------------

  void SiPixelPhase1bumpbonding::bookMEs(DQMStore::IBooker& iBooker) {

    iBooker.cd();
    iBooker.setCurrentFolder(topFolderName_+"/PXForward/");

    residuals_["digi_occupancy_per_col_per_row_FPix_BmI_D1_BLD1_PNL1_RNG1_again"] =
      iBooker.book2D("digi_occupancy_per_col_per_row_FPix_BmI_D1_BLD1_PNL1_RNG1_again",
                       "FPix_BmI_D1_BLD1_PNL1_RNG1",416,-0.5,415.5, 160,-0.5,159.5);

    //Reset the iBooker
    iBooker.setCurrentFolder("PixelPhase1/");


      
  }

  //------------------------------------------------------------------
  // Fill the MEs
  //------------------------------------------------------------------

  void SiPixelPhase1bumpbonding::fillMEs(DQMStore::IBooker& iBooker, DQMStore::IGetter& iGetter) {

    MonitorElement* digi_1 =
      iGetter.get("PixelPhase1/Phase1_MechanicalView/PXForward/HalfCylinder_mI/PXRing_1/PXDisk_-1/SignedBlade_1/digi_occupancy_per_col_per_row_FPix_BmI_D1_BLD1_PNL1_RNG1");

    for (int x = 0; x < digi_1->getNbinsX(); x++) {
      for (int y = 0; y < digi_1->getNbinsY(); y++) {
	std::cout << "SiPixelPhase1bumpbonding::~SiPixelPhase1bumpbonding: x: "<< x <<" y: "<< y << " Content: "<< digi_1->getBinContent(x,y) << std::endl;
	residuals_["digi_occupancy_per_col_per_row_FPix_BmI_D1_BLD1_PNL1_RNG1_again"]->setBinContent(x,y,digi_1->getBinContent(x,y));
      }
    }
  }
  
}  //namespace

DEFINE_FWK_MODULE(SiPixelPhase1bumpbonding);
