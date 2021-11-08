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

const static int COL = 416;
const static int ROW = 160;
const static int ROCCOL = 80;
const static int ROCROW = 52;
const static int ROC = 16;

namespace {

  struct module_bin {
    std::string name;
    bool pixel[ROW][COL];
    bool roc[ROC][ROCROW][ROCCOL];
  };
  
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
    void dqmEndLuminosityBlock(DQMStore::IBooker& iBooker,
			       DQMStore::IGetter& iGetter,
			       edm::LuminosityBlock const& lumiSeg,
			       edm::EventSetup const& c) override;

  private:
    std::string topFolderName_;
    int minHits_;
    edm::ParameterSet conf_;
    bool firstLumi;
    std::map<std::string, MonitorElement*> residuals_;
    std::map<std::string, MonitorElement*> DRnR_;
    map<std::string, module_bin> modules;
    map<std::string, module_bin> modulesSurrounded;
    
    //Book Monitoring Elements
    void bookMEs(DQMStore::IBooker& iBooker);

    //
    std::vector<uint8_t> get_roc_col_row(uint16_t x, uint16_t y);
    std::vector<int> get_roc_col_row_py(int x, int y);
    bool is_safe(module_bin& M, int row, int col);
    bool is_surrounded(module_bin& M, int row, int col);
    module_bin sourroundedRegion(module_bin& M);
    module_bin convert_to_bool(MonitorElement* h);
  
    //Fill Monitoring Elements
    void fillMEs(DQMStore::IBooker& iBooker, DQMStore::IGetter& iGetter);

  }; //class

  //---------------------------
  SiPixelPhase1bumpbonding::SiPixelPhase1bumpbonding(const edm::ParameterSet& iConfig)
    : DQMEDHarvester(iConfig), conf_(iConfig), firstLumi(true) {
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
  
  
  void SiPixelPhase1bumpbonding::dqmEndLuminosityBlock(DQMStore::IBooker& iBooker,
						       DQMStore::IGetter& iGetter,
						       const edm::LuminosityBlock& lumiSeg,
						       edm::EventSetup const& c) {
    if(firstLumi){
      bookMEs(iBooker);
      firstLumi=false;
    }
    fillMEs(iBooker, iGetter);
  }

  void SiPixelPhase1bumpbonding::dqmEndJob(DQMStore::IBooker& iBooker, DQMStore::IGetter& iGetter) {

  }

  //-------------------------
  //-----------Custom functions--------------------------------------
  //
  //-----------------------------------------------------------------

  std::vector<uint8_t> SiPixelPhase1bumpbonding::get_roc_col_row(uint16_t x, uint16_t y) {
    uint8_t row=99;
    uint8_t col=99;
    uint8_t roc=99;
    if ( y == 0){
      //lower row of ROCS in DQM plot
      row = y;
      roc = 7 - int(x/52);
      col = x%52;
    }else if ( y > 0 and y <= 79 ){
      //lower row of ROCS in DQM plot
      row = y;
      roc = 7 - int(x/52);
      col = x%52;
    }else if ( y >= 80 and y <= 160 ){
      // upper row of ROCS in DQM plot
      row = 159 - y;
      roc = 8 + (int(x/52));
      col = x%52;
    }
    return {roc,col,row};
  }
  //
  std::vector<int> SiPixelPhase1bumpbonding::get_roc_col_row_py(int x, int y) {
    y=y+1;
    int row = 99;
    int col = 99;
    int roc = 99;
    if ( y>=1 and y<=80 ){
      //lower row of ROCS in DQM plot
      row = y-1;
      if(x%52 == 0){
	roc = 8 - int(x/52);
      }else{
	roc = 7 - (int(x/52));
      }
      col = int((416 - x)%52);  
    } else if (y>=81 and y<=160){
      //upper row of ROCS in DQM plot
      row = 160 - y;
      if(x%52 == 0){
	roc = 7 + int(x/52);
      }else{
	roc = 8 + (int(x/52));
      }
      col = int((x-1)%52);
    }
    return {roc,col,row};
  }
  //
  bool SiPixelPhase1bumpbonding::is_safe(module_bin& M, int row, int col){
    int rowM = sizeof(M.pixel) / sizeof(M.pixel[0]);
    int colM = sizeof(M.pixel[0]);
    if(row<0 or col<0){
      return false ;
    }else if(row+1 >= rowM or col+1 >= colM){
      return false;
    }else{
      return true;
    }
  }
  //
  bool SiPixelPhase1bumpbonding::is_surrounded(module_bin& M, int row, int col){
    const static int rowNbr[8] = {-1, -1, -1, 0, 0, 1, 1, 1};
    const static int colNbr[8] = {-1, 0, 1, -1, 1, -1, 0, 1};
    int i_surrounded=0;
    for (int k = 0; k<8; k++){
      if (is_safe(M, row + rowNbr[k], col + colNbr[k])==true && M.pixel[row + rowNbr[k]][col + colNbr[k]]==false){
	i_surrounded = i_surrounded + 1;
      }else if (is_safe(M, row + rowNbr[k], col + colNbr[k]==false)){
	i_surrounded = i_surrounded + 1;
      }else {
	return false;
      }
    }
    //std::cout << " i_surrounded: " << i_surrounded << std::endl;
    if(i_surrounded==8){
      return true;
    }
    return false;
  }
  //
  module_bin SiPixelPhase1bumpbonding::sourroundedRegion(module_bin& M) {
    module_bin surrounded;
    for (int i = 0; i < ROW; ++i) {
      for (int j = 0; j < COL; ++j) {
	if(M.pixel[i][j]==false && is_surrounded(M,i,j)==false){
	  surrounded.pixel[i][j] = true;
	}
	//std::cout << " i: " << i << " j: " << j << " M.pixel[i][j]: " << M.pixel[i][j] << " is_surrounded(M,i,j): " << is_surrounded(M,i,j) << " surrounded.pixel[i][j]: " << surrounded.pixel[i][j] << std::endl;
      }
    }
    return surrounded;
  }
  //
  //
  module_bin SiPixelPhase1bumpbonding::convert_to_bool(MonitorElement* h) {
    module_bin M;
    for (int x = 0; x < h->getTH2F()->GetNbinsX(); x++) {
      for (int y = 0; y < h->getTH2F()->GetNbinsY(); y++) {
	std::cout << " x: "<< x <<" y: "<< y << " Content: "<< h->getTH2F()->GetBinContent(x+1, y+1) << std::endl;
	std::vector<uint8_t> roc_roc_x_rocy = get_roc_col_row(uint16_t(x), uint16_t(y));
	std::vector<int> roc_roc_x_rocy_py = get_roc_col_row_py(int(x), int(y));
	//std::cout << "roc_roc_x_rocy[0]: "<< unsigned(roc_roc_x_rocy[0]) <<" roc_roc_x_rocy[1]: "<< unsigned(roc_roc_x_rocy[1]) << " roc_roc_x_rocy[2]: "<< unsigned(roc_roc_x_rocy[2]) << std::endl;
	//std::cout << "roc_roc_x_rocy_py[0]: "<< roc_roc_x_rocy_py[0] <<" roc_roc_x_rocy_py[1]: "<< roc_roc_x_rocy_py[1] << " roc_roc_x_rocy_py[2]: "<< roc_roc_x_rocy_py[2] << std::endl;
	if (h->getTH2F()->GetBinContent(x + 1, y + 1)) {
	  M.roc[roc_roc_x_rocy[0]][roc_roc_x_rocy[1]][roc_roc_x_rocy[2]] = true;
	  M.pixel[y][x] = true;
	} else {
	  M.roc[roc_roc_x_rocy[0]][roc_roc_x_rocy[1]][roc_roc_x_rocy[2]] = false;
	  M.pixel[y][x] = false;
	}
	//std::cout << " M.roc: " << M.roc[roc_roc_x_rocy[0]][roc_roc_x_rocy[1]][roc_roc_x_rocy[2]] << std::endl;
	//std::cout << " M.pixel: " << M.pixel[y][x] << std::endl;
      }
    }
    return M;
  }
  //
  //
  //

  //------------------------------------------------------------------
  // Used to book the MEs
  //------------------------------------------------------------------

  void SiPixelPhase1bumpbonding::bookMEs(DQMStore::IBooker& iBooker) {

    iBooker.cd();
    iBooker.setCurrentFolder(topFolderName_+"/PXForward/HalfCylinder_mI/PXRing_1/PXDisk_-1/SignedBlade_1/");

    residuals_["digi_occupancy_per_col_per_row_FPix_BmI_D1_BLD1_PNL1_RNG1_again"] =
      iBooker.book2D("digi_occupancy_per_col_per_row_FPix_BmI_D1_BLD1_PNL1_RNG1_again",
                       "FPix_BmI_D1_BLD1_PNL1_RNG1",416,-0.5,415.5, 160,-0.5,159.5);

    residuals_["digi_occupancy_per_col_per_row_FPix_BmI_D1_BLD1_PNL1_RNG1_empty"] =
      iBooker.book2D("digi_occupancy_per_col_per_row_FPix_BmI_D1_BLD1_PNL1_RNG1_empty",
                       "FPix_BmI_D1_BLD1_PNL1_RNG1_empty",416,-0.5,415.5, 160,-0.5,159.5);

    //Reset the iBooker
    //iBooker.setCurrentFolder("PixelPhase1/");


      
  }
  
  //------------------------------------------------------------------
  // Fill the MEs
  //------------------------------------------------------------------

  void SiPixelPhase1bumpbonding::fillMEs(DQMStore::IBooker& iBooker, DQMStore::IGetter& iGetter) {

    MonitorElement* digi_1 =
      iGetter.get("PixelPhase1/Phase1_MechanicalView/PXForward/HalfCylinder_mI/PXRing_1/PXDisk_-1/SignedBlade_1/digi_occupancy_per_col_per_row_FPix_BmI_D1_BLD1_PNL1_RNG1");
    //
    modules["digi_occupancy_per_col_per_row_FPix_BmI_D1_BLD1_PNL1_RNG1"] = convert_to_bool(digi_1);
    modules["digi_occupancy_per_col_per_row_FPix_BmI_D1_BLD1_PNL1_RNG1"].name = "digi_occupancy_per_col_per_row_FPix_BmI_D1_BLD1_PNL1_RNG1";
    modulesSurrounded["digi_occupancy_per_col_per_row_FPix_BmI_D1_BLD1_PNL1_RNG1"] = sourroundedRegion(modules["digi_occupancy_per_col_per_row_FPix_BmI_D1_BLD1_PNL1_RNG1"]);

    //MonitorElement* h_ret = digi_1->cloneMEData()
    //for (int x = 0; x < h_ret->getNbinsX(); x++) {
    // for (int y = 0; y < h_ret->getNbinsY(); y++) {
    //        h_ret->setBinContent(x+1,y+1,modulesSurrounded["digi_occupancy_per_col_per_row_FPix_BmI_D1_BLD1_PNL1_RNG1"].pixel[y][x])
    // }
    // }
    //
    
    //
    //module_bin M;
    //
    for (int x = 0; x < digi_1->getNbinsX(); x++) {
      for (int y = 0; y < digi_1->getNbinsY(); y++) {
	/*
	std::cout << "SiPixelPhase1bumpbonding::~SiPixelPhase1bumpbonding: x: "<< x <<" y: "<< y << " Content: "<< digi_1->getBinContent(x+1, y+1) << std::endl;
	std::vector<uint8_t> roc_roc_x_rocy = get_roc_col_row(uint16_t(x), uint16_t(y));
	std::cout << "roc_roc_x_rocy[0]: "<< unsigned(roc_roc_x_rocy[0]) <<" roc_roc_x_rocy[1]: "<< unsigned(roc_roc_x_rocy[1]) << " roc_roc_x_rocy[2]: "<< unsigned(roc_roc_x_rocy[2]) << std::endl;

	if (digi_1->getBinContent(x+1, y+1)) {
	  M.roc[roc_roc_x_rocy[0]][roc_roc_x_rocy[1]][roc_roc_x_rocy[2]] = true;
	} else {
	  M.roc[roc_roc_x_rocy[0]][roc_roc_x_rocy[1]][roc_roc_x_rocy[2]] = false;
	}
	std::cout << " M: " << M.roc[roc_roc_x_rocy[0]][roc_roc_x_rocy[1]][roc_roc_x_rocy[2]] << std::endl;
	*/
	residuals_["digi_occupancy_per_col_per_row_FPix_BmI_D1_BLD1_PNL1_RNG1_again"]->setBinContent(x,y,digi_1->getBinContent(x,y));
	residuals_["digi_occupancy_per_col_per_row_FPix_BmI_D1_BLD1_PNL1_RNG1_empty"]->setBinContent(x+1,y+1,modulesSurrounded["digi_occupancy_per_col_per_row_FPix_BmI_D1_BLD1_PNL1_RNG1"].pixel[y][x]);
	
      }
    }
    //
  }
  
}  //namespace

DEFINE_FWK_MODULE(SiPixelPhase1bumpbonding);
