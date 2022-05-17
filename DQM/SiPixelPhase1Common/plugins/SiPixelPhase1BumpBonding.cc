// -*- C++ -*-
//
// Package:    SiPixelPhase1BumpBonding
// Class:      SiPixelPhase1BumpBonding
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
  
  class SiPixelPhase1BumpBonding : public DQMEDHarvester {
  public:
    explicit SiPixelPhase1BumpBonding(const edm::ParameterSet& conf);
    ~SiPixelPhase1BumpBonding() override;


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
    std::map<std::string, MonitorElement*> badbonding_;
    std::map<std::string, MonitorElement*> deadsize_;
    std::map<std::string, MonitorElement*> deadroc_;
    map<std::string, module_bin> modules;
    map<std::string, module_bin> modulesSurrounded;
    
    //Book Monitoring Elements
    void bookMEs(DQMStore::IBooker& iBooker);

    //
    std::vector<std::string> getDirList();
    std::vector<uint8_t> get_roc_col_row(uint16_t x, uint16_t y);
    std::vector<int> get_roc_col_row_py(int x, int y);
    bool is_safe(module_bin& M, int row, int col);
    bool is_surrounded(module_bin& M, int row, int col);
    module_bin sourroundedRegion(module_bin& M);
    module_bin convert_to_bool(MonitorElement* h);
    int dfs(module_bin& M, int row, int col, module_bin& visited, int count);
    std::vector<std::vector<int>> largestRegion(module_bin& M);
    
  
    //Fill Monitoring Elements
    void fillMEs(DQMStore::IBooker& iBooker, DQMStore::IGetter& iGetter);

  }; //class

  //---------------------------
  SiPixelPhase1BumpBonding::SiPixelPhase1BumpBonding(const edm::ParameterSet& iConfig)
    : DQMEDHarvester(iConfig), conf_(iConfig), firstLumi(true) {
    LogInfo("PixelDQM") << "SiPixelPhase1BumpBonding::SiPixelPhase1BumpBonding: Got DQM BackEnd interface" << endl;
    topFolderName_ = conf_.getParameter<std::string>("TopFolderName");
    minHits_ = conf_.getParameter<int>("MinHits");
  }

  SiPixelPhase1BumpBonding::~SiPixelPhase1BumpBonding() {
    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)
    LogInfo("PixelDQM") << "SiPixelPhase1BumpBonding::~SiPixelPhase1BumpBonding: Destructor" << endl;
  }
  
  void SiPixelPhase1BumpBonding::beginRun(edm::Run const& run, edm::EventSetup const& eSetup) {}
  
  
  void SiPixelPhase1BumpBonding::dqmEndLuminosityBlock(DQMStore::IBooker& iBooker,
						       DQMStore::IGetter& iGetter,
						       const edm::LuminosityBlock& lumiSeg,
						       edm::EventSetup const& c) {
    if(firstLumi){
      bookMEs(iBooker);
      firstLumi=false;
    }
    fillMEs(iBooker, iGetter);
  }

  void SiPixelPhase1BumpBonding::dqmEndJob(DQMStore::IBooker& iBooker, DQMStore::IGetter& iGetter) {

  }

  //-------------------------
  //-----------Custom functions--------------------------------------
  //
  //-----------------------------------------------------------------
  
  std::vector<std::string> SiPixelPhase1BumpBonding::getDirList() {

    std::string histName;
    std::string dirpath;
    std::vector<std::string> alldirs = {};
    
    for (auto itcy : {"mI", "mO", "pI", "pO"}) {  //PXForward
      for (auto itr : {"1", "2"}) {
	for (auto itd : {"1", "2", "3"}) {
	  for (auto itsb : {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11"}) {
	    for (auto itnl : {"1", "2"}) {
	      if (std::string(itcy) == "mI"){
		dirpath = "PixelPhase1/Phase1_MechanicalView/PXForward/HalfCylinder_"+std::string(itcy)+"/PXRing_"+std::string(itr)+"/PXDisk_-"+std::string(itd)+
		  "/SignedBlade_"+std::string(itsb);
		histName = "digi_occupancy_per_col_per_row_FPix_B"+std::string(itcy)+"_D"+std::string(itd)+"_BLD"+std::string(itsb)+"_PNL"+std::string(itnl)+
		  "_RNG"+std::string(itr);
	      }
	      if (std::string(itcy) == "mO"){
		dirpath = "PixelPhase1/Phase1_MechanicalView/PXForward/HalfCylinder_"+std::string(itcy)+"/PXRing_"+std::string(itr)+"/PXDisk_-"+std::string(itd)+
		  "/SignedBlade_-"+std::string(itsb);
		histName = "digi_occupancy_per_col_per_row_FPix_B"+std::string(itcy)+"_D"+std::string(itd)+"_BLD"+std::string(itsb)+"_PNL"+std::string(itnl)+
		  "_RNG"+std::string(itr);
	      }
	      if (std::string(itcy) == "pI"){
		dirpath = "PixelPhase1/Phase1_MechanicalView/PXForward/HalfCylinder_"+std::string(itcy)+"/PXRing_"+std::string(itr)+"/PXDisk_+"+std::string(itd)+
		  "/SignedBlade_"+std::string(itsb);
		histName = "digi_occupancy_per_col_per_row_FPix_B"+std::string(itcy)+"_D"+std::string(itd)+"_BLD"+std::string(itsb)+"_PNL"+std::string(itnl)+
		  "_RNG"+std::string(itr);
	      }
	      if (std::string(itcy) == "pO"){
		dirpath = "PixelPhase1/Phase1_MechanicalView/PXForward/HalfCylinder_"+std::string(itcy)+"/PXRing_"+std::string(itr)+"/PXDisk_+"+std::string(itd)+
		  "/SignedBlade_-"+std::string(itsb);
		histName = "digi_occupancy_per_col_per_row_FPix_B"+std::string(itcy)+"_D"+std::string(itd)+"_BLD"+std::string(itsb)+"_PNL"+std::string(itnl)+
		  "_RNG"+std::string(itr);
	      }
	      //
	      alldirs.push_back(dirpath + "/" + histName);
	      //
	    }
	  }
	}
      }
    }
    return alldirs;
  }
  std::vector<uint8_t> SiPixelPhase1BumpBonding::get_roc_col_row(uint16_t x, uint16_t y) {
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
  std::vector<int> SiPixelPhase1BumpBonding::get_roc_col_row_py(int x, int y) {
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
  bool SiPixelPhase1BumpBonding::is_safe(module_bin& M, int row, int col){
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
  bool SiPixelPhase1BumpBonding::is_surrounded(module_bin& M, int row, int col){
    const static int rowNbr[8] = {-1, -1, -1, 0, 0, 1, 1, 1};
    const static int colNbr[8] = {-1, 0, 1, -1, 1, -1, 0, 1};
    int i_surrounded=0;
    //
    //
    for (int k = 0; k < 8; k++){
      if (is_safe(M, row + rowNbr[k], col + colNbr[k])==true && M.pixel[row + rowNbr[k]][col + colNbr[k]]==false){
	i_surrounded = i_surrounded + 1;
      }else if (is_safe(M, row + rowNbr[k], col + colNbr[k])==false){
	i_surrounded = i_surrounded + 1;
      }else {
	return false;
      }
    }
    if(i_surrounded==8){
      return true;
    }
    return false;
  }
  //
  module_bin SiPixelPhase1BumpBonding::sourroundedRegion(module_bin& M) {
    module_bin surrounded;
    for (int i = 0; i < ROW; ++i) {
      for (int j = 0; j < COL; ++j) {
	if(M.pixel[i][j]==false && is_surrounded(M,i,j)==true){
	  surrounded.pixel[i][j] = true;
	}
	if(M.pixel[i][j]==true && is_surrounded(M,i,j)==true){
	  surrounded.pixel[i][j] = false;
	}
	if(M.pixel[i][j]==false && is_surrounded(M,i,j)==false){
	  surrounded.pixel[i][j] = false;
	}
	if(M.pixel[i][j]==true && is_surrounded(M,i,j)==false){
	  surrounded.pixel[i][j] = false;
	}
      }
    }
    return surrounded;
  }
  //
  //
  module_bin SiPixelPhase1BumpBonding::convert_to_bool(MonitorElement* h) {
    module_bin M;
    for (int x = 0; x < h->getTH2F()->GetNbinsX(); x++) {
      for (int y = 0; y < h->getTH2F()->GetNbinsY(); y++) {
	std::vector<uint8_t> roc_roc_x_rocy = get_roc_col_row(uint16_t(x), uint16_t(y));
	std::vector<int> roc_roc_x_rocy_py = get_roc_col_row_py(int(x), int(y));
	//
	if (h->getTH2F()->GetBinContent(x + 1, y + 1)) {
	  M.roc[roc_roc_x_rocy[0]][roc_roc_x_rocy[1]][roc_roc_x_rocy[2]] = true;
	  M.pixel[y][x] = true;
	} else {
	  M.roc[roc_roc_x_rocy[0]][roc_roc_x_rocy[1]][roc_roc_x_rocy[2]] = false;
	  M.pixel[y][x] = false;
	}
	//
      }
    }
    return M;
  }
  //
  //
  int SiPixelPhase1BumpBonding::dfs(module_bin& M, int row, int col, module_bin& visited, int count){
    
    const static int rowNbr[8] = {-1, -1, -1, 0, 0, 1, 1, 1};
    const static int colNbr[8] = {-1, 0, 1, -1, 1, -1, 0, 1};
    //
    visited.pixel[row][col] = true; //current_label;
    for (int k = 0; k < 8; k++){
      //
      //
      if (is_safe(M, row + rowNbr[k], col + colNbr[k]) == true && M.pixel[row + rowNbr[k]][col + colNbr[k]] == true && visited.pixel[row + rowNbr[k]][col + colNbr[k]] == false){
	count++;
	count = dfs(M, row + rowNbr[k], col + colNbr[k], visited , count);
      } else {
	continue;
      }
    }
    return count;
  }
  //
  //
  std::vector<std::vector<int>> SiPixelPhase1BumpBonding::largestRegion(module_bin& M) {
    module_bin visited;
    //
    for (int i = 0; i < ROW; ++i) {
      for (int j = 0; j < COL; ++j) {
	//
	visited.pixel[i][j] = false;
      }
    }
    //
    vector<vector<int> > results;
    int theshold = 50;
    for (int i = 0; i < ROW; ++i) {
      for (int j = 0; j < COL; ++j) {
	//
	//
	if(M.pixel[i][j]==true && visited.pixel[i][j]==false){
	  int count = 1;
	  count = dfs(M, i, j, visited , count);
	  if (count > theshold){
	    vector<int> v1 = {count,i,j};
	    results.push_back(v1);
	  }
	}
	
      }
    }
    return results;
  }
  //

  //------------------------------------------------------------------
  // Used to book the MEs
  //------------------------------------------------------------------

  void SiPixelPhase1BumpBonding::bookMEs(DQMStore::IBooker& iBooker) {
    //
    std::string digihistName;
    std::string digipathName;
    std::vector<std::string> fullpathNames = getDirList();
    std::cout << " fullpathNames size inside bookME: " << fullpathNames.size() << std::endl;

    for (unsigned int pathIt = 0; pathIt < fullpathNames.size(); pathIt++) {

      digihistName = fullpathNames[pathIt].substr(fullpathNames[pathIt].rfind("/") + 1, fullpathNames[pathIt].length());
      digipathName = fullpathNames[pathIt].substr(0, fullpathNames[pathIt].find("/digi"));

      //std::cout << "digihistName: " << fullpathNames[pathIt].substr(fullpathNames[pathIt].rfind("/") + 1, fullpathNames[pathIt].length()) << std::endl;
      //std::cout << "digipathName: " << fullpathNames[pathIt].substr(0, fullpathNames[pathIt].find("/digi")) << std::endl; 

      iBooker.cd();
      iBooker.setCurrentFolder(digipathName+"/");
      badbonding_[digihistName+"_deadRegion"] = iBooker.book2D(digihistName+"_deadRegion", digihistName+"_deadRegion;col;row", 416,-0.5,415.5, 160,-0.5,159.5);
      //
      deadsize_[digihistName+"_deadRegionSize"] = iBooker.book1D(digihistName+"_deadRegionSize", digihistName+"_deadRegionSize;Size Empty region", 100,0,4100);

      for (int r = 0; r < 16; r++){
	deadroc_[digihistName+"_deadRegion_ROC_"+std::to_string(r)] = iBooker.book1D(digihistName+"_deadRegion_ROC_"+std::to_string(r), digihistName+"_deadRegion_ROC_"+std::to_string(r)+
										     ";Empty Area size/pixel", 100, 0, 4100);
      }
      //
    }
    //
    //Reset the iBooker
    //iBooker.setCurrentFolder("PixelPhase1/");


      
  }
  
  //------------------------------------------------------------------
  // Fill the MEs
  //------------------------------------------------------------------

  void SiPixelPhase1BumpBonding::fillMEs(DQMStore::IBooker& iBooker, DQMStore::IGetter& iGetter) {

    //
    std::vector<std::string> fullpathNames = getDirList();

    std::cout << " fullpathNames size inside fillMEs: " << fullpathNames.size() << std::endl;
    
    for (unsigned int pathIt = 0; pathIt < fullpathNames.size(); pathIt++) {

      std::string digihistName = fullpathNames[pathIt].substr(fullpathNames[pathIt].rfind("/") + 1, fullpathNames[pathIt].length());
      MonitorElement* h_digi = iGetter.get(fullpathNames[pathIt]);
      //
      //
      modules[digihistName] = convert_to_bool(h_digi);
      modules[digihistName].name = digihistName;
      modulesSurrounded[digihistName] = sourroundedRegion(modules[digihistName]);
      //
      std::cout << " digihistName: " << digihistName << std::endl;
      std::vector<std::vector<int>> res = largestRegion(modulesSurrounded[digihistName]);
      //
      for (unsigned int i = 0; i < res.size(); i++) {
	for (unsigned int j = 0; j < res[i].size(); j++)
	  std::cout << res[i][j] << " ";
	std::cout << std::endl;
	std::cout << "size:" << res[i][0] << std::endl;
	if (res[i][0] < 52*80 ) {
	  deadsize_[digihistName+"_deadRegionSize"]->Fill(res[i][0]);
	  std::vector<int> rocdetails = get_roc_col_row_py(res[i][2], res[i][1]);
	  //std::cout << "roc: " << rocdetails[0] << " col X_roc: " << rocdetails[1] << " row Y_roc: " << rocdetails[2] << std::endl ;
	  for (int r = 0; r < 16; r++){
	    if (rocdetails[0] == r) 
	      deadroc_[digihistName+"_deadRegion_ROC_"+std::to_string(r)]->Fill(res[i][0]);
	  }
	}
	  
      }
      //
      for (int x = 0; x < h_digi->getNbinsX(); x++) {
	for (int y = 0; y < h_digi->getNbinsY(); y++) {
	  badbonding_[digihistName+"_deadRegion"]->setBinContent(x+1,y+1,modulesSurrounded[digihistName].pixel[y][x]);
	  
	}
      }
    }
      //
  }
  
}  //namespace

DEFINE_FWK_MODULE(SiPixelPhase1BumpBonding);
