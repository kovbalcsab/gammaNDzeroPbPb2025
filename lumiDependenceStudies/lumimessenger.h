#pragma once

// Reads the "Tree" produced by addlumi.cc: all branches of its input forest
// (e.g. /eos/cms/store/group/phys_heavyions/wangj/Forest2025PbPb/Dzero_*.root)
// plus the instLumi* branches (instantaneous luminosity + lumileveling flag,
// from instlumi.py) that addlumi.cc appends. Exposes them all as plain
// member variables.

#include <vector>

#include <TFile.h>
#include <TTree.h>

class LumiMessenger
{
public:
  TTree* Tree = nullptr;

  // event / vertex
  Int_t Run = 0;
  Long64_t Event = 0;
  Int_t Lumi = 0;
  Int_t ProcessID = 0;
  Int_t clusComp_nPixHits = 0;
  Double_t clusComp_quality = 0.;
  Float_t VX = 0., VY = 0., VZ = 0.;
  Float_t VXError = 0., VYError = 0., VZError = 0.;
  Int_t nVtx = 0;

  // trigger / filter flags
  Bool_t isL1ZDCOr = false, isL1ZDCOr_Min400_Max10000 = false, isL1ZDCOr_Max400_Pixel = false, isL1ZDCOr_Max10000 = false;
  Bool_t isL1ZDCXORJet8 = false, isL1ZDCXORJet12 = false, isL1ZDCXORJet16 = false;
  Bool_t isZeroBias = false, isZeroBias_Min400_Max10000 = false, isZeroBias_Max400_Pixel = false, isZeroBias_Max10000 = false;
  Bool_t isNotBptxOR = false, isUnpairedBunchBptxMinus = false, isUnpairedBunchBptxPlus = false;
  Bool_t selectedBkgFilter = false, selectedVtxFilter = false, ClusterCompatibilityFilter = false, cscTightHalo2015Filter = false;

  // rapidity gap / photonuclear tagging
  Bool_t ZDCgammaN = false, ZDCNgamma = false, gapgammaN = false, gapNgamma = false;
  std::vector<bool>* gammaN = nullptr;
  std::vector<bool>* Ngamma = nullptr;

  // ZDC / HF energies
  Float_t ZDCsumPlus = 0., ZDCsumMinus = 0.;
  Float_t HFEMaxPlus = 0., HFEMaxPlus_forest = 0., HFEMaxPlus_eta5 = 0., HFEMaxPlus_pt0p1 = 0.;
  Float_t HFEMaxMinus = 0., HFEMaxMinus_forest = 0., HFEMaxMinus_eta5 = 0., HFEMaxMinus_pt0p1 = 0.;

  Int_t nTrackInAcceptanceHP = 0;

  // Dzero candidates
  Int_t Dsize = 0;
  std::vector<float>* Dpt = nullptr;
  std::vector<float>* Dy = nullptr;
  std::vector<float>* Dmass = nullptr;

  std::vector<float>* Dtrk1Pt = nullptr;
  std::vector<float>* Dtrk1PtErr = nullptr;
  std::vector<float>* Dtrk1Eta = nullptr;
  std::vector<float>* Dtrk1dedx = nullptr;
  std::vector<float>* Dtrk1MassHypo = nullptr;
  std::vector<float>* Dtrk1PixelHit = nullptr;
  std::vector<float>* Dtrk1StripHit = nullptr;
  std::vector<float>* Dtrk1P = nullptr;
  std::vector<float>* Dtrk1PionScore = nullptr;
  std::vector<float>* Dtrk1KaonScore = nullptr;
  std::vector<float>* Dtrk1ProtScore = nullptr;

  std::vector<float>* Dtrk2Pt = nullptr;
  std::vector<float>* Dtrk2PtErr = nullptr;
  std::vector<float>* Dtrk2Eta = nullptr;
  std::vector<float>* Dtrk2dedx = nullptr;
  std::vector<float>* Dtrk2MassHypo = nullptr;
  std::vector<float>* Dtrk2PixelHit = nullptr;
  std::vector<float>* Dtrk2StripHit = nullptr;
  std::vector<float>* Dtrk2P = nullptr;
  std::vector<float>* Dtrk2PionScore = nullptr;
  std::vector<float>* Dtrk2KaonScore = nullptr;
  std::vector<float>* Dtrk2ProtScore = nullptr;

  std::vector<float>* Dchi2cl = nullptr;
  std::vector<float>* DsvpvDistance = nullptr;
  std::vector<float>* DsvpvDisErr = nullptr;
  std::vector<float>* DsvpvDistance_2D = nullptr;
  std::vector<float>* DsvpvDisErr_2D = nullptr;
  std::vector<float>* Dip3D = nullptr;
  std::vector<float>* Dip3derr = nullptr;
  std::vector<float>* Dalpha = nullptr;
  std::vector<float>* Ddtheta = nullptr;

  std::vector<bool>* DpassCut23PAS = nullptr;
  std::vector<bool>* DpassCut23LowPt = nullptr;
  std::vector<bool>* DpassCut23PASSystDsvpvSig = nullptr;
  std::vector<bool>* DpassCut23PASSystDtrkPt = nullptr;
  std::vector<bool>* DpassCut23PASSystDalpha = nullptr;
  std::vector<bool>* DpassCut23PASSystDchi2cl = nullptr;
  std::vector<bool>* DpassCutNominal = nullptr;
  std::vector<bool>* DpassCutLoose = nullptr;
  std::vector<bool>* DpassCutSystDsvpvSig = nullptr;
  std::vector<bool>* DpassCutSystDtrkPt = nullptr;
  std::vector<bool>* DpassCutSystDalpha = nullptr;
  std::vector<bool>* DpassCutSystDalphaDdtheta = nullptr;
  std::vector<bool>* DpassCutSystDdtheta = nullptr;
  std::vector<bool>* DpassCutSystDchi2cl = nullptr;

  std::vector<float>* Dmva_BDT = nullptr;
  std::vector<float>* Dmva_BDTG = nullptr;

  std::vector<int>* Dgen = nullptr;
  std::vector<bool>* DisSignalCalc = nullptr;
  std::vector<bool>* DisSignalCalcPrompt = nullptr;
  std::vector<bool>* DisSignalCalcFeeddown = nullptr;

  // gen-level (MC only)
  Int_t Gsize = 0;
  std::vector<float>* Gpt = nullptr;
  std::vector<float>* Gy = nullptr;
  std::vector<bool>* GisSignalCalc = nullptr;
  std::vector<bool>* GisSignalCalcPrompt = nullptr;
  std::vector<bool>* GisSignalCalcFeeddown = nullptr;

  // instantaneous luminosity, added by addlumi.cc (from instlumi.py); all in 10^33 cm^-2 s^-1
  Double_t instLumiInit = 0.;
  Double_t instLumiEnd = 0.;
  Double_t instLumi = 0.;
  Bool_t instLumiLeveled = false; // true if instlumi.py classified this lumisection's lumi trend as "leveled"

public:
  LumiMessenger(TFile& File, std::string TreeName = "Tree");
  LumiMessenger(TFile* File, std::string TreeName = "Tree");
  LumiMessenger(TTree* InputTree = nullptr);
  bool Initialize(TTree* InputTree);
  bool Initialize();
  bool GetEntry(Long64_t iEntry);
  Long64_t GetEntries() const;
};

LumiMessenger::LumiMessenger(TFile& File, std::string TreeName)
{
  Tree = (TTree*)File.Get(TreeName.c_str());
  Initialize();
}

LumiMessenger::LumiMessenger(TFile* File, std::string TreeName)
{
  Tree = (File != nullptr) ? (TTree*)File->Get(TreeName.c_str()) : nullptr;
  Initialize();
}

LumiMessenger::LumiMessenger(TTree* InputTree)
{
  Initialize(InputTree);
}

bool LumiMessenger::Initialize(TTree* InputTree)
{
  Tree = InputTree;
  return Initialize();
}

bool LumiMessenger::Initialize()
{
  if (Tree == nullptr) return false;

  Tree->SetBranchAddress("Run", &Run);
  Tree->SetBranchAddress("Event", &Event);
  Tree->SetBranchAddress("Lumi", &Lumi);
  Tree->SetBranchAddress("ProcessID", &ProcessID);
  Tree->SetBranchAddress("clusComp_nPixHits", &clusComp_nPixHits);
  Tree->SetBranchAddress("clusComp_quality", &clusComp_quality);
  Tree->SetBranchAddress("VX", &VX);
  Tree->SetBranchAddress("VY", &VY);
  Tree->SetBranchAddress("VZ", &VZ);
  Tree->SetBranchAddress("VXError", &VXError);
  Tree->SetBranchAddress("VYError", &VYError);
  Tree->SetBranchAddress("VZError", &VZError);
  Tree->SetBranchAddress("nVtx", &nVtx);

  Tree->SetBranchAddress("isL1ZDCOr", &isL1ZDCOr);
  Tree->SetBranchAddress("isL1ZDCOr_Min400_Max10000", &isL1ZDCOr_Min400_Max10000);
  Tree->SetBranchAddress("isL1ZDCOr_Max400_Pixel", &isL1ZDCOr_Max400_Pixel);
  Tree->SetBranchAddress("isL1ZDCOr_Max10000", &isL1ZDCOr_Max10000);
  Tree->SetBranchAddress("isL1ZDCXORJet8", &isL1ZDCXORJet8);
  Tree->SetBranchAddress("isL1ZDCXORJet12", &isL1ZDCXORJet12);
  Tree->SetBranchAddress("isL1ZDCXORJet16", &isL1ZDCXORJet16);
  Tree->SetBranchAddress("isZeroBias", &isZeroBias);
  Tree->SetBranchAddress("isZeroBias_Min400_Max10000", &isZeroBias_Min400_Max10000);
  Tree->SetBranchAddress("isZeroBias_Max400_Pixel", &isZeroBias_Max400_Pixel);
  Tree->SetBranchAddress("isZeroBias_Max10000", &isZeroBias_Max10000);
  Tree->SetBranchAddress("isNotBptxOR", &isNotBptxOR);
  Tree->SetBranchAddress("isUnpairedBunchBptxMinus", &isUnpairedBunchBptxMinus);
  Tree->SetBranchAddress("isUnpairedBunchBptxPlus", &isUnpairedBunchBptxPlus);
  Tree->SetBranchAddress("selectedBkgFilter", &selectedBkgFilter);
  Tree->SetBranchAddress("selectedVtxFilter", &selectedVtxFilter);
  Tree->SetBranchAddress("ClusterCompatibilityFilter", &ClusterCompatibilityFilter);
  Tree->SetBranchAddress("cscTightHalo2015Filter", &cscTightHalo2015Filter);

  Tree->SetBranchAddress("ZDCgammaN", &ZDCgammaN);
  Tree->SetBranchAddress("ZDCNgamma", &ZDCNgamma);
  Tree->SetBranchAddress("gapgammaN", &gapgammaN);
  Tree->SetBranchAddress("gapNgamma", &gapNgamma);
  Tree->SetBranchAddress("gammaN", &gammaN);
  Tree->SetBranchAddress("Ngamma", &Ngamma);

  Tree->SetBranchAddress("ZDCsumPlus", &ZDCsumPlus);
  Tree->SetBranchAddress("ZDCsumMinus", &ZDCsumMinus);
  Tree->SetBranchAddress("HFEMaxPlus", &HFEMaxPlus);
  Tree->SetBranchAddress("HFEMaxPlus_forest", &HFEMaxPlus_forest);
  Tree->SetBranchAddress("HFEMaxPlus_eta5", &HFEMaxPlus_eta5);
  Tree->SetBranchAddress("HFEMaxPlus_pt0p1", &HFEMaxPlus_pt0p1);
  Tree->SetBranchAddress("HFEMaxMinus", &HFEMaxMinus);
  Tree->SetBranchAddress("HFEMaxMinus_forest", &HFEMaxMinus_forest);
  Tree->SetBranchAddress("HFEMaxMinus_eta5", &HFEMaxMinus_eta5);
  Tree->SetBranchAddress("HFEMaxMinus_pt0p1", &HFEMaxMinus_pt0p1);

  Tree->SetBranchAddress("nTrackInAcceptanceHP", &nTrackInAcceptanceHP);

  Tree->SetBranchAddress("Dsize", &Dsize);
  Tree->SetBranchAddress("Dpt", &Dpt);
  Tree->SetBranchAddress("Dy", &Dy);
  Tree->SetBranchAddress("Dmass", &Dmass);

  Tree->SetBranchAddress("Dtrk1Pt", &Dtrk1Pt);
  Tree->SetBranchAddress("Dtrk1PtErr", &Dtrk1PtErr);
  Tree->SetBranchAddress("Dtrk1Eta", &Dtrk1Eta);
  Tree->SetBranchAddress("Dtrk1dedx", &Dtrk1dedx);
  Tree->SetBranchAddress("Dtrk1MassHypo", &Dtrk1MassHypo);
  Tree->SetBranchAddress("Dtrk1PixelHit", &Dtrk1PixelHit);
  Tree->SetBranchAddress("Dtrk1StripHit", &Dtrk1StripHit);
  Tree->SetBranchAddress("Dtrk1P", &Dtrk1P);
  Tree->SetBranchAddress("Dtrk1PionScore", &Dtrk1PionScore);
  Tree->SetBranchAddress("Dtrk1KaonScore", &Dtrk1KaonScore);
  Tree->SetBranchAddress("Dtrk1ProtScore", &Dtrk1ProtScore);

  Tree->SetBranchAddress("Dtrk2Pt", &Dtrk2Pt);
  Tree->SetBranchAddress("Dtrk2PtErr", &Dtrk2PtErr);
  Tree->SetBranchAddress("Dtrk2Eta", &Dtrk2Eta);
  Tree->SetBranchAddress("Dtrk2dedx", &Dtrk2dedx);
  Tree->SetBranchAddress("Dtrk2MassHypo", &Dtrk2MassHypo);
  Tree->SetBranchAddress("Dtrk2PixelHit", &Dtrk2PixelHit);
  Tree->SetBranchAddress("Dtrk2StripHit", &Dtrk2StripHit);
  Tree->SetBranchAddress("Dtrk2P", &Dtrk2P);
  Tree->SetBranchAddress("Dtrk2PionScore", &Dtrk2PionScore);
  Tree->SetBranchAddress("Dtrk2KaonScore", &Dtrk2KaonScore);
  Tree->SetBranchAddress("Dtrk2ProtScore", &Dtrk2ProtScore);

  Tree->SetBranchAddress("Dchi2cl", &Dchi2cl);
  Tree->SetBranchAddress("DsvpvDistance", &DsvpvDistance);
  Tree->SetBranchAddress("DsvpvDisErr", &DsvpvDisErr);
  Tree->SetBranchAddress("DsvpvDistance_2D", &DsvpvDistance_2D);
  Tree->SetBranchAddress("DsvpvDisErr_2D", &DsvpvDisErr_2D);
  Tree->SetBranchAddress("Dip3D", &Dip3D);
  Tree->SetBranchAddress("Dip3derr", &Dip3derr);
  Tree->SetBranchAddress("Dalpha", &Dalpha);
  Tree->SetBranchAddress("Ddtheta", &Ddtheta);

  Tree->SetBranchAddress("DpassCut23PAS", &DpassCut23PAS);
  Tree->SetBranchAddress("DpassCut23LowPt", &DpassCut23LowPt);
  Tree->SetBranchAddress("DpassCut23PASSystDsvpvSig", &DpassCut23PASSystDsvpvSig);
  Tree->SetBranchAddress("DpassCut23PASSystDtrkPt", &DpassCut23PASSystDtrkPt);
  Tree->SetBranchAddress("DpassCut23PASSystDalpha", &DpassCut23PASSystDalpha);
  Tree->SetBranchAddress("DpassCut23PASSystDchi2cl", &DpassCut23PASSystDchi2cl);
  Tree->SetBranchAddress("DpassCutNominal", &DpassCutNominal);
  Tree->SetBranchAddress("DpassCutLoose", &DpassCutLoose);
  Tree->SetBranchAddress("DpassCutSystDsvpvSig", &DpassCutSystDsvpvSig);
  Tree->SetBranchAddress("DpassCutSystDtrkPt", &DpassCutSystDtrkPt);
  Tree->SetBranchAddress("DpassCutSystDalpha", &DpassCutSystDalpha);
  Tree->SetBranchAddress("DpassCutSystDalphaDdtheta", &DpassCutSystDalphaDdtheta);
  Tree->SetBranchAddress("DpassCutSystDdtheta", &DpassCutSystDdtheta);
  Tree->SetBranchAddress("DpassCutSystDchi2cl", &DpassCutSystDchi2cl);

  Tree->SetBranchAddress("Dmva_BDT", &Dmva_BDT);
  Tree->SetBranchAddress("Dmva_BDTG", &Dmva_BDTG);

  Tree->SetBranchAddress("Dgen", &Dgen);
  Tree->SetBranchAddress("DisSignalCalc", &DisSignalCalc);
  Tree->SetBranchAddress("DisSignalCalcPrompt", &DisSignalCalcPrompt);
  Tree->SetBranchAddress("DisSignalCalcFeeddown", &DisSignalCalcFeeddown);

  Tree->SetBranchAddress("Gsize", &Gsize);
  Tree->SetBranchAddress("Gpt", &Gpt);
  Tree->SetBranchAddress("Gy", &Gy);
  Tree->SetBranchAddress("GisSignalCalc", &GisSignalCalc);
  Tree->SetBranchAddress("GisSignalCalcPrompt", &GisSignalCalcPrompt);
  Tree->SetBranchAddress("GisSignalCalcFeeddown", &GisSignalCalcFeeddown);

  Tree->SetBranchAddress("instLumiInit", &instLumiInit);
  Tree->SetBranchAddress("instLumiEnd", &instLumiEnd);
  Tree->SetBranchAddress("instLumi", &instLumi);
  Tree->SetBranchAddress("instLumiLeveled", &instLumiLeveled);

  return true;
}

bool LumiMessenger::GetEntry(Long64_t iEntry)
{
  if (Tree == nullptr) return false;
  Tree->GetEntry(iEntry);
  return true;
}

Long64_t LumiMessenger::GetEntries() const
{
  return Tree ? Tree->GetEntries() : 0;
}
