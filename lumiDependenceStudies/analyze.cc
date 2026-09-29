#include <TNtuple.h>

#include "xjjanauti.h"
#include "lumimessenger.h"

int macro(const std::string& inputname, const std::string& outputname) {
  // Read in input arguments and initialize the LumiMessenger
  auto* inf = TFile::Open(inputname.c_str());
  if (!inf || inf->IsZombie()) { __XJJLOG << "!! failed to open input file, abort." << std::endl; return 1; }
  LumiMessenger msg(inf);
  if (!msg.Tree) { __XJJLOG << "!! no Tree found in input file, abort." << std::endl; return 2; }
  __XJJLOG << ">> " << msg.GetEntries() << " entries loaded from " << inputname << std::endl;

  auto* outf = TFile::Open(outputname.c_str(), "RECREATE");
  if (!outf || outf->IsZombie()) { __XJJLOG << "!! failed to open output file, abort." << std::endl; return 3; }
  outf->cd();

  auto* ntuple = new TNtuple("ntuple", "ntuple", "instLumi:lumiLeveled:ZDCsumPlus:ZDCsumMinus:HFEMaxPlus_eta5:HFEMaxMinus_eta5:nVtx:nTrackInAcceptanceHP");
  // anaylsis loop
  const auto nentries = msg.GetEntries();
  for (Long64_t iEvt = 0; iEvt < nentries; ++iEvt) {
    msg.GetEntry(iEvt);

    if (iEvt % 10000 == 0) {
      __XJJLOG << ">> processing entry " << iEvt << " / " << nentries << std::endl;
    }

    if ((msg.isZeroBias || msg.isZeroBias_Min400_Max10000 || msg.isZeroBias_Max400_Pixel || msg.isZeroBias_Max10000) == false)
        continue;

    ntuple->Fill(msg.instLumi, msg.instLumiLeveled, msg.ZDCsumPlus, msg.ZDCsumMinus, msg.HFEMaxPlus_eta5, msg.HFEMaxMinus_eta5, msg.nVtx, msg.nTrackInAcceptanceHP);

  }
  __XJJLOG << ">> " << ntuple->GetEntries() << " entries passed selection, written to " << outputname << std::endl;

  // make output file, save histograms, etc.
  ntuple->Write("", TObject::kOverwrite);
  outf->Close();
  inf->Close();

  return 0;
}

int main(int argc, char* argv[]) {
  if (argc == 3) {
    return macro(argv[1], argv[2]);
  }
  __XJJLOG << "!! usage: " << argv[0] << " <inputname.root (output of addlumi.cc)> <outputname.root>" << std::endl;
  return 1;
}