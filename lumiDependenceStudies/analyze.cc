#include <TNtuple.h>

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>

#include "xjjanauti.h"
#include "lumimessenger.h"

// triggerSelection:
//   0: any ZB UPC path
//  -1: isNotBptxOR or isUnpairedBunchBptxMinus or isUnpairedBunchBptxPlus
//  -2: isUnpairedBunchBptxMinus
//  -3: isUnpairedBunchBptxPlus
//  -4: isNotBptxOR
static bool validTriggerSelection(int sel) { return sel <= 0 && sel >= -4; }

// good run / lumi section ranges from a CMS certification json: {"<run>": [[first, last], ...], ...}
// (minimal parser for exactly this shape, avoids a json library dependency)
using GoodLumiMap = std::map<int, std::vector<std::pair<int, int>>>;

static bool loadGoodJson(const std::string& path, GoodLumiMap& good) {
  std::ifstream in(path);
  if (!in) return false;
  std::stringstream ss;
  ss << in.rdbuf();
  const std::string txt = ss.str();

  int run = -1;
  std::vector<int> nums;
  auto flush = [&]() {
    if (run >= 0)
      for (size_t i = 0; i + 1 < nums.size(); i += 2) good[run].emplace_back(nums[i], nums[i + 1]);
    nums.clear();
  };
  for (size_t i = 0; i < txt.size();) {
    const char ch = txt[i];
    if (ch == '"') { // run key
      const size_t j = txt.find('"', i + 1);
      if (j == std::string::npos) return false;
      flush();
      run = std::atoi(txt.substr(i + 1, j - i - 1).c_str());
      i = j + 1;
    } else if (ch >= '0' && ch <= '9') {
      char* end = nullptr;
      nums.push_back((int)std::strtol(txt.c_str() + i, &end, 10));
      i = end - txt.c_str();
    } else {
      ++i;
    }
  }
  flush();
  return !good.empty();
}

static bool isGoodLumi(const GoodLumiMap& good, int run, int lumi) {
  const auto it = good.find(run);
  if (it == good.end()) return false;
  for (const auto& r : it->second)
    if (lumi >= r.first && lumi <= r.second) return true;
  return false;
}

int macro(const std::string& inputname, const std::string& outputname, int triggerSelection = 0, const std::string& jsonname = "") {
  if (!validTriggerSelection(triggerSelection)) {
    __XJJLOG << "!! invalid triggerSelection " << triggerSelection << " (must be 0, -1, -2, -3 or -4), abort." << std::endl;
    return 4;
  }
  if (inputname == outputname) { __XJJLOG << "!! input and output are the same file, abort." << std::endl; return 5; }

  GoodLumiMap goodLumis;
  const bool useJson = !jsonname.empty();
  if (useJson) {
    if (!loadGoodJson(jsonname, goodLumis)) { __XJJLOG << "!! failed to read/parse json " << jsonname << ", abort." << std::endl; return 6; }
    __XJJLOG << ">> good-lumi json loaded from " << jsonname << " (" << goodLumis.size() << " runs)" << std::endl;
  }

  // Read in input arguments and initialize the LumiMessenger
  auto* inf = TFile::Open(inputname.c_str());
  if (!inf || inf->IsZombie()) { __XJJLOG << "!! failed to open input file, abort." << std::endl; return 1; }
  LumiMessenger msg(inf);
  if (!msg.Tree) { __XJJLOG << "!! no Tree found in input file, abort." << std::endl; return 2; }
  __XJJLOG << ">> " << msg.GetEntries() << " entries loaded from " << inputname << ", triggerSelection = " << triggerSelection << std::endl;

  auto* outf = TFile::Open(outputname.c_str(), "RECREATE");
  if (!outf || outf->IsZombie()) { __XJJLOG << "!! failed to open output file, abort." << std::endl; return 3; }
  outf->cd();

  auto* ntuple = new TNtuple("ntuple", "ntuple", "instLumi:lumiLeveled:ZDCsumPlus:ZDCsumMinus:HFEMaxPlus_eta5:HFEMaxMinus_eta5:nVtx:nTrackInAcceptanceHP:Run:Lumi");
  // anaylsis loop
  const auto nentries = msg.GetEntries();
  for (Long64_t iEvt = 0; iEvt < nentries; ++iEvt) {
    msg.GetEntry(iEvt);

    if (iEvt % 10000 == 0) {
      __XJJLOG << ">> processing entry " << iEvt << " / " << nentries << std::endl;
    }

    if (useJson && !isGoodLumi(goodLumis, msg.Run, msg.Lumi)) continue;

    bool pass = false;
    switch (triggerSelection) {
      case 0:  pass = msg.isZeroBias || msg.isZeroBias_Min400_Max10000 || msg.isZeroBias_Max400_Pixel || msg.isZeroBias_Max10000; break;
      case -1: pass = msg.isNotBptxOR || msg.isUnpairedBunchBptxMinus || msg.isUnpairedBunchBptxPlus; break;
      case -2: pass = msg.isUnpairedBunchBptxMinus; break;
      case -3: pass = msg.isUnpairedBunchBptxPlus; break;
      case -4: pass = msg.isNotBptxOR; break;
    }
    if (!pass) continue;

    ntuple->Fill(msg.instLumi, msg.instLumiLeveled, msg.ZDCsumPlus, msg.ZDCsumMinus, msg.HFEMaxPlus_eta5, msg.HFEMaxMinus_eta5, msg.nVtx, msg.nTrackInAcceptanceHP, msg.Run, msg.Lumi);

  }
  __XJJLOG << ">> " << ntuple->GetEntries() << " / " << nentries << " entries passed selection, written to " << outputname << std::endl;
  if (ntuple->GetEntries() == 0) __XJJLOG << "!! warning: no entries passed the selection." << std::endl;

  // make output file, save histograms, etc.
  ntuple->Write("", TObject::kOverwrite);
  outf->Close();
  inf->Close();

  return 0;
}

static void usage(const char* prog) {
  std::cerr << "usage: " << prog << " <inputname.root (output of addlumi.cc)> <outputname.root> [triggerSelection=0] [goodlumi.json]\n"
            << "  triggerSelection:\n"
            << "     0  any ZB UPC path (default)\n"
            << "    -1  isNotBptxOR || isUnpairedBunchBptxMinus || isUnpairedBunchBptxPlus\n"
            << "    -2  isUnpairedBunchBptxMinus\n"
            << "    -3  isUnpairedBunchBptxPlus\n"
            << "    -4  isNotBptxOR\n"
            << "  goodlumi.json (optional): CMS certification json; if given, only events whose (run, lumi section)\n"
            << "    fall in the good ranges are kept\n";
}

int main(int argc, char* argv[]) {
  if (argc == 2 && (std::string(argv[1]) == "-h" || std::string(argv[1]) == "--help")) { usage(argv[0]); return 0; }
  if (argc >= 3 && argc <= 5) {
    int triggerSelection = 0;
    if (argc >= 4) {
      char* end = nullptr;
      long v = std::strtol(argv[3], &end, 10);
      if (end == argv[3] || *end != '\0') {
        __XJJLOG << "!! triggerSelection must be an integer, got '" << argv[3] << "'" << std::endl;
        usage(argv[0]);
        return 1;
      }
      triggerSelection = (int)v;
    }
    return macro(argv[1], argv[2], triggerSelection, argc == 5 ? argv[4] : "");
  }
  usage(argv[0]);
  return 1;
}
