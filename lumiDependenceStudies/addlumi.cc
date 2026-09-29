#include <fstream>
#include <filesystem>
#include <unordered_map>
#include <unordered_set>

#include "xjjanauti.h"

namespace {
  // init/end/avg are instantaneous luminosity in units of 10^33 cm^-2 s^-1 (OMS lumisections meta), from instlumi.py
  // leveled is true when instlumi.py's lumi_trend column classified the lumisection as "leveled"
  struct lumival { double init = 0.; double end = 0.; double avg = 0.; bool leveled = false; };

  double lumitoken_to_double(const std::string& tok) {
    if (tok.empty() || tok == "null") return 0.; // mystr() in omstools prints falsy (== 0.) values as "null"
    return std::atof(tok.c_str());
  }

  bool load_instlumi(const std::string& csvname, std::unordered_map<long long, lumival>& lumimap) {
    std::ifstream csv(csvname);
    if (!csv.is_open()) { __XJJLOG << "!! failed to open inst-lumi csv, abort." << std::endl; return false; }
    std::string line;
    bool is_header = true;
    while (std::getline(csv, line)) {
      if (line.empty()) continue;
      if (is_header) { is_header = false; continue; }
      const auto tokens = xjjc::str_divide_trim(line, ",");
      if (tokens.size() < 9) continue;
      const auto run = std::atoi(tokens[0].c_str());
      const auto ls = std::atoi(tokens[1].c_str());
      lumival v;
      v.init = lumitoken_to_double(tokens[2]);
      v.end = lumitoken_to_double(tokens[3]);
      v.avg = lumitoken_to_double(tokens[4]);
      v.leveled = (tokens[8] == "leveled");
      lumimap[(long long)run * 10000000LL + ls] = v;
    }
    return true;
  }
}

int macro(const std::string& inputname, const std::string& instlumicsv, const std::string& outputname, int overwrite = 0) {
  namespace fs = std::filesystem;

  // guard against ever touching the input: output must be a distinct path,
  // and the copy below is the only thing ever opened for writing.
  if (fs::weakly_canonical(inputname) == fs::weakly_canonical(outputname)) {
    __XJJLOG << "!! output must differ from input, abort." << std::endl; return 1;
  }
  if (!overwrite && fs::exists(outputname)) {
    __XJJLOG << "!! output already exists, abort (pass overwrite = 1 to replace)." << std::endl; return 1;
  }

  std::unordered_map<long long, lumival> lumimap;
  if (!load_instlumi(instlumicsv, lumimap)) return 2;
  __XJJLOG << ">> " << lumimap.size() << " run-lumisection entries loaded from " << instlumicsv << std::endl;

  __XJJLOG << ">> copying " << inputname << " -> " << outputname << " ..." << std::endl;
  std::error_code ec;
  fs::copy_file(inputname, outputname, fs::copy_options::overwrite_existing, ec);
  if (ec) { __XJJLOG << "!! copy failed: " << ec.message() << ", abort." << std::endl; return 3; }
  __XJJLOG << ">> copy done." << std::endl;

  auto* outf = TFile::Open(outputname.c_str(), "UPDATE");
  if (!outf || outf->IsZombie()) { __XJJLOG << "!! failed to reopen copied file, abort." << std::endl; return 4; }
  auto* trs = (TTree*)outf->Get("Tree");
  if (!trs) { __XJJLOG << "!! no Tree found in copied file, abort." << std::endl; return 5; }

  trs->SetBranchStatus("*", 0);
  trs->SetBranchStatus("Run", 1);
  trs->SetBranchStatus("Lumi", 1);
  int run, lumi;
  trs->SetBranchAddress("Run", &run);
  trs->SetBranchAddress("Lumi", &lumi);

  double instlumi_init = 0., instlumi_end = 0., instlumi_avg = 0.; // all in 10^33 cm^-2 s^-1
  auto* br_init = trs->Branch("instLumiInit", &instlumi_init, "instLumiInit/D");
  auto* br_end = trs->Branch("instLumiEnd", &instlumi_end, "instLumiEnd/D");
  auto* br_avg = trs->Branch("instLumi", &instlumi_avg, "instLumi/D");
  bool instlumi_leveled = false; // true if instlumi.py classified this lumisection's lumi trend as "leveled"
  auto* br_leveled = trs->Branch("instLumiLeveled", &instlumi_leveled, "instLumiLeveled/O");

  const auto nentries = trs->GetEntries();
  Long64_t nmissing = 0;
  std::unordered_set<long long> warned_runlumi; // warn once per distinct missing (run, lumi), not once per event
  for (Long64_t i = 0; i < nentries; i++) {
    trs->GetEntry(i);
    const auto key = (long long)run * 10000000LL + lumi;
    const auto it = lumimap.find(key);
    if (it != lumimap.end()) {
      instlumi_init = it->second.init;
      instlumi_end = it->second.end;
      instlumi_avg = it->second.avg;
      instlumi_leveled = it->second.leveled;
    } else {
      instlumi_init = instlumi_end = instlumi_avg = -1.; // no match in instlumicsv
      instlumi_leveled = false;
      nmissing++;
      if (warned_runlumi.insert(key).second) {
        __XJJLOG << "?? event " << i << " : run " << run << ", lumi " << lumi << " not found in " << instlumicsv << ", instLumi set to -1." << std::endl;
      }
    }
    br_init->Fill();
    br_end->Fill();
    br_avg->Fill();
    br_leveled->Fill();
    if (i % 2000000 == 0) __XJJLOG << ">> " << i << "/" << nentries << " entries processed..." << std::endl;
  }
  if (nmissing) { __XJJLOG << "?? " << nmissing << "/" << nentries << " entries had no matching run-lumisection in " << instlumicsv << std::endl; }

  trs->SetBranchStatus("*", 1);
  trs->Write("", TObject::kOverwrite);
  outf->Close();

  __XJJLOG << ">> written " << outputname << std::endl;

  return 0;
}

int main(int argc, char* argv[]) {
  if (argc == 4) {
    return macro(argv[1], argv[2], argv[3]);
  }
  if (argc == 5) {
    return macro(argv[1], argv[2], argv[3], atoi(argv[4]));
  }
  __XJJLOG << "!! usage: " << argv[0] << " <inputname.root> <instlumicsv> <outputname.root> [overwrite = 0]" << std::endl;
  return 1;
}
