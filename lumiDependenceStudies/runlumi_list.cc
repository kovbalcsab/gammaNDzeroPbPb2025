#include <fstream>
#include <set>

#include "xjjanauti.h"

int macro(const std::string& inputname, const std::string& outputname) {
  auto* trs = xjjana::chain_files(xjjc::str_divide_trim(inputname, ","), "Tree");
  if (!trs) { __XJJLOG << "!! bad input file, abort." << std::endl; return 1; }

  trs->SetBranchStatus("*", 0);
  trs->SetBranchStatus("Run", 1);
  trs->SetBranchStatus("Lumi", 1);
  int run, lumi;
  trs->SetBranchAddress("Run", &run);
  trs->SetBranchAddress("Lumi", &lumi);

  std::set<std::pair<int, int>> runlumis;
  const auto nentries = trs->GetEntries();
  for (Long64_t i = 0; i < nentries; i++) {
    trs->GetEntry(i);
    runlumis.insert({run, lumi});
  }
  __XJJLOG << ">> " << nentries << " entries scanned, " << runlumis.size() << " unique run-lumi pairs found." << std::endl;

  std::ofstream outf(outputname);
  if (!outf.is_open()) { __XJJLOG << "!! failed to open output file, abort." << std::endl; return 2; }
  for (const auto& [run_number, lumi_number] : runlumis) {
    outf << run_number << " " << lumi_number << std::endl;
  }
  outf.close();
  __XJJLOG << ">> written to " << outputname << std::endl;

  return 0;
}

int main(int argc, char* argv[]) {
  if (argc == 3) {
    return macro(argv[1], argv[2]);
  }
  __XJJLOG << "!! usage: " << argv[0] << " <inputname[,inputname2,...]> <outputname.txt>" << std::endl;
  return 1;
}
