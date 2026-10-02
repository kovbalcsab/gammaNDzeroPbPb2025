#include <filesystem>
#include <algorithm>
#include <vector>
#include <array>
#include <fstream>
#include <iomanip>
#include <cmath>

#include <TF1.h>
#include <TGraphErrors.h>

#include "xjjanauti.h"
#include "../params.h"

// Derives the lumi-dependent HFEMax recalibration (see hfLumiCorrection.h) from the UNCORRECTED ntuple.
// Only the 0n side is calibrated:  HFEMaxPlus_eta5 when ZDCsumPlus < ZDC_0N_CUT (RAW sum),
//                                  HFEMaxMinus_eta5 when ZDCsumMinus < ZDC_0N_CUT.
// Per lumi bin and per side it measures the quantile positions x_p(L) of a fixed grid of cumulative levels p
// (HF_MAX_KNOT is the upper edge of the calibrated region), fits each x_p vs <instLumi>/1e-6 with a polynomial
// (pol2, or pol1 with polyOrder = 1) and writes
//   outputs/calib_TAG/Calib.root, calib_<side>.pdf (selected levels), calib_all_<side>.pdf (all levels, multi-page) and hfLumiCorrParams.h.
// The reference (what every lumi is mapped to) is the fit evaluated at the mean lumi of lumibin refBin; lumi bins below refBin are
// not used in the fits (the lowest-lumi bins have an anomalous high tail), and instLumi below the first used bin is clamped to it.
// Only non-leveled events (lumiLeveled == 0) are used.

namespace {
  const double ZDC_0N_CUT = 1500.;
  const int NBINS_HF = 400; // 0.25 GeV bins
  const double HF_MIN = 0., HF_MAX = 100.;
  const double HF_MAX_KNOT = 55.;    // highest quantile knot position (any lumi); the map is the identity from HF_END
  const double HF_END = 60.;
  const double LUNIT = 1.e-6;
  const std::array<std::string, 2> sides = { "Plus", "Minus" };

  int find_lumibin(double lumi) {
    const auto& edges = params::lumibins;
    const int nbins = (int)edges.size() - 1;
    if (lumi < edges.front() || lumi > edges.back()) return -1;
    auto it = std::upper_bound(edges.begin(), edges.end(), lumi);
    return std::min(std::max(int(it - edges.begin()) - 1, 0), nbins - 1);
  }

  // x at which the cumulative fraction of h reaches p (linear interpolation inside the bin)
  double quantile(const TH1D* h, double p) {
    const double t = h->Integral(0, h->GetNbinsX() + 1);
    double c = h->GetBinContent(0) / t;
    if (c >= p) return h->GetXaxis()->GetXmin();
    for (int b = 1; b <= h->GetNbinsX(); b++) {
      const double pc = c;
      c += h->GetBinContent(b) / t;
      if (c >= p) return h->GetBinLowEdge(b) + h->GetBinWidth(b) * (p - pc) / (c - pc);
    }
    return h->GetXaxis()->GetXmax();
  }

  // cumulative fraction at x (linear interpolation inside the bin)
  double cdf(const TH1D* h, double x) {
    const double t = h->Integral(0, h->GetNbinsX() + 1);
    const int b = h->GetXaxis()->FindBin(x);
    double c = h->Integral(0, b - 1);
    c += h->GetBinContent(b) * (x - h->GetBinLowEdge(b)) / h->GetBinWidth(b);
    return c / t;
  }

  // candidate cumulative levels: 0.02 steps to 0.90, 0.01 steps to 0.99, 0.002 steps to 0.998
  std::vector<double> candidate_levels() {
    std::vector<double> p;
    for (int i = 1; i <= 45; i++) p.push_back(0.02 * i);
    for (int i = 91; i <= 99; i++) p.push_back(0.01 * i);
    for (int i = 992; i <= 998; i += 2) p.push_back(0.001 * i);
    return p;
  }
}

// parity: 0 = all events, 1 = even entries only, 2 = odd entries only (for split-sample closure tests)
int macro(const std::string& inputname, const std::string& tag, int refBin = 0, int polyOrder = 2, int parity = 0) {
  namespace fs = std::filesystem;
  const std::string outdir = "outputs/calib_" + tag;
  fs::create_directories(outdir);

  auto* inf = TFile::Open(inputname.c_str());
  if (!inf || inf->IsZombie()) { __XJJLOG << "!! failed to open input file, abort." << std::endl; return 1; }
  auto* nt = (TTree*)inf->Get("ntuple");
  if (!nt) { __XJJLOG << "!! no ntuple found in input file, abort." << std::endl; return 2; }

  Float_t instLumi = 0., ZDCsumPlus = 0., ZDCsumMinus = 0., HFEMaxPlus = 0., HFEMaxMinus = 0., lumiLeveled = 0.;
  nt->SetBranchStatus("*", 0);
  for (const char* b : { "instLumi", "ZDCsumPlus", "ZDCsumMinus", "HFEMaxPlus_eta5", "HFEMaxMinus_eta5", "lumiLeveled" }) nt->SetBranchStatus(b, 1);
  nt->SetBranchAddress("instLumi", &instLumi);
  nt->SetBranchAddress("ZDCsumPlus", &ZDCsumPlus);
  nt->SetBranchAddress("ZDCsumMinus", &ZDCsumMinus);
  nt->SetBranchAddress("HFEMaxPlus_eta5", &HFEMaxPlus);
  nt->SetBranchAddress("HFEMaxMinus_eta5", &HFEMaxMinus);
  nt->SetBranchAddress("lumiLeveled", &lumiLeveled);

  const int nb = (int)params::lumibins.size() - 1;
  // h[s][i]: HFEMax of side s (0n side only) in lumi bin i
  std::array<std::vector<TH1D*>, 2> h;
  for (int s = 0; s < 2; s++)
    for (int i = 0; i < nb; i++) {
      h[s].push_back(new TH1D(Form("hHF0n_%s_%d", sides[s].c_str(), i), ";HFEMax;Events", NBINS_HF, HF_MIN, HF_MAX));
      h[s].back()->Sumw2();
    }
  std::array<std::vector<double>, 2> sumLumi, nLumi;
  for (int s = 0; s < 2; s++) { sumLumi[s].assign(nb, 0.); nLumi[s].assign(nb, 0.); }

  const auto n = nt->GetEntries();
  for (Long64_t e = 0; e < n; e++) {
    if (parity == 1 && (e % 2) != 0) continue;
    if (parity == 2 && (e % 2) != 1) continue;
    nt->GetEntry(e);
    if (e % 10000000 == 0) __XJJLOG << ">> processing entry " << e << " / " << n << std::endl;
    if ((int)lumiLeveled != 0) continue;
    const int i = find_lumibin(instLumi);
    if (i < 0) continue;
    if (ZDCsumPlus < ZDC_0N_CUT)  { h[0][i]->Fill(HFEMaxPlus);  sumLumi[0][i] += instLumi; nLumi[0][i] += 1.; } // Plus is 0n
    if (ZDCsumMinus < ZDC_0N_CUT) { h[1][i]->Fill(HFEMaxMinus); sumLumi[1][i] += instLumi; nLumi[1][i] += 1.; } // Minus is 0n
  }

  // lumi bins used in the fits, and their mean lumi in units of LUNIT (per side: the 0n samples differ slightly)
  std::array<std::vector<int>, 2> used;
  std::array<std::vector<double>, 2> u;
  for (int s = 0; s < 2; s++) {
    u[s].assign(nb, 0.);
    for (int i = 0; i < nb; i++) {
      if (nLumi[s][i] > 0) u[s][i] = sumLumi[s][i] / nLumi[s][i] / LUNIT;
      if (nLumi[s][i] > 0 && i >= refBin) used[s].push_back(i);
      __XJJLOG << ">> " << sides[s] << " bin " << i << " <L>=" << u[s][i] << "e-6  n(0n)=" << nLumi[s][i]
               << (i < refBin ? "  (below refBin, not used)" : nLumi[s][i] > 0 ? "" : "  (empty, not used)") << std::endl;
    }
    if (std::find(used[s].begin(), used[s].end(), refBin) == used[s].end()) {
      __XJJLOG << "!! reference bin " << refBin << " has no events for side " << sides[s] << ", abort." << std::endl;
      return 4;
    }
  }
  // keep the levels whose position stays below HF_MAX_KNOT in every used bin of both sides
  std::vector<double> levels;
  for (double p : candidate_levels()) {
    bool ok = true;
    for (int s = 0; s < 2; s++)
      for (int i : used[s]) ok = ok && quantile(h[s][i], p) < HF_MAX_KNOT;
    if (ok) levels.push_back(p);
  }
  const int NK = (int)levels.size();
  __XJJLOG << ">> " << NK << " quantile levels, " << levels.front() << " ... " << levels.back() << std::endl;

  // x_p(L_i) and its statistical error  sqrt(p(1-p)/n) / density(x_p)
  std::array<std::vector<std::vector<double>>, 2> val, err; // [s][k][used index]
  std::array<std::vector<double>, 2> uu, uErr;
  for (int s = 0; s < 2; s++) {
    val[s].assign(NK, {}); err[s].assign(NK, {});
    for (int i : used[s]) {
      uu[s].push_back(u[s][i]); uErr[s].push_back(0.);
      for (int k = 0; k < NK; k++) {
        const double x = quantile(h[s][i], levels[k]);
        const double d = 1.; // GeV, half window of the density estimate
        const double dens = std::max((cdf(h[s][i], x + d) - cdf(h[s][i], std::max(x - d, 0.))) / (x + d - std::max(x - d, 0.)), 1.e-6);
        val[s][k].push_back(x);
        err[s][k].push_back(std::sqrt(levels[k] * (1. - levels[k]) / nLumi[s][i]) / dens);
      }
    }
  }

  // fit every level vs u with pol<polyOrder>
  std::array<std::vector<TF1*>, 2> fits;
  std::array<std::vector<TGraphErrors*>, 2> graphs;
  std::ofstream hdr("hfLumiCorrParams.h");
  hdr << "// generated by calibHF.cc (tag " << tag << ", non-leveled events, parity " << parity << ", reference lumibin " << refBin << "); do not edit by hand\n"
      << "// quantile level k of the 0n-side HFEMax spectrum sits at c0 + c1*u + c2*u^2 with u = instLumi / kLumiUnit; reference = same poly at u = kURef[side]\n"
      << "#pragma once\nnamespace hfcorr {\n"
      << "  const double kLumiUnit = " << LUNIT << ";\n  const int kNKnots = " << NK << ";\n"
      << "  const double kHFEnd = " << HF_END << "; // identity from here on, the last quantile knot is joined to (kHFEnd, kHFEnd)\n"
      << "  const double kULow[2] = { " << u[0][used[0].front()] << ", " << u[1][used[1].front()] << " }, kUHigh[2] = { " << u[0][used[0].back()] << ", " << u[1][used[1].back()] << " }; // instLumi is clamped to this range (in units of kLumiUnit)\n";
  hdr << "  const double kURef[2] = { " << u[0][refBin] << ", " << u[1][refBin] << " };\n";
  hdr << "  const double kLevels[" << NK << "] = { ";
  for (int k = 0; k < NK; k++) hdr << levels[k] << (k < NK - 1 ? ", " : " };\n");
  hdr << "  // [side: 0=Plus 1=Minus][level k][c0 c1 c2]\n  const double kKnotPoly[2][" << NK << "][3] = {\n";
  for (int s = 0; s < 2; s++) {
    hdr << "    {\n";
    for (int k = 0; k < NK; k++) {
      auto* g = new TGraphErrors((int)uu[s].size(), uu[s].data(), val[s][k].data(), uErr[s].data(), err[s][k].data());
      g->SetName(Form("g_p%03d_%s", (int)std::lround(levels[k] * 1000.), sides[s].c_str()));
      auto* f = new TF1(Form("f_p%03d_%s", (int)std::lround(levels[k] * 1000.), sides[s].c_str()), polyOrder == 1 ? "pol1" : "pol2", uu[s].front(), uu[s].back());
      g->Fit(f, "QN");
      graphs[s].push_back(g); fits[s].push_back(f);
      hdr << "      { " << std::setprecision(10) << f->GetParameter(0) << ", " << f->GetParameter(1) << ", " << (polyOrder == 1 ? 0. : f->GetParameter(2)) << " }"
          << (k < NK - 1 ? "," : "") << "  // p=" << levels[k] << "\n";
      double maxRes = 0.;
      for (size_t j = 0; j < uu[s].size(); j++) maxRes = std::max(maxRes, std::fabs(val[s][k][j] - f->Eval(uu[s][j])));
      if (k % 5 == 0 || k == NK - 1)
        __XJJLOG << ">> fit " << sides[s] << " p=" << levels[k] << ": ref=" << f->Eval(u[s][refBin]) << "  max |residual|=" << maxRes << " GeV" << std::endl;
    }
    hdr << "    }" << (s == 0 ? "," : "") << "\n";
  }
  hdr << "  };\n};\n";
  hdr.close();

  // diagnostics: x_p vs lumi for a selection of levels, plus all graphs/fits/histograms in Calib.root
  auto* outf = TFile::Open((outdir + "/Calib.root").c_str(), "RECREATE");
  xjjroot::setgstyle(1);
  const std::vector<double> showLevels = { 0.1, 0.3, 0.5, 0.7, 0.9, 0.96, 0.98, 0.99 };
  for (int s = 0; s < 2; s++) {
    auto* c = new TCanvas(Form("c_%d", s), "", 1600, 800);
    c->Divide(4, 2);
    int pad = 1;
    for (double ps : showLevels) {
      int kk = 0;
      for (int k = 0; k < NK; k++) if (std::fabs(levels[k] - ps) < std::fabs(levels[kk] - ps)) kk = k;
      c->cd(pad++);
      graphs[s][kk]->SetTitle(Form("%s p=%g;instLumi / 1e-6;HFEMax (GeV)", sides[s].c_str(), levels[kk]));
      graphs[s][kk]->SetMarkerStyle(20);
      graphs[s][kk]->Draw("ape");
      fits[s][kk]->SetLineColor(kRed + 1);
      fits[s][kk]->Draw("same");
    }
    c->SaveAs((outdir + "/calib_" + sides[s] + ".pdf").c_str());

    // every quantile level: x_p vs lumi with its fit, 12 per page, multi-page calib_all_<side>.pdf
    const int perPage = 12, nPages = (NK + perPage - 1) / perPage;
    const std::string allpdf = outdir + "/calib_all_" + sides[s] + ".pdf";
    auto* call = new TCanvas(Form("call_%d", s), "", 1600, 1200);
    for (int pg = 0; pg < nPages; pg++) {
      call->Clear();
      call->Divide(4, 3);
      for (int j = 0; j < perPage && pg * perPage + j < NK; j++) {
        const int k = pg * perPage + j;
        call->cd(j + 1);
        graphs[s][k]->SetTitle(Form("%s p=%g;instLumi / 1e-6;HFEMax (GeV)", sides[s].c_str(), levels[k]));
        graphs[s][k]->SetMarkerStyle(20);
        graphs[s][k]->Draw("ape");
        fits[s][k]->SetLineColor(kRed + 1);
        fits[s][k]->Draw("same");
      }
      call->Print((allpdf + (pg == 0 ? (nPages == 1 ? "" : "(") : pg == nPages - 1 ? ")" : "")).c_str());
    }
    delete call;
    for (int k = 0; k < NK; k++) { graphs[s][k]->Write(); fits[s][k]->Write(); }
    for (int i = 0; i < nb; i++) h[s][i]->Write();
  }
  outf->Close();

  __XJJLOG << ">> wrote hfLumiCorrParams.h and " << outdir << "/ ; keep the header next to hfLumiCorrection.h" << std::endl;
  inf->Close();
  return 0;
}

int main(int argc, char* argv[]) {
  if (argc >= 3 && argc <= 6) {
    return macro(argv[1], argv[2], argc > 3 ? std::atoi(argv[3]) : 0, argc > 4 ? std::atoi(argv[4]) : 2, argc > 5 ? std::atoi(argv[5]) : 0);
  }
  __XJJLOG << "!! usage: " << argv[0] << " <inputname.root (ntuple from ../analyze.cc, UNCORRECTED)> <TAG (-> outputs/calib_TAG/, ./hfLumiCorrParams.h)> [refBin=0] [polyOrder=2|1] [parity=0 all|1 even|2 odd entries]" << std::endl;
  return 1;
}
