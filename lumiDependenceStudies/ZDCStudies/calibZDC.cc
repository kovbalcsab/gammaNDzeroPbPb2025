#include <filesystem>
#include <algorithm>
#include <array>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>

#include <TF1.h>
#include <TGraphErrors.h>

#include "xjjanauti.h"
#include "../params.h"

// Derives the lumi-dependent ZDC recalibration (see zdcLumiCorrection.h) from the UNCORRECTED ntuple.
// Per lumi bin and per ZDC side it measures
//   - noise core:  q90, q95 of the 0n-side sample (ZDCsum < ZDC_0N_CUT), 5 ADC bins
//   - 1n peak:     Gaussian mean m1 and sigma s1 (An side, other side < ZDC_0N_CUT)
//   - 2n peak:     Gaussian mean m2 (gaus + expo bkg)
// then fits each knot position vs lumi with a pol2 and writes
//   outputs/calib_TAG/Calib.root, outputs/calib_TAG/calib_<side>.pdf and zdcLumiCorrParams.h.
// The reference (what every lumi is mapped to) is the pol2 evaluated at the mean lumi of lumibin 0.
// Only non-leveled events (lumiLeveled == 0) are used.

namespace {
  const int NBINS_ZDC = 200;
  const double ZDC_MIN = 0., ZDC_MAX = 10000.;
  const double ZDC_0N_CUT = 1500.;
  const int NBINS_NOISE = 300; // 0 - ZDC_0N_CUT, 5 ADC bins
  const double LUNIT = 1.e-6;
  const int NK = 6; // knots: n90, n95, m1-s1, m1, m1+s1, m2
  const std::array<std::string, NK> knotNames = { "n90", "n95", "m1lo", "m1", "m1hi", "m2" };
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
    const double t = h->Integral();
    double c = 0.;
    for (int b = 1; b <= h->GetNbinsX(); b++) {
      const double pc = c;
      c += h->GetBinContent(b) / t;
      if (c >= p) return h->GetBinLowEdge(b) + h->GetBinWidth(b) * (p - pc) / (c - pc);
    }
    return h->GetXaxis()->GetXmax();
  }
}

int macro(const std::string& inputname, const std::string& tag) {
  namespace fs = std::filesystem;
  const std::string outdir = "outputs/calib_" + tag;
  fs::create_directories(outdir);

  auto* inf = TFile::Open(inputname.c_str());
  if (!inf || inf->IsZombie()) { __XJJLOG << "!! failed to open input file, abort." << std::endl; return 1; }
  auto* nt = (TTree*)inf->Get("ntuple");
  if (!nt) { __XJJLOG << "!! no ntuple found in input file, abort." << std::endl; return 2; }

  Float_t instLumi = 0., ZDCsumPlus = 0., ZDCsumMinus = 0., lumiLeveled = 0.;
  nt->SetBranchStatus("*", 0);
  for (const char* b : { "instLumi", "ZDCsumPlus", "ZDCsumMinus", "lumiLeveled" }) nt->SetBranchStatus(b, 1);
  nt->SetBranchAddress("instLumi", &instLumi);
  nt->SetBranchAddress("ZDCsumPlus", &ZDCsumPlus);
  nt->SetBranchAddress("ZDCsumMinus", &ZDCsumMinus);
  nt->SetBranchAddress("lumiLeveled", &lumiLeveled);

  const int nb = (int)params::lumibins.size() - 1;
  // hAn[s][i]: ZDC spectrum of side s when the OTHER side is 0n.  hNoise[s][i]: side s when side s itself is 0n.
  std::array<std::vector<TH1D*>, 2> hAn, hNoise;
  for (int s = 0; s < 2; s++)
    for (int i = 0; i < nb; i++) {
      hAn[s].push_back(new TH1D(Form("hAn_%s_%d", sides[s].c_str(), i), "", NBINS_ZDC, ZDC_MIN, ZDC_MAX));
      hNoise[s].push_back(new TH1D(Form("hNoise_%s_%d", sides[s].c_str(), i), "", NBINS_NOISE, 0., ZDC_0N_CUT));
    }
  std::vector<double> sumLumi(nb, 0.), nLumi(nb, 0.);

  const auto n = nt->GetEntries();
  for (Long64_t e = 0; e < n; e++) {
    nt->GetEntry(e);
    if (e % 10000000 == 0) __XJJLOG << ">> processing entry " << e << " / " << n << std::endl;
    if ((int)lumiLeveled != 0) continue;
    const int i = find_lumibin(instLumi);
    if (i < 0) continue;
    sumLumi[i] += instLumi; nLumi[i] += 1.;
    if (ZDCsumPlus < ZDC_0N_CUT) { hNoise[0][i]->Fill(ZDCsumPlus); hAn[1][i]->Fill(ZDCsumMinus); } // Plus is 0n
    if (ZDCsumMinus < ZDC_0N_CUT) { hNoise[1][i]->Fill(ZDCsumMinus); hAn[0][i]->Fill(ZDCsumPlus); } // Minus is 0n
  }

  std::vector<double> lumiMean(nb);
  for (int i = 0; i < nb; i++) lumiMean[i] = sumLumi[i] / nLumi[i];

  // per side, per knot: measured value and error per lumi bin
  std::array<std::array<std::vector<double>, NK>, 2> val, err;
  std::array<std::vector<double>, 2> sig1, sig1Err;
  for (int s = 0; s < 2; s++) {
    for (int k = 0; k < NK; k++) { val[s][k].assign(nb, 0.); err[s][k].assign(nb, 0.); }
    sig1[s].assign(nb, 0.); sig1Err[s].assign(nb, 0.);
    for (int i = 0; i < nb; i++) {
      auto* h = hAn[s][i];
      // 1n peak: gaussian on 1500-3500, refit +-1000 around the mean (same recipe as plotAndFit.cc)
      TF1 g1("g1", "gaus", 1500., 3500.);
      h->Fit(&g1, "RLQN");
      TF1 g1b("g1b", "gaus", g1.GetParameter(1) - 1000., g1.GetParameter(1) + 1000.);
      h->Fit(&g1b, "RLQN");
      const double m1 = g1b.GetParameter(1), s1 = std::fabs(g1b.GetParameter(2));
      // 2n peak: gaus + exp bkg
      const double lo = 1.65 * m1, hi = 2.45 * m1;
      double ymax = 0.;
      for (int b = h->FindBin(lo); b <= h->FindBin(hi); b++) ymax = std::max(ymax, h->GetBinContent(b));
      TF1 g2("g2", "gaus(0)+expo(3)", lo, hi);
      g2.SetParameters(ymax * 0.5, 2. * m1, 700., std::log(ymax * 0.5), 0.);
      g2.SetParLimits(1, 1.7 * m1, 2.4 * m1);
      g2.SetParLimits(2, 300., 1500.);
      h->Fit(&g2, "RQN");
      // noise
      auto* hn = hNoise[s][i];
      val[s][0][i] = quantile(hn, 0.90); err[s][0][i] = 1.;
      val[s][1][i] = quantile(hn, 0.95); err[s][1][i] = 1.;
      val[s][2][i] = m1 - s1; err[s][2][i] = std::hypot(g1b.GetParError(1), g1b.GetParError(2));
      val[s][3][i] = m1;      err[s][3][i] = g1b.GetParError(1);
      val[s][4][i] = m1 + s1; err[s][4][i] = std::hypot(g1b.GetParError(1), g1b.GetParError(2));
      val[s][5][i] = g2.GetParameter(1); err[s][5][i] = g2.GetParError(1);
      sig1[s][i] = s1; sig1Err[s][i] = g1b.GetParError(2);
      __XJJLOG << ">> " << sides[s] << " bin " << i << " <L>=" << lumiMean[i] / LUNIT << "e-6  n90=" << val[s][0][i] << " n95=" << val[s][1][i]
               << " m1=" << m1 << " s1=" << s1 << " m2=" << val[s][5][i] << " (m2/m1=" << val[s][5][i] / m1 << ")" << std::endl;
    }
  }

  // fit every knot vs u = lumi / LUNIT with a pol2
  std::array<std::array<TF1*, NK>, 2> fits;
  std::array<std::array<TGraphErrors*, NK>, 2> graphs;
  std::vector<double> u(nb), uErr(nb, 0.);
  for (int i = 0; i < nb; i++) u[i] = lumiMean[i] / LUNIT;

  std::ofstream hdr("zdcLumiCorrParams.h");
  hdr << "// generated by calibZDC.cc (tag " << tag << ", non-leveled events); do not edit by hand\n"
      << "// knot input position = c0 + c1*u + c2*u^2 with u = instLumi / kLumiUnit; reference = same pol2 at u = kURef (lumibin 0)\n"
      << "#pragma once\nnamespace zdccorr {\n"
      << "  const double kLumiUnit = " << LUNIT << ";\n  const int kNKnots = " << NK << ";\n"
      << "  const double kULow = " << u.front() << ", kUHigh = " << u.back() << "; // instLumi is clamped to this range (in units of kLumiUnit)\n"
      << "  const double kURef = " << u.front() << ";\n"
      << "  // [side: 0=Plus 1=Minus][knot: n90 n95 m1-s1 m1 m1+s1 m2][c0 c1 c2]\n"
      << "  const double kKnotPoly[2][" << NK << "][3] = {\n";
  std::ostringstream refs;
  for (int s = 0; s < 2; s++) {
    hdr << "    {\n";
    for (int k = 0; k < NK; k++) {
      auto* g = new TGraphErrors(nb, u.data(), val[s][k].data(), uErr.data(), err[s][k].data());
      g->SetName(Form("g_%s_%s", knotNames[k].c_str(), sides[s].c_str()));
      auto* f = new TF1(Form("f_%s_%s", knotNames[k].c_str(), sides[s].c_str()), "pol2", u.front(), u.back());
      g->Fit(f, "QN");
      graphs[s][k] = g; fits[s][k] = f;
      hdr << "      { " << std::setprecision(10) << f->GetParameter(0) << ", " << f->GetParameter(1) << ", " << f->GetParameter(2) << " }" << (k < NK - 1 ? "," : "") << "  // " << knotNames[k] << "\n";
      double maxRes = 0.;
      for (int i = 0; i < nb; i++) maxRes = std::max(maxRes, std::fabs(val[s][k][i] - f->Eval(u[i])));
      __XJJLOG << ">> fit " << sides[s] << " " << knotNames[k] << ": ref(u=" << u.front() << ")=" << f->Eval(u.front()) << "  max |residual|=" << maxRes << std::endl;
    }
    hdr << "    }" << (s == 0 ? "," : "") << "\n";
  }
  hdr << "  };\n};\n";
  hdr.close();

  // diagnostics: graphs + fits
  auto* outf = TFile::Open((outdir + "/Calib.root").c_str(), "RECREATE");
  xjjroot::setgstyle(1);
  for (int s = 0; s < 2; s++) {
    auto* c = new TCanvas(Form("c_%d", s), "", 1500, 800);
    c->Divide(3, 2);
    for (int k = 0; k < NK; k++) {
      c->cd(k + 1);
      graphs[s][k]->SetTitle(Form("%s %s;instLumi / 1e-6;ADC", sides[s].c_str(), knotNames[k].c_str()));
      graphs[s][k]->SetMarkerStyle(20);
      graphs[s][k]->Draw("ape");
      fits[s][k]->SetLineColor(kRed + 1);
      fits[s][k]->Draw("same");
      graphs[s][k]->Write();
      fits[s][k]->Write();
    }
    c->SaveAs((outdir + "/calib_" + sides[s] + ".pdf").c_str());
    auto* gs = new TGraphErrors(nb, u.data(), sig1[s].data(), uErr.data(), sig1Err[s].data());
    gs->SetName(Form("g_s1_%s", sides[s].c_str()));
    gs->Write();
  }
  for (int s = 0; s < 2; s++)
    for (int i = 0; i < nb; i++) { hAn[s][i]->Write(); hNoise[s][i]->Write(); }
  outf->Close();

  __XJJLOG << ">> wrote zdcLumiCorrParams.h and " << outdir << "/ ; copy/keep the header next to zdcLumiCorrection.h" << std::endl;
  inf->Close();
  return 0;
}

int main(int argc, char* argv[]) {
  if (argc == 3) return macro(argv[1], argv[2]);
  __XJJLOG << "!! usage: " << argv[0] << " <inputname.root (ntuple from ../analyze.cc, UNCORRECTED)> <TAG (-> outputs/calib_TAG/, ./zdcLumiCorrParams.h)>" << std::endl;
  return 1;
}
