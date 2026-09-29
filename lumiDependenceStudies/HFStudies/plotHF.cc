#include <filesystem>
#include <algorithm>
#include <array>
#include <map>
#include <cmath>

#include <TMultiGraph.h>

#include "xjjanauti.h"
#include "../params.h"

namespace {
  const double ZDC_0N_CUT = 1500.; // ZDCsum below this = "0n" (no neutron) on that side

  // HF-side binning: the 0n side has little forward activity (narrow range),
  // the An (neutron-unconstrained) side can have larger forward activity (wider range)
  const int NBINS_HF_0N = 35;
  const double HF_0N_MIN = 0., HF_0N_MAX = 35.;
  const int NBINS_HF_AN = 100;
  const double HF_AN_MIN = 0., HF_AN_MAX = 100.;

  // the two histogram keys whose side is unconstrained ("An"): the other side carries the 0n cut
  const std::array<std::string, 2> AN_SIDE_KEYS = { "Minus_An0n", "Plus_0nAn" };
  // the two histogram keys that are the 0n (cut) side: complement of AN_SIDE_KEYS
  const std::array<std::string, 2> ZN_SIDE_KEYS = { "Plus_An0n", "Minus_0nAn" };

  bool is_an_side(const std::string& key) {
    return std::find(AN_SIDE_KEYS.begin(), AN_SIDE_KEYS.end(), key) != AN_SIDE_KEYS.end();
  }

  const double HF_0N_FRAC_THRESHOLD = 16.; // fraction of 0n-side events with HFEMax above this value, vs lumi
  const std::array<Color_t, 2> side_colors = { kBlue + 1, kRed + 1 };
  const std::array<Style_t, 2> side_markers = { 20, 21 };

  // index of the params::lumibins bin containing lumi, or -1 if out of range
  int find_lumibin(double lumi) {
    const auto& edges = params::lumibins;
    const int nbins = (int)edges.size() - 1;
    if (lumi < edges.front() || lumi > edges.back()) return -1;
    auto it = std::upper_bound(edges.begin(), edges.end(), lumi);
    int idx = int(it - edges.begin()) - 1;
    if (idx < 0) idx = 0;
    if (idx > nbins - 1) idx = nbins - 1;
    return idx;
  }
}

int macro(const std::string& inputname, const std::string& tag, int leveled = -1) {
  // leveled: -1 = no cut on lumiLeveled (default), 0 = only non-leveled events, 1 = only leveled events
  namespace fs = std::filesystem;

  // all output (histogram file + plots) goes under outputs/TAG[_leveled|_notleveled]/
  const std::string outdir = "outputs/" + tag + (leveled == 1 ? "_leveled" : leveled == 0 ? "_notleveled" : "");
  fs::create_directories(outdir);
  const std::string outputname = outdir + "/Main.root";

  // get input file (ntuple produced by ../analyze.cc)
  auto* inf = TFile::Open(inputname.c_str());
  if (!inf || inf->IsZombie()) { __XJJLOG << "!! failed to open input file, abort." << std::endl; return 1; }
  auto* nt = (TTree*)inf->Get("ntuple");
  if (!nt) { __XJJLOG << "!! no ntuple found in input file, abort." << std::endl; return 2; }

  Float_t instLumi = 0., ZDCsumPlus = 0., ZDCsumMinus = 0., HFEMaxPlus = 0., HFEMaxMinus = 0., lumiLeveled = 0.;
  nt->SetBranchStatus("*", 0);
  nt->SetBranchStatus("instLumi", 1);
  nt->SetBranchStatus("ZDCsumPlus", 1);
  nt->SetBranchStatus("ZDCsumMinus", 1);
  nt->SetBranchStatus("HFEMaxPlus_eta5", 1);
  nt->SetBranchStatus("HFEMaxMinus_eta5", 1);
  nt->SetBranchAddress("instLumi", &instLumi);
  nt->SetBranchAddress("ZDCsumPlus", &ZDCsumPlus);
  nt->SetBranchAddress("ZDCsumMinus", &ZDCsumMinus);
  nt->SetBranchAddress("HFEMaxPlus_eta5", &HFEMaxPlus);
  nt->SetBranchAddress("HFEMaxMinus_eta5", &HFEMaxMinus);
  if (leveled >= 0) {
    nt->SetBranchStatus("lumiLeveled", 1);
    nt->SetBranchAddress("lumiLeveled", &lumiLeveled);
  }

  const auto nentries = nt->GetEntries();
  __XJJLOG << ">> " << nentries << " entries loaded from " << inputname << std::endl;

  const int nlumibins = (int)params::lumibins.size() - 1;

  // make 4 sets of nlumibins histograms for HF values, 2 per side, 2 per 0nAn or An0n cut, one per lumi bin
  const std::array<std::string, 2> sides = { "Plus", "Minus" };
  const std::array<std::string, 2> cuts  = { "An0n", "0nAn" };
  // An0n: ZDCsumPlus < 1500 (tags Plus as 0n; Minus is the unconstrained "An" side)
  // 0nAn: ZDCsumMinus < 1500 (tags Minus as 0n; Plus is the unconstrained "An" side)

  std::map<std::string, std::vector<TH1D*>> hists; // key = "<side>_<cut>"
  for (const auto& side : sides) {
    for (const auto& cut : cuts) {
      const std::string key = side + "_" + cut;
      const bool anSide = is_an_side(key);
      const int nbins = anSide ? NBINS_HF_AN : NBINS_HF_0N;
      const double hmin = anSide ? HF_AN_MIN : HF_0N_MIN;
      const double hmax = anSide ? HF_AN_MAX : HF_0N_MAX;
      auto& hs = hists[key];
      hs.reserve(nlumibins);
      for (int i = 0; i < nlumibins; i++) {
        auto* h = new TH1D(Form("h_HF%s_%s_lumibin%d", side.c_str(), cut.c_str(), i),
                            Form(";HFEMax%s_eta5;Events", side.c_str()),
                            nbins, hmin, hmax);
        h->Sumw2();
        hs.push_back(h);
      }
    }
  }

  // fill histograms with HF values for each lumi bin, for each side, for each 0nAn or An0n cut
  for (Long64_t i = 0; i < nentries; i++) {
    nt->GetEntry(i);
    if (i % 1000000 == 0) __XJJLOG << ">> processing entry " << i << " / " << nentries << std::endl;

    if (leveled >= 0 && (int)lumiLeveled != leveled) continue;

    const auto lumibin = find_lumibin(instLumi);
    if (lumibin < 0) continue;

    if (ZDCsumPlus < ZDC_0N_CUT) { // An0n
      hists.at("Plus_An0n")[lumibin]->Fill(HFEMaxPlus);
      hists.at("Minus_An0n")[lumibin]->Fill(HFEMaxMinus);
    }
    if (ZDCsumMinus < ZDC_0N_CUT) { // 0nAn
      hists.at("Plus_0nAn")[lumibin]->Fill(HFEMaxPlus);
      hists.at("Minus_0nAn")[lumibin]->Fill(HFEMaxMinus);
    }
  }

  // for each 0n-side histogram, compute the fraction of events with HFEMax above HF_0N_FRAC_THRESHOLD, per lumi bin
  // (binomial proportion error: sqrt(p*(1-p)/n))
  std::map<std::string, std::vector<double>> frac_above, frac_aboveErr; // keyed by side ("Plus"/"Minus")
  for (const auto& side : sides) {
    frac_above[side].assign(nlumibins, 0.);
    frac_aboveErr[side].assign(nlumibins, 0.);
  }
  for (const auto& key : ZN_SIDE_KEYS) {
    const std::string side = key.substr(0, key.find('_')); // "Plus" or "Minus"
    auto& hs = hists.at(key);
    for (int i = 0; i < nlumibins; i++) {
      auto* h = hs[i];
      const double n = h->Integral(0, h->GetNbinsX() + 1);
      const double nabove = h->Integral(h->GetXaxis()->FindBin(HF_0N_FRAC_THRESHOLD), h->GetNbinsX() + 1);
      const double p = (n > 0) ? nabove / n : 0.;
      frac_above[side][i] = p;
      frac_aboveErr[side][i] = (n > 0) ? std::sqrt(p * (1. - p) / n) : 0.;
    }
  }

  // build lumi-bin centers/half-widths (x, x-error) and the fraction-vs-lumi graphs
  std::vector<double> lumiCenter(nlumibins), lumiHalfWidth(nlumibins);
  for (int i = 0; i < nlumibins; i++) {
    lumiCenter[i] = 0.5 * (params::lumibins[i] + params::lumibins[i + 1]);
    lumiHalfWidth[i] = 0.5 * (params::lumibins[i + 1] - params::lumibins[i]);
  }
  std::map<std::string, TGraphErrors*> g_frac_above;
  for (const auto& side : sides) {
    auto* g = new TGraphErrors(nlumibins, lumiCenter.data(), frac_above.at(side).data(),
                              lumiHalfWidth.data(), frac_aboveErr.at(side).data());
    g->SetName(Form("g_frac_above%g_vs_lumi_%s", HF_0N_FRAC_THRESHOLD, side.c_str()));
    g->SetTitle(Form(";instLumi;Fraction HFEMax%s > %g (0n side)", side.c_str(), HF_0N_FRAC_THRESHOLD));
    g_frac_above[side] = g;
  }

  // save histograms and fraction-vs-lumi graphs to output file
  auto* outf = TFile::Open(outputname.c_str(), "RECREATE");
  if (!outf || outf->IsZombie()) { __XJJLOG << "!! failed to open output file, abort." << std::endl; return 3; }
  outf->cd();
  for (auto& [key, hs] : hists)
    for (auto* h : hs) h->Write();
  for (auto& [side, g] : g_frac_above) g->Write();
  __XJJLOG << ">> histograms and fraction-vs-lumi graphs written to " << outputname << std::endl;

  // plot each histogram under outputs/TAG/.....pdf
  xjjroot::setgstyle(1);
  auto* c = new TCanvas("c", "", 800, 600);
  c->SetLogy(1);
  for (auto& [key, hs] : hists) {
    for (auto* h : hs) {
      c->cd();
      h->Draw("hist e");
      c->SaveAs((outdir + "/" + h->GetName() + ".pdf").c_str());
    }
  }

  // for each side/role separately (never mix 0n and An sides, they use different binning ranges),
  // overlay that role's HF spectrum for all lumi bins on one canvas: black (lowest lumi) fading to
  // gray (highest lumi), same "black and gray systematic" style as main/hist_fit.cc.
  // each histogram is self-normalized (unit integral) first, so the shapes are comparable across
  // lumi bins regardless of differing statistics; originals in `hists` are left untouched.
  for (auto& [key, hs] : hists) {
    std::vector<TH1D*> hnorm;
    hnorm.reserve(hs.size());
    for (auto* h : hs) {
      auto* hn = (TH1D*)h->Clone(Form("%s_norm", h->GetName()));
      const double integral = hn->Integral();
      if (integral > 0) hn->Scale(1. / integral);
      hn->GetYaxis()->SetTitle("Fraction of events / bin");
      hnorm.push_back(hn);
    }

    xjjana::sethsmax(hnorm, 1.5); // common headroom across the series, log-y min stays automatic
    auto cc = xjjroot::grayscales_alpha((int)hnorm.size());
    auto* legAll = new TLegend(0.70, 0.86 - 0.025 * 1.25 * hnorm.size(), 0.88, 0.86);
    xjjroot::setleg(legAll, 0.025);
    for (int j = 0; j < (int)hnorm.size(); j++) {
      xjjroot::sethempty(hnorm[j], 0, 0.5);
      const auto mstyle = xjjroot::markerlist_solid[j % xjjroot::markerlist_solid.size()];
      xjjroot::setthgrstyle(hnorm[j], kBlack, mstyle, 1.0, kBlack, 1, 1, -1, -1, -1, cc[j], cc[j]);
      legAll->AddEntry(hnorm[j], Form("%.2g-%.2g", params::lumibins[j], params::lumibins[j + 1]), "p");
    }
    c->cd();
    hnorm.front()->Draw("axis");
    for (auto* h : hnorm) h->Draw("pe1 same");
    legAll->Draw();
    c->SaveAs((outdir + "/h_HF" + key + "_alllumibins.pdf").c_str());
  }

  // plot the fraction-above-threshold vs lumi, Plus and Minus overlaid
  c->SetLogy(0);
  auto* mg = new TMultiGraph();
  auto* leg = new TLegend(0.68, 0.74, 0.88, 0.88);
  xjjroot::setleg(leg, 0.04);
  for (int i = 0; i < (int)sides.size(); i++) {
    auto* g = g_frac_above.at(sides[i]);
    xjjroot::setthgrstyle(g, side_colors[i], side_markers[i], 1.3, side_colors[i]);
    mg->Add(g, "pe");
    leg->AddEntry(g, sides[i].c_str(), "pe");
  }
  c->cd();
  mg->Draw("ape");
  mg->GetXaxis()->SetTitle("instLumi");
  mg->GetYaxis()->SetTitle(Form("Fraction HFEMax > %g (0n side)", HF_0N_FRAC_THRESHOLD));
  leg->Draw();
  c->SaveAs((outdir + Form("/frac_above%g_vs_lumi.pdf", HF_0N_FRAC_THRESHOLD)).c_str());

  __XJJLOG << ">> plots saved under " << outdir << "/" << std::endl;

  outf->Close();
  inf->Close();

  return 0;
}

int main(int argc, char* argv[]) {
  if (argc == 4) {
    return macro(argv[1], argv[2], std::atoi(argv[3]));
  }
  if (argc == 3) {
    return macro(argv[1], argv[2]);
  }
  __XJJLOG << "!! usage: " << argv[0] << " <inputname.root (ntuple from ../analyze.cc)> <TAG (-> outputs/TAG[_leveled|_notleveled]/Main.root + *.pdf)> [leveled: -1 (default, no cut) | 0 (non-leveled only, tag gets _notleveled) | 1 (leveled only, tag gets _leveled)]" << std::endl;
  return 1;
}
