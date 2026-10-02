#include <filesystem>
#include <algorithm>
#include <array>
#include <map>
#include <cmath>

#include <TMultiGraph.h>

#include "xjjanauti.h"
#include "../params.h"
#include "hfLumiCorrection.h"

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

  // CMS blue, used for the lumi-inclusive distribution (evaluated lazily: TColor needs ROOT initialized)
  Color_t cms_blue() { return TColor::GetColor("#5790fc"); }

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

int macro(const std::string& inputname, const std::string& tag, int leveled = -1, bool noZDCcut = false, int applyCorr = 0, int parity = 0) {
  // leveled: -1 = no cut on lumiLeveled (default), 0 = only non-leveled events, 1 = only leveled events
  // noZDCcut: if true, no ZDC value cut is applied: every event fills both An0n and 0nAn histograms (tag gets _noZDCcut)
  // applyCorr: 1 = apply the lumi-dependent HFEMax correction (hfLumiCorrection.h, calibHF.cc) to HFEMax on both the 0n and An side (tag gets _corr)
  // parity: 0 = all entries, 1 = even entries only, 2 = odd entries only (split-sample closure; tag gets _even / _odd)
  namespace fs = std::filesystem;

  // all output (histogram file + plots) goes under outputs/TAG[_leveled|_notleveled][_noZDCcut]/
  const std::string outdir = "outputs/" + tag + (leveled == 1 ? "_leveled" : leveled == 0 ? "_notleveled" : "") + (noZDCcut ? "_noZDCcut" : "") + (applyCorr ? "_corr" : "") + (parity == 1 ? "_even" : parity == 2 ? "_odd" : "");
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
  std::map<std::string, TH1D*> hincl;              // same key, no lumi cut (lumi inclusive)
  for (const auto& side : sides) {
    for (const auto& cut : cuts) {
      const std::string key = side + "_" + cut;
      const bool anSide = noZDCcut || is_an_side(key); // no ZDC cut: no side is 0n-constrained, use the wide range
      const int nbins = anSide ? NBINS_HF_AN : NBINS_HF_0N;
      const double hmin = anSide ? HF_AN_MIN : HF_0N_MIN;
      const double hmax = anSide ? HF_AN_MAX : HF_0N_MAX;
      auto& hs = hists[key];
      hs.reserve(nlumibins);
      hincl[key] = new TH1D(Form("h_HF%s_%s_lumiincl", side.c_str(), cut.c_str()),
                            Form(";HFEMax%s_eta5;Events", side.c_str()),
                            nbins, hmin, hmax);
      hincl[key]->Sumw2();
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

    if (parity == 1 && (i % 2) != 0) continue;
    if (parity == 2 && (i % 2) != 1) continue;
    if (leveled >= 0 && (int)lumiLeveled != leveled) continue;

    const auto lumibin = find_lumibin(instLumi); // -1: outside lumibins, still counted in the lumi-inclusive histograms

    // lumi-dependent correction of HFEMax on both the 0n and the An side (the map was calibrated on the 0n side only,
    // per HF side; the ZDC 0n tag always uses the raw ZDC sum)
    const double hfPlus  = applyCorr ? hfcorr::correct(HFEMaxPlus,  instLumi, hfcorr::kPlus)  : HFEMaxPlus;
    const double hfMinus = applyCorr ? hfcorr::correct(HFEMaxMinus, instLumi, hfcorr::kMinus) : HFEMaxMinus;

    if (noZDCcut || ZDCsumPlus < ZDC_0N_CUT) { // An0n: Plus is the 0n side
      hincl.at("Plus_An0n")->Fill(hfPlus);
      hincl.at("Minus_An0n")->Fill(hfMinus);
      if (lumibin >= 0) {
        hists.at("Plus_An0n")[lumibin]->Fill(hfPlus);
        hists.at("Minus_An0n")[lumibin]->Fill(hfMinus);
      }
    }
    if (noZDCcut || ZDCsumMinus < ZDC_0N_CUT) { // 0nAn: Minus is the 0n side
      hincl.at("Plus_0nAn")->Fill(hfPlus);
      hincl.at("Minus_0nAn")->Fill(hfMinus);
      if (lumibin >= 0) {
        hists.at("Plus_0nAn")[lumibin]->Fill(hfPlus);
        hists.at("Minus_0nAn")[lumibin]->Fill(hfMinus);
      }
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

  // fraction above several thresholds (0n side), per lumi bin, as a ratio to the lumi-inclusive fraction: flat = no lumi dependence
  const std::array<double, 6> thresholds = { 8., 10., 12., 16., 20., 30. };
  std::map<std::string, std::vector<TGraphErrors*>> g_thr; // keyed by side
  for (const auto& key : ZN_SIDE_KEYS) {
    const std::string side = key.substr(0, key.find('_'));
    auto* hi = hincl.at(key);
    for (double thr : thresholds) {
      const double ni = hi->Integral(0, hi->GetNbinsX() + 1);
      const double pi = hi->Integral(hi->GetXaxis()->FindBin(thr), hi->GetNbinsX() + 1) / ni;
      std::vector<double> y(nlumibins), ye(nlumibins);
      for (int i = 0; i < nlumibins; i++) {
        auto* h = hists.at(key)[i];
        const double n = h->Integral(0, h->GetNbinsX() + 1);
        const double p = h->Integral(h->GetXaxis()->FindBin(thr), h->GetNbinsX() + 1) / n;
        y[i] = p / pi;
        ye[i] = std::sqrt(p * (1. - p) / n) / pi;
      }
      auto* g = new TGraphErrors(nlumibins, lumiCenter.data(), y.data(), lumiHalfWidth.data(), ye.data());
      g->SetName(Form("g_fracratio_above%g_vs_lumi_%s", thr, side.c_str()));
      g_thr[side].push_back(g);
    }
  }

  // save histograms and fraction-vs-lumi graphs to output file
  auto* outf = TFile::Open(outputname.c_str(), "RECREATE");
  if (!outf || outf->IsZombie()) { __XJJLOG << "!! failed to open output file, abort." << std::endl; return 3; }
  outf->cd();
  for (auto& [key, hs] : hists)
    for (auto* h : hs) h->Write();
  for (auto& [key, h] : hincl) h->Write();
  for (auto& [side, g] : g_frac_above) g->Write();
  for (auto& [side, gs] : g_thr) for (auto* g : gs) g->Write();
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

    // lumi-inclusive distribution (no lumi cut), normalized the same way
    auto* hinclNorm = (TH1D*)hincl.at(key)->Clone(Form("%s_norm", hincl.at(key)->GetName()));
    {
      const double integral = hinclNorm->Integral();
      if (integral > 0) hinclNorm->Scale(1. / integral);
      hinclNorm->GetYaxis()->SetTitle("Fraction of events / bin");
    }

    auto cc = xjjroot::grayscales_alpha((int)hnorm.size());
    auto* legAll = new TLegend(0.70, 0.86 - 0.025 * 1.25 * (hnorm.size() + 1), 0.88, 0.86);
    xjjroot::setleg(legAll, 0.025);
    for (int j = 0; j < (int)hnorm.size(); j++) {
      xjjroot::sethempty(hnorm[j], 0, 0.5);
      const auto mstyle = xjjroot::markerlist_solid[j % xjjroot::markerlist_solid.size()];
      xjjroot::setthgrstyle(hnorm[j], kBlack, mstyle, 1.0, kBlack, 1, 1, -1, -1, -1, cc[j], cc[j]);
      legAll->AddEntry(hnorm[j], Form("%.2g-%.2g", params::lumibins[j], params::lumibins[j + 1]), "p");
    }
    xjjroot::sethempty(hinclNorm, 0, 0.5);
    xjjroot::setthgrstyle(hinclNorm, cms_blue(), 20, 1.2, cms_blue());
    legAll->AddEntry(hinclNorm, "Lumi inclusive", "p");

    // fixed y range for the overlay: 1 down to 1e-7
    for (auto* h : hnorm) { h->SetMinimum(1.e-7); h->SetMaximum(1.); }
    hinclNorm->SetMinimum(1.e-7);
    hinclNorm->SetMaximum(1.);

    // top pad: overlay, bottom pad: ratio to the lumi-inclusive distribution
    auto* cov = new TCanvas(Form("cov_%s", key.c_str()), "", 800, 800);
    auto* pad1 = new TPad(Form("pad1_%s", key.c_str()), "", 0., 0.3, 1., 1.);
    auto* pad2 = new TPad(Form("pad2_%s", key.c_str()), "", 0., 0., 1., 0.3);
    pad1->SetBottomMargin(0.03); pad1->SetLogy(1);
    pad2->SetTopMargin(0.03);    pad2->SetBottomMargin(0.35);
    cov->cd(); pad1->Draw(); pad2->Draw();

    pad1->cd();
    hnorm.front()->GetXaxis()->SetLabelSize(0);
    hnorm.front()->GetXaxis()->SetTitleSize(0);
    hnorm.front()->Draw("axis");
    for (auto* h : hnorm) h->Draw("pe1 same");
    hinclNorm->Draw("pe1 same");
    legAll->Draw();

    pad2->cd();
    std::vector<TH1D*> hratio;
    for (int j = 0; j < (int)hnorm.size(); j++) {
      auto* hr = (TH1D*)hnorm[j]->Clone(Form("%s_ratio", hnorm[j]->GetName()));
      hr->Divide(hinclNorm);
      for (int b = 1; b <= hr->GetNbinsX(); b++) {
        // empty bins: park the point far outside the pad range instead of drawing a ratio of 0
        if (hnorm[j]->GetBinContent(b) <= 0 || hinclNorm->GetBinContent(b) <= 0) {
          hr->SetBinContent(b, -999.);
          hr->SetBinError(b, 0.);
        }
      }
      hr->SetMinimum(0.);
      hr->SetMaximum(2.);
      hratio.push_back(hr);
    }
    auto* hr0 = hratio.front();
    hr0->GetYaxis()->SetTitle("Ratio to lumi incl.");
    hr0->GetYaxis()->SetNdivisions(505);
    hr0->GetYaxis()->SetTitleSize(0.11);
    hr0->GetYaxis()->SetLabelSize(0.10);
    hr0->GetYaxis()->SetTitleOffset(0.5);
    hr0->GetXaxis()->SetTitle(hnorm.front()->GetXaxis()->GetTitle());
    hr0->GetXaxis()->SetTitleSize(0.13);
    hr0->GetXaxis()->SetLabelSize(0.11);
    hr0->GetXaxis()->SetTitleOffset(1.1);
    hr0->Draw("axis");
    for (auto* hr : hratio) hr->Draw("pe1 same");
    TLine lone(hr0->GetXaxis()->GetXmin(), 1., hr0->GetXaxis()->GetXmax(), 1.);
    lone.SetLineColor(cms_blue());
    lone.SetLineStyle(2);
    lone.Draw();

    cov->SaveAs((outdir + "/h_HF" + key + "_alllumibins.pdf").c_str());
    delete cov;
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

  // fraction above each threshold / lumi-inclusive fraction vs lumi, one canvas per side
  for (const auto& side : sides) {
    auto* cth = new TCanvas(Form("cth_%s", side.c_str()), "", 800, 600);
    auto* mgt = new TMultiGraph();
    auto* legt = new TLegend(0.70, 0.62, 0.88, 0.88);
    xjjroot::setleg(legt, 0.035);
    for (int t = 0; t < (int)thresholds.size(); t++) {
      auto* g = g_thr.at(side)[t];
      const Color_t col = xjjroot::colorlist_middle[t % xjjroot::colorlist_middle.size()];
      xjjroot::setthgrstyle(g, col, xjjroot::markerlist_solid[t % xjjroot::markerlist_solid.size()], 1.0, col);
      mgt->Add(g, "pe");
      legt->AddEntry(g, Form("> %g GeV", thresholds[t]), "pe");
    }
    cth->cd();
    mgt->Draw("ape");
    mgt->GetXaxis()->SetTitle("instLumi");
    mgt->GetYaxis()->SetTitle(Form("Fraction HFEMax%s > thr (0n side) / lumi-inclusive", side.c_str()));
    mgt->SetMinimum(0.7); mgt->SetMaximum(1.5);
    legt->Draw();
    TLine l1(params::lumibins.front(), 1., params::lumibins.back(), 1.);
    l1.SetLineStyle(2); l1.Draw();
    cth->SaveAs((outdir + "/fracratio_thresholds_vs_lumi_" + side + ".pdf").c_str());
    delete cth;
  }

  __XJJLOG << ">> plots saved under " << outdir << "/" << std::endl;

  outf->Close();
  inf->Close();

  return 0;
}

int main(int argc, char* argv[]) {
  if (argc >= 3 && argc <= 7) {
    return macro(argv[1], argv[2], argc > 3 ? std::atoi(argv[3]) : -1, argc > 4 && std::atoi(argv[4]) != 0,
                 argc > 5 ? std::atoi(argv[5]) : 0, argc > 6 ? std::atoi(argv[6]) : 0);
  }
  __XJJLOG << "!! usage: " << argv[0] << " <inputname.root (ntuple from ../analyze.cc)> <TAG (-> outputs/TAG[_leveled|_notleveled][_noZDCcut][_corr][_even|_odd]/Main.root + *.pdf)> [leveled: -1 (default, no cut) | 0 (non-leveled only, tag gets _notleveled) | 1 (leveled only, tag gets _leveled)] [noZDCcut: 0 (default, apply ZDC 0n cuts) | 1 (no ZDC value cut, tag gets _noZDCcut)] [applyCorr: 0 (default) | 1 (lumi-dependent HFEMax correction on 0n and An side, tag gets _corr)] [parity: 0 (default, all entries) | 1 (even entries) | 2 (odd entries)]" << std::endl;
  return 1;
}
