#include <filesystem>
#include <algorithm>
#include <array>
#include <map>
#include <cmath>

#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooFitResult.h"
#include "RooPlot.h"

#include <TMultiGraph.h>

#include "xjjanauti.h"
#include "../params.h"

namespace {
  const int NBINS_ZDC = 200;
  const double ZDC_MIN = 0., ZDC_MAX = 10000.;
  const double ZDC_0N_CUT = 1500.; // ZDCsum below this = "0n" (no neutron) on that side
  const double PEAK_FIT_MIN = 1500., PEAK_FIT_MAX = 3500.; // fit range for the first (1n) peak
  const double PEAK_FIT_SIGMA_INIT = 250.;
  // if true, after the first fit (range PEAK_FIT_MIN-PEAK_FIT_MAX), refit with range = fitted mean +- REFIT_WINDOW_HALFWIDTH
  const bool REFIT_AROUND_MEAN = true;
  const double REFIT_WINDOW_HALFWIDTH = 1000.;
  // the two histogram keys whose side is unconstrained ("An"): the other side carries the 0n cut
  const std::array<std::string, 2> AN_SIDE_KEYS = { "Minus_An0n", "Plus_0nAn" };

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
  const std::string levelLabel = leveled == 1 ? "leveled" : leveled == 0 ? "not leveled" : "all";
  fs::create_directories(outdir);
  const std::string outputname = outdir + "/Main.root";

  // get input file (ntuple produced by ../analyze.cc)
  auto* inf = TFile::Open(inputname.c_str());
  if (!inf || inf->IsZombie()) { __XJJLOG << "!! failed to open input file, abort." << std::endl; return 1; }
  auto* nt = (TTree*)inf->Get("ntuple");
  if (!nt) { __XJJLOG << "!! no ntuple found in input file, abort." << std::endl; return 2; }

  Float_t instLumi = 0., ZDCsumPlus = 0., ZDCsumMinus = 0., lumiLeveled = 0.;
  nt->SetBranchStatus("*", 0);
  nt->SetBranchStatus("instLumi", 1);
  nt->SetBranchStatus("ZDCsumPlus", 1);
  nt->SetBranchStatus("ZDCsumMinus", 1);
  nt->SetBranchAddress("instLumi", &instLumi);
  nt->SetBranchAddress("ZDCsumPlus", &ZDCsumPlus);
  nt->SetBranchAddress("ZDCsumMinus", &ZDCsumMinus);
  if (leveled >= 0) {
    nt->SetBranchStatus("lumiLeveled", 1);
    nt->SetBranchAddress("lumiLeveled", &lumiLeveled);
  }

  const auto nentries = nt->GetEntries();
  __XJJLOG << ">> " << nentries << " entries loaded from " << inputname << std::endl;

  const int nlumibins = (int)params::lumibins.size() - 1;

  // make 4 sets of 10 histograms for ZDC values, 2 per ZDC side, 2 per 0nAn or An0n cut, 10 per lumi bins defined in params.h
  const std::array<std::string, 2> sides = { "Plus", "Minus" };
  const std::array<std::string, 2> cuts  = { "An0n", "0nAn" };
  // An0n: ZDCsumPlus < 1500
  // 0nAn: ZDCsumMinus < 1500

  std::map<std::string, std::vector<TH1D*>> hists; // key = "<side>_<cut>"
  for (const auto& side : sides) {
    for (const auto& cut : cuts) {
      auto& hs = hists[side + "_" + cut];
      hs.reserve(nlumibins);
      for (int i = 0; i < nlumibins; i++) {
        auto* h = new TH1D(Form("h_ZDC%s_%s_lumibin%d", side.c_str(), cut.c_str(), i),
                            Form(";ZDCsum%s;Events", side.c_str()),
                            NBINS_ZDC, ZDC_MIN, ZDC_MAX);
        h->Sumw2();
        hs.push_back(h);
      }
    }
  }

  // number of events passing each cut topology (An0n/0nAn), per lumi bin (independent of side,
  // since both sides of a given cut are filled together below from the same selected events)
  std::map<std::string, std::vector<Long64_t>> nevents; // key = cut ("An0n"/"0nAn")
  for (const auto& cut : cuts) nevents[cut].assign(nlumibins, 0);

  // fill histograms with ZDC values for each lumi bin, for each ZDC side, for each 0nAn or An0n cut
  for (Long64_t i = 0; i < nentries; i++) {
    nt->GetEntry(i);
    if (i % 1000000 == 0) __XJJLOG << ">> processing entry " << i << " / " << nentries << std::endl;

    if (leveled >= 0 && (int)lumiLeveled != leveled) continue;

    const auto lumibin = find_lumibin(instLumi);
    if (lumibin < 0) continue;

    if (ZDCsumPlus < ZDC_0N_CUT) { // An0n
      hists.at("Plus_An0n")[lumibin]->Fill(ZDCsumPlus);
      hists.at("Minus_An0n")[lumibin]->Fill(ZDCsumMinus);
      nevents["An0n"][lumibin]++;
    }
    if (ZDCsumMinus < ZDC_0N_CUT) { // 0nAn
      hists.at("Plus_0nAn")[lumibin]->Fill(ZDCsumPlus);
      hists.at("Minus_0nAn")[lumibin]->Fill(ZDCsumMinus);
      nevents["0nAn"][lumibin]++;
    }
  }


  // fit the first (1n) peak of each "An"-side histogram with a RooFit Gaussian, in PEAK_FIT_MIN-PEAK_FIT_MAX
  // fit_mean/fit_sigma (+ errors) are keyed by side ("Plus"/"Minus"), one entry per lumi bin
  std::map<std::string, std::vector<double>> fit_mean, fit_meanErr, fit_sigma, fit_sigmaErr;
  for (const auto& side : sides) {
    fit_mean[side].assign(nlumibins, 0.);
    fit_meanErr[side].assign(nlumibins, 0.);
    fit_sigma[side].assign(nlumibins, 0.);
    fit_sigmaErr[side].assign(nlumibins, 0.);
  }

  xjjroot::setgstyle(1);
  auto* c_fit = new TCanvas("c_fit", "", 800, 600);
  for (const auto& key : AN_SIDE_KEYS) {
    const std::string side = key.substr(0, key.find('_')); // "Plus" or "Minus"
    auto& hs = hists.at(key);
    for (int i = 0; i < nlumibins; i++) {
      auto* h = hs[i];

      // initial mean = position of the histogram max within the fit range
      double meanInit = h->GetXaxis()->GetBinCenter(h->GetMaximumBin());
      double maxContent = -1.;
      for (int b = h->GetXaxis()->FindBin(PEAK_FIT_MIN); h->GetXaxis()->GetBinCenter(b) < PEAK_FIT_MAX; b++) {
        if (h->GetBinContent(b) > maxContent) { maxContent = h->GetBinContent(b); meanInit = h->GetXaxis()->GetBinCenter(b); }
      }

      RooRealVar x(Form("x_%s", h->GetName()), h->GetXaxis()->GetTitle(), ZDC_MIN, ZDC_MAX);
      x.setBins(NBINS_ZDC);
      x.setRange("fitRange", PEAK_FIT_MIN, PEAK_FIT_MAX);
      RooDataHist dh(Form("dh_%s", h->GetName()), "", RooArgList(x), RooFit::Import(*h));

      RooRealVar mean(Form("mean_%s", h->GetName()), "mean", meanInit, PEAK_FIT_MIN, PEAK_FIT_MAX);
      RooRealVar sigma(Form("sigma_%s", h->GetName()), "sigma", PEAK_FIT_SIGMA_INIT, 50., 1500.);
      RooGaussian gauss(Form("gauss_%s", h->GetName()), "gauss", x, mean, sigma);

      auto* fitr = gauss.fitTo(dh, RooFit::Range("fitRange"), RooFit::Save(true), RooFit::PrintLevel(-1));

      std::string fitRangeName = "fitRange";
      if (REFIT_AROUND_MEAN) {
        const double refitMin = std::max(ZDC_MIN, mean.getVal() - REFIT_WINDOW_HALFWIDTH);
        const double refitMax = std::min(ZDC_MAX, mean.getVal() + REFIT_WINDOW_HALFWIDTH);
        x.setRange("fitRangeRefit", refitMin, refitMax);
        delete fitr;
        fitr = gauss.fitTo(dh, RooFit::Range("fitRangeRefit"), RooFit::Save(true), RooFit::PrintLevel(-1));
        fitRangeName = "fitRangeRefit";
      }

      __XJJLOG << ">> fit " << h->GetName() << ": status=" << fitr->status()
               << " mean=" << mean.getVal() << "+-" << mean.getError()
               << " sigma=" << sigma.getVal() << "+-" << sigma.getError() << std::endl;

      fit_mean[side][i] = mean.getVal();
      fit_meanErr[side][i] = mean.getError();
      fit_sigma[side][i] = sigma.getVal();
      fit_sigmaErr[side][i] = sigma.getError();

      c_fit->cd();
      c_fit->SetLogy(1);
      auto* frame = x.frame(RooFit::Title(h->GetName()));
      dh.plotOn(frame);
      gauss.plotOn(frame, RooFit::Range(fitRangeName.c_str()), RooFit::NormRange(fitRangeName.c_str()), RooFit::LineColor(kRed+1));
      frame->SetMinimum(0.5);
      frame->Draw();
      c_fit->SaveAs((outdir + "/" + h->GetName() + "_fit.pdf").c_str());
      delete frame;
      delete fitr;
    }
  }

  // build mean/sigma vs lumi graphs (one point per lumi bin, x error = half bin width)
  std::vector<double> lumiCenter(nlumibins), lumiHalfWidth(nlumibins);
  for (int i = 0; i < nlumibins; i++) {
    lumiCenter[i] = 0.5 * (params::lumibins[i] + params::lumibins[i + 1]);
    lumiHalfWidth[i] = 0.5 * (params::lumibins[i + 1] - params::lumibins[i]);
  }
  std::vector<TGraphErrors*> fitGraphs;
  for (const auto& side : sides) {
    auto* g_mean = new TGraphErrors(nlumibins, lumiCenter.data(), fit_mean.at(side).data(),
                                    lumiHalfWidth.data(), fit_meanErr.at(side).data());
    g_mean->SetName(Form("g_mean_vs_lumi_%s", side.c_str()));
    g_mean->SetTitle(Form(";instLumi;Gaussian mean (ZDCsum%s)", side.c_str()));
    fitGraphs.push_back(g_mean);

    auto* g_sigma = new TGraphErrors(nlumibins, lumiCenter.data(), fit_sigma.at(side).data(),
                                     lumiHalfWidth.data(), fit_sigmaErr.at(side).data());
    g_sigma->SetName(Form("g_sigma_vs_lumi_%s", side.c_str()));
    g_sigma->SetTitle(Form(";instLumi;Gaussian sigma (ZDCsum%s)", side.c_str()));
    fitGraphs.push_back(g_sigma);
  }

  // number-of-events-vs-lumi graphs, one per cut topology (An0n/0nAn); y-error = Poisson sqrt(n)
  std::map<std::string, TGraphErrors*> g_nevents;
  for (const auto& cut : cuts) {
    std::vector<double> n(nlumibins), nErr(nlumibins);
    for (int i = 0; i < nlumibins; i++) {
      n[i] = (double)nevents.at(cut)[i];
      nErr[i] = std::sqrt(n[i]);
    }
    auto* g = new TGraphErrors(nlumibins, lumiCenter.data(), n.data(), lumiHalfWidth.data(), nErr.data());
    g->SetName(Form("g_nevents_vs_lumi_%s", cut.c_str()));
    g->SetTitle(Form(";instLumi;N events (%s cut, %s)", cut.c_str(), levelLabel.c_str()));
    g_nevents[cut] = g;
    fitGraphs.push_back(g);
  }

  // save histograms and fit-result graphs to output file
  auto* outf = TFile::Open(outputname.c_str(), "RECREATE");
  if (!outf || outf->IsZombie()) { __XJJLOG << "!! failed to open output file, abort." << std::endl; return 3; }
  outf->cd();
  for (auto& [key, hs] : hists)
    for (auto* h : hs) h->Write();
  for (auto* g : fitGraphs) g->Write();
  __XJJLOG << ">> histograms and fit-result graphs written to " << outputname << std::endl;

  // plot each histogram under outputs/TAG/.....pdf
  auto* c = new TCanvas("c", "", 800, 600);
  c->SetLogy(1);
  for (auto& [key, hs] : hists) {
    for (auto* h : hs) {
      c->cd();
      h->Draw("hist e");
      c->SaveAs((outdir + "/" + h->GetName() + ".pdf").c_str());
    }
  }

  // for each ZDC "An"-side histogram, overlay all lumi bins on one canvas: black (lowest lumi) fading
  // to gray (highest lumi), distinct marker shape per lumi bin, same style as main/hist_fit.cc.
  // each histogram is self-normalized (unit integral) first, so the shapes are comparable across
  // lumi bins regardless of differing statistics; originals in `hists` are left untouched.
  for (const auto& key : AN_SIDE_KEYS) {
    auto& hs = hists.at(key);
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
    c->SaveAs((outdir + "/h_ZDC" + key + "_alllumibins.pdf").c_str());
  }

  // number of events passing each cut topology (An0n/0nAn) vs lumi, overlaid; reflects the
  // currently selected lumiLeveled choice (levelLabel), since nevents was counted after that cut
  {
    const std::array<Color_t, 2> cutColors = { kBlue + 1, kRed + 1 };
    const std::array<Style_t, 2> cutMarkers = { 20, 21 };
    c->SetLogy(1);
    auto* mg = new TMultiGraph();
    auto* leg = new TLegend(0.68, 0.76, 0.88, 0.88);
    xjjroot::setleg(leg, 0.04);
    double ymin = 1.e+300, ymax = -1.e+300;
    for (int i = 0; i < (int)cuts.size(); i++) {
      auto* g = g_nevents.at(cuts[i]);
      xjjroot::setthgrstyle(g, cutColors[i], cutMarkers[i], 1.3, cutColors[i]);
      mg->Add(g, "pe");
      leg->AddEntry(g, cuts[i].c_str(), "pe");
      for (int k = 0; k < g->GetN(); k++) {
        ymin = std::min(ymin, g->GetY()[k] - g->GetEY()[k]);
        ymax = std::max(ymax, g->GetY()[k] + g->GetEY()[k]);
      }
    }
    mg->SetMinimum((ymin > 0 ? ymin : ymax * 1.e-3) / 2.); // headroom so points/error bars aren't clipped at the frame edge
    mg->SetMaximum(ymax * 2.);
    c->cd();
    mg->Draw("ape");
    mg->GetXaxis()->SetTitle("instLumi");
    mg->GetYaxis()->SetTitle(Form("N events (%s)", levelLabel.c_str()));
    leg->Draw();
    c->SaveAs((outdir + "/nevents_vs_lumi.pdf").c_str());
  }

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
