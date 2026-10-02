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
#include "zdcLumiCorrection.h" // lumi-dependent ZDC recalibration (parameters from calibZDC.exe)

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

  // CMS blue, used for the lumi-inclusive distribution (evaluated lazily: TColor needs ROOT initialized)
  Color_t cms_blue() { return TColor::GetColor("#5790fc"); }

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

int macro(const std::string& inputname, const std::string& tag, int leveled = -1, bool noZDCcut = false, int applyCorr = 0) {
  // leveled: -1 = no cut on lumiLeveled (default), 0 = only non-leveled events, 1 = only leveled events
  // noZDCcut: if true, no ZDC value cut is applied: every event fills both An0n and 0nAn histograms (tag gets _noZDCcut)
  // applyCorr: 0 = off (default); 1 = ZDCsumPlus/Minus replaced by the smooth (PCHIP) zdccorr::correct(ZDC, instLumi, side) before
  //            the filling (the 0n cut always uses the uncorrected sums) (tag gets _corr); 2 = same with the legacy piecewise-linear map (tag gets _corrlin)
  namespace fs = std::filesystem;

  // all output (histogram file + plots) goes under outputs/TAG[_leveled|_notleveled][_noZDCcut][_corr]/
  const std::string outdir = "outputs/" + tag + (leveled == 1 ? "_leveled" : leveled == 0 ? "_notleveled" : "") + (noZDCcut ? "_noZDCcut" : "") + (applyCorr == 1 ? "_corr" : applyCorr == 2 ? "_corrlin" : "");
  const bool smoothCorr = applyCorr != 2;
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

  // applied correction maps ZDC' = f_L(ZDC): one curve per lumi bin (evaluated at the bin centre) and side.
  // Drawn as corrected vs original ZDC sum on the full range and zoomed on the noise region, plus df/dx; the f(x) TGraphs are written to Main.root below
  std::vector<TGraph*> mapGraphs;
  if (applyCorr) {
    xjjroot::setgstyle(1);
    gStyle->SetPalette(kViridis);
    auto* cm = new TCanvas("c_map", "", 800, 600);
    const std::array<double, 3> xmax = { 8000., 300., 8000. };
    const std::array<std::pair<double, double>, 3> yr = { std::make_pair(0., 8000.), std::make_pair(0., 300.), std::make_pair(0., 2.0) };
        for (int s = 0; s < 2; s++) {
      const std::string sideName = s == 0 ? "Plus" : "Minus";
      for (int r = 0; r < 3; r++) {
        // legend strictly inside the frame: pad margins give the frame edges. r = 0,1 (corrected vs original): top-left, two
        // columns (the curves lie below the diagonal); r = 2 (slope): full-width strip at the top
        const double fx1 = cm->GetLeftMargin() + 0.015, fx2 = 1. - cm->GetRightMargin() - 0.015, fy2 = 1. - cm->GetTopMargin() - 0.015;
        auto* leg = r == 2 ? new TLegend(fx1, fy2 - 0.26, fx2, fy2) : new TLegend(fx1, fy2 - 0.40, fx1 + 0.40, fy2);
        leg->SetNColumns(r == 2 ? 4 : 2);
        leg->SetHeader(Form("ZDC %s, lumi bins [10^{-6}]", sideName.c_str()));
        xjjroot::setleg(leg, r == 2 ? 0.03 : 0.028);
        auto* frame = new TH2D(Form("frame_%s_%d", sideName.c_str(), r), Form(";original ZDCsum%s;%s", sideName.c_str(), r == 2 ? "df/dx of the map" : Form("corrected ZDCsum%s", sideName.c_str())),
                               10, 0., xmax[r], 10, yr[r].first, yr[r].second);
        frame->SetStats(0);
        cm->cd();
        frame->Draw("axis");
        for (int i = 0; i < nlumibins; i++) {
          const double lumiC = 0.5 * (params::lumibins[i] + params::lumibins[i + 1]);
          auto* g = new TGraph();
          g->SetName(Form("g_map_%s_lumibin%d", sideName.c_str(), i));
          g->SetTitle(Form(";ZDCsum%s;corrected ZDCsum%s", sideName.c_str(), sideName.c_str()));
          auto* gd = new TGraph(); // shift, for drawing only
          const double step = ZDC_MAX / 2000.;
          for (double x = 0.; x <= ZDC_MAX; x += step) {
            double slope = 1.;
            const double fx = zdccorr::correct(x, lumiC, s, smoothCorr, &slope);
            g->SetPoint(g->GetN(), x, fx);
            gd->SetPoint(gd->GetN(), x, r == 2 ? slope : fx);
          }
          const int col = TColor::GetColorPalette(i * (TColor::GetNumberOfColors() - 1) / std::max(1, nlumibins - 1));
          gd->SetLineColor(col); gd->SetLineWidth(2);
          gd->Draw("l same");
          leg->AddEntry(gd, Form("%.2g-%.2g", params::lumibins[i] / 1.e-6, params::lumibins[i + 1] / 1.e-6), "l");
          if (r == 0) mapGraphs.push_back(g);
        }
        if (r < 2) { // identity reference
          auto* id = new TLine(0., 0., xmax[r], xmax[r]);
          id->SetLineStyle(2); id->SetLineColor(kGray + 1);
          id->Draw();
        }
        leg->Draw();
        cm->SaveAs((outdir + (r == 2 ? "/slope_" : "/map_") + sideName + (r == 1 ? "_zoom" : "") + ".pdf").c_str());
      }
    }
    delete cm;
  }

  // make 4 sets of 10 histograms for ZDC values, 2 per ZDC side, 2 per 0nAn or An0n cut, 10 per lumi bins defined in params.h
  const std::array<std::string, 2> sides = { "Plus", "Minus" };
  const std::array<std::string, 2> cuts  = { "An0n", "0nAn" };
  // An0n: uncorrected ZDCsumPlus < 1500
  // 0nAn: uncorrected ZDCsumMinus < 1500

  std::map<std::string, std::vector<TH1D*>> hists; // key = "<side>_<cut>"
  std::map<std::string, TH1D*> hincl;              // same key, no lumi cut (lumi inclusive)
  for (const auto& side : sides) {
    for (const auto& cut : cuts) {
      auto& hs = hists[side + "_" + cut];
      hs.reserve(nlumibins);
      auto*& hi = hincl[side + "_" + cut];
      hi = new TH1D(Form("h_ZDC%s_%s_lumiincl", side.c_str(), cut.c_str()),
                    Form(";ZDCsum%s;Events", side.c_str()),
                    NBINS_ZDC, ZDC_MIN, ZDC_MAX);
      hi->Sumw2();
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

    const auto lumibin = find_lumibin(instLumi); // -1: outside lumibins, still counted in the lumi-inclusive histograms

    // the 0n cut is always applied on the uncorrected ZDC sum: it is motivated by the sample missing events with both sides above 1500
    const Float_t rawPlus = ZDCsumPlus, rawMinus = ZDCsumMinus;
    if (applyCorr) {
      ZDCsumPlus = zdccorr::correct(ZDCsumPlus, instLumi, zdccorr::kPlus, smoothCorr);
      ZDCsumMinus = zdccorr::correct(ZDCsumMinus, instLumi, zdccorr::kMinus, smoothCorr);
    }

    if (noZDCcut || rawPlus < ZDC_0N_CUT) { // An0n
      hincl.at("Plus_An0n")->Fill(ZDCsumPlus);
      hincl.at("Minus_An0n")->Fill(ZDCsumMinus);
      if (lumibin >= 0) {
        hists.at("Plus_An0n")[lumibin]->Fill(ZDCsumPlus);
        hists.at("Minus_An0n")[lumibin]->Fill(ZDCsumMinus);
        nevents["An0n"][lumibin]++;
      }
    }
    if (noZDCcut || rawMinus < ZDC_0N_CUT) { // 0nAn
      hincl.at("Plus_0nAn")->Fill(ZDCsumPlus);
      hincl.at("Minus_0nAn")->Fill(ZDCsumMinus);
      if (lumibin >= 0) {
        hists.at("Plus_0nAn")[lumibin]->Fill(ZDCsumPlus);
        hists.at("Minus_0nAn")[lumibin]->Fill(ZDCsumMinus);
        nevents["0nAn"][lumibin]++;
      }
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
  for (auto& [key, h] : hincl) h->Write();
  for (auto* g : fitGraphs) g->Write();
  for (auto* g : mapGraphs) g->Write();
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

    // lumi-inclusive distribution (no lumi cut), normalized the same way, in CMS blue to stand out
    auto* hinclNorm = (TH1D*)hincl.at(key)->Clone(Form("%s_norm", hincl.at(key)->GetName()));
    {
      const double integral = hinclNorm->Integral();
      if (integral > 0) hinclNorm->Scale(1. / integral);
      hinclNorm->GetYaxis()->SetTitle("Fraction of events / bin");
    }

    auto hall = hnorm; // lumi bins + inclusive, for the common headroom
    hall.push_back(hinclNorm);
    xjjana::sethsmax(hall, 1.5); // common headroom across the series, log-y min stays automatic
    auto cc = xjjroot::grayscales_alpha((int)hnorm.size());
    auto* legAll = new TLegend(0.70, 0.86 - 0.025 * 1.25 * hall.size(), 0.88, 0.86);
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
    c->cd();
    hnorm.front()->Draw("axis");
    for (auto* h : hnorm) h->Draw("pe1 same");
    hinclNorm->Draw("pe1 same");
    legAll->Draw();
    c->SaveAs((outdir + "/h_ZDC" + key + "_alllumibins.pdf").c_str());

    // ratio of each normalized lumi-bin spectrum to the lowest-lumi one: 1 everywhere = perfect lumi independence;
    // kinks/bumps at the map knots (noise edge, 1n -/+ sigma, 2n) show up here
    {
      gStyle->SetPalette(kViridis);
      std::vector<TH1D*> hratio;
      auto* legR = new TLegend(c->GetLeftMargin() + 0.015, 1. - c->GetTopMargin() - 0.015 - 0.20, 1. - c->GetRightMargin() - 0.015, 1. - c->GetTopMargin() - 0.015);
      legR->SetNColumns(4);
      xjjroot::setleg(legR, 0.025);
      for (int j = 1; j < (int)hnorm.size(); j++) {
        auto* hr = (TH1D*)hs[j]->Clone(Form("%s_ratio", hs[j]->GetName()));
        hr->Scale(1. / hs[j]->Integral());
        auto* h0 = (TH1D*)hs[0]->Clone(Form("%s_ref", hs[0]->GetName()));
        h0->Scale(1. / hs[0]->Integral());
        hr->Divide(h0);
        const int col = TColor::GetColorPalette((j - 1) * (TColor::GetNumberOfColors() - 1) / std::max(1, (int)hnorm.size() - 2));
        hr->SetLineColor(col); hr->SetMarkerColor(col); hr->SetMarkerStyle(1); hr->SetLineWidth(2);
        hr->GetYaxis()->SetTitle("spectrum / spectrum (lowest lumi bin)");
        hratio.push_back(hr);
        legR->AddEntry(hr, Form("%.2g-%.2g", params::lumibins[j] / 1.e-6, params::lumibins[j + 1] / 1.e-6), "l");
      }
      auto* frameR = new TH2D(Form("frameR_%s", key.c_str()), Form(";ZDCsum%s;spectrum / spectrum (lowest lumi bin)", key.substr(0, key.find('_')).c_str()), 10, 0., ZDC_MAX, 10, 0., 3.6);
      frameR->SetStats(0);
      c->cd();
      c->SetLogy(0);
      frameR->Draw("axis");
      for (auto* hr : hratio) hr->Draw("hist same");
      legR->Draw();
      c->SaveAs((outdir + "/h_ZDC" + key + "_ratio_to_bin0.pdf").c_str());
      c->SetLogy(1);
    }
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
  if (argc == 6) {
    return macro(argv[1], argv[2], std::atoi(argv[3]), std::atoi(argv[4]) != 0, std::atoi(argv[5]));
  }
  if (argc == 5) {
    return macro(argv[1], argv[2], std::atoi(argv[3]), std::atoi(argv[4]) != 0);
  }
  if (argc == 4) {
    return macro(argv[1], argv[2], std::atoi(argv[3]));
  }
  if (argc == 3) {
    return macro(argv[1], argv[2]);
  }
  __XJJLOG << "!! usage: " << argv[0] << " <inputname.root (ntuple from ../analyze.cc)> <TAG (-> outputs/TAG[_leveled|_notleveled][_noZDCcut]/Main.root + *.pdf)> [leveled: -1 (default, no cut) | 0 (non-leveled only, tag gets _notleveled) | 1 (leveled only, tag gets _leveled)] [noZDCcut: 0 (default, apply ZDC 0n cuts) | 1 (no ZDC value cut, tag gets _noZDCcut)] [applyCorr: 0 (default) | 1 (smooth zdcLumiCorrection.h map, tag gets _corr) | 2 (legacy piecewise-linear map, tag gets _corrlin)]" << std::endl;
  return 1;
}
