#include <filesystem>
#include <algorithm>
#include <array>
#include <cmath>

#include <TGraph.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TLine.h>

#include "xjjanauti.h"
#include "../params.h"
#include "hfLumiCorrection.h"

// Plots the lumi-dependent HFEMax correction map f_L(x) (hfLumiCorrection.h, parameters from calibHF.cc) against the original HFEMax.
// No data input: the curves are evaluated directly at the middle of every lumi bin of params.h (lumi below the first fitted bin is
// clamped by the correction, i.e. those curves are the identity). Writes under outputs/TAG/
//   hf_map_<side>.pdf : left f_L(x) vs x, right f_L(x) - x, one curve per lumi bin (dashed = identity)
//   Main.root         : the graphs

namespace {
  const double HF_MIN = 0., HF_MAX = 60.;
  const int NX = 600;
  const std::array<std::string, 2> sides = { "Plus", "Minus" };
}

int macro(const std::string& tag) {
  namespace fs = std::filesystem;
  const std::string outdir = "outputs/" + tag;
  fs::create_directories(outdir);
  auto* outf = TFile::Open((outdir + "/Main.root").c_str(), "RECREATE");
  xjjroot::setgstyle(1);
  gStyle->SetPalette(kBird);

  const auto& edges = params::lumibins;
  const int nb = (int)edges.size() - 1;
  for (int s = 0; s < 2; s++) {
    auto* cm = new TCanvas(Form("cm_%d", s), "", 1600, 700);
    cm->Divide(2, 1);
    auto* leg = new TLegend(0.14, 0.45, 0.5, 0.88, "<L> (10^{33} cm^{-2}s^{-1}) of bin");
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.03); leg->SetNColumns(2);
    auto* mg = new TMultiGraph(Form("mg_map_%s", sides[s].c_str()), Form("%s;HFEMax original (GeV);HFEMax corrected (GeV)", sides[s].c_str()));
    auto* md = new TMultiGraph(Form("mg_diff_%s", sides[s].c_str()), Form("%s;HFEMax original (GeV);corrected - original (GeV)", sides[s].c_str()));
    for (int i = 0; i < nb; i++) {
      const double lumi = 0.5 * (edges[i] + edges[i + 1]);
      auto* g = new TGraph(); auto* gd = new TGraph();
      for (int j = 0; j <= NX; j++) {
        const double x = HF_MIN + (HF_MAX - HF_MIN) * j / NX;
        const double y = hfcorr::correct(x, lumi, s == 0 ? hfcorr::kPlus : hfcorr::kMinus);
        g->SetPoint(j, x, y); gd->SetPoint(j, x, y - x);
      }
      const Color_t col = TColor::GetColorPalette(int(254. * i / std::max(nb - 1, 1)));
      for (auto* gg : { g, gd }) { gg->SetLineColor(col); gg->SetLineWidth(2); }
      g->SetName(Form("g_map_%s_bin%d", sides[s].c_str(), i));
      gd->SetName(Form("g_diff_%s_bin%d", sides[s].c_str(), i));
      mg->Add(g, "l"); md->Add(gd, "l");
      g->Write(); gd->Write();
      leg->AddEntry(g, Form("%d: %.2f", i, lumi / 1.e-6), "l");
    }
    cm->cd(1); mg->Draw("a");
    mg->GetXaxis()->SetLimits(HF_MIN, HF_MAX); mg->GetYaxis()->SetRangeUser(HF_MIN, HF_MAX);
    TLine d1(HF_MIN, HF_MIN, HF_MAX, HF_MAX); d1.SetLineStyle(2); d1.SetLineColor(kGray + 2); d1.Draw();
    leg->Draw();
    cm->cd(2); md->Draw("a");
    md->GetXaxis()->SetLimits(HF_MIN, HF_MAX);
    TLine d0(HF_MIN, 0., HF_MAX, 0.); d0.SetLineStyle(2); d0.SetLineColor(kGray + 2); d0.Draw();
    cm->SaveAs((outdir + "/hf_map_" + sides[s] + ".pdf").c_str());
  }
  outf->Close();
  __XJJLOG << ">> wrote " << outdir << "/" << std::endl;
  return 0;
}

int main(int argc, char* argv[]) {
  if (argc == 2) return macro(argv[1]);
  __XJJLOG << "!! usage: " << argv[0] << " <TAG (-> outputs/TAG/hf_map_<side>.pdf)>" << std::endl;
  return 1;
}
