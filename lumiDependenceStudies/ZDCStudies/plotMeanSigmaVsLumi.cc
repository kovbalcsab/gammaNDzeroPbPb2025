#include <TMultiGraph.h>

#include "xjjanauti.h"

namespace {
  const std::array<std::string, 2> sides = { "Plus", "Minus" };
  const std::array<Color_t, 2> side_colors = { kBlue + 1, kRed + 1 };
  const std::array<Style_t, 2> side_markers = { 20, 21 };
}

// reads outputs/TAG/Main.root (produced by plotAndFit.cc) and draws, for the two ZDC sides overlaid:
//  - outputs/TAG/mean_vs_lumi.pdf  : fitted Gaussian mean vs instLumi
//  - outputs/TAG/sigma_vs_lumi.pdf : fitted Gaussian sigma vs instLumi
int macro(const std::string& tag) {
  const std::string outdir = "outputs/" + tag;
  const std::string inputname = outdir + "/Main.root";

  auto* inf = TFile::Open(inputname.c_str());
  if (!inf || inf->IsZombie()) { __XJJLOG << "!! failed to open " << inputname << ", abort." << std::endl; return 1; }

  std::map<std::string, TGraphErrors*> g_mean, g_sigma;
  for (const auto& side : sides) {
    auto* gm = (TGraphErrors*)inf->Get(Form("g_mean_vs_lumi_%s", side.c_str()));
    auto* gs = (TGraphErrors*)inf->Get(Form("g_sigma_vs_lumi_%s", side.c_str()));
    if (!gm || !gs) { __XJJLOG << "!! missing fit-result graph for side " << side << ", abort." << std::endl; return 2; }
    g_mean[side] = gm;
    g_sigma[side] = gs;
  }

  xjjroot::setgstyle(1);
  auto* c = new TCanvas("c", "", 800, 600);

  auto draw_overlay = [&](std::map<std::string, TGraphErrors*>& graphs, const std::string& ytitle, const std::string& outname,
                          double ymin, double ymax) {
    c->cd();
    auto* mg = new TMultiGraph();
    auto* leg = new TLegend(0.68, 0.74, 0.88, 0.88);
    xjjroot::setleg(leg, 0.04);
    for (int i = 0; i < (int)sides.size(); i++) {
      auto* g = graphs.at(sides[i]);
      xjjroot::setthgrstyle(g, side_colors[i], side_markers[i], 1.3, side_colors[i]);
      mg->Add(g, "pe");
      leg->AddEntry(g, sides[i].c_str(), "pe");
    }
    mg->SetMinimum(ymin);
    mg->SetMaximum(ymax);
    mg->Draw("ape");
    mg->GetXaxis()->SetTitle("instLumi");
    mg->GetYaxis()->SetTitle(ytitle.c_str());
    leg->Draw();
    c->SaveAs((outdir + "/" + outname).c_str());
  };

  draw_overlay(g_mean, "Gaussian mean", "mean_vs_lumi.pdf", 2000., 3400.);
  draw_overlay(g_sigma, "Gaussian sigma", "sigma_vs_lumi.pdf", 100., 1000.);

  __XJJLOG << ">> plots saved under " << outdir << "/" << std::endl;

  inf->Close();

  return 0;
}

int main(int argc, char* argv[]) {
  if (argc == 2) {
    return macro(argv[1]);
  }
  __XJJLOG << "!! usage: " << argv[0] << " <TAG (reads outputs/TAG/Main.root, writes outputs/TAG/mean_vs_lumi.pdf + sigma_vs_lumi.pdf)>" << std::endl;
  return 1;
}
