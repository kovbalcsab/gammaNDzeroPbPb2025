// plotRunLumi: occupancy of (Run, Lumi section) for the selected events, to hunt for anomalous run/lumisections.
//
// Input: ntuple from ../analyze.cc (needs the Run and Lumi branches).
// Selections (the unconstrained "An" side carries the ZDC window, the other side must be 0n):
//   0nAn: ZDCsumMinus < ZDC_0N_CUT  and  ZDC_0NAN_MIN < ZDCsumPlus  < ZDC_0NAN_MAX
//   An0n: ZDCsumPlus  < ZDC_0N_CUT  and  ZDC_AN0N_MIN < ZDCsumMinus < ZDC_AN0N_MAX
// x axis = every distinct (Run, Lumi) pair seen in either selection, sorted by (Run, Lumi), one bin each
// (so equal lumi sections of different runs never overlap); y axis = lumi bin of params.h; colour = N events.
//
// Output (outputs/TAG_runLumi/): RunLumi.pdf (0nAn on top, An0n below, dashed lines + labels at run boundaries), RunLumi.root,
//   ZDCAn_perRun/ZDCAn_run<RUN>.pdf: one plot per run of the An-side ZDC sum (0n on the other side, no An window; red = window),
//   with the lumi bins of params.h overlaid, each normalised to unit area

#include <filesystem>
#include <algorithm>
#include <utility>
#include <vector>
#include <map>
#include <set>
#include <array>

#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>

#include "xjjanauti.h"
#include "../params.h"

namespace {
  const double ZDC_0N_CUT = 1500.;
  const double ZDC_0NAN_MIN = 400., ZDC_0NAN_MAX = 800.;   // ZDCsumPlus window for 0nAn
  const double ZDC_AN0N_MIN = 500., ZDC_AN0N_MAX = 1000.;  // ZDCsumMinus window for An0n
  const int ZDC_AN_NBINS = 150;                            // per-run An-side ZDC sum histograms
  const double ZDC_AN_XMIN = 0., ZDC_AN_XMAX = 3000.;
  const double MIN_ENTRIES_OVERLAY = 100.;                 // lumi bins with fewer entries in a run are not overlaid

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

  struct Hit { int run, lumi, lb; };
}

int macro(const std::string& inputname, const std::string& tag, int leveled = -1) {
  // leveled: -1 = no cut on lumiLeveled (default), 0 = only non-leveled events, 1 = only leveled events
  namespace fs = std::filesystem;

  const std::string outdir = "outputs/" + tag + (leveled == 1 ? "_leveled" : leveled == 0 ? "_notleveled" : "") + "_runLumi";
  fs::create_directories(outdir);

  auto* inf = TFile::Open(inputname.c_str());
  if (!inf || inf->IsZombie()) { __XJJLOG << "!! failed to open input file, abort." << std::endl; return 1; }
  auto* nt = (TTree*)inf->Get("ntuple");
  if (!nt) { __XJJLOG << "!! no ntuple found in input file, abort." << std::endl; return 2; }
  if (!nt->GetBranch("Run") || !nt->GetBranch("Lumi")) {
    __XJJLOG << "!! ntuple has no Run/Lumi branches, rerun ../analyze.exe, abort." << std::endl;
    return 3;
  }

  Float_t instLumi = 0., ZDCsumPlus = 0., ZDCsumMinus = 0., lumiLeveled = 0., Run = 0., Lumi = 0.;
  nt->SetBranchStatus("*", 0);
  for (const char* b : { "instLumi", "ZDCsumPlus", "ZDCsumMinus", "Run", "Lumi" }) nt->SetBranchStatus(b, 1);
  nt->SetBranchAddress("instLumi", &instLumi);
  nt->SetBranchAddress("ZDCsumPlus", &ZDCsumPlus);
  nt->SetBranchAddress("ZDCsumMinus", &ZDCsumMinus);
  nt->SetBranchAddress("Run", &Run);
  nt->SetBranchAddress("Lumi", &Lumi);
  if (leveled >= 0) {
    nt->SetBranchStatus("lumiLeveled", 1);
    nt->SetBranchAddress("lumiLeveled", &lumiLeveled);
  }

  const auto nentries = nt->GetEntries();
  __XJJLOG << ">> " << nentries << " entries loaded from " << inputname << std::endl;

  const int nlumibins = (int)params::lumibins.size() - 1;
  const std::array<std::string, 2> cuts = { "0nAn", "An0n" };
  std::array<std::vector<Hit>, 2> hits;
  std::set<std::pair<int, int>> keys; // distinct (run, lumi), sorted
  // run -> [0nAn: ZDCsumPlus, An0n: ZDCsumMinus] -> one An-side ZDC sum histogram per lumi bin
  std::map<int, std::array<std::vector<TH1F*>, 2>> anSum;

  for (Long64_t i = 0; i < nentries; ++i) {
    nt->GetEntry(i);
    if (i % 1000000 == 0) __XJJLOG << ">> processing entry " << i << " / " << nentries << std::endl;
    if (leveled >= 0 && (int)lumiLeveled != leveled) continue;
    const int lb = find_lumibin(instLumi);
    if (lb < 0) continue;
    // per-run An-side ZDC sum: 0n required on the opposite side, no window on the An side
    auto& hrun = anSum[(int)Run];
    if (hrun[0].empty()) {
      for (int c = 0; c < 2; ++c) {
        for (int b = 0; b < nlumibins; ++b) {
          auto* hh = new TH1F(Form("hZDCAn_%s_run%d_lumibin%d", cuts[c].c_str(), (int)Run, b), "", ZDC_AN_NBINS, ZDC_AN_XMIN, ZDC_AN_XMAX);
          hh->SetDirectory(nullptr);
          hrun[c].push_back(hh);
        }
      }
    }
    if (ZDCsumMinus < ZDC_0N_CUT) hrun[0][lb]->Fill(ZDCsumPlus);
    if (ZDCsumPlus < ZDC_0N_CUT) hrun[1][lb]->Fill(ZDCsumMinus);

    const bool is0nAn = ZDCsumMinus < ZDC_0N_CUT && ZDCsumPlus > ZDC_0NAN_MIN && ZDCsumPlus < ZDC_0NAN_MAX;
    const bool isAn0n = ZDCsumPlus < ZDC_0N_CUT && ZDCsumMinus > ZDC_AN0N_MIN && ZDCsumMinus < ZDC_AN0N_MAX;
    if (!is0nAn && !isAn0n) continue;
    const Hit h = { (int)Run, (int)Lumi, lb };
    keys.insert({ h.run, h.lumi });
    if (is0nAn) hits[0].push_back(h);
    if (isAn0n) hits[1].push_back(h);
  }
  __XJJLOG << ">> selected: 0nAn " << hits[0].size() << ", An0n " << hits[1].size() << " events in " << keys.size() << " (run,lumi) pairs" << std::endl;
  if (keys.empty()) { __XJJLOG << "!! no events selected, abort." << std::endl; return 4; }

  // x bin of each (run,lumi); first x bin of each run
  std::map<std::pair<int, int>, int> xbin;
  std::vector<std::pair<int, int>> runStart; // (run, first x index)
  for (const auto& k : keys) {
    const int ix = (int)xbin.size();
    xbin[k] = ix;
    if (runStart.empty() || runStart.back().first != k.first) runStart.emplace_back(k.first, ix);
  }
  const int nx = (int)keys.size();

  std::array<TH2F*, 2> h2;
  for (int c = 0; c < 2; ++c) {
    h2[c] = new TH2F(("hRunLumi_" + cuts[c]).c_str(), "", nx, 0, nx, nlumibins, 0, nlumibins);
    for (const auto& h : hits[c]) h2[c]->Fill(xbin[{ h.run, h.lumi }] + 0.5, h.lb + 0.5);
  }

  gStyle->SetOptStat(0);
  gStyle->SetPalette(kViridis);
  TCanvas cv("cv", "cv", 2800, 1400);
  cv.Divide(1, 2, 0., 0.);
  for (int c = 0; c < 2; ++c) {
    cv.cd(c + 1);
    gPad->SetLeftMargin(0.05); gPad->SetRightMargin(0.07);
    gPad->SetTopMargin(c == 0 ? 0.10 : 0.); gPad->SetBottomMargin(c == 0 ? 0. : 0.16);
    auto* h = h2[c];
    h->SetMinimum(0.5);
    h->GetYaxis()->SetTitle("lumi bin");
    h->GetYaxis()->SetNdivisions(nlumibins, false);
    h->GetXaxis()->SetLabelSize(0);
    h->GetXaxis()->SetTickLength(0);
    h->GetYaxis()->SetTitleSize(0.05); h->GetYaxis()->SetTitleOffset(0.4); h->GetYaxis()->SetLabelSize(0.045);
    if (c == 1) h->GetXaxis()->SetTitle("(Run, Lumi) pairs in time order, grouped by run");
    h->Draw("COLZ");
    gPad->Update();

    const double ymax = nlumibins;
    TLatex t; t.SetNDC(); t.SetTextSize(0.05);
    const std::string win = c == 0 ? Form("%g<ZDCsumPlus<%g, ZDCsumMinus<%g", ZDC_0NAN_MIN, ZDC_0NAN_MAX, ZDC_0N_CUT)
                                   : Form("%g<ZDCsumMinus<%g, ZDCsumPlus<%g", ZDC_AN0N_MIN, ZDC_AN0N_MAX, ZDC_0N_CUT);
    t.DrawLatex(0.06, c == 0 ? 0.92 : 0.86, Form("%s: %s, N = %zu", cuts[c].c_str(), win.c_str(), hits[c].size()));
    for (const auto& rs : runStart) {
      TLine l; l.SetLineStyle(3); l.SetLineColor(kGray + 1);
      l.DrawLine(rs.second, 0, rs.second, ymax);
    }
    if (c == 1) { // run labels under the x axis
      TLatex r; r.SetTextAngle(90); r.SetTextSize(0.028); r.SetTextAlign(32);
      const double yl = -0.02 * ymax;
      for (const auto& rs : runStart) r.DrawLatex(rs.second + 0.3 * nx / 2000., yl, Form("%d", rs.first));
    }
  }
  cv.SaveAs((outdir + "/RunLumi.pdf").c_str());

  // one plot per run: An-side ZDC sum for 0nAn (left) and An0n (right)
  const std::string rundir = outdir + "/ZDCAn_perRun";
  fs::create_directories(rundir);
  const std::array<double, 2> wlo = { ZDC_0NAN_MIN, ZDC_AN0N_MIN }, whi = { ZDC_0NAN_MAX, ZDC_AN0N_MAX };
  const std::array<std::string, 2> xt = { "ZDCsumPlus", "ZDCsumMinus" };
  std::vector<TH1F*> keep;
  for (const auto& [run, hs] : anSum) {
    TCanvas cr(Form("cr_%d", run), "cr", 1600, 700);
    cr.Divide(2, 1);
    for (int c = 0; c < 2; ++c) {
      cr.cd(c + 1);
      gPad->SetLogy();
      gPad->SetLeftMargin(0.13); gPad->SetRightMargin(0.04); gPad->SetTopMargin(0.10);
      // overlay of the lumi bins, each normalised to unit area (shape comparison); empty bins skipped
      TLegend leg(0.62, 0.45, 0.95, 0.89);
      leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.03);
      TH1F* first = nullptr;
      double ymax = 0.;
      std::vector<TH1F*> drawn;
      for (int b = 0; b < nlumibins; ++b) {
        const TH1F* raw = hs[c][b];
        if (raw->GetEntries() < MIN_ENTRIES_OVERLAY) continue;
        auto* h = (TH1F*)raw->Clone(Form("%s_norm", raw->GetName()));
        h->SetDirectory(nullptr);
        h->Scale(1. / h->Integral());
        const int col = gStyle->GetColorPalette(nlumibins > 1 ? b * (gStyle->GetNumberOfColors() - 1) / (nlumibins - 1) : 0);
        h->SetLineColor(col); h->SetLineWidth(2);
        ymax = std::max(ymax, h->GetMaximum());
        leg.AddEntry(h, Form("lumi %g-%g (N=%.0f)", params::lumibins[b], params::lumibins[b + 1], raw->GetEntries()), "l");
        drawn.push_back(h);
      }
      if (drawn.empty()) continue;
      for (auto* h : drawn) {
        if (!first) {
          first = h;
          h->GetXaxis()->SetTitle((xt[c] + " (An side)").c_str());
          h->GetYaxis()->SetTitle("fraction of events");
          h->SetMaximum(ymax * 3.);
          h->SetMinimum(ymax * 1e-5);
          h->Draw("HIST");
        } else h->Draw("HIST SAME");
      }
      leg.Draw();
      for (double x : { wlo[c], whi[c] }) {
        TLine l; l.SetLineStyle(2); l.SetLineColor(kRed + 1);
        l.DrawLine(x, first->GetMinimum(), x, ymax);
      }
      TLatex t; t.SetNDC(); t.SetTextSize(0.045);
      t.DrawLatex(0.14, 0.92, Form("Run %d, %s", run, cuts[c].c_str()));
      gPad->Modified();
      cr.cd(c + 1)->Update();
      // histograms must outlive the SaveAs below
      for (auto* h : drawn) keep.push_back(h);
    }
    cr.SaveAs(Form("%s/ZDCAn_run%d.pdf", rundir.c_str(), run));
    for (auto* h : keep) delete h;
    keep.clear();
  }
  __XJJLOG << ">> " << anSum.size() << " per-run plots written to " << rundir << std::endl;

  auto* outf = TFile::Open((outdir + "/RunLumi.root").c_str(), "RECREATE");
  for (auto* h : h2) h->Write();
  for (auto& [run, hs] : anSum) for (auto& v : hs) for (auto* h : v) h->Write();
  outf->Close();
  inf->Close();
  __XJJLOG << ">> plot written to " << outdir << "/RunLumi.pdf" << std::endl;
  return 0;
}

int main(int argc, char* argv[]) {
  if (argc >= 3 && argc <= 4) return macro(argv[1], argv[2], argc == 4 ? std::atoi(argv[3]) : -1);
  std::cerr << "usage: " << argv[0] << " <ntuple.root (output of ../analyze.exe)> <tag> [leveled=-1 (-1 all, 0 not leveled, 1 leveled)]\n";
  return 1;
}
