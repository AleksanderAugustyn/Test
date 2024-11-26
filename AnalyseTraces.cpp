#include <iostream>

#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TString.h>
#include <TSystem.h>
#include <TPaveText.h>

#include "main.h"

/**
 * Creates and returns a decorated scatter plot
 * @param Name Histogram name
 * @param Title Histogram title
 * @param XTitle X-axis title
 * @param YTitle Y-axis title
 * @param NBinsX Number of X bins
 * @param XMin X minimum
 * @param XMax X maximum
 * @param NBinsY Number of Y bins
 * @param YMin Y minimum
 * @param YMax Y maximum
 * @return Pointer to created histogram
 */
TH2D *CreateScatterPlot(const char *Name, const char *Title,
                        const char *XTitle, const char *YTitle,
                        const Int_t NBinsX, const Double_t XMin, const Double_t XMax,
                        const Int_t NBinsY, const Double_t YMin, const Double_t YMax)
{
    auto *Hist = new TH2D(Name, Title,
                          NBinsX, XMin, XMax,
                          NBinsY, YMin, YMax);

    Hist->GetXaxis()->SetTitle(XTitle);
    Hist->GetYaxis()->SetTitle(YTitle);
    Hist->SetMarkerStyle(20);
    Hist->SetMarkerSize(0.5);

    return Hist;
}

AnalysisHistograms InitializeHistograms(const std::vector<std::string>& Channels)
{
    AnalysisHistograms Histograms;

    for (const auto& Channel : Channels)
    {
        // Initialize scatter plots for decay and rise time vs position
        Histograms.ScatterPlots[Channel].push_back(CreateScatterPlot(
            Form("%s_decay_vs_x", Channel.c_str()),
            Form("%s Decay Time vs X Position", Channel.c_str()),
            "X Position", "Decay Time [ns]",
            500, 0.0, 0.5, 200, 0, 200));

        Histograms.ScatterPlots[Channel].push_back(CreateScatterPlot(
            Form("%s_rise_vs_x", Channel.c_str()),
            Form("%s Rise Time vs X Position", Channel.c_str()),
            "X Position", "Rise Time [ns]",
            500, 0.0, 0.5, 100, 0, 100));

        Histograms.ScatterPlots[Channel].push_back(CreateScatterPlot(
            Form("%s_decay_vs_y", Channel.c_str()),
            Form("%s Decay Time vs Y Position", Channel.c_str()),
            "Y Position", "Decay Time [ns]",
            500, 0.0, 0.5, 200, 0, 200));

        Histograms.ScatterPlots[Channel].push_back(CreateScatterPlot(
            Form("%s_rise_vs_y", Channel.c_str()),
            Form("%s Rise Time vs Y Position", Channel.c_str()),
            "Y Position", "Rise Time [ns]",
            500, 0.0, 0.5, 100, 0, 100));

        // Create corresponding profile histograms
        for (const auto* Scatter : Histograms.ScatterPlots[Channel])
        {
            auto* Profile = new TProfile(Form("%s_prof", Scatter->GetName()),
                                       Form("%s Profile", Scatter->GetTitle()),
                                       100, Scatter->GetXaxis()->GetXmin(),
                                       Scatter->GetXaxis()->GetXmax());
            Profile->GetXaxis()->SetTitle(Scatter->GetXaxis()->GetTitle());
            Profile->GetYaxis()->SetTitle(Scatter->GetYaxis()->GetTitle());
            Histograms.Profiles[Channel].push_back(Profile);
        }

        // Create count histogram - using TH2D instead of TProfile2D for actual counts
        Histograms.CountHist[Channel] = new TH2D(
            Form("%s_counts_2d", Channel.c_str()),
            Form("%s Count Distribution", Channel.c_str()),
            500, 0.0, 0.5, 500, 0.0, 0.5);
        Histograms.CountHist[Channel]->GetXaxis()->SetTitle("X Position");
        Histograms.CountHist[Channel]->GetYaxis()->SetTitle("Y Position");
        Histograms.CountHist[Channel]->GetZaxis()->SetTitle("Counts");
        Histograms.CountHist[Channel]->SetTitleSize(0.05);
        Histograms.CountHist[Channel]->GetXaxis()->SetTitleSize(0.05);
        Histograms.CountHist[Channel]->GetYaxis()->SetTitleSize(0.05);
        Histograms.CountHist[Channel]->GetZaxis()->SetTitleSize(0.05);

        // Create 2D profiles
        Histograms.Profile2Ds[Channel + "_decay"] = new TProfile2D(
            Form("%s_decay_2d_prof", Channel.c_str()),
            Form("%s Decay Time vs Position", Channel.c_str()),
            500, 0.0, 0.5, 500, 0.0, 0.5);
        Histograms.Profile2Ds[Channel + "_decay"]->GetXaxis()->SetTitle("X Position");
        Histograms.Profile2Ds[Channel + "_decay"]->GetYaxis()->SetTitle("Y Position");
        Histograms.Profile2Ds[Channel + "_decay"]->GetZaxis()->SetTitle("Decay Time [ns]");
        Histograms.Profile2Ds[Channel + "_decay"]->SetTitleSize(0.05);
        Histograms.Profile2Ds[Channel + "_decay"]->GetXaxis()->SetTitleSize(0.05);
        Histograms.Profile2Ds[Channel + "_decay"]->GetYaxis()->SetTitleSize(0.05);
        Histograms.Profile2Ds[Channel + "_decay"]->GetZaxis()->SetTitleSize(0.05);

        Histograms.Profile2Ds[Channel + "_rise"] = new TProfile2D(
            Form("%s_rise_2d_prof", Channel.c_str()),
            Form("%s Rise Time vs Position", Channel.c_str()),
            500, 0.0, 0.5, 500, 0.0, 0.5);
        Histograms.Profile2Ds[Channel + "_rise"]->GetXaxis()->SetTitle("X Position");
        Histograms.Profile2Ds[Channel + "_rise"]->GetYaxis()->SetTitle("Y Position");
        Histograms.Profile2Ds[Channel + "_rise"]->GetZaxis()->SetTitle("Rise Time [ns]");
        Histograms.Profile2Ds[Channel + "_rise"]->SetTitleSize(0.05);
        Histograms.Profile2Ds[Channel + "_rise"]->GetXaxis()->SetTitleSize(0.05);
        Histograms.Profile2Ds[Channel + "_rise"]->GetYaxis()->SetTitleSize(0.05);
        Histograms.Profile2Ds[Channel + "_rise"]->GetZaxis()->SetTitleSize(0.05);

        Histograms.Profile2Ds[Channel + "_power"] = new TProfile2D(
            Form("%s_power_2d_prof", Channel.c_str()),
            Form("%s Rise Power vs Position", Channel.c_str()),
            500, 0.0, 0.5, 500, 0.0, 0.5);
        Histograms.Profile2Ds[Channel + "_power"]->GetXaxis()->SetTitle("X Position");
        Histograms.Profile2Ds[Channel + "_power"]->GetYaxis()->SetTitle("Y Position");
        Histograms.Profile2Ds[Channel + "_power"]->GetZaxis()->SetTitle("Rise Power");
        Histograms.Profile2Ds[Channel + "_power"]->SetTitleSize(0.05);
        Histograms.Profile2Ds[Channel + "_power"]->GetXaxis()->SetTitleSize(0.05);
        Histograms.Profile2Ds[Channel + "_power"]->GetYaxis()->SetTitleSize(0.05);
        Histograms.Profile2Ds[Channel + "_power"]->GetZaxis()->SetTitleSize(0.05);

        // Create Rise Power histogram
        Histograms.RisePowerHists[Channel] = new TH1D(
            Form("%s_rise_power_hist", Channel.c_str()),
            Form("%s Rise Power Distribution;Rise Power;Counts", Channel.c_str()),
            3000, 1.0, 4.0);

        // NEW: Create Rise Power vs Position profile for filtered events
        Histograms.FilteredRisePowerScatter[Channel] = new TProfile2D(
            Form("%s_rise_power_vs_pos_filtered", Channel.c_str()),
            Form("%s Rise Power vs Position (0.1-0.4 region);X Position;Y Position", Channel.c_str()),
            300, 0.1, 0.4, 300, 0.1, 0.4);
        Histograms.FilteredRisePowerScatter[Channel]->GetZaxis()->SetTitle("Rise Power");
        Histograms.FilteredRisePowerScatter[Channel]->SetTitleSize(0.05);
        Histograms.FilteredRisePowerScatter[Channel]->GetXaxis()->SetTitleSize(0.05);
        Histograms.FilteredRisePowerScatter[Channel]->GetYaxis()->SetTitleSize(0.05);
        Histograms.FilteredRisePowerScatter[Channel]->GetZaxis()->SetTitleSize(0.05);
    }

    return Histograms;
}

void ProcessInputFiles(const std::vector<std::pair<Int_t, Int_t>>& RunsToAnalyze,
                      AnalysisHistograms& Histograms,
                      const std::vector<std::string>& Channels)
{
    for (const auto& [RunNumber, SubRunNumber] : RunsToAnalyze)
    {
        TString FileName = Form("analysis_%03d_%02d.root", RunNumber, SubRunNumber);
        std::cout << "Processing " << FileName.Data() << std::endl;

        TFile* InputFile = TFile::Open(FileName.Data());
        if (!InputFile || InputFile->IsZombie())
        {
            std::cerr << "Could not open file: " << FileName.Data() << std::endl;
            continue;
        }

        const auto Tree = dynamic_cast<TTree*>(InputFile->Get("analysis"));
        if (!Tree)
        {
            std::cerr << "Could not find analysis tree in " << FileName.Data() << std::endl;
            InputFile->Close();
            delete InputFile;
            continue;
        }

        // Set up branch variables
        Double_t PosX, PosY;
        std::map<std::string, ChannelFitData> FitData;

        // Set branch addresses
        Tree->SetBranchAddress("pos_x", &PosX);
        Tree->SetBranchAddress("pos_y", &PosY);

        for (const auto& Channel : Channels)
        {
            Tree->SetBranchAddress(Form("%s_decay_constant", Channel.c_str()),
                                 &FitData[Channel].DecayTime);
            Tree->SetBranchAddress(Form("%s_rise_time", Channel.c_str()),
                                 &FitData[Channel].RiseTime);
            Tree->SetBranchAddress(Form("%s_rise_power", Channel.c_str()),
                                 &FitData[Channel].RisePower);
        }

        // Process all entries
        const Long64_t Entries = Tree->GetEntries();
        for (Long64_t Entry = 0; Entry < Entries; Entry++)
        {
            Tree->GetEntry(Entry);

            for (const auto& Channel : Channels)
            {
                // Fill histograms
                Histograms.ScatterPlots[Channel][0]->Fill(PosX, FitData[Channel].DecayTime);
                Histograms.ScatterPlots[Channel][1]->Fill(PosX, FitData[Channel].RiseTime);
                Histograms.ScatterPlots[Channel][2]->Fill(PosY, FitData[Channel].DecayTime);
                Histograms.ScatterPlots[Channel][3]->Fill(PosY, FitData[Channel].RiseTime);

                Histograms.Profiles[Channel][0]->Fill(PosX, FitData[Channel].DecayTime);
                Histograms.Profiles[Channel][1]->Fill(PosX, FitData[Channel].RiseTime);
                Histograms.Profiles[Channel][2]->Fill(PosY, FitData[Channel].DecayTime);
                Histograms.Profiles[Channel][3]->Fill(PosY, FitData[Channel].RiseTime);

                Histograms.Profile2Ds[Channel + "_decay"]->Fill(PosX, PosY, FitData[Channel].DecayTime);
                Histograms.Profile2Ds[Channel + "_rise"]->Fill(PosX, PosY, FitData[Channel].RiseTime);
                Histograms.Profile2Ds[Channel + "_power"]->Fill(PosX, PosY, FitData[Channel].RisePower);

                // Fill count histogram with actual event counts
                Histograms.CountHist[Channel]->Fill(PosX, PosY);

                Histograms.RisePowerHists[Channel]->Fill(FitData[Channel].RisePower);

                // NEW: Fill filtered rise power scatter plot
                if (PosX >= 0.1 && PosX <= 0.4 && PosY >= 0.1 && PosY <= 0.4) {
                    Histograms.FilteredRisePowerScatter[Channel]->Fill(PosX, PosY, FitData[Channel].RisePower);
                }
            }
        }

        InputFile->Close();
        delete InputFile;
    }
}

void CreateAndSaveChannelPlots(const AnalysisHistograms& Histograms,
                              const std::vector<std::string>& Channels,
                              const char* OutputDirectory)
{
    const TString OutputFileName = Form("%s/position_fit_analysis.root", OutputDirectory);
    TFile OutputFile(OutputFileName, "RECREATE");

    // Set up style for plots
    gStyle->SetOptStat(111111);
    gStyle->SetPalette(1);
    gStyle->SetOptFit(1);
    gStyle->SetTextSize(12);

    for (const auto& Channel : Channels)
    {
        // Create 1D correlation canvas
        auto* Canvas = new TCanvas(Form("%s_canvas", Channel.c_str()),
                                 Form("%s Correlations", Channel.c_str()),
                                 1600, 1200);
        Canvas->Divide(2, 2);

        // Draw scatter plots with overlaid profiles
        for (size_t i = 0; i < Histograms.ScatterPlots.at(Channel).size(); ++i)
        {
            Canvas->cd(static_cast<Int_t>(i + 1));
            gPad->SetGridx();
            gPad->SetGridy();

            Histograms.ScatterPlots.at(Channel)[i]->Draw();
            Histograms.Profiles.at(Channel)[i]->SetLineColor(kRed);
            Histograms.Profiles.at(Channel)[i]->SetLineWidth(2);
            Histograms.Profiles.at(Channel)[i]->Draw("SAME");

            auto* Legend = new TLegend(0.65, 0.75, 0.85, 0.85);
            Legend->AddEntry(Histograms.ScatterPlots.at(Channel)[i], "Events", "p");
            Legend->AddEntry(Histograms.Profiles.at(Channel)[i], "Profile", "l");
            Legend->Draw();
        }

        Canvas->Write();
        Canvas->SaveAs(Form("%s/%s_correlations.png", OutputDirectory, Channel.c_str()));

        // Create 2D profile canvas
        auto* Canvas2D = new TCanvas(Form("%s_2d_canvas", Channel.c_str()),
                                   Form("%s Time Constants vs Position", Channel.c_str()),
                                   3000, 600);
        Canvas2D->Divide(4, 1);  // Now 4 plots instead of 3

        // Plot profiles
        Canvas2D->cd(1);
        gPad->SetGridx();
        gPad->SetGridy();
        auto* decayHist = Histograms.Profile2Ds.at(Channel + "_decay");
        decayHist->SetMinimum(0);  // Set minimum of Z axis
        decayHist->SetMaximum(200);  // Set maximum of Z axis
        decayHist->Draw("COLZ");
        gPad->SetRightMargin(0.15);

        Canvas2D->cd(2);
        gPad->SetGridx();
        gPad->SetGridy();
        auto* riseHist = Histograms.Profile2Ds.at(Channel + "_rise");
        riseHist->SetMinimum(0);  // Set minimum of Z axis
        riseHist->SetMaximum(100);  // Set maximum of Z axis
        riseHist->Draw("COLZ");
        gPad->SetRightMargin(0.15);

        Canvas2D->cd(3);
        gPad->SetGridx();
        gPad->SetGridy();
        auto* powerHist = Histograms.Profile2Ds.at(Channel + "_power");
        powerHist->SetMinimum(1.0);  // Set minimum of Z axis
        powerHist->SetMaximum(4.0);  // Set maximum of Z axis
        powerHist->Draw("COLZ");
        gPad->SetRightMargin(0.15);

        // Add counts plot
        Canvas2D->cd(4);
        gPad->SetGridx();
        gPad->SetGridy();
        gPad->SetRightMargin(0.15);

        auto* countsHist = Histograms.CountHist.at(Channel);
        countsHist->SetMinimum(0);  // Set minimum of Z axis
        countsHist->Draw("COLZ");

        // Calculate events in the 0.1-0.4 region
        const Int_t bin_0_1_x = countsHist->GetXaxis()->FindBin(0.1);
        const Int_t bin_0_4_x = countsHist->GetXaxis()->FindBin(0.4);
        const Int_t bin_0_1_y = countsHist->GetYaxis()->FindBin(0.1);
        const Int_t bin_0_4_y = countsHist->GetYaxis()->FindBin(0.4);

        Double_t totalEvents = 0;
        Double_t filteredEvents = 0;

        for (Int_t i = 1; i <= countsHist->GetNbinsX(); i++) {
            for (Int_t j = 1; j <= countsHist->GetNbinsY(); j++) {
                const Double_t binContent = countsHist->GetBinContent(i, j);
                if (binContent > 0) {  // Only count non-empty bins
                    totalEvents += binContent;
                    if (i >= bin_0_1_x && i <= bin_0_4_x &&
                        j >= bin_0_1_y && j <= bin_0_4_y) {
                        filteredEvents += binContent;
                    }
                }
            }
        }

        // Add statistics box
        const auto statsBox = new TPaveText(0.25, 0.80, 0.75, 0.90, "NDC");
        statsBox->SetFillColor(0);
        statsBox->SetBorderSize(1);
        statsBox->SetTextAlign(12);
        statsBox->SetTextSize(0.035);

        statsBox->AddText(Form("Total Events: %.0f", totalEvents));
        statsBox->AddText(Form("Events in [0.1-0.4]: %.0f (%.1f%%)",
                         filteredEvents, 100.0 * filteredEvents / totalEvents));
        statsBox->Draw();

        Canvas2D->Write();
        Canvas2D->SaveAs(Form("%s/%s_parameters_2d.png", OutputDirectory, Channel.c_str()));

        // NEW: Create filtered Rise Power vs Position plot with X projection
        auto* FilteredCanvas = new TCanvas(Form("%s_filtered_rise_power_canvas", Channel.c_str()),
                                         Form("%s Rise Power vs Position (0.1-0.4 region)", Channel.c_str()),
                                         1200, 600);
        FilteredCanvas->Divide(2, 1);  // Divide canvas into two pads

        // Left pad: 2D plot
        FilteredCanvas->cd(1);
        gPad->SetRightMargin(0.15);
        gPad->SetGridx();
        gPad->SetGridy();

        auto* filteredPowerHist = Histograms.FilteredRisePowerScatter.at(Channel);
        filteredPowerHist->SetMinimum(1.0);
        filteredPowerHist->SetMaximum(4.0);
        filteredPowerHist->Draw("COLZ");

        // Right pad: X projection
        FilteredCanvas->cd(2);
        gPad->SetGridx();
        gPad->SetGridy();

        // Create projection onto X axis using ProfileX
        auto* xProj = filteredPowerHist->ProfileX(
            Form("%s_x_proj", filteredPowerHist->GetName()));

        xProj->SetTitle(Form("%s X Projection;X Position;Average Rise Power", Channel.c_str()));
        xProj->SetStats(1);  // Enable statistics box for fit parameters
        xProj->SetMinimum(1.0);
        xProj->SetMaximum(4.0);
        xProj->SetMarkerStyle(20);
        xProj->SetMarkerSize(0.5);

        // Create and perform parabolic fit
        auto parabola = new TF1("parabola", "[0] + [1]*(x-[5]) + [2]*(x-[5])^2 + [3]*(x-[5])^3 + [4]*(x-[5])^4", 0.1, 0.4);
        parabola->SetParameters(5.0, -20.0, 25.0, -20.0, 15.0, 0.25);
        parabola->SetParNames("Offset", "Linear", "Quadratic", "Cubic", "Quartic", "Center");
        parabola->SetLineColor(kRed);
        xProj->Fit(parabola, "R");  // R = fit in specified range

        xProj->Draw("E1");
        parabola->Draw("same");

        FilteredCanvas->Write();
        FilteredCanvas->SaveAs(Form("%s/%s_rise_power_filtered.png", OutputDirectory, Channel.c_str()));
        delete parabola;
        delete FilteredCanvas;
        delete xProj;

        // Create Rise Power distribution plot
        auto* RisePowerCanvas = new TCanvas(Form("%s_rise_power_canvas", Channel.c_str()),
                                          Form("%s Rise Power Distribution", Channel.c_str()),
                                          800, 600);
        RisePowerCanvas->SetLogy();
        Histograms.RisePowerHists.at(Channel)->Draw();

        // Add rise power statistics
        const Double_t TotalEvents = Histograms.RisePowerHists.at(Channel)->GetEntries();
        auto* Hist = Histograms.RisePowerHists.at(Channel);

        const Int_t bin1_1 = Hist->FindBin(1.0);
        const Int_t bin1_2 = Hist->FindBin(1.1);
        const Int_t bin1_5 = Hist->FindBin(1.5);
        const Int_t bin2 = Hist->FindBin(2.0);
        const Int_t bin3 = Hist->FindBin(3.0);
        const Int_t bin4 = Hist->FindBin(4.0);

        auto* StatsBox = new TPaveText(0.45, 0.60, 0.70, 0.85, "NDC");
        StatsBox->SetFillColor(0);
        StatsBox->SetBorderSize(1);
        StatsBox->SetTextAlign(12);
        StatsBox->SetTextSize(0.035);

        StatsBox->AddText("Range statistics:");
        for (const auto& Range : std::vector<std::pair<std::pair<Int_t, Int_t>, const char*>>{
            {{bin1_1, bin1_2 - 1}, "1.0-1.1"},
            {{bin1_1, bin1_5 - 1}, "1.0-1.5"},
            {{bin1_1, bin2 - 1}, "1.0-2.0"},
            {{bin2, bin3 - 1}, "2.0-3.0"},
            {{bin3, bin4 - 1}, "3.0-4.0"}})
        {
            const Double_t Count = Hist->Integral(Range.first.first, Range.first.second);
            StatsBox->AddText(Form("%s: %.0f (%.1f%%)", Range.second,
                                 Count, 100.0 * Count / TotalEvents));
        }

        StatsBox->Draw();
        RisePowerCanvas->Write();
        RisePowerCanvas->SaveAs(Form("%s/%s_rise_power_dist.png", OutputDirectory, Channel.c_str()));

        // Clean up canvases
        delete StatsBox;
        delete RisePowerCanvas;
        delete Canvas;
        delete Canvas2D;

        // Write histograms
        for (auto* Histo : Histograms.ScatterPlots.at(Channel)) { Histo->Write(); }
        for (auto* Prof : Histograms.Profiles.at(Channel)) { Prof->Write(); }
        for (const auto& ProfileType : {"_decay", "_rise", "_power"})
        {
            Histograms.Profile2Ds.at(Channel + ProfileType)->Write();
        }
        Histograms.RisePowerHists.at(Channel)->Write();
        Histograms.FilteredRisePowerScatter.at(Channel)->Write();  // NEW: Write filtered histogram
    }

    OutputFile.Close();
}

void CleanupHistograms(AnalysisHistograms& Histograms)
{
    for (auto& [Channel, Hists] : Histograms.ScatterPlots)
    {
        for (const auto* Hist : Hists) { delete Hist; }
    }
    for (auto& [Channel, Profs] : Histograms.Profiles)
    {
        for (const auto* Prof : Profs) { delete Prof; }
    }
    for (auto& [Key, Prof] : Histograms.Profile2Ds)
    {
        delete Prof;
    }
    for (auto& [Channel, Hist] : Histograms.RisePowerHists)
    {
        delete Hist;
    }
    for (auto& [Channel, Hist] : Histograms.CountHist)
    {
        delete Hist;
    }
    // NEW: Cleanup filtered rise power scatter plots
    for (auto& [Channel, Hist] : Histograms.FilteredRisePowerScatter)
    {
        delete Hist;
    }
}

// Main analysis function that uses all the above components
void AnalyzePositionVsFitParameters(const std::vector<std::pair<Int_t, Int_t>>& RunsToAnalyze,
                                  const char* OutputDirectory)
{
    std::cout << "\n[AnalyzePositionVsFitParameters] Starting analysis of "
              << RunsToAnalyze.size() << " runs" << std::endl;

    // Create output directory if it doesn't exist
    gSystem->mkdir(OutputDirectory, kTRUE);

    // Define channels to analyze
    const std::vector<std::string> Channels = {"xa", "xb", "ya", "yb"};

    // Initialize histograms
    auto Histograms = InitializeHistograms(Channels);

    // Process all input files
    ProcessInputFiles(RunsToAnalyze, Histograms, Channels);

    // Create and save plots
    CreateAndSaveChannelPlots(Histograms, Channels, OutputDirectory);

    // Clean up
    CleanupHistograms(Histograms);
}