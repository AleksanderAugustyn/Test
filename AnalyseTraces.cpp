#include <iostream>

#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TString.h>
#include <TSystem.h>

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
TH2D* CreateScatterPlot(const char* Name, const char* Title,
                        const char* XTitle, const char* YTitle,
                        const Int_t NBinsX, const Double_t XMin, const Double_t XMax,
                        const Int_t NBinsY, const Double_t YMin, const Double_t YMax)
{
    auto* Hist = new TH2D(Name, Title,
                         NBinsX, XMin, XMax,
                         NBinsY, YMin, YMax);

    Hist->GetXaxis()->SetTitle(XTitle);
    Hist->GetYaxis()->SetTitle(YTitle);
    Hist->SetMarkerStyle(20);
    Hist->SetMarkerSize(0.5);

    return Hist;
}

/**
 * Analyzes the relationship between position and fit parameters across multiple runs
 * @param RunsToAnalyze Vector of run numbers to analyze
 * @param OutputDirectory Directory to save plots
 */
void AnalyzePositionVsFitParameters(const std::vector<std::pair<Int_t, Int_t>>& RunsToAnalyze,
                                  const char* OutputDirectory)
{
    std::cout << "\n[AnalyzePositionVsFitParameters] Starting analysis of "
              << RunsToAnalyze.size() << " runs" << std::endl;

    // Create output directory if it doesn't exist
    gSystem->mkdir(OutputDirectory, kTRUE);

    // Set up histograms for each anode channel
    const std::vector<std::string> Channels = {"xa", "xb", "ya", "yb"};

    // Map to store histograms for each channel
    std::map<std::string, std::vector<TH2D*>> ScatterPlots;
    std::map<std::string, std::vector<TProfile*>> Profiles;
    std::map<std::string, TProfile2D*> Profile2Ds;

    // Initialize histograms for each channel
    for (const auto& Channel : Channels)
    {
        // Decay time (p[2]) vs X position
        ScatterPlots[Channel].push_back(CreateScatterPlot(
            Form("%s_decay_vs_x", Channel.c_str()),
            Form("%s Decay Time vs X Position", Channel.c_str()),
            "X Position", "Decay Time [ns]",
            100, 0.0, 0.5, 100, 0, 200));

        // Rise time (p[3]) vs X position
        ScatterPlots[Channel].push_back(CreateScatterPlot(
            Form("%s_rise_vs_x", Channel.c_str()),
            Form("%s Rise Time vs X Position", Channel.c_str()),
            "X Position", "Rise Time [ns]",
            100, 0.0, 0.5, 100, 0, 100));

        // Decay time vs Y position
        ScatterPlots[Channel].push_back(CreateScatterPlot(
            Form("%s_decay_vs_y", Channel.c_str()),
            Form("%s Decay Time vs Y Position", Channel.c_str()),
            "Y Position", "Decay Time [ns]",
            100, 0.0, 0.5, 100, 0, 200));

        // Rise time vs Y position
        ScatterPlots[Channel].push_back(CreateScatterPlot(
            Form("%s_rise_vs_y", Channel.c_str()),
            Form("%s Rise Time vs Y Position", Channel.c_str()),
            "Y Position", "Rise Time [ns]",
            100, 0.0, 0.5, 100, 0, 100));

        // Create corresponding profile histograms
        for (const auto* Scatter : ScatterPlots[Channel])
        {
            auto* Profile = new TProfile(Form("%s_prof", Scatter->GetName()),
                                       Form("%s Profile", Scatter->GetTitle()),
                                       100, Scatter->GetXaxis()->GetXmin(),
                                       Scatter->GetXaxis()->GetXmax());
            Profile->GetXaxis()->SetTitle(Scatter->GetXaxis()->GetTitle());
            Profile->GetYaxis()->SetTitle(Scatter->GetYaxis()->GetTitle());
            Profiles[Channel].push_back(Profile);
        }

        // Create 2D profile of decay time vs position
        Profile2Ds[Channel] = new TProfile2D(
            Form("%s_decay_2d_prof", Channel.c_str()),
            Form("%s Decay Time vs Position;X Position;Y Position;Decay Time [ns]", Channel.c_str()),
            50, 0.0, 0.5, 50, 0.0, 0.5);
    }

    // Process each run
    for (const auto& [RunNumber, SubRunNumber] : RunsToAnalyze)
    {
        // Construct the analysis file name
        TString FileName = Form("analysis_%03d_%02d.root", RunNumber, SubRunNumber);

        std::cout << "Processing " << FileName.Data() << std::endl;

        TFile* InputFile = TFile::Open(FileName.Data());
        if (!InputFile || InputFile->IsZombie())
        {
            std::cerr << "Could not open file: " << FileName.Data() << std::endl;
            continue;
        }

        auto Tree = dynamic_cast<TTree*>(InputFile->Get("analysis"));
        if (!Tree)
        {
            std::cerr << "Could not find analysis tree in " << FileName.Data() << std::endl;
            InputFile->Close();
            delete InputFile;
            continue;
        }

        // Set up branch variables
        Double_t PosX, PosY;
        struct ChannelFitData
        {
            Double_t DecayTime;
            Double_t RiseTime;
        };
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
        }

        // Process all entries
        const Long64_t Entries = Tree->GetEntries();
        for (Long64_t Entry = 0; Entry < Entries; Entry++)
        {
            Tree->GetEntry(Entry);

            for (const auto& Channel : Channels)
            {
                // Fill scatter plots and profiles
                ScatterPlots[Channel][0]->Fill(PosX, FitData[Channel].DecayTime);
                ScatterPlots[Channel][1]->Fill(PosX, FitData[Channel].RiseTime);
                ScatterPlots[Channel][2]->Fill(PosY, FitData[Channel].DecayTime);
                ScatterPlots[Channel][3]->Fill(PosY, FitData[Channel].RiseTime);

                // Fill profile histograms
                Profiles[Channel][0]->Fill(PosX, FitData[Channel].DecayTime);
                Profiles[Channel][1]->Fill(PosX, FitData[Channel].RiseTime);
                Profiles[Channel][2]->Fill(PosY, FitData[Channel].DecayTime);
                Profiles[Channel][3]->Fill(PosY, FitData[Channel].RiseTime);

                // Fill 2D profile
                Profile2Ds[Channel]->Fill(PosX, PosY, FitData[Channel].DecayTime);
            }
        }

        InputFile->Close();
        delete InputFile;
    }

    // Create output file
    TString OutputFileName = Form("%s/position_fit_analysis.root", OutputDirectory);
    TFile OutputFile(OutputFileName, "RECREATE");

    // Set up style for plots
    gStyle->SetOptStat(111111);
    gStyle->SetPalette(1);
    gStyle->SetOptFit(1);

    // Create canvases and save plots
    for (const auto& Channel : Channels)
    {
        // Create canvas for 1D correlations
        auto* Canvas = new TCanvas(Form("%s_canvas", Channel.c_str()),
                                 Form("%s Correlations", Channel.c_str()),
                                 1600, 1200);
        Canvas->Divide(2, 2);

        // Draw scatter plots with overlaid profiles
        for (size_t i = 0; i < ScatterPlots[Channel].size(); ++i)
        {
            Canvas->cd(static_cast<Int_t>(i + 1));
            gPad->SetGridx();
            gPad->SetGridy();

            ScatterPlots[Channel][i]->Draw();
            Profiles[Channel][i]->SetLineColor(kRed);
            Profiles[Channel][i]->SetLineWidth(2);
            Profiles[Channel][i]->Draw("SAME");

            auto* Legend = new TLegend(0.65, 0.75, 0.85, 0.85);
            Legend->AddEntry(ScatterPlots[Channel][i], "Events", "p");
            Legend->AddEntry(Profiles[Channel][i], "Profile", "l");
            Legend->Draw();
        }

        // Save canvas
        Canvas->Write();
        Canvas->SaveAs(Form("%s/%s_correlations.png", OutputDirectory, Channel.c_str()));

        // Create and save 2D profile canvas
        auto* Canvas2D = new TCanvas(Form("%s_2d_canvas", Channel.c_str()),
                                   Form("%s 2D Profile", Channel.c_str()),
                                   800, 600);
        Profile2Ds[Channel]->Draw("COLZ");
        Canvas2D->Write();
        Canvas2D->SaveAs(Form("%s/%s_2d_profile.png", OutputDirectory, Channel.c_str()));

        // Clean up
        delete Canvas;
        delete Canvas2D;

        // Write histograms
        for (auto* Hist : ScatterPlots[Channel]) { Hist->Write(); }
        for (auto* Prof : Profiles[Channel]) { Prof->Write(); }
        Profile2Ds[Channel]->Write();
    }

    OutputFile.Close();

    // Clean up histograms
    for (auto& [Channel, Hists] : ScatterPlots)
    {
        for (auto* Hist : Hists) { delete Hist; }
    }
    for (auto& [Channel, Profs] : Profiles)
    {
        for (auto* Prof : Profs) { delete Prof; }
    }
    for (auto& [Channel, Prof] : Profile2Ds)
    {
        delete Prof;
    }
}