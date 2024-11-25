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

/**
 * Analyzes the relationship between position and fit parameters across multiple runs
 * @param RunsToAnalyze Vector of run numbers to analyze
 * @param OutputDirectory Directory to save plots
 */
void AnalyzePositionVsFitParameters(const std::vector<std::pair<Int_t, Int_t> > &RunsToAnalyze,
                                    const char *OutputDirectory)
{
    std::cout << "\n[AnalyzePositionVsFitParameters] Starting analysis of "
            << RunsToAnalyze.size() << " runs" << std::endl;

    // Create output directory if it doesn't exist
    gSystem->mkdir(OutputDirectory, kTRUE);

    // Set up histograms for each anode channel
    const std::vector<std::string> Channels = {"xa", "xb", "ya", "yb"};

    // Map to store histograms for each channel
    std::map<std::string, std::vector<TH2D *> > ScatterPlots;
    std::map<std::string, std::vector<TProfile *> > Profiles;
    std::map<std::string, TProfile2D *> Profile2Ds;
    std::map<std::string, TH1D *> RisePowerHists;

    // Initialize histograms for each channel
    for (const auto &Channel: Channels)
    {
        // Decay time (p[2]) vs X position
        ScatterPlots[Channel].push_back(CreateScatterPlot(
            Form("%s_decay_vs_x", Channel.c_str()),
            Form("%s Decay Time vs X Position", Channel.c_str()),
            "X Position", "Decay Time [ns]",
            500, 0.0, 0.5, 200, 0, 200));

        // Rise time (p[3]) vs X position
        ScatterPlots[Channel].push_back(CreateScatterPlot(
            Form("%s_rise_vs_x", Channel.c_str()),
            Form("%s Rise Time vs X Position", Channel.c_str()),
            "X Position", "Rise Time [ns]",
            500, 0.0, 0.5, 100, 0, 100));

        // Decay time vs Y position
        ScatterPlots[Channel].push_back(CreateScatterPlot(
            Form("%s_decay_vs_y", Channel.c_str()),
            Form("%s Decay Time vs Y Position", Channel.c_str()),
            "Y Position", "Decay Time [ns]",
            500, 0.0, 0.5, 200, 0, 200));

        // Rise time vs Y position
        ScatterPlots[Channel].push_back(CreateScatterPlot(
            Form("%s_rise_vs_y", Channel.c_str()),
            Form("%s Rise Time vs Y Position", Channel.c_str()),
            "Y Position", "Rise Time [ns]",
            500, 0.0, 0.5, 100, 0, 100));

        // Create corresponding profile histograms
        for (const auto *Scatter: ScatterPlots[Channel])
        {
            auto *Profile = new TProfile(Form("%s_prof", Scatter->GetName()),
                                         Form("%s Profile", Scatter->GetTitle()),
                                         100, Scatter->GetXaxis()->GetXmin(),
                                         Scatter->GetXaxis()->GetXmax());
            Profile->GetXaxis()->SetTitle(Scatter->GetXaxis()->GetTitle());
            Profile->GetYaxis()->SetTitle(Scatter->GetYaxis()->GetTitle());
            Profiles[Channel].push_back(Profile);
        }

        // Create 2D profiles
        Profile2Ds[Channel + "_decay"] = new TProfile2D(
            Form("%s_decay_2d_prof", Channel.c_str()),
            Form("%s Decay Time vs Position;X Position;Y Position;Decay Time [ns]", Channel.c_str()),
            500, 0.0, 0.5, 500, 0.0, 0.5);

        Profile2Ds[Channel + "_rise"] = new TProfile2D(
            Form("%s_rise_2d_prof", Channel.c_str()),
            Form("%s Rise Time vs Position;X Position;Y Position;Rise Time [ns]", Channel.c_str()),
            500, 0.0, 0.5, 500, 0.0, 0.5);

        Profile2Ds[Channel + "_power"] = new TProfile2D(
            Form("%s_power_2d_prof", Channel.c_str()),
            Form("%s Rise Power vs Position;X Position;Y Position;Rise Power", Channel.c_str()),
            500, 0.0, 0.5, 500, 0.0, 0.5);

        // Create Rise Power histogram
        RisePowerHists[Channel] = new TH1D(
            Form("%s_rise_power_hist", Channel.c_str()),
            Form("%s Rise Power Distribution;Rise Power;Counts", Channel.c_str()),
            4000, 1.0, 4.0); // Range based on typical values
    }

    // Process each run
    for (const auto &[RunNumber, SubRunNumber]: RunsToAnalyze)
    {
        // Construct the analysis file name
        TString FileName = Form("analysis_%03d_%02d.root", RunNumber, SubRunNumber);

        std::cout << "Processing " << FileName.Data() << std::endl;

        TFile *InputFile = TFile::Open(FileName.Data());
        if (!InputFile || InputFile->IsZombie())
        {
            std::cerr << "Could not open file: " << FileName.Data() << std::endl;
            continue;
        }

        auto Tree = dynamic_cast<TTree *>(InputFile->Get("analysis"));
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
            Double_t RisePower;
        };
        std::map<std::string, ChannelFitData> FitData;

        // Set branch addresses
        Tree->SetBranchAddress("pos_x", &PosX);
        Tree->SetBranchAddress("pos_y", &PosY);

        for (const auto &Channel: Channels)
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

            for (const auto &Channel: Channels)
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

                // Fill 2D profiles
                Profile2Ds[Channel + "_decay"]->Fill(PosX, PosY, FitData[Channel].DecayTime);
                Profile2Ds[Channel + "_rise"]->Fill(PosX, PosY, FitData[Channel].RiseTime);
                Profile2Ds[Channel + "_power"]->Fill(PosX, PosY, FitData[Channel].RisePower);

                // Fill Rise Power histogram
                RisePowerHists[Channel]->Fill(FitData[Channel].RisePower);
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
    for (const auto &Channel: Channels)
    {
        // Create canvas for 1D correlations
        auto *Canvas = new TCanvas(Form("%s_canvas", Channel.c_str()),
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

            auto *Legend = new TLegend(0.65, 0.75, 0.85, 0.85);
            Legend->AddEntry(ScatterPlots[Channel][i], "Events", "p");
            Legend->AddEntry(Profiles[Channel][i], "Profile", "l");
            Legend->Draw();
        }

        // Save canvas
        Canvas->Write();
        Canvas->SaveAs(Form("%s/%s_correlations.png", OutputDirectory, Channel.c_str()));

        // Create and save 2D profile canvases
        auto *Canvas2D = new TCanvas(Form("%s_2d_canvas", Channel.c_str()),
                                     Form("%s Time Constants vs Position", Channel.c_str()),
                                     2400, 600);
        Canvas2D->Divide(3, 1);

        // Plot decay time profile
        Canvas2D->cd(1);
        gPad->SetGridx();
        gPad->SetGridy();
        Profile2Ds[Channel + "_decay"]->Draw("COLZ");
        gPad->SetRightMargin(0.15);

        // Plot rise time profile
        Canvas2D->cd(2);
        gPad->SetGridx();
        gPad->SetGridy();
        Profile2Ds[Channel + "_rise"]->Draw("COLZ");
        gPad->SetRightMargin(0.15);

        // Plot rise power profile
        Canvas2D->cd(3);
        gPad->SetGridx();
        gPad->SetGridy();
        Profile2Ds[Channel + "_power"]->Draw("COLZ");
        gPad->SetRightMargin(0.15);

        // Save canvas
        Canvas2D->Write();
        Canvas2D->SaveAs(Form("%s/%s_parameters_2d.png", OutputDirectory, Channel.c_str()));

        // Create and save Rise Power histogram
        auto *RisePowerCanvas = new TCanvas(Form("%s_rise_power_canvas", Channel.c_str()),
                                            Form("%s Rise Power Distribution", Channel.c_str()),
                                            800, 600);
        RisePowerCanvas->SetLogy(); // Use log scale for y-axis
        RisePowerHists[Channel]->Draw();

        // Calculate statistics for specific ranges
        const Double_t TotalEvents = RisePowerHists[Channel]->GetEntries();
        const Int_t bin1_1 = RisePowerHists[Channel]->FindBin(1.0);
        const Int_t bin1_2 = RisePowerHists[Channel]->FindBin(1.1);
        const Int_t bin1_1_5 = RisePowerHists[Channel]->FindBin(1.5);
        const Int_t bin2 = RisePowerHists[Channel]->FindBin(2.0);
        const Int_t bin3 = RisePowerHists[Channel]->FindBin(3.0);
        const Int_t bin4 = RisePowerHists[Channel]->FindBin(4.0);

        const Double_t Count1_1_1 = RisePowerHists[Channel]->Integral(bin1_1, bin1_2 - 1);
        const Double_t Count1_1_5 = RisePowerHists[Channel]->Integral(bin1_1, bin1_1_5 - 1);
        const Double_t Count1_2 = RisePowerHists[Channel]->Integral(bin1_1, bin2 - 1);
        const Double_t Count2_3 = RisePowerHists[Channel]->Integral(bin2, bin3 - 1);
        const Double_t Count3_4 = RisePowerHists[Channel]->Integral(bin3, bin4 - 1);

        // Create statistics text box
        auto StatsBox = new TPaveText(0.45, 0.60, 0.70, 0.85, "NDC");
        StatsBox->SetFillColor(0);
        StatsBox->SetBorderSize(1);
        StatsBox->SetTextAlign(12); // Left align
        StatsBox->SetTextSize(0.035);

        StatsBox->AddText("Range statistics:");
        StatsBox->AddText(Form("1.0-1.1: %.0f (%.1f%%)",
                               Count1_1_1, 100.0 * Count1_1_1 / TotalEvents));
        StatsBox->AddText(Form("1.0-1.5: %.0f (%.1f%%)",
                               Count1_1_5, 100.0 * Count1_1_5 / TotalEvents));
        StatsBox->AddText(Form("1.0-2.0: %.0f (%.1f%%)",
                               Count1_2, 100.0 * Count1_2 / TotalEvents));
        StatsBox->AddText(Form("2.0-3.0: %.0f (%.1f%%)",
                               Count2_3, 100.0 * Count2_3 / TotalEvents));
        StatsBox->AddText(Form("3.0-4.0: %.0f (%.1f%%)",
                               Count3_4, 100.0 * Count3_4 / TotalEvents));

        StatsBox->Draw();

        RisePowerCanvas->Write();
        RisePowerCanvas->SaveAs(Form("%s/%s_rise_power_dist.png", OutputDirectory, Channel.c_str()));

        // Clean up
        delete StatsBox;
        delete RisePowerCanvas;

        // Clean up
        delete Canvas;
        delete Canvas2D;

        // Write histograms
        for (auto *Hist: ScatterPlots[Channel]) { Hist->Write(); }
        for (auto *Prof: Profiles[Channel]) { Prof->Write(); }
        for (const auto &ProfileType: {"_decay", "_rise", "_power"})
        {
            Profile2Ds[Channel + ProfileType]->Write();
        }
        RisePowerHists[Channel]->Write();
    }

    OutputFile.Close();

    // Clean up histograms
    for (auto &[Channel, Hists]: ScatterPlots)
    {
        for (auto *Hist: Hists) { delete Hist; }
    }
    for (auto &[Channel, Profs]: Profiles)
    {
        for (auto *Prof: Profs) { delete Prof; }
    }
    for (auto &[Key, Prof]: Profile2Ds)
    {
        delete Prof;
    }
    for (auto &[Channel, Hist]: RisePowerHists)
    {
        delete Hist;
    }
}
