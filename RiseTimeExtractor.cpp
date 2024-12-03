#include <TFile.h>
#include <TProfile2D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <memory>

class RiseTimeMapExtractor {
private:
    std::vector<std::string> Channels = {"xa", "xb", "ya", "yb"};
    const Int_t NBinsTarget = 50;
    std::map<std::string, TH2D*> RiseTimeMaps;

    void ProcessChannels(TFile* InputFile, TFile* OutputFile) {
        for (const auto& Channel : Channels) {
            std::cout << "Processing channel: " << Channel << std::endl;

            const std::string ProfileName = Channel + "_rise_2d_prof";
            auto* OriginalProfile = dynamic_cast<TProfile2D*>(InputFile->Get(ProfileName.c_str()));
            if (!OriginalProfile) {
                std::cerr << "Failed to get profile for channel " << Channel << std::endl;
                continue;
            }

            auto* RebinnedMap = CreateRebinnedMap(OriginalProfile, Channel);
            if (!RebinnedMap) {
                std::cerr << "Failed to create rebinned map for channel " << Channel << std::endl;
                continue;
            }

            OutputFile->cd();
            RebinnedMap->Write();
            RiseTimeMaps[Channel] = RebinnedMap;

            std::cout << "Completed processing for channel " << Channel << std::endl;
            PrintMapStatistics(RebinnedMap, Channel);
        }
    }

    TH2D* CreateRebinnedMap(TProfile2D* OriginalProfile, const std::string& Channel) {
        const Int_t OriginalBinsX = OriginalProfile->GetNbinsX();
        const Int_t OriginalBinsY = OriginalProfile->GetNbinsY();

        const Int_t RebinFactorX = OriginalBinsX / NBinsTarget;
        const Int_t RebinFactorY = OriginalBinsY / NBinsTarget;

        auto* RebinnedMap = new TH2D(
            (Channel + "_rise_time_map").c_str(),
            (Channel + " Rise Time Map;X Position;Y Position;Rise Time [ns]").c_str(),
            NBinsTarget, 0.0, 0.50,
            NBinsTarget, 0.0, 0.50
        );

        for (Int_t i = 1; i <= NBinsTarget; ++i) {
            for (Int_t j = 1; j <= NBinsTarget; ++j) {
                Double_t SumContent = 0.0;
                Double_t SumWeight = 0.0;

                for (Int_t ix = (i-1)*RebinFactorX + 1; ix <= i*RebinFactorX; ++ix) {
                    for (Int_t iy = (j-1)*RebinFactorY + 1; iy <= j*RebinFactorY; ++iy) {
                        const Double_t BinContent = OriginalProfile->GetBinContent(ix, iy);
                        const Double_t BinEntries = OriginalProfile->GetBinEntries(
                            OriginalProfile->GetBin(ix, iy));

                        if (BinEntries > 0) {
                            SumContent += BinContent * BinEntries;
                            SumWeight += BinEntries;
                        }
                    }
                }

                if (SumWeight > 0) {
                    RebinnedMap->SetBinContent(i, j, SumContent / SumWeight);
                }
            }
        }

        return RebinnedMap;
    }

    void PrintMapStatistics(TH2D* Map, const std::string& Channel) {
        std::cout << "\nStatistics for " << Channel << " rise time map:" << std::endl;
        std::cout << "Mean rise time: " << Map->GetMean(3) << " ns" << std::endl;
        std::cout << "RMS rise time: " << Map->GetRMS(3) << " ns" << std::endl;
        std::cout << "Min rise time: " << Map->GetMinimum() << " ns" << std::endl;
        std::cout << "Max rise time: " << Map->GetMaximum() << " ns" << std::endl;
    }

    void CreateVisualization() {
        gStyle->SetOptStat(1);
        gStyle->SetPalette(kBird);

        auto Canvas = std::make_unique<TCanvas>("RiseTimeMaps", "Rise Time Maps", 2000, 1000);
        Canvas->Divide(2, 2, 0.01, 0.01);

        Int_t Pad = 1;
        for (const auto& Channel : Channels) {
            Canvas->cd(Pad++);

            auto* Map = RiseTimeMaps[Channel];
            if (!Map) continue;

            Map->GetXaxis()->SetTitleSize(0.05);
            Map->GetYaxis()->SetTitleSize(0.05);
            Map->GetZaxis()->SetTitleSize(0.05);
            Map->GetXaxis()->SetTitleOffset(1.2);
            Map->GetYaxis()->SetTitleOffset(1.2);
            Map->GetZaxis()->SetTitleOffset(1.2);

            Map->GetXaxis()->SetLabelSize(0.04);
            Map->GetYaxis()->SetLabelSize(0.04);
            Map->GetZaxis()->SetLabelSize(0.04);

            gPad->SetRightMargin(0.15);
            gPad->SetLeftMargin(0.15);
            gPad->SetTopMargin(0.1);
            gPad->SetBottomMargin(0.15);

            Map->Draw("COLZ");
        }

        Canvas->SaveAs("rise_time_maps.png");
    }

public:
    ~RiseTimeMapExtractor() {
        // ROOT will handle cleanup of histograms written to file
        RiseTimeMaps.clear();
    }

    void Process(const char* InputFileName) {
        std::unique_ptr<TFile> InputFile(TFile::Open(InputFileName, "READ"));
        if (!InputFile || InputFile->IsZombie()) {
            throw std::runtime_error("Failed to open input file");
        }

        std::unique_ptr<TFile> OutputFile(TFile::Open("rise_time_maps.root", "RECREATE"));
        if (!OutputFile || OutputFile->IsZombie()) {
            throw std::runtime_error("Failed to create output file");
        }

        ProcessChannels(InputFile.get(), OutputFile.get());
        OutputFile->Write();

        CreateVisualization();
    }
};