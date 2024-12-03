#include <iomanip>
#include <vector>
#include <map>
#include <iostream>

#include <TROOT.h>
#include <TFile.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TAxis.h>
#include <TF1.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TText.h>

#include "PaassRootStruct.hpp"

#include "main.h"

void SaveTraceGraphs(TTree *TreeInput, const Long64_t Entry, const char *ImagePath)
{
    if (!MeetsSelectionCriteria(TreeInput, Entry))
    {
        return;
    }

    // Define the channels we want to plot and their titles
    const std::map<Int_t, std::pair<std::string, std::string> > ChannelMap = {
        {4, {"xa", "X Anode A Signal"}},
        {7, {"xb", "X Anode B Signal"}},
        {6, {"ya", "Y Anode A Signal"}},
        {5, {"yb", "Y Anode B Signal"}}
    };

    TTreeReader Reader;
    Reader.SetTree(TreeInput);
    TTreeReaderArray<processor_struct::ROOTDEV> RootDevVector = {Reader, "rootdev_vec_"};
    Reader.SetEntry(Entry);

    // Store graphs for each channel type
    std::map<std::string, TGraph *> TraceGraphs;

    // Set global style parameters
    gStyle->SetTextSize(0.2); // Increase default text size
    gStyle->SetLabelSize(0.1); // Increase axis label size
    gStyle->SetTitleSize(2.5); // Increase title size
    //gStyle-> SetTitleFontSize(18); // Increase title font size

    if (RootDevVector.GetSize() > 0)
    {
        for (UInt_t DeviceIndex = 0; DeviceIndex < RootDevVector.GetSize(); DeviceIndex++)
        {
            const auto &Device = RootDevVector.At(DeviceIndex);

            if (!Device.hasValidTimingAnalysis || !Device.hasValidWaveformAnalysis)
            {
                continue;
            }

            // Process dynode_high
            if (Device.subtype == "dynode_high")
            {
                const auto TraceGraph = CreateTraceGraph(Device, "Dynode High Signal", static_cast<Int_t>(DeviceIndex));
                TraceGraphs["dynode"] = TraceGraph;
                continue;
            }

            // Process anodes
            if (Device.subtype == "anode_high")
            {
                if (auto ChannelIter = ChannelMap.find(Device.chanNum); ChannelIter != ChannelMap.end())
                {
                    const auto TraceGraph = CreateTraceGraph(Device,
                                                             ChannelIter->second.second, static_cast<Int_t>(DeviceIndex));
                    TraceGraphs[ChannelIter->second.first] = TraceGraph;
                }
            }
        }
    }

    // Only proceed if we have all required traces
    if (TraceGraphs.size() == 5) // 4 anodes + 1 dynode
    {
        // Create larger canvas with higher DPI
        // Set higher resolution for saved images
        gStyle->SetImageScaling(3.0); // Increase image resolution

        const auto CombinedCanvas = new TCanvas("AllTraces", "All Traces", 1600, 1000);
        CombinedCanvas->Divide(1, 5);

        // Plot traces in specific order
        const std::vector<std::string> PlotOrder = {"xa", "xb", "ya", "yb", "dynode"};

        for (size_t i = 0; i < PlotOrder.size(); i++)
        {
            CombinedCanvas->cd(static_cast<Int_t>(i + 1));
            gPad->SetGrid();
            TraceGraphs[PlotOrder[i]]->Draw("ALP");
        }

        // Save combined view
        char CombinedPngName[100];
        sprintf(CombinedPngName, "%s/event_%lld_traces.png", ImagePath, Entry);
        CombinedCanvas->SaveAs(CombinedPngName);

        // Cleanup
        delete CombinedCanvas;
        for (auto &[Key, Graph]: TraceGraphs)
        {
            delete Graph;
        }
    }
}

TGraph *CreateTraceGraph(const processor_struct::ROOTDEV &Device,
                         const std::string &Title, const Int_t DeviceIndex)
{
    const auto PointCount = static_cast<Int_t>(Device.trace.size());
    const auto TraceGraph = new TGraph(PointCount);

    char GraphName[100];
    sprintf(GraphName, "TraceGraph_%d", DeviceIndex);
    TraceGraph->SetName(GraphName);
    TraceGraph->SetTitle(Title.c_str());

    // Enhanced appearance
    TraceGraph->SetMarkerStyle(8);
    TraceGraph->SetMarkerSize(1.0); // Increased marker size
    TraceGraph->SetMarkerColor(kBlue);
    TraceGraph->SetLineColor(kBlue);
    TraceGraph->SetLineWidth(2); // Thicker lines

    // Set axis titles with larger font
    TraceGraph->GetXaxis()->SetTitle("Time [ns]");
    TraceGraph->GetYaxis()->SetTitle("Amplitude");
    TraceGraph->GetXaxis()->SetTitleSize(0.1);
    TraceGraph->GetYaxis()->SetTitleSize(0.1);
    TraceGraph->GetXaxis()->SetLabelSize(0.1);
    TraceGraph->GetYaxis()->SetLabelSize(0.1);

    for (Int_t Point = 0; Point < PointCount; Point++)
    {
        TraceGraph->SetPoint(Point, Point, Device.trace.at(Point));
    }

    return TraceGraph;
}

/**
 * Updates the SaveTraceGraphs function to include peak fitting
 * @param TreeInput Pointer to the input tree
 * @param Entry Entry number to process
 * @param ImagePath Path to save the output images
 */
void SaveTraceGraphsWithFit(TTree *TreeInput, const Long64_t Entry, const char *ImagePath)
{
    if (!MeetsSelectionCriteria(TreeInput, Entry))
    {
        return;
    }

    const std::map<Int_t, std::pair<std::string, std::string>> ChannelMap = {
        {4, {"xa", "X Anode A Signal"}},
        {7, {"xb", "X Anode B Signal"}},
        {6, {"ya", "Y Anode A Signal"}},
        {5, {"yb", "Y Anode B Signal"}}
    };

    TTreeReader Reader;
    Reader.SetTree(TreeInput);
    TTreeReaderValue<Double_t> HighGainPosX(Reader, "high_gain_.pos_x_");
    TTreeReaderValue<Double_t> HighGainPosY(Reader, "high_gain_.pos_y_");
    TTreeReaderArray<processor_struct::ROOTDEV> RootDevVector(Reader, "rootdev_vec_");
    Reader.SetEntry(Entry);

    std::map<std::string, TGraph*> TraceGraphs;
    std::vector<TF1*> FitFunctions;

    Double_t PositionX = *HighGainPosX;
    [[maybe_unused]] Double_t PositionY = *HighGainPosY;

    std::cout << "Position: " << PositionX << ", " << PositionY << std::endl;

    // Global style settings
    gStyle->SetOptTitle(1);
    gStyle->SetTitleSize(0.08, "t");
    gStyle->SetTitleSize(0.06, "xy");
    gStyle->SetLabelSize(0.05, "xy");
    gStyle->SetTitleOffset(0.8, "y");
    gStyle->SetTitleOffset(0.9, "x");
    gStyle->SetTitleFontSize(0.08);
    gStyle->SetGridWidth(2);
    gStyle->SetLineWidth(1);
    gStyle->SetOptFit(1);
    gStyle->SetFuncWidth(4);

    if (RootDevVector.GetSize() > 0)
    {
        for (UInt_t DeviceIndex = 0; DeviceIndex < RootDevVector.GetSize(); DeviceIndex++)
        {
            const auto& Device = RootDevVector.At(DeviceIndex);

            if (!Device.hasValidTimingAnalysis || !Device.hasValidWaveformAnalysis)
            {
                continue;
            }

            if (Device.subtype == "dynode_high")
            {
                const auto TraceGraph = CreateTraceGraph(Device, "Dynode High Signal",
                                                       static_cast<Int_t>(DeviceIndex));
                TraceGraphs["dynode"] = TraceGraph;

                try
                {
                    TF1* DynodeFitResult = FitDynodePeak(TraceGraph, 0, TraceGraph->GetN());
                    DynodeFitResult->SetLineColor(kRed);
                    DynodeFitResult->SetLineWidth(3);
                    DynodeFitResult->SetNpx(2000);
                    FitFunctions.push_back(DynodeFitResult);
                }
                catch (const std::exception& Error)
                {
                    std::cerr << "Dynode fitting error: " << Error.what() << std::endl;
                }
            }
            else if (Device.subtype == "anode_high")
            {
                if (auto ChannelIter = ChannelMap.find(Device.chanNum);
                    ChannelIter != ChannelMap.end())
                {
                    const auto TraceGraph = CreateTraceGraph(Device,
                                                           ChannelIter->second.second,
                                                           static_cast<Int_t>(DeviceIndex));
                    TraceGraphs[ChannelIter->second.first] = TraceGraph;

                    try
                    {
                        // Always use X position for rise power calculation
                        TF1* FitResult = FitPeakToTrace(TraceGraph, 0, TraceGraph->GetN(),
                                                      ChannelIter->second.first, PositionX, PositionY);
                        FitResult->SetLineColor(kRed);
                        FitResult->SetLineWidth(5);
                        FitResult->SetNpx(2000);
                        FitFunctions.push_back(FitResult);
                    }
                    catch (const std::exception& Error)
                    {
                        std::cerr << "Fitting error for " << ChannelIter->second.first
                                 << ": " << Error.what() << std::endl;
                    }
                }
            }
        }
    }

    // Create visualization
    if (TraceGraphs.size() == 5)
    {
        const auto CombinedCanvas = new TCanvas("AllTraces", "All Traces", 2000, 1600);
        CombinedCanvas->SetWindowSize(2000, 1600);
        CombinedCanvas->Divide(1, 5, 0.001, 0.005);

        // Print the X and Y positions on the canvas, at the top
        char PositionText[100];
        sprintf(PositionText, "X = %.5f, Y = %.5f", PositionX, PositionY);
        auto PositionLabel = new TText(0.5, 0.95, PositionText);
        PositionLabel->SetNDC();
        PositionLabel->SetTextSize(0.05);
        PositionLabel->Draw();

        const std::vector<std::string> PlotOrder = {"xa", "xb", "ya", "yb", "dynode"};

        for (size_t i = 0; i < PlotOrder.size(); i++)
        {
            CombinedCanvas->cd(static_cast<Int_t>(i + 1));
            gPad->SetGrid(1, 1);
            gPad->SetLeftMargin(0.05);
            gPad->SetRightMargin(0.05);
            gPad->SetBottomMargin(0.20);
            gPad->SetTopMargin(0.15);

            const auto graph = TraceGraphs[PlotOrder[i]];
            graph->Draw("ALP");

            graph->GetXaxis()->SetTitle("Time [ns]");
            graph->GetYaxis()->SetTitle("Amplitude");
            graph->GetXaxis()->CenterTitle();
            graph->GetYaxis()->CenterTitle();

            gPad->Update();

            if (PlotOrder[i] != "dynode")
            {
                for (const auto FitFunc: FitFunctions)
                {
                    if (TString(FitFunc->GetName()).Contains(PlotOrder[i]))
                    {
                        FitFunc->Draw("same C");
                        gPad->Modified();
                        gPad->Update();
                    }
                }
            }
        }

        char CombinedPngName[100];
        sprintf(CombinedPngName, "%s/Event_%lld_Traces_Fit.png", ImagePath, Entry);
        CombinedCanvas->SaveAs(CombinedPngName);

        delete CombinedCanvas;
        for (auto &[Key, Graph]: TraceGraphs)
        {
            delete Graph;
        }
        for (const auto &FitFunc: FitFunctions)
        {
            delete FitFunc;
        }
    }
}

void GraphFirstNEvents(TTree *TreeInput, const std::vector<Long64_t> &QualifyingEvents,
                       const Long64_t NumberOfEvents, const char *OutputPath)
{
    if (!TreeInput)
    {
        throw std::runtime_error("Invalid tree pointer");
    }

    if (NumberOfEvents <= 0)
    {
        throw std::runtime_error("Number of events must be positive");
    }

    if (QualifyingEvents.empty())
    {
        throw std::runtime_error("No qualifying events found");
    }

    const Long64_t EventsToProcess = std::min(static_cast<Long64_t>(QualifyingEvents.size()),
                                              NumberOfEvents);

    // std::cout << "Processing first " << EventsToProcess << " qualifying events..." << std::endl;

    for (Long64_t i = 0; i < EventsToProcess; i++)
    {
        std::cout << "Processing event " << QualifyingEvents[i] << " ("
                << i + 1 << "/" << EventsToProcess << ")" << std::endl;
        SaveTraceGraphsWithFit(TreeInput, QualifyingEvents[i], OutputPath);
    }
}

/**
 * Extracts fit parameters from an anode fit function
 * @param FitFunc Pointer to the fitted function
 * @return Optional channel fit parameters
 */
std::optional<AnalysisResults::ChannelFit> ExtractAnodeFitParameters(const TF1* FitFunc)
{
    if (!FitFunc)
    {
        return std::nullopt;
    }

    AnalysisResults::ChannelFit Params;
    Params.Amplitude = FitFunc->GetParameter(0);
    Params.PeakPosition = FitFunc->GetParameter(1);
    Params.DecayConstant = FitFunc->GetParameter(2);
    Params.RiseTimeConstant = FitFunc->GetParameter(3);
    Params.RisePower = FitFunc->GetParameter(4);
    Params.Baseline = FitFunc->GetParameter(5);

    return Params;
}

/**
 * Extracts fit parameters from a dynode fit function
 * @param FitFunc Pointer to the fitted function
 * @return Optional dynode fit parameters
 */
std::optional<AnalysisResults::DynodeFit> ExtractDynodeFitParameters(const TF1* FitFunc)
{
    if (!FitFunc)
    {
        return std::nullopt;
    }

    AnalysisResults::DynodeFit Params;
    Params.Amplitude = FitFunc->GetParameter(0);
    Params.PeakPosition = FitFunc->GetParameter(1);
    Params.FastDecay = FitFunc->GetParameter(2);
    Params.SlowDecay = FitFunc->GetParameter(3);
    Params.RiseTime = FitFunc->GetParameter(4);
    Params.UndershootAmp = FitFunc->GetParameter(5);
    Params.UndershootRecovery = FitFunc->GetParameter(6);
    Params.FastFraction = FitFunc->GetParameter(7);
    Params.Baseline = FitFunc->GetParameter(8);

    return Params;
}

/**
 * Modified version of SaveTraceGraphsWithFit that returns analysis results
 * @param TreeInput Pointer to the input tree
 * @param Entry Entry number to process
 * @return Optional analysis results containing fit parameters and positions
 */
std::optional<AnalysisResults> GetEventFitParameters(TTree* TreeInput, const Long64_t Entry)
{
    if (!MeetsSelectionCriteria(TreeInput, Entry))
    {
        return std::nullopt;
    }

    const std::map<Int_t, std::pair<std::string, std::string>> ChannelMap = {
        {4, {"xa", "X Anode A Signal"}},
        {7, {"xb", "X Anode B Signal"}},
        {6, {"ya", "Y Anode A Signal"}},
        {5, {"yb", "Y Anode B Signal"}}
    };

    AnalysisResults Results;
    Results.EventNumber = Entry;

    TTreeReader Reader;
    Reader.SetTree(TreeInput);
    TTreeReaderValue<Double_t> HighGainPosX(Reader, "high_gain_.pos_x_");
    TTreeReaderValue<Double_t> HighGainPosY(Reader, "high_gain_.pos_y_");
    TTreeReaderArray<processor_struct::ROOTDEV> RootDevVector(Reader, "rootdev_vec_");
    Reader.SetEntry(Entry);

    Results.PosX = *HighGainPosX;
    Results.PosY = *HighGainPosY;

    Bool_t ValidFits = true;

    if (RootDevVector.GetSize() > 0)
    {
        for (UInt_t DeviceIndex = 0; DeviceIndex < RootDevVector.GetSize(); DeviceIndex++)
        {
            const auto& Device = RootDevVector.At(DeviceIndex);

            if (!Device.hasValidTimingAnalysis || !Device.hasValidWaveformAnalysis)
            {
                continue;
            }

            auto* TraceGraph = new TGraph(static_cast<Int_t>(Device.trace.size()));
            for (size_t i = 0; i < Device.trace.size(); i++)
            {
                TraceGraph->SetPoint(static_cast<Int_t>(i),
                                   static_cast<Double_t>(i), Device.trace[i]);
            }

            if (Device.subtype == "dynode_high")
            {
                try
                {
                    TF1* DynodeFitResult = FitDynodePeak(TraceGraph, 0, TraceGraph->GetN());
                    auto DynodeParams = ExtractDynodeFitParameters(DynodeFitResult);
                    if (DynodeParams)
                    {
                        Results.DynodeFitParams = *DynodeParams;
                    }
                    else
                    {
                        ValidFits = false;
                    }
                    delete DynodeFitResult;
                }
                catch (const std::exception& Error)
                {
                    std::cerr << "Dynode fitting error for event " << Entry
                             << ": " << Error.what() << std::endl;
                    ValidFits = false;
                }
            }
            else if (Device.subtype == "anode_high")
            {
                if (auto ChannelIter = ChannelMap.find(Device.chanNum);
                    ChannelIter != ChannelMap.end())
                {
                    try
                    {
                        //std::cout << "Position: " << Results.PosX << std::endl;
                        // Always use X position for rise power calculation
                        TF1* FitResult = FitPeakToTrace(TraceGraph, 0.0, TraceGraph->GetN(),
                                                      ChannelIter->second.first, Results.PosX, Results.PosY);
                        auto AnodeParams = ExtractAnodeFitParameters(FitResult);
                        if (AnodeParams)
                        {
                            Results.AnodeFits[ChannelIter->second.first] = *AnodeParams;
                        }
                        else
                        {
                            ValidFits = false;
                        }
                        delete FitResult;
                    }
                    catch (const std::exception& Error)
                    {
                        throw std::runtime_error(Error.what());
                    }
                }
            }
            delete TraceGraph;
        }
    }

    return ValidFits ? std::optional(Results) : std::nullopt;
}

void SaveAnalysisResults(const std::vector<AnalysisResults>& Results, const Int_t RunNumber, const Int_t SubRunNumber)
{
    std::cout << "\n[SaveAnalysisResults] Run " << std::setfill('0') << std::setw(3) << RunNumber
              << "_" << std::setfill('0') << std::setw(2) << SubRunNumber
              << ": Saving " << Results.size() << " events\n" << std::endl;

    // Create output filename based on run numbers
    std::ostringstream OutputFileName;
    OutputFileName << "analysis_"
                  << std::setfill('0') << std::setw(3) << RunNumber
                  << "_"
                  << std::setfill('0') << std::setw(2) << SubRunNumber
                  << ".root";

    TFile OutputFile(OutputFileName.str().c_str(), "RECREATE");

    // Create a tree to store the results
    TTree ResultTree("analysis", "Analysis Results");

    // Variables for tree branches
    AnalysisResults CurrentEvent;

    // Set up branches
    ResultTree.Branch("event_number", &CurrentEvent.EventNumber);
    ResultTree.Branch("pos_x", &CurrentEvent.PosX);
    ResultTree.Branch("pos_y", &CurrentEvent.PosY);

    // Branches for anode fit parameters
    for (const auto& Channel : {"xa", "xb", "ya", "yb"})
    {
        std::string Prefix = std::string(Channel) + "_";
        ResultTree.Branch((Prefix + "amplitude").c_str(), &CurrentEvent.AnodeFits[Channel].Amplitude);
        ResultTree.Branch((Prefix + "peak_position").c_str(), &CurrentEvent.AnodeFits[Channel].PeakPosition);
        ResultTree.Branch((Prefix + "decay_constant").c_str(), &CurrentEvent.AnodeFits[Channel].DecayConstant);
        ResultTree.Branch((Prefix + "rise_time").c_str(), &CurrentEvent.AnodeFits[Channel].RiseTimeConstant);
        ResultTree.Branch((Prefix + "rise_power").c_str(), &CurrentEvent.AnodeFits[Channel].RisePower);
        ResultTree.Branch((Prefix + "baseline").c_str(), &CurrentEvent.AnodeFits[Channel].Baseline);
    }

    // Branches for dynode fit parameters
    ResultTree.Branch("dynode_amplitude", &CurrentEvent.DynodeFitParams.Amplitude);
    ResultTree.Branch("dynode_peak_position", &CurrentEvent.DynodeFitParams.PeakPosition);
    ResultTree.Branch("dynode_fast_decay", &CurrentEvent.DynodeFitParams.FastDecay);
    ResultTree.Branch("dynode_slow_decay", &CurrentEvent.DynodeFitParams.SlowDecay);
    ResultTree.Branch("dynode_rise_time", &CurrentEvent.DynodeFitParams.RiseTime);
    ResultTree.Branch("dynode_undershoot_amp", &CurrentEvent.DynodeFitParams.UndershootAmp);
    ResultTree.Branch("dynode_undershoot_recovery", &CurrentEvent.DynodeFitParams.UndershootRecovery);
    ResultTree.Branch("dynode_fast_fraction", &CurrentEvent.DynodeFitParams.FastFraction);
    ResultTree.Branch("dynode_baseline", &CurrentEvent.DynodeFitParams.Baseline);

    // Fill tree with results
    for (const auto& Result : Results)
    {
        CurrentEvent = Result;
        ResultTree.Fill();
    }

    ResultTree.Write();
    OutputFile.Close();
}