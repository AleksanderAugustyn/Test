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

    const std::map<Int_t, std::pair<std::string, std::string> > ChannelMap = {
        {4, {"xa", "X Anode A Signal"}},
        {7, {"xb", "X Anode B Signal"}},
        {6, {"ya", "Y Anode A Signal"}},
        {5, {"yb", "Y Anode B Signal"}}
    };

    TTreeReader Reader;
    Reader.SetTree(TreeInput);
    TTreeReaderValue<Double_t> HighGainPosX(Reader, "high_gain_.pos_x_");
    TTreeReaderValue<Double_t> HighGainPosY(Reader, "high_gain_.pos_y_");
    TTreeReaderArray<processor_struct::ROOTDEV> RootDevVector = {Reader, "rootdev_vec_"};
    Reader.SetEntry(Entry);

    /*// Print position information
    std::cout << "\nEvent " << Entry << " Position Information:" << std::endl;
    std::cout << "High Gain Position: X = " << std::fixed << std::setprecision(5) << *HighGainPosX << ", Y = " << *HighGainPosY << std::endl;*/

    std::map<std::string, TGraph *> TraceGraphs;
    std::vector<TF1 *> FitFunctions;

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
        /*std::map<std::string, std::vector<Double_t> > FitParameters;
        std::cout << "\nFit Parameters:" << std::endl;
        std::cout << "Channel\tAmplitude\tPosition\tDecay(τ1)\tRise(τ2)\tPower\tBaseline" << std::endl;
        std::cout << "------------------------------------------------------------------------" << std::endl;*/

        for (UInt_t DeviceIndex = 0; DeviceIndex < RootDevVector.GetSize(); DeviceIndex++)
        {
            const auto &Device = RootDevVector.At(DeviceIndex);

            if (!Device.hasValidTimingAnalysis || !Device.hasValidWaveformAnalysis)
            {
                continue;
            }

            if (Device.subtype == "dynode_high")
            {
                const auto TraceGraph = CreateTraceGraph(Device, "Dynode High Signal", static_cast<Int_t>(DeviceIndex));
                TraceGraphs["dynode"] = TraceGraph;

                try
                {
                    /*// Print header for dynode fit parameters
                    std::cout << "\nDynode Signal Fit Parameters:" << std::endl;
                    std::cout << "-------------------------------------------------------------" << std::endl;
                    std::cout << "Parameter           Value           Error" << std::endl;
                    std::cout << "-------------------------------------------------------------" << std::endl;*/

                    TF1 *DynodeFitResult = FitDynodePeak(TraceGraph, 0, TraceGraph->GetN());
                    DynodeFitResult->SetLineColor(kRed);
                    DynodeFitResult->SetLineWidth(3);
                    DynodeFitResult->SetNpx(2000);
                    FitFunctions.push_back(DynodeFitResult);

                    /*// Store and print fit parameters with formatting
                    const char *ParamNames[] = {
                        "Amplitude", "Peak Position", "Fast Decay τ1", "Slow Decay τ2",
                        "Rise Time τr", "Undershoot Amp", "Undershoot τu", "Fast Fraction",
                        "Baseline"
                    };

                    for (Int_t i = 0; i < DynodeFitResult->GetNpar(); i++)
                    {
                        const Double_t Value = DynodeFitResult->GetParameter(i);
                        const Double_t Error = DynodeFitResult->GetParError(i);

                        std::cout << std::left << std::setw(18) << ParamNames[i]
                                << std::fixed << std::setprecision(3)
                                << std::setw(15) << Value
                                << " ± " << std::setw(10) << Error << std::endl;
                    }*/

                    /*// Calculate and print fit quality metrics
                    const Double_t ChiSquare = DynodeFitResult->GetChisquare();
                    const Int_t NDF = DynodeFitResult->GetNDF();
                    const Double_t RedChiSquare = ChiSquare / NDF;

                    std::cout << "\nFit Quality Metrics:" << std::endl;
                    std::cout << "Chi-square: " << std::fixed << std::setprecision(2) << ChiSquare << std::endl;
                    std::cout << "NDF: " << NDF << std::endl;
                    std::cout << "Reduced chi-square: " << std::fixed << std::setprecision(3) << RedChiSquare << std::endl;*/

                    // Add statistics box to the plot
                    gStyle->SetOptFit(1); // Show fit statistics on plot
                    gStyle->SetStatX(0.9);
                    gStyle->SetStatY(0.9);
                }
                catch (const std::exception &Error)
                {
                    std::cerr << "Dynode fitting error: " << Error.what() << std::endl;
                }
            }

            if (Device.subtype == "anode_high")
            {
                if (auto ChannelIter = ChannelMap.find(Device.chanNum); ChannelIter != ChannelMap.end())
                {
                    const auto TraceGraph = CreateTraceGraph(Device, ChannelIter->second.second, static_cast<Int_t>(DeviceIndex));
                    TraceGraphs[ChannelIter->second.first] = TraceGraph;

                    try
                    {
                        TF1 *FitResult = FitPeakToTrace(TraceGraph, 0, TraceGraph->GetN());
                        FitResult->SetLineColor(kRed);
                        FitResult->SetLineWidth(5);
                        FitResult->SetNpx(2000);
                        FitFunctions.push_back(FitResult);

                        /*// Store fit parameters
                        std::vector<Double_t> Params;
                        Params.reserve(FitResult->GetNpar());
                        for (Int_t i = 0; i < FitResult->GetNpar(); i++)
                        {
                            Params.push_back(FitResult->GetParameter(i));
                        }
                        FitParameters[ChannelIter->second.first] = Params;

                        // Print fit parameters
                        std::cout << std::fixed << std::setprecision(5)
                                << ChannelIter->second.first << "\t"
                                << FitResult->GetParameter(0) << "\t\t" // Amplitude
                                << FitResult->GetParameter(1) << "\t\t" // Peak Position
                                << FitResult->GetParameter(2) << "\t\t" // Decay Time (τ1)
                                << FitResult->GetParameter(3) << "\t\t" // Rise Time (τ2)
                                << FitResult->GetParameter(4) << "\t" // Rise Time Power
                                << FitResult->GetParameter(5) << std::endl; // Baseline*/
                    }
                    catch (const std::exception &Error)
                    {
                        std::cerr << "Fitting error for " << ChannelIter->second.first
                                << ": " << Error.what() << std::endl;
                    }
                }
            }
        }
    }

    // Create visualization (rest of the function remains the same)
    if (TraceGraphs.size() == 5)
    {
        const auto CombinedCanvas = new TCanvas("AllTraces", "All Traces", 2000, 1600);
        CombinedCanvas->SetWindowSize(2000, 1600);
        CombinedCanvas->Divide(1, 5, 0.001, 0.005);

        // Print the X and Y positions on the canvas, at the top
        char PositionText[100];
        sprintf(PositionText, "X = %.5f, Y = %.5f", *HighGainPosX, *HighGainPosY);
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
