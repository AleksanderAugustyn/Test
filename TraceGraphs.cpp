#include <vector>
#include <map>
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

#include "PaassRootStruct.hpp"

#include "main.h"

void SaveTraceGraphs(TTree* TreeInput, const Long64_t Entry, const char* ImagePath)
{
    if (!MeetsSelectionCriteria(TreeInput, Entry))
    {
        return;
    }

    // Define the channels we want to plot and their titles
    const std::map<Int_t, std::pair<std::string, std::string>> ChannelMap = {
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
    std::map<std::string, TGraph*> TraceGraphs;

    // Set global style parameters
    gStyle->SetTextSize(0.2);    // Increase default text size
    gStyle->SetLabelSize(0.1);   // Increase axis label size
    gStyle->SetTitleSize(2.5);   // Increase title size
    //gStyle-> SetTitleFontSize(18); // Increase title font size

    if (RootDevVector.GetSize() > 0)
    {
        for (UInt_t DeviceIndex = 0; DeviceIndex < RootDevVector.GetSize(); DeviceIndex++)
        {
            const auto& Device = RootDevVector.At(DeviceIndex);

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
    if (TraceGraphs.size() == 5)  // 4 anodes + 1 dynode
    {
        // Create larger canvas with higher DPI
        // Set higher resolution for saved images
        gStyle->SetImageScaling(3.0);  // Increase image resolution

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
        for (auto& [Key, Graph] : TraceGraphs)
        {
            delete Graph;
        }
    }
}

TGraph* CreateTraceGraph(const processor_struct::ROOTDEV& Device,
    const std::string& Title, const Int_t DeviceIndex)
{
    const auto PointCount = static_cast<Int_t>(Device.trace.size());
    const auto TraceGraph = new TGraph(PointCount);

    char GraphName[100];
    sprintf(GraphName, "TraceGraph_%d", DeviceIndex);
    TraceGraph->SetName(GraphName);
    TraceGraph->SetTitle(Title.c_str());

    // Enhanced appearance
    TraceGraph->SetMarkerStyle(8);
    TraceGraph->SetMarkerSize(1.0);    // Increased marker size
    TraceGraph->SetMarkerColor(kBlue);
    TraceGraph->SetLineColor(kBlue);
    TraceGraph->SetLineWidth(2);        // Thicker lines

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
void SaveTraceGraphsWithFit(TTree* TreeInput, const Long64_t Entry, const char* ImagePath)
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
    TTreeReaderArray<processor_struct::ROOTDEV> RootDevVector = {Reader, "rootdev_vec_"};
    Reader.SetEntry(Entry);

    std::map<std::string, TGraph*> TraceGraphs;
    std::vector<TF1*> FitFunctions;

    // Global style settings
    gStyle->SetOptTitle(1);        // Show titles
    gStyle->SetTitleSize(0.08, "t");  // Title size for plot titles
    gStyle->SetTitleSize(0.06, "xy"); // Title size for axes
    gStyle->SetLabelSize(0.05, "xy"); // Label size for axes
    gStyle->SetTitleOffset(0.8, "y"); // Adjust y-title position
    gStyle->SetTitleOffset(0.9, "x"); // Adjust x-title position
    gStyle->SetTitleFontSize(0.08);   // Global title font size
    gStyle->SetGridWidth(2);          // Make grid lines thicker
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
                const auto TraceGraph = CreateTraceGraph(Device, "Dynode High Signal", static_cast<Int_t>(DeviceIndex));
                TraceGraphs["dynode"] = TraceGraph;
                continue;
            }

            if (Device.subtype == "anode_high")
            {
                if (auto ChannelIter = ChannelMap.find(Device.chanNum); ChannelIter != ChannelMap.end())
                {
                    const auto TraceGraph = CreateTraceGraph(Device, ChannelIter->second.second, static_cast<Int_t>(DeviceIndex));
                    TraceGraphs[ChannelIter->second.first] = TraceGraph;

                    try
                    {
                        TF1* FitResult = FitPeakToTrace(TraceGraph, 0, TraceGraph->GetN());
                        // Set fit function attributes - much thicker line
                        FitResult->SetLineColor(kRed);
                        FitResult->SetLineWidth(5);    // Significantly increased line width
                        FitResult->SetNpx(2000);       // More points for smoother curve
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

    if (TraceGraphs.size() == 5)  // 4 anodes + 1 dynode
    {
        // Create canvas with higher pixel density
        const auto CombinedCanvas = new TCanvas("AllTraces", "All Traces", 2000, 1600);
        CombinedCanvas->SetWindowSize(2000, 1600);
        CombinedCanvas->Divide(1, 5, 0.001, 0.005);  // Minimal spacing between pads

        const std::vector<std::string> PlotOrder = {"xa", "xb", "ya", "yb", "dynode"};

        for (size_t i = 0; i < PlotOrder.size(); i++)
        {
            CombinedCanvas->cd(static_cast<Int_t>(i + 1));

            // Set pad properties
            gPad->SetGrid(1, 1);  // Enable both x and y grid
            gPad->SetLeftMargin(0.05);
            gPad->SetRightMargin(0.05);
            gPad->SetBottomMargin(0.20);  // Increased for bigger labels
            gPad->SetTopMargin(0.15);     // Increased for bigger title

            const auto graph = TraceGraphs[PlotOrder[i]];
            graph->Draw("ALP");

            // Customize graph appearance
            graph->GetXaxis()->SetTitle("Time [ns]");
            graph->GetYaxis()->SetTitle("Amplitude");
            graph->GetXaxis()->CenterTitle();
            graph->GetYaxis()->CenterTitle();

            // Force update for proper layering
            gPad->Update();

            if (PlotOrder[i] != "dynode")
            {
                for (const auto FitFunc : FitFunctions)
                {
                    if (TString(FitFunc->GetName()).Contains(PlotOrder[i]))
                    {
                        // Redraw with same attributes to ensure visibility
                        FitFunc->SetLineColor(kRed);
                        FitFunc->SetLineWidth(5);
                        FitFunc->SetNpx(2000);
                        FitFunc->Draw("same C");

                        // Force immediate update
                        gPad->Modified();
                        gPad->Update();
                    }
                }
            }
        }

        // Save with high quality
        char CombinedPngName[100];
        sprintf(CombinedPngName, "%s/event_%lld_traces_fitted.png", ImagePath, Entry);
        CombinedCanvas->SaveAs(CombinedPngName);

        delete CombinedCanvas;
        for (auto& [Key, Graph] : TraceGraphs)
        {
            delete Graph;
        }
        for (const auto& FitFunc : FitFunctions)
        {
            delete FitFunc;
        }
    }
}