#include <vector>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TLegend.h>
#include <TAxis.h>

#include "PaassRootStruct.hpp"

// Function to check if an event meets our selection criteria
Bool_t MeetsSelectionCriteria(TTree* TreeInput, Long64_t Entry)
{
    // Setup branch variables
    Int_t HighGainValid;
    Int_t LowGainValid;
    Double_t HighGainQdc;

    // Set branch addresses
    TreeInput->SetBranchAddress("high_gain_.valid_", &HighGainValid);
    TreeInput->SetBranchAddress("low_gain_.valid_", &LowGainValid);
    TreeInput->SetBranchAddress("high_gain_.qdc_", &HighGainQdc);

    // Get the values for this entry
    TreeInput->GetEntry(Entry);

    // Check conditions
    if (HighGainValid != 1)
    {
        return false;
    }

    if (LowGainValid != 0)
    {
        return false;
    }

    if (HighGainQdc <= 10000 || HighGainQdc >= 50000)
    {
        return false;
    }

    return true;
}

void SaveTraceGraphs(TTree* TreeInput, Long64_t Entry, const char* ImagePath = "./")
{
    // First check if this event meets our criteria
    if (!MeetsSelectionCriteria(TreeInput, Entry))
    {
        return;
    }

    TTreeReader Reader;
    Reader.SetTree(TreeInput);
    TTreeReaderArray<processor_struct::ROOTDEV> RootDevVector = {Reader, "rootdev_vec_"};
    Reader.SetEntry(Entry);

    Int_t TraceCount = 0;
    std::vector<TGraph*> TraceGraphVector;

    if (RootDevVector.GetSize() > 0)
    {
        std::cout << "Found trace data: " << RootDevVector.GetSize() << std::endl;

        for (UInt_t DeviceIndex = 0; DeviceIndex < RootDevVector.GetSize(); DeviceIndex++)
        {
            if (RootDevVector.At(DeviceIndex).hasValidTimingAnalysis)
            {
                Int_t PointCount = RootDevVector.At(DeviceIndex).trace.size();
                std::cout << "Number of Points: " << PointCount << std::endl;

                TGraph* TraceGraph = new TGraph(PointCount);
                char GraphName[50];
                sprintf(GraphName, "TraceGraph_%d", DeviceIndex);
                TraceGraph->SetName(GraphName);
                TraceGraph->SetTitle(GraphName);

                // Enhance graph appearance
                TraceGraph->SetMarkerStyle(8);
                TraceGraph->SetMarkerSize(0.5);
                TraceGraph->SetMarkerColor(kBlue);
                TraceGraph->SetLineColor(kBlue);
                TraceGraph->SetLineWidth(1);

                // Set axis titles
                TraceGraph->GetXaxis()->SetTitle("Sample Number");
                TraceGraph->GetYaxis()->SetTitle("Amplitude");

                for (Int_t Point = 0; Point < PointCount; Point++)
                {
                    TraceGraph->SetPoint(Point, Point,
                        RootDevVector.At(DeviceIndex).trace.at(Point));
                }

                TraceGraphVector.push_back(TraceGraph);
                TraceCount++;
            }
        }
    }

    std::cout << "Valid trace data: " << TraceCount << std::endl;

    if (TraceCount > 0)
    {
        // Create combined canvas
        TCanvas* CombinedCanvas = new TCanvas("AllTraces", "All Traces", 800,
            200 * TraceCount);
        CombinedCanvas->Divide(1, TraceCount);

        for (Int_t GraphIndex = 0; GraphIndex < TraceCount; GraphIndex++)
        {
            CombinedCanvas->cd(GraphIndex + 1);
            gPad->SetGrid();
            TraceGraphVector[GraphIndex]->Draw("ALP");
        }

        // Save combined view with event number in filename
        char CombinedPngName[100];
        sprintf(CombinedPngName, "%s/event_%lld_traces.png", ImagePath, Entry);
        CombinedCanvas->SaveAs(CombinedPngName);

        // Cleanup
        delete CombinedCanvas;

        for (auto Graph : TraceGraphVector)
        {
            delete Graph;
        }
    }
}

int main()
{
    // Load required libraries
    gSystem->Load("libPaassRootStruct.so");
    gSystem->Load("libyaml-cpp.so");
    gSystem->Load("libTraceAnalyzerLib.so");
    gSystem->Load("libmerger_data_dic.so");
    gSystem->Load("libMergerLib.so");

    TFile* InputFile = TFile::Open("pixie_bigrips_traces_119_35.root");
    TTree* Tree = (TTree*)InputFile->Get("pspmt");

    Long64_t Entries = Tree->GetEntries();
    std::cout << "Total number of entries: " << Entries << std::endl;

    // Create output directory if it doesn't exist
    gSystem->Exec("mkdir -p traces");

    // Set maximum number of events to save
    const Long64_t MaxEventsToSave = 100;
    Long64_t SavedEvents = 0;
    std::vector<Long64_t> SelectedEventNumbers;

    // First pass: scan all events to count total qualifying events
    std::cout << "Scanning all events to count those meeting criteria..." << std::endl;
    Long64_t TotalQualifyingEvents = 0;

    for (Long64_t Event = 0; Event < Entries; Event++)
    {
        if (Event % 10000 == 0)
        {
            std::cout << "Scanning event " << Event << "/" << Entries << "\r" << std::flush;
        }

        if (MeetsSelectionCriteria(Tree, Event))
        {
            TotalQualifyingEvents++;
            if (SavedEvents < MaxEventsToSave)
            {
                SelectedEventNumbers.push_back(Event);
                SavedEvents++;
            }
        }
    }
    std::cout << std::endl;  // New line after progress indicator

    std::cout << "Total events meeting criteria: " << TotalQualifyingEvents << std::endl;
    std::cout << "Selected first " << SavedEvents << " events for plotting" << std::endl;

    // Second pass: save traces for selected events
    std::cout << "Saving traces for selected events..." << std::endl;
    for (const auto& EventNumber : SelectedEventNumbers)
    {
        std::cout << "Processing event " << EventNumber << std::endl;
        SaveTraceGraphs(Tree, EventNumber, "traces");
    }

    std::cout << "Finished saving " << SavedEvents << " events" << std::endl;

    InputFile->Close();
    delete InputFile;

    return 0;
}