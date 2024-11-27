#include <vector>
#include <iostream>
#include <iomanip>
#include <sstream>

#include <TFile.h>
#include <TTreeReaderValue.h>

#include "main.h"


int main()
{
    try
    {
        LoadRequiredLibraries();

        // Analyze only past runs
        if (0)
        {
            // Perform position analysis on all processed runs
            std::vector<std::pair<Int_t, Int_t> > RunsToAnalyze = {
                {55, 20}, {55, 21},
                {55, 22}, {55, 23},
                {55, 24}, {55, 25},
                {55, 26}, {55, 27},
                {55, 28}, {55, 29},
                {54, 20}, {54, 21},
                {54, 22}, {54, 23},
                {54, 24}, {54, 25},
                {54, 26}, {54, 27},
                {54, 28}, {54, 29},
                {53, 20}, {53, 21},
                {53, 22}, {53, 23},
                {53, 24}, {53, 25},
                {53, 26}, {53, 27},
                {53, 28}, {53, 29},
                {119, 31}, {119, 32},
                {119, 33}, {119, 34},
                {119, 35}, {119, 36},
                {119, 37}, {119, 38},
                {119, 39}, {119, 40},
                {118, 31}, {118, 32},
                {118, 33}, {118, 34},
                {118, 35}, {118, 36},
                {118, 37}, {118, 38},
                {118, 39}, {118, 40},
                {117, 31}, {117, 32},
                {117, 33}, {117, 34},
                {117, 35}, {117, 36},
                {117, 37}, {117, 38},
                {117, 39}, {117, 40}
            };
            AnalyzePositionVsFitParameters(RunsToAnalyze, "position_analysis");

            return 0;
        }

        // Define run ranges to process
        const std::vector<std::pair<Int_t, Int_t> > RunsToProcess = {
            {119, 31}, {119, 32},
            {119, 33}, {119, 34},
            {119, 35}, {119, 36},
            {119, 37}, {119, 38},
            {119, 39}, {119, 40},
        };

        Int_t ProcessedFiles = 0;
        constexpr Int_t MaxFilesToProcess = 100;

        for (const auto &[RunNumber, SubRunNumber]: RunsToProcess)
        {
            if (ProcessedFiles >= MaxFilesToProcess)
            {
                std::cout << "Reached maximum number of files to process ("
                        << MaxFilesToProcess << ")" << std::endl;
                break;
            }

            // Construct input filename
            std::ostringstream InputFileName;
            InputFileName << "pixie_bigrips_traces_"
                    << std::setfill('0') << std::setw(3) << RunNumber
                    << "_"
                    << std::setfill('0') << std::setw(2) << SubRunNumber
                    << ".root";

            std::cout << "\nProcessing file: " << InputFileName.str() << std::endl;

            // Open input file and get tree
            TFile *InputFile = OpenRootFile(InputFileName.str().c_str());
            TTree *Tree = GetTree(InputFile, "pspmt");

            // Create output directory for trace images
            const auto OutputDirectory = CreateTraceDirectory({RunNumber, SubRunNumber});

            // Get qualifying events
            const std::vector<Long64_t> QualifyingEvents = GetAllQualifyingEvents(Tree);

            if (QualifyingEvents.empty())
            {
                std::cout << "No qualifying events found in "
                        << InputFileName.str() << std::endl;
                InputFile->Close();
                delete InputFile;
                continue;
            }

            // Store analysis results
            std::vector<AnalysisResults> Results;
            Results.reserve(QualifyingEvents.size());

            // Process all qualifying events
            std::cout << "Processing " << QualifyingEvents.size() << " qualifying events..." << std::endl;

            int EventCounter = 0;

            for (Long64_t i = 0; i < static_cast<Long64_t>(QualifyingEvents.size()); i++) // static_cast<Long64_t>(QualifyingEvents.size())
            {
                EventCounter++;

                const Long64_t EventNumber = QualifyingEvents[i];

                if (EventCounter % 1000 == 0)
                {
                    std::cout << "Processing event " << EventCounter << " of "
                            << QualifyingEvents.size() << "..." << std::endl;
                }

                auto EventResults = GetEventFitParameters(Tree, EventNumber);
                if (EventResults)
                {
                    Results.push_back(*EventResults);
                }

                // Break after processing 1000 events, for testing purposes
                /*if (EventCounter >= 1000)
                {
                    std::cout << "Finished processing event " << EventCounter << std::endl;
                    break;
                }*/
            }
            std::cout << "\nFinished processing events." << std::endl;

            // Create graphs for a subset of events
            constexpr Long64_t EventsToGraph = 100;
            std::cout << "Graphing first " << EventsToGraph << " qualifying events..." << std::endl;
            GraphFirstNEvents(Tree, QualifyingEvents, EventsToGraph, OutputDirectory.c_str());

            // Save results to output ROOT file
            SaveAnalysisResults(Results, RunNumber, SubRunNumber);

            InputFile->Close();
            delete InputFile;

            ProcessedFiles++;
        }

        // Perform position analysis on all processed runs
        std::vector<std::pair<Int_t, Int_t> > RunsToAnalyze = RunsToProcess;
        AnalyzePositionVsFitParameters(RunsToAnalyze, "position_analysis");

        return 0;
    }
    catch (const std::exception &Error)
    {
        std::cerr << "Error: " << Error.what() << std::endl;
        return 1;
    }
}
