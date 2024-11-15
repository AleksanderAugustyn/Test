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

        // Define run ranges to process
        /*const std::vector<std::pair<Int_t, Int_t>> RunsToProcess = {
            {55, 24}, {55, 25}, {119, 35}, {119, 36}
        };*/
        const std::vector<std::pair<Int_t, Int_t>> RunsToProcess = {
            {119, 35}
        };


        Int_t ProcessedFiles = 0;
        constexpr Int_t MaxFilesToProcess = 10;

        for (const auto& [RunNumber, SubRunNumber] : RunsToProcess)
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
            TFile* InputFile = OpenRootFile(InputFileName.str().c_str());
            TTree* Tree = GetTree(InputFile, "pspmt");

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
        std::vector<std::pair<Int_t, Int_t>> RunsToAnalyze = RunsToProcess;
        AnalyzePositionVsFitParameters(RunsToAnalyze, "position_analysis");

        return 0;
    }
    catch (const std::exception& Error)
    {
        std::cerr << "Error: " << Error.what() << std::endl;
        return 1;
    }
}