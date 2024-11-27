#include <TTreeReader.h>

#include "PaassRootStruct.hpp"

#include "main.h"

/**
 * Scans all events in a tree and returns selected event numbers
 * @param Tree Pointer to the input tree
 * @param MaxEventsToSave Maximum number of events to save
 * @return Vector of selected event numbers
 */
std::vector<Long64_t> ScanEvents(TTree* Tree, const Long64_t MaxEventsToSave)
{
    if (!Tree)
    {
        throw std::runtime_error("Invalid tree pointer");
    }

    const Long64_t Entries = Tree->GetEntries();
    std::vector<Long64_t> SelectedEventNumbers;
    Long64_t SavedEvents = 0;
    Long64_t TotalQualifyingEvents = 0;

    std::cout << "Scanning " << Entries << " events..." << std::endl;

    for (Long64_t Event = 0; Event < Entries; Event++)
    {
        /*if (Event % 10000 == 0)
        {
            std::cout << "Scanning event " << Event << "/" << Entries << "\r" << std::flush;
        }*/

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

    std::cout << "\nTotal events meeting criteria: " << TotalQualifyingEvents << std::endl;
    std::cout << "Selected " << SavedEvents << " events for processing" << std::endl;

    return SelectedEventNumbers;
}

Bool_t MeetsSelectionCriteria(TTree* TreeInput, const Long64_t Entry)
{
    if (!TreeInput)
    {
        throw std::runtime_error("Invalid tree pointer");
    }

    // Setup branch variables for previous criteria
    Int_t HighGainValid;
    Int_t LowGainValid;
    Double_t HighGainQdc;
    Double_t PosX;
    Double_t PosY;

    // Set branch addresses
    TreeInput->SetBranchAddress("high_gain_.valid_", &HighGainValid);
    TreeInput->SetBranchAddress("low_gain_.valid_", &LowGainValid);
    TreeInput->SetBranchAddress("high_gain_.qdc_", &HighGainQdc);
    TreeInput->SetBranchAddress("high_gain_.pos_x_", &PosX);
    TreeInput->SetBranchAddress("high_gain_.pos_y_", &PosY);

    // Get the values for this entry
    TreeInput->GetEntry(Entry);

    // Check previous conditions first
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

    // Check position range
    constexpr Double_t MinPos = 0.1;
    constexpr Double_t MaxPos = 0.4;

    if (PosX < MinPos || PosX > MaxPos || PosY < MinPos || PosY > MaxPos)
    {
        return false;
    }

    // Set up reader for rootdev_vec_ branch
    TTreeReader Reader;
    Reader.SetTree(TreeInput);
    TTreeReaderValue<std::vector<processor_struct::ROOTDEV>> RootDevVector(Reader, "rootdev_vec_");

    // Get entry
    Reader.SetEntry(Entry);

    // Variables to track if we found all required channels
    Bool_t FoundXa = false;
    Bool_t FoundXb = false;
    Bool_t FoundYa = false;
    Bool_t FoundYb = false;
    Bool_t FoundDynodeHigh = false;

    // Check each device in the vector
    for (const auto& Device : *RootDevVector)
    {
        // Check for valid timing and waveform analysis
        if (!Device.hasValidTimingAnalysis || !Device.hasValidWaveformAnalysis)
        {
            continue;
        }

        // Check for anode channels
        if (Device.subtype == "anode_high")
        {
            switch (Device.chanNum)
            {
                case 4:  // xa
                    FoundXa = true;
                    break;
                case 7:  // xb
                    FoundXb = true;
                    break;
                case 6:  // ya
                    FoundYa = true;
                    break;
                case 5:  // yb
                    FoundYb = true;
                    break;
                default:
                    break;
            }
        }
        // Check for dynode
        else if (Device.subtype == "dynode_high")
        {
            FoundDynodeHigh = true;
        }
    }

    // Event is selected only if all required channels are present
    return FoundXa && FoundXb && FoundYa && FoundYb && FoundDynodeHigh;
}

std::vector<Long64_t> GetAllQualifyingEvents(TTree* TreeInput)
{
    if (!TreeInput)
    {
        throw std::runtime_error("Invalid tree pointer");
    }

    const Long64_t Entries = TreeInput->GetEntries();
    std::vector<Long64_t> QualifyingEvents;

    std::cout << "Scanning " << Entries << " events for qualification..." << std::endl;

    int EventCounter = 1;

    for (Long64_t Event = 0; Event < Entries; Event++)
    {


        /*if (EventCounter % 10000 == 0)
        {
            std::cout << "Processing event " << Event << "/" << Entries << "\r" << std::flush;
        }*/

        if (MeetsSelectionCriteria(TreeInput, Event))
        {
            QualifyingEvents.push_back(Event);
        }

        EventCounter++;
    }

    std::cout << "\nFound " << QualifyingEvents.size() << " qualifying events" << std::endl;
    return QualifyingEvents;
}

