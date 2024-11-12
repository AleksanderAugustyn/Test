#pragma once


#include <TFile.h>
#include <TTree.h>
#include <string>
#include <memory>

class ROOTFileReader
{
public:
    // Constructor
    ROOTFileReader();

    // Destructor
    ~ROOTFileReader();

    // Open ROOT file
    bool OpenFile(const std::string &filename);

    // Close ROOT file
    void CloseFile();

    // Get pointer to the TTree named "pspmt"
    [[nodiscard]] TTree *GetPSPMTTree() const;

    // Check if file is open
    [[nodiscard]] bool IsOpen() const;

private:
    static const std::string DATA_PATH;
    std::unique_ptr<TFile> ROOTFile;
    TTree *PSPMTTree;
};
