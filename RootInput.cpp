#include <iomanip>
#include <regex>

#include <TROOT.h>
#include <TFile.h>
#include <TTreeReader.h>
#include <TSystem.h>

#include "main.h"

/**
 * Loads required ROOT libraries for the analysis
 * Throws runtime_error if any library fails to load
 */
void LoadRequiredLibraries()
{
    const char* RequiredLibraries[] = {
        "libPaassRootStruct.so",
        "libyaml-cpp.so",
        "libTraceAnalyzerLib.so",
        "libmerger_data_dic.so",
        "libMergerLib.so"
    };

    for (const auto& Library : RequiredLibraries)
    {
        if (gSystem->Load(Library) < 0)
        {
            throw std::runtime_error("Failed to load library: " + std::string(Library));
        }
    }
}

/**
 * Opens a ROOT file and returns the file pointer
 * @param FileName Path to the ROOT file
 * @return Pointer to the opened TFile
 * @throws runtime_error if file cannot be opened
 */
TFile* OpenRootFile(const char* FileName)
{
    // /home/aaugustyn/data/FileName
    auto* path = "/home/aaugustyn/data/";
    TFile* InputFile = TFile::Open((std::string(path) + FileName).c_str());

    if (!InputFile || InputFile->IsZombie())
    {
        throw std::runtime_error("Failed to open file: " + std::string(FileName));
    }

    return InputFile;
}

/**
 * Retrieves a tree from an open ROOT file
 * @param InputFile Pointer to the open ROOT file
 * @param TreeName Name of the tree to retrieve
 * @return Pointer to the retrieved TTree
 * @throws runtime_error if tree cannot be found
 */
TTree* GetTree(TFile* InputFile, const char* TreeName)
{
    if (!InputFile)
    {
        throw std::runtime_error("Invalid file pointer");
    }

    const auto Tree = dynamic_cast<TTree*>(InputFile->Get(TreeName));

    if (!Tree)
    {
        throw std::runtime_error("Failed to retrieve tree: " + std::string(TreeName));
    }

    return Tree;
}

/**
 * Extracts run numbers from a ROOT file name
 * Expected format: pixie_bigrips_traces_XXX_YY.root
 * where XXX is the main run number and YY is the sub-run number
 *
 * @param FileName Name of the ROOT file
 * @return std::pair<Int_t, Int_t> where first is main run number and second is sub-run number
 * @throws std::runtime_error if filename doesn't match expected format
 */
std::pair<Int_t, Int_t> ExtractRunNumbers(const std::string& FileName)
{
    // Regular expression to match the expected filename pattern
    const std::regex Pattern(R"(pixie_bigrips_traces_(\d{3})_(\d{2})\.root)");
    std::smatch Matches;

    if (!std::regex_match(FileName, Matches, Pattern))
    {
        throw std::runtime_error("Invalid filename format. Expected: pixie_bigrips_traces_XXX_YY.root");
    }

    // Matches[0] is the full match
    // Matches[1] is the first capture group (main run number)
    // Matches[2] is the second capture group (sub-run number)
    try
    {
        const Int_t MainRun = std::stoi(Matches[1].str());
        const Int_t SubRun = std::stoi(Matches[2].str());

        // Validate the extracted numbers
        if (MainRun < 0 || SubRun < 0)
        {
            throw std::runtime_error("Extracted negative run numbers");
        }

        return std::make_pair(MainRun, SubRun);
    }
    catch (const std::exception& Error)
    {
        throw std::runtime_error("Failed to convert extracted numbers to integers");
    }
}

/**
 * Creates a directory for storing trace graphs based on run numbers
 * Directory format: Traces_XXX_YY where XXX is main run and YY is sub-run
 *
 * @param RunNumbers Pair of run numbers (MainRun, SubRun)
 * @return std::string Path to the created directory
 * @throws std::runtime_error if directory creation fails
 */
std::string CreateTraceDirectory(const std::pair<Int_t, Int_t>& RunNumbers)
{
    std::ostringstream DirectoryName;
    DirectoryName << "Traces_"
                 << std::setfill('0') << std::setw(3) << RunNumbers.first  // Main run with leading zeros
                 << "_"
                 << std::setfill('0') << std::setw(2) << RunNumbers.second;  // Sub run with leading zeros

    // Create directory using ROOT's gSystem (works across platforms)
    if (gSystem->mkdir(DirectoryName.str().c_str(), kTRUE) != 0)  // kTRUE enables recursive creation
    {
        // Check if directory already exists
        if (gSystem->AccessPathName(DirectoryName.str().c_str()) != 0)  // Returns 0 if path exists
        {
            throw std::runtime_error("Failed to create directory: " + DirectoryName.str());
        }
    }

    return DirectoryName.str();
}