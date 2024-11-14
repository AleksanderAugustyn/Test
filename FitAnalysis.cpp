#include <TF1.h>
#include <TGraph.h>
#include <iostream>

/**
 * Function defining the peak shape for fitting
 * p[0] = amplitude
 * p[1] = peak position
 * p[2] = decay constant (tau1)
 * p[3] = rise time constant (tau2)
 * p[4] = rise time power
 * p[5] = baseline
 */
Double_t PeakFunction(const Double_t* X, const Double_t* Parameters)
{
    // Before peak, return baseline
    if (X[0] <= Parameters[1])
    {
        return Parameters[5];
    }

    const Double_t TimeOffset = X[0] - Parameters[1];

    // Improved rise time calculation with variable power
    const Double_t Rise = 1.0 - TMath::Exp(-TMath::Power(TimeOffset / Parameters[3], Parameters[4]));

    // Exponential decay
    const Double_t Decay = TMath::Exp(-TimeOffset / Parameters[2]);

    // Combine components
    return Parameters[5] + Parameters[0] * Rise * Decay;
}

/**
 * Fits the peak function to a trace and adds the fit to the graph
 * @param TraceGraph Pointer to the graph containing the trace
 * @param FitRangeStart Start of the fit range
 * @param FitRangeEnd End of the fit range
 * @return Pointer to the fitted function
 */
TF1* FitPeakToTrace(TGraph* TraceGraph, const Double_t FitRangeStart, const Double_t FitRangeEnd)
{
    if (!TraceGraph)
    {
        throw std::runtime_error("Invalid graph pointer");
    }

    // std::cout << "\nFitting trace: " << TraceGraph->GetTitle() << std::endl;
    // std::cout << "Fit range: [" << FitRangeStart << ", " << FitRangeEnd << "]" << std::endl;

    // Create the fit function
    static int FitCounter = 0;
    const TString FitName = TString::Format("PeakFit_%d", FitCounter++);
    const auto FitFunc = new TF1(FitName, PeakFunction, FitRangeStart, FitRangeEnd, 6);

    // Set parameter names
    FitFunc->SetParName(0, "Amplitude");
    FitFunc->SetParName(1, "PeakPosition");
    FitFunc->SetParName(2, "DecayConstant");
    FitFunc->SetParName(3, "RiseTimeConstant");
    FitFunc->SetParName(4, "RiseTimePower");
    FitFunc->SetParName(5, "Baseline");

    // Find approximate parameters
    Double_t MaxY = -1e9;
    Double_t MaxX = 0;
    Double_t BaselineValue = 0;
    Int_t BaselineSamples = 0;

    // Use first ~20 points for baseline estimation
    constexpr Int_t BaselinePoints = 20;
    for (Int_t i = 0; i < TMath::Min(BaselinePoints, TraceGraph->GetN()); i++)
    {
        Double_t X, Y;
        TraceGraph->GetPoint(i, X, Y);
        BaselineValue += Y;
        BaselineSamples++;
    }
    BaselineValue /= BaselineSamples;

    // Calculate baseline RMS for better baseline constraints
    Double_t BaselineRMS = 0;
    for (Int_t i = 0; i < BaselineSamples; i++)
    {
        Double_t X, Y;
        TraceGraph->GetPoint(i, X, Y);
        BaselineRMS += (Y - BaselineValue) * (Y - BaselineValue);
    }
    BaselineRMS = TMath::Sqrt(BaselineRMS / BaselineSamples);

    // Find peak and estimate rise time
    Double_t RiseStartX = 0;
    Bool_t RiseStartFound = false;
    const Double_t RiseThreshold = BaselineValue + 5 * BaselineRMS;

    for (Int_t i = 0; i < TraceGraph->GetN(); i++)
    {
        Double_t X, Y;
        TraceGraph->GetPoint(i, X, Y);

        // Find rise start
        if (!RiseStartFound && Y > RiseThreshold)
        {
            RiseStartX = X;
            RiseStartFound = true;
        }

        // Track maximum
        if (Y > MaxY)
        {
            MaxY = Y;
            MaxX = X;
        }
    }

    // Estimate decay time constant
    Double_t DecayEndX = MaxX;
    for (auto i = static_cast<Int_t>(MaxX); i < TraceGraph->GetN(); i++)
    {
        Double_t X, Y;
        TraceGraph->GetPoint(i, X, Y);
        if (Y <= BaselineValue + (MaxY - BaselineValue) * 1/M_E)  // 1/e decay point
        {
            DecayEndX = X;
            break;
        }
    }
    const Double_t EstimatedDecayConstant = DecayEndX - MaxX;
    const Double_t EstimatedRiseTime = MaxX - RiseStartX;

    // std::cout << "Initial parameter estimation:" << std::endl;
    // std::cout << "  Max amplitude: " << MaxY - BaselineValue << " at position: " << MaxX << std::endl;
    // std::cout << "  Baseline: " << BaselineValue << " ± " << BaselineRMS << std::endl;
    // std::cout << "  Estimated rise time: " << EstimatedRiseTime << std::endl;
    // std::cout << "  Estimated decay constant: " << EstimatedDecayConstant << std::endl;

    // Set initial parameters with improved estimates
    FitFunc->SetParameter(0, MaxY - BaselineValue);  // Amplitude
    FitFunc->SetParameter(1, MaxX);                  // Peak position
    FitFunc->SetParameter(2, EstimatedDecayConstant);// Decay constant
    FitFunc->SetParameter(3, EstimatedRiseTime);     // Rise time constant
    FitFunc->SetParameter(4, 2.0);                   // Rise time power (start between linear and quadratic)
    FitFunc->SetParameter(5, BaselineValue);         // Baseline

    // Set parameter limits
    FitFunc->SetParLimits(0, 0.5 * (MaxY - BaselineValue), 1.5 * (MaxY - BaselineValue));  // Amplitude
    FitFunc->SetParLimits(1, MaxX - 25, MaxX + 25);        // Peak position
    FitFunc->SetParLimits(2, 1, 100);                     // Decay constant (reasonable range for PMT)
    FitFunc->SetParLimits(3, 1, 50);                       // Rise time constant
    FitFunc->SetParLimits(4, 1.0, 4.0);                    // Rise time power
    FitFunc->SetParLimits(5, BaselineValue - 5*BaselineRMS, BaselineValue + 5*BaselineRMS);  // Baseline

    // Perform the fit
    [[maybe_unused]] const Int_t FitStatus = TraceGraph->Fit(FitFunc, "QR");  // Q: quiet, R: use range

    /*// Print fit results
     std::cout << "Fit status: " << FitStatus << " (0 = success)" << std::endl;
    std::cout << "Fitted parameters:" << std::endl;
    for (Int_t i = 0; i < FitFunc->GetNpar(); i++)
    {
        std::cout << "  " << FitFunc->GetParName(i) << ": "
                 << FitFunc->GetParameter(i) << " +/- "
                 << FitFunc->GetParError(i) << std::endl;
    }

    // Calculate fit quality metrics
    const Double_t ChiSquare = FitFunc->GetChisquare();
    const Int_t NDF = FitFunc->GetNDF();
    const Double_t RedChiSquare = ChiSquare / NDF;
    std::cout << "Fit quality:" << std::endl;
    std::cout << " Chi-square: " << ChiSquare << std::endl;
    std::cout << " NDF: " << NDF << std::endl;
    std::cout << " Reduced chi-square: " << RedChiSquare << std::endl;*/

    // Set fit function style
    FitFunc->SetLineColor(kRed);
    FitFunc->SetLineWidth(100);
    FitFunc->SetNpx(2000);

    return FitFunc;
}