#include <TF1.h>
#include <TGraph.h>
#include <iostream>

// Define polynomial coefficients struct for rise power functions
struct RisePowerCoefficients
{
    Double_t Offset; // p[0]
    Double_t Linear; // p[1]
    Double_t Quadratic; // p[2]
    Double_t Cubic; // p[3]
    Double_t Quartic; // p[4]
    Double_t Center; // p[5]
};

// Channel-specific rise power polynomial coefficients
const std::map<std::string, RisePowerCoefficients> ChannelRisePowerFits = {
    {"xa", {1.178, 0.1166, 7.657, 31.42, 1222.0, 0.2509}},
    {"xb", {1.177, 0.1157, 8.094, 33.55, 1201.0, 0.2511}},
    {"ya", {1.178, 0.1174, 7.860, 32.89, 1212.0, 0.2510}},
    {"yb", {1.178, 0.1177, 7.531, 33.51, 1235.0, 0.2512}}
};

Double_t CalculateRisePower(const std::string& Channel, const Double_t Position)
{
    if (ChannelRisePowerFits.find(Channel) == ChannelRisePowerFits.end())
    {
        throw std::runtime_error("Invalid channel name");
    }

    // Validate X position range
    if (Position < 0.1 || Position > 0.4)
    {
        throw std::runtime_error("Invalid X position: " + std::to_string(Position));
    }

    const auto&[Offset, Linear, Quadratic, Cubic, Quartic, Center] = ChannelRisePowerFits.at(Channel);
    const Double_t X = Position - Center;

    return Offset +
           Linear * X +
           Quadratic * std::pow(X, 2) +
           Cubic * std::pow(X, 3) +
           Quartic * std::pow(X, 4);
}

/**
 * Function defining the peak shape for fitting the anode trace
 * p[0] = amplitude
 * p[1] = peak position
 * p[2] = decay constant (tau1)
 * p[3] = rise time constant (tau2)
 * p[4] = rise time power
 * p[5] = baseline
 */
Double_t AnodePeakFunction(const Double_t *X, const Double_t *Parameters)
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
 * Modified version of FitPeakToTrace that uses position information to constrain rise power
 * @param TraceGraph Pointer to the graph containing the trace
 * @param FitRangeStart Start of the fit range
 * @param FitRangeEnd End of the fit range
 * @param Channel Channel identifier (xa, xb, ya, yb)
 * @param Position X or Y position depending on the channel
 * @return Pointer to the fitted function
 */
TF1* FitPeakToTrace(TGraph* TraceGraph, const Double_t FitRangeStart,
                    const Double_t FitRangeEnd, const std::string& Channel,
                    const Double_t Position)
{
    if (!TraceGraph)
    {
        throw std::runtime_error("Invalid graph pointer");
    }

    // Create the fit function
    static int FitCounter = 0;
    const TString FitName = TString::Format("PeakFit_%d", FitCounter++);
    const auto FitFunc = new TF1(FitName, AnodePeakFunction, FitRangeStart, FitRangeEnd, 6);

    // Set parameter names
    FitFunc->SetParName(0, "Amplitude");
    FitFunc->SetParName(1, "PeakPosition");
    FitFunc->SetParName(2, "DecayConstant");
    FitFunc->SetParName(3, "RiseTimeConstant");
    FitFunc->SetParName(4, "RisePower");
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

    // Calculate baseline RMS for better constraints
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
    const Double_t RiseThreshold = BaselineValue + 10 * BaselineRMS;

    for (Int_t i = 0; i < TraceGraph->GetN(); i++)
    {
        Double_t X, Y;
        TraceGraph->GetPoint(i, X, Y);

        if (!RiseStartFound && Y > RiseThreshold)
        {
            RiseStartX = X;
            RiseStartFound = true;
        }

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

    // Calculate expected rise power from position
    const Double_t ExpectedRisePower = CalculateRisePower(Channel, Position);
    constexpr Double_t RisePowerTolerance = 0.05; // Allow some variation around the expected value

    // Set initial parameters with improved estimates
    FitFunc->SetParameter(0, MaxY - BaselineValue);  // Amplitude
    FitFunc->SetParameter(1, MaxX);                  // Peak position
    FitFunc->SetParameter(2, EstimatedDecayConstant);// Decay constant
    FitFunc->SetParameter(3, EstimatedRiseTime);     // Rise time constant
    FitFunc->SetParameter(4, ExpectedRisePower);     // Rise time power (from position)
    FitFunc->SetParameter(5, BaselineValue);         // Baseline

    // Set parameter limits
    FitFunc->SetParLimits(0, 0.5 * (MaxY - BaselineValue), 1.2 * MaxY);
    FitFunc->SetParLimits(1, MaxX - 30, MaxX + 30);
    FitFunc->SetParLimits(2, 1, 200);
    FitFunc->SetParLimits(3, 1, 100);
    // Constrain rise power around the expected value from position
    FitFunc->SetParLimits(4, ExpectedRisePower - RisePowerTolerance,
                            ExpectedRisePower + RisePowerTolerance);
    FitFunc->SetParLimits(5, BaselineValue - 5*BaselineRMS, BaselineValue + 5*BaselineRMS);

    // Perform the fit
    [[maybe_unused]] const Int_t FitStatus = TraceGraph->Fit(FitFunc, "QR");

    // Set fit function style
    FitFunc->SetLineColor(kRed);
    FitFunc->SetLineWidth(100);
    FitFunc->SetNpx(2000);

    return FitFunc;
}

/**
 * Function to fit dynode PMT signals with double exponential decay and undershoot
 * Parameters:
 * p[0] = amplitude
 * p[1] = peak position (t0)
 * p[2] = fast decay time constant (τ1)
 * p[3] = slow decay time constant (τ2)
 * p[4] = rise time constant (τr)
 * p[5] = undershoot amplitude
 * p[6] = undershoot recovery time constant (τu)
 * p[7] = relative weight of fast component
 * p[8] = baseline
 */
Double_t DynodePeakFunction(const Double_t *X, const Double_t *Parameters)
{
    const Double_t T = X[0] - Parameters[1]; // Time relative to peak

    // Before peak, return baseline
    if (T < 0)
    {
        return Parameters[8];
    }

    // Rise time component (error function for Gaussian response)
    const Double_t Rise = 1.0 - TMath::Exp(-T / Parameters[4]);

    // Double exponential decay
    const Double_t FastDecay = Parameters[7] * TMath::Exp(-T / Parameters[2]);
    const Double_t SlowDecay = (1.0 - Parameters[7]) * TMath::Exp(-T / Parameters[3]);
    const Double_t Decay = FastDecay + SlowDecay;

    // Undershoot and recovery
    const Double_t Undershoot = Parameters[5] * (1.0 - TMath::Exp(-T / Parameters[6]));

    // Combine all components
    return Parameters[8] + Parameters[0] * Rise * Decay - Undershoot;
}

/**
 * Fits the dynode peak function to a trace
 * @param TraceGraph Pointer to the graph containing the trace
 * @param FitRangeStart Start of the fit range
 * @param FitRangeEnd End of the fit range
 * @return Pointer to the fitted function
 */
TF1 *FitDynodePeak(TGraph *TraceGraph, const Double_t FitRangeStart, const Double_t FitRangeEnd)
{
    if (!TraceGraph)
    {
        throw std::runtime_error("Invalid graph pointer");
    }

    // Create the fit function
    static Int_t FitCounter = 0;
    const TString FitName = TString::Format("DynodeFit_%d", FitCounter++);
    const auto FitFunc = new TF1(FitName, DynodePeakFunction, FitRangeStart, FitRangeEnd, 9);

    // Set parameter names
    FitFunc->SetParName(0, "Amplitude");
    FitFunc->SetParName(1, "PeakPosition");
    FitFunc->SetParName(2, "FastDecay");
    FitFunc->SetParName(3, "SlowDecay");
    FitFunc->SetParName(4, "RiseTime");
    FitFunc->SetParName(5, "UndershootAmp");
    FitFunc->SetParName(6, "UndershootRecovery");
    FitFunc->SetParName(7, "FastFraction");
    FitFunc->SetParName(8, "Baseline");

    // Find initial parameter estimates
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

    // Find peak
    for (Int_t i = 0; i < TraceGraph->GetN(); i++)
    {
        Double_t X, Y;
        TraceGraph->GetPoint(i, X, Y);
        if (Y > MaxY)
        {
            MaxY = Y;
            MaxX = X;
        }
    }

    // Estimate undershoot
    Double_t MinAfterPeak = 1e9;
    [[maybe_unused]] Double_t MinAfterPeakX = 0;
    for (auto i = static_cast<Int_t>(MaxX + 100); i < TraceGraph->GetN(); i++)
    {
        Double_t X, Y;
        TraceGraph->GetPoint(i, X, Y);
        if (Y < MinAfterPeak)
        {
            MinAfterPeak = Y;
            MinAfterPeakX = X;
        }
    }

    // Set initial parameters
    FitFunc->SetParameter(0, MaxY - BaselineValue); // Amplitude
    FitFunc->SetParameter(1, MaxX); // Peak position
    FitFunc->SetParameter(2, 20.0); // Fast decay (~10ns)
    FitFunc->SetParameter(3, 40.0); // Slow decay (~40ns)
    FitFunc->SetParameter(4, 3.0); // Rise time
    FitFunc->SetParameter(5, BaselineValue - MinAfterPeak); // Undershoot amplitude
    FitFunc->SetParameter(6, 500.0); // Undershoot recovery
    FitFunc->SetParameter(7, 2.0); // Fast fraction
    FitFunc->SetParameter(8, BaselineValue); // Baseline

    // Set parameter limits
    FitFunc->SetParLimits(0, 0.5 * (MaxY - BaselineValue), 2.5 * (MaxY - BaselineValue));
    FitFunc->SetParLimits(1, MaxX - 50, MaxX + 50);
    FitFunc->SetParLimits(2, 1.0, 100.0); // Fast decay
    FitFunc->SetParLimits(3, 10.0, 200.0); // Slow decay
    FitFunc->SetParLimits(4, 0.5, 20.0); // Rise time
    FitFunc->SetParLimits(5, 0.0, BaselineValue - MinAfterPeak + 300);
    FitFunc->SetParLimits(6, 50.0, 1000.0); // Undershoot recovery
    FitFunc->SetParLimits(7, 0.0, 50.0); // Fast fraction
    FitFunc->SetParLimits(8, BaselineValue - 100, BaselineValue + 100);

    // Perform the fit
    TraceGraph->Fit(FitFunc, "QR"); // Q: quiet, R: use range

    return FitFunc;
}
