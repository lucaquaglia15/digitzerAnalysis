// PeakFinder.h
#ifndef PEAKFINDER_H
#define PEAKFINDER_H

#include <vector>
#include <utility>
#include <tuple>
#include "TMath.h"
#include "TLine.h"

namespace PeakFinder {
    const float EPS = 2.2204e-16f;

    /*
        Inputs
        x0: input signal        
        extrema: 1 if maxima are desired, -1 if minima are desired
        includeEndpoints - If true the endpoints will be included as possible extrema otherwise they will not be included
        Output
        peakInds: Indices of peaks in x0
    */
    //void findPeaks(std::vector<float> x0, std::vector<int>& peakInds, bool includeEndpoints=true, float extrema=1);
    int findPeaks2(std::vector<float> x0, std::vector<float>& avgNoise, std::vector<float>& avg, std::vector<int>& peakInds, bool includeEndpoints=true, int threshold=5);
    int findPeaks(std::vector<float> x0, std::vector<int>& peakInds, bool includeEndpoints=true, int threshold=5);
    //auto findPeaks2(std::vector<float> x0, std::vector<int>& peakInds, bool includeEndpoints=true, int threshold=5);
    int findPeaksSmoothed(std::vector<float> x0, std::vector<float> meanNoise, std::vector<float> mean, std::vector<int>& peakInds, bool includeEndpoints=true, int threshold=5);
    //double FindIntegrationInterval(double avgNoise, int signalLength, int muonStart, vector<float> datamV);
    double FindIntegrationInterval(int event, int strip, double avgNoise, int signalLength, int muonStart, vector<float> datamV);
    //double FindIntegrationInterval2(int event, int strip, double avgNoise, int signalLength, int muonStart, vector<float> datamV);
    //pair<int,int> FindIntegrationInterval2(int event, int strip, double avgNoise, int signalLength, int muonStart, vector<float> datamV);
    //pair<int,int> FindIntegrationInterval2(int event, int strip, double avgNoise, int signalLength, int muonStart, vector<float> datamV, vector<TLine*> &testMark);
    tuple<vector<int>,vector<int>,int,double> FindIntegrationInterval2(int event, int strip, double avgNoise, int signalLength, int muonStart, vector<float> datamV);
    double IntegrateMuonPeak(int signStart, int signEnd, vector<float> datamV);
}

#endif
