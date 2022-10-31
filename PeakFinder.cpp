// PeakFinder.cpp
#include "PeakFinder.h"

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <numeric>

//My implementation
int PeakFinder::findPeaks(std::vector<float> x0, std::vector<int>& peakInds, bool includeEndpoints, int threshold) //old version returning only number of peaks
{
	int peakStart = 0, peakEnd = 0, numPeaks = 0; //Beginning and ending of the peak (in terms of vector indeces)
	bool mightBePeak = false;
	float meanNoise = 0, sq_sum = 0, stdev = 0;
	unsigned int noiseMin = 0, noiseMax = 150;
	int noiseWindow = 150;

	//checks when the values go over and when they go under the threshold value
	unsigned int inputVectSize = x0.size();
	for (unsigned int i = 0; i < inputVectSize; i+= 1024) { //Cycle through input vector
		//cout << "Event " << i/1024 << endl;
		//cout << x0.at(noiseMin+i) << "\t" << x0.at(noiseMax+i) << endl;
		meanNoise = std::accumulate(&x0.at(noiseMin+i), &x0.at(noiseMax+i), 0.0) / (noiseWindow);
		sq_sum = inner_product(&x0.at(noiseMin+i), &x0.at(noiseMax+i), &x0.at(noiseMin+i),0.0);
		stdev = std::sqrt(sq_sum/(noiseWindow) - meanNoise*meanNoise);
		//cout << "Mean noise " << meanNoise << "\t sq_sum " << sq_sum << "\t stdev " << stdev << endl;

		for (int j = 0; j < 1024; j++) {

			if ((x0[j+i] > meanNoise + threshold*stdev) && mightBePeak == false) { //Element of the vector is above threshold
				mightBePeak = true;
				peakStart = j; //assign peak beg index
				//cout << "peak start " << peakStart+i << endl;
			}

			else if (mightBePeak == true && (x0[j+i] < meanNoise + threshold*stdev)) { //End of the peak
				mightBePeak = false;
				numPeaks++;
				//cout << "picco!" << endl;
				peakEnd = j;
				//cout << "peak end " << peakEnd+i << endl;
				//cout << std::max_element(&x0.at(peakStart), &x0.at(peakEnd)) - &x0.at(peakStart) + i << endl;
				peakInds.push_back(std::max_element(&x0.at(peakStart), &x0.at(peakEnd)) - &x0.at(peakStart) + i);
				//peakInds.push_back(*max_element(x0[peakStart],x0[peakEnd]));
			}

			if (j == 1023)
				mightBePeak = false; //between one event and the other the mean might change, so we force no peaks
		}
		meanNoise = 0;
		sq_sum = 0;
		stdev = 0;
	}	
	return numPeaks;
}


//int PeakFinder::findPeaks2(std::vector<float> x0, std::vector<int>& peakInds, bool includeEndpoints, int threshold) //old version returning only number of peaks
int PeakFinder::findPeaks2(std::vector<float> x0, std::vector<float>& avgNoise, std::vector<float>& avg, std::vector<int>& peakInds, bool includeEndpoints, int threshold) //new version w multiple returns
{
	int peakStart = 0, peakEnd = 0, numPeaks = 0; //Beginning and ending of the peak (in terms of vector indeces)
	bool mightBePeak = false;
	float meanNoise = 0, sq_sum = 0, stdev = 0;
	unsigned int noiseMin = 0, noiseMax = 150;
	int noiseWindow = 150;

	//checks when the values go over and when they go under the threshold value
	unsigned int inputVectSize = x0.size();
	for (unsigned int i = 0; i < inputVectSize; i+= 1024) { //Cycle through input vector
		//cout << "Event " << i/1024 << endl;
		//cout << x0.at(noiseMin+i) << "\t" << x0.at(noiseMax+i) << endl;
		meanNoise = std::accumulate(&x0.at(noiseMin+i), &x0.at(noiseMax+i), 0.0) / (noiseWindow);
		sq_sum = inner_product(&x0.at(noiseMin+i), &x0.at(noiseMax+i), &x0.at(noiseMin+i),0.0);
		stdev = std::sqrt(sq_sum/(noiseWindow) - meanNoise*meanNoise);
		avg.push_back(meanNoise); //push back mean oscillation
		avgNoise.push_back(stdev); //push back average noise fluctuation
		//cout << "Mean noise " << meanNoise << "\t sq_sum " << sq_sum << "\t stdev " << stdev << endl;

		for (int j = 0; j < 1024; j++) {

			if ((x0[j+i] > meanNoise + threshold*stdev) && mightBePeak == false) { //Element of the vector is above threshold
				mightBePeak = true;
				peakStart = j; //assign peak beg index
				//cout << "peak start " << peakStart+i << endl;
			}

			else if (mightBePeak == true && (x0[j+i] < meanNoise + threshold*stdev)) { //End of the peak
				mightBePeak = false;
				numPeaks++;
				//cout << "picco!" << endl;
				peakEnd = j;
				//cout << "peak end " << peakEnd+i << endl;
				//cout << std::max_element(&x0.at(peakStart), &x0.at(peakEnd)) - &x0.at(peakStart) + i << endl;
				peakInds.push_back(std::max_element(&x0.at(peakStart), &x0.at(peakEnd)) - &x0.at(peakStart) + i);
				//peakInds.push_back(*max_element(x0[peakStart],x0[peakEnd]));
			}

			if (j == 1023)
				mightBePeak = false; //between one run and the other the mean might change, so we force no peaks
		}
		meanNoise = 0;
		sq_sum = 0;
		stdev = 0;
	}	
	return numPeaks;
	//return retVal{noise,mean,numPeaks};
}

//peak finder in smoothed data, using the same thresholds as non-smoothed data
int PeakFinder::findPeaksSmoothed(std::vector<float> x0, std::vector<float> meanNoise, std::vector<float> mean, std::vector<int>& peakInds, bool includeEndpoints, int threshold)
{
	int peakStart = 0, peakEnd = 0, numPeaks = 0, eventCounter = 0; //Beginning and ending of the peak (in terms of vector indeces)
	bool mightBePeak = false;
	float noise = 0, oscillation = 0;

	//checks when the values go over and when they go under the threshold value
	unsigned int inputVectSize = x0.size();
	for (unsigned int i = 0; i < inputVectSize; i+= 1024) { //Cycle through input vector
		//cout << "Event " << i/1024 << endl;
		//cout << x0.at(noiseMin+i) << "\t" << x0.at(noiseMax+i) << endl;
		noise = meanNoise[eventCounter]; 
		oscillation = mean[eventCounter];
		//cout << "Mean noise " << meanNoise << "\t sq_sum " << sq_sum << "\t stdev " << stdev << endl;

		for (int j = 0; j < 1024; j++) {

			if ((x0[j+i] > oscillation + threshold*noise) && mightBePeak == false) { //Element of the vector is above threshold
				mightBePeak = true;
				peakStart = j; //assign peak beg index
				//cout << "peak start " << peakStart+i << endl;
			}

			else if (mightBePeak == true && (x0[j+i] < oscillation + threshold*noise)) { //End of the peak
				mightBePeak = false;
				numPeaks++;
				//cout << "picco!" << endl;
				peakEnd = j;
				//cout << "peak end " << peakEnd+i << endl;
				//cout << std::max_element(&x0.at(peakStart), &x0.at(peakEnd)) - &x0.at(peakStart) + i << endl;
				peakInds.push_back(std::max_element(&x0.at(peakStart), &x0.at(peakEnd)) - &x0.at(peakStart) + i);
				//peakInds.push_back(*max_element(x0[peakStart],x0[peakEnd]));
			}

			if (j == 1023)
				mightBePeak = false; //between one run and the other the mean might change, so we force no peaks
		}
	}	
	return numPeaks;
}

double PeakFinder::IntegrateMuonPeak(int signStart, int signEnd, vector<float> datamV) { //Calculate integral in range
	double integral = 0.;
	
	for (int i = signStart; i < signEnd; i++) { //Compute integra with trapezoid method
		//Trapzoid area: 0.5*(B+b)*h -> in our 

		//y                 //In our case the base (on x axys) is always 1 -> because of the sampling rate (1024 samples @ 1 Gs/s)
		//|                 //The two bases are two conseutive values of the points y[i] and y[i+1]
		//|	      /|y[i+1]  // In order to calculate the integral I su over all these contributions
		//|  y[i]| | 
		//|    b | | B
		//|______|_|_____x
		//        h

		integral += (datamV[i+1] + datamV[i]);
	}

	return (0.5*integral)/50; //due to the trapezoid rule, I have to sum ((b+B)*h)/2 -> it's the same as dividing the final sum by 0.5 and h = 1
	//Divided by 50 due to the fact that we read the signal on a 50 ohm resistor and also Q = I*t and I = V/R -> Q = (V/R)*t in the integral 
	//pC unit of measure
}

bool IsGreaterThanOne (int i) {
        return (i >= 1);
}

//Go through the datamV vector

//This function is called everytime a peak is found
//double PeakFinder::FindIntegrationInterval(double avgNoise, int signalLength, int muonStart, vector<float> datamV) {
double PeakFinder::FindIntegrationInterval(int event, int strip, double avgNoise, int signalLength, int muonStart, vector<float> datamV) {
	//avgNoise = noise before the peak (5*rms)
	//muonStart = start of muon window, signalLength = 1024 (elements per event)
	//Logic: calulcate integral from the first point where signal/noise > 1 until the las point where this is true
	//For muons we keep a smaller window, usually the signals are around 150-250 ns
	//If signal is short (order of 20 ns) -> we need to implement something else, if the signal is longer this method is quite ok
	int signStart = 0, signEnd = 0; //Start and end of signal/noise > 1
	float signBefore = 0, signAfter = 0; //Keeps track of the ratio of values in order to get the proper integration interval

	//Calculate signal/noise ratio -> absolute value because signals can be negative as well
	vector<float> signalToNoise;
	signalToNoise.reserve(1024);

	for (int i = 0; i < signalLength; i++) {
		signalToNoise.push_back(TMath::Abs(datamV[i]/avgNoise)); 
	}

    cout << endl << "Event # " << event << " stip " << strip << endl;

    vector<float>::iterator itFirst = find_if(signalToNoise.begin()+muonStart, signalToNoise.end(), IsGreaterThanOne); //Find iterator to first S/N > 1
    cout << "The index of first element is: " << std::distance(signalToNoise.begin(), itFirst)   << '\n';
    vector<float>::reverse_iterator itLast = find_if(signalToNoise.rbegin(), signalToNoise.rend(), IsGreaterThanOne); //Find iterator to last S/N > 1
    cout << "The index of last element is: " << std::distance(itLast, signalToNoise.rend()) - 1  << '\n';
    
    //signStart = std::distance(signalToNoise.begin()+muonStart, itFirst);
    signStart = std::distance(signalToNoise.begin(), itFirst);
    signEnd = std::distance(itLast, signalToNoise.rend()) - 1;

    cout << "Int start before " << signStart << " int end before " << signEnd << endl;

    for (int i = signStart; i > 0; i--) { //Find true start of the signal
    	signBefore = datamV[i]/datamV[i-1];
    	cout << "signBefore " << signBefore << endl;
    	if (signBefore >= 0) continue;
    	else {
    		//signStart = datamV[i];
    		signStart = i;
    		break;
    	}
    }

    cout << "------------" << endl;

    /*for (unsigned int i = signEnd; i < datamV.size()-1; i++) { //Find true end of the signal
    	signAfter = datamV[i]/datamV[i+1];
    	cout << "signAfter " << signAfter << endl;
    	if (signAfter >= 0) continue;
    	else {
    		signEnd = i-1;
    		break;
    	}

    }*/

    signEnd = signStart + 200;

    cout << "Int start after " << signStart << " int end after " << signEnd << endl;
 	cout << "Signal integral " << IntegrateMuonPeak(signStart,signEnd,datamV) << endl;	
    return IntegrateMuonPeak(signStart,signEnd,datamV);
}

tuple<vector<int>,vector<int>,int,double> PeakFinder::FindIntegrationInterval2(int event, int strip, double avgNoise, int signalLength, int muonStart, vector<float> datamV) {
//pair<int,int> PeakFinder::FindIntegrationInterval2(int event, int strip, double avgNoise, int signalLength, int muonStart, vector<float> datamV) {
//pair<int,int> PeakFinder::FindIntegrationInterval2(int event, int strip, double avgNoise, int signalLength, int muonStart, vector<float> datamV, vector<TLine*> &testMark) {

	//avgNoise = noise before the first peak (5*rms)
	//muonStart = start of muon window, signalLength = 1024 (elements per event)
	//Logic: calulcate integral from the first point where signal/noise > 1 until the las point where this is true
	//For muons we keep a smaller window, usually the signals are around 150-250 ns
	//If signal is short (order of 20 ns) -> we need to implement something else, if the signal is longer this method is quite ok

	bool mightBePeak = false, firstPeak = false, positivePeak = false, negativePeak = false, firstPeakDone = false, wasPeak = false;
	bool rising = false, falling = false, printOut = false;
	int peakStart = 0, peakEnd = 0, peakEndTemp = 0, peakStartTemp = 0, timeDiff = 0, timeDiffEndPeak = 0, timeDiffStartPeak = 0, peakNum = 1, trigStart = 680, trigEnd = 785;
	double signBefore = 0, ratio = 0, signAfter = 0, avgPeaks = 0, sq_sum_gamma = 0, avgValue = 0.;
	double m = 0, q = 0, intersectionStart = 0, intersectionEnd;
	vector<int> peakValues, peakSign, inizioPicco, finePicco; 
	vector<double> vIntStart, vIntEnd;
	 
	if(printOut) 
		cout << "----" << endl;
	if(printOut) 
		cout << "Event # " << event << " strip " << strip << endl;
	if(printOut) 
		cout << "----" << endl;

	if (printOut) cout << "+5 rms " << avgNoise << " -5 rms " << -avgNoise << endl;

	avgValue = accumulate(&datamV[0],&datamV[muonStart],0.0)/muonStart; //average value of the waveform in noise window
	if(printOut) cout << "Avg value " << avgValue << endl;

	for (unsigned int i = 150; i < datamV.size(); i++) { //Go through the mV vector

		//I distinguish in particular the first peak because with the algorithm I find the value where the signal crosses 5 sigma

		// y
		// |         *
		// |        *  *
		// |      |*|   |*|
		// |------------------- 5 sigma
		// |      *        *
		// |     *          /*/**
		// |*_*/*/_______________ x    |*| point found by the algorithm /*/ real start, reason why I distinguis the first peak and last peak

		//----First peak----//
		if (firstPeakDone == false) {

			if (i > 1000) {
				if (mightBePeak == true && positivePeak == true) {//Data is still above 5 sigma
					peakEnd = 1000;
					finePicco.push_back(peakEnd); //end of peak assumed to be at 1000 ns
					vIntEnd.push_back(peakEnd); //end of time over threshold as well
				}

				else if (mightBePeak == true && negativePeak == true) {//Data is still below -5 sigma
					peakEnd = 1000;
					finePicco.push_back(peakEnd); //end of peak assumed to be at 1000 ns
					vIntEnd.push_back(peakEnd); //end of time over threshold as well
				}

				else if (mightBePeak == false) //Data is not above 5 sigma or below -5 sigma
					peakEnd=peakEndTemp;
			break;
			}

			else if (datamV[i] > avgNoise && firstPeak == false && mightBePeak == false) { //Positive peak
				mightBePeak = true;
				positivePeak = true;
				firstPeak = true;
				if (printOut) cout << "First point over 5rms " << i << endl;
				//Calculate time over threshold starting point
						//----//
				m = (datamV[i]-datamV[i-1])/(i-(i-1)); //Slope of line
				q = datamV[i]-m*i;
				intersectionStart = (avgNoise-q)/m;
				cout << "Intersection " << intersectionStart << endl;
				vIntStart.push_back(intersectionStart);
						//----//
				for (int j = i; j > muonStart; j--) { //Calculate derivative to find real signal start
					signBefore = (datamV[j]-datamV[j-1])/(j-(j-1));
			    	if (printOut) cout << "Point " << j << " signBefore " << signBefore << endl;
			    	if (signBefore >= 0) continue;
	    			else {
	    				peakStart = j;
	    				break;
	    			}
				}
				if (printOut) cout << "First real point for first peak " << peakStart << endl;
				inizioPicco.push_back(peakStart);
				continue;
			}

			else if (datamV[i] < -avgNoise && firstPeak == false && mightBePeak == false) { //Negative peak
				mightBePeak = true;
				negativePeak = true;
				firstPeak = true;
				if (printOut) cout << "First point below -5rms " << i << endl;
				//Calculate time over threshold starting point
						//----//
				m = (datamV[i]-datamV[i-1])/(i-(i-1)); //Slope of line
				q = datamV[i]-m*i;
				intersectionStart = (-avgNoise-q)/m;
				cout << "Intersection " << intersectionStart << endl;
				vIntStart.push_back(intersectionStart);
						//----//
				for (int j = i; j > muonStart; j--) { //Calculate derivative to find real signal start
					signBefore = (datamV[j]-datamV[j-1])/(j-(j-1));
			    	if (printOut) cout << "Point " << j << " signBefore " << signBefore << endl;
			    	if (signBefore >= 0) continue;
	    			else {
	    				peakStart = j;
	    				break;
	    			}
				}
				if (printOut) cout << "First real point for first peak " << peakStart << endl;
				inizioPicco.push_back(peakStart);
				continue;
			}

			//----End of first peak----//
			else if (datamV[i] < avgNoise && firstPeak == true && mightBePeak == true && positivePeak == true) { //Positive
				peakSign.push_back(1);
				positivePeak = false;
				mightBePeak = false;
				firstPeakDone = true;
				//Calculate time over threshold ending point
						//----//
				m = (datamV[i]-datamV[i-1])/(i-(i-1)); //Slope of line
				q = datamV[i]-m*i;
				intersectionEnd = (avgNoise-q)/m;
				cout << "Intersection " << intersectionEnd << endl;
				vIntEnd.push_back(intersectionEnd);
						//----//
				for (int j = i; j < datamV.size(); j++) { //calculate derivative to find real peak end
					signAfter = (datamV[j+1]-datamV[j])/(j+1-j);
					if (printOut) cout << "Point " << j << "signAfter " << signAfter << endl;
					if (signAfter <= 0) continue;
					else {
						peakEndTemp = j;
						break;
					}  
				}
				peakValues.push_back(std::distance(&datamV[0],max_element(&datamV[peakStart],&datamV[peakEndTemp])));
				if (printOut) cout << "End of the first positive peak " << peakEndTemp << endl;
				if (printOut) cout << "Positive first peak value " << *max_element(&datamV[peakStartTemp],&datamV[peakEndTemp]) << endl;
				finePicco.push_back(peakEndTemp);
				continue;
			}

			else if (datamV[i] > -avgNoise && firstPeak == true && mightBePeak == true && negativePeak == true) { //Negative
				peakSign.push_back(-1);
				negativePeak = false;
				mightBePeak = false;
				firstPeakDone = true;
				//Calculate time over threshold ending point
						//----//
				m = (datamV[i]-datamV[i-1])/(i-(i-1)); //Slope of line
				q = datamV[i]-m*i;
				intersectionEnd = (-avgNoise-q)/m;
				cout << "Intersection " << intersectionEnd << endl;
				vIntEnd.push_back(intersectionEnd);
						//----//
				for (int j = i; j < datamV.size(); j++) { //calculate derivative to find real peak end
					signAfter = (datamV[j+1]-datamV[j])/(j+1-j);
					if (printOut) cout << "Point " << j << "signAfter " << signAfter << endl;
					if (signAfter >= 0) continue;
					else {
						peakEndTemp = j;
						break;
					}  
				}
				peakValues.push_back(std::distance(&datamV[0],min_element(&datamV[peakStart],&datamV[peakEndTemp])));
				if (printOut) cout << "End of the first negative peak " << peakEndTemp << endl;
				if (printOut) cout << "Negative first peak value " << *min_element(&datamV[peakStartTemp],&datamV[peakEndTemp]) << endl;
				finePicco.push_back(peakEndTemp);
				continue;
			}
		}

		//-----Analysis of following peaks----//
		else if (firstPeakDone == true) {

			if (i > 1000) {
				if (mightBePeak == true && positivePeak == true) {//Data is still above 5 sigma
					peakEnd = 1000;
					finePicco.push_back(peakEndTemp);
					vIntEnd.push_back(peakEndTemp);
				}

				else if (mightBePeak == true && negativePeak == true) {//Data is still below -5 sigma
					peakEnd = 1000;
					finePicco.push_back(peakEndTemp);
					vIntEnd.push_back(peakEndTemp);
				}

				else if (mightBePeak == false) //Data is not above 5 sigma or below -5 sigma
					peakEnd=peakEndTemp;
			break;
			}

			else if (i >= trigStart && i <= trigEnd) {
				if (printOut) cout << "I'm in the trigger zone" << endl;
				continue;
			}

			else if (datamV[i] > avgNoise && firstPeak == true && mightBePeak == false) { //Positive
				peakStartTemp= i-1;
				mightBePeak = true;
				positivePeak = true;
				if (printOut) cout << "Positive peak start " << i << endl;
				//Calculate time over threshold starting point
						//----//
				m = (datamV[i]-datamV[i-1])/(i-(i-1)); //Slope of line
				q = datamV[i]-m*i;
				intersectionStart = (avgNoise-q)/m;
				cout << "Intersection " << intersectionStart << endl;
				vIntStart.push_back(intersectionStart);
						//----//
				inizioPicco.push_back(peakStartTemp);
				continue;
			}

			else if (datamV[i] < -avgNoise && firstPeak == true && mightBePeak == false) { //Negative
				peakStartTemp = i-1;
				mightBePeak = true;
				negativePeak = true;
				if (printOut) cout << "Negative peak start " << i << endl;
				//Calculate time over threshold starting point
						//----//
				m = (datamV[i]-datamV[i-1])/(i-(i-1)); //Slope of line
				q = datamV[i]-m*i;
				intersectionStart = (-avgNoise-q)/m;
				cout << "Intersection " << intersectionStart << endl;
				vIntStart.push_back(intersectionStart);
						//----//
				inizioPicco.push_back(peakStartTemp);
				continue;
			}
		
			//-----End of following peaks----//
			else if (datamV[i] < avgNoise && firstPeak == true && mightBePeak == true && positivePeak == true) { //End of positive peak
				peakEndTemp = i;
				if ((peakEndTemp - peakStartTemp) > 4) { //At least three points in the peak -> push back the peak
					peakValues.push_back(std::distance(&datamV[0],max_element(&datamV[peakStartTemp],&datamV[peakEndTemp])));
					peakSign.push_back(1);
					mightBePeak = false;
					positivePeak = false;
					//Calculate time over threshold ending point
						//----//
					m = (datamV[i]-datamV[i-1])/(i-(i-1)); //Slope of line
					q = datamV[i]-m*i;
					intersectionEnd = (avgNoise-q)/m;
					cout << "Intersection " << intersectionEnd << endl;
					vIntEnd.push_back(intersectionEnd);
						//----//
					if (printOut) cout << "Positive peak end " << i << endl;
					if (printOut) cout << "Positive peak value " << *max_element(&datamV[peakStartTemp],&datamV[peakEndTemp]) << endl;
					finePicco.push_back(peakEndTemp);
					continue;
				}
				else if ((peakEndTemp - peakStartTemp) <= 4) { //Not a peak -> Delete last value in peak start and peak value 
					mightBePeak = false;
					positivePeak = false;
					inizioPicco.pop_back();
					vIntStart.pop_back();
					peakEndTemp = finePicco.back();
					if (printOut) cout << "Positive peak too short -> deleted" << endl;
					continue;
				}
			}

			else if (datamV[i] > -avgNoise && firstPeak == true && mightBePeak == true && negativePeak == true) { //Negative
				peakEndTemp = i;
				if ((peakEndTemp - peakStartTemp) > 4) { //At least three points in the peak -> push back the peak
					peakValues.push_back(std::distance(&datamV[0],min_element(&datamV[peakStartTemp],&datamV[peakEndTemp])));
					peakSign.push_back(-1);
					mightBePeak = false;
					negativePeak = false;
					//Calculate time over threshold ending point
						//----//
					m = (datamV[i]-datamV[i-1])/(i-(i-1)); //Slope of line
					q = datamV[i]-m*i;
					intersectionEnd = (-avgNoise-q)/m;
					cout << "Intersection " << intersectionEnd << endl;
					vIntEnd.push_back(intersectionEnd);
						//----//
					if (printOut) cout << "Negative peak end " << i << endl;
					if (printOut) cout << "Negative peak value " << *min_element(&datamV[peakStartTemp],&datamV[peakEndTemp]) << endl;
					finePicco.push_back(peakEndTemp);
					continue;
				}
				else if ((peakEndTemp - peakStartTemp) <= 4) { //Not a peak -> Delete last value in peak start and peak value 
					mightBePeak = false;
					negativePeak = false;
					inizioPicco.pop_back();
					vIntStart.pop_back();
					//vIntEnd.pop_back();
					peakEndTemp = finePicco.back();
					if (printOut) cout << "Negative peak too short -> deleted" << endl;
				}
			}
		}
	}

	if (printOut) cout << "Done with integration interval search....start is: " << peakStart << " peak end is: " << peakEnd << endl; 
	if (printOut) cout << "There are " << peakValues.size() << " # of peaks and their values are" << endl;
	if (printOut) cout << "Inizio picco size " << inizioPicco.size() << " fine picco size " << finePicco.size() << endl;
	if (printOut) {
		for (unsigned int i = 0; i < peakValues.size(); i++) {
			cout << "Start " << inizioPicco[i] << " peak " << peakValues[i] << " end " << finePicco[i] << " and it is " << peakSign[i] << endl;
			//testMark.push_back(new TLine(peakValues[i],-20,peakValues[i],20));
		}
	}

	if(printOut) cout << "peakValues size " << peakValues.size() << endl;

	for (unsigned int i = 0; i < peakValues.size(); i++) { //compute difference between peaks to check if it's a ripple or two distinct peaks
		
		if (peakValues[i] == peakValues.back() && peakValues.size() > 1) { //last point of the vector, we must be more careful
			if(printOut) cout << "Peak number i " << i << endl;
			if(printOut) cout << "Peak value i " << peakValues[i] << endl;
			if(printOut) cout << "Sign of peak i " << peakSign[i] << endl;
			if(printOut) cout << "i-1 " << i-1 << endl;
			timeDiff = peakValues[i] - peakValues[i-1];
			timeDiffEndPeak = peakValues[i] - finePicco[i-1];
			timeDiffStartPeak = inizioPicco[i] - finePicco[i-1];
			avgPeaks = accumulate(&datamV.at(peakValues[i]-10), &datamV.at(peakValues[i]+10) , 0.0)/(20);
		}

		else if (peakValues[i] == peakValues.back() && peakValues.size() == 1) {
			if(printOut) cout << "Size of peak size 1 " << peakValues.size() << endl;
			if(printOut) cout << "Peak number i " << i << endl;
			if(printOut) cout << "Peak value i " << peakValues[i] << endl;
			if(printOut) cout << "Sign of peak i " << peakSign[i] << endl;
			break;
		}

		else {
			if(printOut) cout << "Peak number i " << i << endl;
			if(printOut) cout << "Peak value i " << peakValues[i] << endl;
			if(printOut) cout << "Sign of peak i " << peakSign[i] << endl;
			if(printOut) cout << "Peak number i+1 " << i+1 << endl;
			if(printOut) cout << "Peak value i+1 " << peakValues[i+1] << endl;
			if(printOut) cout << "Sign of peak i+1 " << peakSign[i+1] << endl;
			timeDiff = peakValues[i+1] - peakValues[i]; //Time difference between peaks
			timeDiffEndPeak = peakValues[i+1] - finePicco[i]; //Time difference from last point above threhsold and current peak
			timeDiffStartPeak = inizioPicco[i+1] - finePicco[i];
			avgPeaks = accumulate(&datamV.at(peakValues[i+1]-10), &datamV.at(peakValues[i+1]+10) , 0.0)/(20);
		}

		if ((peakSign[i] == 1 && peakSign[i+1] == 1) || (peakSign[i] == -1 && peakSign[i+1] == -1))  {//Both are positive or negative
			if(printOut) cout << "Time Diff " << timeDiff << endl;
			if(printOut) cout << "Time Diff End Peak " << timeDiffEndPeak << endl;
			if(printOut) cout << "Time Diff Start Peak " << timeDiffStartPeak << endl;

			if (timeDiff <= 40 || timeDiffEndPeak <= 40 || timeDiffStartPeak <= 40) { //Two same sign peaks with time difference <= 40 ns
				if (printOut) cout << "Two same sign peaks with time difference <= 40 ns" << endl;
				continue;
			}

			else { //Two same sign peaks with time difference > 40 ns
				if (avgPeaks >= 9*avgValue || avgPeaks <= -9*avgValue) {
					if (printOut) cout << "Peak that is > 40 ns and above noise " << peakValues[i+1] << endl;
					if (printOut) cout << "Time difference > 40 ns but average above noise, hence not really peaks";

					if (peakValues[i+1] == peakValues.back()) {
						inizioPicco.pop_back();
						finePicco.pop_back();
						peakValues.pop_back();
						peakSign.pop_back();
						vIntStart.pop_back();
						vIntEnd.pop_back();
					}
					
					else {
						inizioPicco.erase(inizioPicco.begin()+i+1); //Delete start of peak
						finePicco.erase(finePicco.begin()+i+1); //Delete end of peak
						peakValues.erase(peakValues.begin()+i+1); //Delete peak value
						vIntStart.erase(vIntStart.begin()+i+1); //Delete start of ToT
						vIntEnd.erase(vIntEnd.begin()+i+1); //Delete end of ToT
						if (printOut) cout << "Peak after deletion " << peakValues[i+1] << "\t" << peakValues[i] << endl;
						i--;
					}
					if(printOut) cout << "-> peak deleted" << endl;
					continue;
				}

				else { //Time difference 
					peakNum++;
					if (printOut) cout << "Same sign peaks with time difference > 40 ns " << endl;
					continue;
				}
			}
		}

		else if ((peakSign[i] == 1 && peakSign[i+1] == -1) || (peakSign[i] == -1 && peakSign[i+1] == 1))  {//One is positive and the other is negative or vice-versa
			if(printOut) cout << "Time Diff " << timeDiff << endl;
			if(printOut) cout << "Time Diff End Peak " << timeDiffEndPeak << endl;
			if(printOut) cout << "Time Diff Start Peak " << timeDiffStartPeak << endl;

			if (timeDiff <= 40 || timeDiffEndPeak <= 40 || timeDiffStartPeak <= 40) { //Time difference < 40 ns 
				if (printOut) cout << "Two different sign peaks with time difference <= 40 ns " << endl;
				continue;
			}

			else {//Two different sign peaks with > 40 ns time difference
				if (printOut) cout << "avgPeaks " << avgPeaks << endl;

				if (avgPeaks >= 9*avgValue || avgPeaks <= -9*avgValue) {
					if (printOut) cout << "Peak that is > 40 ns and above noise " << peakValues[i+1] << endl;
					if (printOut) cout << "Time difference > 40 ns but average above noise, hence not really peaks";

					if (peakValues[i+1] == peakValues.back()) {
						inizioPicco.pop_back();
						finePicco.pop_back();
						peakValues.pop_back();
						peakSign.pop_back();
						vIntStart.pop_back();
						vIntEnd.pop_back();
					}
				
					else {
						inizioPicco.erase(inizioPicco.begin()+i+1); //Delete start of peak
						finePicco.erase(finePicco.begin()+i+1); //Delete end of peak
						peakValues.erase(peakValues.begin()+i+1); //Delete peak value
						vIntStart.erase(vIntStart.begin()+i+1); //Delete start of ToT
						vIntEnd.erase(vIntEnd.begin()+i+1); //Delete end of ToT
						if (printOut) cout << "Peak after deletion " << peakValues[i+1] << "\t" << peakValues[i] << endl;
						i--;
					}

					if(printOut) cout << "-> peak deleted" << endl;
					continue;
				}

				else {
					peakNum++;
					if (printOut) cout << "Two different sign peaks with time difference > 40 ns " << endl;
					continue;
				}
			}
		}
		//cout << "-------" << endl;
	}

	if (printOut) cout << "size init before " << inizioPicco.size() << " size end before " << finePicco.size() << endl;

	for (unsigned int i = 0; i < inizioPicco.size(); i++) {
		if (printOut) cout << "Start " << inizioPicco[i] << " picco " << peakValues[i] << " fine picco " << finePicco[i] << endl;
		if (printOut) cout << i << endl;

		if (finePicco[i] != finePicco.back()) {
			if (peakValues[i+1] - peakValues[i] <= 40) {
				finePicco.erase(finePicco.begin()+i);
				inizioPicco.erase(inizioPicco.begin()+i+1);
			}
		}
		else 
			break;
	}

	//Calculate Time Over Threshold by summing the different ToTs if a signal has ripples;
	if (printOut) cout << "size init after " << inizioPicco.size() << " size end after " << finePicco.size() << endl;
	if (printOut) cout << "size ToT after " << vIntStart.size() << " size ToT after " << vIntEnd.size() << endl;

	double ToT = 0;
	for (unsigned int i = 0; i < vIntStart.size(); i++) {
		ToT += (vIntEnd[i] - vIntStart[i]);
	}	

	cout << "Time over thrsehold " << ToT << endl;

	//cout << "Inizio picco size dopo processsamento " << inizioPicco.size() << "\t" << "fine picco size dopo processamento " << finePicco.size() << endl;
	if(peakNum > 1) cout << "Event " << event << " strip " << strip << " number of real peaks " << peakNum << endl;
	//cout << "Event " << event << " strip " << strip << " number of real peaks " << peakNum << endl;
	return make_tuple(inizioPicco,finePicco,peakNum,ToT);
}