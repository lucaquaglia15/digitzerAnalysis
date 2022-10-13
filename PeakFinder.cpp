// PeakFinder.cpp
#include "PeakFinder.h"

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <numeric>

void diff(std::vector<float> in, std::vector<float>& out)
{
	out = std::vector<float>(in.size()-1);

	for(unsigned int i=1; i<in.size(); ++i)
		out[i-1] = in[i] - in[i-1];
}

void vectorElementsProduct(std::vector<float> a, std::vector<float> b, std::vector<float>& out)
{
	out = std::vector<float>(a.size());

	for(unsigned int i=0; i<a.size(); ++i)
		out[i] = a[i] * b[i];
}

void findIndicesLessThan(std::vector<float> in, float threshold, std::vector<int>& indices)
{
	for(unsigned int i=0; i<in.size(); ++i)
		if(in[i]<threshold)
			indices.push_back(i+1);
}

void selectElementsFromIndices(std::vector<float> in, std::vector<int> indices, std::vector<float>& out)
{
	for(unsigned int i=0; i<indices.size(); ++i)
		out.push_back(in[indices[i]]);
}

void selectElementsFromIndices(std::vector<int> in, std::vector<int> indices, std::vector<int>& out)
{
	for(unsigned int i=0; i<indices.size(); ++i)
		out.push_back(in[indices[i]]);
}

void signVector(std::vector<float> in, std::vector<int>& out)
{
	out = std::vector<int>(in.size());

	for(unsigned int i=0; i<in.size(); ++i)
	{
		if(in[i]>0)
			out[i]=1;
		else if(in[i]<0)
			out[i]=-1;
		else
			out[i]=0;
	}
}

void scalarProduct(float scalar, std::vector<float> in, std::vector<float>& out)
{
	out = std::vector<float>(in.size());

	for(unsigned int i=0; i<in.size(); ++i)
		out[i] = scalar * in[i];
}

/*void PeakFinder::findPeaks(std::vector<float> x0, std::vector<int>& peakInds, bool includeEndpoints, float extrema)
{
	int minIdx = distance(x0.begin(), min_element(x0.begin(), x0.end()));
	int maxIdx = distance(x0.begin(), max_element(x0.begin(), x0.end()));

	float sel = (x0[maxIdx]-x0[minIdx])/4.0;
	int len0 = x0.size();

	scalarProduct(extrema, x0, x0);

	std::vector<float> dx;
	diff(x0, dx);
	replace(dx.begin(), dx.end(), 0.0f, -PeakFinder::EPS);
	std::vector<float> dx0(dx.begin(), dx.end()-1);
	std::vector<float> dx0_1(dx.begin()+1, dx.end());
	std::vector<float> dx0_2;

	vectorElementsProduct(dx0, dx0_1, dx0_2);

	std::vector<int> ind;
	findIndicesLessThan(dx0_2, 0, ind); // Find where the derivative changes sign	
	std::vector<float> x;
	float leftMin;
	int minMagIdx;
	float minMag;
	
	if(includeEndpoints)
	{
		//x = [x0(1);x0(ind);x0(end)];	
		selectElementsFromIndices(x0, ind, x);		
		x.insert(x.begin(), x0[0]);
		x.insert(x.end(), x0[x0.size()-1]);
		//ind = [1;ind;len0];
		ind.insert(ind.begin(), 1);
		ind.insert(ind.end(), len0);
		minMagIdx = distance(x.begin(), std::min_element(x.begin(), x.end()));
		minMag = x[minMagIdx];		
		//std::cout<<"Hola"<<std::endl;
		leftMin = minMag;
	}
	else
	{
		selectElementsFromIndices(x0, ind, x);
		if(x.size()>2)
		{
			minMagIdx = distance(x.begin(), std::min_element(x.begin(), x.end()));		
			minMag = x[minMagIdx];				
			leftMin = x[0]<x0[0]?x[0]:x0[0];
		}
	}

	int len = x.size();

	if(len>2)
	{
		float tempMag = minMag;
    	bool foundPeak = false;
    	int ii;

		if(includeEndpoints)
		{
    		// Deal with first point a little differently since tacked it on
        	// Calculate the sign of the derivative since we tacked the first
        	//  point on it does not neccessarily alternate like the rest.
    		std::vector<float> xSub0(x.begin(), x.begin()+3);//tener cuidado subvector
    		std::vector<float> xDiff;//tener cuidado subvector
    		diff(xSub0, xDiff);

    		std::vector<int> signDx;
    		signVector(xDiff, signDx);

        	if (signDx[0] <= 0) // The first point is larger or equal to the second
        	{
				if (signDx[0] == signDx[1]) // Want alternating signs
				{
					x.erase(x.begin()+1);
					ind.erase(ind.begin()+1);
					len = len-1;
				}
        	}
        	else // First point is smaller than the second
        	{
				if (signDx[0] == signDx[1]) // Want alternating signs
				{
					x.erase(x.begin());
					ind.erase(ind.begin());
					len = len-1;
				}
        	}
		}

		//Skip the first point if it is smaller so we always start on maxima
		if ( x[0] >= x[1] )
			ii = 0;
		else
			ii = 1;

		//Preallocate max number of maxima
		float maxPeaks = ceil((float)len/2.0);
		std::vector<int> peakLoc(maxPeaks,0);
		std::vector<float> peakMag(maxPeaks,0.0);
		int cInd = 1;
		int tempLoc;		
    
    	while(ii < len)
    	{
        	ii = ii+1;//This is a peak
        	//Reset peak finding if we had a peak and the next peak is bigger
        	//than the last or the left min was small enough to reset.
        	if(foundPeak)
        	{
            	tempMag = minMag;
            	foundPeak = false;
            }
        
        	//Found new peak that was lager than temp mag and selectivity larger
        	//than the minimum to its left.
        
        	if( x[ii-1] > tempMag && x[ii-1] > leftMin + sel )
        	{
            	tempLoc = ii-1;
            	tempMag = x[ii-1];
        	}

        	//Make sure we don't iterate past the length of our vector
        	if(ii == len)
            	break; //We assign the last point differently out of the loop

        	ii = ii+1; // Move onto the valley
        	
        	//Come down at least sel from peak
        	if(!foundPeak && tempMag > sel + x[ii-1])
            {            	
	            foundPeak = true; //We have found a peak
	            leftMin = x[ii-1];
	            peakLoc[cInd-1] = tempLoc; // Add peak to index
	            peakMag[cInd-1] = tempMag;
	            cInd = cInd+1;
	        }
        	else if(x[ii-1] < leftMin) // New left minima
            	leftMin = x[ii-1];
            
        }

		// Check end point
		if(includeEndpoints)
		{        
			if ( x[x.size()-1] > tempMag && x[x.size()-1] > leftMin + sel )
			{
				peakLoc[cInd-1] = len-1;
				peakMag[cInd-1] = x[x.size()-1];
				cInd = cInd + 1;
			}
			else if( !foundPeak && tempMag > minMag )// Check if we still need to add the last point
			{
				peakLoc[cInd-1] = tempLoc;
				peakMag[cInd-1] = tempMag;
				cInd = cInd + 1;
			}
		}
		else if(!foundPeak)
		{
			float minAux = x0[x0.size()-1]<x[x.size()-1]?x0[x0.size()-1]:x[x.size()-1];
			if ( x[x.size()-1] > tempMag && x[x.size()-1] > leftMin + sel )
			{
				peakLoc[cInd-1] = len-1;
				peakMag[cInd-1] = x[x.size()-1];
				cInd = cInd + 1;
			}
			else if( !tempMag >  minAux + sel)// Check if we still need to add the last point
			{
				peakLoc[cInd-1] = tempLoc;
				peakMag[cInd-1] = tempMag;
				cInd = cInd + 1;
			}
		}

		//Create output
    	if( cInd > 0 )
    	{        	
        	std::vector<int> peakLocTmp(peakLoc.begin(), peakLoc.begin()+cInd-1);
			selectElementsFromIndices(ind, peakLocTmp, peakInds);        	
        }		

	}
	//else
	//{
		//input signal length <= 2
	//}
}*/

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