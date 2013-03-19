#include "StdAfx.h"
#include "TNM_TRC.h"

using namespace std;

TNM_TRC::TNM_TRC()
{
}

TNM_TRC::~TNM_TRC(void)
{
}

void TNM_TRC::Initialize(vector<LineProps> lines, string distType, bool info) 
{
	DistType = distType;
	transform(DistType.begin(), DistType.end(), DistType.begin(), ::tolower);
    Info = info;
    Lines = lines;

	MapDist();
	
	//// Derive the parameters of the particular distribution from the mean and variance
	//switch (DistTypes[DistType])
	//{
	//case Expon :
	//	for (int i = 0; i < NumOfLines; i++)
	//	{
	//		ExponRate.push_back(1 / HeadwayMean[i]);
	//	}
	//	break;
	//case Erlang :
	//	for (int i = 0; i < NumOfLines; i++)
	//	{
	//		ErlangRate.push_back(1 / HeadwayMean[i]);
	//		ErlangShape.push_back(int(HeadwayMean[i] * HeadwayMean[i] / HeadwayVar[i] + 0.5));
	//	}
	//	break;
	//case Deterministic :
	//	for (int i = 0; i < NumOfLines; i++)
	//	{
	//		DeterminRate.push_back(1 / HeadwayMean[i]);
	//	}
	//	break;
	//default :
	//	printf("Error: the distribution type is not supported.\n");
	//}
}

double TNM_TRC::WaitingTimePDF(double x, int lineIndex)
{
	double f = 0;
	double lambda;
	double k;

	switch (DistTypes[DistType])
	{
	case Expon :
		lambda = ExponRate[lineIndex];
		f = lambda * exp(- lambda * x);
		break;
	case Erlang :
		lambda = ErlangRate[lineIndex];
		k = ErlangShape[lineIndex];
		for (int i = 0; i < k; i++)
		{
			f += pow(k * lambda * x, i) / factorial(i);
		}
		f *= lambda * exp(- k * lambda * x); 
		break;
	case Deterministic :
		if (x >= 0 && x <= HeadwayMean[lineIndex])
			f = DeterminRate[lineIndex];
		else
			f = 0;
		break;
	default :
		printf("Error: the distribution type is not supported.\n");
	}

	return f;
}

double TNM_TRC::WaitingTimeCDF(double x, int lineIndex)
{
	double f = 0;
	double lambda;
	double k;

	switch (DistTypes[DistType])
	{
	case Expon :
		f = 1 - exp(- ExponRate[lineIndex] * x);
		break;
	case Erlang :
		lambda = ErlangRate[lineIndex];
		k = ErlangShape[lineIndex];
		for (int i = 0; i < k; i++)
		{
			f += (1.0 - i / k) * pow(k * lambda * x, i) / factorial(i);
		}
		f *= exp(- k * lambda * x);
		f = 1 - f;
		break;
	case Deterministic :
		if (x >= 0 && x <= HeadwayMean[lineIndex])
			f = DeterminRate[lineIndex] * x;
		else
			f = 1;
		break;
	default :
		printf("Error: the distribution type is not supported.\n");
	}

	return f;
}

double TNM_TRC::ConditionalWaitingTimePDF(double x, void *params)
{
	int lineIndex = *(int *) params;
	double f = 1;

	for (int i = 0; i < NumOfLines; i++)
	{
		if (i != lineIndex)
			f *= (1 - WaitingTimeCDF(x, i));
		else
			f *= WaitingTimePDF(x, i);
	}

	return f;
}

double TNM_TRC::IntFuncOfEWT(double x, void *params)
{
	double f = 0;

	for (int i = 0; i < NumOfLines; i++)
	{
		f += ConditionalWaitingTimePDF(x, &i);
	}
	f *= x;

	return f;
}

// Notes: to get the route choice probability, the expected waiting time and other related stuff,
// we could use the numerical integration to integrate the corresponding integrand from 0 to infinity
// for all types of distributions. For special types of distributions, for example, exponential distribution,
// we can directly use the analytical form, and for deterministic case, we don't have to integrate the 
// functions from 0 to infinity because we know that the waiting time could not exceed a certain value. If we
// take advantage of this, it may save some computation time.


vector<double> TNM_TRC::GetProb()
{
	vector<double> pVector;                     // a vector to store the probability of takin each line

	// declare some variables for the numerical integration
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
	double result, error;
    
    // the sum of the rate (frequency) for exponential distribution 
	double sumRate = 0;                                      

	gsl_function F;
	F.function = &TNM_TRC::ConditionalWaitingTimePDFWrapper;
    CCallbackHolder holder;
	holder.cls = this;

	switch (DistTypes[DistType])
	{
	// analytical solution for exponential distribution
	case Expon :                       
		for (int i = 0; i < NumOfLines; i++)
		{
			sumRate += ExponRate[i];
		}
		for (int i = 0; i < NumOfLines; i++)
		{
			pVector.push_back(ExponRate[i] / sumRate);
		}
		break;
	// only integrate to the constant headway of each line for deterministic headway
	// because the line waiting time cannot be greater than the corresponding constant headway
	case Deterministic :
		for (int i = 0; i < NumOfLines; i++)
		{
			holder.params = &i;
			F.params = &holder;
			gsl_integration_qag(&F, 0, HeadwayMean[i], 0, 1e-6, 1000, 6, w, &result, &error);
			pVector.push_back(result);
		}
		break;
	// numerical integration from 0 to infinity for all other types of distribution
	// note that this also works for exponentail distribution headway and deterministic headway
	default :
		for (int i = 0; i < NumOfLines; i++)
		{
			holder.params = &i;
			F.params = &holder;
			gsl_integration_qagiu(&F, 0, 0, 1e-6, 1000, w, &result, &error);
			pVector.push_back(result);
		}
	}

	gsl_integration_workspace_free (w);

	Prob = pVector;
	return pVector;
}

double TNM_TRC::GetExpectedWaitingTime()
{
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
	double result, error;
	double sumRate = 0;
	double minHeadway = HeadwayMean[0];

	gsl_function F;
	F.function = &TNM_TRC::IntFuncOfEWTWrapper;
    CCallbackHolder holder;
	holder.cls = this;
	F.params = &holder;

	switch (DistTypes[DistType])
	{
	// analytical solution for exponential distribution
	case Expon :
		for (int i = 0; i < NumOfLines; i++)
		{
			sumRate += ExponRate[i];
		}
		result = 1 / sumRate;
		break;
    // only integrate to the minimum headway among all lines for deterministic headway
	// because the waiting time cannot be greater than the minimum headway
	// in fact, integrating to any value that is greater than the minimum headway is fine,
	// the result should not be affected
	case Deterministic :
		for (int i = 0; i < NumOfLines; i++)
		{
			if (HeadwayMean[i] < minHeadway)
				minHeadway = HeadwayMean[i];
		}
		gsl_integration_qag(&F, 0, minHeadway, 0, 1e-6, 1000, 6, w, &result, &error);
		break;
	// numerical integration from 0 to infinity for all other types of distribution
	// note that this also works for exponentail distribution headway and deterministic headway
	default :
		gsl_integration_qagiu(&F, 0, 0, 1e-6, 1000, w, &result, &error);
	}

	gsl_integration_workspace_free (w);

	return result;
}

double TNM_TRC::GetExpectedTravelTimeAfterBoarding()
{
	double expTTAB = 0;
	vector<double> pVector;

	pVector = GetProb();
	for (int i = 0; i < NumOfLines; i++)
	{
		expTTAB += ExpectedTT[i] * pVector[i];
	}

	return expTTAB;
}

double TNM_TRC::GetExpectedTotalTravelTime()
{
	return GetExpectedWaitingTime() + GetExpectedTravelTimeAfterBoarding(); 
}

void TNM_TRC::AddOneLine(LineProps line)
{
	double mean;
	double var;
	double expRTT;
	mean = line.headwayMean;
	var = line.headwayVar;
	expRTT = line.expRTT;

	HeadwayMean.push_back(mean);
	HeadwayVar.push_back(var);
	ExpectedTT.push_back(expRTT);

	// Update the distribution parameters
	switch (DistTypes[DistType])
	{
	case Expon :
		ExponRate.push_back(1 / mean);
		break;
	case Erlang :
		ErlangRate.push_back(1 / mean);
		ErlangShape.push_back(int(mean * mean / var + 0.5));
		break;
	case Deterministic :
		DeterminRate.push_back(1 / mean);
		break;
	default :
		printf("Error: the distribution type is not supported.\n");
	}
	
	// Update NumOfLines
	NumOfLines++;
}

void TNM_TRC::UpdateAttractiveSet()
{
	double expTotalTT;

	// sort the lines by the expected remaining travel time
	sort(Lines.begin(), Lines.end(), SortByERTT()); 
    
	NumOfLines = 0;
	AddOneLine(Lines[0]);
    expTotalTT = GetExpectedTotalTravelTime();
    MinExpTotalTT = expTotalTT;
	AttractiveSet.push_back(Lines[0]);
	AttractiveSet[0].prob = Prob[0];

	for (int i = 1; i < Lines.size(); i++)
	{
		AddOneLine(Lines[i]);
		expTotalTT = GetExpectedTotalTravelTime();
		if (expTotalTT < MinExpTotalTT)
		{
			MinExpTotalTT = expTotalTT;
			AttractiveSet.push_back(Lines[i]);
			// update the probability
			for (int j = 0; j < i+1; j++)
			{
				AttractiveSet[j].prob = Prob[j];
			}
		}
		else
		{
			break;
		}
	}
}