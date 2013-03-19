#pragma once

#include "stdafx.h"
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <gsl/gsl_integration.h>

using namespace std;

class TNM_TRC;

struct CCallbackHolder
{
	TNM_TRC *cls;
	void *params;
};

struct LineProps
{
	void *linePointer;
	double headwayMean;
	double headwayVar;
	double expRTT;                       // expected remaining travel time after boarding this line
	double prob;
};

struct SortByERTT
{
	inline bool operator() (const LineProps& line1, const LineProps& line2)
	{
		return (line1.expRTT < line2.expRTT);
	}
};

class TNM_TRC
{
private:
	string DistType;                     // headway distribution type
	bool Info;                           // True/Flase: with/without online information at stops
	int NumOfLines;                      // number of lines in the choice set
	vector<LineProps> Lines;             // a vector to store line properties for each line
    vector<double> HeadwayMean;          // a list to store the mean headway of each line
	vector<double> HeadwayVar;           // a list to store the headway variance of each line
	vector<double> ExpectedTT;           // a list to store the expected travel time of each line
    vector<LineProps> AttractiveSet;    
    vector<double> Prob;
	double MinExpTotalTT; 

	enum Distribution {Error, Expon, Erlang, Deterministic};
	map<string, Distribution> DistTypes;

	void MapDist()
	{
		DistTypes["expon"] = Expon;
		DistTypes["erlang"] = Erlang;
		DistTypes["deterministic"] = Deterministic;
	}

	vector<double> ErlangRate;           // parameter for erlang distribution
	vector<int> ErlangShape;             // parameter for erlang distribution

	vector<double> ExponRate;            // parameter for exponential distribution

	vector<double> DeterminRate;         // parameter for deterministic distribution

	// PDF of waiting time for a specific line
	double WaitingTimePDF(double x, int lineIndex);

	// CDF of waiting time for a specific line
	double WaitingTimeCDF(double x, int lineIndex);

	// PDF of waiting time at the stop conditional to boarding a carrier of line i
	// params store the line identifier i
	double ConditionalWaitingTimePDF(double x, void * params);
    // a helper function for ConditionalWaitingTimePDF to be used in gsl integration functions
	static double ConditionalWaitingTimePDFWrapper(double x, void * holder)
	{
		CCallbackHolder *h = static_cast<CCallbackHolder*>(holder);
		return h->cls->ConditionalWaitingTimePDF(x, h->params);
	}

	// integrand of expected waiting time
	double IntFuncOfEWT(double x, void * params);
	// a helper function for IntFuncOfEWT to be used in gsl integration functions
	static double IntFuncOfEWTWrapper(double x, void * holder)
	{
		CCallbackHolder *h = static_cast<CCallbackHolder*>(holder);
		return h->cls->IntFuncOfEWT(x, h->params);
	}

	// calculate the factorial of an integer
	// the return type has to be double, because the return value gets very large as n increases
	double factorial(int n)
	{
		return (n == 1 || n == 0) ? 1 : n * factorial(n - 1);
	}

	void AddOneLine(LineProps line);
public:
	TNM_TRC();
	~TNM_TRC(void);
	void Initialize(vector<LineProps> lines, string distType = "Expon", bool info = false);

	vector<double> GetProb();
	double GetExpectedWaitingTime();
	double GetExpectedTravelTimeAfterBoarding();
	double GetExpectedTotalTravelTime();
	void UpdateAttractiveSet();

	double GetMinExpTotalTravelTime()
	{
		return MinExpTotalTT;
	}

	vector<LineProps> GetAttractiveSet()
	{
		return AttractiveSet;
	}
};

