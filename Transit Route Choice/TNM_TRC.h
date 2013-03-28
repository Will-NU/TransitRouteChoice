//////////////////////// Transit route choice model
//////////////// 
//written by Peng Cheng and Marco Nie, March 2013. 

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

struct TRCCallbackHolder
{
	TNM_TRC *cls;
	void *params;
};


//Base class for LineProps
//class TNMDER_API LineProps
class LineProps
{
public:
    LineProps()
    {
        linePointer = NULL;
        headwayMean = 1.0;
        headwayVar  = 0.0;
        expRTT      = 0.0;
        prob        = 0.0;
        id = 0;
        par.clear();
    };
    enum Distribution {Unknown, Expon, Erlang, Deterministic};
    virtual ~LineProps() {}
    // a helper function for IntFuncOfEWT to be used in gsl integration functions
	static double ExpWaitingTimeWrapper(double x, void * holder)
	{
		return x*((LineProps*)holder)->WaitingTimePDF(x);
	}
	// calculate the factorial of an integer
	// the return type has to be double, because the return value gets very large as n increases
	inline double factorial(int n)
	{
		return (n == 1 || n == 0) ? 1 : n * factorial(n - 1);
	}        
	void *linePointer;
	double headwayMean;
	double headwayVar;
	double expRTT;                       // expected remaining travel time after boarding this line
	double prob;
    int    id; //this is the position of the line in the Lines object;
    void   Print();
    vector<double> par; //parameters for the distribution. For Erlang the first is rate, the second is scale
                        //for expon, the first is rate;
                        //for determinisitc, the first is rate;
    virtual void   InitializePar()
    {
      par.clear();  
      par.push_back(1 / headwayMean);
    }
    virtual double WaitingTimePDF(double x) = 0; //empty virtual, pdf
    virtual double WaitingTimeCDF(double x) = 0;  //empty virtual function, cdf. 
    virtual Distribution GetType() = 0;
    virtual string       GetTypeStr() = 0;
    virtual bool   UpdateProb(TNM_TRC* trc, gsl_integration_workspace * w = NULL, int size = 1000);
    virtual double GetExpectedWaitingTime(gsl_integration_workspace *w = NULL, int size = 1000); //expected waiting just for this line.
};

class LinePropsExp : public LineProps
{
public:
    LinePropsExp() : LineProps() {;}
    virtual ~LinePropsExp() {;}
    virtual double WaitingTimePDF(double x);
    virtual double WaitingTimeCDF(double x);
    virtual Distribution GetType() {return Expon;}
    virtual string       GetTypeStr()  {return "Exponential";}
    virtual double GetExpectedWaitingTime(gsl_integration_workspace *w = NULL, int size = 1000)
    {return headwayMean;}
};

class LinePropsErlang : public LineProps
{
public:
    LinePropsErlang() : LineProps() {;}
    virtual ~LinePropsErlang() {;}
    virtual double WaitingTimePDF(double x);
    virtual double WaitingTimeCDF(double x);
    virtual Distribution GetType() {return Erlang;}
    virtual string       GetTypeStr()  {return "Erlang";}
    virtual void   InitializePar();

};

class LinePropsDert : public LineProps
{
public:
    LinePropsDert() : LineProps() {;}
    virtual ~LinePropsDert() {;}
    virtual double WaitingTimePDF(double x);
    virtual double WaitingTimeCDF(double x);
    virtual Distribution GetType() {return Deterministic;}
    virtual string       GetTypeStr()  {return "Deterministic";}
    virtual bool   UpdateProb(TNM_TRC* trc, gsl_integration_workspace * w = NULL, int size = 1000);
     virtual double GetExpectedWaitingTime(gsl_integration_workspace *w = NULL, int size = 1000)
    {return 0.5*headwayMean;}
};


//Base class TNM_TRC
//class  TNMDER_API TNM_TRC
class TNM_TRC
{
protected:
	bool Info;                           // True/Flase: with/without online information at stops
	vector<LineProps*> Lines;             // a vector to store line properties for each line
    vector<int>  AttractiveSet; //Attractive set, save as ids in Lines.    
	double MinExpTotalTT;  //updated minimum expected total travel time;
    double ExpWaitTime;    //update expected waiting time;
    LineProps::Distribution  DefaultDist; //default distribution type for the class 
    bool      SameDist;
    bool      NeedUpdate;
protected:
	double   ConditionalWaitingTimePDF(double x, void * params);
	// integrand of expected waiting time
	double   IntFuncOfEWT(double x, void * params);
    double   GetExpectedTravelTimeAfterBoarding();
    void     UpdateExpectedWaitingTime(); //main internal operation, calculated expected waiting time.
    double   GetExpectedTotalTravelTime(); //main internal operation, call UpdateProb and UpdateExpectedWaitingTime to calculate the total 
    void     UpdateProb(); ////main internal operation, update choice probability. 
    void     CleanLines();
    void     CheckAllSubSet(int *arr, int size, int left, int index, std::vector<int> &l);
public:
	TNM_TRC();
	virtual ~TNM_TRC(void);
    //initialize functions
    bool Initialize(LineProps::Distribution dist = LineProps::Expon, bool info = false);
    void AddLine(void *pointer, double mean, double var, double exRtt, LineProps::Distribution dist = LineProps::Unknown);
//////////////////////////////////////////////////////////////////////
    ///////////// properties//////////////////////////////////////////////
    bool     IsSameDist() {return SameDist;}
    bool     IsInfoAvailable() {return Info;}
    LineProps::Distribution GetDefaultDist() {return DefaultDist;}
    double   GetExpectedWaitingTime() {return ExpWaitTime;}	
    double   GetTotalFreq(bool attracitveonly = true);
    double   GetMinHeadway(bool attractiveonly = true);
    double   GetMinExpTotalTravelTime()
	{
		return MinExpTotalTT;
	}
	bool GetAttractiveSet(vector<void*> &linePtrs)
	{
        if(IsUpdated())
        {
            linePtrs.clear();
            for(int i = 0;i<AttractiveSet.size();i++)
            {
                linePtrs.push_back(Lines[AttractiveSet[i]]->linePointer);
            }
            return true;
        }
        else 
        {
            cout<<"\tAttractiveSet not updated yet."<<endl;
            return false;
        }
	}
    int      GetNumberOfLines() {return Lines.size();}
    int      GetSizeOfAttractiveSet() {return AttractiveSet.size();}
///////////////////////////////////////////////////////////////////////////////////
    ///Main operations.
    void     UpdateGreedy(bool rankwithlinetimeonly = true); //default is the greedy method.         
    void     UpdateEnum(); //update based on enumeration method. caution, limmited to lines < 10;
    bool     IsUpdated() {return !NeedUpdate;}
    void     Print(bool attractiveonly = true); //io functions. 

    ///////////////////utility functions, not intended to be used by regular users.
  // a helper function for ConditionalWaitingTimePDF to be used in gsl integration functions
    static double ConditionalWaitingTimePDFWrapper(double x, void * holder)
	{
		TRCCallbackHolder *h = static_cast<TRCCallbackHolder*>(holder);
		return h->cls->ConditionalWaitingTimePDF(x, h->params);
	}
// a helper function for IntFuncOfEWT to be used in gsl integration functions
	static double IntFuncOfEWTWrapper(double x, void * holder)
	{
		TRCCallbackHolder *h = static_cast<TRCCallbackHolder*>(holder);
		return h->cls->IntFuncOfEWT(x, h->params);
	}
};

