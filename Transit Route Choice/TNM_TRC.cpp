#include "StdAfx.h"
#include "TNM_TRC.h"

using namespace std;

void LineProps::Print()
{
    cout<<setw(12)<<GetTypeStr()
        <<setw(12)<<this->headwayMean<<setw(12)<<this->headwayVar
        <<setw(12)<<this->expRTT<<setw(12)<<this->prob<<" ";
    for(int i = 0;i<par.size();i++)
        cout<<par[i]<<",";
    cout<<endl;
}

double LineProps::GetExpectedWaitingTime(gsl_integration_workspace *w, int size)
{
    bool neww = false;
    if(w == NULL)
    {
        neww = true;
        if(size > 0) w = gsl_integration_workspace_alloc (size);
        else         
        {
            w = gsl_integration_workspace_alloc (1000);
            size = 1000;
        }
    }
    else
    {
        if(size <=0) return false;
    }
    gsl_function F;
    F.function = &LineProps::ExpWaitingTimeWrapper;
    F.params = this;    
    double error, expwait;
    gsl_integration_qagiu(&F, 0, 0, 1e-6, size, w, &expwait, &error);
   // gsl_integration_qag(&F, 0, headwayMean, 0, 1e-6, size, 6, w, &exptwait, &error);
    if(neww)
    {
        gsl_integration_workspace_free (w);
    }
    return true;
}
bool LineProps::UpdateProb(TNM_TRC *trc, gsl_integration_workspace *w, int size)
{

    bool neww = false;
    if(w == NULL)
    {
        neww = true;
        if(size > 0) w = gsl_integration_workspace_alloc (size);
        else         
        {
            w = gsl_integration_workspace_alloc (1000);
            size = 1000;
        }
    }
    else
    {
        if(size <=0) return false;
    }
    gsl_function F;
	F.function = &TNM_TRC::ConditionalWaitingTimePDFWrapper;
    TRCCallbackHolder holder;
	holder.cls = trc;
	holder.params = &id;
    F.params = &holder;    
    double error;
    gsl_integration_qagiu(&F, 0, 0, 1e-6, size, w, &prob, &error);
    if(neww)
    {
        gsl_integration_workspace_free (w);
    }
    return true;
    
}
double LinePropsExp::WaitingTimePDF(double x)
{
    return par[0] * exp(- par[0] * x);
}
double LinePropsExp::WaitingTimeCDF(double x)
{
    return 1 - exp(- par[0] * x);
}


double LinePropsDert::WaitingTimePDF(double x)
{
    if (x >= 0 && x <= headwayMean) return par[0];//HeadwayMean[lineIndex])
    else return 0.0;
}

double LinePropsDert::WaitingTimeCDF(double x)
{
    if (x >= 0 && x <= headwayMean) return par[0] *x;
    else return 1.0;
}

bool LinePropsDert::UpdateProb(TNM_TRC *trc, gsl_integration_workspace *w, int size)
{
    bool neww = false;
    if(w == NULL)
    {
        neww = true;
        if(size > 0) w = gsl_integration_workspace_alloc (size);
        else         
        {
            w = gsl_integration_workspace_alloc (1000);
            size = 1000;
        }
    }
    else
    {
        if(size <=0) return false;
    }
    gsl_function F;
	F.function = &TNM_TRC::ConditionalWaitingTimePDFWrapper;
    TRCCallbackHolder holder;
	holder.cls = trc;
	holder.params = &id;
    F.params = &holder;    
    double error;
    gsl_integration_qag(&F, 0, headwayMean, 0, 1e-6, size, 6, w, &prob, &error);
    if(neww)
    {
        gsl_integration_workspace_free (w);
    }
    return true;
    
}
void LinePropsErlang::InitializePar()
{
    par.clear();
    par.push_back(1 / headwayMean);
	par.push_back(int(headwayMean * headwayMean / headwayVar + 0.5));
}
double LinePropsErlang::WaitingTimePDF(double x)
{
    double f = 0.0;
	double lambda = par[0];
	double k      = par[1];
	for (int i = 0; i < k; i++)
	{
		f += pow(k * lambda * x, i) / factorial(i);
	}
	return  f*lambda * exp(- k * lambda * x); 
    
}
double LinePropsErlang::WaitingTimeCDF(double x)
{
    double f = 0.0;
    double lambda = par[0];
	double k      = par[1];
	for (int i = 0; i < k; i++)
	{
		f += (1.0 - i / k) * pow(k * lambda * x, i) / factorial(i);
	}
	f *= exp(- k * lambda * x);
	f = 1 - f;
    return f;
}

TNM_TRC::TNM_TRC()
{
    SameDist = true;
    NeedUpdate = true;
    DefaultDist = LineProps::Unknown;
}

TNM_TRC::~TNM_TRC(void)
{
    CleanLines();
}

void TNM_TRC::CleanLines()
{
    for(int i = 0;i<Lines.size();i++)
    {
        delete Lines[i];
    }
    Lines.clear();
    AttractiveSet.clear();
}

void TNM_TRC::AddLine(void *pointer, double mean, double var, double exRtt, LineProps::Distribution dist)
{
    LineProps::Distribution pdist;
    if(dist != LineProps::Unknown) pdist = dist;
    else                pdist = DefaultDist;
    if(pdist != DefaultDist) SameDist = false;
    LineProps *lp;
    switch(pdist)
    {
    case LineProps::Expon:
        lp = new LinePropsExp;
        break;
    case LineProps::Erlang:
        lp = new LinePropsErlang;        
        break;
    case LineProps::Deterministic:
        lp = new LinePropsDert;
        break;
    }
    lp->linePointer = pointer;
    lp->headwayMean = mean;
    lp->headwayVar  = var;
    lp->expRTT = exRtt;
    lp->prob   = 0.0;
    lp->id     = Lines.size();
    lp->InitializePar(); //call this function untill all parameters are set. 
    NeedUpdate = true;
    Lines.push_back(lp);
}
bool TNM_TRC::Initialize(LineProps::Distribution dist, bool info) 
{
    if(dist == LineProps::Unknown) 
    {
        cout<<"\tUnsuppoted distribution type"<<endl;
        return false;
    }
    else   DefaultDist = dist;
    Info = info;
    CleanLines();
    SameDist      = true;
    NeedUpdate    = true;
    ExpWaitTime   = -1;
    MinExpTotalTT = -1;
    return true;
}

double TNM_TRC::ConditionalWaitingTimePDF(double x, void *params)
{
	int lineIndex = *(int *) params;
	double f = 1;    
    for (vector<int>::iterator pv = AttractiveSet.begin(); pv!=AttractiveSet.end();pv++) 
	{
        if (*pv != lineIndex)
            f *= (1 - Lines[*pv]->WaitingTimeCDF(x));
		else
            f *= Lines[*pv]->WaitingTimePDF(x);
	}

	return f;
}

double TNM_TRC::IntFuncOfEWT(double x, void *params)
{
	double f = 0;	
    for (vector<int>::iterator  pv = AttractiveSet.begin(); pv!=AttractiveSet.end();pv++) 
	{
        int id = *pv;
        f += ConditionalWaitingTimePDF(x, &id);
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


void TNM_TRC::UpdateProb()
{
	//vector<double> pVector;                     // a vector to store the probability of takin each line

	// declare some variables for the numerical integration
	
    
    if(SameDist && DefaultDist == LineProps::Expon)//in this case, no integration is needed. 
    {	                  
        double sumRate = this->GetTotalFreq();	
        for (vector<int>::iterator pv = AttractiveSet.begin(); pv!=AttractiveSet.end();pv++) 
		{
            LineProps *prop = Lines[*pv];
            prop->prob = prop->par[0] / sumRate;
		}
    }
    else
    {
        gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);      
        for (vector<int>::iterator pv = AttractiveSet.begin(); pv!=AttractiveSet.end();pv++) 
	    {
            Lines[*pv]->UpdateProb(this, w, 1000);
        }
       	gsl_integration_workspace_free (w);

    }

	//Prob = pVector;
	//return pVector;
}

double  TNM_TRC::GetTotalFreq(bool attractiveonly )
{
    double sumRate = 0.0;
    if(attractiveonly)
    {       
        for (vector<int>::iterator pv = AttractiveSet.begin(); pv!=AttractiveSet.end();pv++) 
		{
            sumRate += 1/Lines[*pv]->headwayMean;//ExponRate[i];
		}
    }
    else
    {
        for (vector<LineProps*>::iterator pv = Lines.begin(); pv!=Lines.end();pv++) 
		{
            sumRate += 1/(*pv)->headwayMean;//ExponRate[i];
		}
    }
    return sumRate;
}

double TNM_TRC::GetMinHeadway(bool attractiveonly)
{
    double minHeadway = 1e20;
    if(attractiveonly)
    {        
        for (vector<int>::iterator pv = AttractiveSet.begin(); pv!=AttractiveSet.end();pv++) 
		{
            LineProps *prop = Lines[*pv];
            if(minHeadway > prop->headwayMean)//ExponRate[i];
            minHeadway = prop->headwayMean;
		}
    }
    else
    {
        for (vector<LineProps*>::iterator pv = Lines.begin(); pv!=Lines.end();pv++) 
		{           
            if(minHeadway > (*pv)->headwayMean)//ExponRate[i];
            minHeadway = (*pv)->headwayMean;
		}
    }
    return minHeadway;
}
void TNM_TRC::UpdateExpectedWaitingTime()
{
    
    if(SameDist && DefaultDist == LineProps::Expon)//in this case, no integration is needed. 
    {	  
       ExpWaitTime = 1/GetTotalFreq();       
    }
        // only integrate to the minimum headway among all lines for deterministic headway
	// because the waiting time cannot be greater than the minimum headway
	// in fact, integrating to any value that is greater than the minimum headway is fine,
	// the result should not be affected
    else if(SameDist && DefaultDist == LineProps::Deterministic)
    {
        gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
	    double result, error;
        gsl_function F;
	    F.function = &TNM_TRC::IntFuncOfEWTWrapper;
        TRCCallbackHolder holder;
	    holder.cls = this;
	    F.params = &holder;
        double minHeadway = GetMinHeadway();
        gsl_integration_qag(&F, 0, minHeadway, 0, 1e-6, 1000, 6, w, &ExpWaitTime, &error);
        gsl_integration_workspace_free (w);
    }
    else //all other cases numerical integration from 0 to infinity for all other types of distribution
	// note that this also works for exponentail distribution headway and deterministic headway
    {
        gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
	    double result, error;
        gsl_function F;
	    F.function = &TNM_TRC::IntFuncOfEWTWrapper;
        TRCCallbackHolder holder;
	    holder.cls = this;
	    F.params = &holder;
        gsl_integration_qagiu(&F, 0, 0, 1e-6, 1000, w, &ExpWaitTime, &error);
        gsl_integration_workspace_free (w);
    }

	
}

double TNM_TRC::GetExpectedTravelTimeAfterBoarding()
{
	double expTTAB = 0;  
    
    for (vector<int>::iterator pv = AttractiveSet.begin(); pv!=AttractiveSet.end();pv++) 	
    {
        LineProps *prop  = Lines[*pv];
        expTTAB += prop->expRTT * prop->prob;//pVector[i];
	}
    return expTTAB;
}

double TNM_TRC::GetExpectedTotalTravelTime()
{
    UpdateProb();
    UpdateExpectedWaitingTime();
	return GetExpectedWaitingTime() + GetExpectedTravelTimeAfterBoarding(); 
}


//implementation of greedy method.
void TNM_TRC::UpdateGreedy(bool rankwithlttonly)
{
    if(IsUpdated()) 
    {
        AttractiveSet.clear();
        NeedUpdate = true;
    }
    multimap<double, LineProps*, less<double> > plines;
    //rank Lines.
    if(rankwithlttonly) //rank based on line travel time only
    {
        for(int i = 0;i<Lines.size();i++)
        {
            LineProps* prop = Lines[i];          
            plines.insert(pair<double, LineProps*>(prop->expRTT, prop));
        }
    }
    else //rank based on the sum of line travel time and expected line waiting time. 
    {
        gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
        for(int i = 0;i<Lines.size();i++)
        {
            LineProps* prop = Lines[i];
            double wt = prop->GetExpectedWaitingTime(w, 1000);
            plines.insert(pair<double, LineProps*>(prop->expRTT + wt, prop));
        }
        gsl_integration_workspace_free (w);
    }
    //double curTT = 1e20;
    if(DefaultDist == LineProps::Expon && IsSameDist()) //special case for 
    {
        multimap<double, LineProps*, less<double> >::iterator pv= plines.begin();       
        do 
        {
            AttractiveSet.push_back(pv->second->id);
            MinExpTotalTT = GetExpectedTotalTravelTime();
            pv++;
        }while(pv!=plines.end() && pv->second->expRTT < MinExpTotalTT); //if line travel time > expectation total tt, then add this line will gurantee increase the travel time. 
    }
    else
    {
        double preMinETT;
        MinExpTotalTT= 1e20;
        multimap<double, LineProps*, less<double> >::iterator pv= plines.begin();
        do
        {
            preMinETT = MinExpTotalTT;
          //  cout<<"current total TT "<<preMinETT<<endl;
            AttractiveSet.push_back(pv->second->id);
            MinExpTotalTT = GetExpectedTotalTravelTime();
            pv++;
        }while(preMinETT > MinExpTotalTT && pv!=plines.end());
       // cout<<"\tafter loop, expeted total = "<<MinExpTotalTT <<endl;
        if(preMinETT < MinExpTotalTT)
        {
            AttractiveSet.pop_back();
            MinExpTotalTT = GetExpectedTotalTravelTime();//take the last one.
        }
    }
    NeedUpdate = false;
}

void TNM_TRC::CheckAllSubSet(int *arr, int size, int left, int index, std::vector<int> &l)
{
    if(left==0){
        std::vector<int> pa = AttractiveSet;
        AttractiveSet = l;
        double ett = GetExpectedTotalTravelTime();
        if(ett < MinExpTotalTT)
        {
            MinExpTotalTT = ett;
        }
        else
        {
            AttractiveSet = pa; ///otherwise reset Attractive set;
            MinExpTotalTT = GetExpectedTotalTravelTime();
        }
        return;
    }
    for(int i=index; i<size;i++){
        l.push_back(arr[i]);
        CheckAllSubSet(arr,size,left-1,i+1,l);
        l.pop_back();
    }
}
void TNM_TRC::UpdateEnum()
{
    if(Lines.size() > 10) 
    {
        cout<<"\tEnumeration method is not allowed becasue the number of lines > 10"<<endl;
        return; 
    }
    if(IsUpdated()) 
    {
        AttractiveSet.clear();
        NeedUpdate = true;
    }
    int nl = Lines.size();
    int *parray = new int[nl];//{1,2,3,4,5};
    for(int i = 0;i<nl;i++) parray[i] = i;
    std::vector<int> lt;   
    MinExpTotalTT = 1e20;
    for(int i = 1;i<=nl;i++)
    {
        CheckAllSubSet(parray,nl,i,0,lt); //check all subsets of size i. 
    }
    NeedUpdate = false;
}
void TNM_TRC::Print(bool attractive)
{
   
    if(attractive)
    {
        cout<<"Size of attractive set: "<<setw(4)<<AttractiveSet.size()<<" Online informaiton "<<(Info?"Yes":"No")<<endl;
        if(IsUpdated())
        { 
            cout<<"Minimum expected total travel time = "<<setw(12)<<MinExpTotalTT
                <<", of which expected waiting time = "<<setw(12)<<ExpWaitTime<<endl;
        cout<<setw(12)<<"DistType"
        <<setw(12)<<"MeanHeadway"<<setw(12)<<"VarHeadway"<<setw(12)<<"RunTime"
        <<setw(12)<<"ChoiceProb"<<" Parameters"<<endl;
            for(int i = 0 ;i<AttractiveSet.size();i++)
            {
                Lines[AttractiveSet[i]]->Print();
            }
        }
        else
        {
            cout<<"\tAttractive set is not updated."<<endl;
        }
    }
    else
    {
         cout<<"Number of lines = "<<Lines.size()<<endl;
         cout<<setw(12)<<"DistType"
        <<setw(12)<<"MeanHeadway"<<setw(12)<<"VarHeadway"<<setw(12)<<"RunTime"
        <<setw(12)<<"ChoiceProb"<<" Parameters"<<endl;
            for(int i = 0 ;i<Lines.size();i++)
            {
                Lines[i]->Print();
            }
        
    }
}
