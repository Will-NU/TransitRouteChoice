// Transit Route Choice.cpp : Defines the entry point for the console application.
//
#include "StdAfx.h"
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <string>
#include "TNM_TRC.h"
using namespace std;
////////////////////////////////////////////////////////////


void RunTRC(TNM_TRC* trc)
{
    //print all transit lines;
     trc->Print(false);
    ////Main operation, calculate choice probabilit, waiting time and minimum expected time.
    //cout<<"\tFind attractive set using the greedy method, based on line travel time only.."<<endl;
    //trc->UpdateGreedy();
    ////print attractive transit lines;
    //trc->Print();
    //cout<<"\tFind attractive set using the greedy method, based on total line travel time.."<<endl;
    //trc->UpdateGreedy(false);
    ////print attractive transit lines;
    //trc->Print();
    //cout<<"\tFind attractive set using the enumeration method, based on total line travel time.."<<endl;
    //trc->UpdateEnum();
    ////print attractive transit lines;
    //trc->Print();
    trc->UpdateAttractiveSet();
	trc->Print();
}
int main(int argc, char **argv)
{  
	void *p;
	p = NULL;
    //declear TRC object;
	TNM_TRC Obj;
    ////initialization;
    //cout<<"\tBasic exponential distribution."<<endl; //Table 6 in Gentile et al. 2005. TS.
    //Obj.Initialize(LineProps::Expon);
    //Obj.AddLine(p, 30, 900, 27);
    //Obj.AddLine(p, 50, 50*50, 38);
    //Obj.AddLine(p, 5, 25, 40);
    //RunTRC(&Obj);
  
    //
    //cout<<"\n\tAnother exponential distribution."<<endl; //this example shows that ranking based on total line travel time is also problematic. 
    //Obj.Initialize(LineProps::Expon);
    //Obj.AddLine(p, 20, 400, 30);
    //Obj.AddLine(p, 15, 225, 40);
    //Obj.AddLine(p, 10, 100, 44.3);
    //RunTRC(&Obj);

    cout<<"\n\tErlang distribution."<<endl; // partial information
    Obj.Initialize(LineProps::Erlang);
    Obj.AddLine(p, 20, 400/4, 30, true);
    Obj.AddLine(p, 15, 225/3, 40, true);
    Obj.AddLine(p, 10, 100/2, 45);
	Obj.AddLine(p, 8, 64/2, 35);
    RunTRC(&Obj);

    //cout<<"\n\tDeterministic distribution."<<endl; //Table 7, in Gentile et al. 2005, TC
    //Obj.Initialize(LineProps::Deterministic);
    //Obj.AddLine(p, 30, 0, 27, true);
    //Obj.AddLine(p, 50, 0, 38, true);
    //Obj.AddLine(p, 5,  0, 40, true);
    //RunTRC(&Obj);

    cout<<"\n\tExponential distribution."<<endl; // no information 
    Obj.Initialize(LineProps::Expon);
    Obj.AddLine(p, 20, 400, 30);
    Obj.AddLine(p, 15, 225, 40);
    Obj.AddLine(p, 10, 100, 44.3);
	Obj.AddLine(p, 8, 64, 20);
    RunTRC(&Obj);

    cout<<"\n\tExponential distribution."<<endl; // partial information case 1
    Obj.Initialize(LineProps::Expon);
    Obj.AddLine(p, 20, 400, 30, true);
    Obj.AddLine(p, 15, 225, 40, true);
    Obj.AddLine(p, 10, 100, 44.3);
	Obj.AddLine(p, 8, 64, 20);
    RunTRC(&Obj);

	cout<<"\n\tExponential distribution."<<endl; // partial information case 2
    Obj.Initialize(LineProps::Expon);
    Obj.AddLine(p, 20, 400, 30, true);
    Obj.AddLine(p, 15, 225, 40);
    Obj.AddLine(p, 10, 100, 44.3);
	Obj.AddLine(p, 8, 64, 20, true);
    RunTRC(&Obj);

	cout<<"\n\tExponential distribution."<<endl; // complete information 
    Obj.Initialize(LineProps::Expon);
    Obj.AddLine(p, 20, 400, 30, true);
    Obj.AddLine(p, 15, 225, 40, true);
    Obj.AddLine(p, 10, 100, 44.3, true);
	Obj.AddLine(p, 8, 64, 20, true);
    RunTRC(&Obj);

 //   int n = 2;
	//int m = 3;
	//double **x;
	//x = new double*[n];
	//for (int i = 0; i < n; i++)
	//	x[i] = new double[m];
	//x[1][2] = 4;
	//cout << x[1][2] << endl;
	//delete [] x;
	getchar();
	return 0;
}


