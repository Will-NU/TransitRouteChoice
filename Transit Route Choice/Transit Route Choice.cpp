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
    //Main operation, calculate choice probabilit, waiting time and minimum expected time.
    cout<<"\tFind attractive set using the greedy method, based on line travel time only.."<<endl;
    trc->UpdateGreedy();
    //print attractive transit lines;
    trc->Print();
    cout<<"\tFind attractive set using the greedy method, based on total line travel time.."<<endl;
    trc->UpdateGreedy(false);
    //print attractive transit lines;
    trc->Print();
    cout<<"\tFind attractive set using the enumeration method, based on total line travel time.."<<endl;
    trc->UpdateEnum();
    //print attractive transit lines;
    trc->Print();
}
int main(int argc, char **argv)
{  
	void *p;
	p = NULL;
    //declear TRC object;
	TNM_TRC Obj;
    //initialization;
    cout<<"\tBasic exponential distribution."<<endl; //Table 6 in Gentile et al. 2005. TS.
    Obj.Initialize(LineProps::Expon);
    Obj.AddLine(p, 30, 900, 27);
    Obj.AddLine(p, 50, 50*50, 38);
    Obj.AddLine(p, 5, 25, 40);
    RunTRC(&Obj);
  
    
    cout<<"\n\tAnother exponential distribution."<<endl; //this example shows that ranking based on total line travel time is also problematic. 
    Obj.Initialize(LineProps::Expon);
    Obj.AddLine(p, 20, 400, 30);
    Obj.AddLine(p, 15, 225, 40);
    Obj.AddLine(p, 10, 100, 44.3);
    RunTRC(&Obj);

    cout<<"\n\tErlang distribution."<<endl; 
    Obj.Initialize(LineProps::Erlang);
    Obj.AddLine(p, 20, 400/4, 30);
    Obj.AddLine(p, 15, 225/3, 40);
    Obj.AddLine(p, 10, 100/2, 45);
    RunTRC(&Obj);

    cout<<"\n\tDeterministic distribution."<<endl; //Table 7, in Gentile et al. 2005, TC
    Obj.Initialize(LineProps::Deterministic);
    Obj.AddLine(p, 30, 0, 27);
    Obj.AddLine(p, 50, 0, 38);
    Obj.AddLine(p, 5,  0, 40);
    RunTRC(&Obj);
	getchar();
	return 0;
}


