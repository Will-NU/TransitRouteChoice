// Transit Route Choice.cpp : Defines the entry point for the console application.
//
#include "stdafx.h"
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <string>
#include <gsl/gsl_integration.h>
#include "TNM_TRC.h"

using namespace std;
////////////////////////////////////////////////////////////


int main(int argc, char **argv)
{
    printf("=========== test ===========\n");
	string distType; 
	distType = "Expon";            // could be "Erlang", "Expon", or "Deterministic", case insensitive
	void *p;
	p = NULL;
	LineProps line1 = {p, 20, 400, 30, 0};
	LineProps line2 = {p, 15, 225, 40, 0};
	LineProps line3 = {p, 10, 100, 45, 0};
	vector<LineProps> lines;
	lines.push_back(line1);
	lines.push_back(line2);
	lines.push_back(line3);

	vector<LineProps> AttrSet;
	double minETTT;

	TNM_TRC Obj;
	Obj.Initialize(lines);
	Obj.UpdateAttractiveSet();
	AttrSet = Obj.GetAttractiveSet();
	minETTT = Obj.GetMinExpTotalTravelTime();

	cout << "Number of lines in the available line set: " << lines.size() << endl;
	cout << "Number of lines in the attractive set: " << AttrSet.size() << endl;

	for (int i = 0; i < AttrSet.size(); i++)
	{
		cout << "The probability of taking line " << i + 1 << ": " << AttrSet[i].prob << endl;
	}
	//vector<double> hMean;
	//hMean.push_back(30.0);
	//hMean.push_back(50.0);
	//hMean.push_back(5.0);
	//vector<double> hVar;
	//hVar.push_back(400.0 / 2);
	//hVar.push_back(225.0 / 2);
	//hVar.push_back(25.0 / 2);
	//vector<double> eTT;
	//eTT.push_back(27.0);
	//eTT.push_back(38.0);
	//eTT.push_back(40.0);
	//vector<double> pVec;
	//TNM_TRC lines(distType, hMean, hVar, eTT);
	//pVec = lines.GetProb();

	//cout << "The headway distribution is : " << distType << endl;
	//for (int i = 0; i < pVec.size(); i++)
	//{
	//	cout << "Probability of taking line " << i+1 << ": " << pVec[i] << endl;
	//}

	//double expWT = lines.GetExpectedWaitingTime();
	//double expTTAB = lines.GetExpectedTravelTimeAfterBoarding();

	//cout << "expectedWaitingTime: " << expWT << endl;
	//cout << "expectedTravelTimeAfterBoarding: " << expTTAB << endl;
	//cout << "expectedTotalTravelTime: " << expWT + expTTAB << endl;
	getchar();
	return 0;
}


