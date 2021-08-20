#include "CorrectedPval.H"
#include <vector>
#include <algorithm>
#include <iostream>

CorrectedPval::CorrectedPval()
{
}

CorrectedPval::~CorrectedPval()
{
}

int
CorrectedPval::estimateQvalues(map<string,double>& pvalue_map,map<string,double>& qvalue_map)
{
	vector<Pval_T> pvals;
	int cnt = 0;
	for (map<string,double>::iterator piter=pvalue_map.begin();piter!=pvalue_map.end();piter++)
	{
		Pval_T p;
		p.name = piter->first;
		p.pval = piter->second;
		p.qval = 0;
		p.index = cnt;
		cnt ++;
		pvals.push_back(p);
	}
	sort(pvals.begin(), pvals.end(), customLess());
	double m=(double)pvals.size();
	for(int k=0;k<pvals.size();k++)
	{
		double cpval=(pvals[k].pval*m)/((double)(k+1));
		pvals[k].qval = cpval;
	}
	double minfdr=pvals[pvals.size()-1].qval;
	for(int l=pvals.size()-1;l>=0;l--)
	{
		double q = pvals[l].qval;
		if(q>minfdr)
		{
			pvals[l].qval=minfdr;
		}
		else 
		{
			minfdr=q;
		}
	}
	for(int k=0;k<pvals.size();k++)
	{
		qvalue_map[pvals[k].name] = pvals[k].qval;
	}
	/*
	for (int k=0;k<pvals.size();k++)
	{
		cout << pvals[k].index << "\t" << pvals[k].name << "\t" << pvals[k].pval << "\t" << pvals[k].qval << endl;
	}
	*/
	return 0;
}
