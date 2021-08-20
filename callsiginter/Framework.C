#include <iostream>
#include <fstream>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <sys/timeb.h>
#include <sys/time.h>
#include "gsl/gsl_randist.h"
#include "gsl/gsl_rng.h"
#include "gsl/gsl_cdf.h"
#include "gsl/gsl_statistics_double.h"
#include "HyperGeomPval.H"
#include "CorrectedPval.H"
#include "Framework.H"

Framework::Framework()
{
	totalCounts=0;
}

Framework::~Framework()
{
}

int
Framework::readPairFile(const char* aFName, int binSize, int radius,int scale)
{
	ifstream inFile(aFName);
	char buffer[1024];
	int totalPairs=0;
	if(scale<1)
	{
		scale=1;
	}
	while(inFile.good())
	{
		inFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		if(strchr(buffer,'-')==NULL && strchr(buffer,'\t')==NULL)
		{
			continue;
		}
		char* tok=strchr(buffer,'-');
		if(tok!=NULL)
		{
			*tok='\0';
		}
		char* e=buffer;
		char* p=tok+1;
		char* tok2=strchr(p,'\t');
		double paircnt=0;
		if(tok2!=NULL)
		{
			*tok2='\0';
			char* cntval=tok2+1;
			char* end=strchr(cntval,'\t');
			if(end!=NULL)
			{
				*end='\0';
			}
			paircnt=atof(cntval);
		}
		string eKey(e);
		if(regionSet.find(eKey)==regionSet.end())
		{
			Framework::Region* eregion=new Framework::Region;
			regionSet[eKey]=eregion;	
			string chrom;	
			int start;
			int end;
			tok=strtok(e,"_");
			int tokCnt=0;
			while(tok!=NULL)
			{
				if(tokCnt==0)
				{
					eregion->chromosome.append(tok);
				}
				else if(tokCnt==1)
				{	
					eregion->begin=atoi(tok);
				}
				else if(tokCnt==2)
				{
					eregion->end=atoi(tok);
				}
				tok=strtok(NULL,"_");
				tokCnt++;
			}
		}
		string pKey(p);
		if(regionSet.find(pKey)==regionSet.end())
		{
			Framework::Region* pregion=new Framework::Region;
			regionSet[pKey]=pregion;	
			tok=strtok(p,"_");
			int tokCnt=0;
			while(tok!=NULL)
			{
				if(tokCnt==0)
				{
					pregion->chromosome.append(tok);
				}
				else if(tokCnt==1)
				{	
					pregion->begin=atoi(tok);
				}
				else if(tokCnt==2)
				{
					pregion->end=atoi(tok);
				}
				tok=strtok(NULL,"_");
				tokCnt++;
			}
			regionSet[pKey]=pregion;
		}
		/*else
		{
			Region* pregion=regionSet[pKey];
			pregion->cnt=pregion->cnt+paircnt;	
		}*/
		int dist=getDistance(eKey,pKey);
		//cout << "dist=" << dist << endl;
		if(dist>radius)
		{
			continue;
		}
		if(dist<20000)
		{
			//continue;
		}
		Framework::Pair* pair=new Framework::Pair;
		pair->e.append(eKey.c_str());
		pair->p.append(pKey.c_str());
		pair->cnt=paircnt;
		//This is a hack to make binomial distribution give more significant pvalues..
		//pair->cnt=paircnt*scale;
		pair->cnt=paircnt+log(scale);		
		pair->dist=dist;
		int distBin=dist/binSize;
		vector<Pair*>* pairSet=NULL;
		if(binnedPairs.find(distBin)==binnedPairs.end())
		{
			pairSet=new vector<Pair*>;
			binnedPairs[distBin]=pairSet;
		}
		else
		{
			pairSet=binnedPairs[distBin];
		}
		pairSet->push_back(pair);
		totalPairs++;
		totalCounts=totalCounts+paircnt;
		if((totalPairs%1000000)==0)
		{
			cout <<"Found " << totalPairs << " so far" <<endl;
		}
	}
	inFile.close();
	//cout <<"Found " << totalPairs << " pairs at dist " << radius<< endl;
	return 0;
}

int
Framework::getSignificantInteractions(const char* aFName)
{
	ofstream oFile(aFName);
	map<string,double> allPvals;
	map<string,double> allPvals_HGP;
	map<string,double> allPvalsG;
	map<int,vector<string>*> indexSet;
	int iter=0;
	for(map<int,vector<Pair*>*>::iterator dIter=binnedPairs.begin();dIter!=binnedPairs.end();dIter++)	
	{
		cout <<"Calling significant interactions for bin " << dIter->first << endl;
		int countForDist=0;
		double countForDistL=0;
		vector <Pair*>* pairSet=dIter->second;
		cout << "Number of interactions is " << pairSet->size() << endl;
		vector<string>* indices=new vector<string>;
		indexSet[dIter->first]=indices;
		double binP=1/((double)pairSet->size());
		double* dataL = new double[pairSet->size()];
		for(int i=0;i<pairSet->size();i++)
		{
			Pair* p=(*pairSet)[i];	
			countForDist=countForDist+(int)round(exp(p->cnt)-1);
			countForDistL+=p->cnt;
			dataL[i]=p->cnt;
		}
		double aveL=(double)countForDistL/(double)pairSet->size();
		double std=gsl_stats_sd_m(dataL,1,pairSet->size(),aveL);		
		delete [] dataL;
		double binP_Fithic=(((double)countForDistL/(double)pairSet->size())/totalCounts);
		cout << "Ave=" << aveL << "\t" << "std=" << std << endl;
		//continue;
		//double binP_Fithic=(((double)countForDist)/totalCounts);
		//Now go over all the pair entries and get a binomial pvalue. WE can do the hypergeom later
		for(int i=0;i<pairSet->size();i++)
		{
			Pair* p=(*pairSet)[i];
			int k=round(exp(p->cnt)-1);
			double kL=p->cnt;
			/*Region* r1=regionSet[p->e];
			Region* r2=regionSet[p->p];
			if(r1->cnt.find(dIter->first)==r1->cnt.end())
			{
				cout <<"Logical error! No cnts at distance bin " << dIter->first << " for region"  << r1->chromosome <<"_" << r1->begin<<"_" << r1->end<< endl;
				exit(0);
			}
			if(r2->cnt.find(dIter->first)==r2->cnt.end())
			{
				cout <<"Logical error! No cnts at distance bin " << dIter->first << " for region" << r2->chromosome<<"_"<< r2->begin<<"_"<< r2->end << endl;
				exit(0);
			}
			int t1=r1->cnt[dIter->first];
			int t2=r2->cnt[dIter->first];*/
			double binom_P=gsl_cdf_binomial_P(k, binP_Fithic, countForDist);
			double binom_P_Duan=gsl_cdf_binomial_P(k, binP, countForDist);
			double gaussian_P=1-gsl_cdf_gaussian_P(kL-aveL,std);
			binom_P=1-binom_P;
			binom_P_Duan=1-binom_P_Duan;
			//HyperGeomPval hgp;
			//double hyperGeomP=hgp.getOverRepPval(k,t1,t2,countForDist-t2);
			char key[256];
			sprintf(key,"%d",iter);
			string keys(key);
			indices->push_back(keys);
			//allPvals[keys]=binom_P;
			allPvals[keys]=binom_P_Duan;
			allPvalsG[keys]=gaussian_P;
			//allPvals_HGP[keys]=hyperGeomP;
			iter++;
		}
	}
	CorrectedPval cp;
	map<string,double> allQvals;
	cp.estimateQvalues(allPvals,allQvals);
	map<string,double> allQvalsG;
	cp.estimateQvalues(allPvalsG,allQvalsG);
	//map<string,double> allQvals_HGP;
	//cp.estimateQvalues(allPvals_HGP,allQvals_HGP);
	//oFile <<"Pair\tCnt\tBinomQ\tBinomP\tHGPQ\tHGPP"<<endl;
	oFile <<"Pair\tCnt\tDist\tBinomQ\tBinomP\tGaussianQ\tGaussianP"<<endl;
	for(map<int,vector<string>*>::iterator kIter=indexSet.begin();kIter!=indexSet.end();kIter++)
	{
		vector<string>* indices=kIter->second;
		vector<Pair*>* pset=binnedPairs[kIter->first];
		for(int i=0;i<indices->size();i++)
		{
			Pair* p=(*pset)[i];
			oFile << p->e << "-" << p->p<< "\t" << p->cnt << "\t" << p->dist << "\t" << allQvals[(*indices)[i]] << "\t" << allPvals[(*indices)[i]] << "\t" << allQvalsG[(*indices)[i]] << "\t" << allPvalsG[(*indices)[i]]
									//<< "\t" << allQvals_HGP[(*indices)[i]] <<"\t"<< allPvals_HGP[(*indices)[i]] 
									<< endl; 
		}
	}
	oFile.close();
	return 0;
}

int
Framework::getDistance(string& ekey, string& pkey)
{
	Region* e=regionSet[ekey];
	Region* p=regionSet[pkey];
	int dist=0;
	int first1=e->begin;
	int first2=e->end;
	int second1=p->begin;
	int second2=p->end;
	if(first1>second1)
	{
		first1=p->begin;
		first2=p->end;
		second1=e->begin;
		second2=e->end;
	}
	if(second1<first2) 
	{
		//there is overlap
		dist=0;
	}
	else
	{
		dist=second1-first2;
	}
	return dist;
}

int
main(int argc, const char** argv)
{
	if(argc<6)
	{	
		cout <<"Usage: getSignificantIntr inputfile distbin radius scale outputfile"<<endl;
		return 0;
	}
	Framework fw;
	fw.readPairFile(argv[1],atoi(argv[2]),atoi(argv[3]),atoi(argv[4]));
	fw.getSignificantInteractions(argv[5]);
	
	return 0;
}
