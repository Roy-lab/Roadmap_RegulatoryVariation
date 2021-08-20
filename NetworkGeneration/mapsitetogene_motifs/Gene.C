#include <stdio.h>
#include "Gene.H"

Gene::Gene()
{
	upstreamNeighbor=NULL;
	downstreamNeighbor=NULL;
}

Gene::~Gene()
{
}

int 
Gene::setName(const char* anid)
{
	name.append(anid);
	return 0;
}

int 
Gene::setGeneCoordinate(int begin, int end)
{
	geneStart=begin;
	geneEnd=end;
	return 0;
}

int 
Gene::set5utr(int begin, int end)
{
	char keybuff[256];
	sprintf(keybuff,"%d-%d",begin,end);
	string key(keybuff);
	if(fiveUtrSet.find(key)==fiveUtrSet.end())
	{
		Coordinate* coord=new Coordinate;
		coord->begin=begin;
		coord->end=end;
		fiveUtrSet[key]=coord;
	}
	return 0;
}

int 
Gene::setCDS(int begin,int end)
{
	char keybuff[256];
	sprintf(keybuff,"%d-%d",begin,end);
	string key(keybuff);
	if(cdsSet.find(key)==cdsSet.end())
	{
		Coordinate* coord=new Coordinate;
		coord->begin=begin;
		coord->end=end;
		cdsSet[key]=coord;
	}
	return 0;
}

int 
Gene::setStrand(char astrand)
{
	strand=astrand;
	return 0;
}

int
Gene::setChromosome(const char* achrom)
{
	chromosome.append(achrom);
	return 0;
}

char
Gene::getStrand()
{
	return strand;
}

const char*
Gene::getChromosome()
{
	return chromosome.c_str();
}

int 
Gene::inCDS(int start, int end, int margin)
{
	int cdsCnt=0;
	return cdsCnt;
}	

int 
Gene::inGene(int start,int end, int margin)
{
	int lowerbound=geneStart-margin;
	int upperbound=geneStart+margin;
	if(start>=lowerbound && end<=upperbound)
	{
		return 1;
	}
	return 0;
}

int 
Gene::in5utr(int start,int end, int margin)
{
	int hits=0;
	for(map<string,Coordinate*>::iterator aIter=fiveUtrSet.begin();aIter!=fiveUtrSet.end();aIter++)
	{
		int utrbegin=aIter->second->begin;
		int utrend=aIter->second->end;
		int lowerbound=utrbegin-margin;
		int upperbound=utrbegin+margin;
		if(start>=lowerbound && end<=upperbound)
		{
			hits++;
		}
	}
	return hits;
}

int
Gene::get5utrCnts()
{
	int utrCnt=fiveUtrSet.size();
	return utrCnt;
}
