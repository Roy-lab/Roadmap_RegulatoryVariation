#include <iostream>
#include <fstream>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "Framework.H"


Framework::Framework()
{
}

Framework::~Framework()
{
}

//Don't change
int 
Framework::readInteractions(const char* epFName)
{
	ifstream inFile(epFName);
	//cout<<epFName<<endl;
	char buffer[1024];
	int total=0;
	int mapped=0;
	while(inFile.good())
	{
		inFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}	
		char* tok=strtok(buffer,"\t");
		string chrom1;
		int begin1=0;
		int end1=0;
		string chrom2;
		int begin2=0;
		int end2=0;
		int tokCnt=0;
		string gene;
		//each part of the line
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				getChromInfo(tok,chrom1,begin1,end1);
			}
			else if(tokCnt==1)
			{
				getChromInfo(tok,chrom2,begin2,end2);
			}
			else if(tokCnt==2)
			{
				gene.append(tok);
				//cout << gene << endl;
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
		if(gene.length()==0)
		{
			continue;
		}
		if(strcmp(chrom1.c_str(),chrom2.c_str())!=0)	
		{
			continue;
		}
		//pair is two regions and a name
		Pair* p=new Pair;
		//cout<<gene<<endl;
		//non-promoter is always the first bin
		p->nonpromoter.begin=begin1;
		p->nonpromoter.end=end1;
		char ekeychr[1024]; //represents the name of the non-promoter
		sprintf(ekeychr,"%s_%d_%d",chrom1.c_str(),begin1,end1);
		p->nonpromoter.name.append(ekeychr);
		//do the same for the non-promoter.
		p->promoter.begin=begin2;
		p->promoter.end=end2;
		p->gene.append(gene.c_str());
		p->chrom.append(chrom1.c_str());
		char pkeychr[1024];
		sprintf(pkeychr,"%s_%d_%d",chrom1.c_str(),begin2,end2);
		p->promoter.name.append(pkeychr);
		char keychr[1024];
		sprintf(keychr,"%s_%d_%d-%s_%d_%d",chrom1.c_str(),begin1,end1,chrom2.c_str(),begin2,end2);
		string key(keychr);
		interactions[keychr]=p;
		cout <<p->gene<<endl;
		total++;
		// num enhancer per gene?? 
	}
	inFile.close();
	return 0;
}

//dont change
int
Framework::getChromInfo(char* tok,string& chr, int& start, int& end)
{
	char* begin=tok;
	int tokCnt=0;
	while(begin!=NULL)
	{
		char* pos=strchr(begin,'_');
		if(pos!=NULL)
		{
			*pos='\0';
		}
		if(tokCnt==0)
		{
			chr.append(begin);
		}
		else if(tokCnt==1)
		{
			start=atoi(begin);
		}
		else if(tokCnt==2)
		{
			end=atoi(begin);
		}
		if(pos!=NULL)
		{
			begin=pos+1;
		}
		else
		{
			begin=NULL;
		}
		tokCnt++;
	}
	return 0;
}

int
Framework::readPIQMotifInstances(const char* aFName)
{
	ifstream inFile(aFName);
	char buffer[1024];
	while(inFile.good())
	{
		inFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{		
			continue;
		}
		if(strchr(buffer,'#')!=NULL)
		{
			continue;
		}
		//reading in the file line by line
		string chrName;
		string motifName;
		string strand;
		double pval;
		int position;
		int tokCnt=0;
		char* tok=strtok(buffer,"\t");
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				chrName.append(tok);
			}
			else if (tokCnt==1)
			{
				position=atof(tok);
			}
			else if(tokCnt==2)
			{
				pval=atof(tok);
			}
			else if(tokCnt==3)
			{
				strand.append(tok);
			}
			else if(tokCnt==4)
			{
				motifName.append(tok);
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
		//if pval >> 0.9 then I want the chrName, position and peak name stored in a structure
		if (pval > 0.9)
		{
			Motif* m=new Motif;
			m->chrom=chrName;
			m->pos=position;
			m->motifID=motifName;
			char keychr[1024];
			sprintf(keychr,"%s_%d",chrName.c_str(),position);
			string key(keychr);
			sigMotifs[keychr]=m;
		}
	}
	inFile.close();
	return 0;
}

int
Framework::generateDistalNetwork(const char* aFName)
{
	ofstream oFile(aFName);
	//interating over the set of all interactions
	for(map<string,Pair*>::iterator pIter=interactions.begin();pIter!=interactions.end();pIter++)
	{
		Pair* p=pIter->second;
		//cout << p->chrom <<endl;
		for(map<string,Motif*>::iterator aIter=sigMotifs.begin();aIter!=sigMotifs.end();aIter++)
		{
			Motif* m=aIter->second;
			if (p->chrom==m->chrom)
			{
				//cout << m->chrom <<endl;
				string motifName;
				if(getOverlap(&p->nonpromoter,m,motifName))
				{
					oFile << m->motifID <<"\t"<< p->gene << endl;
				}
				else
				{
					continue;
				}
			}
			else
			{
				continue;
			}
		}
	}
	oFile.close();
	return 0;
}


bool
Framework::getOverlap(Region* r,  Motif* m, string& nameme)
{
	bool found=false;

	int b1,e1,pos;
	b1=r->begin;
	pos=m->pos;
	e1=r->end;

	if(pos>b1 & pos<e1)
	{
		nameme.clear();
		nameme.append(m->motifID);
		found=true;
	}
	return found;
}

int
main(int argc, const char** argv)
{
	Framework fw;
	fw.readInteractions(argv[1]);
	fw.readPIQMotifInstances(argv[2]);
	fw.generateDistalNetwork(argv[3]);
	return 0;
}
