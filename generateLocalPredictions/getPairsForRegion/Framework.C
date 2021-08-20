#include <iostream>
#include <fstream>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <sys/timeb.h>
#include <sys/time.h>
#include "gsl/gsl_randist.h"
#include "Framework.H"

Framework::Framework()
{
}

Framework::~Framework()
{
}


int
Framework::setMaxDist(int mDist)
{
	maxDist=mDist;
	return 0;
}
	

//Read the enhancer-promoter interactions,  use the count, put the features in a block
int 
Framework::readPairs(const char* aFName, int start, int end, int countField)
{
	ifstream inFile(aFName);
	char buffer[1024];
	int totalPairs=0;
	int lineCnt=0;
	int regionID=0;
	while(inFile.good())
	{
		inFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		if(lineCnt==0)
		{
			char* tok=strtok(buffer,"\t");
			int tokCnt=0;
			while(tok!=NULL)
			{
				string cName;
				cName.append(tok);
				colNames.push_back(cName);
				tok=strtok(NULL,"\t");	
				tokCnt++;
			}
			lineCnt++;
			continue;
		}
		if(strchr(buffer,'-')==NULL)
		{
			continue;
		}
		string pairKey;
		char e[1024];
		char p[1024];
		double true_paircnt=0;
		double pred_paircnt=0;
		char* tok=strtok(buffer,"\t");
		int tokCnt=0;
		vector<double>* featVect=new vector<double>;
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				char* pos=strchr(tok,'-');
				if(pos==NULL)
				{
					cout <<"Bad format, expecting - "<<endl;
					exit(0);
				}
				*pos='\0';
				strcpy(e,tok);
				strcpy(p,pos+1);
			}
			else
			{
				featVect->push_back(atof(tok));	
				if(tokCnt==countField)
				{
					true_paircnt=atof(tok);
				}
			}
			/*else if(tokCnt==3)
			{
				pred_paircnt=atof(tok);
			}*/
			tok=strtok(NULL,"\t");
			tokCnt++;
		}

		//featVect->push_back(true_paircnt);
		string eKey(e);
		int enhID=0;
		int proID=0;
		Framework::Region* eregion=NULL;
		if(regionSet.find(eKey)==regionSet.end())
		{
			eregion=new Framework::Region;
			regionID++;
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
					int rID=eregion->begin/5000;
					if(maxRegionID<rID)
					{
						maxRegionID=rID;
					}
				}
				else if(tokCnt==2)
				{
					eregion->end=atoi(tok);
				}
				tok=strtok(NULL,"_");
				tokCnt++;
			}
			enhID=eregion->begin/5000;
		}
		else
		{
			eregion=regionSet[eKey];
			enhID=eregion->begin/5000;
		}
		string pKey(p);
		Framework::Region* pregion=NULL;
		if(regionSet.find(pKey)==regionSet.end())
		{
			pregion=new Framework::Region;
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
					int rID=pregion->begin/5000;
					if(maxRegionID<rID)
					{
						maxRegionID=rID;
					}
				}
				else if(tokCnt==2)
				{
					pregion->end=atoi(tok);
				}
				tok=strtok(NULL,"_");
				tokCnt++;
			}
			regionSet[pKey]=pregion;
			proID=pregion->begin/5000;
		}
		else
		{
			pregion=regionSet[pKey];
			proID=pregion->begin/5000;
		}
		
		//checks if 1 mb
		int dist=getDistance(eKey,pKey);
		//if(dist>maxDist)
		//{
		//	continue;
		//}
		//Now do the selection of the regions in the window
		
		//this can no longer be the case
		
		if(eregion->begin>=start && eregion->end<=end || pregion->begin>=start && pregion->end<=end)
		{

		pairKey.append(eKey);
		pairKey.append("-");
		pairKey.append(pKey);
		allPairs[pairKey]=totalPairs;

		//found th regions here
		if(allPairsIDToName.find(totalPairs)!=allPairsIDToName.end())
		{
			cout <<"ID " << totalPairs << " already associated with " << allPairsIDToName[totalPairs] << endl;
			exit(0);
		}


		allPairsIDToName[totalPairs]=pairKey;
		allPairs_Features[pairKey]=featVect;
		totalPairs++;
		PairIDs* pid=new PairIDs;
		pid->e=enhID;
		pid->p=proID;
		allPairs_RegionPairs[pairKey]=pid;
		map<string,double>* promSet_True=NULL;
		if(pairSet_True.find(eKey)==pairSet_True.end())
		{
			promSet_True=new map<string,double>;
			pairSet_True[eKey]=promSet_True;
		}	
		else
		{
			promSet_True=pairSet_True[eKey];
		}

		(*promSet_True)[pKey]=(true_paircnt);
		//We will make this bi-directional just so to keep the logic simple
		/*map<string,double>* enhSet_True=NULL;
		map<string,double>* enhSet_Pred=NULL;
		if(pairSet_True.find(pKey)==pairSet_True.end())
		{
			enhSet_True=new map<string,double>;
			pairSet_True[pKey]=enhSet_True;
			enhSet_Pred=new map<string,double>;
			pairSet_Pred[pKey]=enhSet_Pred;
		}
		else
		{
			enhSet_True=pairSet_True[pKey];
			enhSet_Pred=pairSet_Pred[pKey];
		}
		(*enhSet_True)[eKey]=true_paircnt;
		(*enhSet_Pred)[eKey]=pred_paircnt;*/
	}

	}
	inFile.close();
	//cout <<"Found " << totalPairs << " pairs at dist " << maxDist << endl;
	return 0;
}

int
Framework::setRegionNeighborSetSize(int size)
{
	regionDistNeighborSize=size;
	return 0;
}

//Now we add neighbor region pairs
/*int
Framework::createNeighborGraph()
{
	int pairID=0;
	for(map<string,map<string,double>*>::iterator rIter=pairSet_True.begin();rIter!=pairSet_True.end();rIter++)
	{
		
		Region* row=regionSet[rIter->first];
		int rowID=row->begin/5000;
		map<string,double>* colSet=rIter->second;
		for(map<string,double>::iterator cIter=colSet->begin();cIter!=colSet->end();cIter++)
		{
			Region* col=regionSet[cIter->first];
			int colID=col->begin/5000;
			int rowbegin=rowID-regionDistNeighborSize;
			if(rowbegin<0)
			{
				rowbegin=0;
			}
			int rowend=rowID+regionDistNeighborSize;
			if(rowend>maxRegionID)
			{
				rowend=maxRegionID;
			}
			int colbegin=colID-regionDistNeighborSize;
			if(colbegin<0)
			{
				colbegin=0;
			}
			int colend=colID+regionDistNeighborSize;
			if(colend>maxRegionID)
			{
				colend=maxRegionID;
			}
			//myname
			string thisnode;
			char buffer[1024];
			sprintf(buffer,"%s-%s",rIter->first.c_str(),cIter->first.c_str());
			thisnode.append(buffer);
			int thisnodeID=allPairs[thisnode];
			map<int,int>* neighbors=NULL;
			if(regionPairGraph.find(thisnodeID)==regionPairGraph.end())
			{
				neighbors=new map<int,int>;
				regionPairGraph[thisnodeID]=neighbors;
			}
			else
			{
				neighbors=regionPairGraph[thisnodeID];
			}
			for(int r=rowbegin;r<rowend;r++)
			{
				for(int c=colbegin;c<colend;c++)
				{
					char buffer[1024];
					sprintf(buffer,"%s_%d_%d-%s_%d_%d",row->chromosome.c_str(),r*5000,(r+1)*5000,col->chromosome.c_str(),c*5000,(c+1)*5000);
					string buff(buffer);
					if(allPairs.find(buff)==allPairs.end())
					{
						//skippity skip
						continue;
					}
					int othernodeID=allPairs[buff];
					if(thisnodeID==othernodeID)
					{
						continue;
					}
					//Eventually we can also compute the similarity of region pairs based on the chromatin signals.
					(*neighbors)[othernodeID]=0;
				}
			}
		}
	}
	return 0;
}*/

int
Framework::generateNeighborGraph(const char* outSuff)
{
	//createNeighborGraph();
	char fName[1024];
	//sprintf(fName,"%s/graph.txt",outSuff);
	//ofstream oFile(fName);
	
	sprintf(fName,"%s/featurefile.txt",outSuff);
	ofstream pFile(fName);
	
	for(int i=0;i<colNames.size();i++)
	{
		if(i==0)
		{
			pFile <<colNames[i]<<"\tID\tRegion1\tRegion2";
		}
		else
		{
			pFile <<"\t"<<colNames[i];
		}
	}
	pFile <<endl;
	for(map<string,vector<double>*>::iterator pIter=allPairs_Features.begin();pIter!=allPairs_Features.end();pIter++)
	{
		int ID=allPairs[pIter->first];
		pFile <<pIter->first<<"\t"<< ID;
		PairIDs* regpairid=allPairs_RegionPairs[pIter->first];
		pFile <<"\t"<<regpairid->e<<"\t"<<regpairid->p;
		vector<double>* feats=allPairs_Features[pIter->first];
		for(int i=0;i<feats->size();i++)
		{
			pFile <<"\t" <<(*feats)[i];
		}
		pFile <<endl;
		/*int pid=allPairs[pIter->first];
		if(regionPairGraph.find(pid)==regionPairGraph.end())
		{
			cout <<"No graph for "<< pid << " "<< pIter->first <<endl;
			continue;
		}
		map<int,int>* neighbors=regionPairGraph[pid];
		for(map<int,int>::iterator nIter=neighbors->begin();nIter!=neighbors->end();nIter++)
		{
			oFile <<pid<<"\t" << nIter->first<<"\t"<<1 << endl;
		}*/
	}
	//oFile.close();
	pFile.close();
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
	if(argc!=8)
	{
		cout <<"Usage: graphGen featurefile start end  count field dist neighborsize outputsuffix" << endl;
		return 0;
	}
	Framework fw;
	fw.setMaxDist(atoi(argv[5]));
	fw.readPairs(argv[1],atoi(argv[2]),atoi(argv[3]),atoi(argv[4]));
	fw.setRegionNeighborSetSize(atoi(argv[6]));
	fw.generateNeighborGraph(argv[7]);
	return 0;
}
