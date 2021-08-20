#include <iostream>
#include <algorithm>
#include <fstream>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>

#include "Randomizer.H"
#include "Kmeans.H"


Kmeans::Kmeans()
{
	withinDist=-1;
	acrossDist=-1;
}

Kmeans::~Kmeans()
{
}

int 
Kmeans::setConvergenceThreshold(double aThresh)
{
	threshold=aThresh;
	return 0;
}

int 
Kmeans::setClusterCnt(int aCnt)
{
	clusterCnt=aCnt;
	return 0;
}

int 
Kmeans::setDistanceType(Kmeans::DistanceType type)
{
	dtype=type;	
	return 0;
}

int 
Kmeans::cluster(KMEANS_SP_DATA* dataPtr)
{	
	dataSet=dataPtr;
	bool notDone=true;
	setInitialCenters_DataPts();
	formClusters();
	withinDist=computeWithinDistance();
	cout <<withinDist<< endl;
	int currIter=0;
	while((notDone) && (currIter<maxIter))
	{
		resetCenters();
		formClusters();
		double currDist=computeWithinDistance();
		cout << currDist << endl;
		double diff=fabs(withinDist-currDist);
		if(diff<=threshold)
		{
			notDone=false;
		}
		withinDist=currDist;
		currIter++;
	}
	cout <<"Found " << clusters.size() << endl;
	acrossDist=computeAcrossDistance();
	sortClusters();
	return 0;
}

//added by BB to sort based on the mean value
int 
Kmeans::sortClusters()
{

//for each of the clusters
vector<double> sortedMeans;
vector<double> unsortedMeans;
vector<int> clusterIDs;
//CLUSTERS sortedClusters;

for(CLUSTERS_ITER aIter=clusters.begin();aIter!=clusters.end();aIter++)
{
	int clusterID = aIter->first;
	vector<int> members = (*aIter->second);

	//loop over the members of the cluster
	double totalSum = 0.0;
	for(vector<int>::iterator it = members.begin(); it != members.end(); ++it) 
	{
    	INTDBLMAP* dPt=(*dataSet)[*it];
    	double memberValue =(*dPt)[0];
    	totalSum = totalSum+memberValue;
	}

	double clusterMean = totalSum/members.size();
	sortedMeans.push_back(clusterMean);
	unsortedMeans.push_back(clusterMean);
	clusterIDs.push_back(clusterID);
}

//sort the means
sort(sortedMeans.begin(), sortedMeans.end());

//loop over each member of the sorted clusterMean, find the value of unsorted (clusterMeans2) equal to it and use
int sortedClusterID = 1;
for(vector<double>::iterator it = sortedMeans.begin(); it != sortedMeans.end(); ++it)
{
	//the actual value of the iterator - sorted mean
	double value = *it; 
	//want to find the sorted mean in the usorted means
	int n=0;
	for(vector<double>::iterator jt = unsortedMeans.begin(); jt != unsortedMeans.end(); ++jt, ++n)
	{
		if (*jt == value)
		{
			int clusterID = clusterIDs[n];
			sortedClusters[sortedClusterID] = clusters[clusterID];
			sortedClusterID = sortedClusterID + 1;
			break;
		}
	} 
	//do we write over cluster?? 
}

//need to sort clusters from smallest to largest mean
return 0;
}

int
Kmeans::showClusters(const char* oFName)
{
	ofstream oFile(oFName);
	for(CLUSTERS_ITER cIter=clusters.begin();cIter!=clusters.end();cIter++)
	{
		oFile << "ClusterId " << cIter->first << endl;
		MEMBER* m=cIter->second;
		for(unsigned int i=0;i<m->size();i++)
		{
			int dId=(*m)[i];
			oFile <<" "<< dId;
		}
		oFile << endl;
	}
	oFile.close();
	return 0;
}

CLUSTERS&
Kmeans::getClusters()
{
	return sortedClusters;
}

double
Kmeans::getWithinClusterDist()
{
	return withinDist;
}

double
Kmeans::getAcrossClusterDist()
{
	return acrossDist;
}

int 
Kmeans::setInitialCenters()
{
	//Get the dimension of data
	//for each dimension find the max and min
	for(int d=0;d<maxDim;d++)
	{
		double max=0;
		double min=10000;
		for(KMEANS_SP_DATA_ITER kdIter=dataSet->begin();kdIter!=dataSet->end();kdIter++)
		{
			if(kdIter->second->find(d)==kdIter->second->end())
			{
				continue;
			}
			double dVal=(*kdIter->second)[d];
			if(dVal>max)
			{
				max=dVal;
			}
			if(dVal<min)
			{
				min=dVal;
			}
		}
		double incr=(max-min)/(double)clusterCnt;
		int cId=0;
		while(cId<clusterCnt)
		{
			double dVal=min+(cId*incr)+(0.1*incr);
			if(centers.find(cId)==centers.end())
			{
				INTDBLMAP* centerCoord=new INTDBLMAP;
				(*centerCoord)[d]=dVal;
				centers[cId]=centerCoord;
			}
			else
			{
				INTDBLMAP* centerCoord=centers[cId];
				(*centerCoord)[d]=dVal;
			}
			cId++;
		}
			
	}
	
	showCenters();
	return 0;
}


int 
Kmeans::setInitialCenters_DataPts()
{
	//Here we will use the datapoints as the centers of the clusters.
	int cId=0;
	map<int,int> centroidIndex;
	Randomizer rgen;
	rgen.initialize(0,dataSet->size()-1);

	while(cId<clusterCnt)
	{
		int rNo=rgen.getRandomNumber();
		while(centroidIndex.find(rNo)!=centroidIndex.end())
		{
			rNo=rgen.getRandomNumber();
		}
		INTDBLMAP* dPt=(*dataSet)[rNo];

		INTDBLMAP* centerCoord=new INTDBLMAP;
		for(INTDBLMAP_ITER idIter=dPt->begin();idIter!=dPt->end();idIter++)
		{
			(*centerCoord)[idIter->first]=idIter->second;
		}
		centers[cId]=centerCoord;
		cId++;
		cout <<"Setting cluster center " <<cId << " to datapoint " << rNo << endl;
	}
//	showCenters();		
	return 0;
}


//Go over the data 
int 
Kmeans::formClusters()
{
	for(CLUSTERS_ITER cIter=clusters.begin();cIter!=clusters.end();cIter++)
	{
		cIter->second->clear();
	}
	clusters.clear();

	//for each data point find the center to which it is closest
	for(KMEANS_SP_DATA_ITER kdIter=dataSet->begin();kdIter!=dataSet->end();kdIter++)
	{
		int dId=kdIter->first;
		double minDist=-1;
		int closestCenter=-1;
		INTDBLMAP* dataPt=kdIter->second;
		map<int,double> clusterDist;
		for(KMEANS_SP_DATA_ITER aIter=centers.begin();aIter!=centers.end();aIter++)
		{
			int cId=aIter->first;
			double pwDist=getPairwiseDist(aIter->second,dataPt);
			if(minDist==-1)
			{
				minDist=pwDist;
				closestCenter=cId;
			}
			else
			{
				if(dtype==Kmeans::EUCLID)
				{
					if(pwDist<=minDist)
					{
						minDist=pwDist;
						closestCenter=cId;
					}
				}
				else if(dtype==Kmeans::MI)
				{
					if(pwDist>=minDist)
					{
						minDist=pwDist;
						closestCenter=cId;
					}
				}
			}
			clusterDist[closestCenter]=minDist;
		}
		if(closestCenter==-1)
		{
			cout << "Did not find any center for a data point "<< dId << endl;
			exit(0);
		}
		
		if(clusters.find(closestCenter)==clusters.end())
		{
			MEMBER* m=new MEMBER;
			m->push_back(dId);
			clusters[closestCenter]=m;
			//clusters[smallestClusterID]=m;
		}
		else
		{
			//MEMBER* m=clusters[smallestClusterID];
			MEMBER* m=clusters[closestCenter];
			m->push_back(dId);
		}
	}
	cout <<"Found " << clusters.size() << " clusters " <<  endl;
	return 0;
}

int 
Kmeans::resetCenters()
{
	centers.clear();
	//cout <<"Resetting centers" <<endl;
	for(CLUSTERS_ITER cIter=clusters.begin();cIter!=clusters.end();cIter++)
	{
		MEMBER* m=cIter->second;
		INTDBLMAP* aCenter=new INTDBLMAP;
		for(unsigned int i=0;i<m->size();i++)
		{
			INTDBLMAP* dataPt=(*dataSet)[(*m)[i]];
			for(INTDBLMAP_ITER dIter=dataPt->begin();dIter!=dataPt->end();dIter++)
			{
				if(aCenter->find(dIter->first)==aCenter->end())
				{
					(*aCenter)[dIter->first]=dIter->second;
				}
				else
				{
					(*aCenter)[dIter->first]=(*aCenter)[dIter->first]+dIter->second;
				}
			}
		}
		for(INTDBLMAP_ITER dIter=aCenter->begin();dIter!=aCenter->end();dIter++)
		{
			dIter->second=dIter->second/m->size();
			if(dtype==Kmeans::MI)
			{
				//Project to discrete space
				//Hardcoding this for now
				double step=2.0/3.0;
				double rval=dIter->second;
				double projectedval=floor(rval/step);
			}
		}
		centers[cIter->first]=aCenter;
	}
	//showCenters();
	return 0;
}

int
Kmeans::showCenters()
{
	for(KMEANS_SP_DATA_ITER aIter=centers.begin();aIter!=centers.end();aIter++)
	{
		cout <<"Center" << aIter->first;
		INTDBLMAP * dMap=aIter->second;
		for(INTDBLMAP_ITER idIter=dMap->begin();idIter!=dMap->end();idIter++)
		{
			cout <<  " "<< idIter->first<<":"<<idIter->second;
		}
		cout << endl;
	}
	return 0;
}

double 
Kmeans::computeWithinDistance()
{
	double totalDist=0;
	for(CLUSTERS_ITER cIter=clusters.begin();cIter!=clusters.end();cIter++)
	{
		MEMBER* m=cIter->second;
		INTDBLMAP* aCenter=centers[cIter->first];
		for(unsigned int i=0;i<m->size();i++)
		{
			int dId=(*m)[i];
			INTDBLMAP* dataPt=(*dataSet)[dId];
			double dist=getPairwiseDist(dataPt,aCenter);
			totalDist=totalDist+dist;
		}
	}

	return totalDist;
}

double
Kmeans::computeAcrossDistance()
{
	double totalDist=0;
	for(CLUSTERS_ITER cIter=clusters.begin();cIter!=clusters.end();cIter++)
	{
		MEMBER* m=cIter->second;
		for(KMEANS_SP_DATA_ITER idIter=centers.begin();idIter!=centers.end();idIter++)
		{
			if(idIter->first==cIter->first)
			{
				continue;
			}
			INTDBLMAP* aCenter=centers[idIter->first];
			for(unsigned int i=0;i<m->size();i++)
			{
				int dId=(*m)[i];
				INTDBLMAP* dataPt=(*dataSet)[dId];
				double dist=getPairwiseDist(dataPt,aCenter);
				totalDist=totalDist+dist;
			}
		}
	}
	totalDist=totalDist/((double) clusters.size());
	return totalDist;
}


int
Kmeans::clear()
{
	for(KMEANS_SP_DATA_ITER ksIter=centers.begin();ksIter!=centers.end();ksIter++)
	{
		ksIter->second->clear();
		delete ksIter->second;
	}
	centers.clear();
	for(CLUSTERS_ITER cIter=clusters.begin();cIter!=clusters.end();cIter++)
	{
		MEMBER* m=cIter->second;
		m->clear();
		delete m;
	}
	clusters.clear();
	return 0;
}

int
Kmeans::setMaxIter(int amax)
{
	maxIter=amax;
	return 0;
}

//This is where I can use a Distance class and use many different types of distance metrics
double 
Kmeans::getPairwiseDist(DBLVECT* aVect,DBLVECT* bVect)
{
	double dist=0;
	double sum1=0;
	double sum2=0;
	for(unsigned int i=0;i<aVect->size();i++)
	{
		double diff=fabs((*aVect)[i]-(*bVect)[i]);
		sum1=sum1+((*aVect)[i]*(*aVect)[i]);
		sum2=sum2+((*bVect)[i]*(*bVect)[i]);
		dist=dist+(diff*diff);
		//dist=dist+((*aVect)[i]*(*bVect)[i]);
	}
	//dist=dist/(sqrt(sum1)*sqrt(sum2));
	return dist;
}

double
Kmeans::getPairwiseDist(INTDBLMAP* aMap, INTDBLMAP* bMap)
{
	double distance=-1;
	switch(dtype)
	{
		case Kmeans::EUCLID:
		{
			distance=getPairwiseDist_Euclid(aMap,bMap);
			break;
		}
		case Kmeans::MI:
		{
			//distance=getPairwiseDist_MI(aMap,bMap);
			distance=getPairwiseDist_CC(aMap,bMap);
			break;
		}
		default:
		{
			cout <<"No Distance type set!!"<< endl;
			exit(0);
		}
	}
	return distance;
}

double
Kmeans::getPairwiseDist_Euclid(INTDBLMAP* aMap, INTDBLMAP* bMap)
{
	double dist=0;
	for(INTDBLMAP_ITER aIter=aMap->begin();aIter!=aMap->end();aIter++)
	{
		double diff;
		if(bMap->find(aIter->first)!=bMap->end())
		{
			diff=aIter->second-(*bMap)[aIter->first];
		}
		else
		{
			diff=aIter->second;
		}
		dist=dist+(diff*diff);
	}

	
	for(INTDBLMAP_ITER aIter=bMap->begin();aIter!=bMap->end();aIter++)
	{
		//If this dimension is occurs in the first map, we have already taken that into account
		if(aMap->find(aIter->first)!=aMap->end())
		{
			continue;
		}
		double diff=aIter->second;
		dist=dist+(diff*diff);
	}
	return dist;
}


double
Kmeans::getPairwiseDist_CC(INTDBLMAP* aMap, INTDBLMAP* bMap)
{
	double cc=0;
	double m1=0;
	for(INTDBLMAP_ITER aIter=aMap->begin();aIter!=aMap->end();aIter++)
	{
		m1=m1+aIter->second;
	}
	m1=m1/aMap->size();

	double m2=0;
	for(INTDBLMAP_ITER aIter=bMap->begin();aIter!=bMap->end();aIter++)
	{
		m2=m2+aIter->second;
	}
	m2=m2/bMap->size();
	
	double xx=0;
	double yy=0;
	double xy=0;
	for(INTDBLMAP_ITER aIter=aMap->begin();aIter!=aMap->end();aIter++)
	{
		double diff1=aIter->second-m1;
		xx=xx + (diff1 * diff1);
		INTDBLMAP_ITER bIter=bMap->find(aIter->first);
		double diff2=bIter->second-m2;
		yy=yy+ (diff2 * diff2);
		xy=xy+(diff1*diff2);  	
	}
	cc=sqrt((xy*xy)/(xx*yy));
	return fabs(cc);
}

//Assume that this data is discrete and the DBL values are actually integers
double
Kmeans::getPairwiseDist_MI(INTDBLMAP* aMap,INTDBLMAP* bMap)
{
	map<string,double> jointConf;
	map<int,double> marginal1;
	map<int,double> marginal2;
	for(INTDBLMAP_ITER dIter=aMap->begin();dIter!=aMap->end();dIter++)
	{
		int v1=(int) ceil(dIter->second);
		int v2=(int) ceil((*bMap)[dIter->first]);
		char key[256];
		sprintf(key,"%d-%d",v1,v2);
		string keystr(key);
		if(jointConf.find(keystr)==jointConf.end())
		{
			jointConf[keystr]=1.0;
		}
		else
		{
			jointConf[keystr]=jointConf[keystr]+1.0;
		}
		if(marginal1.find(v1)==marginal1.end())
		{
			marginal1[v1]=1.0;
		}
		else
		{
			marginal1[v1]=marginal1[v1]+1.0;
		}

		if(marginal2.find(v2)==marginal2.end())
		{
			marginal2[v2]=1.0;
		}
		else
		{
			marginal2[v2]=marginal2[v2]+1.0;
		}
	}
	double total=(double) aMap->size();
	double jtEntropy=0;
	for(map<string,double>::iterator aIter=jointConf.begin();aIter!=jointConf.end();aIter++)
	{
		double pval=aIter->second/total;
		jtEntropy=jtEntropy+(pval*log(pval));
	}
	double margEntropy1=0;
	double margEntropy2=0;
	for(map<int,double>::iterator aIter=marginal1.begin();aIter!=marginal1.end();aIter++)
	{
		double pval=aIter->second/total;
		margEntropy1=margEntropy1+(pval*log(pval));
	}	

	for(map<int,double>::iterator aIter=marginal2.begin();aIter!=marginal2.end();aIter++)
	{
		double pval=aIter->second/total;
		margEntropy2=margEntropy2+(pval*log(pval));
	}
	//This formula is actually H(X) + H(Y) - H(X,Y). But since we have not
	//negated the entropies the formula should be what is below
	double mi=jtEntropy-margEntropy1-margEntropy2;
	jointConf.clear();
	marginal1.clear();
	marginal2.clear();
	return mi;
}