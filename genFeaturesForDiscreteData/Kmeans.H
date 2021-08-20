/*This class implements the functionality of the K means clustering in 1d.
 * The user of this class provides a 1d data and gets back k clusters of the
 * data*/
#include <vector>
#include <map>

using namespace std;

typedef vector<int> INTVECT;
typedef vector<string> STRVECT;
typedef vector<double> DBLVECT;

typedef map<string,int> STRINTMAP;
typedef map<string,int>::iterator STRINTMAP_ITER;

typedef map<int,double> INTDBLMAP;
typedef map<int,double>::iterator INTDBLMAP_ITER;

typedef map<string,double> STRDBLMAP;
typedef map<string,double>::iterator STRDBLMAP_ITER;

typedef map<int,double> INTDBLMAP;
typedef map<int,double>::iterator INTDBLMAP_ITER;

typedef map<int,int> INTINTMAP;
typedef map<int,int>::iterator INTINTMAP_ITER;
#define CONSTR_LEN 256

typedef map<int,INTDBLMAP*> KMEANS_SP_DATA;
typedef map<int,INTDBLMAP*>::iterator KMEANS_SP_DATA_ITER;


typedef vector<int> MEMBER;
typedef map<int,MEMBER*> CLUSTERS;
typedef map<int,MEMBER*>::iterator CLUSTERS_ITER;



class Kmeans
{
	public:
		Kmeans();
		~Kmeans();
		int setConvergenceThreshold(double);
		int setClusterCnt(int);
		int setMaxDim(int);
	
		typedef enum
		{
			EUCLID,
			MI
		}DistanceType;
		int setDistanceType(DistanceType);
		int cluster(KMEANS_SP_DATA*);
		int showClusters(const char*);
		int sortClusters();
		CLUSTERS& getClusters();
		//This is the sum of all within cluster simillarity
		double getWithinClusterDist();
		double getAcrossClusterDist();
		double computeWithinDistance();
		double computeAcrossDistance();
		int clear();
		int setMaxIter(int);
	private:
		int setInitialCenters();
		int setInitialCenters_DataPts();
		int formClusters();
		int resetCenters();
	
		int showCenters();
		
		double getPairwiseDist(DBLVECT*,DBLVECT*);
		double getPairwiseDist(INTDBLMAP*,INTDBLMAP*);
		double getPairwiseDist_Euclid(INTDBLMAP*, INTDBLMAP*);
		double getPairwiseDist_CC(INTDBLMAP*, INTDBLMAP*);
		double getPairwiseDist_MI(INTDBLMAP*, INTDBLMAP*);
	
		
		int clusterCnt;
		double threshold;
		double withinDist;
		double acrossDist;
		int maxDim;
		int maxIter;
			
		KMEANS_SP_DATA centers;
		CLUSTERS clusters;
		CLUSTERS sortedClusters;
		KMEANS_SP_DATA* dataSet;
		DistanceType dtype;
};
