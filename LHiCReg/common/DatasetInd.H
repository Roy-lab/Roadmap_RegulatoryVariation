#ifndef _DATASET_IND
#define _DATASET_IND
#include <map>
#include <vector>
using namespace std;
class DatasetInd
{
	public:
	DatasetInd();
	~DatasetInd();
	int setDatasetOrder(int,vector<int>&);
	vector<int>* getDatasetOrder(int);
	private:
	//We will store all the sorted dataset IDs per attribute here
	map<int,vector<int>*> sortOrders;
};
#endif
