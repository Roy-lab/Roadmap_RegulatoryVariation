#include <iostream>
#include <fstream>
#include <string.h>
#include "EdgeList.H"

EdgeList::EdgeList(AddType a)
{
	elist.clear();
	ecnt.clear();
	aType = a;
}

EdgeList::~EdgeList()
{
	clear();
}

int
EdgeList::clear()
{
	for(map<string,map<string,double>* >::iterator eiter=elist.begin(); eiter!=elist.end(); eiter++)
	{
		map<string,double>* vec = eiter->second;
		vec->clear();
		delete vec;
	}
	elist.clear();
	for(map<string,map<string,double>* >::iterator eiter=ecnt.begin(); eiter!=ecnt.end(); eiter++)
	{
		map<string,double>* vec = eiter->second;
		vec->clear();
		delete vec;
	}
	ecnt.clear();
	return 0;
}

int
EdgeList::addEdge(string tf, string tg, double v)
{
	map<string,double>* es;
	map<string,double>* ec;
	if(elist.find(tf) == elist.end())
	{
		es = new map<string,double>;
		elist[tf] = es;
		ec = new map<string,double>;
		ecnt[tf]  = ec;
	}
	else
	{
		es = elist[tf];
		ec = ecnt[tf];
	}
	if (es->find(tg) == es->end())
	{
		(*es)[tg] = v;
		(*ec)[tg] = 1;
	}
	else
	{
		double cur_c = (*ec)[tg];
		cur_c += 1;
		(*ec)[tg] = cur_c;

		double cur_v = (*es)[tg];
		if (aType == T_MAX)
		{
			if(cur_v < v)
			{
				(*es)[tg] = v;
			}
		}
		else if (aType == T_MIN)
		{
			if(cur_v > v)
			{
				(*es)[tg] = v;
			}
		}
		else if (aType == T_SUM || aType == T_AVG)
		{
			(*es)[tg] = v+cur_v;
		}
	}
	return 0;
}

int
EdgeList::writeToFile(const char* outFName)
{
	ofstream oFile(outFName);
	for(map<string,map<string,double>* >::iterator eiter=elist.begin(); eiter!=elist.end(); eiter++)
	{
		string tf = eiter->first;
		map<string,double>* es = eiter->second;
		map<string,double>* ec = ecnt[tf];
		for (map<string,double>::iterator titr=es->begin(); titr!=es->end(); titr++)
		{
			string tg = titr->first;
			double v  = titr->second;
			if (aType == T_AVG)
			{
				double c = (*ec)[tg];
				v = v/c;
			}
			oFile << tf << "\t" << tg << "\t" << v << endl;
		}
	}
	oFile.close();
	return 0;
}
