#include <iostream>
#include <fstream>
#include <string.h>
#include "MotifNameMapper.H"

MotifNameMapper::MotifNameMapper()
{
	mot2tf.clear();
}

MotifNameMapper::~MotifNameMapper()
{
	clear();
}

int
MotifNameMapper::clear()
{
	for(map<string,vector<string>* >::iterator miter=mot2tf.begin(); miter!=mot2tf.end(); miter++)
	{
		vector<string>* vec = miter->second;
		vec->clear();
		delete vec;
	}
	mot2tf.clear();
	return 0;
}

vector<string>*
MotifNameMapper::splitTFs(char* tfnames)
{
	vector<string>* vec = new vector<string>;
	char* begin=tfnames;
	char* end=tfnames;

	while(begin!=NULL)
	{
		end=strchr(begin,':');
		if(end!=NULL)
		{
			*end='\0';
		}
		string name(begin);
		vec->push_back(name);
		if(end!=NULL)
		{
			begin=end+2;
		}
		else
		{
			begin=NULL;
		}
	}
	return vec;
}

int
MotifNameMapper::readMap(const char* aFName)
{
	clear();
	ifstream inFile(aFName);
	char buffer[4096];

	while(inFile.good())
	{
		inFile.getline(buffer,4095);
		if(strlen(buffer)<=0)
		{
			continue;
		}
			
		int tokCnt=0;
		string motname;
		vector<string>* vec = NULL;
		char* begin=buffer;
		char* end=buffer;

		while(begin!=NULL)
		{
			end=strchr(begin,'\t');
			if(end!=NULL)
			{
				*end='\0';
			}
			if(tokCnt==0)
			{
				motname.append(begin);
			}
			else if(tokCnt==1)
			{
				vec = splitTFs(begin);
			}
			tokCnt++;
			if(end!=NULL)
			{
				begin=end+1;
			}
			else
			{
				begin=NULL;
			}
		}
		mot2tf[motname]=vec;
	}
	inFile.close();
	return 0;
}


vector<string>*
MotifNameMapper::getTFName(const char* mot)
{
	string key(mot);
	if(mot2tf.find(key)==mot2tf.end())
	{
		vector<string>* vec = new vector<string>;
		vec->push_back(key);
		mot2tf[key] = vec;
		//return NULL;
	}
	return mot2tf[key];
}

/*
int
main(int argc, char** argv)
{
	MotifNameMapper* motmap = new MotifNameMapper;
	motmap->readMap(argv[1]);
	vector<string>* vec = motmap->getTFName(argv[2]);
	if (vec == NULL)
	{
		cout << "Motif name does not exist: " << argv[2] << endl;
	}
	else
	{
		for (int i=0;i<vec->size();i++)
		{
			cout << vec->at(i) << endl;
		}
	}
	return 0;
}
*/
