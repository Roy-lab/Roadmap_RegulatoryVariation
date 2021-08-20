
#include <iostream>
#include <fstream>
#include <string.h>
#include "GeneNameMapper.H"


GeneNameMapper::GeneNameMapper()
{
}

GeneNameMapper::~GeneNameMapper()
{
}

int
GeneNameMapper::readGeneNames(const char* aFName)
{
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
		string commName;
		string orfName;
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
				commName.append(begin);
				char cName[100];
				int len=strlen(begin);
				int i=0;
				while(i<len)
				{
					cName[i]=tolower(begin[i]);
					i++;
				}
				cName[len]='\0';
			}
			else if(tokCnt==1)
			{
				orfName.append(begin);
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
		
		orfToCommon[orfName]=commName;		
	}
	inFile.close();
	return 0;
}


const char*
GeneNameMapper::getCommonName(const char* aGene)
{
	string key(aGene);
	if(orfToCommon.find(key)==orfToCommon.end())
	{
		return aGene;
	}
	return orfToCommon[key].c_str();
}

