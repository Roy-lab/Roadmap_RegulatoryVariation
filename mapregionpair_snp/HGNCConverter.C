#include <iostream>
#include <fstream>
#include <string.h>

#include "HGNCConverter.H"

HGNCConverter::HGNCConverter()
{
}

HGNCConverter::~HGNCConverter()
{
}

int
HGNCConverter::readGeneNames(const char* aFName)
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
		char* tok=buffer;
		int tokCnt=0;
		string geneName;
		string ensembleName1;
		string ensembleName2;
		while(tok!=NULL)
		{
			char* end=strchr(tok,'\t');
			if(end!=NULL)
			{
				*end='\0';
			}
			if(tokCnt==1)
			{
				geneName.append(tok);
			}
			else if(tokCnt==18)
			{
				ensembleName1.append(tok);
			}
			else if(tokCnt==36)
			{
				ensembleName2.append(tok);
			}
			if(end!=NULL)
			{
				tok=end+1;
			}
			else
			{
				tok=NULL;
			}
			tokCnt++;
		}
		if(ensembleName1.length()<=0 && ensembleName2.length()<=0)
		{
			continue;
		}
		if(ensembleName1.length()>0)
		{
			if(ensembleToCommon.find(ensembleName1)==ensembleToCommon.end())
			{
				ensembleToCommon[ensembleName1]=geneName;
			}
			else
			{
				string& old=ensembleToCommon[ensembleName1];
				if(strcmp(old.c_str(),geneName.c_str())!=0)
				{
					//cout <<"Mismatch of genenames1 "<< old.c_str() << "\t" << geneName.c_str() << endl;
				}
			}
		}
		if(ensembleName2.length()>0)
		{
			if(ensembleToCommon.find(ensembleName2)==ensembleToCommon.end())
			{
				ensembleToCommon[ensembleName2]=geneName;
			}
			else
			{
				string& old=ensembleToCommon[ensembleName2];
				if(strcmp(old.c_str(),geneName.c_str())!=0)
				{
					//cout <<"Mismatch of genenames2 "<< old.c_str() << "\t" << geneName.c_str() << endl;
				}
			}
		}
	}
	inFile.close();
	return 0;
}

const char*
HGNCConverter::getCommonName(string& key)
{
	if(ensembleToCommon.find(key)==ensembleToCommon.end())
	{
		return key.c_str();
	//	return NULL;
	}
	return ensembleToCommon[key].c_str();
}
