#include <iostream>
#include <fstream>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include "Gene.H"
#include "GFFReader.H"

GFFReader::GFFReader(TSSType t)
{
	tType = t;
}

GFFReader::~GFFReader()
{
}

int 
GFFReader::readGFFFile(const char* aFName)
{
	//cout << "Touch GFFReader::readGFFFile\n";
	ifstream inFile(aFName);
	string buffstr;
	char* buffer=new char[4096];
	int bufflen=4096;
	while(inFile.good())
	{
		//cout << "In While\n";
		getline(inFile,buffstr);
		if(buffstr.length()<=0)
		{
			continue;
		}
		if(bufflen<=buffstr.length())
		{
			delete[] buffer;
			bufflen=buffstr.length()+1;
			buffer=new char[bufflen];
		}
		strcpy(buffer,buffstr.c_str());
		if(buffer==NULL)cout << "Buffer is NULL\n";
		if(buffstr.c_str()==NULL)cout << "Buffer string is NULL\n";
		//if(strstr(buffer,"KNOWN")==NULL)cout << "Warning: strstr KNOWN is returning NULL!\n";
		if(strstr(buffer,"KNOWN")==NULL)
		{
			//continue;
		}
		//cout << "Touch GFFReader::readGFFFile II\n";
		int tokCnt=0;
		char* tok=strtok(buffer,"\t");
		int startCoordinate;
		int endCoordinate;
		string id;
		vector<string> all_gene_ids;
		string featureType;
		string chromosome;
		char strand;
		string src;
		bool known=false;
		bool pseudo=false;
		string gene_name;
		string gene_id;
		string transcript_id;
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				chromosome.append(tok);
				////cout << chromosome << endl;
			}
			else if(tokCnt==1)
			{
				src.append(tok);		
				//cout <<"src variable " << src << "\n";
			}
			else if(tokCnt==2)
			{
				featureType.append(tok);
				//cout << "featureType " << featureType << "\n";
			}
			else if(tokCnt==3)
			{
				startCoordinate=atoi(tok);
				//cout << "Start coordinate " << startCoordinate << "\n";
			}
			else if(tokCnt==4)
			{
				endCoordinate=atoi(tok);
			}
			//was ==6 for previous gff, it's also position 6 in the present file
			else if(tokCnt==6)
			{
				strand=tok[0];
				//cout << "Strand "<< strand << "\n";
			}
			else if(tokCnt==8)
			{
				if(strstr(tok,"KNOWN")!=NULL)
				{
					known=true;
				}
				if(strstr(tok,"NOVEL")!=NULL || strstr(tok,"PUTATIVE")!=NULL)
				{
					pseudo=true;
				}
				//if(pseudo) cout << "Warning dubious marker selection\n";
				//char* startpos=strstr(tok,"trlist");
				//char* stoppos=strstr(tok,"trbiotlist");
				//char* stoppos=strstr(tok,"transcript_id");
				//cout << "Touch ID assessment\n";
				
				////gene_id ENSMUSG00000056763 trlist ENSMUST00000143606, 
				if(tType == T_GENE)
				{  
					char* startpos=strstr(tok,"gene_id");
					char* stoppos=strstr(tok,"trlist");
					id = string(startpos+8,stoppos-1);
					all_gene_ids.push_back(id);
					//id.append(startpos+9,stoppos-3);
					////cout << "id set to " << id << "\n";
				}
				if(tType == T_TRANS)
				{
					char* startpos=strstr(tok,"trlist");
					char* stoppos=strstr(startpos+7," ");
					char* bpos=startpos+7;
					char* epos=strstr(bpos,",");
					while(epos != NULL)
					{
						id = string(bpos,epos);
						all_gene_ids.push_back(id);
						bpos=epos+1;
						epos=strstr(bpos,",");
						if (bpos>=stoppos)
						{
							break;
						}
					}
					//id.append(startpos+9,stoppos-3);
					////cout << "id set to " << id << "\n";
				}
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
		Gene* gene=NULL;
		if(strand=='-')
		{
			int temp=startCoordinate;
			startCoordinate=endCoordinate;
			endCoordinate=temp;	
		}
		//if(1)//strcmp(featureType.c_str(),"transcript")==0)
		for(int iid=0; iid<all_gene_ids.size(); iid++)
		{
			id = all_gene_ids[iid];
			if(geneSet.find(id)!=geneSet.end())
			{
				////cout <<"Gene " << id << " already exists" << endl;
				if(strand=='+')
				{
					int cStart=geneSet.find(id)->second->geneStart;
					int cEnd=geneSet.find(id)->second->geneEnd;
					int nStart;
					int nEnd;
					////cout << "Potential new coordinates are " << startCoordinate << " and " << endCoordinate << endl;
					////cout << "Current Coordinates are " << cStart << " and " << cEnd << endl;
					if(startCoordinate<cStart)nStart=startCoordinate;
					else nStart=cStart;
					if(endCoordinate>cEnd)nEnd=endCoordinate;
					else nEnd=cEnd;
					if(!nStart || !nEnd)cout << "Warning new range values not defined" <<   endl;
					geneSet.find(id)->second->setGeneCoordinate(nStart,nEnd);
					////cout << "Updated Coordinates are " << geneSet.find(id)->second->geneStart << " and " << geneSet.find(id)->second->geneEnd << endl;
				}
				else
				{
					int cStart=geneSet.find(id)->second->geneStart;
					int cEnd=geneSet.find(id)->second->geneEnd;
					int nStart;
					int nEnd;
					////cout << "Potential new coordinates are " << startCoordinate << " and " << endCoordinate << endl;
					////cout << "Current Coordinates are " << cStart << " and " << cEnd << endl;
					if(startCoordinate>cStart)nStart=startCoordinate;
					else nStart=cStart;
					if(endCoordinate<cEnd)nEnd=endCoordinate;          
					else nEnd=cEnd;
					if(!nStart || !nEnd)cout << "Warning new range values not defined" <<   endl;
					geneSet.find(id)->second->setGeneCoordinate(nStart,nEnd);          
					////cout << "Updated Coordinates are " << geneSet.find(id)->second->geneStart << " and " << geneSet.find(id)->second->geneEnd << endl;
				}
				continue;
			}
			////cout << "new Gene " << id << " and geneSet.size()is " << geneSet.size() << endl;
			gene=new Gene;
			gene->setName(id.c_str());
			geneSet[id]=gene;
			gene->setStrand(strand);
			gene->setGeneCoordinate(startCoordinate,endCoordinate);
			gene->setChromosome(chromosome.c_str());
			map<string,Gene*>* chromGenes=NULL;
			if(geneSetPerChrom.find(chromosome)==geneSetPerChrom.end())
			{
				chromGenes=new map<string,Gene*>;
				geneSetPerChrom[chromosome]=chromGenes;
			}
			else
			{
				chromGenes=geneSetPerChrom[chromosome];
			}
			(*chromGenes)[id]=gene;
			char strandChrom[256];
			sprintf(strandChrom,"%s_%c",chromosome.c_str(),strand);
			string key(strandChrom);
			map<int,Gene*>* chromstrandGenes=NULL;
			if(geneSetPerChromStrand.find(key)==geneSetPerChromStrand.end())
			{
				chromstrandGenes=new map<int,Gene*>;
				geneSetPerChromStrand[key]=chromstrandGenes;
			}
			else
			{
				chromstrandGenes=geneSetPerChromStrand[key];
			}
			(*chromstrandGenes)[startCoordinate]=gene;
		}
	}
	delete [] buffer;
	inFile.close();
	return 0;
}

map<string,map<string,Gene*>*>&
GFFReader::getGeneSet()
{
	return geneSetPerChrom;
}

int
GFFReader::setNeighbors()
{
	for(map<string,map<int,Gene*>*>::iterator cIter=geneSetPerChromStrand.begin();cIter!=geneSetPerChromStrand.end();cIter++)
	{
		map<int,Gene*>* geneSet=cIter->second;
		////cout <<"Genes in " << cIter->first << " " << geneSet->size() << endl;
		map<int,Gene*>::iterator gIter=geneSet->begin();
		while(gIter!=geneSet->end())
		{
			map<int,Gene*>::iterator fIter=gIter;
			fIter++;
			Gene* g1=gIter->second;
			if(fIter!=geneSet->end())
			{
				Gene* g2=fIter->second;
				g1->downstreamNeighbor=g2;
				g2->upstreamNeighbor=g1;
			}
			gIter=fIter;
		}
	}
}


map<string,Gene*>* 
GFFReader::getGeneSetForChromosome(string& chr)
{
	map<string,Gene*>* geneSetOnChrom=NULL;
	if(geneSetPerChrom.find(chr)!=geneSetPerChrom.end())
	{
		geneSetOnChrom=geneSetPerChrom[chr];
	}
	return geneSetOnChrom;	
}

