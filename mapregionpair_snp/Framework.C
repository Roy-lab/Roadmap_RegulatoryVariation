#include <iostream>
#include <fstream>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "SNP.H"
#include "Node.H"
#include "HGNCConverter.H"
#include "Framework.H"

Framework::Framework()
{
}

Framework::~Framework()
{
}

int 
Framework::readGWAS(const char* aFName)
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
		string chromosome;
		string snpname;
		int position;
		while(tok!=NULL)
		{
			char* end=strchr(tok,'\t');
			if(end!=NULL)
			{
				*end='\0';
			}
			if(tokCnt==0)
			{
				snpname.append(tok);
			}
			if(tokCnt==1)
			{
				char chromKey[1024];
				if(strcmp(tok,"23")==0)
				{
					sprintf(chromKey,"chrX");
				}
				else
				{
					sprintf(chromKey,"chr%s",tok);
				}
				chromosome.append(chromKey);
			}
			else if(tokCnt==2)
			{
				position=atoi(tok);
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
		SNP* snp=new SNP;
		snp->name.append(snpname.c_str());
		snp->start=position;
		map<int,SNP*>* loci=NULL;
		if(snpSet.find(chromosome)==snpSet.end())
		{
			loci=new map<int,SNP*>;
			snpSet[chromosome]=loci;
		}
		else
		{
			loci=snpSet[chromosome];
		}
		(*loci)[position]=snp;
	}
	inFile.close();
	return 0;
}

int
Framework::readGeneNameMap(const char* aFName)
{
	hgncconverter.readGeneNames(aFName);
	return 0;
}

//This file is formatted as follows: region1\tregion2\tpredicted_cnt\tdistance\tqvalue\tgene1(several)\tgene2(several)
int 
Framework::readRegionPairs(const char* aFName)
{
	ifstream inFile(aFName);
	char* buffer=NULL;
	int buffLen=0;
	string buffstr;
	while(inFile.good())
	{
		getline(inFile,buffstr);
		if(buffstr.length()<=0)
		{
			continue;
		}
		if(buffstr.length()>=buffLen)
		{
			if(buffer==NULL)
			{
				delete[] buffer;
			}
			buffLen=buffstr.length()+1;
			buffer=new char[buffLen];
		}
		strcpy(buffer,buffstr.c_str());
		
		double count=0;
		double qvalue=0;
		int dist=0;
		Node* r1=NULL;
		Node* r2=NULL;
		string regName1;
		string regName2;
		char* tok=strtok(buffer,"\t");
		int tokCnt=0;
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				regName1.append(tok);
				if(nodeSet.find(regName1)==nodeSet.end())
				{
					r1=new Node;
					nodeSet[regName1]=r1;
					r1->setName(tok);
				}
				else
				{
					r1=nodeSet[regName1];
				}
			}
			else if(tokCnt==1)
			{
				regName2.append(tok);
				if(nodeSet.find(regName2)==nodeSet.end())
				{
					r2=new Node;
					nodeSet[regName2]=r2;
					r2->setName(tok);
				}
				else
				{
					r2=nodeSet[regName2];
				}
				
			}
			else if(tokCnt==2)
			{
				count=atof(tok);
			}	
			else if(tokCnt==3)
			{
				dist=atoi(tok);
			}
			else if(tokCnt==4)
			{
				qvalue=atof(tok);
			}
			else if(tokCnt==5 && strcmp(tok,"<nodata>")!=0)
			{
				r1->addGenes(tok);
			}
			else if(tokCnt==6 && strcmp(tok,"<nodata>")!=0)
			{
				r2->addGenes(tok);
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
		r1->setTarget(r2,count,dist,qvalue);
		r2->setTarget(r1,count,dist,qvalue);
		//Now store by chrom.
		map<int,Node*>* regsOnChr=NULL;
		if(nodeSetByChrom.find(r1->chrom)==nodeSetByChrom.end())
		{
			regsOnChr=new map<int,Node*>;
			nodeSetByChrom[r1->chrom]=regsOnChr;
		}
		else
		{
			regsOnChr=nodeSetByChrom[r1->chrom];
		}
		(*regsOnChr)[r1->start]=r1;
		regsOnChr=NULL;
		if(nodeSetByChrom.find(r2->chrom)==nodeSetByChrom.end())
		{
			regsOnChr=new map<int,Node*>;
			nodeSetByChrom[r2->chrom]=regsOnChr;
		}
		else
		{
			regsOnChr=nodeSetByChrom[r2->chrom];
		}
		(*regsOnChr)[r2->start]=r2;
	}
	inFile.close();
	return 0;
}

int 
Framework::mapSNPsToRegions(const char* aFSuff)
{
	char aFName[1024];
	sprintf(aFName,"%s_region_snp.txt",aFSuff);
	ofstream oFile(aFName);
	int totalRegionGeneAssocs_Pred=0;
	int hitRegionGeneAssocs=0;
	int totalRegionsWithSnps=0;
	map<string,int> snpsFound;
	for(map<string,map<int,Node*>*>::iterator cIter=nodeSetByChrom.begin();cIter!=nodeSetByChrom.end();cIter++)
	{
		map<int,Node*>* regionSet=cIter->second;
		if(snpSet.find(cIter->first)==snpSet.end())
		{
			continue;
		}
		for(map<int,Node*>::iterator nIter=regionSet->begin();nIter!=regionSet->end();nIter++)
		{
			Node* n=nIter->second;
			map<int,int> mappedSNP;
			map<int,SNP*>* loci=snpSet[cIter->first];
			int hit=0;
			for(map<int,SNP*>::iterator gIter=loci->begin();gIter!=loci->end();gIter++)
			{
				int dist=0;
				SNP* snp=gIter->second;
				int ig;
				//if(gIter->first>=(n->start-2000) && gIter->first <=(n->end+2000))
				if(gIter->first>=(n->start) && gIter->first <=(n->end))

				{
					n->setSNP(snp);
					hit++;
					snpsFound[snp->name]=0;
				}
			}
			//For every SNP mapped to this region, show all the pairs that get connected
			map<int,SNP*>& snpset=n->snpset;
			if(snpset.size()==0)
			{
				continue;
			}
			map<string,Node*>& targets=n->targets;
			totalRegionsWithSnps++;
			map<string,int>& genes=n->geneset;
			map<string,int> showGenes;

			for(map<int,SNP*>::iterator sIter=snpset.begin();sIter!=snpset.end();sIter++)
			{	
				SNP* snp=sIter->second;
				for(map<string,Node*>::iterator tIter=targets.begin();tIter!=targets.end();tIter++)
				{
					Node* m=tIter->second;
					double count=n->counts[tIter->first];
					double qval=n->qvals[tIter->first];
					oFile<<snp->name<<"\t"<< snp->start<<"\t"<<n->name<<"\t"<<m->name<<"\t"<< count <<"\t"<< qval<< "\t";
					for(map<string,int>::iterator geneIter=genes.begin();geneIter!=genes.end();geneIter++)
					{
						if(geneIter!=genes.begin())
						{
							oFile<<","<<hgncconverter.getCommonName((string&)geneIter->first);
						}
						else
						{
							oFile<<hgncconverter.getCommonName((string&)geneIter->first);
						}
						showGenes[geneIter->first]=0;
					}
					//Now get the neighbor region's gene set;
					map<string,int>& genes_m=m->geneset;
					for(map<string,int>::iterator geneIter=genes_m.begin();geneIter!=genes_m.end();geneIter++)
					{
						if(showGenes.size()==0)
						{
							oFile<<hgncconverter.getCommonName((string&)geneIter->first);
						}
						else if(showGenes.find(geneIter->first)==showGenes.end())
						{
							oFile<<","<<hgncconverter.getCommonName((string&)geneIter->first);
						}
						showGenes[geneIter->first]=0;
					}
					oFile<<endl;
					//Now showGenes has all the genes associated with this SNP
					for(map<string,int>::iterator geneIter=showGenes.begin();geneIter!=showGenes.end();geneIter++)
					{
						map<int,SNPEffect*>* snpEffectSet=NULL;
						string& gname=(string&)geneIter->first;
						if(geneSNPEffect.find(gname)==geneSNPEffect.end())
						{
							snpEffectSet=new map<int,SNPEffect*>;
							geneSNPEffect[gname]=snpEffectSet;
						}
						else
						{
							snpEffectSet=geneSNPEffect[gname];
						}
						SNPEffect* effect=new SNPEffect;
						effect->name.append(snp->name);
						effect->cnt=count;
						effect->qval=qval;
						(*snpEffectSet)[snp->start]=effect;
					}
					showGenes.clear();
				}
			}
		}
	}
	cout <<"Total regions with SNPs " << totalRegionsWithSnps <<endl;
	cout <<"Total SNPs connected " << snpsFound.size() <<endl;
	oFile.close();
	return 0;
}

//We will write two files: file1 will have snp-gene pair with count and qvalue, and file2 will have the avg qvalue because of multiple snps/interactions possibly
int
Framework::aggregateQvalues(const char* outSuffix)
{
	char oFName[1024];
	sprintf(oFName,"%s_snpgene_pair.txt",outSuffix);
	ofstream oFile1(oFName);
	sprintf(oFName,"%s_genescores.txt",outSuffix);
	ofstream oFile2(oFName);

	for(map<string,map<int,SNPEffect*>*>::iterator gIter=geneSNPEffect.begin();gIter!=geneSNPEffect.end();gIter++)
	{
		double avgLogQval=0;
		double avgCnt=0;
		map<int,SNPEffect*>* snps=gIter->second;
		for(map<int,SNPEffect*>::iterator sIter=snps->begin();sIter!=snps->end();sIter++)
		{
			SNPEffect* se=sIter->second;
			oFile1<<hgncconverter.getCommonName((string&)gIter->first)<<"\t"<< sIter->first<<"\t"<< se->cnt << "\t"<< se->qval<< endl;
			avgCnt=avgCnt+se->cnt;
			double qval_thresh=se->qval;
			if(se->qval<1e-100)
			{
				qval_thresh=1e-100;
			}
			double lqval=-1*log(qval_thresh);
			avgLogQval=avgLogQval+lqval;
		}
		oFile2<<hgncconverter.getCommonName((string&)gIter->first) << "\t" <<avgLogQval/((double)snps->size()) << "\t"<< avgCnt/((double)snps->size()) <<endl;
	}
	
	oFile1.close();
	oFile2.close();
	return 0;
}

int
main(int argc, const char** argv)
{
	if(argc!=5)
	{
		cout <<"Usage: mapSNPToEnhancers snp genenamemapping enhancerprom outputsuff" << endl;
		return 0;
	}
	Framework fw;
	fw.readGWAS(argv[1]);
	fw.readGeneNameMap(argv[2]);
	fw.readRegionPairs(argv[3]);
	fw.mapSNPsToRegions(argv[4]);
	fw.aggregateQvalues(argv[4]);
	return 0;
}
