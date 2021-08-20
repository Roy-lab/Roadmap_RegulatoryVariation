#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include "Gene.H"
#include "GFFReader.H"
#include "GeneNameMapper.H"
#include "MotifNameMapper.H"
#include "EdgeList.H"
#include "Framework.H"

Framework::Framework()
{
}

Framework::~Framework()
{
}

int 
Framework::readChromMap(const char* aFName)
{	
	ifstream inFile(aFName);
	char buffer[1024];
	while(inFile.good())
	{
		inFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		string srcChrom;
		string targetChrom;
		char*  tok=strtok(buffer,"\t");
		int tokCnt=0;
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				srcChrom.append(tok);
			}
			else if(tokCnt==1)
			{
				targetChrom.append(tok);
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
		srcTargetChromMap[srcChrom]=targetChrom;
	}
	inFile.close();
	return 0;
}

int
Framework::readSites()
{
	/*
	if(strcmp(type,"single")==0)
	{
		readSitesSingle(aFName);
	}
	else if(strcmp(type,"double")==0)
	{
		readSitesDouble(aFName);
	}
	*/
	readSitesBed(motname);
	return 0;
}

int 
Framework::readSitesSingle(const char* aFName)
{
	ifstream inFile(aFName);
	char buffer[1024];
	int lineCnt=0;
	while(inFile.good())
	{
		inFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		string chr;
		int pos=-1;
		char* tok=strtok(buffer,"\t");
		int tokCnt=0;
		Framework::SiteInfo*  site=new Framework::SiteInfo;
		while(tok!=NULL)
		{
			switch(tokCnt)
			{
				case 0:
				{
					string srcKey(tok);
					if(srcTargetChromMap.find(srcKey)==srcTargetChromMap.end())
					{
						cout <<"No chromosome found for " << srcKey<< endl;
						srcTargetChromMap[srcKey] = srcKey;
					}	
					site->chrom.append(srcTargetChromMap[srcKey]);
					break;
				}
				case 1:
				{
					site->pos=atoi(tok);
					break;
				}
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
		siteSet.push_back(site);
	}
	inFile.close();
	return 0;
}

int 
Framework::readSitesDouble(const char* aFName)
{
	ifstream inFile(aFName);
	char buffer[1024];
	int lineCnt=0;
	while(inFile.good())
	{
		inFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		string chr;
		int pos=-1;
		int begin=0;
		int end=0;
		string name;
		char* tok=strtok(buffer," \t");
		int tokCnt=0;
		Framework::SiteInfo*  site=new Framework::SiteInfo;
		while(tok!=NULL)
		{
			char* tempTok=tok;
			char* pos=strchr(tok,'"');
			if(pos!=NULL)
			{
				tempTok=tok+1;
			}
			char* endpos=strchr(tempTok,'"');
			if(endpos!=NULL)
			{
				*endpos='\0';
			}
			switch(tokCnt)
			{
				case 0:
				{
					string srcKey(tempTok);
					if(srcTargetChromMap.find(srcKey)==srcTargetChromMap.end())
					{
						cout <<"No chromosome found for " << srcKey<< endl;
					}	
					else
					{
						site->chrom.append(srcTargetChromMap[srcKey]);
					}
					break;
				}
				case 1:
				{
					begin=atoi(tempTok);
					break;
				}
				case 2:
				{
					end=atoi(tempTok);
					break;
				}
				case 4:
				{
					name.append(tok);
					break;
				}		
			}
			tok=strtok(NULL," \t");
			tokCnt++;
		}
		if(site->chrom.length()==0)
		{
			delete site;
			continue;
		}
		pos=begin+(end-begin)/2;
		site->pos=pos;
		int length=(end-begin);
		site->length=length;
		site->name=name;
		siteSet.push_back(site);
	}
	inFile.close();
	return 0;
}

int
Framework::readBedConfig(const char* aFName)
{
	mySiteConf.name_index = -1;
	mySiteConf.chr_index = -1;
	mySiteConf.beg_index = -1;
	mySiteConf.end_index = -1;
	mySiteConf.score_index = -1;
	mySiteConf.strand_index = -1;

	ifstream inFile(aFName);
	char buffer[1024];
	while(inFile.good())
	{
		inFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		if(buffer[0]=='#')
		{
			continue;
		}
		string name;
		int index;
		char* tok=strtok(buffer," \t");
		int tokCnt=0;
		while(tok!=NULL)
		{
			switch(tokCnt)
			{
				case 0:
				{
					name = string(tok);
					break;
				}
				case 1:
				{
					index=atoi(tok);
					break;
				}
			}
			tok=strtok(NULL," \t");
			tokCnt++;
		}
		if(name == "name")
		{
			mySiteConf.name_index = index;
		}
		if(name == "chr")
		{
			mySiteConf.chr_index = index;
		}
		if(name == "beg")
		{
			mySiteConf.beg_index = index;
		}
		if(name == "end")
		{
			mySiteConf.end_index = index;
		}
		if(name == "score")
		{
			mySiteConf.score_index = index;
		}
		if(name == "strand")
		{
			mySiteConf.strand_index = index;
		}
	}
	inFile.close();
	return 0;
}

int 
Framework::readSitesBed(const char* aFName)
{
	//0		1		2		3	4			5
	//chr1	3004298	3004313	+	8.830266	M6519_1.02
	ifstream inFile(aFName);
	char buffer[1024];
	int lineCnt=0;
	while(inFile.good())
	{
		inFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		if(buffer[0]=='#')
		{
			continue;
		}
		string chr;
		int pos=-1;
		int begin=-1;
		int end=-1;
		double score=-1;
		string name;
		string strand;
		char* tok=strtok(buffer," \t");
		int tokCnt=0;
		Framework::SiteInfo*  site=new Framework::SiteInfo;
		while(tok!=NULL)
		{
			if(tokCnt == mySiteConf.chr_index)
			{
				string srcKey(tok);
				if(srcTargetChromMap.find(srcKey)==srcTargetChromMap.end())
				{
					////cout <<"No chromosome found for " << srcKey << endl;
					srcTargetChromMap[srcKey]=srcKey;
				}
				site->chrom.append(srcTargetChromMap[srcKey]);
			}
			if(tokCnt == mySiteConf.beg_index)
			{
				begin=atoi(tok);
			}
			if(tokCnt == mySiteConf.end_index)
			{
				end=atoi(tok);
			}
			if(tokCnt == mySiteConf.score_index)
			{
				score = atof(tok);
			}
			if(tokCnt == mySiteConf.name_index)
			{
				name.append(tok);
			}
			if(tokCnt == mySiteConf.strand_index)
			{
				strand.append(tok);
			}		
			tok=strtok(NULL," \t");
			tokCnt++;
		}
		if(site->chrom.length()==0)
		{
			delete site;
			continue;
		}
		if(end==-1)
		{
			end = begin;
		}
		if(begin==-1)
		{
			cerr << "begin position missing!" << endl;
			cerr << "Exiting!" << endl;
			exit(0);
		}
		pos=begin+(end-begin)/2;
		site->pos=pos;
		int length=(end-begin);
		site->length=length;
		site->name=name;
		site->score = score;
		site->strand = strand;
		siteSet.push_back(site);
	}
	inFile.close();
	return 0;
}

int 
Framework::readGFFFile()
{
	gffreader.readGFFFile(tssname);
	gffreader.setNeighbors();
	return 0;
}

int
Framework::readRefFasta(const char* aFName, const char* chrnamemapping)
{
	char buffer[1024];
	ifstream cmapFile(chrnamemapping);
	map<string,string> nameIDMap;
	while(cmapFile.good())
	{
		cmapFile.getline(buffer,1023);
		char* tok=strtok(buffer,"\t");
		int tokcnt=0;
		string key1;
		string key2;
		while(tok!=NULL)
		{
			if(tokcnt==0)
			{
				key1.append(tok);		
			}
			else if(tokcnt==1)
			{
				key2.append(tok);
			}
			tok=strtok(NULL,"\t");
			tokcnt++;
		}
		nameIDMap[key1]=key2;
	}
	string bases;
	string chromosome;
	ifstream inFile(aFName);
	while(inFile.good())
	{
		inFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		if(strstr(buffer,">")!=NULL)
		{
			if(bases.length()>0)
			{
				refchromLength[chromosome]=bases.length();
			}
			bases.clear();
			chromosome.clear();
			char* tok=strchr(buffer,' ');
			string key;
			if(tok!=NULL)
			{
				*tok='\0';
			}
			key.append(buffer+1);
			if(nameIDMap.find(key)==nameIDMap.end())
			{
				cout <<"No chromosome " << key << endl;
				exit(0);
			}
			chromosome.append(nameIDMap[key]);
		}
		else
		{
			bases.append(buffer);
		}
	}
	if(bases.length()>0)
	{
		refchromLength[chromosome]=bases.length();
	}
	return 0;
}

int
Framework::getChromosomeLength(string& key)
{
	if(refchromLength.find(key)==refchromLength.end())
	{
		cout <<"No chromosome by name " << key << endl;
		//exit(0);
	}
	int len=refchromLength[key];
	return len;
}

int
Framework::readCommonNameMap(const char* genenameMap)
{
	gnm.readGeneNames(genenameMap);
	return 0;
}

int
Framework::readMotifNameMap(const char* fName)
{
	mnm.readMap(fName);
	return 0;
}

int 
Framework::mapSitesToGenes()
{
	int windowSize = wsize;
	char pFName[1024];
	sprintf(pFName,"%s_mapped.txt",outname);
	char outFName[1024];
	sprintf(outFName,"%s_edges.txt",outname);

	//ofstream oFile(outFName);
	ofstream pFile(pFName);
	//Will take this as a parameter too.
	//int windowSize=2000;
	for(int s=0;s<siteSet.size();s++)
	{
		Framework::SiteInfo* site=siteSet[s];
		map<string,Gene*>* genesOnChrom=gffreader.getGeneSetForChromosome(site->chrom);
		//cout << "test " << genesOnChrom << "\n";
		if(genesOnChrom==NULL)
		{
		        cout <<"No chromosome " << site->chrom<< " in site "<< site->pos<< endl;
			exit(0);
		}
		Gene* closestGene=NULL;
		int minDist=100000;
		char possibleHitType='\0';
		int position=site->pos;
		int hits=0;
		for(map<string,Gene*>::iterator gIter=genesOnChrom->begin();gIter!=genesOnChrom->end();gIter++)
		{
			Gene* gene=gIter->second;
			int start=gIter->second->geneStart;
			int end=gIter->second->geneEnd;
			bool hitReg=false;
			bool hitCoding=false;
			int minRegDist=0;
			int minCodingDist=0;
			int distBegin=0;
			int distEnd=0;
			int regdistBegin=0;
			int regdistEnd=0;
			char putativeType;
			int codingLen=0;
			int regulatoryLen=0;
			if(gene->strand=='-')
			{
				//Regulatory region is +n bp or -n bp of the ATG
				int beginRegRegion=start+windowSize;
				//int endRegRegion=start-windowSize;
				int endRegRegion=start-windowSize;
				//Now check limits
//				if(gene->downstreamNeighbor==NULL)
//				{
//					//This is the last gene
//					int chromosomeEnd=getChromosomeLength(gene->chromosome);
//					if(beginRegRegion>chromosomeEnd)
//					{	
//						beginRegRegion=chromosomeEnd;
//					}
//				}
//				else
//				{
//					Gene* downstreamNeighbor=gene->downstreamNeighbor;
//					if(beginRegRegion>=downstreamNeighbor->geneEnd)
//					{
//						beginRegRegion=downstreamNeighbor->geneEnd-1;
//					}
//				}
				if(position>=endRegRegion && position<=beginRegRegion)
				{
					hitReg=true;
					regdistBegin=beginRegRegion-position;
					regdistEnd=endRegRegion-position;
				}
				else
				{
					if(position>=end && position<=start)
					{
						hitCoding=true;
						distEnd=position-end;
						distBegin=start-position;
					}
				}
				if(!hitReg && !hitCoding)
				{
					if(position>beginRegRegion)
					{
						minRegDist=position-beginRegRegion;	
						putativeType='R';
					}
					else if(position<end)
					{
						minCodingDist=end-position;	
						putativeType='C';
					}
				}
				codingLen=start-end;
				regulatoryLen=beginRegRegion-endRegRegion;
			}
			else 
			{
				int beginRegRegion=start-windowSize;
				//int endRegRegion=start+windowSize;
				int endRegRegion=start+windowSize;
				//Now check limits
//				if(gene->upstreamNeighbor==NULL)
//				{
//					//This is the last gene
//					if(beginRegRegion<0)
//					{	
//						beginRegRegion=1;
//					}
//				}
//				else
//				{
//					Gene* upstreamNeighbor=gene->upstreamNeighbor;
//					if(beginRegRegion<=upstreamNeighbor->geneEnd)
//					{
//						beginRegRegion=upstreamNeighbor->geneEnd+1;
//					}
//				}

				if(position>=beginRegRegion && position<=endRegRegion)
				{
					hitReg=true;
					regdistBegin=position-beginRegRegion;
					regdistEnd=endRegRegion-position;
				}
				else
				{
					if(position>=start && position<=end)
					{
						hitCoding=true;
						distEnd=end-position;
						distBegin=position-start;
					}
				}
				if(!hitReg && !hitCoding)
				{
					if(position<beginRegRegion)
					{
						minRegDist=beginRegRegion-position;	
						putativeType='R';
					}
					else if(position>end)
					{
						minCodingDist=position-end;	
						putativeType='C';
					}
				}
				codingLen=end-start;
				regulatoryLen=endRegRegion-beginRegRegion;
			}
			//if(hitReg || hitCoding)
			if(hitReg )
			{
				if(hitReg && hitCoding)
				{
					//cout <<"Reg and coding hits for snp "<< s << "!"<< endl;
				}
				map<int,int>* siteForGene=NULL;
				char hittype='C';
				if(hitReg)
				{
					hittype='R';
				}
				if(geneSiteMap.find(gene->name)==geneSiteMap.end())
				{
					siteForGene=new map<int,int>;
					geneSiteMap[gene->name]=siteForGene;
				}
				else
				{
					siteForGene=geneSiteMap[gene->name];
				}
				(*siteForGene)[site->pos]=0;
				hits++;
				//cout <<"Mapping " << site->chrom <<":" << site->pos << " to " << gnm.getCommonName(gene->name.c_str()) << endl;
				//cout <<"Mapping " << site->chrom <<":" << site->pos << " to " << gene->name.c_str() << " " <<  gnm.getCommonName(gene->name.c_str()) << endl;
				//oFile << gnm.getCommonName(gene->name.c_str()) << "\t" << site->name<< endl;
				vector<string>* tfnames = mnm.getTFName(site->name.c_str());
				string cname = gnm.getCommonName(gene->name.c_str());
				for (int titr=0;titr<tfnames->size();titr++)
				{
					elist.addEdge(tfnames->at(titr),cname,site->score);
				}
			}
			else
			{
				int mindistOfgene=minRegDist;
				if(putativeType=='C')
				{
					mindistOfgene=minCodingDist;
				}
				if(minDist>mindistOfgene)
				{
					minDist=mindistOfgene;
					closestGene=gene;
					possibleHitType=putativeType;
				}
			}
		}
		if(hits==0)
		{
			if(closestGene!=NULL)
			{
				//cout <<"No hits found for site " << s<<  " at " << site->pos << " closest_gene " << closestGene->name << "\tDist="  << minDist
				//<< closestGene->strand << "\t" << closestGene->geneStart << "\t" << closestGene->geneEnd<<  "\t" << possibleHitType << endl;
			}
			else
			{
				//cout <<"No hits found for site " <<  s<< " closest gene null " << endl;
			}
		}
		if(hits)
		{
			pFile << site->chrom << "\t" << site->pos-(site->length/2) << "\t" << site->pos+(site->length/2) << endl;
		}
	}
	/*ofstream oFile(outFName);
	for(map<string,map<int,int>*>::iterator gIter=geneSiteMap.begin();gIter!=geneSiteMap.end();gIter++)
	{
		//oFile << gnm.getCommonName(gIter->first.c_str()) << "\t" << gIter->first<<"\t" << gIter->second->size() << endl;
	  oFile << (string)gIter->first.c_str() << "\t"  << gIter->second->size() << endl;
	}*/
	//oFile.close();
	elist.writeToFile(outFName);
	pFile.close();
}

int
Framework::printHelp()
{
	cout << "Mapping motif instances to TSS of targets" << endl;
	cout << "usage:" << endl;
	cout << "\t-s\t\tTSS file" << endl;
	cout << "\t-x\t\tGet TSS per GENE or TRANSCRIPT (default GENE)" << endl;
	cout << "\t-a\t\tFor duplicate edges, get MAX or MIN or SUM or AVG (default MAX)" << endl;
	cout << "\t-m\t\tMotif instances" << endl;
	cout << "\t-c\t\tConfiguration for motif instances (optional)" << endl;
	cout << "\t-r\t\tMapping for chromosome names (optional)" << endl;
	cout << "\t-g\t\tMapping for gene names (optional)" << endl;
	cout << "\t-t\t\tMapping for motif names (optional)" << endl;
	cout << "\t-w\t\tWindow size around TSS (default: 2000)" << endl;
	cout << "\t-o\t\tOutput file prefix" << endl;
	cout << "\t-h\t\tPrint this message" << endl;
	return 0;
}

int
Framework::init(int argc, char** argv)
{
	int c;
	int condCnt = 1;
	//char outFilePrefix[256];
	bool tssDefault=true;//-s
	bool motifDefault=true;//-m
	bool configDefault=true;//-c
	bool chrmapDefault=true;//-r
	bool genemapDefault=true;//-g
	bool motmapDefault=true;//-t
	bool windowDefault=true;//-w
	bool outDefault=true;//-o
	bool addDefault=true;//-o
	bool trnsDefault=true;//-o

	char addname[1024];
	char trnsname[1024];
	//bool mappedregDefault=true;//-
	int optret='-';
	opterr=1;
	int oldoptind=optind;
 
	while(optret=getopt(argc,argv,"s:x:a:m:c:r:g:t:w:o:h")!=-1)
	{
		if(optret=='?')
		{
			cout <<"Option error " << optopt << endl;
			return 0;
		}
		char c;
		char* my_optarg=NULL;
		c=*(argv[oldoptind]+1);
		if(optind-oldoptind ==2)
		{
			my_optarg=argv[oldoptind+1];	
		}
		else
		{
			my_optarg=argv[oldoptind]+2;
		}
		switch(c)
		{
	while(optret=getopt(argc,argv,"s:m:c:r:g:t:w:h")!=-1)
			case 's'://TSS
			{
				tssDefault = false;
				strcpy(tssname,optarg);
				break;
			}
			case 'x'://TSS Type
			{
				trnsDefault = false;
				strcpy(trnsname,optarg);
				break;
			}
			case 'a'://Add Type
			{
				addDefault = false;
				strcpy(addname,optarg);
				break;
			}
			case 'm'://Motif bed
			{
				motifDefault = false;
				strcpy(motname,optarg);
				break;
			}
			case 'c'://Bed config file
			{
				configDefault = false;
				strcpy(cnfname,optarg);
				break;
			}
			case 'r'://chr2chr map
			{
				chrmapDefault = false;
				strcpy(chrmapname,optarg);
				break;
			}
			case 'g'://gene map
			{
				genemapDefault = false;
				strcpy(genemapname,optarg);
				break;
			}
			case 't'://motif tf map
			{
				motmapDefault = false;
				strcpy(motmapname,optarg);
				break;
			}
			case 'w'://window size
			{
				windowDefault = false;
				wsize = atoi(optarg);
				break;
			}
			case 'o'://outname
			{
				outDefault=false;
				strcpy(outname,optarg);
				break;
			}
			case 'h':
			{
				printHelp();
				return 0;
			}
			default:
			{
				cerr << "Unhandled option " << c << endl;
				return 0;
			}
		}
		oldoptind=optind;
	}
	if(tssDefault)
	{
		cerr << "You should provide TSS file" << endl;
		printHelp();
		return 0;
	}
	if(motifDefault)
	{
		cerr << "You should provide Motif instances" << endl;
		printHelp();
		return 0;
	}
	if(outDefault)
	{
		cerr << "You should provide output name" << endl;
		printHelp();
		return 0;
	}
	if(trnsDefault)
	{
		cerr << "Target type was not specified." << endl;
		cerr << "Reading TSS per gene." << endl;
		gffreader = GFFReader(T_GENE);
	}
	else
	{
		if(strcmp(trnsname,"GENE")==0)
		{
			cerr << "Reading TSS per gene." << endl;
			gffreader = GFFReader(T_GENE);
		}
		else if (strcmp(trnsname,"TRANSCRIPT")==0)
		{
			cerr << "Reading TSS per transcript." << endl;
			gffreader = GFFReader(T_TRANS);
		}
		else
		{
			cerr << "Unknown option for target type: " << trnsname << endl;
			cerr << "Reading TSS per gene." << endl;
			gffreader = GFFReader(T_GENE);
		}
	}
	if(addDefault)
	{
		cerr << "The option for handling duplicates is not set." << endl;
		cerr << "Taking the max of duplicate edges." << endl;
		elist = EdgeList(T_MAX);
	}
	else
	{
		if(strcmp(addname,"MAX")==0)
		{
			cerr << "Taking the max of duplicate edges." << endl;
			elist = EdgeList(T_MAX);
		}
		else if(strcmp(addname,"MIN")==0)
		{
			cerr << "Taking the min of duplicate edges." << endl;
			elist = EdgeList(T_MIN);
		}
		else if(strcmp(addname,"SUM")==0)
		{
			cerr << "Taking the sum of duplicate edges." << endl;
			elist = EdgeList(T_SUM);
		}
		else if(strcmp(addname,"AVG")==0)
		{
			cerr << "Taking the average of duplicate edges." << endl;
			elist = EdgeList(T_AVG);
		}
		else
		{
			cerr << "Unknown option for handling duplicates." << endl;
			cerr << "Taking the max of duplicate edges." << endl;
			elist = EdgeList(T_MAX);
		}
	}
	if(configDefault)
	{
		cerr << "Bed config file was not provided!" << endl;
		cerr << "We assume FIMO format:" << endl;
		cerr << "name:0, chr:1, beg:2, end:3, strand:4, score:5" << endl;
		//LM13	chr12	66272519	66272531	+	6.23383	5.1e-06		ATGCTTATGAGCA
		mySiteConf.name_index=0;
		mySiteConf.chr_index=1;
		mySiteConf.beg_index=2;
		mySiteConf.end_index=3;
		mySiteConf.strand_index=4;
		mySiteConf.score_index=5;
	}
	else
	{
		readBedConfig(cnfname);
	}
	if(chrmapDefault)
	{
		cerr << "Chromosome map was not provided!" << endl;
		cerr << "We keep motif chromosome names as are." << endl;
	}
	else
	{
		readChromMap(motmapname);
	}
	if(genemapDefault)
	{
		cerr << "Gene map was not provided!" << endl;
		cerr << "We keep TSS gene names as they are." << endl;
	}
	else
	{
		readCommonNameMap(genemapname);
	}
	if(motmapDefault)
	{
		cerr << "Motif map was not provided!" << endl;
		cerr << "We keep motif instances names as they are." << endl;
	}
	else
	{
		readMotifNameMap(motmapname);
	}
	if(windowDefault)
	{
		cerr << "Window size was not provided!" << endl;
		cerr << "Set to 2000." << endl;
		wsize = 2000;
	}
	return 1;
}

int
main(int argc, char** argv)
{
	//if(argc!=9)
	//{
	//	cout <<"Usage: ./mapSitesToGenes gfffile chromnamemap_chip sites single|double genenamemap output windowSize pFName"<< endl;
	//	return 0;
	//}
	Framework fw;
	if(fw.init(argc,argv))
	{
		fw.readGFFFile();
		fw.readSites();
		fw.mapSitesToGenes();
	}
	return 0;
}
