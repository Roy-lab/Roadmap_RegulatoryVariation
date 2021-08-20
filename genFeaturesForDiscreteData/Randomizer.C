#include <iostream>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/timeb.h>
#include <time.h>

#include "Randomizer.H"

Randomizer::Randomizer()
{
}

Randomizer::~Randomizer()
{
}

void
Randomizer::initialize(int anUpper, int aLower)
{
	upperBound=anUpper;
	lowerBound=aLower;
	ftime(&atime);
	int pid=getpid();
	 
	srand((unsigned int) atime.millitm);
	//srand(pid);
}

int
Randomizer::getRandomNumber()
{
	int randomNumCnt=upperBound-lowerBound;
//	sleep(1);
	int randNum=rand();
	int randIndex=randNum%(randomNumCnt-1);
	return randIndex;
	
}

#ifdef TEST

int main()
{
	Randomizer randomize;
	randomize.initialize(0,6000);
	for(int i=0;i <100; i++)
	{
		int no=randomize.getRandomNumber();
		cout << no  << endl;	
	}
	return 0;
}
#endif
