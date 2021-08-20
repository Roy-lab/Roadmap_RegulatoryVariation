#ifndef _RANDOMIZER
#define _RANDOMIZER

#include <unistd.h>
#include <sys/types.h>
#include <sys/timeb.h>
#include <time.h>

class Randomizer
{
	public:
		Randomizer();
		~Randomizer();
		void initialize(int, int);
		int getRandomNumber();

	private:

		int lowerBound;
		int upperBound;
		struct timeb atime;
};
#endif
