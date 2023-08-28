#include "tictoc.h"
#include <iostream>

TicTocClock::TicTocClock()
{
	time(&starttime);
}

TicTocClock::~TicTocClock()
{

}

void TicTocClock::toc(std::string message)
{
#ifdef OXLEY_ENABLE_PROFILE_TIMERS
	time_t currenttime;
	time(&currenttime);

	std::cout << "\033[1;31m[" << difftime(currenttime,starttime) << " seconds]\033[0m: " << message << std::endl;
#endif
}

void TicTocClock::setTime(time_t newtime)
{
#ifdef OXLEY_ENABLE_PROFILE_TIMERS
	starttime=newtime;
#endif
}

time_t TicTocClock::getTime()
{
#ifdef OXLEY_ENABLE_PROFILE_TIMERS
	return starttime;
#endif
}

incrementClock::incrementClock(int n)
{
#ifdef OXLEY_ENABLE_PROFILE_INCREMENT_TIMERS
	time(&starttime);
	time(&currenttime);
	time(&lasttime);
	numberOfIncrements = n;
	timer.resize(n);
#endif
}

incrementClock::~incrementClock()
{

}

void incrementClock::toc(int n)
{
#ifdef OXLEY_ENABLE_PROFILE_INCREMENT_TIMERS
	if(n >= numberOfIncrements)
		std::cout << "n is out of bounds " << std::endl;

	time(&currenttime);
	timer[n-1]+=difftime(currenttime,lasttime);
	time(&lasttime);
#endif
}

void incrementClock::report()
{
#ifdef OXLEY_ENABLE_PROFILE_INCREMENT_TIMERS
	for(int i = 0; i < numberOfIncrements; i++)
	{
		std::cout << timer[i] << " -- ";
	}
	std::cout << std::endl;
#endif
}
