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

	std::cout << "[" << difftime(currenttime,starttime) << " seconds]: " << message << std::endl;
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