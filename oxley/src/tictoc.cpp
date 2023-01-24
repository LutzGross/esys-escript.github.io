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