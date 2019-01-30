#include "common.h"
#include "semel.h"

#ifdef __APPLE__
#include <mach/clock.h>
#include <mach/mach.h>
#else
#include <time.h>
#endif

void semel_get_monotonic_time(
    struct timespec * ts)
{
#ifdef __APPLE__
  clock_serv_t cclock;
  mach_timespec_t mts;
  host_get_clock_service(mach_host_self(), SYSTEM_CLOCK, &cclock);
  clock_get_time(cclock, &mts);
  mach_port_deallocate(mach_task_self(), cclock);
  ts->tv_sec = mts.tv_sec;
  ts->tv_nsec = mts.tv_nsec;
#else
  //clock_gettime(CLOCK_MONOTONIC, ts);
	timespec_get(ts, TIME_UTC);
#endif
}

uint64_t semel_lap(
		void)
{
	uint64_t nsec = 0;
	struct timespec ts;
	semel_get_monotonic_time(&ts);
	nsec = (uint64_t) ((C_1_BLN * ts.tv_sec) + ts.tv_nsec);
	return nsec;
}
