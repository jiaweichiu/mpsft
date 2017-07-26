#ifndef TIMER_H
#define TIMER_H

//#define TIMER_TYPE CLOCK_REALTIME
#define TIMER_TYPE CLOCK_PROCESS_CPUTIME_ID

void reset_timer();
double get_time();

#endif
