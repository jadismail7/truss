
#ifndef _TIMER_H_
#define _TIMER_H_

#include <stdio.h>
#include <time.h>

typedef struct
{
        clock_t startTime;
        clock_t endTime;
} Timer;

static void startTime(Timer *timer)
{
        timer->startTime = clock();
}

static void stopTime(Timer *timer)
{
        timer->endTime = clock();
}

static void printElapsedTime(Timer timer, const char *s)
{
        float t = (float)(timer.endTime - timer.startTime) / CLOCKS_PER_SEC;
        printf("%s: %f s\n", s, t);
}

#endif
