/***************************************************************************
 *
 *   File        : main.c
 *   Student Id  : 756560
 *   Name        : Jiawei Liao
 *
 ***************************************************************************/



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <string.h>
#include "tasks.h"
#include <assert.h>



void print_tasktiming(struct timeval start, struct timeval stop, int task_num);

int main(int argc, char *argv[]) {

	/* TODO: Parse Command Line Arguments */
	char* flow_file = NULL;
	int resolution = 0;

	/* TODO: Add timing for each task and output running time in ms */

	/*structs to store the start and finish time of the program*/
	struct timeval start;
	struct timeval stop;
	/*begin timing*/
	gettimeofday(&start, NULL);

	/*check if input has sufficient number of components*/
	if (argc != 3) {
       printf("USAGE: %s <in file name> <resolution>\n", argv[0]);
       exit(EXIT_FAILURE);
  }
	/*the flow_data file input*/
	flow_file = argv[1];
	/*the resolution of task2*/
	resolution = atoi(argv[2]);

	/* Task 1: Find the maximum velocity difference */
	maxveldiff(flow_file);

	gettimeofday(&stop, NULL);
	print_tasktiming(start, stop, 1);

	gettimeofday(&start, NULL);
	/* Task 2: Coarser Grid */
	coarsegrid(flow_file, resolution);

	gettimeofday(&stop, NULL);
	print_tasktiming(start, stop, 2);

	gettimeofday(&start, NULL);
	/* Task 3: Statistics */
	velstat(flow_file);

	gettimeofday(&stop, NULL);
	print_tasktiming(start, stop, 3);

	gettimeofday(&start, NULL);
	/* Task 4: Wake height and visualisation */
	wakevis(flow_file);

	gettimeofday(&stop, NULL);
	print_tasktiming(start, stop, 4);

	return (EXIT_SUCCESS);
}


void print_tasktiming(struct timeval start, struct timeval stop, int task_num)
{
		double elapsed_ms = (stop.tv_sec - start.tv_sec) * 1000.0;
		elapsed_ms += (stop.tv_usec - start.tv_usec) / 1000.0;
		printf("TASK %d:  %.2f milliseconds\n",task_num, elapsed_ms);
}
