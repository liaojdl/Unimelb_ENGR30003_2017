/***************************************************************************
 *
 *   File        : tasks.c
 *   Student Id  : 756560
 *   Name        : Jiawei liao
 *
 ***************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <string.h>
#include <assert.h>
#include "tasks.h"

/*constants*/
/*file names*/
#define TASK1HEADER "x,y,u,v" /*the header of input file and task 1 output*/
#define TASK2HEADER "x,y,u,v,S"/*the header of input file and task 2 output*/
#define TASK3HEADER "threshold,points,percentage" /*task 3 output header*/
#define TASK4HEADER "x,y_h" /*task 4 output header*/
#define TASK1FILE "task1.csv" /*file name for task1 output*/
#define TASK2FILE "task2.csv" /*file name for task1 output*/
#define TASK3FILE "task3.csv" /*file name for task3 output*/
#define TASK4FILE "task4_1.csv" /*file name for task4 output*/

/*constans numbers*/
#define INIT_FLOWDATASIZE 1 /*default number of the flowpoints expected*/
#define TASK1_XMIN 20.000000 /*minimum value of x for task1*/

/*domain of task 2 analysis*/
#define X_BOUND_L 10.000000
#define X_BOUND_H 70.000000
#define Y_BOUND_L -20.000000
#define Y_BOUND_H 20.000000

/*task 3 threshold and threshold increments*/
#define TASK3_THRESHOLD_MIN 0.500000
#define TASK3_THRESHOLD_INC 0.100000

/*task 4 threshold and threshold increments*/
#define TASK4_WAKE_INC 5.000000
#define TASK4_X_BOUND 0.050000

/*boolean types casted as int, to determine whether a point has been processed
for certain tasks*/
#define HASPOINT 1
#define NOPOINT -1






/*implements task1, find the pair of points with the maximum difference
in u velocity and the pair of points with the maximum difference in v
velocity and ouputs the four points into the file task1.csv*/
void maxveldiff(const char* flow_file)
{

    /*temporaily stores the data of one point*/
    flowpoint* next = (flowpoint*)safe_malloc(sizeof(flowpoint));
    /*stores the point with maximum u velocity*/
    flowpoint* u_max = (flowpoint*)safe_malloc(sizeof(flowpoint));
    /*stores the point with minium u velocity*/
    flowpoint* u_min = (flowpoint*)safe_malloc(sizeof(flowpoint));
    /*stores the point with maximum v velocity*/
    flowpoint* v_max = (flowpoint*)safe_malloc(sizeof(flowpoint));
    /*stores the point with minium v velocity*/
    flowpoint* v_min = (flowpoint*)safe_malloc(sizeof(flowpoint));
    /*whether anypoint is being read at all, preset as not*/
    int has_points = NOPOINT;


    /*open file for input*/
    FILE* fp = safe_fopen(flow_file, "r");
    /*check the header is before the data points*/
    check_header(fp);
    /*take inputs*/
    while (fscanf(fp, "%f,%f,%f,%f", &next->x, &next->y,&next->u,&next->v) == 4)
    {

        /*only consider points where x>20*/
        if (next->x > TASK1_XMIN)
        {
            /*presets u_max,u_min,v_max,v_min as the first point read in*/
            if (has_points == NOPOINT)
            {
              assign_point(u_max,next);
              assign_point(u_min,next);
              assign_point(v_max,next);
              assign_point(v_min,next);
              has_points = HASPOINT;
            /*update u_max,u_min,v_max,v_min accordingly*/
            }else {
              if(next->u > u_max->u) {assign_point(u_max,next);}
              if(next->u < u_min->u) {assign_point(u_min,next);}
              if(next->v > v_max->v) {assign_point(v_max,next);}
              if(next->v < v_min->v) {assign_point(v_min,next);}
            }
        }
    }

    fclose(fp); /*close input file*/

    /*outputing*/
    FILE* f1 = safe_fopen(TASK1FILE, "w");
    print_header(f1, TASK1HEADER);
    print_point(f1,u_max);
    print_point(f1,u_min);
    print_point(f1,v_max);
    print_point(f1,v_min);

    /*cleanup*/
    free(next);
    free(u_max);
    free(u_min);
    free(v_max);
    free(v_min);
    fclose(f1);/*close output file*/
}


/*implements task2*/
void coarsegrid(const char* flow_file, int resolution)
{
    /*the step size of x_grid and y_grid for the cells*/
    double x_grid_size = (X_BOUND_H-X_BOUND_L)/resolution;
    double y_grid_size = (Y_BOUND_H-Y_BOUND_L)/resolution;

    /*intialize resolution^2 number of cells, with their boundaries
    calculated and stored, all other values zeroed*/
    cell* cells = (cell*)safe_malloc(pow(resolution,2)*sizeof(cell));
    for (int i = 0; i < pow(resolution,2); i++)
    {
      new_cell(&cells[i], i, resolution, x_grid_size, y_grid_size);
    }

    /*temporaily stores the data of one point*/
    flowpoint* next = (flowpoint*)safe_malloc(sizeof(flowpoint));

    /*the x_grid index*/
    int x_grid;
    /*the x_grid index*/
    int y_grid;
    /*the actual index in the 1d array,
    represented as x_grid*resolution+y_grid*/
    int grid_index;


    /*open file for input*/
    FILE* fp = safe_fopen(flow_file, "r");
    /*check the header is before the data points*/
    check_header(fp);
    /*inputs*/
    while (fscanf(fp, "%f,%f,%f,%f", &next->x, &next->y,&next->u,&next->v) == 4)
    {
        /*check in bound of domain*/
        if((next->x >= X_BOUND_L) && (next->x <= X_BOUND_H) &&
            (next->y >= Y_BOUND_L) && (next->y <= Y_BOUND_H))
        {
            /*calculate x degree index*/
            x_grid = (int)((next->x - (double)X_BOUND_L)/x_grid_size);
            /*change x from resolution to resolution-1 in case of error*/
            if (x_grid == resolution)
            {
              x_grid --;
            }
            /*calculate y_degree index*/
            y_grid = (int)((next->y - (double)Y_BOUND_L)/y_grid_size);
            /*change y from resolution to resolution-1 in case of error*/
            if (y_grid == resolution)
            {
              y_grid --;
            }
            /*calculate actual 1d index*/
            grid_index = x_grid*resolution + y_grid;

            /*check for cell[grid_index] against the point*/
            if ((next->x >= cells[grid_index].xlow) &&
                    (next->x <= cells[grid_index].xhigh) &&
                          (next->y >= cells[grid_index].ylow) &&
                                    (next->y <= cells[grid_index].yhigh))
            {
                /*add the point if in bound*/
                add_point(&cells[grid_index],next);

                /*check neighbour cells*/

                /*check left boundary for neighbour cell*/
                if (next->x == cells[grid_index].xlow && (x_grid > 0))
                {
                  add_point(&cells[grid_index-1],next);
                }
                /*check right boundary for neighbour cell*/
                if (next->x == cells[grid_index].xhigh &&
                                                (x_grid < (resolution-1)))
                {
                  add_point(&cells[grid_index+1],next);
                }
                /*check bottom boundary for neighbour cell*/
                if (next->y == cells[grid_index].ylow && (y_grid > 0))
                {
                  add_point(&cells[(y_grid-1)*resolution+x_grid],next);
                }
                /*check top boundary for neighbour cell*/
                if (next->y == cells[grid_index].yhigh &&
                                                (y_grid < (resolution-1)))
                {
                  add_point(&cells[(y_grid+1)*resolution+x_grid],next);
                }
            }
        }
    }
    /*close input and cleanup of next*/
    fclose(fp);
    free(next);


    /*calculate the required average values and S of each cell*/
    for (int i = 0; i < pow(resolution,2); i++)
    {
      calculate_point(&cells[i]);
    }

    /*sort with respect to S, in descending order*/
    qsort(cells, pow(resolution,2), sizeof(cell), compare_s);


    /*outputing*/
    FILE* f2 = safe_fopen(TASK2FILE, "w");
    print_header(f2, TASK1HEADER);
    for (int i = 0; i < pow(resolution,2); i++)
    {
      print_cell(f2,&cells[i]);
    }

    /*close ouput and cleanup cells*/
    fclose(f2);
    free(cells);


    //exit(EXIT_FAILURE);
}



/*implements Task3, calculate the number and the percentage of points
below each threshold of magnitude of u velocity, starting from 0.5,
and increases by 0.1 everytime,until all the points are included*/
void velstat(const char* flow_file)
{

    /*intialize the dynamic array of only size 1*/
    int fsize = INIT_FLOWDATASIZE;
    /*the malloc array to store all the points in the data*/
    flowpoint* flowdata = (flowpoint*)safe_malloc(fsize*sizeof(flowpoint));
    /*temporaily stores the data of one point*/
    flowpoint* next = (flowpoint*)safe_malloc(sizeof(flowpoint));
    /*keeps track of the points read in total*/
    int points_read = 0;



    /*open file for input*/
    FILE* fp = safe_fopen(flow_file, "r");
    /*check the header is before the data points*/
    check_header(fp);

    /*take inputs*/
    while (fscanf(fp, "%f,%f,%f,%f", &next->x, &next->y,&next->u,&next->v) == 4)
    {
        /*double the flowdata size if required, max log(n) times*/
        if (points_read >= fsize)
        {
          fsize *= 2;
          flowdata = (flowpoint*)realloc(flowdata, fsize*sizeof(flowpoint));
          assert (flowdata != NULL);
        }
        /*check within domain*/
        if ((next->x>=X_BOUND_L)&&(next->x<=X_BOUND_H)&&
              (next->y>=Y_BOUND_L)&&(next->y<=Y_BOUND_H))
        {
          /*add the point into flowdata*/
          assign_point(&flowdata[points_read],next);
          points_read++;
        }
    }
    /*close input and clean up next*/
    fclose(fp);
    free(next);

    /*quicksort the flowdata in terms of magnitude of u velocity in ascending
    order*/
    qsort (flowdata, points_read, sizeof(flowpoint), compare_u_abs);


    /*outputing*/
    FILE* f3 = safe_fopen(TASK3FILE, "w");
    print_header(f3, TASK3HEADER);

    /*tracks the number of points considered*/
    int points_count = 0;
    /*records the threshold*/
    double threshold = TASK3_THRESHOLD_MIN;
    /*records the percentage below each threshold, initiiaze as zero*/
    float percentage = 0;
    /*track the number of thresholds used,initiiaze as zero*/
    int iterations = 0;
    /*dont stop till all points  are considered*/
    while (points_count < points_read)
    {
        /*stop if the points's u magnitude exceed the current threshold
        and when all points are accounted for*/
        while (points_count < points_read &&
                fabsf(flowdata[points_count].u) < threshold)
        {
            /*adds one point*/
            points_count ++;
        }
        /*calculate the percentage of points below each threshold
        the points of the last threshold will not be counted again*/
        percentage = 100*(float)points_count/points_read;
        if (fprintf(f3,"%.6f,%d,%.6f\n",threshold,points_count,percentage) < 0)
        {
            perror("file write error");
        }
        /*increase one more iteration*/
        iterations ++;
        /*change the threshold according to the number of iterations,
        instead of using ++ everytime to avoid added error*/
        threshold = TASK3_THRESHOLD_INC*iterations + TASK3_THRESHOLD_MIN;
    }

    /*close ouput and cleanup flowdata*/
    fclose(f3);
    free(flowdata);
    //exit(EXIT_FAILURE);
}



/*implements Task4, generate the wake_data at locations 10, 15, 20 etc all the
way to 65, and prints out both in csv format and the part 2 prints out a txt
visualization file*/
void wakevis(const char* flow_file)
{
    /*part 1*/

    /*the number of wake_data points, which should be 12*/
    int num_wakes = (int)(X_BOUND_H-X_BOUND_L)/TASK4_WAKE_INC;
    /*temporaily stores the data of one point*/
    flowpoint* next = (flowpoint*)safe_malloc(sizeof(flowpoint));
    /*the malloc array to store the wake_points data*/
    wakepoint* wakedata = (wakepoint*)safe_malloc(num_wakes*sizeof(wakepoint));
    /*preset the wake_points in wakedata*/
    for (int i = 0; i < num_wakes; i++)
    {
        /*preset as no points processed*/
        wakedata[i].isempty = NOPOINT;
        /*the attencipated x value for each wake_point, used
        to compare against the points read*/
        wakedata[i].x_exact = X_BOUND_L + i*TASK4_WAKE_INC;
    }


    /*open the input file*/
    FILE* fp = safe_fopen(flow_file, "r");
    /*check the header is before the data points*/
    check_header(fp);
    /*tracks the index wakepoint in wake_data for comparing against*/
    int wakeindex;
    /*tracks the point's x value minus the inital x value, 10*/
    float x_abs;
    /*the dx for checking which x-value to check closely*/
    float dx = TASK4_WAKE_INC/2;
    /*abosulute difference between a point's x value and the prefered*/
    float difference;
    /*take in points*/
    while (fscanf(fp, "%f,%f,%f,%f", &next->x, &next->y,&next->u,&next->v) == 4)
    {
        /*calculates the x- X_BOUND_L value*/
        x_abs = (next->x)-X_BOUND_L;
        /*work out which x-value in wake_profile this point is closest to
        for close checking*/
        wakeindex = (int)(x_abs+dx/2)/TASK4_WAKE_INC;

        /*process a point each time ,and a point is processed only once*/
        if (wakeindex < num_wakes)
        {
            difference = fabsf(next->x - wakedata[wakeindex].x_exact);
            /*check if this point is whithin lower and upperbound*/
            if (difference <= TASK4_X_BOUND)
            {
                /*preset the wake_point as the first point with x-value within
                bound*/
                if (wakedata[wakeindex].isempty == NOPOINT)
                {
                  wakedata[wakeindex].x_offset = difference;
                  assign_wakepoint(&wakedata[wakeindex],next);
                  wakedata[wakeindex].isempty = HASPOINT;
                }
                /*each wake point is preseted already*/
                else
                {
                    /*if a point's x value is closer to the x_exact value for
                    this wake_point, replace the original*/
                    if (difference < wakedata[wakeindex].x_offset)
                    {
                        wakedata[wakeindex].x_offset = difference;
                        assign_wakepoint(&wakedata[wakeindex],next);
                    }
                    /*compare u value
                    if a point has exact x value as the previous*/
                    else if (next->x == wakedata[wakeindex].x)
                    {
                        /*find the maximum u-velocity point*/
                        if ((next->u) > (wakedata[wakeindex].u))
                        {
                            wakedata[wakeindex].u = next->u;
                            wakedata[wakeindex].y_abs = fabsf(next->y);
                        }
                        /*if a few points have same u velcity, choose
                        the one with the smallest y value*/
                        else if ((next->u) == (wakedata[wakeindex].u) &&
                                (next->y) < (wakedata[wakeindex].y))
                        {
                            wakedata[wakeindex].y_abs = fabsf(next->y);
                        }
                    }
                }
            }
        }
    }
    /*close input and cleanup next*/
    fclose(fp);
    free(next);


    /*outputing*/
    FILE* f4 = safe_fopen(TASK4FILE, "w");
    print_header(f4, TASK4HEADER);
    for (int i = 0; i < num_wakes; i++) {
        print_wakepoint(f4, &wakedata[i]);
    }
    /*close output*/
    fclose(f4);


    /*Task4 part 2, visualization of wakedata*/

    int i,j;
    int n = 12; // Location in x for wake visualization
    float* yheight;
    yheight = (float*) calloc(n,sizeof(float));

    /* Task 4: Part 2, nothing is to be changed here
       Remember to output the spacing into the array yheight
       for this to work. You also need to initializ?
       e i,j and
       yheight so the skeleton as it stands will not compile */
    /*calculate the spacing for part 2*/
    for (i = 0; i < n; i++) {
        yheight[i] = (float)ceil(10 * wakedata[i].y_abs);
    }

    FILE *ft42;
    ft42 = fopen("task4_2.txt","w");
    for (j = 11; j>=0; j--){
	     for (i=0;i<yheight[j]-yheight[0]+4;i++){
 	       fprintf(ft42, " ");
	      }
    	fprintf(ft42, "*\n");
    }
    for (i=0;i<5; i++){
    	fprintf(ft42, "III\n");
    }
    for(j = 0; j<12; j++ ){
    	for (i=0;i<yheight[j]-yheight[0]+4;i++){
    	    fprintf(ft42, " ");
    	}
    	fprintf(ft42, "*\n");
    }
    fclose(ft42);

    /* Cleanup */

    free(yheight);
    free(wakedata);

    //exit(EXIT_FAILURE);
}


/*assign the pointa with the attributes of point b*/
void assign_point(flowpoint *pointa, flowpoint *pointb) {
    pointa->x = pointb->x;
    pointa->y = pointb->y;
    pointa->u = pointb->u;
    pointa->v = pointb->v;
}


/*writes the data of a single point*/
void print_point(FILE* fp,flowpoint* point)
{
    if (fprintf(fp,"%.6f,%.6f,%.6f,%.6f\n",
    point->x, point->y,point->u,point->v) < 0)
    {
      perror("file write error");
    }
}


/*create a new cell result for task2, index is the index of the cell in the
cells data, x_grid_size and y_grid_size are the dimensions ot the cell*/
void new_cell(cell* newcell, int index, int resolution,
                  double x_grid_size, double y_grid_size)
{
    newcell->num_points = 0;
    newcell->avx = 0;
    newcell->avy = 0;
    newcell->avu = 0;
		newcell->avv = 0;
    newcell->s = 0;
    newcell->xlow = X_BOUND_L + (index/resolution)*x_grid_size;
    newcell->xhigh = X_BOUND_L + (index/resolution+1)*x_grid_size;
    newcell->ylow = Y_BOUND_L + (index%resolution)*y_grid_size;
    newcell->yhigh = Y_BOUND_L + (index%resolution+1)*y_grid_size;
}


/*add a point's values to a cell in task2*/
void add_point(cell* cell, flowpoint* point)
{
    cell->avx += point->x;
    cell->avy += point->y;
    cell->avu += point->u;
    cell->avv += point->v;
    cell->num_points ++;
}


/*calculate the values in the cell in task2*/
void calculate_point(cell* cell)
{
    cell->avx /= cell->num_points;
    cell->avy /= cell->num_points;
    cell->avu /= cell->num_points;
    cell->avv /= cell->num_points;
    cell->s = 100*sqrt(pow(cell->avu,2) + pow(cell->avv,2))/
                  sqrt(pow(cell->avx,2) + pow(cell->avy,2));
}

/*assign the relevant attributes of a flowpoint to a wakepoint*/
void assign_wakepoint(wakepoint* wake, flowpoint* flow)
{
  wake->x = flow->x;
  wake->y = flow->y;
  wake->u = flow->u;
  wake->y_abs = fabsf(flow->y);
}

/*open a file with exception handelling, taken from L2/file_io.c*/
FILE* safe_fopen(const char* path, const char* mode)
{
    FILE* fp = fopen(path, mode);
    if (fp == NULL) {
        perror("file open error.");
        exit(EXIT_FAILURE);
    }
    return fp;
}

/*malloc with expception handelling, taken from L2/file_io.c*/
void* safe_malloc(size_t num_bytes)
{
		void* ptr = malloc(num_bytes);
		if (ptr == NULL) {
			printf("ERROR: malloc(%lu)",num_bytes);
			exit(EXIT_FAILURE);
		}
		return ptr;
}

/*check the header of the file is present.*/
void check_header(FILE* fp)
{
    char x,y,u,v;
    if (fscanf(fp, "%c,%c,%c,%c", &x, &y,&u,&v) != 4) {
      perror("Error confirming header\n");
    }
}



/*writes the data of a single cell*/
void print_cell(FILE* fp,cell* cell)
{
    if (fprintf(fp,"%.6f,%.6f,%.6f,%.6f,%.6f\n",
    cell->avx, cell->avy, cell->avu, cell->avv, cell->s) < 0)
    {
      perror("file write error");
    }
}

/*writes the data of a single wakepoint*/
void print_wakepoint(FILE* fp,wakepoint* wakepoint)
{
    if (fprintf(fp,"%.6f,%.6f\n",
    wakepoint->x, wakepoint->y_abs) < 0)
    {
      perror("file write error");
    }
}

/*writes the header of each tasks*/
void print_header(FILE* fp, char* header)
{
    if (fprintf(fp,"%s\n", header) < 0)
    {
      perror("file write error");
    }
}

/*a compare function for magnitude of u, ascending order
returns 1 if a>b, -1 if b<a, 0 if a=b*/
int compare_u_abs(const void * a, const void * b)
{
  flowpoint* pointa = (flowpoint*)a;
  flowpoint* pointb = (flowpoint*)b;
  int agreater = (fabsf(pointa->u) > fabsf(pointb->u));
  int bgreater = (fabsf(pointa->u) < fabsf(pointb->u));
  return (agreater-bgreater);
}

/*a compare function for s for qsort, descending order*
returns 1 if a<b, -1 if a>b ,0 if a=b*/
int compare_s(const void * a, const void * b)
{
  cell* cella = (cell*)a;
  cell* cellb = (cell*)b;
  int agreater = cella->s > cellb->s;
  int bgreater = cella->s < cellb->s;
  return (bgreater - agreater);
}
