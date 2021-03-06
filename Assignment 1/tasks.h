/***************************************************************************
 *
 *   File        : tasks.h
 *   Student Id  : 756560
 *   Name        : Jiawei Liao
 *
 ***************************************************************************/

#ifndef TASKS_H
#include <stdio.h>


/*data type to store the information of single flowpoint*/
typedef struct
{
  	float x;
  	float y;
  	float u;
  	float v;
}	flowpoint;


/*data type to store the information of a cell in Task2*/
typedef struct
{
    /*number of points in the cell*/
    int num_points;
    /*average values of attributes of points in cell*/
  	float avx;
    float avy;
    float avu;
    float avv;
    float s;
    /*boundaries of the cell*/
    double xlow;
    double xhigh;
    double ylow;
    double yhigh;
}	cell;


/*data type to store the information of a wakepoint in Task4*/
typedef struct
{
    float x;
    float y;
    float y_abs; /*magnitude of y*/
    float u;
    float x_exact;  /*the exact x value expected*/
    float x_offset; /*the difference between the actual x and the x expected*/
    int isempty;  /*whether a flowpoint is assigned to the wakepoint*/
} wakepoint;


/*taks functions prototypes*/
void maxveldiff(const char* flow_file);

void coarsegrid(const char* flow_file, int resolution);

void velstat(const char* flow_file);

void wakevis(const char* flow_file);


/*subfunction prototypes*/



/*READ IN SPECIFIC*/
/*----------------------------------------------------------------------------------*/

/*open a file with exception handelling, taken from L2/file_io.c*/
FILE* safe_fopen(const char* path, const char* mode);
/*malloc with expception handelling, taken from L2/file_io.c*/
void* safe_malloc(size_t num_bytes);
/*check the header of the file is present.*/
void check_header(FILE* fp);

/*----------------------------------------------------------------------------------*/




/*TASK 1 SPECIFIC*/
/*----------------------------------------------------------------------------------*/

/*assign the pointa with the attributes of point b*/
void assign_point(flowpoint *pointa, flowpoint *pointb);
/*writes the data of a single point*/
void print_point(FILE* fp,flowpoint* point);

/*----------------------------------------------------------------------------------*/




/*TASK 2 SPECIFIC*/
/*----------------------------------------------------------------------------------*/

/*create a new cell result for task2, index is the index of the cell in the
cells data, x_grid_size and y_grid_size are the dimensions ot the cell*/
void new_cell(cell* newcell, int index,
                      int resolution, double x_grid_size, double y_grid_size);
/*add a point's values to a cell in task2*/
void add_point(cell* cell, flowpoint* point);
/*calculate the values in the cell in task2*/
void calculate_point(cell* cell);
/*writes the data of a single cell*/
void print_cell(FILE* fp,cell* cell);
/*a compare function for s for qsort, descending order*
returns 1 if a<b, -1 if a>b ,0 if a=b*/
int compare_s(const void * a, const void * b);

/*----------------------------------------------------------------------------------*/




/*TASK 3 SPECIFIC*/
/*----------------------------------------------------------------------------------*/

/*a compare function for magnitude of u, ascending order
returns 1 if a>b, -1 if b<a, 0 if a=b*/
int compare_u_abs(const void * a, const void * b);

/*----------------------------------------------------------------------------------*/




/*TASK 4 SPECIFIC*/
/*----------------------------------------------------------------------------------*/

/*writes the data of a single wakepoint*/
void print_wakepoint(FILE* fp,wakepoint* wakepoint);
/*assign the relevant attributes of a flowpoint to a wakepoint*/
void assign_wakepoint(wakepoint* wake, flowpoint* flow);

/*----------------------------------------------------------------------------------*/




/*OUTPUT SPECIFIC*/
/*----------------------------------------------------------------------------------*/

/*writes the header of each tasks*/
void print_header(FILE* fp, char* header);

/*----------------------------------------------------------------------------------*/


#endif
