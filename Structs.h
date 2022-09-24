

#ifndef STRUCTS_H
#define STRUCTS_H



struct point2i //(or pixel).
{
	int x,y;
};

struct point2d
{
	double x,y;
};

struct point1d1i
{
	double data;
	int cnt;
};

struct point3d
{
	double x,y;
	double r;
};

struct point3i
{
	int x,y;
	int z;
};

struct point2d1i
{
	double x,y;
	int z;
};

struct  point5d
{
	double x,y;
	double a,b;
	double phi;
};

/*----------------------------------------------------------------------------*/
/** Rectangle structure: line segment with width.
 */
struct rect
{
  double x1,y1,x2,y2;  /* first and second point of the line segment */
  double width;        /* rectangle width */
  double x,y;          /* center of the rectangle */
  double theta;        /* angle */
  double dx,dy;        /* (dx,dy) is vector oriented as the line segment,dx = cos(theta), dy = sin(theta) */
  int   polarity;     /* if the arc direction is the same as the edge direction, polarity = 1, else if opposite ,polarity = -1.*/
  double prec;         /* tolerance angle */
  double p;            /* probability of a point with angle within 'prec' */
};

typedef struct
{
  double vx[4];  /* rectangle's corner X coordinates in circular order */
  double vy[4];  /* rectangle's corner Y coordinates in circular order */
  double ys,ye;  /* start and end Y values of current 'column' */
  int x,y;       /* coordinates of currently explored pixel */
} rect_iter;

typedef struct image_double_s
{
  double * data;
  int xsize,ysize;
} * image_double;

/*----------------------------------------------------------------------------*/
/** Chained list of coordinates.
 */
struct coorlist
{
  int x,y;
  struct coorlist * next;
};
typedef struct ntuple_list_s
{
  int size;
  int max_size;
  int dim;
  double * values;
} * ntuple_list;



//================================Generate Ellipse Candidates=========================================
//匹配组对，组对的索引参数，椭圆参数
// (match group pair, index parameter of group pair, ellipse parameter)
typedef struct PairGroup_s
{
	point2i pairGroupInd;
	point2d center;  //(x0,y0)
	point2d axis;    //(a,b)
	double  phi;     //angle of orientation  
}PairGroup;

//匹配组对节点 (Match group pair nodes)
typedef struct PairGroupNode_s
{
	point2i pairGroupInd;
	point2d center;  //(x0,y0)
	point2d axis;    //(a,b)
	double  phi;     //angle of orientation  
	PairGroupNode_s* next;
}PairGroupNode;

typedef struct  PairGroupList_s
{
	int length;
	PairGroup * pairGroup;
}PairGroupList;

typedef struct Point2dNode_s
{
	point2d point;
	Point2dNode_s * next;
}Point2dNode;

typedef struct Point3dNode_s
{
	point3d point;
	Point3dNode_s * next;
}Point3dNode;

typedef struct Point5dNode_s
{
	point2d center;
	point2d axis;
	double  phi;
	Point5dNode_s * next;
}Point5dNode;

typedef struct Point1dNode_s
{
	double data;
	Point1dNode_s * next;
}Point1dNode;

PairGroupList * pairGroupListInit( int length)
{
	if(length <= 0)
		error("paired groups length less equal than 0");
	PairGroupList * pairGroupList = (PairGroupList*)malloc(sizeof(PairGroupList));
	pairGroupList->length = length;
	pairGroupList->pairGroup = (PairGroup*)malloc(sizeof(PairGroup)*length);
	if(pairGroupList->pairGroup == NULL)
		error("pairGroupListInit,not enough memory");
	return pairGroupList;
}

void freePairGroupList( PairGroupList * list)
{
	if(list == NULL || list->pairGroup == NULL)
		error("freePairGroupList,invalidate free");
	free(list->pairGroup);
    list->pairGroup = NULL;
	free(list);
	list = NULL;
}

#endif //_STRUCTS_H
