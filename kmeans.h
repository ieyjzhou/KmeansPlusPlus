#ifndef _KMEANS_H
#define _KMEANS_H

//#define  USING_OMP

#include<iostream> 
#include<list>
typedef struct Point
{
	Point(float x, float y, int g = 0)
	{
		_x = x;
		_y = y;
		_group = g;
	}
	float _x;
	float _y;
	int _group;
};
inline float Distance(Point &p1, Point &p2)
{
	return std::sqrtf((p1._x - p2._x)*(p1._x - p2._x) + (p1._y - p2._y)*(p1._y - p2._y));
}
class Kmeans
{
public:
	Kmeans(int k, int pointnumber);
	Kmeans(int k);
	void InitPoints();
	void InitPoints(std::list<Point> & pointlist);
	void InitPoints(std::string fileName);
	void InitCenters();
	void InitSpecifiedCenters();
	int NearestCenter(Point &p, float &minDistance);
	void RunKmean();
	void Cluster();
	void Center();
	void SaveEPS(std::string fileName);
	void PrintPointLis(std::list<Point> & pointList);
public:
	std::list<Point> _Points;
	std::list<Point> _Centers;
	int _MaxIteration;
	int _K;
	int _PointNumber;
	
} ;
#endif
