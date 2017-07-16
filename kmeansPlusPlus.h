#ifndef KMEANS_PLUS_PLUS_H
#define KMEANS_PLUS_PLUS_H
#include "kmeans.h"
class KmeansPlusPlus: public Kmeans
{
public:
	KmeansPlusPlus(int k, int pointnumber);
	KmeansPlusPlus(int k);;
	
	void InitCenters();
	int NearestCenter(Point &p,int alreadyInitCenterNumber, float &minDistance);
private:

};




#endif
