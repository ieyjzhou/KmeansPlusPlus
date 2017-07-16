#include <chrono>
#include <random>
#include "kmeansPlusPlus.h"

KmeansPlusPlus::KmeansPlusPlus(int k, int pointnumber) :Kmeans(k, pointnumber)
{
	
}
KmeansPlusPlus::KmeansPlusPlus(int k) : Kmeans(k)
{

}
int KmeansPlusPlus::NearestCenter(Point &p, int alreadyInitCenterNumber, float &minDistance)
{
	minDistance = std::numeric_limits<float>::max();
	int k_id = -1;
	float dis;
	std::list<Point>::iterator centersIter = _Centers.begin();
#ifdef USING_OMP
#pragma ompparallel for 
#endif
	for (int k = 0; k <= alreadyInitCenterNumber; centersIter++, k++)
	{
		dis = Distance(p, *centersIter);
		if (dis < minDistance)
		{
			minDistance = dis;
			k_id = k;
		}
	}
	return k_id;
}
void KmeansPlusPlus::InitCenters()
{
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_int_distribution<int> distribution(0, _PointNumber - 1);
	int id = distribution(gen);
	std::list<Point>::iterator it = _Points.begin();
	std::advance(it, id);
#ifdef USING_OMP
#pragma ompparallel for 
#endif
	for (int i = 0; i < _K; i++)
	{
		_Centers.push_back(*it);
	}

	float sum,min_distance;
	std::list<Point>::iterator centersIter = _Centers.begin();
	std::list<Point>::iterator pointIter = _Points.begin();
	std::list<float> nearestDis(_PointNumber,0);
	std::list<float>::iterator floatIt = nearestDis.begin();
	
	for (int k = 1; k < _K; centersIter++, k++)
	{
		sum = 0;
		pointIter = _Points.begin();
		floatIt = nearestDis.begin();
#ifdef USING_OMP
#pragma ompparallel for 
#endif
		for (int p = 0; p < _PointNumber; pointIter++, p++)
		{
			NearestCenter(*pointIter,k, min_distance);
			*floatIt = min_distance;
			sum += min_distance;
			floatIt++;
		}
		std::random_device rd;  //Will be used to obtain a seed for the random number engine
		std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
		std::uniform_real_distribution<float> distribution(0.0, 1.0);
		float probability = distribution(gen);
		//std::cout << "orignional sum " << sum << std::endl;
		sum = sum*probability;
		//std::cout << "orignional sum " << sum << std::endl;
		pointIter = _Points.begin();
		floatIt = nearestDis.begin();
#ifdef USING_OMP
#pragma ompparallel for 
#endif
		for (int p = 0; p < _PointNumber; pointIter++, floatIt++, p++)
		{
			 sum =sum- *floatIt;
			 //std::cout <<"p " << p<< " sum "<<sum << std::endl;
			if (sum >0)
				continue;
			centersIter->_x = pointIter->_x;
			centersIter->_y = pointIter->_y;
			break;
		}
		 
	}
	/*float *d = malloc(sizeof(double)* len);

	point p, c;
	cent[0] = pts[rand() % len];
	for (n_cluster = 1; n_cluster < n_cent; n_cluster++) {
		sum = 0;
		for (j = 0, p = pts; j < len; j++, p++)
		{
			nearest(p, cent, n_cluster, d + j);
			sum += d[j];
		}
		sum = randf(sum);
		for (j = 0, p = pts; j < len; j++, p++)
		{
			if ((sum -= d[j]) > 0) continue;
			cent[n_cluster] = pts[j];
			break;
		}
	}
	for (j = 0, p = pts; j < len; j++, p++)
		p->group = nearest(p, cent, n_cluster, 0);*/
	
}