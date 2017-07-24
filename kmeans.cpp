#include<iostream> 
#include<list>
#include <chrono>
#include <random>
#include <fstream>
#include <iomanip>
#include"kmeans.h"
#include "simple_svg.hpp"
#include <map>


#ifdef USING_OMP
#include<omp.h>
#endif

Kmeans::Kmeans(int k)
{
	_K = k;

}
Kmeans::Kmeans(int k,int pointnumber)
{
	_K=k;
	_PointNumber = pointnumber;
}
void Kmeans::InitPoints()
{
	
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_real_distribution<float> distribution(0.0, 100.0);
	for (int i = 0; i < _PointNumber; i++)
  {
	  
		Point p(distribution(gen), distribution(gen));
	 // std::cout << p._x<<"\t" << p._y << std::endl;
	  _Points.push_back(p);
  }
}
// this function is used for compareing calculation time (Using omp)
void Kmeans::InitSpecifiedCenters()
{
#ifdef USING_OMP
#pragma ompparallel for 
#endif
	for (int i = 0; i < _K; i++)
	{


		std::list<Point>::iterator it = _Points.begin();
		std::advance(it, i);		
		_Centers.push_back(*it);
	}

}
void  Kmeans::InitCenters()
{
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_int_distribution<int> distribution(0, _PointNumber - 1);
#ifdef USING_OMP
#pragma ompparallel for 
#endif
	std::map<int, int> uniqueMap;
	while (uniqueMap.size() < _K)
	{
		int id = distribution(gen);
		uniqueMap.insert(std::pair<int, int>(id, id));
	}
	std::map<int, int>::iterator itMap = uniqueMap.begin();
	for (int i = 0; i < _K; i++)
	{

		int id = itMap->first;
		std::list<Point>::iterator it = _Points.begin();
		//std::advance(it, id);

		int count = 0;
		while (count != id)
		{
			it++;
			count++;
		}
		_Centers.push_back(*it);

		itMap++;
	}
}
void Kmeans::InitPoints(std::list<Point> & pointlist)
{
	_Points.assign(pointlist.begin(), pointlist.end());
	_PointNumber = pointlist.size();
}
int   Kmeans::NearestCenter(Point &p)
{
	
	float minDistance = std::numeric_limits<float>::max();
	int k_id = -1;
	float dis;
	std::list<Point>::iterator centersIter = _Centers.begin();
#ifdef USING_OMP
#pragma ompparallel for 
#endif
	for (int k = 0; k < _K; centersIter++, k++)
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
void  Kmeans::Cluster()
{
	
	std::list<Point>::iterator pointsIter = _Points.begin();

#ifdef USING_OMP
#pragma ompparallel for 
#endif
	float minDistance;
	for (int p = 0; p < _PointNumber; p++, pointsIter++)
	{
		pointsIter->_group = NearestCenter(*pointsIter);
				 
	}
}
void  Kmeans::Center()
{
	Point zeroPoint(0,0);
	std::vector<Point> center(_K, zeroPoint);	
	std::vector<int> count(_K, 0);
	std::list<Point>::iterator pointsIter = _Points.begin();
	std::list<Point>::iterator centerIter = _Centers.begin();
#ifdef USING_OMP
#pragma ompparallel for 
#endif
	for (int p = 0; p < _PointNumber; p++, pointsIter++)
	{
		center[pointsIter->_group]._x += pointsIter->_x;
		center[pointsIter->_group]._y += pointsIter->_y;
		count[pointsIter->_group]++;
	}
 
	for (int i = 0; i < center.size();i++)
	{
		center[i]._x = center[i]._x / (1.0*count[i]);
		center[i]._y = center[i]._y / (1.0*count[i]);
		centerIter->_x = center[i]._x;
		centerIter->_y = center[i]._y;
		centerIter++;
	}
	 
	 
}
void Kmeans::RunKmean()
{
	std::list<Point> oldCenter;
	_MaxIteration = 100;
	for (int iteration = 0; iteration < _MaxIteration; iteration++)
	{
		oldCenter = _Centers;
		//PrintPointLis(oldCenter);
		Cluster();
		Center();
		//std::cout << "-------------------------\n";
		//PrintPointLis(_Centers);
		float sum = 0;
		std::list<Point>::iterator currentIter = _Centers.begin();
		std::list<Point>::iterator oldIter = oldCenter.begin();
		for (int k = 0; k < _K;k++)
		{
			sum += Distance(*oldIter, *currentIter);
			oldIter++;
			currentIter++;
		}
		std::cout << "iteration "<< iteration<<" sum " << sum << std::endl;
		if (sum < 0.0001)
		{

			break;
		}
	}
}
 
void Kmeans::SaveEPS(std::string fileName)
{
	std::ofstream fs;
	fs.open(fileName.c_str(), std::ofstream::out);
	if (fs.is_open())
	{
		fs << std::fixed << std::setprecision(2);
		int W = 1000;
		int H = 1000;
		float min_x, max_x, min_y, max_y, scale, cx, cy;
		float *colors = new float [_K* 3];

		for (int i = 0; i < _K; i++)
		{
			colors[3 * i + 0] = (3 * (i + 1) % 11) / 11.;
			colors[3 * i + 1] = (7 * i % 11) / 11.;
			colors[3 * i + 2] = (9 * i % 11) / 11.;
		}

		max_x = max_y = -(min_x = min_y = HUGE_VAL);
		std::list<Point>::iterator pointsIter = _Points.begin();
		for (int j = 0; j < _PointNumber; j++, pointsIter++)
		{
			if (max_x < pointsIter->_x) max_x = pointsIter->_x;
			if (min_x > pointsIter->_x) min_x = pointsIter->_x;
			if (max_y < pointsIter->_y) max_y = pointsIter->_y;
			if (min_y > pointsIter->_y) min_y = pointsIter->_y;
		}
		scale = W / (max_x - min_x);
		if (scale > H / (max_y - min_y)) 
			scale = H / (max_y - min_y);
		cx = (max_x + min_x) / 2.0;
		cy = (max_y + min_y) / 2.0;
		 
		fs << "%%!PS-Adobe-3.0\n%%%%BoundingBox: -5 -5 " << W + 10 <<" "<< H + 10 << std::endl;
		
		fs << "/l {rlineto} def /m {rmoveto} def\n";
		fs << "/c { .25 sub exch .25 sub exch .5 0 360 arc fill } def\n";
		fs << "/s { moveto -2 0 m 2 2 l 2 -2 l -2 -2 l closepath ";
		fs << "	gsave 1 setgray fill grestore gsave 3 setlinewidth";
		fs << " 1 setgray stroke grestore 0 setgray stroke } def\n";
		std::list<Point>::iterator centersIter = _Centers.begin();
		for (int k = 0; k < _K; k++, centersIter++)
		{
			fs << colors[3 * k] << " " << colors[3 * k + 1] << " " << colors[3 * k + 2] << " setrgbcolor\n";
			pointsIter = _Points.begin();
			for (int j = 0; j < _PointNumber; j++, pointsIter++)
			{
				if (pointsIter->_group != k) 
					continue;
				fs << (pointsIter->_x - cx) * scale + W / 2 << " " << (pointsIter->_y - cy) * scale + H / 2 << " c\n";
					
			}
			fs << "\n0 setgray " << (centersIter->_x - cx) * scale + W / 2 << " " << (centersIter->_y - cy) * scale + H / 2 << " s\n";
		}
		fs<<"\n%%%%EOF";
		delete colors;

		std::cout << "Operation successfully performed\n";
		fs.close();
	}
	else
	{
		std::cout <<fileName.c_str() << " error opening file\n";
	}
}
void Kmeans::SaveSVG(std::string fileName)
{
	int W = 1000;
	int H = 1000;
	svg::Dimensions dimensions(W, H);
	svg::Document doc(fileName.c_str(), svg::Layout(dimensions, svg::Layout::BottomLeft));

	svg::Color pointColor(svg::Color::Red);
	svg::Color centerColor(svg::Color::Green);
	float min_x, max_x, min_y, max_y, scale, cx, cy;
	float *colors = new float[_K * 3];

	for (int i = 0; i < _K; i++)
	{
		colors[3 * i + 0] = (3 * (i + 1) % 11) / 11.;
		colors[3 * i + 1] = (7 * i % 11) / 11.;
		colors[3 * i + 2] = (9 * i % 11) / 11.;
	}
	max_x = max_y = -(min_x = min_y = HUGE_VAL);
	std::list<Point>::iterator pointsIter = _Points.begin();
	for (int j = 0; j < _PointNumber; j++, pointsIter++)
	{
		if (max_x < pointsIter->_x) max_x = pointsIter->_x;
		if (min_x > pointsIter->_x) min_x = pointsIter->_x;
		if (max_y < pointsIter->_y) max_y = pointsIter->_y;
		if (min_y > pointsIter->_y) min_y = pointsIter->_y;
	}
	scale = W / (max_x - min_x);
	if (scale > H / (max_y - min_y))
		scale = H / (max_y - min_y);
	cx = (max_x + min_x) / 2.0;
	cy = (max_y + min_y) / 2.0;
	
	
	std::list<Point>::iterator centerIter = _Centers.begin();
	for (int k = 0; k < _K; k++, centerIter++)
	{
		int r = colors[3 * k]*255;
		int g = colors[3 * k + 1] * 255;
		int b = colors[3 * k + 2] * 255;
		svg::Color c(r, g, b);
		
		pointsIter = _Points.begin();
		for (int j = 0; j < _PointNumber; j++, pointsIter++)
		{
			if (pointsIter->_group != k)
				continue;
			doc << svg::Circle(svg::Point((pointsIter->_x - cx) * scale + W / 2, (pointsIter->_y - cy) * scale + H / 2), 1, svg::Fill(svg::Color::White), svg::Stroke(1, c));

		}
		doc << svg::Rectangle(svg::Point((centerIter->_x - cx) * scale + W / 2, (centerIter->_y - cy) * scale + H / 2), 5,5, svg::Fill(svg::Color::White), svg::Stroke(2, c));
	}

	
	
	
	




	doc.save();
}
void Kmeans::PrintPointLis(std::list<Point> & pointList)
{
	std::list<Point>::iterator pointsIter = pointList.begin();
	for (; pointsIter != pointList.end();pointsIter++)
	{
		std::cout << pointsIter->_x << "\t" << pointsIter->_y << std::endl;
	}
}
void Kmeans::InitPoints(std::string fileName)
{
	std::ifstream fs;
	fs.open(fileName.c_str(), std::ofstream::in);
	if (fs.is_open())
	{
		float x, y;
		while (!fs.eof())
		{
			fs >> x >> y;
			_Points.push_back(Point(x, y));
		}
		_PointNumber = _Points.size();
		fs.close();
	}
	else
	{
		std::cout << fileName.c_str() << " error opening file\n";
	}
}