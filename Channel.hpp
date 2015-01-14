#ifndef _CHANNEL_HPP
#define _CHANNEL_HPP

#include <assert.h>
#include <algorithm>
#include <iostream>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>

#include <wand/magick_wand.h>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/density.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/moment.hpp>
#include <boost/accumulators/statistics/stats.hpp>

 
using namespace boost;
using namespace boost::accumulators;
 
typedef accumulator_set<double, features<tag::density, stats<tag::mean, tag::moment<2> > > > acc;
typedef accumulator_set<double, stats<tag::mean, tag::moment<2> > > accstat;
typedef iterator_range<std::vector<std::pair<double, double> >::iterator > histogram_type; 

class Channel
{
private:
    int rows;
    int cols;
    double *I;
    double max, min, avg, std;
public:
    Channel(int _rows, int _cols);
    ~Channel();

    int getRows(){ return rows; };
    int getCols(){ return cols; };

    void transferData( uint32_t *_I );
    void operator+=(const Channel &B); 
    void operator-=(const Channel &B); 
    void operator*=(const Channel &B); 
    void operator/=(const Channel &B); 
    void operator=(const Channel &B);

    void drawChannel( const char *filename, int(* transferFunc)(double, void *args), void *funcArgs );
    void toMask(double quantile);
    void scaleMean(Channel *mask);
    void extractHistogram(double **x, double **pdf, double **cdf, double *mu, double *sigma, int bins);
    void extractStats( double *stats );

    double getMax() { return max; };
    double getMin() { return min; };
    double getAvg() { return avg; };
    double getStd() { return std; };
};

void importData(const char*filename, std::vector<Channel *> *channels);

#endif
