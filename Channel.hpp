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
    void operator+=(Channel &B); 
    void operator-=(Channel &B); 
    void operator*=(Channel &B); 
    void operator/=(Channel &B); 

    void printParams() { std::cout << rows << ", " << cols << std::endl; } ;
    void drawChannel( const char *filename, int(* transferFunc)(double, void *args), void *funcArgs );
    void toMask(double quantile);
    void scaleMean(Channel *mask);

    double getMax() { return max; };
    double getMin() { return min; };
    double getAvg() { return avg; };
    double getStd() { return std; };
};

void importData(const char*filename, std::vector<Channel *> *channels);

#endif
