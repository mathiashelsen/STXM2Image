#ifndef _CHANNEL_HPP
#define _CHANNEL_HPP

#include <assert.h>
#include <iostream>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>

class Channel
{
private:
    int rows;
    int cols;
    double *I;
public:
    Channel(int _rows, int _cols);
    ~Channel();

    void transferData( uint32_t *_I );
    void operator+=(Channel &B); 
    void operator-=(Channel &B); 
    void operator*=(Channel &B); 
    void operator/=(Channel &B); 

    void printParams() { std::cout << rows << ", " << cols << std::endl; } ;
    void drawChannel( const char *filename, int(* transferFunc)(double, void *args), void *funcArgs );
};

void importData(const char*filename, std::vector<Channel *> *channels);

#endif
