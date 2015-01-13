#include "Channel.hpp"

#include <stdlib.h>
#include <vector>

int linearMapping( double x, void *args)
{
    double *ptr = (double *)args;
    double min = ptr[0], max = ptr[1];
    return (int)(255.0 * (x-min)/(max-min));
};

int main(int argc, char **argv)
{
    assert(argc > 1);
    std::vector<Channel *> channels;
    importData(argv[1], &channels);

    Channel *a = channels.at(0);
    Channel staticImage(a->getRows(), a->getCols());

    for(int i = 0; i < channels.size(); i++)
    {
	a = channels.at(i);
	staticImage += *a;
    }

    double args[2];
    args[0] = staticImage.getMin();
    args[1] = staticImage.getMax();
    staticImage.drawChannel("static_image.gif", linearMapping, (void*)args);
    staticImage.toMask(atof(argv[2]));
    args[0] = staticImage.getMin();
    args[1] = staticImage.getMax();
    staticImage.drawChannel("mask.gif", linearMapping, (void*)args);

    a = channels.at(0);
    a->scaleMean(&staticImage);
    args[0] = a->getMin();
    args[1] = a->getMax(); 
    a->drawChannel("test.gif", linearMapping, (void*)args);
 
    return 0;
};
