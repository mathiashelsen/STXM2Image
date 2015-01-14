#include "Channel.hpp"

#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include <boost/math/special_functions/erf.hpp>

#define NBINS 100
#define INVSQRT2 0.7071067811865475244

int linearMapping( double x, void *args)
{
    Channel *ptr = (Channel *)args;
    double min = ptr->getMin(), max = ptr->getMax();
    return (int)(255.0 * (x-min)/(max-min));
};

int erfMapping( double x, void *args)
{
    double *ptr = (double *)args;
    double mu = ptr[0], sigma = ptr[1];
    double tmp = INVSQRT2*(x-mu)/sigma;
    //template double boost::math::erf<>(double);
    return (int)(125.5 * boost::math::erf<>(-tmp));
    return 1;
};

int main(int argc, char **argv)
{
    assert(argc > 1);
    std::vector<Channel *> channels;
    importData(argv[1], &channels);

    Channel *a = channels.at(0);
    Channel staticImage(a->getRows(), a->getCols());
    Channel mask(a->getRows(), a->getCols());

    for(unsigned int i = 0; i < channels.size(); i++)
    {
	a = channels.at(i);
	staticImage += *a;
    }

    staticImage.drawChannel("static_image.gif", linearMapping, &staticImage);
    mask = staticImage;
    mask.toMask(0.6);
    mask.drawChannel("mask.gif", linearMapping, &mask);

    staticImage.scaleMean(&mask);

    double params[2];
    char filename[1024];
    for(unsigned int i = 0; i < channels.size(); i++)
    {
	a = channels.at(i);
	a->scaleMean(&mask);
	*a -= staticImage;
	a->extractStats( params );
	//params[1]*=0.8;
	sprintf(filename, "image_%03d.gif", i);
	a->drawChannel(filename, erfMapping, (void *)params);

    }
 
    return 0;
};
