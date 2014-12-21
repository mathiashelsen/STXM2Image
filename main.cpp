#include "Channel.hpp"

#include <stdlib.h>
#include <vector>

int main(int argc, char **argv)
{
    assert(argc > 1);
    std::vector<Channel *> channels;
    importData(argv[argc-1], &channels);

    for(int i = 0; i < channels.size(); i++)
    {
	Channel *a = channels.at(i);
    }
 
    return 0;
};
