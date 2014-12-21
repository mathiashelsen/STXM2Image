#include "Channel.hpp"

Channel::Channel(int _rows, int _cols)
{
    rows = _rows;
    cols = _cols;
    I = new double[rows*cols];
};

Channel::~Channel()
{
    delete[] I;
};

void Channel::transferData(uint32_t *_I)
{
    //memcpy(I, _I, sizeof(double)*rows*cols);
    for(int i = 0; i < (rows*cols); i++)
    {
	I[i] = (double)_I[i];
    }
};

void Channel::operator+=(Channel &B)
{
    assert(cols == B.cols);
    assert(rows == B.rows);
    for(int i = 0; i < cols; i++ )
    {
	for(int j = 0; j < rows; j++)
	{
	    I[i*rows + j] += B.I[i*rows + j];
	}
    }
};

void Channel::operator-=(Channel &B)
{
    assert(cols == B.cols);
    assert(rows == B.rows);
    for(int i = 0; i < cols; i++ )
    {
	for(int j = 0; j < rows; j++)
	{
	    I[i*rows + j] -= B.I[i*rows + j];
	}
    }
};

void Channel::operator*=(Channel &B)
{
    assert(cols == B.cols);
    assert(rows == B.rows);
    for(int i = 0; i < cols; i++ )
    {
	for(int j = 0; j < rows; j++)
	{
	    I[i*rows + j] *= B.I[i*rows + j];
	}
    }
};

void Channel::operator/=(Channel &B)
{
    assert(cols == B.cols);
    assert(rows == B.rows);
    for(int i = 0; i < cols; i++ )
    {
	for(int j = 0; j < rows; j++)
	{
	    I[i*rows + j] /= B.I[i*rows + j];
	}
    }
};

void importData(const char*filename, std::vector<Channel *> *channels)
{
    FILE *datafile = fopen(filename, "r");
    assert(datafile != NULL);
    uint8_t bufferSize1[4] = {0, 0, 0, 0};
    uint8_t bufferSize2[4] = {0, 0, 0, 0};
    uint8_t bufferSize3[4] = {0, 0, 0, 0};

    rewind(datafile);
    fread( bufferSize1, 1, sizeof(bufferSize1), datafile);
    fread( bufferSize2, 1, sizeof(bufferSize1), datafile);
    fread( bufferSize3, 1, sizeof(bufferSize1), datafile);

    int size1, size2, size3, totalsize;
    int multiplier[4] = {16777216, 65536, 256, 1};
    size1 = 0;
    size2 = 0;
    size3 = 0;

    for(int i = 0 ; i < 4; i++ )
    {
       size1 += bufferSize1[i]*multiplier[i]; 
       size2 += bufferSize2[i]*multiplier[i]; 
       size3 += bufferSize3[i]*multiplier[i]; 
    }

    totalsize = 4*size1*size2*size3;
    uint8_t *buffer;
    //printf("Require %d bytes\n", (size_t)totalsize * sizeof(uint8_t));
    buffer = new uint8_t[totalsize];
    if( buffer == NULL )
    {
        fprintf( stdout, "Could not allocate memory.\n" );
        return;
    }

    fread( buffer, 1, totalsize , datafile );

    uint32_t *buffer2;
    buffer2 = new uint32_t[totalsize/4];
    int j = 0;
    for(int i = 0; i < totalsize/4; i++ )
    {
        buffer2[i] = buffer[j] * multiplier[0];
        buffer2[i] += buffer[j+1] * multiplier[1];
        buffer2[i] += buffer[j+2] * multiplier[2];
        buffer2[i] += buffer[j+3] * multiplier[3];
        buffer2[i] = buffer2[i] & (multiplier[0] - 1);
        j += 4;
        if( j > totalsize )
        {
            fprintf( stdout, "Array index out of bounds.\n" );
            return;
        }
    }
    delete[] buffer;

    for(int i = 0; i < size1; i++)
    {
	Channel *newChannel = new Channel(size2, size3);
	newChannel->transferData( &(buffer2[i*size2*size3]) );
	channels->push_back(newChannel);
    }
    delete[] buffer2;
};

void Channel::drawChannel( const char *filename, int(* transferFunc)(double, void *args), void *funcArgs )
{

};
