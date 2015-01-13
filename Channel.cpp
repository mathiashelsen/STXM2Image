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
    avg = 0.0;
    max = (double) _I[0];
    min = (double) _I[0];
    for(int i = 0; i < (rows*cols); i++)
    {
	I[i] = (double)_I[i];
	max = (max > I[i]) ? max : I[i];
	min = (min < I[i]) ? min: I[i];
	avg += I[i];
    }
    avg /= (double) (rows*cols);
    std = 0.0;
    for(int i = 0; i < (rows*cols); i++)
    {
	std += (I[i] - avg)*(I[i] - avg);
    }
    std /= (double) (rows*cols - 1);
    std = sqrt(std);
};

void Channel::operator+=(Channel &B)
{
    assert(cols == B.cols);
    assert(rows == B.rows);
    max = I[0]+B.I[0];
    min = I[0]+B.I[0];
    for(int i = 0; i < cols; i++ )
    {
	for(int j = 0; j < rows; j++)
	{
	    I[i*rows + j] += B.I[i*rows + j];
	    max = (max > I[i*rows + j]) ? max : I[i*rows + j];
	    min = (min < I[i*rows + j]) ? min: I[i*rows + j];
	}
    }
};

void Channel::operator-=(Channel &B)
{
    assert(cols == B.cols);
    assert(rows == B.rows);
    max = I[0]+B.I[0];
    min = I[0]+B.I[0];
    for(int i = 0; i < cols; i++ )
    {
	for(int j = 0; j < rows; j++)
	{
	    I[i*rows + j] -= B.I[i*rows + j];
	    max = (max > I[i*rows + j]) ? max : I[i*rows + j];
	    min = (min < I[i*rows + j]) ? min: I[i*rows + j];
	}
    }
};

void Channel::operator*=(Channel &B)
{
    assert(cols == B.cols);
    assert(rows == B.rows);
    max = I[0]+B.I[0];
    min = I[0]+B.I[0];
    for(int i = 0; i < cols; i++ )
    {
	for(int j = 0; j < rows; j++)
	{
	    I[i*rows + j] *= B.I[i*rows + j];
	    max = (max > I[i*rows + j]) ? max : I[i*rows + j];
	    min = (min < I[i*rows + j]) ? min: I[i*rows + j];
	}
    }
};

void Channel::operator/=(Channel &B)
{
    assert(cols == B.cols);
    assert(rows == B.rows);
    max = I[0]+B.I[0];
    min = I[0]+B.I[0];
    for(int i = 0; i < cols; i++ )
    {
	for(int j = 0; j < rows; j++)
	{
	    I[i*rows + j] /= B.I[i*rows + j];
	    max = (max > I[i*rows + j]) ? max : I[i*rows + j];
	    min = (min < I[i*rows + j]) ? min: I[i*rows + j];
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
    MagickWand *m_wand = NULL;
    PixelWand *p_wand = NULL;
    PixelIterator *iterator = NULL;
    PixelWand **pixels = NULL;
    int gray;
    char hex[128];

    MagickWandGenesis();

    p_wand = NewPixelWand();
    PixelSetColor(p_wand, "white");
    m_wand = NewMagickWand();
    MagickNewImage(m_wand, cols, rows, p_wand);
    iterator = NewPixelIterator(m_wand);
    
    for(int y = 0; y < cols; y++)
    {
	int x = 0;
        pixels=PixelGetNextIteratorRow(iterator,(size_t*) &x);
        for(x = 0; x < rows; x++)
        {
            gray = transferFunc( I[y*rows + x], funcArgs );
            sprintf(hex, "#%02x%02x%02x",gray,gray,gray);
            PixelSetColor(pixels[x], hex);
        }
        PixelSyncIterator(iterator);
    }
    MagickWriteImage(m_wand, filename);

    iterator = DestroyPixelIterator(iterator);
    DestroyMagickWand(m_wand);
    MagickWandTerminus();
};

void Channel::toMask(double quantile)
{
    double *pixeltmp = new double[rows*cols];
    for(int i = 0; i < (rows*cols); i++)
    {
	pixeltmp[i] = I[i];
    }
    int quantile_pos = (int)((double)(rows*cols) * quantile);
    std::nth_element(pixeltmp, &(pixeltmp[quantile_pos]), &(pixeltmp[rows*cols-1]));
    double quant_pixel_val = pixeltmp[quantile_pos];
    delete[] pixeltmp;
    for(int i = 0; i < (rows*cols); i++)
    {
	I[i] = (I[i] > quant_pixel_val ? 1.0 : 0.0);
    }
    min = 0.0;
    max = 1.0;

}

void Channel::scaleMean(Channel *mask)
{
    double sum = 0.0, norm = 0.0;
    int Nvals = 0;  
    for(int i = 0; i < (rows*cols); i++)
    {
	if( mask->I[i] == 1.0 )
	{
	    sum += I[i];
	    Nvals++;
	}
    }

    norm = (double)Nvals/sum;
    std::cout << max << ", " << min << ", " << norm << std::endl;
    max *= norm;
    min *= norm;
    for(int i = 0; i < (rows*cols); i++)
    {
	I[i] *= norm;
    }
}
