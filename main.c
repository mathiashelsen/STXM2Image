#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include <wand/magick_wand.h>

#include <gsl/gsl_histogram.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort.h>

typedef struct
{
    double *** dataArray;
    int nChannels, sizeX, sizeY;
} rawData;

typedef struct
{
    int nChannels, sizeX, sizeY, normChannel;
    double *** normalizedDataArray; 
    gsl_histogram ** h;
    gsl_histogram ** p;
    double * maxPixelVal, * minPixelVal, * average, * standardDeviation;
} normalizedData;

typedef struct
{
    double time;
    int index;
} timeline;

int importData( FILE *inputFile, rawData * data);
int subtractChannel( rawData * data, int normChannel );
int subtractAvgPP( rawData * data );
int subtractAvgG( rawData * data );
int subtractChannelAvg( rawData * data );
int shiftModalPoint( rawData * data );
int shiftQuantile( rawData * data );
int shiftQuantileEdge( rawData * data, double *** mask );
int shiftMeanEdge( rawData * data, double *** mask );
int medianFilter( rawData * data, int nPixels );
int averageGaussian( rawData * data, double sigma, int pixels);
int createHistograms( rawData * data, normalizedData * result, int normChannel );
//int sort( rawData * result, double frequency, timeline *** sortingGuide );
int sort( rawData * data, double nChannels, double magicNumber, timeline *** sortingGuide );
int sortFunction( const void * a, const void * b );
int mergeChannels( rawData * data, rawData * result, timeline *** sortingGuide );
int drawChannel( normalizedData * data, int channelIndex );
int drawChannelFit( normalizedData * data, int channelIndex, double offset, double scale );
int drawChannelFitGlobal( normalizedData * data, int channelIndex, double mu, double sigma,
        double offset, double scale );
int drawMask( double *** mask, int sizeX, int sizeY );
int extractMask( rawData * data, double *** mask );
int calcFlatParams( normalizedData * data, double *mu, double *sigma);
int statisticalAnalysis( normalizedData * data );
void plotStatic( normalizedData *data);

int freeRawData( rawData * data );
int freeNormalizedData( normalizedData * data );

int main(int argc, char **argv)
{
    int i;
    FILE *inputFile;

    inputFile = fopen( argv[1], "rb" );
    if( inputFile == NULL )
    {
        fprintf( stdout, "File error.\n" );
        exit(1);
    }
    rawData raw;
    normalizedData norm;
    timeline *** t = malloc( sizeof(timeline **) );
    double scale = 1.0, offset = 0.0, sigmaSmooth = 1.0;
    double magicNumber = 1.0;
    int pixels = 0, medianPixels = 0;
    if( argc == 8 )
    {
        scale = atof( argv[2] );
        offset = atof( argv[3] );
        sigmaSmooth = atof( argv[4] );
        pixels = atoi( argv[5] );
        magicNumber = atof( argv[6] );
        medianPixels = atoi( argv[7] );
    }
    else
    {
        fprintf( stdout, "Usage: AnalyseSoftware [inputfile] [scale] [offset] [sigmaFilter] [pixelsFilter] [magic number] [medianPixels]\n" );
        return 0;
    }

    importData( inputFile, &raw );
    medianFilter( &raw, medianPixels );
    fclose( inputFile );

    double ** mask = malloc( sizeof(double **) * raw.sizeX );
    for(i = 0 ; i < raw.sizeX; i++ )
    {
        mask[i] = malloc( sizeof(double) * raw.sizeY );
    }
    extractMask( &raw, &mask );
    drawMask( &mask, raw.sizeX, raw.sizeY );
    shiftMeanEdge( &raw, &mask );
    //shiftQuantileEdge( &raw, &mask );

    //subtractAvgPP( &raw );
    if( pixels > 0 )
    {
        averageGaussian(&raw, sigmaSmooth, pixels);
    }
    sort( &raw, (double)(raw.nChannels), magicNumber, t );

    FILE * pixelFile = fopen( "pixels.txt", "w" );
    for( i = 0 ; i < raw.nChannels; i++ )
    {
        fprintf(pixelFile, "%d\t%f\n", i, raw.dataArray[i][0][0] );
    }
    fclose(pixelFile);

    int normChannel = raw.nChannels+2;
    createHistograms( &raw, &norm, normChannel );
    double mu, sigma;
    calcFlatParams( &norm, &mu, &sigma);
    statisticalAnalysis( &norm );

    FILE * histoFile = fopen("histograms.txt", "w");
    FILE * statFile = fopen("statistics.txt", "w");
    for( i = 0; i < norm.nChannels; i++ )
    {
        if( i != norm.normChannel )
        {
            fprintf( statFile, "%d\t%f\t%f\n", i, gsl_histogram_mean( norm.h[i] ), 
                    gsl_histogram_sigma( norm.h[i] ) );
            int j;
            double center, val, upper, lower;

            for( j = 0; j < gsl_histogram_bins( norm.h[i] ) ;  j++ )
            {
               gsl_histogram_get_range( norm.h[i], j, &lower, &upper );
               center = 0.5 * ( upper + lower );
               val = gsl_histogram_get( norm.h[i], j);
               fprintf( histoFile, "%f\t%f\n", center, val );
            }
            fprintf( histoFile, "\n\n" );
            //drawChannelFit( &norm, i, offset, scale );
            drawChannelFitGlobal( &norm, i, mu, sigma, offset, scale );
        }
    }
    fclose(histoFile);
    fclose(statFile);
    freeRawData( &raw );
    freeNormalizedData( &norm );
    return 0;
}

int importData( FILE *inputFile, rawData * data)
{
    uint8_t bufferSize1[4] = {0, 0, 0, 0};
    uint8_t bufferSize2[4] = {0, 0, 0, 0};
    uint8_t bufferSize3[4] = {0, 0, 0, 0};

    rewind(inputFile);
    printf("%d bytes read for bufferSize1\n", fread( bufferSize1, 1, sizeof(bufferSize1), inputFile));
    printf("%d bytes read for bufferSize1\n", fread( bufferSize2, 1, sizeof(bufferSize1), inputFile));
    printf("%d bytes read for bufferSize1\n", fread( bufferSize3, 1, sizeof(bufferSize1), inputFile));
    printf( "size buffer 1 = %d, buffer 2 = %d and buffer 3 = %d\n", bufferSize1[0], bufferSize2[0], bufferSize3[0] );
    printf( "size buffer 1 = %d, buffer 2 = %d and buffer 3 = %d\n", bufferSize1[1], bufferSize2[1], bufferSize3[1] );
    printf( "size buffer 1 = %d, buffer 2 = %d and buffer 3 = %d\n", bufferSize1[2], bufferSize2[2], bufferSize3[2] );
    printf( "size buffer 1 = %d, buffer 2 = %d and buffer 3 = %d\n", bufferSize1[3], bufferSize2[3], bufferSize3[3] );

    int i, j, k, l = 0;
    int size1, size2, size3, totalsize;
    int multiplier[4] = {16777216, 65536, 256, 1};
    size1 = 0;
    size2 = 0;
    size3 = 0;

    for( i = 0 ; i < 4; i++ )
    {
       size1 += bufferSize1[i]*multiplier[i]; 
       size2 += bufferSize2[i]*multiplier[i]; 
       size3 += bufferSize3[i]*multiplier[i]; 
    }

    totalsize = 4*size1*size2*size3;
    uint8_t *buffer;
    printf("Require %d bytes\n", (size_t)totalsize * sizeof(uint8_t));
    buffer = malloc( (size_t) totalsize * sizeof(uint8_t) );
    if( buffer == NULL )
    {
        fprintf( stdout, "Could not allocate memory.\n" );
        exit(1);
    }

    printf("# The array is %d x %d x %d, total bytes to be read: %d\n", size1, size2, size3, totalsize);
    fread( buffer, 1, totalsize , inputFile );

    uint32_t *buffer2;
    buffer2 = malloc( totalsize/4 * sizeof(uint32_t) );
    j = 0;
    for( i = 0; i < totalsize/4; i++ )
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
            exit(1);
        }
    }
    free( buffer );

    data->dataArray = malloc( size1 * sizeof(double **) );
    for( i = 0; i < size1; i++ )
    {
        (data->dataArray)[i] = malloc( size2 * sizeof(double *) );
        for( j = 0; j < size2; j++ )
        {
            (data->dataArray)[i][j] = malloc( size3 * sizeof(double) );
            for( k = 0; k < size3; k++ )
            {
                (data->dataArray)[i][j][k] = (double) buffer2[l];
                l++;
                if( l > totalsize/4 )
                {
                    fprintf( stdout, "Array index out of bounds.\n" );
                    exit(1);
                }
            }
        }

    }
    free( buffer2 );

    data->sizeX = size2;
    data->sizeY = size3;
    data->nChannels = size1;
    return 0;
}

int subtractChannel( rawData * data, int normChannel )
{
    int i, j, k;
    for( i = 0 ; i < data->nChannels; i++ )
    {
        for( j = 0; j < data->sizeX ; j++ )
        {
            for( k = 0 ; k < data->sizeY ; k++ )
            {
                data->dataArray[i][j][k] -= data->dataArray[normChannel][j][k];
            }
        }
    }
    return 0;
}

int subtractAvgPP( rawData * data)
{
    int i, j, k;
    double pixelValue;
    for( i = 0; i < data->sizeX; i++ )
    {
        for( j = 0 ; j < data->sizeY; j++ )
        {
            pixelValue = 0.0;
            for( k = 0 ; k < data->nChannels ; k++ )
            {
                pixelValue += (double) data->dataArray[k][i][j];

            }
            pixelValue /= (double) data->nChannels;
            for( k = 0; k < data->nChannels ; k++ )
            {
                data->dataArray[k][i][j] -= pixelValue;
            }
        }
    }
    return 0;
}

int subtractAvgG( rawData * data )
{
    int i, j, k;
    double pixelValue;
    for( i = 0; i < data->sizeX; i++ )
    {
        for( j = 0 ; j < data->sizeY; j++ )
        {
            pixelValue = 0.0;
            for( k = 0 ; k < data->nChannels ; k++ )
            {
                pixelValue += (double) data->dataArray[k][i][j];

            }
        }
    }
    pixelValue /= (double) data->nChannels;
    pixelValue /= (double) data->sizeX;
    pixelValue /= (double) data->sizeY;
    for( i = 0; i < data->sizeX; i++ )
    {
        for( j = 0 ; j < data->sizeY; j++ )
        {
            pixelValue = 0.0;
            for( k = 0 ; k < data->nChannels ; k++ )
            {
                data->dataArray[k][i][j] -= pixelValue;

            }
        }
    }
    return 0;
}

int subtractChannelAvg( rawData * data )
{
    int i, j, k;
    for( i = 0 ; i < data->nChannels; i++ )
    {
        double channelAvg = 0.0;
        for( j = 0 ; j < data->sizeX; j++ )
        {
            for(k = 0 ; k < data->sizeY; k++ )
            {
                channelAvg += data->dataArray[i][j][k];
            }
        }
        channelAvg /= (double) data->sizeX * (double) data->sizeY;
        for( j = 0 ; j < data->sizeX; j++ )
        {
            for(k = 0 ; k < data->sizeY; k++ )
            {
                data->dataArray[i][j][k] -= channelAvg;
            }
        }
    }

    return 0;
}

int averageGaussian( rawData * data, double sigma, int pixels )
{
    int lowerX, upperX, lowerY, upperY;
    lowerX = pixels;
    upperX = data->sizeX - pixels;
    lowerY = pixels;
    upperY = data->sizeY - pixels;
    int i, j, k;
    double *** dataArray = malloc( sizeof(double **) * data->nChannels );
    for( i = 0 ; i < data->nChannels ; i++ )
    {
        dataArray[i] = malloc( sizeof(double *) * data->sizeX );
        for( j = 0; j < data->sizeX; j++ )
        {
            dataArray[i][j] = malloc( sizeof(double) * data->sizeY );
            for( k = 0; k < data->sizeY; k++ )
            {
                int jj, kk;
                double pixelVal = 0.0;
                lowerX = ( (j - pixels) > 0 ) ? j - pixels : 0;
                upperX = ( (j + pixels) >= data->sizeX ) ? (data->sizeX - 1) : j + pixels ;
                lowerY = ( (k - pixels) > 0 ) ? k - pixels : 0;
                upperY = ( (k + pixels) >= data->sizeY ) ? (data->sizeY - 1) : k + pixels;

                for( jj = lowerX; jj <= upperX; jj++ )
                {
                    if( jj >= data->sizeX )
                    {
                        printf( "jj out of bounds\n" );
                    }
                    for( kk = lowerY; kk <= upperY; kk++ )
                    {
                        if( kk >= data->sizeY )
                        {
                            printf( "kk out of bounds\n" );
                        }
                        double r = sqrt( (kk-k)*(kk-k) + (jj-j)*(jj-j) );
                        pixelVal += data->dataArray[i][jj][kk] * exp( -0.5 * (r / sigma) * (r / sigma) ) / (sigma * 2.50662827463100050242) ;
                    }
                    dataArray[i][j][k] = pixelVal;
                }
            }
        }
    }
    for( i = 0 ; i < data->nChannels; i++ )
    {
        for( j = 0 ; j < data->sizeX; j++ )
        {
            free( data->dataArray[i][j] );
        }
        free( data->dataArray[i] );
    }
    free( data->dataArray );
    data->dataArray = dataArray;



    return 0;
}

int medianFilter( rawData * data, int nPixels )
{
    int i, j, k;
    double *** newValues = malloc( sizeof(double**) * data->nChannels );
    int lowerX, upperX, lowerY, upperY;
    for ( i = 0 ; i < data->nChannels ; i++ )
    {
        newValues[i] = malloc( sizeof(double*) * data->sizeX );
        for( j = 0 ; j < data->sizeX; j++ )
        {
            newValues[i][j] = malloc( sizeof(double) * data->sizeY );
            for( k = 0 ; k < data->sizeY; k++ )
            {

                lowerX = ( (j - nPixels) > 0) ? (j - nPixels) : 0;
                upperX = ( (j + nPixels + 1) < (data->sizeX) ) ? (j + nPixels + 1) : data->sizeY;
                lowerY = ( (k - nPixels) > 0) ? (k - nPixels) : 0;
                upperY = ( (k + nPixels + 1) < (data->sizeY) ) ? (k + nPixels + 1) : data->sizeY;
                //printf( "%d %d %d \t %d %d %d %d\n", i, j, k, lowerX, upperX, lowerY, upperY );
                int l, m, size = 0;
                double r, *pixels = malloc( sizeof(double) * (upperX - lowerX) * (upperY - lowerY));
                for( l = lowerX ; l < upperX ; l++ )
                {
                    for( m = lowerY ; m < upperY ; m++ )
                    {
                        r = (l - j)*(l - j) + (m - k)*(m - k);
                        if( r <= nPixels*nPixels ){
                            pixels[size] = data->dataArray[i][l][m]; 
                            size++;
                        }
                    }
                }
                gsl_sort( pixels, 1, size );
                newValues[i][j][k] = gsl_stats_median_from_sorted_data( pixels, 1, size );
                if( ( (newValues[i][j][k] - data->dataArray[i][j][k]) / data->dataArray[i][j][k] ) > 0.05 )
                {
                    //printf( "%d %d %d %f %f\n", i, j, k, data->dataArray[i][j][k], newValues[i][j][k]);
                }
            }
        }
    }
    for( i = 0 ; i < data->nChannels ; i++ )
    {
        for( j = 0 ; j < data->sizeX ; j++ )
        {
            free(data->dataArray[i][j]);
        }
        free(data->dataArray[i]);
    }
    free(data->dataArray);
    data->dataArray = newValues;
    return 0;
}

int createHistograms( rawData * data, normalizedData * result, int normChannel )
{
    //printf( "%d\n", data->nChannels );
    result->maxPixelVal = malloc( sizeof(double) * data->nChannels );
    result->minPixelVal = malloc( sizeof(double) * data->nChannels );
    result->average = malloc( sizeof(double) * data->nChannels );
    result->standardDeviation = malloc( sizeof(double) * data->nChannels );
    result->p = malloc( sizeof(gsl_histogram_pdf *) * data->nChannels );
    result->h = malloc( sizeof(gsl_histogram *) * data->nChannels );
    result->sizeX = data->sizeX;
    result->sizeY = data->sizeY;
    result->nChannels = data->nChannels;
    result->normalizedDataArray = data->dataArray;
    int i, j, k;
    for( i = 0; i < data->nChannels; i++ )
    {
        result->minPixelVal[i] = result->normalizedDataArray[i][0][0];
        result->maxPixelVal[i] = result->normalizedDataArray[i][0][0];
        for( j = 0; j < data->sizeX; j++ )
        {
            for( k = 0; k < data->sizeY; k++ )
            {
               result->maxPixelVal[i] = (result->maxPixelVal[i] > result->normalizedDataArray[i][j][k]) ? 
                   result->maxPixelVal[i] : result->normalizedDataArray[i][j][k];
                result->minPixelVal[i] = (result->minPixelVal[i] < result->normalizedDataArray[i][j][k]) ? 
                   result->minPixelVal[i] : result->normalizedDataArray[i][j][k];
            }
        }

        if( i != normChannel )
        {
            result->h[i] = gsl_histogram_alloc( data->sizeX * data->sizeY / 10 );
            gsl_histogram_set_ranges_uniform( result->h[i], (double) result->minPixelVal[i], (double) result->maxPixelVal[i] );
            for( j = 0; j < result->sizeX; j++ )
            {
                for( k = 0; k < result->sizeY; k++ )
                {   
                    gsl_histogram_increment( result->h[i], (double) result->normalizedDataArray[i][j][k] );
                }
            }
            result->p[i] = gsl_histogram_alloc( gsl_histogram_bins( result->h[i] ) );
            gsl_histogram_set_ranges_uniform( result->p[i], (double) result->minPixelVal[i] - 1.0, (double) result->maxPixelVal[i] + 1.0 );
            double sum = gsl_histogram_sum( result->h[i] );
            double upper, lower, center, val = 0;
            for( j = 0; j < gsl_histogram_bins( result->h[i] ); j++ )
            {
                gsl_histogram_get_range( result->h[i], j, &lower, &upper );
                center = 0.5*( lower + upper );
                val += (gsl_histogram_get( result->h[i], j ) / sum) ;
                gsl_histogram_accumulate( result->p[i], center, val );
            }
            result->average[i] = gsl_histogram_mean( result->h[i] );
            result->standardDeviation[i] = gsl_histogram_sigma( result->h[i] );
        }
        else
        {
            result->h[i] = NULL;
            result->p[i] = NULL;
        }
    }
    return 0;
}

int shiftModalPoint( rawData * data )
{
    int i, j, k;
    double * list = malloc(sizeof(double) * data->sizeX * data->sizeY );
    double * median = malloc( sizeof(double) * data->nChannels );
    double avgMedian;

    for( i = 0 ; i < data->nChannels ; i++ )
    {
        for ( j = 0 ; j < data->sizeX ; j++)
        {
            for( k = 0 ; k < data->sizeY; k++ )
            {
                list[ j * data->sizeY + k ] = data->dataArray[i][j][k];
            }
        }
        gsl_sort( list, 1, data->sizeX * data->sizeY );
        median[i] = gsl_stats_median_from_sorted_data( list, 1, data->sizeX * data->sizeY );
    }
    avgMedian = gsl_stats_mean( median, 1, data->nChannels );

    for( i = 0 ; i < data->nChannels ; i++ )
    {
        for ( j = 0 ; j < data->sizeX ; j++)
        {
            for( k = 0 ; k < data->sizeY; k++ )
            {
                data->dataArray[i][j][k] *= avgMedian / median[i];
            }
        }
    }
    free(list);
    free(median);
    return 0;
}

int shiftQuantile( rawData * data )
{
    int i, j, k;
    double * list = malloc(sizeof(double) * data->sizeX * data->sizeY );
    double * upperQuantile = malloc( sizeof(double) * data->nChannels );
    double * lowerQuantile = malloc( sizeof(double) * data->nChannels );
    double avgUpperQuantile, avgLowerQuantile;

    for( i = 0 ; i < data->nChannels ; i++ )
    {
        for ( j = 0 ; j < data->sizeX ; j++)
        {
            for( k = 0 ; k < data->sizeY; k++ )
            {
                list[ j * data->sizeY + k ] = data->dataArray[i][j][k];
            }
        }
        gsl_sort( list, 1, data->sizeX * data->sizeY );
        upperQuantile[i] = gsl_stats_quantile_from_sorted_data( list, 1, data->sizeX*data->sizeY, 0.8 );
        lowerQuantile[i] = gsl_stats_quantile_from_sorted_data( list, 1, data->sizeX*data->sizeY, 0.2 );

    }
    avgUpperQuantile = gsl_stats_mean( upperQuantile, 1, data->nChannels );
    avgLowerQuantile = gsl_stats_mean( lowerQuantile, 1, data->nChannels );

    double theta, phi;

    for( i = 0 ; i < data->nChannels ; i++ )
    {
        theta = (avgUpperQuantile - avgLowerQuantile) / (upperQuantile[i] - lowerQuantile[i]);
        phi = avgUpperQuantile - (upperQuantile[i]*theta);
        for ( j = 0 ; j < data->sizeX ; j++)
        {
            for( k = 0 ; k < data->sizeY; k++ )
            {
                data->dataArray[i][j][k] *= theta;
                data->dataArray[i][j][k] += phi;
            }
        }
    }
    free(list);
    free(upperQuantile);
    free(lowerQuantile);
    return 0;
}

int shiftQuantileEdge( rawData * data, double *** mask )
{
    int i, j, k, size = 0;
    double * list = malloc(sizeof(double) * data->sizeX * data->sizeY );
    double * upperQuantile = malloc( sizeof(double) * data->nChannels );
    double * lowerQuantile = malloc( sizeof(double) * data->nChannels );
    double avgUpperQuantile, avgLowerQuantile;

    for( i = 0 ; i < data->nChannels ; i++ )
    {
        size = 0;
        for ( j = 0 ; j < data->sizeX ; j++)
        {
            for( k = 0 ; k < data->sizeY; k++ )
            {
                if( (*mask)[j][k] == 1.0 )
                {
                    list[ size ] = data->dataArray[i][j][k];
                    size++;
                }
            }
        }
        gsl_sort( list, 1, size );
        upperQuantile[i] = gsl_stats_quantile_from_sorted_data( list, 1, size, 0.8 );
        lowerQuantile[i] = gsl_stats_quantile_from_sorted_data( list, 1, size, 0.2 );

    }
    avgUpperQuantile = gsl_stats_mean( upperQuantile, 1, data->nChannels );
    avgLowerQuantile = gsl_stats_mean( lowerQuantile, 1, data->nChannels );

    double theta, phi;

    for( i = 0 ; i < data->nChannels ; i++ )
    {
        theta = (avgUpperQuantile - avgLowerQuantile) / (upperQuantile[i] - lowerQuantile[i]);
        phi = avgUpperQuantile - (upperQuantile[i]*theta);
        for ( j = 0 ; j < data->sizeX ; j++)
        {
            for( k = 0 ; k < data->sizeY; k++ )
            {
                data->dataArray[i][j][k] *= theta;
                data->dataArray[i][j][k] += phi;
            }
        }
    }
    free(list);
    free(upperQuantile);
    free(lowerQuantile);
    return 0;
}

int shiftMeanEdge( rawData * data, double *** mask )
{
    int i, j, k, size = 0;
    double * list = malloc(sizeof(double) * data->sizeX * data->sizeY );
    double * mean = malloc( sizeof(double) * data->nChannels );
    double avgMean;

    for( i = 0 ; i < data->nChannels ; i++ )
    {
        size = 0;
        for ( j = 0 ; j < data->sizeX ; j++)
        {
            for( k = 0 ; k < data->sizeY; k++ )
            {
                if( (*mask)[j][k] == 1.0 )
                {
                    list[ size ] = data->dataArray[i][j][k];
                    size++;
                }
            }
        }
        mean[i] = gsl_stats_mean( list, 1, size );

    }
    avgMean = gsl_stats_mean( mean, 1, data->nChannels );
    double std = gsl_stats_sd_m( mean, 1, data->nChannels, avgMean );
    printf("# %f\t%f\n", avgMean, std );

    double theta, phi;

    for( i = 0 ; i < data->nChannels ; i++ )
    {
        theta = avgMean / mean[i];
        phi = 0.0;
        for ( j = 0 ; j < data->sizeX ; j++)
        {
            for( k = 0 ; k < data->sizeY; k++ )
            {
                data->dataArray[i][j][k] *= theta;
                data->dataArray[i][j][k] += phi;
            }
        }
    }
    free(list);
    free(mean);
    return 0;
}
/*
int sort( rawData * data, double frequency, timeline *** sortingGuide)
*/
int sort( rawData * data, double nChannels, double magicNumber, timeline *** sortingGuide )
{
    *sortingGuide  = malloc( sizeof(timeline *) * data->nChannels );
    int i = 0;
    for( i = 0 ; i < data->nChannels ; i++ )
    {
        (*sortingGuide)[i] = malloc( sizeof(timeline) );
        (*sortingGuide)[i]->index = i;
        (*sortingGuide)[i]->time = fmod( magicNumber*(double)i, nChannels )/magicNumber;
    }

    qsort( *sortingGuide, data->nChannels , sizeof(timeline *), sortFunction );
    double ***a = malloc( sizeof(double **) * data->nChannels );

  
    for( i = 0; i < data->nChannels; i++ )
    {
        a[i] = data->dataArray[ (*sortingGuide)[i]->index ];
    }

    free( data->dataArray );
    data->dataArray = a;

    return 0;
}

int sortFunction( const void * a, const void * b)
{
    const timeline *ia = *(timeline **) a;
    const timeline *ib = *(timeline **) b;
    if( ia->time > ib->time )
    {
        return 1;
    }
    else if( ia->time < ib->time )
    {
        return -1;
    }
    return 0;
}

int mergeChannels( rawData * data, rawData * result, timeline *** sortingGuide )
{
    int i, j, k, index;
    int * sorting = malloc( sizeof(int) * data->nChannels );
    sorting[0] = 0;
    for(i = 1 ; i < data->nChannels; i++ )
    {
        fprintf( stdout, "%f, %f\n", (*sortingGuide)[i]->time, (*sortingGuide)[i-1]->time );
       if(  (*sortingGuide)[i]->time == (*sortingGuide)[i-1]->time )
       {
            sorting[i] = sorting[i-1];
       }
       else
       {
           sorting[i] = sorting[i-1]+1;
       }
    }
    result->nChannels = sorting[ data->nChannels - 1 ] + 1;
    result->dataArray = malloc( sizeof( double ** ) * result->nChannels );
    result->sizeX = data->sizeX;
    result->sizeY = data->sizeY;
    for( i = 0; i < data->nChannels ; i++ )
    {
        result->dataArray[i] = malloc( sizeof( double *) * result->sizeX );
        index = sorting[i];
        for ( j = 0 ; j < result->sizeX; j++ )
        {
            result->dataArray[i][j] = malloc( sizeof( double ) * result->sizeY );
            result->dataArray[index][j][0] = data->dataArray[i][j][0];
            for ( k = 1 ; k < result->sizeY ; k++ )
            {
                result->dataArray[index][j][k] += data->dataArray[i][j][k];
            }
        }
    }

    return 0;
}

int drawChannel( normalizedData * data, int channelIndex )
{
    MagickWand *m_wand = NULL;
    PixelWand *p_wand = NULL;
    PixelIterator *iterator = NULL;
    PixelWand **pixels = NULL;
    int x,y, gray;
    char hex[128];

     size_t index;

    MagickWandGenesis();

    p_wand = NewPixelWand();
    PixelSetColor(p_wand, "white");
    m_wand = NewMagickWand();
    MagickNewImage(m_wand, data->sizeX, data->sizeY, p_wand);
    iterator = NewPixelIterator(m_wand);
    
    for(y=0; y< data->sizeX; y++)
    {
        pixels=PixelGetNextIteratorRow(iterator,(size_t*) &x);
        for(x = 0; x< data->sizeY; x++)
        {
            gsl_histogram_find( data->p[channelIndex], data->normalizedDataArray[channelIndex][x][y] , &index );

            gray = (int) (gsl_histogram_get( data->p[channelIndex], index ) * 255.0);
            sprintf(hex, "#%02x%02x%02x",gray,gray,gray);
            PixelSetColor(pixels[x], hex);
        }
        PixelSyncIterator(iterator);
    }
    char filename[64];
    sprintf( filename, "image_%d.gif", channelIndex);
    MagickWriteImage(m_wand, filename);

    iterator = DestroyPixelIterator(iterator);
    DestroyMagickWand(m_wand);
    MagickWandTerminus();

    return 0;
}

int drawChannelFit( normalizedData * data, int channelIndex , double offset, double scale)
{
    MagickWand *m_wand = NULL;
    PixelWand *p_wand = NULL;
    PixelIterator *iterator = NULL;
    PixelWand **pixels = NULL;
    int i,j, gray;
    char hex[128];
    double x, mu, sigma;

    MagickWandGenesis();

    p_wand = NewPixelWand();
    PixelSetColor(p_wand, "white");
    m_wand = NewMagickWand();
    MagickNewImage(m_wand, data->sizeX, data->sizeY, p_wand);
    iterator = NewPixelIterator(m_wand);
    
    for(j=0; j< data->sizeX; j++)
    {
        pixels=PixelGetNextIteratorRow(iterator,(size_t*) &i);
        for(i = 0; i< data->sizeY; i++)
        {
            x = data->normalizedDataArray[channelIndex][i][j];
            mu = data->average[channelIndex] * (1.0 + offset);
            sigma = data->standardDeviation[channelIndex] * scale;

            gray = (int) (255.0 * 0.5 * (1.0 + gsl_sf_erf( (x - mu) / (1.4142135623730950488 * sigma ) ) ) );
            sprintf(hex, "#%02x%02x%02x",gray,gray,gray);
            PixelSetColor(pixels[i], hex);
        }
        PixelSyncIterator(iterator);
    }
    char filename[64];
    sprintf( filename, "image_%d.gif", channelIndex);
    MagickWriteImage(m_wand, filename);

    iterator = DestroyPixelIterator(iterator);
    DestroyMagickWand(m_wand);
    MagickWandTerminus();

    return 0;
}

int drawChannelFitGlobal( normalizedData * data, int channelIndex, double mu, double sigma,
        double offset, double scale )
{
    MagickWand *m_wand = NULL;
    PixelWand *p_wand = NULL;
    PixelIterator *iterator = NULL;
    PixelWand **pixels = NULL;
    int i,j, gray;
    char hex[128];
    double x;

    MagickWandGenesis();

    p_wand = NewPixelWand();
    PixelSetColor(p_wand, "white");
    m_wand = NewMagickWand();
    MagickNewImage(m_wand, data->sizeX, data->sizeY, p_wand);
    iterator = NewPixelIterator(m_wand);
    
    for(j=0; j< data->sizeX; j++)
    {
        pixels=PixelGetNextIteratorRow(iterator,(size_t*) &i);
        for(i = 0; i< data->sizeY; i++)
        {
            x = data->normalizedDataArray[channelIndex][i][j];

            gray = (int) (255.0 * 0.5 * (1.0 + gsl_sf_erf( (x - mu - offset) / (1.4142135623730950488 * sigma * scale ) ) ) );
            sprintf(hex, "#%02x%02x%02x",gray,gray,gray);
            PixelSetColor(pixels[i], hex);
        }
        PixelSyncIterator(iterator);
    }
    char filename[64];
    sprintf( filename, "image_%d.gif", channelIndex);
    MagickWriteImage(m_wand, filename);

    iterator = DestroyPixelIterator(iterator);
    DestroyMagickWand(m_wand);
    MagickWandTerminus();

    return 0;
}

int drawMask( double *** mask, int sizeX, int sizeY )
{
    MagickWand *m_wand = NULL;
    PixelWand *p_wand = NULL;
    PixelIterator *iterator = NULL;
    PixelWand **pixels = NULL;
    int x,y, gray;
    char hex[128];

    MagickWandGenesis();

    p_wand = NewPixelWand();
    PixelSetColor(p_wand, "white");
    m_wand = NewMagickWand();
    MagickNewImage(m_wand, sizeX, sizeY, p_wand);
    iterator = NewPixelIterator(m_wand);
    
    for(y=0; y< sizeX; y++)
    {
        pixels=PixelGetNextIteratorRow(iterator,(size_t*) &x);
        for(x = 0; x< sizeY; x++)
        {
            gray = (int) ((*mask)[x][y] * 255.0);
            sprintf(hex, "#%02x%02x%02x",gray,gray,gray);
            PixelSetColor(pixels[x], hex);
        }
        PixelSyncIterator(iterator);
    }
    MagickWriteImage(m_wand, "mask.gif");

    iterator = DestroyPixelIterator(iterator);
    DestroyMagickWand(m_wand);
    MagickWandTerminus();

    return 0;
}

int freeRawData( rawData * data )
{
    int i, j;
    for( i = 0 ; i < data->nChannels ; i++ )
    {
        for( j = 0; j < data->sizeX; j++ )
        {
            free( data->dataArray[i][j] );
        }
        free( data->dataArray[i] );
    }
    free( data->dataArray );
    return 0;
}

int freeNormalizedData( normalizedData * data )
{
    int i;
    for( i = 0 ; i < data->nChannels ; i++ )
    {
        gsl_histogram_free( data->h[i] );
        gsl_histogram_free( data->p[i] );
    }
    free( data->p );
    free( data->h );
    free( data->maxPixelVal );
    free( data->minPixelVal );
    free( data->average );
    free( data->standardDeviation );
    return 0;
}

int extractMask( rawData * data, double *** mask )
{
    int i, j, k;
    double *pixels = malloc(sizeof(double) * data->sizeX * data->sizeY * data->nChannels );
    double **flat = malloc(sizeof(double *) * data->sizeX );
    for( j = 0 ; j < data->sizeX; j++ )
    {
        flat[j] = malloc(sizeof(double) * data->sizeY );
        for( k = 0 ; k < data->sizeY; k++ )
        {
            flat[j][k] = 0.0;
            for( i = 0 ; i < data->nChannels; i++ )
            {
                flat[j][k] += data->dataArray[i][j][k];
                pixels[(i * data->sizeX * data->sizeY) + (j * data->sizeY) + k] = data->dataArray[i][j][k];
            }
            flat[j][k] /= (double) data->nChannels;
        }
    }
    gsl_sort( pixels, 1, data->sizeX*data->sizeY*data->nChannels );
    double center = gsl_stats_quantile_from_sorted_data( pixels, 1, data->sizeX*data->sizeY*data->nChannels, 0.9 );
    for( i = 0 ; i < data->sizeX; i++ )
    {
        for( j = 0 ; j < data->sizeY; j++ )
        {
            if( flat[i][j] > center )
            {
                (*mask)[i][j] = 1.0;
            }
            else
            {
                (*mask)[i][j] = 0.0;
            }
        }
    }
    free( pixels );
    return 0;
}

int calcFlatParams( normalizedData * data, double *mu, double *sigma)
{
    int i, j, k, range;
    range = data->sizeX * data->sizeY;
    double *pixelList = malloc( sizeof(double) * range );
    double *sigmas = malloc( sizeof(double) * data->nChannels );
    double *mus = malloc( sizeof(double) * data->nChannels );
    for( i = 0; i < data->nChannels; i++ )
    {
        for( j = 0 ; j < data->sizeX; j++ )
        {
            for( k = 0 ; k < data->sizeY; k++ )
            {
                pixelList[ (data->sizeY * j) + k ] = data->normalizedDataArray[i][j][k];
            }
        }
        sigmas[i] = gsl_stats_sd( pixelList, 1, range );
        mus[i] = gsl_stats_mean( pixelList, 1, range );
    }
    *sigma = gsl_stats_sd( sigmas, 1, data->nChannels );
    *mu = gsl_stats_mean( mus, 1, data->nChannels );
    /*
    for( j = 0 ; j < data->sizeX; j++ )
    {
        for( k = 0 ; k < data->sizeY; k++)
        {
            pixelList[ (data->sizeY * j) + k ] = 0.0;
            for( i = 0; i < data->nChannels; i++ )
            {
                pixelList[ (data->sizeY * j) + k ] += data->normalizedDataArray[i][j][k];
            }
            pixelList[ (data->sizeY * j) + k ] /= (double) data->nChannels;
        }
    }
    *mu = gsl_stats_mean( pixelList, 1, range );
    *sigma = gsl_stats_sd_m( pixelList, 1, range, *mu );
    */
    return 0;
}

int statisticalAnalysis( normalizedData * data )
{
    FILE * statFile = fopen( "channelStatistics.txt", "w" );
    fprintf( statFile, "#channel\tmean\tsigma\tmedian\t0.25Quantile\t0.75Quantile\n" );
    int i, j, k, range;
    range = data->sizeX * data->sizeY;
    double *pixelList = malloc( sizeof(double) * range );
    double mean, sigma, median, upperQuantile, lowerQuantile;
    //printf( "%d", data->nChannels );
    for( i = 0 ; i < data->nChannels; i++ )
    {
        for( j = 0 ; j < data->sizeX; j++ )
        {
            for( k = 0; k < data->sizeY; k++ )
            {
                pixelList[ (j*data->sizeY) + k ] = data->normalizedDataArray[i][j][k];
            }
        }
        mean = gsl_stats_mean( pixelList, 1, range );
        sigma = gsl_stats_sd_m( pixelList, 1, range, mean );
        gsl_sort( pixelList, 1, range );
        median = gsl_stats_median_from_sorted_data( pixelList, 1, range );
        lowerQuantile = gsl_stats_quantile_from_sorted_data ( pixelList, 1, range, 0.25 );
        upperQuantile = gsl_stats_quantile_from_sorted_data ( pixelList, 1, range, 0.75 );
        fprintf(statFile,  "%d\t%f\t%f\t%f\t%f\t%f\n", i, mean, sigma, median, lowerQuantile, upperQuantile );
    }
    fclose( statFile );

    return 0;
}

void plotStatic( normalizedData *data)
{
    double **avg = malloc(sizeof(double *)*data->sizeX);
    int i =0;
    for( i = 0; i < data->sizeX ; i++ )
    {
	avg[i] = malloc(sizeof(double)*data->sizeY);
    }



    for( i = 0; i < data->sizeX ; i++ )
    {
	free(avg[i]);
    }
}
