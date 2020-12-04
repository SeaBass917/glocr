#include "imageprocessing.h"

#include <iostream>
#include <cmath>

/*
 * MISC
 */

void drawASmile(char const* addrOut){
    
    unsigned const srcRes = 8;
    unsigned const tgtRes = 1024;
    unsigned const scale = tgtRes/srcRes;

    // Lower resolution image template
    std::vector<std::vector<unsigned>> imgTemplate {
        {0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 1, 0, 0, 1, 0, 0},
        {0, 0, 1, 0, 0, 1, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0},
        {0, 1, 0, 0, 0, 0, 1, 0},
        {0, 1, 0, 0, 0, 0, 1, 0},
        {0, 0, 1, 1, 1, 1, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0},
    };
    
    // Create the image
    png::image< png::rgb_pixel > image(tgtRes, tgtRes);
    for (unsigned y = 0; y < tgtRes; y++){
        unsigned ySrc = y/scale;
        for (unsigned x = 0; x < tgtRes; x++){
            unsigned xSrc = x/scale;
            
            if(imgTemplate[ySrc][xSrc]) image[y][x] = png::rgb_pixel(0, 0, 0);
            else image[y][x] = png::rgb_pixel(127, 200, 0);
        }
    }

    // Gaussian blur
    unsigned const masklength = 25u;
    float const sigma = 54.0f;
    
    // Write image out
    image.write(addrOut);
}

void drawASmile(std::string addrOut){
    drawASmile(addrOut.c_str());
}

template<class T>
std::vector<std::vector<T>> gaussiankernel(unsigned const masklength, T const sigma){

    std::vector<std::vector<T>> kernel(masklength);

    // Mask length must be odd
    if(masklength & 1 == 1){
        const int radius = (masklength - 1) / 2;

        // intialising standard deviation
        T var = 2.0 * sigma * sigma;
        T denom = M_PI * var;

        // Fill the kernel row by row
        T sum = 0.0; // sum is for normalization 
        for (int x = -radius; x <= radius; x++) {
            std::vector<T> row(masklength);
            T rx = x * x;
            for (int y = -radius; y <= radius; y++) {
                T Gxy = exp(-(rx + y * y) / var) / denom;
                row[y + radius] = Gxy;
                sum += Gxy; 
            } 
            kernel[x + radius] = row;
        }

        // normalising the Kernel 
        for (int i = 0; i < masklength; ++i){
            for (int j = 0; j < masklength; ++j){
                kernel[i][j] /= sum;
            }
        }
    }
    
    return kernel;
} 

/*
 * Cropping
 */

png::image<png::gray_pixel> crop(png::image<png::gray_pixel> const &img, unsigned const x0, unsigned const x1, unsigned const y0, unsigned const y1){

    // Bounds checking
    unsigned const width = img.get_width();
    unsigned const height = img.get_height();
    if( 0 <= x0 && x0 < x1 && x1 < width && 
        0 <= y0 && y0 < y1 && y1 < height){

        // Prepare a new image for the crop
        unsigned const widthNew = x1 - x0 + 1;
        unsigned const heightNew = y1 - y0 + 1;
        if(0 < widthNew && 0 < heightNew){
            png::image<png::gray_pixel> imgCropped(widthNew, heightNew);

            // Loop through each image and move the data
            for(unsigned y = y0, yy = 0; y <= y1; y++, yy++){
                std::vector<png::gray_pixel> rowTgt = imgCropped[yy];
                std::vector<png::gray_pixel> rowSrc = img[y];
                for(unsigned x = x0, xx = 0; x <= x1; x++, xx++){
                    rowTgt[xx] = rowSrc[x];
                }
                imgCropped[yy] = rowTgt;
            }
            
            return imgCropped;
        }
        else{
            std::cerr << "\tERROR! Attempted to crop 0 width/height image: x0={"<<x0<<"} x1={"<<x1<<"} y0={"<<y0<<"} y1={"<<y1<<"}." << std::endl;
            throw std::exception();
        }
    }
    else{
        std::cerr << "\tERROR! Cannot crop image with dimmensions: (width="<<width<<", height="<<height<<") at x0={"<<x0<<"} x1={"<<x1<<"} y0={"<<y0<<"} y1={"<<y1<<"}." << std::endl;
        throw std::exception();
    }
}

png::image<png::gray_pixel> staggeredYCrop(png::image<png::gray_pixel> const& img, std::vector<unsigned> const& yTopArray, std::vector<unsigned> const& yBotArray, unsigned const yPadding){

    // Img info
    unsigned height = img.get_height();
    unsigned width = img.get_width();

    // Determine the min and max Y values used for image allocation
    unsigned yTopMin = yTopArray[0];
    unsigned yTopMax = yTopMin;
    unsigned yBotMax = yBotArray[0];
    unsigned yBotMin = yBotMax;
    for(unsigned i = 1; i < width; i++){
        int yTop = yTopArray[i];
        int yBot = yBotArray[i];
        if(yTop < yTopMin){
            yTopMin = yTop;
        }
        if(yTopMax < yTop){
            yTopMax = yTop;
        }
        if(yBot < yBotMin){
            yBotMin = yBot;
        }
        if(yBotMax < yBot){
            yBotMax = yBot;
        }
    }
    if(yTopMin < yBotMax && yTopMax < yBotMin){

        // Height of the cropped image
        unsigned heightCropped = yBotMax-yTopMin + 2 * yPadding + 1;

        // Initialize output image with whitespace
        png::image<png::gray_pixel> imgCropped(width, heightCropped);
        for(unsigned y = 0; y < heightCropped; y++){
            std::vector<png::gray_pixel> row = imgCropped[y];
            for(unsigned x = 0; x < width; x++)
                row[x] = 255;
            imgCropped[y] = row;
        }

        // Crop out the pixels and store them in the output image
        // Basically read the image sideways
        // TODO: See if its worth it to transpose first and read rows of memory with memcpys
        //       Dont forget that we want whitespace default
        for(unsigned x = 0; x < width; x++){
            unsigned y0 = yTopArray[x];
            unsigned y1 = yBotArray[x];

            // NOTE: we compute where to start in the output image as 
            // difference between the heighest y0 + padding
            unsigned yyStart = y0 - yTopMin + yPadding;
            for(unsigned y = y0, yy = yyStart; y <= y1; y++, yy++){
                imgCropped[yy][x] = img[y][x];
            }
        }

        return imgCropped;
    }
    else{
        std::cerr << "\tWARNING! Min Array crosses over Max Array staggeredYCrop()." << std::endl;
        throw std::exception();
    }
}

png::image<png::gray_pixel> staggeredXCrop(png::image<png::gray_pixel> const& img, std::vector<unsigned> const& xTopArray, std::vector<unsigned> const& xBotArray, unsigned const xPadding){

    // Img info
    unsigned height = img.get_height();
    unsigned width = img.get_width();

    // Determine the min and max Y values used for image allocation
    unsigned xTopMin = xTopArray[0];
    unsigned xTopMax = xTopMin;
    unsigned xBotMax = xBotArray[0];
    unsigned xBotMin = xBotMax;
    for(unsigned i = 1; i < height; i++){
        int xTop = xTopArray[i];
        int xBot = xBotArray[i];
        if(xTop < xTopMin){
            xTopMin = xTop;
        }
        if(xTopMax < xTop){
            xTopMax = xTop;
        }
        if(xBot < xBotMin){
            xBotMin = xBot;
        }
        if(xBotMax < xBot){
            xBotMax = xBot;
        }
    }
    if(xTopMin < xBotMax && xTopMax < xBotMin){

        // Height of the cropped image
        unsigned widthCropped = xBotMax-xTopMin + 2 * xPadding + 1;

        // Initialize output image with whitespace
        png::image<png::gray_pixel> imgCropped(widthCropped, height);
        for(unsigned y = 0; y < height; y++){
            std::vector<png::gray_pixel> row = imgCropped[y];
            for(unsigned x = 0; x < widthCropped; x++)
                row[x] = 255;
            imgCropped[y] = row;
        }

        // Crop out the pixels and store them in the output image
        for(unsigned y = 0; y < height; y++){
            unsigned x0 = xTopArray[y];
            unsigned x1 = xBotArray[y];
            std::vector<png::gray_pixel> rowSrc = img[y];
            std::vector<png::gray_pixel> rowTgt = imgCropped[y];
            for(unsigned x = 0; x < widthCropped; x++){

                // NOTE: we compute where to start in the output image as 
                // difference between the heighest x0 + padding
                unsigned xxStart = x0 - xTopMin + xPadding;
                for(unsigned x = x0, xx = xxStart; x <= x1; x++, xx++){
                    rowTgt[xx] = rowSrc[x];
                }
            }
            imgCropped[y] = rowTgt;
        }

        return imgCropped;
    }
    else{
        std::cerr << "\tWARNING! Min Array crosses over Max Array staggeredXCrop()." << std::endl;
        throw std::exception();
    }
}

/*
 * Segmentation
 */

png::image<png::gray_pixel> segmentImageThreshold(png::image<png::gray_pixel> const &imgSrc, unsigned const T){

    // Image dimensions
    unsigned height = imgSrc.get_height();
    unsigned width = imgSrc.get_width();

    // Create and return the segmented image
    png::image< png::gray_pixel > imgSeg(width, height);
#pragma omp parallel for
    for(unsigned y = 0; y < height; y++){
        std::vector<png::gray_pixel> rowSrc = imgSrc[y];
        std::vector<png::gray_pixel>& rowTgt = imgSeg[y];
        for(unsigned x = 0; x < width; x++){
            if(T <= rowSrc[x]) rowTgt[x] = 255;
            else rowTgt[x] = 0;
        }
    }

    return imgSeg;
}

unsigned determineKittlerThreshold(png::image<png::gray_pixel> const &imgSrc){

    // Image dimensions
    unsigned height = imgSrc.get_height();
    unsigned width = imgSrc.get_width();
    unsigned numPixels = height * width;

    // Determine the threshold for segmenting the image
    // Using Kittler's Method (iterative)

    // Generate grayscale histogram from image source
    double hist[256] = {0.0};
    for(unsigned y = 0; y < height; y++){
        for(unsigned x = 0; x < width; x++){
            hist[ imgSrc[y][x] ]++;
        }
    }

    // Normalize the histogram
    // DURING THE LOOP: Compute total mu for later
    double mu = 0.0;
    for(unsigned i = 0; i < 256; i++){
        double v = hist[i] / numPixels;
        hist[i] = v;
        mu += i * v;
    }

    // Look for the first and last non-zero bins 
    unsigned i_first = 0u;
    unsigned i_last = 255u;
    for(unsigned i = 0; i < 255; i++){
        if(0 < hist[i]){
            i_first = i;
            break;
        }
    }
    for(unsigned i = 255; 0 < i; i--){
        if(0 < hist[i]){
            i_last = i;
            break;
        }
    }

    // Initlialize our iterative variables 
    double errMin = std::numeric_limits<double>::max(); // The argument to minimize
    unsigned T = 0;                                     // Lowest error thrteshold
    double q1 = hist[i_first];                          // ratio in the bimodal model 1st mode
    double q2 = 1.0 - q1;                               // ratio in the bimodal model 2nd mode
    double mu1 = i_first;                               // Mean of 1st mode
    double mu2 = mu / q2;                               // Mean of 2nd mode
    double var1 = 0.0;                                  // Variance of 1st mode
    double var2 = 0.0;                                  // Variance of 2nd mode
    for(unsigned i = i_first+1; i < i_last+1; i++){
        double d = i - mu2;
        var2 += d*d*hist[i] / q2;
    }

    // Iteratively compute the error at each threshold t
    double mu1prev;
    double mu2prev;
    double var1prev;
    double var2prev;
    double q1prev;
    double q2prev;
    for(unsigned t = i_first+1; t < i_last-1; t++){
        
        // Store the prev computations
        mu1prev = mu1;
        mu2prev = mu2;
        var1prev = var1;
        var2prev = var2;
        q1prev = q1;
        q2prev = q2;
        
        // Update
        double p = hist[t];
        q1 = q1prev + p;
        q2 = 1 - q1;
        mu1 = (q1prev*mu1prev + t*p)/q1;
        mu2 = (mu - q1*mu1)/q2;
        double d1 = mu1prev - mu1;
        double d2 = mu2prev - mu2;
        var1 = (q1prev*(var1prev + (d1)*(d1)) + p*(t - mu1)*(t - mu1))/q1;
        var2 = (q2prev*(var2prev + (d2)*(d2)) - p*(t - mu2)*(t - mu2))/q2;
        
        // Compute error and update the minargs
        double err = (q1*log(var1) + q2*log(var2)) / 2.0 - q1*log(q1) - q2*log(q2);
        if(err < errMin){
            errMin = err;
            T = t;
        }
    }

    return T;
}

/*
 * Histograms
 */

std::vector<unsigned> horizontalProjectionHistogram(png::image<png::gray_pixel> const &img, unsigned const thresh, unsigned const binWidth){

    // Img info
    unsigned height = img.get_height();
    unsigned width = img.get_width();

    std::vector<unsigned> hist(height / binWidth + 1);

    // For each row count pixels that are below the threshold
    for(unsigned r = 0; r < height; r++){
        unsigned rr = r / binWidth;

        // Read row once
        const std::vector<png::byte, std::allocator<png::byte>> row = img[r];

        unsigned sum = 0;
        for(unsigned c = 0; c < width; c++){
            if(row[c] < thresh){
                sum++;
            }
        }

        hist[rr] = sum;
    }

    return hist;
}

template<typename depth>
std::vector<depth> horizontalProjectionHistogramNorm(png::image<png::gray_pixel> const &img, unsigned const thresh, unsigned const binWidth){

    // Img info
    unsigned height = img.get_height();
    unsigned width = img.get_width();

    std::vector<depth> hist(height / binWidth + 1);

    // For normilization
    unsigned totalSum = 0;

    // For each row count pixels that are below the threshold
    for(unsigned r = 0; r < height; r++){
        unsigned rr = r / binWidth;

        // Read row once
        const std::vector<png::gray_pixel> row = img[r];

        unsigned sum = 0;
        for(unsigned c = 0; c < width; c++){
            if(row[c] < thresh){
                sum++;
                totalSum++;
            }
        }

        hist[rr] = sum;
    }

    // Normalize the histogram
    if(0 < totalSum){
        for(unsigned r = 0; r < height; r++)
            hist[r] /= (depth)totalSum;
    }

    return hist;
}

std::vector<unsigned> verticalProjectionHistogram(png::image<png::gray_pixel> const &img, unsigned const thresh, unsigned const binWidth){

    // Img info
    unsigned height = img.get_height();
    unsigned width = img.get_width();

    std::vector<unsigned> hist(width / binWidth + 1);

    // For each col count pixels that are below the threshold
    for(unsigned c = 0; c < width; c++){
        unsigned cc = c / binWidth;

        unsigned sum = 0;
        for(unsigned r = 0; r < height; r++){
            if(img[r][c] < thresh){
                sum++;
            }
        }

        hist[cc] = sum;
    }

    return hist;
}

template<typename depth>
std::vector<depth> verticalProjectionHistogramNorm(png::image<png::gray_pixel> const &img, unsigned const thresh, unsigned const binWidth){
    
    // Img info
    unsigned height = img.get_height();
    unsigned width = img.get_width();

    std::vector<depth> hist(width / binWidth + 1);

    // For normilization
    unsigned totalSum = 0;

    // For each col count pixels that are below the threshold
    for(unsigned c = 0; c < width; c++){
        unsigned cc = c / binWidth;

        unsigned sum = 0;
        for(unsigned r = 0; r < height; r++){
            if(img[r][c] < thresh){
                sum++;
                totalSum++;
            }
        }

        hist[cc] = sum;
    }

    // Normalize the histogram
    if(0 < totalSum){
        for(unsigned r = 0; r < width; r++)
            hist[r] /= (depth)totalSum;
    }

    return hist;
}

std::vector<tuple2<int>> getMidPointBounds(std::vector<unsigned> hist, unsigned const peakThreshold, unsigned const minNoiseFloor, unsigned const minGapThresh){
    
    unsigned const histLen = hist.size();

    // Determine Transition points
    // Look for a peak
    // sweep out until hitting 0
    std::vector<int> transitionPoints;
    unsigned i = 0;
    while(i < histLen){
        if(peakThreshold < hist[i]){
            for(int ii = i; 0 <= ii; ii--){
                if(hist[ii] <= minNoiseFloor || 0 == ii){
                    transitionPoints.push_back(ii);
                    break;
                }
            }
            for(i; i < histLen; i++){
                if(hist[i] <= minNoiseFloor || i == histLen-1){
                    transitionPoints.push_back(i);
                    break;
                }
            }
        }

        i++;
    }
        
    // Check that the number of transition points makes sense
    // Value should be even
    unsigned numTransitions = transitionPoints.size();
    unsigned numMidPoints = numTransitions / 2 - 1;
    if((numTransitions & 1) == 0){

        // Get a list of the midpoints
        std::vector<tuple2<int>> midPointBounds;
        for(unsigned i=2; i < numTransitions; i+=2){
            int tp = transitionPoints[i];
            int tpPrev = transitionPoints[i-1];
            int gapDistance = tp - tpPrev;
            if(minGapThresh <= gapDistance){
                tuple2<int> midPointBound = {tpPrev, tp};
                midPointBounds.push_back(midPointBound);
            }
        }

        return midPointBounds;
    }
    else{
        std::cerr << "\tERROR! getMidPointBounds() Found an odd number of segments. Cannot confidently segment document." << std::endl;
        throw std::exception();
    }
}

/*
 * GRADIENTS
 */

float computeGradientVector_sobel(png::image<png::gray_pixel> const &image, unsigned const r, unsigned const c, unsigned const width, unsigned const height){

    // Initialize variables for the readablity
    int sA, sB, sC, sD, sF, sG, sH, sI;

    // Mirror on the edges
    unsigned rPrev = (r == 0)? r + 1 : r - 1;
    unsigned rNext = (r == height - 1)? r - 1 : r + 1;
    unsigned cPrev = (c == 0)? c + 1 : c - 1;
    unsigned cNext = (c == width - 1)? c - 1 : c + 1;

    // Read the prev row
    const std::vector<png::byte, std::allocator<png::byte>> rowPrev = image[rPrev];
    sA = (int)rowPrev[cPrev];
    sB = (int)rowPrev[c];
    sC = (int)rowPrev[cNext];
    // Read the current row
    const std::vector<png::byte, std::allocator<png::byte>> row = image[r];
    sD = (int)row[cPrev];
    // sE = row[c] (Not used)
    sF = (int)row[cNext];
    // Read the next row
    const std::vector<png::byte, std::allocator<png::byte>> rowNext = image[rNext];
    sG = (int)rowNext[cPrev];
    sH = (int)rowNext[c];
    sI = (int)rowNext[cNext];

    // Compute the gradient components
    float gx = sC + 2.0f*sF + sI - sA - 2.0f*sD - sG;
    float gy = sG + 2.0f*sH + sI - sA - 2.0f*sB - sC;

    // Return the magnitude
    return sqrtf(gx*gx + gy*gy) / 8.0f;
}

tuple2<float> computeDirectionalGradientVector_sobel(png::image<png::gray_pixel> const &image, unsigned const r, unsigned const c, unsigned const width, unsigned const height){

    // Initialize variables for the readablity
    int sA, sB, sC, sD, sF, sG, sH, sI;

    // Mirror on the edges
    unsigned rPrev = (r == 0)? r + 1 : r - 1;
    unsigned rNext = (r == height - 1)? r - 1 : r + 1;
    unsigned cPrev = (c == 0)? c + 1 : c - 1;
    unsigned cNext = (c == width - 1)? c - 1 : c + 1;

    // Read the prev row
    const std::vector<png::byte, std::allocator<png::byte>> rowPrev = image[rPrev];
    sA = (int)rowPrev[cPrev];
    sB = (int)rowPrev[c];
    sC = (int)rowPrev[cNext];
    // Read the current row
    const std::vector<png::byte, std::allocator<png::byte>> row = image[r];
    sD = (int)row[cPrev];
    // sE = row[c] (Not used)
    sF = (int)row[cNext];
    // Read the next row
    const std::vector<png::byte, std::allocator<png::byte>> rowNext = image[rNext];
    sG = (int)rowNext[cPrev];
    sH = (int)rowNext[c];
    sI = (int)rowNext[cNext];

    // Compute the gradient components
    float gx = sC + 2.0f*sF + sI - sA - 2.0f*sD - sG;
    float gy = sG + 2.0f*sH + sI - sA - 2.0f*sB - sC;

    // Package the gradient vector into a tuple and return
    tuple2<float> grad;
    grad._0 = gx;
    grad._1 = gy;

    return grad;
}

float computeGradientVector_median(png::image<png::gray_pixel> const &image, unsigned const r, unsigned const c, unsigned const width, unsigned const height){

    // Initialize variables for the readablity
    int sA, sB, sC, sD, sF, sG, sH, sI;

    // Mirror on the edges
    unsigned rPrev = (r == 0)? r + 1 : r - 1;
    unsigned rNext = (r == height - 1)? r - 1 : r + 1;
    unsigned cPrev = (c == 0)? c + 1 : c - 1;
    unsigned cNext = (c == width - 1)? c - 1 : c + 1;

    // Read the prev row
    const std::vector<png::byte, std::allocator<png::byte>> rowPrev = image[rPrev];
    sA = (int)rowPrev[cPrev];
    sB = (int)rowPrev[c];
    sC = (int)rowPrev[cNext];
    // Read the current row
    const std::vector<png::byte, std::allocator<png::byte>> row = image[r];
    sD = (int)row[cPrev];
    // sE = row[c] (Not used)
    sF = (int)row[cNext];
    // Read the next row
    const std::vector<png::byte, std::allocator<png::byte>> rowNext = image[rNext];
    sG = (int)rowNext[cPrev];
    sH = (int)rowNext[c];
    sI = (int)rowNext[cNext];

    // Fill the sample arrays
    int x0[5] = {
        sA, sB, sD, sG, sH
    };
    int x1[5] = {
        sB, sC, sF, sH, sI
    };
    int y0[5] = {
        sA, sB, sC, sD, sF
    };
    int y1[5] = {
        sD, sF, sG, sH, sI
    };

    // NOTE: the typecast to float so that the gradient datatype is consistent
    float gx = (float)(median(x1, 5) - median(x0, 5));
    float gy = (float)(median(y1, 5) - median(y0, 5));

    // Return the magnitude
    return sqrtf(gx*gx + gy*gy);
}

tuple2<float> computeDirectionalGradientVector_median(png::image<png::gray_pixel> const &image, unsigned const r, unsigned const c, unsigned const width, unsigned const height){

    // Initialize variables for the readablity
    int sA, sB, sC, sD, sF, sG, sH, sI;

    // Mirror on the edges
    unsigned rPrev = (r == 0)? r + 1 : r - 1;
    unsigned rNext = (r == height - 1)? r - 1 : r + 1;
    unsigned cPrev = (c == 0)? c + 1 : c - 1;
    unsigned cNext = (c == width - 1)? c - 1 : c + 1;

    // Read the prev row
    const std::vector<png::byte, std::allocator<png::byte>> rowPrev = image[rPrev];
    sA = (int)rowPrev[cPrev];
    sB = (int)rowPrev[c];
    sC = (int)rowPrev[cNext];
    // Read the current row
    const std::vector<png::byte, std::allocator<png::byte>> row = image[r];
    sD = (int)row[cPrev];
    // sE = row[c] (Not used)
    sF = (int)row[cNext];
    // Read the next row
    const std::vector<png::byte, std::allocator<png::byte>> rowNext = image[rNext];
    sG = (int)rowNext[cPrev];
    sH = (int)rowNext[c];
    sI = (int)rowNext[cNext];

    // Fill the sample arrays
    int x0[5] = {
        sA, sB, sD, sG, sH
    };
    int x1[5] = {
        sB, sC, sF, sH, sI
    };
    int y0[5] = {
        sA, sB, sC, sD, sF
    };
    int y1[5] = {
        sD, sF, sG, sH, sI
    };

    // NOTE: the typecast to float so that the gradient datatype is consistent
    float gx = (float)(median(x1, 5) - median(x0, 5));
    float gy = (float)(median(y1, 5) - median(y0, 5));

    // Package the gradient vector into a tuple and return
    tuple2<float> grad;
    grad._0 = gx / 2.0f;
    grad._1 = gy / 2.0f;
    return grad;
}

std::vector<std::vector<float>> generateGradientMap(png::image<png::gray_pixel> const &image, gradType const gradientType){

    // Image dimensions + Safety check
    unsigned const height = image.get_height();
    unsigned const width = image.get_width();
    if(0 == height || 0 == width){
        std::cerr << "\tWARNING! User tried to create gradient map on empty image." << std::endl;
    }

    std::vector<std::vector<float>> gradMap(height, std::vector<float>(width, 0));

    // Determine what function to use for the gradient calculation
    float (*computeGradientVector)(png::image<png::gray_pixel> const&, unsigned const, unsigned const, unsigned const, unsigned const);
    switch (gradientType){
    case SOBEL:
        computeGradientVector = &computeGradientVector_sobel;
        break;
    case MEDIAN:
        computeGradientVector = &computeGradientVector_median;
        break;
    default:
        std::cerr << "\tERROR! Invalid gradient type: "<<gradientType<<" passed to createGradientMap()." << std::endl;
        throw std::exception();
        break;
    }
        
    // Loop through th image pixels and compute gradient at each point
#pragma omp parallel for
    for (unsigned r = 0; r < height; r++){
        std::vector<float> gradMapRow = gradMap[r];
        for (unsigned c = 0; c < width; c++){
            gradMapRow[c] = computeGradientVector(image, r, c, width, height);
        }
        gradMap[r] = gradMapRow;
    }

    return gradMap;
}

std::vector<std::vector<tuple2<float>>> generateDirectionalGradientMap(png::image<png::gray_pixel> const &image, gradType const gradientType){

    // Image dimensions + Safety check
    unsigned const height = image.get_height();
    unsigned const width = image.get_width();
    if(0 == height || 0 == width){
        std::cerr << "\tWARNING! User tried to create gradient map on empty image." << std::endl;
    }

    std::vector<std::vector<tuple2<float>>> gradMap(height, std::vector<tuple2<float>>(width, {0,0}));

    // Determine what function to use for the gradient calculation
    tuple2<float> (*computeGradientVector)(png::image<png::gray_pixel> const&, unsigned const, unsigned const, unsigned const, unsigned const);
    switch (gradientType){
    case SOBEL:
        computeGradientVector = &computeDirectionalGradientVector_sobel;
        break;
    case MEDIAN:
        computeGradientVector = &computeDirectionalGradientVector_median;
        break;
    default:
        std::cerr << "\tERROR! Invalid gradient type: "<<gradientType<<" passed to createGradientMap()." << std::endl;
        throw std::exception();
        break;
    }
        
    // Loop through th image pixels and compute gradient at each point
#pragma omp parallel for
    for (unsigned r = 0; r < height; r++){
        std::vector<tuple2<float>> gradMapRow = gradMap[r];
        for (unsigned c = 0; c < width; c++){
            gradMapRow[c] = computeGradientVector(image, r, c, width, height);
        }
    }

    return gradMap;
}

/*
 * Document Cleaning
 */

tuple4<unsigned> findIAMTextBounds(png::image<png::gray_pixel> const &img){

    unsigned histThreshold = 180u;      // Minimum intensity for the histogram to count a pixel
    unsigned histLowerThreshold = 3;    // Minimum pixels for a row/col to be marked important
    unsigned cropPaddingWidth = 150;     // How much whitespace padding to add to the crop
    unsigned segmentRadius = 10;                     // Estimated segment Radius
    unsigned segmentDiameter = segmentRadius*2;     // Diameter calculation is also needed for min spacing

    unsigned width = img.get_width();
    unsigned height = img.get_height();

    // Generate a horizontal histogram of the image to find the segments
    std::vector<unsigned> histHorizontal = horizontalProjectionHistogram(img, histThreshold);
    
    // Find the max 3 peaks (these will be the segments)
    // Stored as (value, index)
    tuple2<int> yPeakList[3] = {{-INT32_MAX, -INT32_MAX}, {-INT32_MAX, -INT32_MAX}, {-INT32_MAX, -INT32_MAX}};
    for(int r = 0; r < height-segmentRadius; r++){
        int v = (int)histHorizontal[r];

        // Check that we arent next to a row we already marked
        // If so update the max value of the respective row if its greater
        bool nextToAPeak = false;
        for(unsigned i = 0; i < 3; i++){
            unsigned d = abs(yPeakList[i]._1 - r);
            if(d <= segmentDiameter){
                nextToAPeak = true;
                if(yPeakList[i]._0 < v){
                    yPeakList[i]._0 = v;
                    yPeakList[i]._1 = r;

                    // keep the list in order if we updated it
                    std::sort(yPeakList, yPeakList + 3, [](tuple2<int>a,tuple2<int>b){return a._0 > b._0;});
                }

                break;
            }
        }

        if(!nextToAPeak){
            for(unsigned i = 0; i < 3; i++){
                if(yPeakList[i]._0 < v){

                    // Inserting into ordered list, push everything up
                    for(unsigned ii = 2; i+1 <= ii; ii--){
                        yPeakList[ii] = yPeakList[ii-1];
                    }
                    yPeakList[i]._0 = v;
                    yPeakList[i]._1 = r;

                    break;
                }
            }
        }
    }
    
    // Saftey checks:
    // -- Should be three peaks (No negative values in the array)
    for(unsigned i = 0; i < 3; i++){
        tuple2<int> yPeak = yPeakList[i];
        if(yPeak._0 < 0 || yPeak._1 < 0){
            std::cerr << "\tERROR! findTextBounds() did not find all three segments." << std::endl;
            throw std::exception();
        }
    }

    // Sort the peaks in order of index
    std::sort(yPeakList, yPeakList+3, [](tuple2<int>a,tuple2<int>b){return a._1 < b._1;});

    // Return the bounds of the handwritten text 
    // NOTE: We overestimate segment width to be +/-5
    int y0 = yPeakList[1]._1 + segmentRadius;
    int y1 = yPeakList[2]._1 - segmentRadius;

    // Close the margins further by moving in to the next non-zero row
    // NOTE Crop is done 5 rows above/below the next non-zero bin, check that we went at least 5 rows in
    for(unsigned r = y0; r < y1; r++){
        unsigned v = histHorizontal[r];
        if(histLowerThreshold < v){
            int yy0 = r - cropPaddingWidth;
            if(y0 < yy0){
                y0 = yy0;
            }
            break;
        }
    }
    for(unsigned r = y1; y0 < r; r--){
        unsigned v = histHorizontal[r];
        if(histLowerThreshold < v){
            int yy1 = r + cropPaddingWidth;
            if(yy1 < y1){
                y1 = yy1;
            }
            break;
        }
    }

    // Determine the vertical bounds using the same out-in method above
    int x0 = 0;
    int x1 = width-1;
    std::vector<unsigned> histVertical = verticalProjectionHistogram(img, histThreshold);
    for(int c = x0; c < x1; c++){
        unsigned v = histVertical[c];
        if(histLowerThreshold < v){
            int xx0 = c - cropPaddingWidth;
            if(x0 < xx0){
                x0 = xx0;
            }
            break;
        }
    }
    for(int c = x1; x0 < c; c--){
        unsigned v = histVertical[c];
        if(histLowerThreshold < v){
            int xx1 = c + cropPaddingWidth;
            if(xx1 <= x1){
                x1 = xx1-1;
            }
            break;
        }
    }

    // Package and return the bounds
    tuple4<unsigned> textBounds = {(unsigned)y0, (unsigned)y1, (unsigned)x0, (unsigned)x1};
    return textBounds;
}

png::image<png::gray_pixel> cleanIAMDocument(png::image<png::gray_pixel> const &img){

    unsigned width = img.get_width();
    unsigned height = img.get_height();

    // Filter noise and isolate edges
    png::image<png::gray_pixel> edgeMapImage = edgeMapImg(img);

    // Compute the bounds on the handwritten text
    tuple4<unsigned> textBounds = findIAMTextBounds(edgeMapImage);
    unsigned y0 = textBounds._0;
    unsigned y1 = textBounds._1;
    unsigned x0 = textBounds._2;
    unsigned x1 = textBounds._3;

    // Crop the image into a new cropped image
    png::image<png::gray_pixel> imgYCropped = crop(img, x0, x1, y0, y1);

    return imgYCropped;
}

/*
 * Hit and Miss
 */

/*
 * Image Filtering
 */

png::image<png::gray_pixel> edgeMapImg(png::image<png::gray_pixel> const &img, gradType const gradientType, float const threshold){

    // Image dimensions
    unsigned const height = img.get_height();
    unsigned const width = img.get_width();

    // Fills the map with the (x,y) gradient at each pixel
    std::vector<std::vector<float>> gradMap = generateGradientMap(img, gradientType);

    // Threshold the gradient map
    png::image<png::gray_pixel> imgEdgeMap(width, height);
#pragma omp parallel for
    for(unsigned r = 0; r < height; r++){
        std::vector<float> rowG = gradMap[r];
        std::vector<png::gray_pixel> rowE = imgEdgeMap[r];
        for(unsigned c = 0; c < width; c++){
            if(threshold < rowG[c]) rowE[c] = 0;
            else rowE[c] = 255;
        }
        imgEdgeMap[r] = rowE;
    }

    return imgEdgeMap;
}

png::image<png::gray_pixel> erodeImg(png::image<png::gray_pixel> const& imgDoc, unsigned const kernelRadius){

    unsigned const height = imgDoc.get_height();
    unsigned const width = imgDoc.get_width();

    png::image<png::gray_pixel> imgEroded(width, height);

#pragma omp parallel for
    for(int y = 0; y < height; y++){
        std::vector<png::gray_pixel> rowTgt = imgEroded[y];
        for(int x = 0; x < width; x++){

            // Determine min pixel value in neighbourhood
            unsigned pixelMin = 255;
            for(int n = -(int)kernelRadius; n <= (int)kernelRadius; n++){
                int yy = y + n;
                if(0 <= yy && yy < height){
                    std::vector<png::gray_pixel> rowSrc = imgDoc[yy];
                    for(int m = -(int)kernelRadius; m <= (int)kernelRadius; m++){
                        int xx = x + m;
                        if(0 <= xx && xx < width){
                            unsigned pixel = rowSrc[xx];
                            if(pixel < pixelMin){
                                pixelMin = pixel;
                            }
                        }
                    }
                }
            }

            // Store min in center
            rowTgt[x] = pixelMin;
        }

        // Write row back
        imgEroded[y] = rowTgt;
    }

    return imgEroded;
}

png::image<png::gray_pixel> thinImg(png::image<png::gray_pixel> const& img){

    unsigned const maxRuns = 1000;  // Safety cap for the thinning process

    unsigned const height = img.get_height();
    unsigned const width = img.get_width();

    // Make a copy of the image that we will be editing and returning
    png::image<png::gray_pixel> imgCopy(img);

    // Subtract a hit&miss operation from the image on 8 rotations
    // See: https://homepages.inf.ed.ac.uk/rbf/HIPR2/thin.htm
    // Repeat this untill no changes are made
    bool changes;
    unsigned t = 0;
    do {
        changes = false;

        // For each pixel, check neighbours
#pragma omp parallel for
        for(int y = 1; y < height-1; y++){
            std::vector<png::gray_pixel> rowPrev = imgCopy[y-1];
            std::vector<png::gray_pixel>& row = imgCopy[y];
            std::vector<png::gray_pixel> rowNext = imgCopy[y+1];
            for(int x = 1; x < width-1; x++){
                if(0 == row[x]){
                    png::gray_pixel A = rowPrev[x-1];
                    png::gray_pixel B = rowPrev[x];
                    png::gray_pixel C = rowPrev[x+1];
                    png::gray_pixel D = row[x-1];
                    // png::gray_pixel E = row[x];
                    png::gray_pixel F = row[x+1];
                    png::gray_pixel G = rowNext[x-1];
                    png::gray_pixel H = rowNext[x];
                    png::gray_pixel I = rowNext[x+1];

                    // 0 degree
                    if( 0 < A && 0 < B && 0 < C && 
                        0 == G && 0 == H && 0 == I){
                        row[x] = 255;
                        changes = true;
                        continue;
                    }
                    // 45 degree
                    if(     0 < B && 
                        0 < D && 0 == F &&
                            0 == H){
                        row[x] = 255;
                        changes = true;
                        continue;
                    }
                    // 90 degree
                    if( 0 < A && 0 == C &&
                        0 < D && 0 == F &&
                        0 < G && 0 == I){
                        row[x] = 255;
                        changes = true;
                        continue;
                    }
                    // 135 degree
                    if(     0 == B &&
                        0 < D && 0 == F &&
                            0 < H){
                        row[x] = 255;
                        changes = true;
                        continue;
                    }
                    // 180 degree
                    if( 0 == A && 0 == B && 0 == C && 
                        0 < G && 0 < H && 0 < I){
                        row[x] = 255;
                        changes = true;
                        continue;
                    }
                    // 225 degree
                    if(      0 == B && 
                        0 == D && 0 < F &&
                             0 < H){
                        row[x] = 255;
                        changes = true;
                        continue;
                    }
                    // 270 degree
                    if( 0 == A && 0 < C &&
                        0 == D && 0 < F &&
                        0 == G && 0 < I){
                        row[x] = 255;
                        changes = true;
                        continue;
                    }
                    // 315 degree
                    if(      0 < B &&
                        0 == D && 0 < F &&
                             0 == H){
                        row[x] = 255;
                        changes = true;
                        continue;
                    }
                }
            }
        }

        t++;
    }while (changes && t < maxRuns);

    if(t == maxRuns){
        std::cerr << "\tERROR! Thinning process hit "<<maxRuns<<" itterations." << std::endl;
        throw std::exception();
    }

    return imgCopy;
}
