#include "imageprocessing.h"

#include <iostream>
#include <cmath>

// -----------------
// Image Misc
// -----------------

png::image<png::gray_pixel> loadPNG(char const* addr){
    png::image<png::gray_pixel> img(addr);
    return img;
}

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

// -----------------
// Image Analysis
// -----------------

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

void populateDirectionalGradientMap(tuple2<float>** gradMap, png::image<png::gray_pixel> const &image, gradType gradientType){

    // Image dimensions
    unsigned const height = image.get_height();
    unsigned const width = image.get_width();
    if(0 < height && 0 < width){
        if(gradMap){

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
            for (unsigned r = 0; r < height; r++){
                tuple2<float>* gradMapRow = gradMap[r];
                for (unsigned c = 0; c < width; c++){
                    tuple2<float> grad = computeGradientVector(image, r, c, width, height);
                    gradMapRow[c] = grad;
                }
                gradMap[r] = gradMapRow;
            }
        }
        else{
            std::cerr << "\tERROR! Null gradient map passed to createGradientMap()." << std::endl;
            throw std::exception();
        }
    }
    else{
        std::cerr << "\tWARNING! User tried to create gradient map on empty image." << std::endl;
    }
}

void populateGradientMap(float** gradMap, png::image<png::gray_pixel> const &image, gradType gradientType){

    // Image dimensions
    unsigned const height = image.get_height();
    unsigned const width = image.get_width();
    if(0 < height && 0 < width){
        if(gradMap){

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
                float* gradMapRow = gradMap[r];
                for (unsigned c = 0; c < width; c++){
                    gradMapRow[c] = computeGradientVector(image, r, c, width, height);
                }
                gradMap[r] = gradMapRow;
            }
        }
        else{
            std::cerr << "\tERROR! Null gradient map passed to createGradientMap()." << std::endl;
            throw std::exception();
        }
    }
    else{
        std::cerr << "\tWARNING! User tried to create gradient map on empty image." << std::endl;
    }
}

void populateEdgeMap(bool** edgeMap, png::image<png::gray_pixel> const &image, gradType gradientType, float threshold){
    if(edgeMap){

        // Image dimensions
        unsigned height = image.get_height();
        unsigned width = image.get_width();

        // Allocate gradient map
        float** gradMap = (float**)malloc(sizeof(float*) * height);
        if(gradMap){
            for(unsigned i = 0; i < height; i++){
                float* row = (float*)malloc(sizeof(float) * width);
                if(row){
                    gradMap[i] = row;
                }
                else{
                    std::cerr << "Failed to allocate memory for gradientmap." << std::endl;
                    throw std::exception();
                }
            }

            // Fills the map with the (x,y) gradient at each pixel
            populateGradientMap(gradMap, image, gradientType);

            // Threshold the gradient map
#pragma omp parallel for
            for(unsigned r = 0; r < height; r++){
                float* rowG = gradMap[r];
                bool* rowE = edgeMap[r];
                for(unsigned c = 0; c < width; c++){
                    if(threshold < rowG[c]) rowE[c] = true;
                    else rowE[c] = false;
                }
            }

            // Clean up the gradient map
            for(unsigned i = 0; i < height; i++)
                free(gradMap[i]);
            free(gradMap);
        }
        else{
            std::cerr << "createEdgeMap() recieved null gradientmap." << std::endl;
            throw std::exception();
        }
    }
    else{
        std::cerr << "createEdgeMap() recieved null edgemap." << std::endl;
        throw std::exception();
    }
}

void horizontalProjectionHistogram(png::image<png::gray_pixel> const &img, unsigned* hist, unsigned thresh){
    if(hist){

        // Img info
        unsigned height = img.get_height();
        unsigned width = img.get_width();

        // For each row count pixels that are below the threshold
        for(unsigned r = 0; r < height; r++){

            // Read row once
            const std::vector<png::byte, std::allocator<png::byte>> row = img[r];

            unsigned sum = 0;
            for(unsigned c = 0; c < width; c++){
                if(row[c] < thresh){
                    sum++;
                }
            }

            hist[r] = sum;
        }
    }
    else{
        std::cout << "\tERROR! horizontalProjectionHistogram() recieved a null histogram. Failed to make a histogram." << std::endl;
    }
}

void horizontalProjectionHistogramNorm(png::image<png::gray_pixel> const &img, double* hist, unsigned thresh){
    if(hist){

        // Img info
        unsigned height = img.get_height();
        unsigned width = img.get_width();

        // For normilization
        unsigned totalSum = 0;

        // For each row count pixels that are below the threshold
        for(unsigned r = 0; r < height; r++){

            // Read row once
            const std::vector<png::byte, std::allocator<png::byte>> row = img[r];

            unsigned sum = 0;
            for(unsigned c = 0; c < width; c++){
                if(row[c] < thresh){
                    sum++;
                    totalSum++;
                }
            }

            hist[r] = sum;
        }

        // Normalize the histogram
        double normScale = 1.0 / (double)totalSum;
        for(unsigned r = 0; r < height; r++)
            hist[r] *= normScale;
    }
    else{
        std::cout << "\tERROR! horizontalProjectionHistogram() recieved a null histogram. Failed to make a histogram." << std::endl;
    }
}

void verticalProjectionHistogram(png::image<png::gray_pixel> const &img, unsigned* hist, unsigned thresh){
    if(hist){

        // Img info
        unsigned height = img.get_height();
        unsigned width = img.get_width();

        // For each col count pixels that are below the threshold
        for(unsigned c = 0; c < width; c++){

            unsigned sum = 0;
            for(unsigned r = 0; r < height; r++){
                if(img[r][c] < thresh){
                    sum++;
                }
            }

            hist[c] = sum;
        }
    }
    else{
        std::cout << "\tERROR! verticalProjectionHistogram() recieved a null histogram. Failed to make a histogram." << std::endl;
    }
}

void verticalProjectionHistogramNorm(png::image<png::gray_pixel> const &img, double* hist, unsigned thresh){
    if(hist){

        // Img info
        unsigned height = img.get_height();
        unsigned width = img.get_width();

        // For normilization
        unsigned totalSum = 0;

        // For each col count pixels that are below the threshold
        for(unsigned c = 0; c < width; c++){

            unsigned sum = 0;
            for(unsigned r = 0; r < height; r++){
                if(img[r][c] < thresh){
                    sum++;
                    totalSum++;
                }
            }

            hist[c] = sum;
        }

        // Normalize the histogram
        double normScale = 1.0 / (double)totalSum;
        for(unsigned r = 0; r < width; r++)
            hist[r] *= normScale;
    }
    else{
        std::cout << "\tERROR! verticalProjectionHistogram() recieved a null histogram. Failed to make a histogram." << std::endl;
    }
}

// -----------------
// Image processing
// -----------------

bool sortTupleValue(tuple2<int> a, tuple2<int> b){
    return a._0 > b._0;
}

bool sortTupleIndex(tuple2<int> a, tuple2<int> b){
    return a._1 < b._1;
}

tuple4<unsigned> findTextBounds(png::image<png::gray_pixel> const &img){

    unsigned histThreshold = 180u;      // Minimum intensity for the histogram to count a pixel
    unsigned histLowerThreshold = 3;    // Minimum pixels for a row/col to be marked important
    unsigned cropPaddingWidth = 150;     // How much whitespace padding to add to the crop
    unsigned segmentRadius = 10;                     // Estimated segment Radius
    unsigned segmentDiameter = segmentRadius*2;     // Diameter calculation is also needed for min spacing

    unsigned width = img.get_width();
    unsigned height = img.get_height();

    // Generate a horizontal histogram of the image to find the segments
    unsigned* histHorizontal = (unsigned*) malloc(height * sizeof(unsigned));
    unsigned* histVertical = (unsigned*) malloc(width * sizeof(unsigned));
    if(histHorizontal && histVertical){
        horizontalProjectionHistogram(img, histHorizontal, histThreshold);

        // DEBUG
        // for(unsigned i = 0; i < height; i++)
        //     std::cout << i << ',' << histHorizontal[i] << std::endl;
        
        // Find the max 3 peaks (these will be the segments)
        // Stored as (value, index)
        tuple2<int> yPeakList[3] = {{-INT32_MAX, -INT32_MAX}, {-INT32_MAX, -INT32_MAX}, {-INT32_MAX, -INT32_MAX}};
        for(int r = 0; r < height; r++){
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
                        std::sort(yPeakList, yPeakList + 3, sortTupleValue);
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
        std::sort(yPeakList, yPeakList+3, sortTupleIndex);

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
        int x1 = width;
        verticalProjectionHistogram(img, histVertical, histThreshold);
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
                if(xx1 < x1){
                    x1 = xx1;
                }
                break;
            }
        }

        // Free and return the bounds 
        free(histHorizontal);
        free(histVertical);
        tuple4<unsigned> textBounds = {(unsigned)y0, (unsigned)y1, (unsigned)x0, (unsigned)x1};
        return textBounds;
    }
    else{
        std::cerr << "\tERROR! Failed to allocated histogram needed in findTextBounds()." << std::endl;
        throw std::exception();
    }
}

png::image<png::gray_pixel> isolateHandwrittenText(png::image<png::gray_pixel> const &img, png::image<png::gray_pixel> const &edgeMapImage){

    unsigned width = img.get_width();
    unsigned height = img.get_height();

    // Compute the bounds on the handwritten text
    tuple4<unsigned> textBounds = findTextBounds(edgeMapImage);
    unsigned y1 = textBounds._0;
    unsigned y2 = textBounds._1;
    unsigned x1 = textBounds._2;
    unsigned x2 = textBounds._3;

    // Crop the image into a new cropped image
    unsigned heightCropped = y2 - y1;
    unsigned widthCropped = x2 - x1;
    png::image<png::gray_pixel> imgYCropped(widthCropped, heightCropped);
    for(unsigned r = y1, rr = 0; r < y2; r++, rr++){
        std::vector<png::gray_pixel> const rowSrc = img[r];
        std::vector<png::gray_pixel> rowTgt = imgYCropped[rr];
        for(unsigned c = x1, cc = 0; c < x2; c++, cc++){
            rowTgt[cc] = rowSrc[c];
        }
        imgYCropped[rr] = rowTgt;
    }

    return imgYCropped;
}

png::image<png::gray_pixel> convertToEdgeMap(png::image<png::gray_pixel> const &img){
    unsigned const height = img.get_height();
    unsigned const width = img.get_width();

    bool** edgeMap = (bool**)malloc(sizeof(bool*) * height);
    if(edgeMap){
        for(unsigned i = 0; i < height; i++){
            bool* row = (bool*)malloc(sizeof(bool) * width);
            if(row){
                edgeMap[i] = row;
            }
            else{
                std::cerr << "Failed to allocate memory for edgemap." << std::endl;
                throw std::exception();
            }
        }

        // Create boolean edge map, convert to png image
        populateEdgeMap(edgeMap, img, SOBEL, 60);
        png::image<png::gray_pixel> imgEdgeMap = imageFromEdgeMap(edgeMap, height, width);

        for(unsigned i = 0; i < height; i++)
            free(edgeMap[i]);
        free(edgeMap);

        return imgEdgeMap;
    }
    else{
        std::cerr << "Failed to allocate memory for edgemap." << std::endl;
        throw std::exception();
    }
}

png::image<png::gray_pixel> imageFromEdgeMap(bool** edgeMap, unsigned const height, unsigned const width){
    png::image<png::gray_pixel> img(width, height);

    for(unsigned r = 0; r < height; r++){
        bool* rowEdge = edgeMap[r];
        std::vector<png::gray_pixel> rowImg = img[r];
        for(unsigned c = 0; c < width; c++){
            if(rowEdge[c]) rowImg[c] = 0;
            else rowImg[c] = 255;
        }
        img[r] = rowImg;
    }

    return img;
}

png::image<png::gray_pixel> preProcessDocument(png::image<png::gray_pixel> const &img){

    // Filter noise and isolate edges
    png::image<png::gray_pixel> edgeMapImage = convertToEdgeMap(img);

    // DEBUG
    // edgeMapImage.write("../edge.png");

    // Isolate the handwritten text
    png::image<png::gray_pixel> imgYCropped = isolateHandwrittenText(img, edgeMapImage);
    
    return imgYCropped;
}

png::image<png::gray_pixel> segmentImageThreshold(png::image<png::gray_pixel> const &imgSrc, unsigned T){

    // Image dimensions
    unsigned height = imgSrc.get_height();
    unsigned width = imgSrc.get_width();

    // Create and return the segmented image
    png::image< png::gray_pixel > imgSeg(width, height);
    for(unsigned y = 0; y < height; y++){
        for(unsigned x = 0; x < width; x++){
            if(T < imgSrc[y][x]) imgSeg[y][x] = 255;
            else imgSeg[y][x] = 0;
        }   
    }

    return imgSeg;
}

png::image<png::gray_pixel> segmentImageKittler(png::image<png::gray_pixel> const &imgSrc){

    // Image dimensions
    unsigned height = imgSrc.get_height();
    unsigned width = imgSrc.get_width();
    unsigned numPixels = height * width;

    // Determine the threshold for segmenting the image
    // Using Kittler's Method

    // Generate grayscale histogram from image source
    double hist[256] = {0};
    for(unsigned y = 0; y < height; y++){
        for(unsigned x = 0; x < width; x++){
            hist[ imgSrc[y][x] ]++;
        }
    }

    // Normalize the histogram
    for(unsigned i = 0; i < 256; i++){
        hist[i] /= numPixels;
    }

    // Compute the error for each threshold
    double errMin = std::numeric_limits<double>::max();
    unsigned T = 0;
    double errConst = 0.5+0.5*log(2.0*M_PI); // Const used in error calculations
    for(unsigned t = 0; t < 256; t++){
        double q1 = 0;
        for(unsigned i = 0; i <= t; i++){
            q1 += hist[i];
        }

        double q2 = 1.0 - q1;

        double m1 = 0;
        for(unsigned i = 0; i <= t; i++){
            m1 += i * hist[i] / q1;
        }

        double m2 = 0;
        for(unsigned i = t+1; i < 256; i++){
            m2 += i * hist[i] / q1;
        }

        double v1 = 0;
        double diff;
        for(unsigned i = 0; i <= t; i++){
            diff = i - m1;
            v1 += diff * diff * hist[i] / q1;
        }

        double v2 = 0;
        for(unsigned i = t+1; i < 256; i++){
            diff = i - m2;
            v2 += diff * diff * hist[i] / q2;
        }

        double err = errConst - q1*log(q1) - q2*log(q2) + q1*log(v1)/2.0 + q2*log(v2)/2.0;

        if(err < errMin){
            errMin = err;
            T = t;
        }
    }

    std::cout << "Threshold: " << T << std::endl;

    // Create and return the segmented image
    png::image< png::gray_pixel > imgSeg(width, height);
    for(unsigned y = 0; y < height; y++){
        for(unsigned x = 0; x < width; x++){
            if(T < imgSrc[y][x]) imgSeg[y][x] = 0;
            else imgSeg[y][x] = 255;
        }   
    }

    return imgSeg;
}

template<class T>
T** gaussiankernel(unsigned const masklength, T const sigma){

    T** kernel;

    // Mask length must be odd
    if(masklength & 1 == 1){

        const int radius = (masklength - 1) / 2;

        // intialising standard deviation
        T var = 2.0 * sigma * sigma;
        T denom = M_PI * var;
    
        // Allocating kernel 
        kernel = (T**)malloc(masklength * sizeof(T*));
        for (unsigned i = 0; i < masklength; i++){
            kernel[i] = (T*)malloc(masklength * sizeof(T));
        }

        T sum = 0.0; // sum is for normalization 
        for (int x = -radius; x <= radius; x++) {
            T rx = x * x;
            for (int y = -radius; y <= radius; y++) {
                T Gxy = exp(-(rx + y * y) / var) / denom;
                kernel[x + radius][y + radius] = Gxy;
                sum += Gxy; 
            } 
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
