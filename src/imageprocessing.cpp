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

tuple2<float> computeGradientVector_sobel(png::image<png::gray_pixel> const &image, unsigned const r, unsigned const c){

    // Image dimensions
    unsigned const height = image.get_height();
    unsigned const width = image.get_width();

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
    float gx = sC + 2.0*sF + sI - sA - 2.0*sD - sG;
    float gy = sG + 2.0*sH + sI - sA - 2.0*sB - sC;

    // Package the gradient vector into a tuple and return
    tuple2<float> grad;
    grad._0 = gx / 8.0f;
    grad._1 = gy / 8.0f;
    return grad;
}

tuple2<float> computeGradientVector_median(png::image<png::gray_pixel> const &image, unsigned const r, unsigned const c){

    // Image dimensions
    unsigned const height = image.get_height();
    unsigned const width = image.get_width();

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
            tuple2<float> (*computeGradientVector)(png::image<png::gray_pixel> const&, unsigned const, unsigned const);
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
            for (unsigned r = 0; r < height; r++){
                for (unsigned c = 0; c < width; c++){
                    tuple2<float> grad = computeGradientVector(image, r, c);
                    gradMap[r][c] = grad;
                }
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
            tuple2<float> (*computeGradientVector)(png::image<png::gray_pixel> const&, unsigned const, unsigned const);
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
            for (unsigned r = 0; r < height; r++){
                for (unsigned c = 0; c < width; c++){
                    tuple2<float> grad = computeGradientVector(image, r, c);
                    float gradX = grad._0;
                    float gradY = grad._1;
                    gradMap[r][c] = sqrtf(gradX*gradX + gradY*gradY);
                }
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

void drawEdgeMap(bool** edgeMap, unsigned height, unsigned width, char const* pathSmileEdge){
    png::image<png::gray_pixel> img(width, height);

    for(unsigned r = 0; r < height; r++){
        for(unsigned c = 0; c < width; c++){
            if(edgeMap[r][c]) img[r][c] = 255;
            else img[r][c] = 0;
        }
    }

    img.write(pathSmileEdge);
}
void drawEdgeMap(bool** edgeMap, unsigned height, unsigned width, std::string pathSmileEdge){
    png::image<png::gray_pixel> img(width, height);

    for(unsigned r = 0; r < height; r++){
        for(unsigned c = 0; c < width; c++){
            if(edgeMap[r][c]) img[r][c] = 255;
            else img[r][c] = 0;
        }
    }

    img.write(pathSmileEdge);
}

void horizontalProjectionHistogram(png::image<png::gray_pixel> const &img, double* hist, bool normalize){
    if(hist){

        // Img info
        unsigned height = img.get_height();
        unsigned width = img.get_width();

        unsigned histLeng = *(&hist + 1) - hist;
        if(height <= histLeng){

            // For normilization
            unsigned totalSum = 0;

            // For each row count pixels that are below the threshold
            unsigned thresh = 64;
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

            // Normalize the histogram if requested
            if(normalize){
                double normScale = 1.0 / (double)totalSum;
                for(unsigned r = 0; r < height; r++)
                    hist[r] *= normScale;
            }
        }
        else{
            std::cout << "\tERROR! horizontalProjectionHistogram() recieved a histogram smaller["<<histLeng<<"] than image height["<<height<<"]. Failed to make a histogram." << std::endl;
            std::cout << sizeof(hist) << std::endl;
            std::cout << sizeof(hist[0]) << std::endl;
        }
    }
    else{
        std::cout << "\tERROR! horizontalProjectionHistogram() recieved a null histogram. Failed to make a histogram." << std::endl;
    }
}

void verticalProjectionHistogram(png::image<png::gray_pixel> const &img, double* hist, bool normalize){
    if(hist){

        // Img info
        unsigned height = img.get_height();
        unsigned width = img.get_width();

        unsigned histLeng = *(&hist + 1) - hist;
        if(width <= histLeng){


            // For normilization
            unsigned totalSum = 0;

            // For each col count pixels that are below the threshold
            unsigned thresh = 127;
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

            // Normalize the histogram if requested
            if(normalize){
                double normScale = 1.0 / (double)totalSum;
                for(unsigned r = 0; r < width; r++)
                    hist[r] *= normScale;
            }
        }
        else{
            std::cout << "\tERROR! verticalProjectionHistogram() recieved a histogram smaller["<<histLeng<<"] than image width["<<width<<"]. Failed to make a histogram." << std::endl;
        }
    }
    else{
        std::cout << "\tERROR! verticalProjectionHistogram() recieved a null histogram. Failed to make a histogram." << std::endl;
    }
}

// -----------------
// Image processing
// -----------------

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
