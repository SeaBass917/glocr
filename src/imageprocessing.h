#ifndef IMAGEPROCESSING_H_DEFINED
#define IMAGEPROCESSING_H_DEFINED

#include "utils.h"

#include "png.hpp"

enum gradType{
    SOBEL,
    MEDIAN
};

// -----------------
// Image Misc
// -----------------

// Loads a png file from the specified location
// Will read as a grayscale 8-bit image
png::image<png::gray_pixel> loadPNG(char const* addr);

// Draws a smile image to the address specified
void drawASmile(char const* addrOut);

// -----------------
// Image Analysis
// -----------------

/*
         [-1][0][1]       [-1][-2][-1]  
    h1 = [-2][0][2]  h2 = [ 0][ 0][ 0]  scaling factor: 1/8
         [-1][0][1]       [ 1][ 2][ 1]

    [A][B][C] 
    [D][E][F] 
    [G][H][I] 

    Returns: 
        Gx, Gy
*/
tuple2<float> computeGradientVector_sobel(png::image<png::gray_pixel> const &image, unsigned const r, unsigned const c);

/*
    𝐼𝑥(𝑟, 𝑐) = median{𝐼(𝑟 − 1, 𝑐),𝐼(𝑟 − 1, 𝑐 + 1),𝐼(𝑟, 𝑐 + 1),𝐼(𝑟 + 1, 𝑐),𝐼(𝑟 + 1, 𝑐 + 1)}
                    − median{𝐼(𝑟 − 1, 𝑐 − 1),𝐼(𝑟 − 1, 𝑐),𝐼(𝑟, 𝑐 − 1),𝐼(𝑟 + 1, 𝑐 − 1),𝐼(𝑟 + 1, 𝑐)}
                    
    𝐼𝑦(𝑟, 𝑐) = median{𝐼(𝑟, 𝑐 − 1),𝐼(𝑟, 𝑐 + 1),𝐼(𝑟 − 1, 𝑐 − 1),𝐼(𝑟 − 1, 𝑐),𝐼(𝑟 − 1, 𝑐 + 1)}
                    − median{𝐼(𝑟, 𝑐 − 1),𝐼(𝑟, 𝑐 + 1),𝐼(𝑟 + 1, 𝑐 − 1),𝐼(𝑟 + 1, 𝑐),𝐼(𝑟 + 1, 𝑐 + 1)}

    [A][B][C]   [A][B]   [B][C]    [A][B][C]          
    [D][E][F]   [D]         [F]    [D]   [F]   [D]   [F]
    [G][H][I]   [G][H]   [H][I]                [G][H][I]
                    x0       x1          y0          y1

    𝐼𝑥(𝑟, 𝑐) = median(x1) - median(x0)
    𝐼𝑦(𝑟, 𝑐) = median(y0) - median(y1)

    Returns: 
        Gx, Gy
*/
tuple2<float> computeGradientVector_median(png::image<png::gray_pixel> const &image, unsigned const r, unsigned const c);

// Populates a directional gradient map based on the provided image
void populateDirectionalGradientMap(tuple2<float>** gradMap, png::image<png::gray_pixel> &image, gradType gradientType);

// Populates a gradient map based on the provided image
void populateGradientMap(tuple2<float>** gradMap, png::image<png::gray_pixel> &image, gradType gradientType);

// Populates a boolean edge map based on the provided image
void populateEdgeMap(bool** edgeMap, png::image<png::gray_pixel> const &image, gradType gradientType, float threshold);

// Draws an edgemap as a local png
void drawEdgeMap(bool** edgeMap, unsigned height, unsigned width, char const* pathSmileEdge);
void drawEdgeMap(bool** edgeMap, unsigned height, unsigned width, std::string pathSmileEdge);

// Create a histogram of the horizontal projection of a grayscale image
void horizontalProjectionHistogram(png::image<png::gray_pixel> const &img, unsigned* hist, unsigned thresh = 64);

// Create a histogram of the horizontal projection of a grayscale image (normalizes the output)
void horizontalProjectionHistogramNorm(png::image<png::gray_pixel> const &img, double* hist, unsigned thresh = 64);

// Create a histogram of the vertical projection of a grayscale image
void verticalProjectionHistogram(png::image<png::gray_pixel> const &img, unsigned* hist, unsigned thresh = 64);

// Create a histogram of the vertical projection of a grayscale image (normalizes the output)
void verticalProjectionHistogramNorm(png::image<png::gray_pixel> const &img, double* hist, unsigned thresh = 64);

// -----------------
// Image processing
// -----------------

// Find the yMin and yMax values for the handwritten part of an IAM dataset text document
// NOTE: Cuttoff is made at the line segments, there are 2 on the top and one on the bottom
//       The text will be located between the top 2 and bottom
//           -------------------
//           Computer text
//           -------------------
//           *Handwritten text*
//           *Handwritten text*
//           *Handwritten text*
//           -------------------
tuple4<unsigned> findTextBounds(png::image<png::gray_pixel> const &img);

// Cleans a document from the IAM dataset
// - Crops out everything but the centered text
png::image<png::gray_pixel> preProcessDocument(png::image<png::gray_pixel> const &img);

// Segment the image into black and white using a fixed threshold
png::image<png::gray_pixel> segmentImageThreshold(png::image<png::gray_pixel> const &imgSrc, unsigned T);

// Segment the image into black and white using Kittlers method
png::image<png::gray_pixel> segmentImageKittler(png::image<png::gray_pixel> const &imgSrc);

// Function to create Gaussian kernel on length n
template<class T>
T** gaussiankernel(unsigned const masklength, T const sigma);

#endif