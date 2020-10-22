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

// Segment the document image into a list of line images
std::vector<png::image<png::gray_pixel>> lineSegmentation(png::image<png::gray_pixel> const &imgDoc);

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
float computeGradientVector_sobel(png::image<png::gray_pixel> const &image, unsigned const r, unsigned const c, unsigned const width, unsigned const height);
tuple2<float> computeDirectionalGradientVector_sobel(png::image<png::gray_pixel> const &image, unsigned const r, unsigned const c, unsigned const width, unsigned const height);

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
float computeGradientVector_median(png::image<png::gray_pixel> const &image, unsigned const r, unsigned const c, unsigned const width, unsigned const height);
tuple2<float> computeDirectionalGradientVector_median(png::image<png::gray_pixel> const &image, unsigned const r, unsigned const c, unsigned const width, unsigned const height);

// Populates a directional gradient map based on the provided image
void populateDirectionalGradientMap(tuple2<float>** gradMap, png::image<png::gray_pixel> const &image, gradType const gradientType);

// Populates a gradient map based on the provided image
void populateGradientMap(tuple2<float>** gradMap, png::image<png::gray_pixel> const &image, gradType const gradientType);

// Populates a boolean edge map based on the provided image
void populateEdgeMap(bool** edgeMap, png::image<png::gray_pixel> const &image, gradType const gradientType, float const threshold);

// Draws an edgemap as a png
png::image<png::gray_pixel> imageFromEdgeMap(bool** edgeMap, unsigned const height, unsigned const width);

// Create a histogram of the horizontal projection of a grayscale image
void horizontalProjectionHistogram(unsigned* hist, png::image<png::gray_pixel> const &img, unsigned const thresh = 64);

// Create a histogram of the horizontal projection of a grayscale image (normalizes the output)
void horizontalProjectionHistogramNorm(double* hist, png::image<png::gray_pixel> const &img, unsigned const thresh = 64);

// Create a histogram of the vertical projection of a grayscale image
void verticalProjectionHistogram(unsigned* hist, png::image<png::gray_pixel> const &img, unsigned const thresh = 64);

// Create a histogram of the vertical projection of a grayscale image (normalizes the output)
void verticalProjectionHistogramNorm(double* hist, png::image<png::gray_pixel> const &img, unsigned const thresh = 64);

// -----------------
// Image processing
// -----------------

// Traces along the given midpoint in the given image and populates the yArray with 
// the y values needed to segment the image at the midpoint
void spaceTrace(unsigned* yArray, png::image<png::gray_pixel> const& img, unsigned const midPoint);

// Crops the image between two staggered y bounds
//  - Fills the empty space with whitespace
// USAGE:
// Requires two arrays that span the width and specify the top and bottom y values respecively
// (Optional) padding field that will add whitespace above and below the highest and lowest y points
png::image<png::gray_pixel> staggeredYCrop(png::image<png::gray_pixel> const& img, unsigned const* yTopArray, unsigned const* yBotArray, unsigned const yPadding = 0);

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
png::image<png::gray_pixel> segmentImageThreshold(png::image<png::gray_pixel> const &imgSrc, unsigned const T);

// Segment the image into black and white using Kittlers method
unsigned determineKittlerThreshold(png::image<png::gray_pixel> const &imgSrc);

// Function to create Gaussian kernel on length n
template<class T>
T** gaussiankernel(unsigned const masklength, T const sigma);

#endif