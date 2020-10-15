#ifndef IMAGE_H_DEFINED
#define IMAGE_H_DEFINED
/**************************************
 * Image handler for the glocr project
 * Currently just works with png.
 * TODO: Focus of this file is unclear
 *       Scope + split things that dont 
 *       need to be here
 */

#include <png.hpp>
#include <iostream>
#include <cmath> 

// Segment the image into black and white
png::image<png::gray_pixel> segmentImage(png::image<png::gray_pixel> const &img);

// Create a histogram of the horizontal projection of a grayscale image
void horizontalProjectionHistogram(png::image<png::gray_pixel> const &img, double* hist, bool normalize = true);

// Create a histogram of the vertical projection of a grayscale image
void verticalProjectionHistogram(png::image<png::gray_pixel> const &img, double* hist, bool normalize = true);

// Function to create Gaussian kernel on length n
template<typename T = float>
T** gaussiankernel(unsigned const masklength, T const sigma);

// Perform a 2d convolution on an image object
template<typename T = float>
png::image<png::gray_pixel> conv2D(png::image<png::rgb_pixel> &image, unsigned const maskHeight, unsigned const maskWidth, T** kernel);

// Draws a smile image to the address specified
void drawASmile(char const* addrOut);

#endif