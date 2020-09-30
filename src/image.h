#ifndef IMAGE_H_DEFINED
#define IMAGE_H_DEFINED

#include <png.hpp>
#include <iostream>
#include <cmath> 

// Create a histogram of the horizontal projection of a grayscale image
void horizontalProjectionHistogram(png::image<png::gray_pixel> const &img, double* hist, bool normalize = true);

// Create a histogram of the vertical projection of a grayscale image
void verticalProjectionHistogram(png::image<png::gray_pixel> const &img, double* hist, bool normalize = true);

template<typename T = float>
T** gaussiankernel(unsigned const masklength, T const sigma);

template<typename T = float>
png::image<png::gray_pixel> conv2D(png::image<png::rgb_pixel> &image, unsigned const maskHeight, unsigned const maskWidth, T** kernel);

void drawASmile(char const* addrOut);

#endif