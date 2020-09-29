#ifndef IMAGE_H_DEFINED
#define IMAGE_H_DEFINED

#include <png.hpp>
#include <iostream>
#include <cmath> 

template<typename T = float>
T** gaussiankernel(unsigned const masklength, T const sigma);

template<typename T = float>
png::image<png::gray_pixel> conv2D(png::image<png::rgb_pixel> &image, unsigned const maskHeight, unsigned const maskWidth, T** kernel);

void drawASmile(char const* addrOut);

#endif