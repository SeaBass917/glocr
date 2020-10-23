#ifndef IMAGESEGMENTATION_H_DEFINED
#define IMAGESEGMENTATION_H_DEFINED

#include "png.hpp"

/*
 * Segmentation
 */

// Includes steps: Sobel edge marking, and image eroding
png::image<png::gray_pixel> preProcessDocumentImage(png::image<png::gray_pixel> const& imgDoc);

// Returns a 2d list of word images from a given document
std::vector<std::vector<png::image<png::gray_pixel>>> wordSegmentation(png::image<png::gray_pixel> const &imgDoc);

/*
 * Tracing
 */

// Traces along the given midpoint in the given image and populates the yArray with 
// the y values needed to segment the image at the midpoint
std::vector<unsigned> spaceTraceHorizontal(png::image<png::gray_pixel> const& img, unsigned const midPoint);

// Traces along the given midpoint in the given image and populates the xArray with 
// the x values needed to segment the image at the midpoint
std::vector<unsigned> spaceTraceVertical(png::image<png::gray_pixel> const& img, unsigned const midPoint);

#endif