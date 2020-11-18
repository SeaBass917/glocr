#ifndef IMAGESEGMENTATION_H_DEFINED
#define IMAGESEGMENTATION_H_DEFINED

#include "png.hpp"
#include <map>

/*
 * IAM Utilities
 */

enum segQuality{
    PRT = 0,    // some lines correctly segmented
    ALL = 1     // all lines correctly segmented
};

std::map<std::string, std::array<unsigned, 7>> loadIAMDocumentsMetadata(std::string metaDataPath);

/*
 * Segmentation
 */

// Includes steps: Sobel edge marking, and image eroding
png::image<png::gray_pixel> preProcessDocumentImage(png::image<png::gray_pixel> const& imgDoc);

// Steps:
// a) edge detection
// b) skeletonization
// c) noise removal
png::image<png::gray_pixel> preProcessWordImage(png::image<png::gray_pixel> const& imgDoc);

// Returns a 2d list of word images from a given document
std::vector<std::vector<png::image<png::gray_pixel>>> wordSegmentation(png::image<png::gray_pixel> const &imgDoc);

// Returns a list of character images from the word image
std::vector<png::image<png::gray_pixel>> charSegmentation(png::image<png::gray_pixel> const &imgword);

/*
 * Tracing
 */

// Traces along the given midpoint in the given image and populates the yArray with 
// the y values needed to segment the image at the midpoint
std::vector<unsigned> spaceTraceHorizontal(png::image<png::gray_pixel> const& img, unsigned const yMin, unsigned const yMax);

// Traces along the given midpoint in the given image and populates the xArray with 
// the x values needed to segment the image at the midpoint
std::vector<unsigned> spaceTraceVertical(png::image<png::gray_pixel> const& img, unsigned const xMin, unsigned const xMax);

#endif