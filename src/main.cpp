#include "image.h"
#include "utils.h"

#include <iostream>

int main(int argc, char const *argv[]){

    char const* addr_test = "data/document.png";

    if(pathExists(addr_test)){
        png::image<png::gray_pixel> img(addr_test);

        unsigned height = img.get_height();
        unsigned width = img.get_width();

        std::cout << "width  " << width  << std::endl;
        std::cout << "height " << height << std::endl;

        double* hist = (double*)malloc(height * sizeof(double));
        horizontalProjectionHistogram(img, hist, false);

        for(unsigned r = 0; r < height; r++)
            std::cout << hist[r] << std::endl;

        free(hist);
    }
    else{
        std::cout << "\tERROR! Failed to read image file \"" << addr_test << "\"." << std::endl;
    }

    return 0;
}
