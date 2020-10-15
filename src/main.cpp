#include "image.h"
#include "utils.h"

#include <iostream>

int main(int argc, char const *argv[]){

    char const* addr_test = "data/document1.png";

    if(fs::exists(addr_test)){
        png::image<png::gray_pixel> img(addr_test);
        
        png::image<png::gray_pixel> imgSeg = segmentImage(img);

        imgSeg.write("data/document1_seg.png");
    }
    else{
        std::cout << "\tERROR! Failed to read image file \"" << addr_test << "\"." << std::endl;
    }

    return 0;
}
