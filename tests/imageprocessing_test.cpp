#include "../src/imageprocessing.h"
#include "../src/utils.h"

#include <iostream>
#include <string>

int main(int argc, char const *argv[]){

    // Create a directory for doing tests
    std::string testDir = "imageprocessing_test/";
    fs::create_directory(testDir);

    // Draw a smile
    std::string pathSmile = testDir+"smile.png";
    drawASmile(pathSmile.c_str());

    std::string testImg = "document";

    // Address for test image edgemaps
    std::string pathEdgeSobel = testDir+testImg+"_edge-sobel.png";
    std::string pathEdgeMedian = testDir+testImg+"_edge-median.png";

    png::image<png::gray_pixel> img(testImg+".png");
    // png::image<png::gray_pixel> img(pathSmile);
    unsigned height = img.get_height();
    unsigned width = img.get_width();

    // Allocate and populate an edgemap
    bool** edgeMap = (bool**)malloc(sizeof(bool*) * height);
    if(edgeMap){
        for(unsigned i = 0; i < height; i++){
            bool* row = (bool*)malloc(sizeof(bool) * width);
            if(row){
                edgeMap[i] = row;
            }
            else{
                std::cerr << "Failed to allocate memory for edgemap." << std::endl;
                throw std::exception();
            }
        }

        populateEdgeMap(edgeMap, img, SOBEL, 40);
        drawEdgeMap(edgeMap, height, width, pathEdgeSobel);

        populateEdgeMap(edgeMap, img, MEDIAN, 40);
        drawEdgeMap(edgeMap, height, width, pathEdgeMedian);

        for(unsigned i = 0; i < height; i++)
            free(edgeMap[i]);
        free(edgeMap);
    }
    else{
        std::cerr << "Failed to allocate memory for edgemap." << std::endl;
        throw std::exception();
    }

    return 0;
}