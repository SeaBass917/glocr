#include "../src/imageprocessing.h"
#include "../src/utils.h"

#include <iostream>
#include <string>

// TODO: Automate test a little bit
void edgeDetection_test(std::string testDir){

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
        png::image<png::gray_pixel> edgeMapImage = imageFromEdgeMap(edgeMap, height, width);
        edgeMapImage.write(pathEdgeSobel);

        populateEdgeMap(edgeMap, img, MEDIAN, 40);
        edgeMapImage = imageFromEdgeMap(edgeMap, height, width);
        edgeMapImage.write(pathEdgeMedian);

        for(unsigned i = 0; i < height; i++)
            free(edgeMap[i]);
        free(edgeMap);
    }
    else{
        std::cerr << "Failed to allocate memory for edgemap." << std::endl;
        throw std::exception();
    }
}

void preProcessing_test(std::string testDir){

    // Check the document preprocessor

    // Test image
    std::string testImg = "../data/IAM/documents/a02-000.png";

    if(fs::exists(testImg)){

        png::image<png::gray_pixel> imgDoc(testImg);

        png::image<png::gray_pixel> imgDocClean = preProcessDocument(imgDoc);

        std::string fileOut = testDir+"document_clean.png";
        imgDocClean.write(fileOut);
    }
    else{
        std::cerr << "\tERROR! preProcessing_test() cannot find \""<<testImg<<"\" needed for test." << std::endl;
    }
}

void lineSegmentation_test(std::string testDir){

    std::string testDocumentPath = "document_cleaned.png";
    if(fs::exists(testDocumentPath)){

        png::image<png::gray_pixel> imgDoc(testDocumentPath);

        std::vector<png::image<png::gray_pixel>> lineImgs = lineSegmentation(imgDoc);

        std::string outputFilename = getPathNoExtension(testDocumentPath) + "_line_";
        unsigned i = 0;
        for(auto& lineImg : lineImgs){
            std::string fileOut = testDir+outputFilename+std::to_string(i)+".png";
            lineImg.write(fileOut);
            i++;
        }
    }
    else{
        std::cerr << "\tERROR! preProcessing_test() cannot find \""<<testDocumentPath<<"\" needed for test." << std::endl;
    }
}

int main(int argc, char const *argv[]){

    // Create a directory for doing tests
    std::string testDir = "imageprocessing_test/";
    if(!fs::exists(testDir)) fs::create_directory(testDir);

    // edgeDetection_test(testDir);
    // preProcessing_test(testDir);
    lineSegmentation_test(testDir);

    return 0;
}