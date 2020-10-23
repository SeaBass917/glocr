#include "../src/imagesegmentation.h"
#include "../src/utils.h"

void wordSegmentation_test(std::string testDir){

    std::string testDocumentPath = "document_cleaned.png";
    if(fs::exists(testDocumentPath)){
        std::string testDocumentPathNoExt = getPathNoExtension(testDocumentPath);

        png::image<png::gray_pixel> imgDoc(testDocumentPath);

        std::vector<std::vector<png::image<png::gray_pixel>>> lineImgs = wordSegmentation(imgDoc);

        // Double loop through the lines and images
        unsigned i = 0;
        for(auto& imgWords : lineImgs){
            unsigned j = 0;
            for(auto& imgWord : imgWords){
                imgWord.write(testDir+testDocumentPathNoExt+"_line_"+std::to_string(i)+"_word_"+std::to_string(j)+".png");

                j++;
            }
            i++;
        }
    }
    else{
        std::cerr << "\tERROR! preProcessing_test() cannot find \""<<testDocumentPath<<"\" needed for test." << std::endl;
    }
}

int main(int argc, char const *argv[]){

    // Create a directory for doing tests
    std::string testDir = "imagesegmentation_test/";
    if(!fs::exists(testDir)) fs::create_directory(testDir);

    // Test word segmentation in its own sub directory
    std::string testDir0 = testDir+"wordSegmentation_test/";
    if(!fs::exists(testDir0)) fs::create_directory(testDir0);
    wordSegmentation_test(testDir0);
    
    return 0;
}