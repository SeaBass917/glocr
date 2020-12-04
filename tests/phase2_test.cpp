#include "../src/imagesegmentation.h"
#include "../src/utils.h"

// Loop through each document
// Make a folder for the document
// Segment the characters into individual word sub directories
// outputdir/doc0/w0/c0.png
// outputdir/doc0/w0/c1.png
// outputdir/doc0/w1/c0.png
// ...
// outputdir/doc1/w0/c0.png
// ...
void charSegmentation_testfull(std::string documentsDir, std::string outputDir, bool clearOutput){
    
    // Create or clear our working directory
    if(!fs::exists(outputDir)) fs::create_directory(outputDir);
    else if(clearOutput){
        fs::remove_all(outputDir);
        fs::create_directory(outputDir);
    }

    // Count the number of documents so we can inform user
    unsigned numDocs = 0;
    for(auto& docDir : fs::directory_iterator(documentsDir))
        if(fs::is_directory(docDir)) numDocs++;

    // Loop through each document
    unsigned iDoc = 1;
    for(auto& docDirPath : fs::directory_iterator(documentsDir)){
        if(fs::is_directory(docDirPath)){
            std::string docname = docDirPath.path().filename();
            std::string documentOutputDir = outputDir + "/" + docname + "/";
            std::cout << "Segmenting characters in document: \""<<docname<<"\". ("<<iDoc<<"/"<<numDocs<<")" << std::endl;
            iDoc++;

            // Check if the output has already been computed (if it has, then skip this document)
            bool isEmpty = true;
            if((fs::exists(documentOutputDir))) {
                unsigned cnt = 0;
                for(auto& f : fs::recursive_directory_iterator(documentOutputDir))
                    if(fs::is_regular_file(f))
                        cnt++;
                isEmpty = 0 == cnt;
            }
            else fs::create_directories(documentOutputDir); // NOTE: Make the directory if its not there
            if(isEmpty){

                // For each segmented word in the document
                unsigned i_word = 0;
                for(auto& wordFilePath : fs::directory_iterator(docDirPath)){
                    std::string wordFilename = wordFilePath.path().filename();
                    std::string wordname = getPathNoExtension(wordFilename);    // Name associated with the word (e.g. line_0-word_1)
                    std::string wordPath = wordFilePath.path().string();
                    std::string wordOutputDir = documentOutputDir + '/' + wordname + '/';
                    std::string ext = toLower(getFileExtension(wordFilename));
                    if(ext == "png"){

                        //Load the word image
                        png::image<png::gray_pixel> img(wordPath);

                        // Segment the word
                        std::vector<png::image<png::gray_pixel>> charImgs = charSegmentation(img);

                        // Store the resulting characters in our output directory
                        unsigned i_char = 0;
                        fs::create_directory(wordOutputDir);
                        for(png::image<png::gray_pixel>& charImg : charImgs){
                            std::string charOutputPath = wordOutputDir+"char_"+std::to_string(i_char)+".png";
                            charImg.write(charOutputPath);
                            i_char++;
                        }

                        i_word++;
                    }
                }
            }
        }
    }

}

void charSegmentation_testsimple(){
    std::string wordPath = "data/word.png";
    std::string outputDir = "Phase2_SimpleTest/";
    std::string perfReportPath = outputDir+"score.txt";
    if(fs::exists(wordPath)){

        // Create output directory
        if(!fs::exists(outputDir)) fs::create_directory(outputDir);

        //Load the word image
        png::image<png::gray_pixel> img(wordPath);

        // Segment the word
        std::vector<png::image<png::gray_pixel>> charImgs = charSegmentation(img);

        // Store the resulting characters in our output directory
        unsigned i_char = 0;
        for(png::image<png::gray_pixel>& charImg : charImgs){
            charImg.write(outputDir+"char_"+std::to_string(i_char)+".png");
            i_char++;
        }
        
        // Scoring
        int numChars = charImgs.size();
        int const numExpectedChars = 11;
        std::cout << "Estimated Accuracy: " << 1.0f - abs(numChars-numExpectedChars) / (float)numExpectedChars << std::endl;
    }
    else{
        std::cerr << "\tERROR! Cannot find word.png file that came with the project." << std::endl;
    }
}

void printHelp(){
    std::cerr << "\tERROR! Bad call to phase 2 tests." << std::endl;
    std::cerr << "Usage:" << std::endl << std::endl;
    std::cerr << "Simple test:" << std::endl;
    std::cerr << "\t./phase2" << std::endl;
    std::cerr << "Full Test:" << std::endl;
    std::cerr << "\t./phase2 SEGMENTEDDOCS_DIR OUTPUT_DIR [--clearOutput]" << std::endl;
    std::cerr << "\t\t--clearOutput : On True: Will wipe the output directory instead of continuing from last run point" << std::endl;
}

int main(int argc, char const *argv[]){

    // Full test
    if(3 == argc || 4 == argc){
        // On True: Will wipe the output directory instead of continuing from last run point
        bool clearOutput = (4 == argc)? std::string("--clearOutput").compare(argv[3]) : false;
        
        std::string wordsDir = argv[1];
        std::string outputDir = argv[2];

        if(fs::exists(wordsDir)){
            charSegmentation_testfull(wordsDir, outputDir, clearOutput);
        }
        else{
            std::cerr << "\tERROR! Cannnot access SEGMENTEDDOCS_DIR \""<<wordsDir<<"\"." << std::endl;
        }
    }
    // Simple Test, using local file
    else if(argc == 1){
        charSegmentation_testsimple();
    }
    else{
        // BAD ARGS
        printHelp();
    }

    return 0;
}