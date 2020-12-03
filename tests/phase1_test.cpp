#include "../src/imagesegmentation.h"
#include "../src/utils.h"

class perfReport{

public:

    // Parse a report from a file
    perfReport(std::ifstream &fp){

        // Parse the report
        std::string line;
        std::getline(fp, line);

        this->writerID = std::stoi(line.substr(9));
        std::getline(fp, line);
        this->isDocClean = (line.substr(11, 4) == "true")? true : false;
        std::getline(fp, line);
        this->numLines = std::stoi(line.substr(9));
        std::getline(fp, line);
        this->numWords = std::stoi(line.substr(9));
        std::getline(fp, line);
        this->numBadLines = std::stoi(line.substr(12));
        std::getline(fp, line);
        this->numBadWords = std::stoi(line.substr(12));
        std::getline(fp, line);
        this->lineDiff = std::stoi(line.substr(9));
        std::getline(fp, line);
        this->wordDiff = std::stoi(line.substr(9));

    }
    // generate a report from stats
    perfReport(unsigned const writerID, segQuality const docSegQuality, 
                int const numLines, int const numWords, int const numLinesGood, 
                int const numWordsGood, int const numLinesSegmented, int const numWordsSegmented){
        this->writerID = writerID;
        this->isDocClean = (docSegQuality == ALL)? true : false;
        this->numLines = numLines;
        this->numWords = numWords;
        this->numBadLines = numLines - numLinesGood;
        this->numBadWords = numWords - numWordsGood;
        this->lineDiff = numLinesSegmented - numLines;
        this->wordDiff = numWordsSegmented - numWords;
    }

    // To String
    // String is in ini file format
    std::string toString(){
        std::string perfString = "";
        perfString += "writerID=" + std::to_string(this->writerID) + "\n";
        perfString += (this->isDocClean)? "isDocClean=true\n" : "isDocClean=false\n";
        perfString += "numLines=" + std::to_string(this->numLines) + "\n";
        perfString += "numWords=" + std::to_string(this->numWords) + "\n";
        perfString += "numBadLines=" + std::to_string(this->numBadLines) + "\n";
        perfString += "numBadWords=" + std::to_string(this->numBadWords) + "\n";
        perfString += "lineDiff=" + std::to_string(this->lineDiff) + "\n";
        perfString += "wordDiff=" + std::to_string(this->wordDiff);// + "\n";
        return perfString;
    }

    unsigned writerID;
    bool isDocClean;
    int numLines;
    int numWords;
    int numBadLines;
    int numBadWords;
    int lineDiff;
    int wordDiff;

private:
};

void scoreWordSegmentation(std::string segmentedWordsDir, std::string outputPath){

    // Metrics to be computed
    // Note some documents are marked as unclean, we give clean docs their own isolated measurement
    unsigned numDocs = 0;
    unsigned numDocsClean = 0;
    unsigned numPerfectLineSegmentations = 0;
    unsigned numPerfectLineSegmentationsClean = 0;
    // unsigned numPerfectWordsSegmentations = 0;  //Cant calculate this, the words segmentation needs manual analysis
    // unsigned numPerfectWordsSegmentationsClean = 0; //Cant calculate this, the words segmentation needs manual analysis

    // build lists of the acc (#correct segment / #total)
    std::vector<float> lineSegmentationAccList;
    std::vector<float> lineSegmentationAccListClean;
    std::vector<float> wordSegmentationAccList;
    std::vector<float> wordSegmentationAccListClean;

    for(auto& docFolder : fs::directory_iterator(segmentedWordsDir)){
        std::string sDocFolder = docFolder.path();
        std::string sPerfReportPath = sDocFolder + "/perfReport.txt";

        // Open the performance report for this document
        if(fs::exists(sPerfReportPath)){
            std::ifstream fp(sPerfReportPath);
            if(fp.is_open()){
                perfReport rep(fp);

                bool isClean = rep.isDocClean;
                int lineDiff = rep.lineDiff;
                int wordDiff = rep.wordDiff;
                int numLines = rep.numLines;
                int numWords = rep.numWords;
                bool correctNumLineSegmentations = lineDiff == 0;

                numDocs++;
                if(isClean) numDocsClean++;
                
                if(correctNumLineSegmentations){
                    numPerfectLineSegmentations++;
                    if(isClean) numPerfectLineSegmentationsClean++;
                }

                float lineAcc = 1.0f - (float)abs(lineDiff) / (float)numLines;
                lineSegmentationAccList.push_back(lineAcc);
                if(isClean) lineSegmentationAccListClean.push_back(lineAcc);

                float wordAcc = 1.0f - (float)abs(wordDiff) / (float)numWords;
                wordSegmentationAccList.push_back(wordAcc);
                if(isClean) wordSegmentationAccListClean.push_back(wordAcc);

            }
            else{
                std::cerr << "\tWarning! failed to open performance report at \""<<sPerfReportPath<<"\"." << std::endl;
            }

        }
        else{
            std::cerr << "\tWarning! Missing a performance report for document \""<<docFolder<<"\"." << std::endl;
        }
    }

    std::string finalReport = "";
    finalReport += "# Docs:                                " + std::to_string(numDocs) + "\n";
    finalReport += "# Docs Clean:                          " + std::to_string(numDocsClean) + "\n";
    finalReport += "# Perfect Line Segmentations:          " + std::to_string(numPerfectLineSegmentations) + "\n";
    finalReport += "# Perfect Line Segmentations Clean:    " + std::to_string(numPerfectLineSegmentationsClean) + "\n";
    finalReport += "% Perfect Line Segmentations:          " + std::to_string((float)numPerfectLineSegmentations / (float) numDocs) + "%\n";
    finalReport += "% Perfect Line Segmentations Clean:    " + std::to_string((float)numPerfectLineSegmentationsClean / (float) numDocsClean) + "%\n\n";

    finalReport += "% Line Segmentation Avg:         " + std::to_string(avgv(lineSegmentationAccList)) + "%\n";
    finalReport += "% Line Segmentation Median:      " + std::to_string(medianv(lineSegmentationAccList)) + "%\n";
    finalReport += "% Line Segmentation Max:         " + std::to_string(maxv(lineSegmentationAccList)) + "%\n";
    finalReport += "% Line Segmentation Min:         " + std::to_string(minv(lineSegmentationAccList)) + "%\n\n";

    finalReport += "% Line Segmentation Clean Avg:    " + std::to_string(avgv(lineSegmentationAccListClean)) + "%\n";
    finalReport += "% Line Segmentation Clean Median: " + std::to_string(medianv(lineSegmentationAccListClean)) + "%\n";
    finalReport += "% Line Segmentation Clean Max:    " + std::to_string(maxv(lineSegmentationAccListClean)) + "%\n";
    finalReport += "% Line Segmentation Clean Min:    " + std::to_string(minv(lineSegmentationAccListClean)) + "%\n\n";

    finalReport += "% Word Segmentation Avg:         " + std::to_string(avgv(wordSegmentationAccList)) + "%\n";
    finalReport += "% Word Segmentation Median:      " + std::to_string(medianv(wordSegmentationAccList)) + "%\n";
    finalReport += "% Word Segmentation Max:         " + std::to_string(maxv(wordSegmentationAccList)) + "%\n";
    finalReport += "% Word Segmentation Min:         " + std::to_string(minv(wordSegmentationAccList)) + "%\n\n";

    finalReport += "% Word Segmentation Clean Avg:    " + std::to_string(avgv(wordSegmentationAccListClean)) + "%\n";
    finalReport += "% Word Segmentation Clean Median: " + std::to_string(medianv(wordSegmentationAccListClean)) + "%\n";
    finalReport += "% Word Segmentation Clean Max:    " + std::to_string(maxv(wordSegmentationAccListClean)) + "%\n";
    finalReport += "% Word Segmentation Clean Min:    " + std::to_string(minv(wordSegmentationAccListClean)) + "%\n\n";

    std::cout << finalReport;

    std::ofstream fp(outputPath);
    if(fp.is_open()){
        fp << finalReport;
        fp.close();
    }
    else{
        std::cerr << "\tError! Cannot access \""<<outputPath<<"\". Could not store final performance report." << std::endl; 
    }
}

void wordSegmentation_testfull(std::string documentsDir, std::string metaDataPath, std::string outputDir, bool clearOutput){

    // Create or clear our working directory
    if(!fs::exists(outputDir)) fs::create_directory(outputDir);
    else if(clearOutput){
        fs::remove_all(outputDir);
        fs::create_directory(outputDir);
    }

    // Count the number of documents so we can inform user
    unsigned numDocs = 0;
    for(auto& docFile : fs::directory_iterator(documentsDir))
        numDocs++;

    // Load the metadata on the documents for validation of output
    std::map<std::string, std::array<unsigned, 7>> metaData = loadIAMDocumentsMetadata(metaDataPath);

    // Loop through each document
    unsigned iDoc = 1;
    for(auto& docFile : fs::directory_iterator(documentsDir)){
        std::string sDocFile = docFile.path().filename();
        std::string docname = getPathNoExtension(sDocFile);
        std::string documentPath = documentsDir + "/" + sDocFile;
        std::string segmentedOutputDir = outputDir + "/" + docname + "/";
        std::string perfReportPath = segmentedOutputDir + "perfReport.txt";

        std::cout << "Segmenting words in document: \""<<docname<<"\". ("<<iDoc<<"/"<<numDocs<<")" << std::endl;
        iDoc++;

        // Read in the info for this document
        if(metaData.contains(docname)){
            std::array<unsigned, 7> info = metaData[docname];
            unsigned const writerID = info[0];
            segQuality const docSegQuality = (segQuality)(info[2]);
            unsigned const numLines = info[3];
            unsigned const numLinesGood = info[4];
            unsigned const numWords = info[5];
            unsigned const numWordsGood = info[6];

            // Create an output directory for this document
            // If it already exists check for a performance report
            // If its there we have alreaady processed this document
            if(!fs::exists(segmentedOutputDir)) fs::create_directory(segmentedOutputDir);
            else if(fs::exists(perfReportPath)) continue;

            // Segment this document
            png::image<png::gray_pixel> imgDoc(documentPath);
            try{
                std::vector<std::vector<png::image<png::gray_pixel>>> wordImgs = wordSegmentation(imgDoc);

                // Save the segmented images in their own directory
                // include a performance report in that directory
                unsigned numWordsSegmented = 0;
                unsigned numLinesSegmented = 0;
                unsigned i = 0;
                for(auto& imgWords : wordImgs){
                    unsigned j = 0;
                    for(auto& imgWord : imgWords){
                        imgWord.write(segmentedOutputDir+"line_"+std::to_string(i)+"-word_"+std::to_string(j)+".png");

                        numWordsSegmented++;
                        j++;
                    }
                    numLinesSegmented++;
                    i++;
                }

                // Save performance report
                std::ofstream fp(perfReportPath);
                if(fp){
                    perfReport rep(writerID, docSegQuality, numLines, numWords, numLinesGood, numWordsGood, numLinesSegmented, numWordsSegmented);
                    std::string perfString = rep.toString();
                    fp << perfString;
                    std::cout << perfString << std::endl;
                    fp.close();
                }
                else{
                    std::cerr << "\tWARNING! Cannot access \""<<perfReportPath<<"\". Could not store performance report." << std::endl; 
                }
            }
            catch(const std::exception& e)
            {
                std::cerr << "Caught an exception while segmenting words." << std::endl;
                std::cerr << e.what() << std::endl;
            }
        }
        else{
            std::cerr << "\tWARNING! Cannot find document \""<<docname<<"\" in metadata file. Skipping." << std::endl; 
        }
    }
}

void wordSegmentation_testsimple(){
    std::string documentPath = "document_cleaned.png";
    std::string outputDir = "Phase1_SimpleTest/";
    std::string perfReportPath = outputDir+"score.txt";
    if(fs::exists(documentPath)){

        // Create output directory
        if(!fs::exists(outputDir)) fs::create_directory(outputDir);

        // Load and segment the test word
        png::image<png::gray_pixel> imgDoc(documentPath);
        try{
            std::cout << "Segmenting words in test document..." << std::endl;
            std::vector<std::vector<png::image<png::gray_pixel>>> wordImgs = wordSegmentation(imgDoc);

            // Save the segmented images in their own directory
            // include a performance report in that directory
            std::cout << "Saving word images locally." << std::endl;
            unsigned numWordsSegmented = 0;
            unsigned numLinesSegmented = 0;
            unsigned i = 0;
            for(auto& imgWords : wordImgs){
                unsigned j = 0;
                for(auto& imgWord : imgWords){
                    imgWord.write(outputDir+"line_"+std::to_string(i)+"-word_"+std::to_string(j)+".png");

                    numWordsSegmented++;
                    j++;
                }
                numLinesSegmented++;
                i++;
            }

            // Save performance report
            std::ofstream fp(perfReportPath);
            if(fp){
                // a01-011u 000 2 all 10 10 68 68
                perfReport rep(0, ALL, 10, 68, 10, 68, numLinesSegmented, numWordsSegmented);
                std::string perfString = rep.toString();
                fp << perfString;
                fp.close();

                std::cout << "Performance:" << std::endl;
                std::cout << "\tLine Accuracy: " << 1.0f - (float)abs(rep.lineDiff) / (float)rep.numLines << std::endl;
                std::cout << "\tWord Accuracy: " << 1.0f - (float)abs(rep.wordDiff) / (float)rep.numWords << std::endl;
            }
            else{
                std::cerr << "\tWARNING! Cannot access \""<<perfReportPath<<"\". Could not store performance report." << std::endl; 
            }
        }
        catch(const std::exception& e)
        {
            std::cerr << "Caught an exception while segmenting words." << std::endl;
            std::cerr << e.what() << std::endl;
        }
    }
    else{
        std::cerr << "\tERROR! Cannot find document_cleaned.png file that came with the project." << std::endl;
    }
}

void printHelp(){
    std::cerr << "\tERROR! Bad call to phase 1 tests." << std::endl;
    std::cerr << "Usage:" << std::endl << std::endl;
    std::cerr << "Simple test:" << std::endl;
    std::cerr << "\t./phase1" << std::endl;
    std::cerr << "Full Test:" << std::endl;
    std::cerr << "\t./phase1 CLEANED_IAMDOCS_DIR IAMDOCS_METADATA_PATH OUTPUT_DIR [--clearOutput]" << std::endl;
    std::cerr << "\t\t--clearOutput : On True: Will wipe the output directory instead of continuing from last run point" << std::endl;
}

int main(int argc, char const *argv[]){

    // Full test
    if(4 == argc || 5 == argc){
        // On True: Will wipe the output directory instead of continuing from last run point
        bool clearOutput = (5 == argc)? "--clearOutput" == argv[4] : false;
        
        std::string documentsDir = argv[1];
        std::string metaDataPath = argv[2];
        std::string outputDir = argv[3];
        std::string outputPath = outputDir+"/full-score.txt";

        if(fs::exists(documentsDir)){
            if(fs::exists(metaDataPath)){

                wordSegmentation_testfull(documentsDir, metaDataPath, outputDir, clearOutput);
                scoreWordSegmentation(outputDir, outputPath);

            }
            else{
                std::cerr << "\tERROR! Cannnot access IAMDOCS_METADATA_PATH \""<<metaDataPath<<"\"." << std::endl;
            }
        }
        else{
            std::cerr << "\tERROR! Cannnot access CLEANED_IAMDOCS_DIR \""<<documentsDir<<"\"." << std::endl;
        }
    }
    // Simple Test, using local file
    else if(argc == 1){
        wordSegmentation_testsimple();
    }
    else{
        // BAD ARGS
        printHelp();
    }

    return 0;
}