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

void wordSegmentation_test(std::string documentsDir, std::string metaDataPath, std::string outputDir, bool clearOutput){

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
                std::vector<std::vector<png::image<png::gray_pixel>>> lineImgs = wordSegmentation(imgDoc);

                // Save the segmented images in their own directory
                // include a performance report in that directory
                unsigned numWordsSegmented = 0;
                unsigned numLinesSegmented = 0;
                unsigned i = 0;
                for(auto& imgWords : lineImgs){
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

// Loop through each document
// Make a folder for the document
// Segment the characters into individual word sub directories
// outputdir/doc0/w0/c0.png
// outputdir/doc0/w0/c1.png
// outputdir/doc0/w1/c0.png
// ...
// outputdir/doc1/w0/c0.png
// ...
void charSegmentation_test(std::string documentsDir, std::string outputDir, bool clearOutput){
    
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
                    if(toLower(getFileExtension(wordFilename)) == "png"){

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
                        }

                        i_word++;
                    }
                }
            }
        }
    }

}

int main(int argc, char const *argv[]){

    // Test word segmentation in its own sub directory
    std::string documentsDir = "/home/seabass/extra/IAM/documentsCleaned_small";
    std::string metaDataPath = "/home/seabass/extra/IAM/metadata/documents.txt";
    std::string output0Dir = "/home/seabass/extra/IAM/documentsSegmented_small";
    bool clearOutput = false;   // On True: Will wipe the output directory instead of continuing from last run point
    // wordSegmentation_test(documentsDir, metaDataPath, output0Dir, clearOutput);
    
    std::string outputPath = "/home/seabass/extra/IAM/documentsSegmented_small/score.txt";
    // scoreWordSegmentation(outputDir, outputPath);

    std::string wordsDir = "/home/seabass/extra/IAM/documentsSegmented_small";
    std::string output1Dir = "/home/seabass/extra/IAM/charsSegmented_small";
    charSegmentation_test(wordsDir, output1Dir, false);

    return 0;
}