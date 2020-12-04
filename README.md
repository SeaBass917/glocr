# glocr
______

## Description
This project will be building some of the necessary image analysis foundations to interpret and transmit grocery lists from a board to your phone. These foundations include the following operations on hand-written phrases: character segmentation, word/phrase segmentation, and character recognition.

##### Phase I (Word Segmentation)
Using the “documents” portion of the IAM database (see Figure 1) an algorithm will be constructed that takes in a handwritten document and outputs the sentences and words as image files. Word segmentation is the critical part here, but for analysis purposes it would be productive to output the line segmentation step. There will also be some pre-processing involved including but not limited to: isolating the written text from the full document and noise reduction.

##### Phase II (Character Segmentation)
The goal of this phase is to segment handwritten words into its character components. This will be tested on the “words” subsection of the IAM database (see Figure 3) and the resulting images from Phase I. Once provided an image of a word this algorithm will produce an ordered list of images for each of the characters in the word.

##### Phase III (Character Recognition)
Character recognition is often done using pretrained and heavy deep neural nets. Such nets are compute intensive and would not fit the scope of this project, since ideally this would be run on an embedded device. For that reason this phase will involved a simple algorithm like the one described in [5]. This would be trained on EMNIST dataset[6] or a derivative of this dataset.

## Requirements
 - [zlib](https://zlib.net/)
 - [libpng](http://www.libpng.org/pub/png/libpng.html)
 - [png++](https://www.nongnu.org/pngpp/)
 - [OpenMP](https://www.openmp.org/) (Optional, codebase is very unoptimized)

## Aquiring the IAM Handwritten Dataset
 - Visit [https://fki.tic.heia-fr.ch/databases/iam-handwriting-database](https://fki.tic.heia-fr.ch/databases/iam-handwriting-database).
 - Data is available here, but they will ask for name/institution/etc.
 - Just follow the instructions on the site.
 - This tool uses a cleaned version of these documents (see cleanIAMDocs tool for getting this data)

## Building (CMake)
Warning current CMake implementation has issues with OpenMP
Recommending vanilla build
 - `mkdir build`
 - `cd build`
 - `cmake ../`
 - `make`

## Building (Vanilla Linux)
 - Phase1_tests
   ```
   g++ -std=c++2a \
   -g tests/phase1_test.cpp src/imagesegmentation.cpp src/imageprocessing.cpp src/utils.cpp \
   -I ./include \
   -L ./lib -lz -lpng -fopenmp \
   -o phase1_tests
   ```
 - Phase2_tests
   ```
   g++ -std=c++2a \
   -g tests/phase2_test.cpp src/imagesegmentation.cpp src/imageprocessing.cpp src/utils.cpp \
   -I ./include \
   -L ./lib -lz -lpng -fopenmp \
   -o phase2_tests
   ```
 - ImageProcessing_tests
   ```
   g++ -std=c++2a \
   -g tests/imageprocessing_test.cpp src/imageprocessing.cpp src/utils.cpp \
   -I ./include \
   -L ./lib -lz -lpng -fopenmp \
   -o imageprocessing_tests
   ```
 - Utils_tests
   ```
   g++ -std=c++2a \
   -g tests/utils_test.cpp src/utils.cpp \
   -o utils_tests
   ```
 - CleanIAMDocs_tool
   ```
   g++ -std=c++2a \
   -g src/cleanIAMDocs.cpp src/imageprocessing.cpp src/utils.cpp \
   -o cleanIAMDocs_tool
   ```

## Usage

### Clean IAM Docs Tool
 - `./cleanIAMDocs_tool` Tool for cleaning the IAM docs

### Testing Phase 1
 - `./phase1_tests` For Simple test using the provided document
 - `./phase1_tests CLEANED_IAMDOCS_DIR IAMDOCS_METADATA_PATH OUTPUT_DIR [--clearOutput]` For Full test using a local copy of cleaned IAM dataset
 - `--clearOutput` If set will clear the output directory before running (otherwise program will continue from where it left off).

### Testing Phase 2
 - `./phase2_tests` For Simple test using the provided document
 - `./phase2_tests SEGMENTEDDOCS_DIR OUTPUT_DIR [--clearOutput]` For Full test using a local copy of Phase1's Full output
 - `--clearOutput` If set will clear the output directory before running (otherwise program will continue from where it left off).

### Testing Misc/Tools
 - `./imageprocessing_tests` Unit tests for Image Processing tools
 - `./utils_tests` Unit tests for Utils
