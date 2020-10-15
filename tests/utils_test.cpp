#include "../src/utils.h"
#include <iostream>
#include <array>
#include <fstream>

void getFileExtension_test(){

    std::cout << "Testing getFileExtension() ..." << std::endl;

    unsigned numTestFiles = 7; // NOTE: Tied to the following list
    char const* testFiles[7][2] = {
        {"", ""},
        {"dir", ""},
        {"file.ext", "ext"},
        {"file.extension", "extension"},
        {"file.d.ext", "ext"},
        {"../file.ext", "ext"},
        {"../dir", ""},
    };

    unsigned numTestDirs = 3; // NOTE: Tied to the following list
    std::array<std::string const, 3> testDirs = {
        "",
        "a/b/c/",
        "a/b/c.d/"
    };

    // Loop through each test case
    unsigned testCasesPassed = 0;
    unsigned testCasesTotal = numTestFiles * numTestDirs * 2;
    for(const auto& testDir : testDirs){
        if(0 < testDir.length()) fs::create_directories(testDir);

        // Load each test
        for(unsigned i = 0; i < numTestFiles; i++){
            char const** test = testFiles[i];
            std::string testFile(test[0]);
            std::string testOutput(test[1]);
            std::string testPath = testDir+testFile;

            // Create the test file in the OS
            // NOTE: path "" is in tests
            std::ofstream fp(testPath);
            if(fp.is_open() || testFile.length() == 0){
                fp.close();

                // Run test on both overloads
                std::string s = getFileExtension(testPath);
                if( testOutput.compare(s) == 0 ) testCasesPassed++;
                else std::cout << "Test failed on \""<<testPath<<"\" for strings. Expected \""<<testOutput<<"\" got \""<<s<<"\"" << std::endl;
                char const* cs = getFileExtension(testPath.c_str());
                if( testOutput.compare(cs) == 0 ) testCasesPassed++;
                else std::cout << "Test failed on \""<<testPath<<"\" for c_strs. Expected \""<<testOutput<<"\" got \""<<s<<"\"" << std::endl;
            }
            else{
                std::cerr << "ERROR! Test failed to open \""<<testPath<<"\"." << std::endl;
            }
        }
    }

    std::cout << testCasesPassed << " / "<< testCasesTotal <<" test cases passed." << std::endl;
    if (testCasesPassed == testCasesTotal) std::cout << "All tests passed." << std::endl;
    std::cout << std::endl;

    // Clean up the directories and files we made
    fs::remove_all("a/");
    for(unsigned i = 0; i < numTestFiles; i++){
        fs::remove(testFiles[i][0]);
    }
}

void getPathNoExtension_test(){

    std::cout << "Testing getPathNoExtension() ..." << std::endl;

    unsigned numTestFiles = 7; // NOTE: Tied to the following list
    char const* testFiles[7][2] = {
        {"", ""},
        {"dir", "dir"},
        {"file.ext", "file"},
        {"file.extension", "file"},
        {"file.d.ext", "file.d"},
        {"../file.ext", "../file"},
        {"../dir", "../dir"},
    };

    unsigned numTestDirs = 3; // NOTE: Tied to the following list
    std::array<std::string const, 3> testDirs = {
        "",
        "a/b/c/",
        "a/b/c.d/"
    };

    // Loop through each test case
    unsigned testCasesPassed = 0;
    unsigned testCasesTotal = numTestFiles * numTestDirs * 2;
    for(const auto& testDir : testDirs){
        if(0 < testDir.length()) fs::create_directories(testDir);

        // Load each test
        for(unsigned i = 0; i < numTestFiles; i++){
            char const** test = testFiles[i];
            std::string testFile(test[0]);
            std::string testOutput(test[1]);
            std::string testPath = testDir+testFile;
            testOutput = testDir+testOutput;

            // Create the test file in the OS
            // NOTE: path "" is in tests
            std::ofstream fp(testPath);
            if(fp.is_open() || testFile.length() == 0){
                fp.close();

                // Run test on both overloads
                std::string s = getPathNoExtension(testPath);
                if( testOutput.compare(s) == 0 ) testCasesPassed++;
                else std::cout << "Test failed on \""<<testPath<<"\" for strings. Expected \""<<testOutput<<"\" got \""<<s<<"\"" << std::endl;
                char const* cs = getPathNoExtension(testPath.c_str());
                if( testOutput.compare(cs) == 0 ) testCasesPassed++;
                else std::cout << "Test failed on \""<<testPath<<"\" for c_strs. Expected \""<<testOutput<<"\" got \""<<s<<"\"" << std::endl;
            }
            else{
                std::cerr << "ERROR! Test failed to open \""<<testPath<<"\"." << std::endl;
            }
        }
    }

    std::cout << testCasesPassed << " / "<< testCasesTotal <<" test cases passed." << std::endl;
    if (testCasesPassed == testCasesTotal) std::cout << "All tests passed." << std::endl;
    std::cout << std::endl;

    // Clean up the directories and files we made
    fs::remove_all("a/");
    for(unsigned i = 0; i < numTestFiles; i++){
        fs::remove(testFiles[i][0]);
    }
}

int main(int argc, char const *argv[]){

    getFileExtension_test();
    getPathNoExtension_test();
    
    return 0;
}