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

void median_test(){

    std::cout << "Testing median() ..." << std::endl;
    
    unsigned testCasesTotal = 0;
    unsigned testCasesPassed = 6;

    char test0[] = "";
    char output0 = '\0';
    char m0 = median(test0, 1);
    if(m0 == output0) testCasesPassed++;
    else std::cout << "Test failed on \"test0\". Expected \""<<output0<<"\" got \""<<m0<<"\"" << std::endl;

    char test1[] = "aebgfdc";
    char output1 = 'd';
    char m1 = median(test1, 8);
    if(m1 == output1) testCasesPassed++;
    else std::cout << "Test failed on \"test1\". Expected \""<<output1<<"\" got \""<<m1<<"\"" << std::endl;

    char test2[] = "addwdddga";
    char output2 = 'd';
    char m2 = median(test2, 10);
    if(m2 == output2) testCasesPassed++;
    else std::cout << "Test failed on \"test2\". Expected \""<<output2<<"\" got \""<<m2<<"\"" << std::endl;

    unsigned test3[] = {3,2,6,4,5,9,1,7,8};
    unsigned output3 = 5;
    unsigned m3 = median(test3, 9);
    if(m3 == output3) testCasesPassed++;
    else std::cout << "Test failed on \"test3\". Expected \""<<output3<<"\" got \""<<m3<<"\"" << std::endl;

    unsigned test4[] = {3,7,4,1,9,2,8,6};
    unsigned output4 = 6;
    unsigned m4 = median(test4, 8);
    if(m4 == output4) testCasesPassed++;
    else std::cout << "Test failed on \"test4\". Expected \""<<output4<<"\" got \""<<m4<<"\"" << std::endl;

    unsigned* test5 = (unsigned*)malloc(sizeof(unsigned)*5);
    if(test5){
        test5[0] = 0;
        test5[1] = 5;
        test5[2] = 1;
        test5[3] = 7;
        test5[4] = 7;
        unsigned output5 = 5;
        unsigned m5 = median(test5, 5);
        free(test5);
        if(m5 == output5) testCasesPassed++;
        else std::cout << "Test failed on \"test5\". Expected \""<<output5<<"\" got \""<<m5<<"\"" << std::endl;
    }
    else{
        std::cerr << "Failed to allocate memory for test5." << std::endl;
        throw std::exception();
    }

    std::cout << testCasesPassed << " / "<< testCasesTotal <<" test cases passed." << std::endl;
    if (testCasesPassed == testCasesTotal) std::cout << "All tests passed." << std::endl;
    std::cout << std::endl;
}

int main(int argc, char const *argv[]){

    getFileExtension_test();
    getPathNoExtension_test();
    median_test();

    return 0;
}