// offlineReverb.cpp : Defines the entry point for the console application.
//

#include <stdio.h>
#include "FDN.h"
#include <ctime>
#include <string>
#include <fstream>
#include <cstdint>
#include "assert.h"

#define FS 44100

using std::cin;
using std::cout;
using std::endl;
using std::fstream;
using std::string;
using namespace std;

int processWaveInput(int type, long samples, std::ifstream* inFile, std::ofstream* outFile)
{
    // open the output file
    std::string outFilename = "./Users/natalieagus/Desktop/offlineReverb/offlineReverb/audioData/audio";
    if (type > 0) outFilename += "H";
    else outFilename += "C";
    outFilename += std::to_string(abs(type));
    outFilename += ".csv";
    outFile->open(outFilename);
    
    // open the input file
    std::string inFilename = "./Users/natalieagus/Desktop/offlineReverb/offlineReverb/audioData/arpSound.wav";
    inFile->open(inFilename, std::ifstream::binary);
    
    FDN reverb = FDN(type);
    reverb.fileResponse(samples, inFile, outFile);
    
    std::cout << "file output saved for type " << type << ".\n";
    
    outFile->close();
    inFile->close();
    
    return 0;
}

double timeReverb(int type, float decayTime, float seconds){
    FDN reverb = FDN(type);
    reverb.setDecayTime(decayTime);
    
    reverb.testProcess(1000);
    std::clock_t start = clock();
    reverb.testProcess(seconds * FS);
    std::clock_t end = clock();
    
    return double(end - start) / CLOCKS_PER_SEC;
}

void saveImpulse(int type, int samples, std::ofstream* of){
    
//    clock_t begin = clock();
    std::string filename = "impulse";
    if (type > 0) filename += "H";
    else filename += "C";
    filename += std::to_string(abs(type));
    filename += ".csv";
    
    
    of->open(filename);
    
    FDN reverb = FDN(type);
    reverb.impulseResponse(samples, of);
    
    std::cout << "impulse saved for type " << type << ".\n";
    of->close();
//    
//    clock_t end = clock();
//    double elapsed_msecs = double(end - begin) / CLOCKS_PER_SEC * 1000.f;
//    
//    printf("Time elapsed: %f ms\n", elapsed_msecs);
}

void saveDensity(int type, int samples, std::ofstream* of){
    std::string filename = "/Users/hans/Desktop/offlineReverb/impulseData/density";
    if (type > 0) filename += "H";
    else filename += "C";
    filename += std::to_string(abs(type));
    filename += ".csv";
    of->open(filename);
    
    FDN reverb = FDN(type);
    reverb.densityResponse(samples, of);
    
    std::cout << "density saved for type " << type << ".\n";
    of->close();
}

//x^y for x and y in int
int powi(int x, int y) {
    int result = 1;
    while (y-- > 0) result *= x;
    return result;
}

int main(int argc, char* argv[])
{
    float impulseLength = 5.0f;
    
    std::ofstream impulse;
    std::ifstream inFile;
    
    //float time = 500.0f; // seconds of reverb
    
    /*
    std::cout << "Timings for hadamard mixing matrices with " << time << " seconds of reverb\n{";
    std::cout.flush();
    for (int i = 0; i <= 5; i++){
        std::cout << "{" << 16 * powi(2,i) << ", " << timeReverb(16*powi(2,i), -1.0f, time) << "}, ";
        std::cout.flush();
    }
    std::cout << "\b\b}\n\n";
     */
    
    /*
    std::cout << "Timings for sparse block circulant mixing matrices with " << time << " seconds of reverb\n{";
    std::cout.flush();
    for (int i = 692; i <= 800; i += 32){
        std::cout << "{" << i << ", " << timeReverb(-i, -1.0f, time) << "}, ";
        std::cout.flush();
    }
    std::cout << "\b\b}\n\n";
    */
    
    /*
     save impules for H types
    for (int i = 0; i < 6; i++){
        int type = 16 * powi(2,i);
        saveImpulse(type, FS*impulseLength, &impulse);
    }
    */
    
//    for (int i = 0; i < 6; i++){
//        int type = 16 * powi(2,i);
//        printf("Tyles : %d \n", type);
////        saveImpulse(type, FS*impulseLength, &impulse);
//    }
    saveImpulse(16, FS*impulseLength, &impulse);
    // save impuses for C types
//    saveImpulse(-20, FS*impulseLength, &impulse);
    //saveImpulse(-52, FS*impulseLength, &impulse);
    //saveImpulse(-116, FS*impulseLength, &impulse);
    //saveImpulse(-248, FS*impulseLength, &impulse);
    //saveImpulse(-456, FS*impulseLength, &impulse);
    //saveImpulse(-848, FS*impulseLength, &impulse);
    
     
    /*
    // save densities for H types
    for (int i = 0; i < 5; i++){
        int type = 16 * powi(2,i);
        saveDensity(type, FS*impulseLength, &impulse);
    }
    
    // save densities for H types
    saveDensity(-20, FS*impulseLength, &impulse);
    saveDensity(-48, FS*impulseLength, &impulse);
    saveDensity(-104, FS*impulseLength, &impulse);
    saveDensity(-224, FS*impulseLength, &impulse);
    saveDensity(-368, FS*impulseLength, &impulse);
     */
    
//    processWaveInput(-128, FS*30, &inFile, &impulse);
//    processWaveInput(128, FS*30, &inFile, &impulse);
    
//    ofstream myfile ("example1234.txt");
//    if (myfile.is_open())
//    {
//
//        myfile << "This is a line.\n";
//        myfile << "This is another line.\n";
//        myfile.close();
//        printf("File is opened");
//
//    }
//    else if  (myfile.fail())
//    {
//        std::cout << "Failed to open outputfile.\n";
//    }
//    else cout << "Unable to open file";

    
    std::cout << "\ndone.\n";
        return 0;
}

