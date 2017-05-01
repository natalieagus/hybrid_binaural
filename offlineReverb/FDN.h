//
//  FDN.h
//
//   A Feedback Delay Network Reverberator
//
//  offlineReverb
//
//  Created by Hans on 19/5/15.
//  Copyright (c) 2015 Hans. All rights reserved.
//

#ifndef __offlineReverb__FDN__
#define __offlineReverb__FDN__

#include <stdio.h>
#include <iostream>
#include <fstream>
#include "BMMultiTapDelay.h"
#include "Cuboid.hpp"
#include "Gains.hpp"
#include "SingleTapDelay.h"
#include "FirstOrderFilter.h"

#endif /* defined(__offlineReverb__FDN__) */


#define NUMDELAYS 16
#define NUMDELAYSEXT 1024
#define AUDIOCHANNELS 1
#define CHANNELS 8


class FDN
{
public:
    // constructor
    FDN();
    FDN(int type);
    void setDecayTime(double rt60);
    ~FDN();
    
protected:
    
    
    float outputsBinaural[TOTALDELAYS];
    
    size_t delayTimesChannel[ER_NUMTAPS];
    
    Vector3D lLocLE = Vector3D();
    Vector3D lLocRE = Vector3D();
    void processTankOut(float fdnTankOut[CHANNELS]);
    
    //To handle direct Rays
    SingleTapDelay directRaysStereo[2];
    void setDirectSingleTapDelay();
    void setDirectGains();
    void setDirectDelayTimes();
    void processDirectRays(float* input, float* directRaysOutput);
    float directDelayTimes[2]; //unit = FREQ * seconds
    void processDirectRays(float* input, float* directRaysOutput);
    
    //Binaural filters
    FirstOrderFilter leftEarFilter[CHANNELS];
    FirstOrderFilter rightEarFilter[CHANNELS];
    FirstOrderFilter directRayFilter[2];
    void setHRTFFilters();
    //Direct ray is inclusive in one of the channels, so there's no need to have another channel angle for this
    float directRayAngles[2];
    //process the fdntankout[channels]
    void filterChannels(float fdnTankOut[CHANNELS], float directRay[2], float fdnTankOutLeft[CHANNELS], float fdnTankOutRight[CHANNELS]);

    //To do channel angle calculations
    void setDelayChannels();
    void setDirectRayAngles();
    size_t determineChannel(float x, float y);
    size_t angleToChannel(float angleInDegrees);
    float channelToAngle(size_t channel);
    float channeltoangleNormal(size_t channel);
    

    void  setMultiTapDelayChannels();
    size_t multiDelayLinePointsChannel[ER_NUMTAPS];
    
    double time_elapsed_msecs = 0.0f;
    bool diffuserSameAsERLength;
    
    int firstState = 0;
    
    float directRayAttenuation;
    
    SingleTapDelay directRays;
    BMMultiTapDelay earlyReflector;
    size_t erDelayTapTimes[ER_NUMTAPS];
    float delayTapGains[ER_NUMTAPS];

    float multiTapOutputs[ER_NUMTAPS];


    float totalRt60Energy;

    BMMultiTapDelay diffuser;
    size_t diffuserTapTimes[DIFFUSER_NUMTAPS];
    
    Vector3D ssLoc;
    Vector3D lLoc;
    Cuboid Room;
    Gains GainValues;
    void setGains();
    
//    float inputGains[MULTITAPDELAYS];
//    float outputGains[MULTITAPDELAYS];
    
    float minDim;
    float maxDim;
    
    float gamma = 1.0f;
    float upsilon = 1.0f;
    
    void processAudioNoGain(float* pInput, float* pOutputL, float* pOutputR);
    
//    float* inputAttenuationArray;
//    float* outputAttenuationArray;
//
    
    float RT60GainsForMultiTapDelay[DIFFUSER_NUMTAPS];
    
    float fdnMultiplexInputAttenuation;
//    float outputAttenuation;
    float matrixAttenuation;
    float wetPct,dryPct;
    float un_matrixAttenuation;
    float tapGains[NUMDELAYSEXT];
    float tapSigns[NUMDELAYSEXT];
    float temp1[NUMDELAYSEXT];
    float temp2[NUMDELAYSEXT];
    float temp3[NUMDELAYSEXT];
    float inputTemp[NUMDELAYSEXT];
    float temp16[16];
    float ppxx_xn, xxpp, pnxx_xn, xxpn;
    
    float* delayBuffers;
    
    // read / write indices
    size_t rwIndices[NUMDELAYSEXT];
    size_t endIndices[NUMDELAYSEXT];
    size_t startIndices[NUMDELAYSEXT];
    size_t totalDelayTime;
    
    // delay times
    size_t delayTimes[NUMDELAYSEXT];
    
    //float maxRand_f;
    
    float outputs[NUMDELAYSEXT];
    
    void resetReadIndices();
    void incrementIndices();
    void randomPermutation(size_t* a, size_t length, size_t audioChannels);
    void randomPermutation1Channel(size_t* a, size_t length, size_t audioChannels);
    void hadamardTransform(float* input, float* output, size_t length);
    void hadamardRecursion(float* input, float* output, size_t length);
    void hadamardIteration(float* input, float* output, size_t length);
    void hadamard16(float* input, float* output);
    void velvetNoiseTapTimes(float minDelaySamples,
                             float maxDelaySamples,
                             size_t* tapTimes,
                             size_t numTaps);
    
    void updateRand();
    int rvType, numDelays;
    int rand, randomSeed;
public:
    void resetDelay();
    void processAudio(float* pInput, float* pOutputL, float* pOutputR);
    void testProcess(long numSamples);
    void impulseResponse(long numSamples, std::ofstream* outputFile);
    void densityResponse(long numSamples, std::ofstream* outputFile);
    void fileResponse(long numSamples, std::ifstream* inputFile, std::ofstream* outputFile);
};



typedef struct  WAV_HEADER
{
    /* RIFF Chunk Descriptor */
    uint8_t         RIFF[4];        // RIFF Header Magic header
    uint32_t        ChunkSize;      // RIFF Chunk Size
    uint8_t         WAVE[4];        // WAVE Header
    /* "fmt" sub-chunk */
    uint8_t         fmt[4];         // FMT header
    uint32_t        Subchunk1Size;  // Size of the fmt chunk
    uint16_t        AudioFormat;    // Audio format 1=PCM,6=mulaw,7=alaw,     257=IBM Mu-Law, 258=IBM A-Law, 259=ADPCM
    uint16_t        NumOfChan;      // Number of channels 1=Mono 2=Sterio
    uint32_t        SamplesPerSec;  // Sampling Frequency in Hz
    uint32_t        bytesPerSec;    // bytes per second
    uint16_t        blockAlign;     // 2=16-bit mono, 4=16-bit stereo
    uint16_t        bitsPerSample;  // Number of bits per sample
    /* "data" sub-chunk */
    uint8_t         Subchunk2ID[4]; // "data"  string
    uint32_t        Subchunk2Size;  // Sampled data length
} wav_hdr;

union byteTo16 {
    int16_t sixteenBit;
    uint8_t bytes[2];
};
