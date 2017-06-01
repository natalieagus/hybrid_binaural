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

#include "ISM.h"
#include "ISMRoom.h"

#include <stdio.h>
#include <iostream>
#include <fstream>
#include "BMMultiTapDelay.h"
#include "Cuboid.hpp"
#include "Gains.hpp"
#include "SingleTapDelay.h"
#include "FirstOrderFilter.h"

#endif /* defined(__offlineReverb__FDN__) */


//#define NUMDELAYS 16
#define NUMDELAYSEXT 1024
#define AUDIOCHANNELS 1
#define CHANNELS 12
#define NUMFDNDELAYS 64


class FDN
{
public:
    // constructor
    FDN();
    FDN(int type);
    
    void setDecayTime(double rt60);
 
    
protected:    
    
    //****ARE-RELATED****//
    size_t ARE_IRLength;
    size_t AREInit(size_t numberOfFirstOrder);




    
    //****ISM-RELATED****//
    ISMVector3D* imageSources;
    size_t totalISM;
    float* reflectionTimesL;
    float* reflectionTimesR;
    float* reflectionGainsL;
    float* reflectionGainsR;
    float* reflectionGains;
    size_t ISMInit(size_t reflectionMaxOrder);
    size_t ISM_IRLength;

    
    //****CHANNEL-RELATED****//
    size_t determineChannel(float x,float y, float orientation);
    size_t delayTimesChannel[NUMFDNDELAYS];//[ER_NUMTAPS];
    void setERDelayOutputChannels();
    void setFDNDelayOutputChannels();
    float determineAngle(float x, float y);
    size_t angleToChannel(float angleInDegrees);
    float channeltoangleNormal(size_t channel);
    float convertAzimuthToLeftEar(float azimuth);
    float convertAzimuthToRightEar(float azimuth);
    
    
    //****ITD-RELATED****//
    int earStatus[CHANNELS];
    void setTempPoints();
    Vector3D tempPoints[CHANNELS];
    void calculateAdditionalDelays();
    float additionalDelays[CHANNELS];
    SingleTapDelay reverbDelays[CHANNELS];
    void addReverbDelay(float* fdnLeft, float*fdnRight);
    
    
    //****HRTF-RELATED****//
    void setupHRTFFunctions();
    void processHRTFER(float ERTankOutL[CHANNELS], float ERTankOutR[CHANNELS], float* outputER);
    FirstOrderFilter leftEarFilter[CHANNELS];
    FirstOrderFilter rightEarFilter[CHANNELS];
    FirstOrderFilter leftEarFilterER[CHANNELS];
    FirstOrderFilter rightEarFilterER[CHANNELS];
    FirstOrderFilter directRayHRTFFilter[2];
    void setHRTFFilters();
    void processHRTFFDN(float* input, float* fdnTankOutLeft, float* fdnTankOutRight);
    
    //****MULTIPLEXER-RELATED****//
    void processTankOutER(float ERTankOut[CHANNELS], float* inputER);
    void processTankOutFDN(float FDNTankOut[CHANNELS], float* inputFDN);
    
    
    
    //****DIRECT-RAYS*****//
    SingleTapDelay directRays;
    float directRayAngles[2];
    SingleTapDelay directRaysStereo[2];
    void setDirectSingleTapDelay();
    void setDirectGains();
    void setDirectDelayTimes();
    void processDirectRays(float input, float* directRaysOutput);
    float directDelayTimes[2]; //unit = FREQ * seconds
    void HRTFFilterDirect(float* directRay, float* directRayAfterHRTF);
    void attenuateDirectRay(float* input);
    void setDirectRayAngles();
    float directRayAttenuation[2];



    
    //****ROOM PARAMETERS****//
    float orientation;
    Vector3D lLocLE = Vector3D();
    Vector3D lLocRE = Vector3D();
    Vector3D ssLoc;
    Vector3D lLoc;
    
    float RT60;
    float roomWidth;
    float roomHeight;
    float roomCeiling;
    Vector3D listenerLoc;
    Vector3D soundSourceLoc;
    
    Cuboid Room;
    Cuboid RoomARE;
    Gains GainValues;
    Gains GainValuesARE;
    void setGains();
    double time_elapsed_msecs = 0.0f;
    float fdnMultiplexInputAttenuation;
    float matrixAttenuation;

    
    //****FDN PROCESSING PARAMS****//
    
    float LRDivisionConstant =  1.f / sqrtf(2.f);
    float* leftEROut;
    float* rightEROut;
    float* leftEROutLastOrder;
    float* rightEROutLastOrder;
    
    float tapGains[NUMDELAYSEXT];
    float tapSigns[NUMDELAYSEXT];
    float temp1[NUMDELAYSEXT];
    float temp2[NUMDELAYSEXT];
    float temp3[NUMDELAYSEXT];
    float temp4[NUMDELAYSEXT];
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
    int rands, randomSeed;
    
public:
    void resetDelay();
    void processAudio(float* pInput, float* pOutputL, float* pOutputR);
    void impulseResponse(long numSamples, std::ofstream* outputFileL, std::ofstream* outputFileR);
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
