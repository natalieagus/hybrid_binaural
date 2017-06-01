//
//  FDN.cpp
//  offlineReverb
//
//  Created by Hans on 19/5/15.
//  Copyright (c) 2015 Hans. All rights reserved.
//

#include "FDN.h"
#include <cstdlib>
#include <cstring>
#include <math.h>
#include <assert.h>
#include <random>

#include "ISM.h"
#include "ISMRoom.h"
#include "ERImpulseResponse.h"

#include "HRTFFilter.h"

#include <stdio.h>
#include <stdlib.h>

#define ISMORDER 3
#define DENSITY_WINDOW_SIZE 882 // 20.0 * (44100.0 / 1000.0); (20 ms)
//#define FDN_RT60 0.39f
//#define GAINRT60 FDN_RT60
#define FDN_SAMPLE_RATE 44100.f
#define AREFIRSTORDER 64

#define FDN_DFSR_MIN_DLY 1000
#define FDN_DFSR_MAX_DLY 4000
#define RADIUSOFHEAD 0.08f
#define DMIN 1.0f

//todo:
//do the ISM without diffuser, just pass to FDN straight
//integrate ARE with the new HRTF method

//***CHOOSE ONE MODE***//
//#define FDN_ER_SECOND_ORDER true
#define USE_ISM  //For ISM. Otherwise use ARE for order 2




inline void FDN::setGains(){
    this->RT60 = 0.36f;
    this->roomWidth = 8.0f; //x
    this->roomHeight = 5.0f; //y
    this->roomCeiling = 3.1f; //z
    this->soundSourceLoc = Vector3D(7.0f, 4.f, 1.2f);
    this->listenerLoc = Vector3D(5.f, 2.55,1.2f);
    
    
    Room = Cuboid(roomWidth, roomHeight, roomCeiling);

    ssLoc = soundSourceLoc;
    lLoc = listenerLoc;

    lLocLE = Vector3D(lLoc.x - RADIUSOFHEAD, lLoc.y, lLoc.z);
    lLocRE = Vector3D(lLoc.x + RADIUSOFHEAD, lLoc.y, lLoc.z);
    
    
#ifdef FDN_ER_SECOND_ORDER
    size_t totalARE = AREInit(AREFIRSTORDER);
    
    
#elif defined (USE_ISM)
    totalISM = ISMInit(ISMORDER);
    
    
#endif
 
    
}


/***********************************
 *    PRODUCE IMPULSE RESPONSE FN   *
 ***********************************/

void FDN::impulseResponse(long numSamples, std::ofstream* outputFileL, std::ofstream* outputFileR){
    
    
    float fdnOutL, fdnOutR, input;
    
    float erHighestOrder;
    float directHRTFOut[2];
    float directRayOut[2];
    float zero = 0.0f;
    float one = 1.0f;

    
    for (int i = 0; i < numSamples; i++){
        erHighestOrder = 0.0f;
        
        fdnOutL = 0.0f;
        fdnOutR = 0.0f;
        
        // set the input to one for the first sample, zero otherwise.
        input = i==0 ? one : zero;

        // direct rays
        processDirectRays(input, directRayOut);
        attenuateDirectRay(directRayOut);
        HRTFFilterDirect(directRayOut, directHRTFOut);
        
#ifdef  USE_ISM

        //multiply by division constant
        if (i<ISM_IRLength){
            erHighestOrder = leftEROutLastOrder[i] * LRDivisionConstant + rightEROutLastOrder[i] * LRDivisionConstant;
        }
        else{
            erHighestOrder = 0.0f;
        }
        
        //FDN Processing
        processAudio(&erHighestOrder, &fdnOutL, &fdnOutR);
        
        if (i<ISM_IRLength){
            *outputFileL << fdnOutL + leftEROut[i] + directHRTFOut[0]<< ",";
            *outputFileR << fdnOutR + rightEROut[i] + directHRTFOut[1]<< ",";
        }
        else{
            *outputFileL << fdnOutL  + directHRTFOut[0] << ",";
            *outputFileR << fdnOutR + directHRTFOut[1] << ",";
        }

#else
        //multiply by division constant
        if (i<ARE_IRLength){
            erHighestOrder = leftEROutLastOrder[i] * LRDivisionConstant + rightEROutLastOrder[i] * LRDivisionConstant;
        }
        else{
            erHighestOrder = 0.0f;
        }
        
        //FDN Processing
        processAudio(&erHighestOrder, &fdnOutL, &fdnOutR);
        
        if (i<ARE_IRLength){
            *outputFileL << fdnOutL + leftEROut[i] + directHRTFOut[0]<< ",";
            *outputFileR << fdnOutR + rightEROut[i] + directHRTFOut[1]<< ",";
        }
        else{
            *outputFileL << fdnOutL  + directHRTFOut[0] << ",";
            *outputFileR << fdnOutR + directHRTFOut[1] << ",";
        }
#endif
        
    }
    
    printf("GAIN OFF \n");

    
    *outputFileL << "\n";
    *outputFileR << "\n";
    

}





/***********************************
 *    CONSTRUCTOR, SETUP VARS   *
 ***********************************/

FDN::FDN(int rvType)
{

    
    clock_t begin = clock();
    
    setGains();
    
    float delay = ssLoc.distance(lLoc) / 340.f;
    directRays.setTimeSafe(delay);
    

    
    /***********************************
     *    setup the FDN               *
     ***********************************/
    
    Room.sliceCube(NUMFDNDELAYS);//FOR SECOND ORDER, change to 9 if want more taps
    GainValues = Gains(DMIN, Room.elements, Room.area, Room.volume, RT60, true);
    GainValues.calculateUpsilon(Room.segmentedSides, lLoc, NUMFDNDELAYS, Room.area);
    Room.getDelayValues(delayTimes, lLocLE, lLocRE, ssLoc, FDN_SAMPLE_RATE, 0);
    
//    for (int i =0; i<NUMFDNDELAYS; i++){
//        printf("Upsilon %f \n", GainValues.upsilon[i]);
//    }

    float meanFreePath = 4*Room.volume / Room.area;
//
//    //Normal Method
//    //Hardcode the length and channel of FDN based on room dimension according to Wendt et al
//    std::default_random_engine generator;
//    std::uniform_real_distribution<float> distribution(-0.1,0.1);
//    
//    std::uniform_real_distribution<float> distributionX(-roomWidth/2.f,roomWidth/2.f);
//    std::uniform_real_distribution<float> distributionY(-roomHeight/2.f,roomHeight/2.f);
//    
////    printf("%f %f %f \n", distributionTwo(generator), distributionTwo(generator),distributionTwo(generator));
//    float averageLength = (Room.xLength + Room.yLength + Room.zLength)/3.0f;
//    
//    delayTimes[0] = Room.xLength/340.f*44100.f + distribution(generator) * averageLength /340.f*44100.f;
////    FDNChannel[0] = determineChannel(Room.floor.getMidpoint().x, Room.floor.getMidpoint().y);
////    randomPointsOnRectangle(Room.side1.corner, Room.side1.S1, Room.side1.S2, &tempPoints[0], 1);
//    tempPoints[0] = Vector3D(Room.side1.getMidpoint().x + distributionX(generator), Room.side1.getMidpoint().y, Room.side1.getMidpoint().z);
//    
//    delayTimes[1] = Room.xLength/340.f*44100.f + distribution(generator) * averageLength /340.f*44100.f;
////    FDNChannel[1] = determineChannel(Room.floor.getMidpoint().x, Room.floor.getMidpoint().y);
//    tempPoints[1] = Vector3D(Room.side3.getMidpoint().x + distributionX(generator), Room.side3.getMidpoint().y, Room.side3.getMidpoint().z);
////    randomPointsOnRectangle(Room.side3.corner, Room.side3.S1, Room.side3.S2, &tempPoints[1], 1);
//    
//    delayTimes[2] = Room.xLength/340.f*44100.f + distribution(generator) * averageLength /340.f*44100.f;
////    FDNChannel[2] = determineChannel(Room.ceiling.getMidpoint().x, Room.ceiling.getMidpoint().y);
//    tempPoints[2] = Vector3D(Room.side3.getMidpoint().x + distributionX(generator), Room.side3.getMidpoint().y, Room.side3.getMidpoint().z);
////    randomPointsOnRectangle(Room.side3.corner, Room.side3.S1, Room.side3.S2, &tempPoints[2], 1);
//    
//    delayTimes[3] = Room.xLength/340.f*44100.f + distribution(generator) * averageLength /340.f*44100.f;
////    FDNChannel[3] = determineChannel(Room.ceiling.getMidpoint().x, Room.ceiling.getMidpoint().y);
//    tempPoints[3] = Vector3D(Room.side1.getMidpoint().x + distributionX(generator), Room.side1.getMidpoint().y, Room.side1.getMidpoint().z);
////    randomPointsOnRectangle(Room.side1.corner, Room.side1.S1, Room.side1.S2, &tempPoints[3], 1);
//    
//    delayTimes[4] = Room.yLength/340.f*44100.f + distribution(generator) * averageLength /340.f*44100.f;
////    FDNChannel[4] = determineChannel(Room.side2.getMidpoint().x, Room.side2.getMidpoint().y);
//    tempPoints[4] = Vector3D(Room.side2.getMidpoint().x , Room.side2.getMidpoint().y + distributionY(generator), Room.side2.getMidpoint().z);
////    randomPointsOnRectangle(Room.side2.corner, Room.side2.S1, Room.side2.S2, &tempPoints[4], 1);
//    
//    delayTimes[5] = Room.yLength/340.f*44100.f + distribution(generator) * averageLength /340.f*44100.f;
////    FDNChannel[5] = determineChannel(Room.side2.getMidpoint().x, Room.side2.getMidpoint().y);
//    tempPoints[5] = Vector3D(Room.side2.getMidpoint().x , Room.side2.getMidpoint().y + distributionY(generator), Room.side2.getMidpoint().z);
//    //    randomPointsOnRectangle(Room.side2.corner, Room.side2.S1, Room.side2.S2, &tempPoints[5], 1);
//    
////    FDNChannel[5] = 2;
//    delayTimes[6] = Room.yLength/340.f*44100.f + distribution(generator) * averageLength /340.f*44100.f;
////    FDNChannel[6] = determineChannel(Room.side4.getMidpoint().x, Room.side4.getMidpoint().y);
//    tempPoints[6] = Vector3D(Room.side4.getMidpoint().x , Room.side4.getMidpoint().y + distributionY(generator), Room.side4.getMidpoint().z);
//    //     randomPointsOnRectangle(Room.side4.corner, Room.side4.S1, Room.side4.S2, &tempPoints[6], 1);
//    
////    FDNChannel[6] = 5;
//    delayTimes[7] = Room.yLength/340.f*44100.f + distribution(generator) * averageLength /340.f*44100.f;
////    FDNChannel[7] = determineChannel(Room.side4.getMidpoint().x, Room.side4.getMidpoint().y);
//    tempPoints[7] = Vector3D(Room.side4.getMidpoint().x , Room.side4.getMidpoint().y + distributionY(generator), Room.side4.getMidpoint().z);
////    randomPointsOnRectangle(Room.side4.corner, Room.side4.S1, Room.side4.S2, &tempPoints[7], 1);
//    
////    FDNChannel[7] = 6;
//    delayTimes[8] = Room.zLength/340.f*44100.f + distribution(generator) * averageLength /340.f*44100.f;
////    FDNChannel[8] = determineChannel(Room.side1.getMidpoint().x, Room.side1.getMidpoint().y);
//    tempPoints[8] = Vector3D(Room.side1.getMidpoint().x + distributionX(generator) , Room.side1.getMidpoint().y , Room.side1.getMidpoint().z);
////    randomPointsOnRectangle(Room.side1.corner, Room.side1.S1, Room.side1.S2, &tempPoints[8], 1);
//    
////    FDNChannel[8] = 3;
//    delayTimes[9] = Room.zLength/340.f*44100.f + distribution(generator) * averageLength /340.f*44100.f;
////    FDNChannel[9] = determineChannel(Room.side1.getMidpoint().x, Room.side1.getMidpoint().y);
//    tempPoints[9] = Vector3D(Room.side1.getMidpoint().x + distributionX(generator) , Room.side1.getMidpoint().y , Room.side1.getMidpoint().z);
////    randomPointsOnRectangle(Room.side1.corner, Room.side1.S1, Room.side1.S2, &tempPoints[9], 1);
//    
//    
////    FDNChannel[9] = 4;
//    delayTimes[10] = Room.zLength/340.f*44100.f + distribution(generator) * averageLength /340.f*44100.f;
////    FDNChannel[10] = determineChannel(Room.side3.getMidpoint().x, Room.side3.getMidpoint().y);
//    tempPoints[10] = Vector3D(Room.side3.getMidpoint().x + distributionX(generator), Room.side3.getMidpoint().y, Room.side3.getMidpoint().z);
////    randomPointsOnRectangle(Room.side1.corner, Room.side1.S1, Room.side1.S2, &tempPoints[10], 1);
//    
//    
////    FDNChannel[10] = 0;
//    delayTimes[11] = Room.zLength/340.f*44100.f + distribution(generator) * averageLength /340.f*44100.f;
////    FDNChannel[11] = determineChannel(Room.side3.getMidpoint().x, Room.side3.getMidpoint().y);
//    tempPoints[11] = Vector3D(Room.side3.getMidpoint().x + distributionX(generator), Room.side3.getMidpoint().y, Room.side3.getMidpoint().z);
////    randomPointsOnRectangle(Room.side3.corner, Room.side3.S1, Room.side3.S2, &tempPoints[11], 1);
//    
//    
////    FDNChannel[11] = 7;
//    delayTimes[12] = Room.yLength/340.f*44100.f + distribution(generator) * averageLength /340.f*44100.f;
////    FDNChannel[12] = determineChannel(Room.side3.getMidpoint().x, Room.ceiling.getMidpoint().y);
//    tempPoints[12] = Vector3D(Room.side2.getMidpoint().x , Room.side2.getMidpoint().y + distributionY(generator), Room.side2.getMidpoint().z);
////    randomPointsOnRectangle(Room.side2.corner, Room.side2.S1, Room.side2.S2, &tempPoints[11], 1);
//    
//    delayTimes[13] = Room.yLength/340.f*44100.f + distribution(generator) * averageLength /340.f*44100.f;
////    FDNChannel[13] = determineChannel(Room.side4.getMidpoint().x, Room.side4.getMidpoint().y);
//    tempPoints[13] = Vector3D(Room.side4.getMidpoint().x , Room.side4.getMidpoint().y + distributionY(generator), Room.side4.getMidpoint().z);
////    randomPointsOnRectangle(Room.side4.corner, Room.side4.S1, Room.side4.S2, &tempPoints[13], 1);
//    
//    
//    delayTimes[14] = Room.zLength/340.f*44100.f + distribution(generator) * averageLength /340.f*44100.f;
////    FDNChannel[14] = determineChannel(Room.side2.getMidpoint().x, Room.side3.getMidpoint().y);
//    tempPoints[14] = Vector3D(Room.side2.getMidpoint().x , Room.side2.getMidpoint().y + distributionY(generator), Room.side2.getMidpoint().z);
////    randomPointsOnRectangle(Room.side2.corner, Room.side2.S1, Room.side2.S2, &tempPoints[14], 1);
//    
//    delayTimes[15] = Room.zLength/340.f*44100.f + distribution(generator) * averageLength /340.f*44100.f;
////    FDNChannel[15] = determineChannel(Room.side1.getMidpoint().x, Room.floor.getMidpoint().y);
//    tempPoints[15] = Vector3D(Room.side4.getMidpoint().x , Room.side4.getMidpoint().y + distributionY(generator), Room.side4.getMidpoint().z);
////    randomPointsOnRectangle(Room.side4.corner, Room.side4.S1, Room.side4.S2, &tempPoints[15], 1);
//    

    // count total delay time
    totalDelayTime = 0;
    for (size_t i = 0; i < NUMFDNDELAYS; i++){
        totalDelayTime += delayTimes[i];
        std::cout << delayTimes[i] << ", ";
    }
    
    //assert its larger than MFP
    assert(totalDelayTime/NUMFDNDELAYS > meanFreePath/340.f * 44100);
    printf("\n Mean free path in Hz %f, average delay length in Hz %f, ratio %f \n", meanFreePath/340.f * 44100, (float)totalDelayTime/(float)NUMFDNDELAYS,(float)totalDelayTime/(float)NUMFDNDELAYS/(meanFreePath/340.f * 44100) );

    
    //Add reverb ITD
    setTempPoints(); //commented for normal method, uncommented for channel method
    calculateAdditionalDelays();
    setupHRTFFunctions();
    
    // randomly shuffle the order of delay times in the array
    randomPermutation1Channel(delayTimes, NUMFDNDELAYS, 1);
    
    // declare and initialise the delay buffers
    delayBuffers = new float[totalDelayTime];
    for (int i = 0; i < totalDelayTime; i++) delayBuffers[i] = 0;
    
    // calculate start and end indices for each delay line
    startIndices[0] = 0;
    endIndices[0] = delayTimes[0];
    for (int i = 1; i < NUMFDNDELAYS; i++){
        startIndices[i] = startIndices[i - 1] + delayTimes[i - 1];
        endIndices[i] = startIndices[i] + delayTimes[i];
    }
    
    // since we copy one input channel into many delay lines, we need to attenuate the input to avoid clipping.
    fdnMultiplexInputAttenuation = 1.0f / sqrt((float)NUMFDNDELAYS);
    
    // we specify our feedback matrix with elements of 1, and -1.  In order to keep it unitary,
    // we need to scale it down.
    // here we calculate the dividend needed to scale the feedback matrix in order to keep it unitary.

    if (rvType > 0) // fast hadamard transform
        matrixAttenuation = 1.0f / sqrt(NUMFDNDELAYS);
    else  // 4x4 circular network
        matrixAttenuation = 1.0f / sqrt(4.0f);

    // randomise output tap signs
    for (int i = 0; i < NUMFDNDELAYS; i++){
        tapSigns[i] = pow(-1.0, i * 31 % 71);
        //printf("tapsigns %d %f \n", i, tapSigns[i]);
    }
    
    hadamardTransform(tapSigns, temp3, NUMFDNDELAYS);
    memcpy(tapSigns, temp3, sizeof(float) * NUMFDNDELAYS);
    
    resetDelay(); //setting jot's gain
    
    clock_t end = clock();
    double elapsed_msecs = double(end - begin) / CLOCKS_PER_SEC * 1000.f;
    printf("\n Setup Time is : %f ms \n", elapsed_msecs);

}


/***********************************
 *   PROCESS FDN FUNCTION        *
 ***********************************/

void FDN::processAudio(float* pInput, float* pOutputL, float* pOutputR)
{
    
    // Read the Input
    *pOutputL = 0;
    *pOutputR = 0;
    
    float FDNTankOut[CHANNELS] = {};
    float FDNTankOutL[CHANNELS] = {};
    float FDNTankOutR[CHANNELS] = {};
    
    int i;
    float xn = *pInput * fdnMultiplexInputAttenuation;
    
    
    // mono output copied to both channels
    for (i = 0; i < NUMFDNDELAYS; i++) {
        // copy delay line signals to output buffer and attenuate
        outputs[i] = delayBuffers[rwIndices[i]] * tapGains[i];
    }
    
    vDSP_vmul(outputs, 1, GainValues.upsilon, 1, temp4, 1, NUMFDNDELAYS);
    processTankOutFDN(FDNTankOut, temp4);
    processHRTFFDN(FDNTankOut, FDNTankOutL, FDNTankOutR);
    
    //ITD
    addReverbDelay(FDNTankOutL, FDNTankOutR); //Temp point delays
    
    vDSP_sve(FDNTankOutL, 1, pOutputL, CHANNELS);
    vDSP_sve(FDNTankOutR, 1, pOutputR, CHANNELS);
    
    // Filling up temp3 with copies of the input sample
    vDSP_vfill(&xn, temp3, 1, NUMFDNDELAYS);
    // Randomizing the sign of the input sample vector
    vDSP_vmul(temp3, 1, tapSigns, 1, temp3, 1, NUMFDNDELAYS);
    // Mixing the input sample vector with the feedback vector
    vDSP_vadd(temp3, 1, outputs, 1, outputs, 1, NUMFDNDELAYS);
    
    // Matrix mixing operation
    hadamardTransform(outputs, temp3, NUMFDNDELAYS);
    
    // Write into the delay line inputs
    for (int i = 0; i < NUMFDNDELAYS; i++) delayBuffers[rwIndices[i]] = temp3[i];
    
    
    incrementIndices();
    
    
}


/***********************************
 *    FDN HELPER FUNCTIONS       *
 ***********************************/
inline void FDN::hadamard16(float* input, float* output){
    // level 1
    // +
    temp16[0] = input[0] + input[8];
    temp16[1] = input[1] + input[9];
    temp16[2] = input[2] + input[10];
    temp16[3] = input[3] + input[11];
    temp16[4] = input[4] + input[12];
    temp16[5] = input[5] + input[13];
    temp16[6] = input[6] + input[14];
    temp16[7] = input[7] + input[15];
    // -
    temp16[8] = input[0] - input[8];
    temp16[9] = input[1] - input[9];
    temp16[10] = input[2] - input[10];
    temp16[11] = input[3] - input[11];
    temp16[12] = input[4] - input[12];
    temp16[13] = input[5] - input[13];
    temp16[14] = input[6] - input[14];
    temp16[15] = input[7] - input[15];
    
    // level 2.1
    //+
    output[0] = temp16[0] + temp16[4];
    output[1] = temp16[1] + temp16[5];
    output[2] = temp16[2] + temp16[6];
    output[3] = temp16[3] + temp16[7];
    //-
    output[4] = temp16[0] - temp16[4];
    output[5] = temp16[1] - temp16[5];
    output[6] = temp16[2] - temp16[6];
    output[7] = temp16[3] - temp16[7];
    
    //+
    output[8] = temp16[8] + temp16[12];
    output[9] = temp16[9] + temp16[13];
    output[10] = temp16[10] + temp16[14];
    output[11] = temp16[11] + temp16[15];
    //-
    output[12] = temp16[8] - temp16[12];
    output[13] = temp16[9] - temp16[13];
    output[14] = temp16[10] - temp16[14];
    output[15] = temp16[11] - temp16[15];
    
    // level 3
    // +
    temp16[0] = output[0] + output[2];
    temp16[1] = output[1] + output[3];
    // -
    temp16[2] = output[0] - output[2];
    temp16[3] = output[1] - output[3];
    
    // +
    temp16[4] = output[4] + output[6];
    temp16[5] = output[5] + output[7];
    // -
    temp16[6] = output[4] - output[6];
    temp16[7] = output[5] - output[7];
    
    // +
    temp16[8] = output[8] + output[10];
    temp16[9] = output[9] + output[11];
    // -
    temp16[10] = output[8] - output[10];
    temp16[11] = output[9] - output[11];
    
    // +
    temp16[12] = output[12] + output[14];
    temp16[13] = output[13] + output[15];
    // -
    temp16[14] = output[12] - output[14];
    temp16[15] = output[13] - output[15];
    
    // level 4
    output[0] = temp16[0] + temp16[1];
    output[1] = temp16[0] - temp16[1];
    output[2] = temp16[2] + temp16[3];
    output[3] = temp16[2] - temp16[3];
    
    output[4] = temp16[4] + temp16[5];
    output[5] = temp16[4] - temp16[5];
    output[6] = temp16[6] + temp16[7];
    output[7] = temp16[6] - temp16[7];
    
    output[8] = temp16[8] + temp16[9];
    output[9] = temp16[8] - temp16[9];
    output[10] = temp16[10] + temp16[11];
    output[11] = temp16[10] - temp16[11];
    
    output[12] = temp16[12] + temp16[13];
    output[13] = temp16[12] - temp16[13];
    output[14] = temp16[14] + temp16[15];
    output[15] = temp16[14] - temp16[15];
}


// input and output must not be the same array
inline void FDN::hadamardTransform(float* input, float* output, size_t length){
    size_t partitionSize = length;
    float* tmpIn = input;
    float* tmpOut = temp1;
    bool tmpState = true;
    
    // iteratively calculated recursion
    while (partitionSize > 16) {
        size_t halfPartitionSize = partitionSize >> 1;
        for (size_t i = 0; i < length; i += partitionSize){
            // copy all the lower terms into place
            memcpy(tmpOut+i,tmpIn+i,sizeof(float)*halfPartitionSize);
            memcpy(tmpOut+i+halfPartitionSize,tmpIn+i,sizeof(float)*halfPartitionSize);
            // sum all the higher terms into place
            for (size_t j=i; j<halfPartitionSize+i; j++) {
                size_t idx2 = j+halfPartitionSize;
                tmpOut[j] += tmpIn[idx2];
                tmpOut[idx2] -= tmpIn[idx2];
            }
        }
        
        // swap temp buffers to avoid using the same memory for reading and writing.
        tmpIn = tmpOut;
        tmpState = !tmpState;
        if (tmpState) tmpOut = temp1;
        else tmpOut = temp2;
        partitionSize >>= 1;
    }
    
    // base case
    for (int i = 0; i < length; i += 16) hadamard16(tmpIn + i, output + i);
    
    // scale the output to make it unitary
    vDSP_vsmul(output,1,&matrixAttenuation,output,1,(size_t) length);
}


void FDN::resetReadIndices(){
    for (size_t i = 0; i < NUMFDNDELAYS; i++){
        rwIndices[i] = startIndices[i];
    }
}








// computes the appropriate feedback gain attenuation
// to get a decay envelope with the specified RT60 time (in seconds)
// from a delay line of the specified length.
//
// This formula comes from solving EQ 11.33 in DESIGNING AUDIO EFFECT PLUG-INS IN C++ by Will Pirkle
// which is attributed to Jot, originally.
double gain(double rt60, double delayLengthInSamples, double sampleRate) {
    //    printf("Delay length in samples: %f\n", delayLengthInSamples);
    return pow(10.f, (-3.0 * delayLengthInSamples) / (rt60 * sampleRate));
}




// larger values decay more slowly
void FDN::setDecayTime(double rt60){
    if (rt60 < 0) {
        for (int i = 0; i < NUMFDNDELAYS; i++){
            tapGains[i] = 1.0f;
        }
    }
    else if (rt60 == 0) {
        for (int i = 0; i < NUMFDNDELAYS; i++){
            tapGains[i] = 0.0f;
        }
    }
    else {
        for (int i = 0; i < NUMFDNDELAYS; i++) {
            tapGains[i] = gain(rt60, delayTimes[i],FDN_SAMPLE_RATE);
//            printf("TapGains: %f delayTime %zu \n", tapGains[i], delayTimes[i]);
        }
    }
}


void FDN::resetDelay()
{
    // clear the buffers
    memset(delayBuffers, 0, sizeof(float)* totalDelayTime);
    memset(outputs, 0, sizeof(float)* NUMDELAYSEXT);
    
    // init read indices
    resetReadIndices();
    
    setDecayTime(RT60);
//    std::srand(0);
}


inline void FDN::incrementIndices(){
    for (int i = 0; i < NUMFDNDELAYS; i++){
        rwIndices[i]++;
        if (rwIndices[i] >= endIndices[i])
            rwIndices[i] = startIndices[i];
    }
}


void swap(size_t* a, size_t b, size_t c){
    size_t tmp = a[b];
    a[b] = a[c];
    a[c] = tmp;
}

// randomly permute the elements of a without swapping any delay times between channels
void FDN::randomPermutation(size_t* a, size_t length, size_t audioChannels){
    for (size_t j = 0; j < audioChannels; j++){
        for (size_t i = j; i < length; i += audioChannels){
            swap(a, i, j + (std::rand() % (length / audioChannels))*audioChannels);
        }
    }
}

void FDN::randomPermutation1Channel(size_t* a, size_t length, size_t audioChannels){
    int j = 0;
    for (int i = j; i < length; i += audioChannels){
        updateRand();
        swap(a, i, j + (rands % (length / audioChannels))*audioChannels);
    }
}

float max(float* window, int length){
    float max = 0.0;
    for (int i = 0; i < length; i++){
        if (fabsf(window[i]) > max) max = fabsf(window[i]);
    }
    return max;
}

int thresholdCount(float* window, int length, float threshold){
    int count = 0;
    for (int i = 0; i < length; i++){
        if (fabsf(window[i]) > threshold) count++;
    }
    return count;
}

//
//// see: http://software.intel.com/en-us/articles/fast-random-number-generator-on-the-intel-pentiumr-4-processor/
//// rather than returning an output, this function updates a class variable so that we only have to generate 1 random number for each sample.
inline void FDN::updateRand()
{
    randomSeed = (214013 * randomSeed + 2531011);
    rands = (randomSeed >> 16) & 0x7FFF;
}


/***********************************
 *    PROCESS AUDIO FUNCTIONS       *
 ***********************************/
//ITD for direct rays
inline void FDN::processDirectRays(float input, float* directRaysOutput){
    directRaysOutput[0] = directRaysStereo[0].process(input);
    directRaysOutput[1] = directRaysStereo[1].process(input);
}

//Multiplexer FDN
void FDN::processTankOutFDN(float FDNTankOut[CHANNELS], float* inputFDN){
    for (size_t i = 0; i < NUMFDNDELAYS; i++){
        size_t channel = delayTimesChannel[i];
        FDNTankOut[channel] += inputFDN[i];
    }
}

//HRTF ER
inline    void FDN::processHRTFER(float ERTankOutL[CHANNELS], float ERTankOutR[CHANNELS], float* inputER){
    for (size_t i = 0; i < CHANNELS; i++){
        ERTankOutL[i] = leftEarFilterER[i].process(inputER[i]);
        ERTankOutR[i] = rightEarFilterER[i].process(inputER[i]);
    }
}

//HRTF FDN
inline void FDN::processHRTFFDN(float* input,  float* fdnTankOutLeft, float* fdnTankOutRight){
    //Filter left & right ear
    for (size_t i = 0; i < CHANNELS; i++){
        fdnTankOutLeft[i] = leftEarFilter[i].process(input[i]);
        fdnTankOutRight[i] = rightEarFilter[i].process(input[i]);
    }
}

//HRTF Direct
inline void FDN::HRTFFilterDirect(float directRay[2], float directRayAfterHRTF[2]){
    //Filter direct rays
    directRayAfterHRTF[0] = directRayHRTFFilter[0].process(directRay[0]);
    directRayAfterHRTF[1] = directRayHRTFFilter[1].process(directRay[1]);
    
}

//Direct Ray Attenuation
void FDN::attenuateDirectRay(float* input){
    input[0] *= directRayAttenuation[0];
    input[1] *= directRayAttenuation[1];
}


/***********************************
 *    setup the HRTF FUNCTIONS    *
 ***********************************/

void FDN::setupHRTFFunctions(){
    
    setFDNDelayOutputChannels();

    //Direct Rays
    setDirectDelayTimes();
    setDirectGains();
    setDirectRayAngles();
    setDirectSingleTapDelay();
    
    //HRTF filters
    setHRTFFilters();
}





/***********************************
 *    DIRECT RAYS SETUP FUNCTIONS    *
 ***********************************/

void FDN::setDirectSingleTapDelay(){
    directRaysStereo[0].setTimeSafe(directDelayTimes[0]);
    directRaysStereo[1].setTimeSafe(directDelayTimes[1]);
}


inline void FDN::setDirectDelayTimes()
{
    
    //Calculate delay from source to receiver
    float directDelayLeft = ssLoc.distance(lLocLE) / SOUNDSPEED;
    float directDelayRight = ssLoc.distance(lLocRE) / SOUNDSPEED;
    printf("direct delay left %f , direct delay right %f \n", directDelayLeft, directDelayRight);
    directDelayTimes[0] = directDelayLeft;
    directDelayTimes[1] = directDelayRight;
}

void FDN::setDirectGains(){
    directRayAttenuation[0] = expf(-ALPHA * ssLoc.distance(lLocLE)) * DMIN/ssLoc.distance(lLocLE);
    directRayAttenuation[1] = expf(-ALPHA * ssLoc.distance(lLocRE)) * DMIN/ssLoc.distance(lLocRE);
    printf("Direct left %f direct right %f \n", directRayAttenuation[0], directRayAttenuation[1]);
}

//Setting direct ray angle with respect to listener location
void FDN::setDirectRayAngles()
{
    float yDiff, xDiff;
    yDiff = ssLoc.y - lLocLE.y;
    xDiff = ssLoc.x - lLocLE.x;
    float angle = atan2(xDiff, yDiff) * 180.0f / M_PI;
    
    //switch to 360.f
    
    //Make it between 0 to 360
    if (fabs(angle - 0.0f)<0.0001){
        angle = 0.0f;
    }
    
    if (angle < 0.0f ){
        angle = 360.f + angle;
    }
    
    directRayAngles[0] = angle;
    
    float yDiff2, xDiff2;
    yDiff2 = ssLoc.y - lLocRE.y;
    xDiff2 = ssLoc.x - lLocRE.x;
    float angle2 = atan2(xDiff2, yDiff2) * 180.0f / M_PI;
    
    //switch to 360.f
    
    //Make it between 0 to 360
    if (fabs(angle2 - 0.0f)<0.0001){
        angle2 = 0.0f;
    }
    
    if (angle2 < 0.0f ){
        angle2 = 360.f + angle2;
    }
    
    directRayAngles[1] = angle2;
}



/***********************************
 *    HRTF SETUP FUNCTIONS    *
 ***********************************/


//Returns angle between 0 to 360
float FDN::channeltoangleNormal(size_t channel){
    float channelWidth = 360.0f / (float)NUMFDNDELAYS;
    //    size_t midChannel = CHANNELS / 2 - 1;
    return (float)channel * channelWidth + (channelWidth / 2.0f);
    //    if (channel <= midChannel){f
    //        return (float)channel * channelWidth + (channelWidth / 2.0f);
    //    }
    //    else{
    //        return (-1.f  * (float) (CHANNELS - channel) * channelWidth )+ (channelWidth / 2.0f);
    //    }
}


size_t FDN::angleToChannel(float angleInDegrees){
    
    //Make it between 0 to 360
    if (fabs(angleInDegrees - 0.0f)<0.0001){
        angleInDegrees = 0.0f;
    }
    
    if (angleInDegrees < 0.0f ){
        angleInDegrees = 360.f + angleInDegrees;
    }
    float channelWidth = 360.0f / CHANNELS;
    
    //Channel 0 to 7 from 12 o'clock, clockwise according from DFX
    size_t channel =(size_t)floorf(angleInDegrees / channelWidth);
    //    printf("From angle: %f to channel %zu \n", angleInDegrees, channel);
    return channel;
}


float FDN::determineAngle(float x,float y){
    float xL = lLoc.x;
    float yL = lLoc.y;
    
    float xDistance = x - xL;
    float yDistance = y - yL;
    
    //Special cases
    if (fabs(yDistance - 0.00001) < 0.0f){
        // zero y, can be at 90 or -90
        if (xDistance >= 0){
            return angleToChannel(90.f);
        }
        else{
            return angleToChannel(-90.f);
        }
    }
    
    if (fabs(xDistance - 0.00001) < 0.0f){
        //zero x, can be at 0 or 180;
        if (yDistance >= 0){
            return angleToChannel(0.f);
        }
        else{
            return angleToChannel(180.f);
        }
    }
    
    //Normal case
    //x over y because we want azimuth
    float angleInDegrees = atan2f(xDistance, yDistance) * 180.0f / M_PI;
    
    //Make it between 0 to 360
    if (fabs(angleInDegrees - 0.0f)<0.0001){
        angleInDegrees = 0.0f;
    }
    
    if (angleInDegrees < 0.0f ){
        angleInDegrees = 360.f + angleInDegrees;
    }
    
    //printf("Angle %f \n", angle2);
    return angleInDegrees;
    
}

void FDN::setHRTFFilters(){
    
    for (int i =0; i<CHANNELS; i++){
        leftEarFilter[i].setAngle(convertAzimuthToLeftEar(determineAngle(tempPoints[i].x, tempPoints[i].y)), FDN_SAMPLE_RATE, false);
        rightEarFilter[i].setAngle(convertAzimuthToRightEar(determineAngle(tempPoints[i].x, tempPoints[i].y)), FDN_SAMPLE_RATE, false);
    }
    
    directRayHRTFFilter[0].setAngle(convertAzimuthToLeftEar(directRayAngles[0]), FDN_SAMPLE_RATE,false);
    directRayHRTFFilter[1].setAngle(convertAzimuthToRightEar(directRayAngles[1]), FDN_SAMPLE_RATE, true);
    
}

//setting the channels for each FDN output
void FDN::setFDNDelayOutputChannels(){
    for(size_t i = 0; i<NUMFDNDELAYS; i++){
        printf("room x %f room y %f Lx %f Ly %f \n", Room.segmentedSides[i].getMidpoint().x, Room.segmentedSides[i].getMidpoint().y, lLoc.x, lLoc.y);
        delayTimesChannel[i] = determineChannel(Room.segmentedSides[i].getMidpoint().x, Room.segmentedSides[i].getMidpoint().y,0);
        printf("i %zu Channels %zu \n", i, delayTimesChannel[i]);
    }
}




float FDN::convertAzimuthToLeftEar(float azimuth){

    //change azimuth to be between 0 to 180 and 0 to -180
    if (azimuth >= 360.f){
        azimuth -= 360.f;
    }
    
    if (azimuth > 180.f){
        azimuth = azimuth - 360.f;
    }
    
    if(azimuth > 90.f and azimuth <= 180.f){
        float leftEarAngle = azimuth + 90.f - 360.f;
        return leftEarAngle;
    }
    
    return azimuth + 90.f;
    
    
}

float FDN::convertAzimuthToRightEar(float azimuth){
    
    //    //azimuth between 0 to 360
    //    float azimuth = channeltoangleNormal(channel);
    //
    //change azimuth to be between 0 to 180 and 0 to -180
    if (azimuth >= 360.f){
        azimuth -= 360.f;
    }
    
    
    if (azimuth > 180.f){
        azimuth = azimuth - 360.f;
    }
    
    if(azimuth >= -180.f and azimuth < -90.f){
        float rightEarAngle = azimuth - 90.f + 360.f;
        return rightEarAngle;
    }
    
    return azimuth - 90.f;
    
}




/***********************************
 *    PROCESS ARE FUNCTION          *
 ***********************************/
size_t FDN::AREInit(size_t numberOfFirstOrder){
    
    
    // find the longest reflection time
    float longestReflectionTime = 0.0f;
    // total number of ARE sources
    size_t numTotalSources = numberOfFirstOrder*numberOfFirstOrder;
    printf("Num totalSources : %zu \n", numTotalSources);
    
    float* leftEarGains = static_cast<float*>(malloc(sizeof(float)*numTotalSources));
    float* rightEarGains = static_cast<float*>(malloc(sizeof(float)*numTotalSources));
    float* leftEarDelayTimes = static_cast<float*>(malloc(sizeof(float)*numTotalSources));
    float* rightEarDelayTimes = static_cast<float*>(malloc(sizeof(float)*numTotalSources));
    
    Vector3D* ARESources = static_cast<Vector3D*>(malloc(sizeof(Vector3D)*numTotalSources));
    
    RoomARE = Cuboid(roomWidth, roomHeight, roomCeiling);
    RoomARE.sliceCube((int)numberOfFirstOrder);

    GainValuesARE = Gains(DMIN, RoomARE.elements, RoomARE.area, RoomARE.volume, RT60, true);
    
    //this comes up with Beta
    GainValuesARE.calculateAllGainsSecondOrder(temp3, RoomARE.segmentedSides, lLocLE, ssLoc, RT60, DMIN, numberOfFirstOrder, leftEarGains, ARESources);
    GainValuesARE.calculateAllGainsSecondOrder(temp3, RoomARE.segmentedSides, lLocRE, ssLoc, RT60, DMIN, numberOfFirstOrder, rightEarGains, ARESources);
    
    //get the second order Delay times
    RoomARE.getSecondOrderDelayValues(leftEarDelayTimes, lLocLE, ssLoc, FDN_SAMPLE_RATE);
    RoomARE.getSecondOrderDelayValues(rightEarDelayTimes, lLocRE, ssLoc, FDN_SAMPLE_RATE);

    
    for(size_t i=0; i<numTotalSources; i++){
//        printf("reflection time L %f reflection time R %f \n", leftEarDelayTimes[i], rightEarDelayTimes[i]);
        if(leftEarDelayTimes[i] > longestReflectionTime)
            longestReflectionTime = leftEarDelayTimes[i];
        if(rightEarDelayTimes[i] > longestReflectionTime)
            longestReflectionTime = rightEarDelayTimes[i];
    }
    
    
    
    // allocate memory for producing the impulse response
    float sampleRate = 44100;
    ARE_IRLength = 1+(longestReflectionTime*sampleRate)+HRTF_IR_LENGTH;
    printf("Longest reflection time :%zu\n", ARE_IRLength);
    
    
    leftEROut  = static_cast<float*>(malloc(sizeof(float)*ARE_IRLength));
    rightEROut = static_cast<float*>(malloc(sizeof(float)*ARE_IRLength));
    leftEROutLastOrder = static_cast<float*>(malloc(sizeof(float)*ARE_IRLength));
    rightEROutLastOrder = static_cast<float*>(malloc(sizeof(float)*ARE_IRLength));
    
    float* anglesLeft = static_cast<float*>(malloc(sizeof(float)*numTotalSources));
    float* anglesRight = static_cast<float*>(malloc(sizeof(float)*numTotalSources));
    
//    ISMVector3D sourcePosition   = (ISMVector3D){ssLoc.x, ssLoc.y, ssLoc.z};
    ISMVector3D listenerLeftEar  = (ISMVector3D){lLoc.x - RADIUSOFHEAD,lLoc.y,lLoc.z};
    ISMVector3D listenerRightEar = (ISMVector3D){lLoc.x + RADIUSOFHEAD,lLoc.y,lLoc.z};

    for (size_t i = 0; i<numTotalSources; i++){
        anglesLeft[i]  = convertAzimuthToLeftEar(angleAzimuth(listenerLeftEar.x, listenerLeftEar.y,ARESources[i].x, ARESources[i].y));
        anglesRight[i] = convertAzimuthToRightEar(angleAzimuth(listenerRightEar.x, listenerRightEar.y,ARESources[i].x, ARESources[i].y));
//        printf("Room wall loc %f %f, listener left ear %f %f Angle left %f listener right ear %f %f angle right %f \n",ARESources[i].x, ARESources[i].y, lLocLE.x, lLocLE.y, anglesLeft[i], lLocRE.x, lLocRE.y, anglesRight[i]);
    }


    
      AREImpulseResponse(leftEarDelayTimes, rightEarDelayTimes,
                            leftEarGains, rightEarGains,
                            anglesLeft, anglesRight,
                            leftEROut, rightEROut,
                            leftEROutLastOrder, rightEROutLastOrder,
                            numTotalSources, ARE_IRLength,
                            numberOfFirstOrder,
                            FDN_SAMPLE_RATE);
    
    
    
    
    return numTotalSources;
    
}



/***********************************
 *    PROCESS ISM FUNCTION          *
 ***********************************/
size_t FDN::ISMInit(size_t reflectionMaxOrder){
    

    size_t* reflectionOrders;
    imageSources = static_cast<ISMVector3D*>(malloc(sizeof(ISMImageSource)*ISM_totalSourcesForOrder(reflectionMaxOrder)));
    reflectionTimesL = static_cast<float*>(malloc(sizeof(float)*ISM_totalSourcesForOrder(reflectionMaxOrder)));
    reflectionTimesR = static_cast<float*>(malloc(sizeof(float)*ISM_totalSourcesForOrder(reflectionMaxOrder)));
    reflectionGainsL = static_cast<float*>(malloc(sizeof(float)*ISM_totalSourcesForOrder(reflectionMaxOrder)));
    reflectionGainsR = static_cast<float*>(malloc(sizeof(float)*ISM_totalSourcesForOrder(reflectionMaxOrder)));
    reflectionOrders = static_cast<unsigned long*>( malloc(sizeof(ISMImageSource)*ISM_totalSourcesForOrder(reflectionMaxOrder)));
    

    // set up inputs for the simulation
    float length = roomHeight;
    float height = roomCeiling;
    float width  = roomWidth;
    float absorptionCoefficient =  sqrtf(1.f-0.1611f * Room.volume / (RT60 * Room.area));
    printf("Abs coeff in pressure %f \n", absorptionCoefficient*absorptionCoefficient);
    float speedOfSound = SOUNDSPEED;
    float micDistance = DMIN; // this is d_M
    ISMVector3D sourcePosition   = (ISMVector3D){ssLoc.x, ssLoc.y, ssLoc.z};
    ISMVector3D listenerLeftEar  = (ISMVector3D){lLoc.x - RADIUSOFHEAD,lLoc.y,lLoc.z};
    ISMVector3D listenerRightEar = (ISMVector3D){lLoc.x + RADIUSOFHEAD,lLoc.y,lLoc.z};
    
    
    // do the simulation
    ISM_simulateRoom(length, height, width,
                     absorptionCoefficient,
                     speedOfSound,
                     micDistance,
                     reflectionMaxOrder,
                     sourcePosition,
                     listenerLeftEar,
                     listenerRightEar,
                     imageSources,           // (output) source positions
                     reflectionTimesL,        // (output)
                     reflectionTimesR,
                     reflectionGainsL,        // (output)
                     reflectionGainsR,
                     reflectionOrders);      // (output) order of each I.S.
    
//    printf("Reflection Gains {");
//    for(size_t i=0; i<ISM_totalSourcesForOrder(reflectionMaxOrder); i++)
//        printf("%f, ",reflectionGainsL[i]);
//    printf("}\n\n");
//    
//    printf("Reflection Times {");
//    for(size_t i=0; i<ISM_totalSourcesForOrder(reflectionMaxOrder); i++)
//        printf("%f, ",reflectionTimesL[i]);
//    printf("}\n\n");
//    
    
    
    // find the longest reflection time
    float longestReflectionTime = 0.0f;
    size_t numSources = ISM_totalSourcesForOrder(reflectionMaxOrder);
    for(size_t i=0; i<numSources; i++){

        if(reflectionTimesL[i] > longestReflectionTime)
            longestReflectionTime = reflectionTimesL[i];
        if(reflectionTimesR[i] > longestReflectionTime)
            longestReflectionTime = reflectionTimesR[i];
    }
    
    
    
    // allocate memory for producing the impulse response
    float sampleRate = 44100;
    ISM_IRLength = 1+(longestReflectionTime*sampleRate)+HRTF_IR_LENGTH;
    leftEROut  = static_cast<float*>(malloc(sizeof(float)*ISM_IRLength));
    rightEROut = static_cast<float*>(malloc(sizeof(float)*ISM_IRLength));
    leftEROutLastOrder = static_cast<float*>(malloc(sizeof(float)*ISM_IRLength));
    rightEROutLastOrder = static_cast<float*>(malloc(sizeof(float)*ISM_IRLength));
    float* anglesLeft = static_cast<float*>(malloc(sizeof(float)*numSources));
    float* anglesRight = static_cast<float*>(malloc(sizeof(float)*numSources));
    
    
    for(size_t i=0; i<numSources; i++){
        //        anglesLeft[i]  = convertAzimuthToLeftEar(angleAzimuth(listenerLeftEar.x, listenerLeftEar.y, imageSources[i].x, imageSources[i].y));
        //        anglesRight[i] = convertAzimuthToRightEar(angleAzimuth(listenerLeftEar.x, listenerLeftEar.y, imageSources[i].x, imageSources[i].y));
        //
        anglesLeft[i]  = convertAzimuthToLeftEar(angleAzimuth(listenerLeftEar.x, listenerLeftEar.y, imageSources[i].x, imageSources[i].y));
        anglesRight[i] = convertAzimuthToRightEar(angleAzimuth(listenerRightEar.x, listenerRightEar.y, imageSources[i].x, imageSources[i].y));

        
        
    }
    
    
    ERImpulseResponse(reflectionMaxOrder, reflectionTimesL, reflectionTimesR,
                      reflectionGainsL, reflectionGainsR,
                      anglesLeft, anglesRight,
                      leftEROut, rightEROut,
                      leftEROutLastOrder, rightEROutLastOrder,
                      numSources, ISM_IRLength,
                      sampleRate);

    
    return ISM_totalSourcesForOrder(reflectionMaxOrder);
    
    
    
    
}


///************************************************////
///*********Reverb ITD Functions******************////
///************************************************////
//ITD for reverb
inline void FDN::addReverbDelay(float* fdnLeft, float*fdnRight){
    for (int  i = 0; i<CHANNELS; i++){
        if (earStatus[i] == 0){
            continue;
        }
        else if (earStatus[i] == 1){
            //add delay to left ear
            fdnLeft[i] = reverbDelays[i].process(fdnLeft[i]);
        }
        else if (earStatus[i]==2){
            //add delay to right ear
            fdnRight[i] = reverbDelays[i].process(fdnRight[i]);
        }
    }
}


void FDN::calculateAdditionalDelays(){
    
    for (int i = 0; i<CHANNELS; i++){
        tempPoints[i].z = lLoc.z; //set to be leveled
        printf("{%f, %f, %f},", tempPoints[i].x, tempPoints[i].y, tempPoints[i].z);
    }
    
    printf("\n Setting Additional Delays\n");
    for (int i = 0; i<CHANNELS; i++){
        tempPoints[i].z = lLoc.z; //set to be leveled
        
        float distanceToLeftEar = tempPoints[i].distance(lLocLE);
        float distanceToRightEar = tempPoints[i].distance(lLocRE);
//        printf("Temp point %f %f %f, distToLeftEar %f distToRightEar %f \n", tempPoints[i].x, tempPoints[i].y, tempPoints[i].z, distanceToLeftEar, distanceToRightEar);
        
       if (fabs(distanceToLeftEar-distanceToRightEar)/340.f * 44100.f< 1.0f){
            earStatus[i] = 0; //means same, no added delay
        }
        else if (distanceToRightEar > distanceToLeftEar){
            printf("i : %d add delay to right ear for %f samples \n", i, fabs((distanceToRightEar-distanceToLeftEar) / SOUNDSPEED*44100.f));
            earStatus[i] = 2;
            additionalDelays[i] = (distanceToRightEar-distanceToLeftEar) / SOUNDSPEED;
            reverbDelays[i].setTimeSafe(additionalDelays[i]);
//            reverbDelays[i].setTimeSafe((float)(i+1)*0.02);
        }
        else if (distanceToRightEar < distanceToLeftEar){
            printf("i : %d add delay to left earfor %f samples \n", i, fabs((distanceToRightEar-distanceToLeftEar) / SOUNDSPEED*44100.f));
            earStatus[i] = 1;
            additionalDelays[i] = (distanceToLeftEar-distanceToRightEar) / SOUNDSPEED;
            reverbDelays[i].setTimeSafe(additionalDelays[i]);
//            reverbDelays[i].setTimeSafe((float)(i+1)*0.02);
        }
    }
}

void FDN::setTempPoints(){

////    If FDN delay length same as proposed method
//    for (int i = 0; i<CHANNELS; i++){
//        tempPoints[i] = Room.segmentedSides[i].getMidpoint();
//    }

    //channel points
    float yBot = 0.0f-lLoc.y;
    float yTop = roomHeight - lLoc.y;
    float xLeft = 0.0f - lLoc.x;
    float xRight = roomWidth - lLoc.x;
    
    float w = lLoc.x;
    float h = lLoc.y;
    
    for (int i = 0; i < CHANNELS/2; i++){
        float angle = channeltoangleNormal(i) ;
        //        printf("Angle %f \n", angle);
        float m = 1.0f/tan(angle * M_PI / 180.f);
        //y = mx + 0
        Vector3D pointArray[4] = {Vector3D(yBot/m, yBot),
            Vector3D(yTop/m, yTop),
            Vector3D(xLeft, m*xLeft),
            Vector3D(xRight, m*xRight)};
        
        for (int j = 0; j< 4; j++){
            float xO = pointArray[j].x + lLoc.x;
            if (xO > roomWidth or xO < 0.0f){
                pointArray[j].mark = false;
                continue;
            }
            float yO = pointArray[j].y + lLoc.y;
            if (yO > roomHeight or yO < 0.0f){
                pointArray[j].mark = false;
                continue;
            }
            if (pointArray[j].mark == true){
                //check for x value
                if (pointArray[j].x >= 0){
                    tempPoints[i].x = pointArray[j].x + w;
                    tempPoints[i].y = pointArray[j].y + h;
                }
                else{
                    tempPoints[i+CHANNELS/2].x = pointArray[j].x + w;
                    tempPoints[i+CHANNELS/2].y = pointArray[j].y + h;
                }
            }
        }
    }
}


//azimuth is clockwise 0 to 180 degree right, 0 to -180 left, p. 149 DAFX
size_t FDN::determineChannel(float x,float y, float orientation){
    float xL = lLoc.x;
    float yL = lLoc.y;
    
    float xDistance = x - xL;
    float yDistance = y - yL;
    
    //Special cases
    if (fabs(yDistance - 0.00001) < 0.0f){
        // zero y, can be at 90 or -90
        if (xDistance >= 0){
            return angleToChannel(90.f);
        }
        else{
            return angleToChannel(-90.f);
        }
    }
    
    if (fabs(xDistance - 0.00001) < 0.0f){
        //zero x, can be at 0 or 180;
        if (yDistance >= 0){
            return angleToChannel(0.f);
        }
        else{
            return angleToChannel(180.f);
        }
    }
    
    //Normal case
    //x over y because we want azimuth
    float angle2 = atan2f(xDistance, yDistance) * 180.0f / M_PI;
    
    //printf("Angle %f \n", angle2);
    
    return angleToChannel(angle2);
    
}


