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
#define DENSITY_WINDOW_SIZE 882 // 20.0 * (44100.0 / 1000.0); (20 ms)
#define FDN_RT60 0.39f
#define GAINRT60 FDN_RT60
#define FDN_SAMPLE_RATE 44100.f
#define FDN_GAIN_MODEL_ON true
#define FDN_ER_SECOND_ORDER true
#define FDN_DFSR_MIN_DLY 1000
#define FDN_DFSR_MAX_DLY 4000
#define RADIUSOFHEAD 0.08f
/************* CHOOSE ONE MODE **************/

#define DIFFUSER_SAME_AS_ER true          //GOOD RESULT, using RT60 formula and solid angles
//#define DIFFUSER_SHORT_TAP true           //GOOD RESULT
//#define DIFFUSER_VELVET true              //NOT GOOD RESULT
//#define FDN_NO_DIFFUSER true              //GOOD RESULT
//#define DIFFUSER_GAMMA true               //GOOD RESULT, based on ARE, using gamma + diffuser that is shorter length, and FDN length with room dimension
//#define FDN_NO_DIFFUSER_WITH_GAMMA true   //NOT GOOD RESULT (THE 'ARE' NEEDS MORE DETAILS ON THE ROOM)



inline void FDN::setGains(){
    Room = Cuboid(8.0f, 5.0f, 3.1f);
    ssLoc = Vector3D(7.0f, 2.65f,1.2f);
    lLoc = Vector3D(7.0f, 4.0f, 1.2f);
    
    lLocLE = Vector3D(lLoc.x - RADIUSOFHEAD, lLoc.y, lLoc.z);
    lLocRE = Vector3D(lLoc.x + RADIUSOFHEAD, lLoc.y, lLoc.z);


#ifdef FDN_ER_SECOND_ORDER
    Room.segmentCube(4);//FOR SECOND ORDER, change to 9 if want more taps
    GainValues = Gains(DMIN, Room.elements, Room.area, Room.volume, FDN_RT60, true);
    printf("SECOND ORDER \n");
#else

//        Room.sliceCube(ER_NUMTAPS); //AN ALTERNATIVE WAY OF DIVIDING THE ROOM
        Room.segmentCube(ER_NUMTAPS/6); //FOR FIRST ORDER

    GainValues = Gains(DMIN, Room.elements, Room.area, Room.volume, FDN_RT60, false);
    printf("FIRST ORDER \n");
#endif
    
    
#ifdef DIFFUSER_SAME_AS_ER
//getting second order delay values
    float max = negativeInfinity;
    float min = positiveInfinity;
    for (int i = 0; i< Room.elements; i++){
        for (int j = 0; j< Room.elements; j++){
            if (i!=j && !(Room.segmentedSides[i].normal.normalize().equal(Room.segmentedSides[j].normal.normalize())) ){
                float distance = Room.segmentedSides[i].getMidpoint().distance(Room.segmentedSides[j].getMidpoint());
                
                if (distance < min){
                    min = distance;
                }
                if (distance > max){
                    max = distance;
                }
            }
        }
    }
    
    minDim = min;
    maxDim = max;
    
    printf("DIFFUSER SAME AS ER,\n minDim %f , maxDim %f \n", min , max);
    diffuserSameAsERLength = true;
    

#elif FDN_NO_DIFFUSER
    //getting second order delay values
    float max = negativeInfinity;
    float min = positiveInfinity;
    for (int i = 0; i< Room.elements; i++){
        for (int j = 0; j< Room.elements; j++){
            if (i!=j && !(Room.segmentedSides[i].normal.normalize().equal(Room.segmentedSides[j].normal.normalize())) ){
                float distance = Room.segmentedSides[i].getMidpoint().distance(Room.segmentedSides[j].getMidpoint());
                
                if (distance < min){
                    min = distance;
                }
                if (distance > max){
                    max = distance;
                }
            }
        }
    }
    
    minDim = min;
    maxDim = max;
    
    printf("FDN NO DIFFUSER,\n minDim %f , maxDim %f \n", min , max);
    diffuserSameAsERLength = false;
    
#elif DIFFUSER_GAMMA
//    //getting second order delay values
//    float max = negativeInfinity;
//    float min = positiveInfinity;
//    for (int i = 0; i< Room.elements; i++){
//        for (int j = 0; j< Room.elements; j++){
//            if (i!=j && !(Room.segmentedSides[i].normal.normalize().equal(Room.segmentedSides[j].normal.normalize())) ){
//                float distance = Room.segmentedSides[i].getMidpoint().distance(Room.segmentedSides[j].getMidpoint());
//                
//                if (distance < min){
//                    min = distance;
//                }
//                if (distance > max){
//                    max = distance;
//                }
//            }
//        }
//    }
//    
//    minDim = min;
//    maxDim = max;
//    

//    diffuserSameAsERLength = true;
    diffuserSameAsERLength = false;
    minDim = 3.5f;
    maxDim = 10.9f;
    printf("DIFFUSER GAMMA USED,\n minDim %f , maxDim %f \n", minDim , maxDim);

    
    
#else
    
    minDim = 3.1f;
    maxDim = 8.0f;


    printf("SET FDN DELAY AS ROOM DIM \n");
    diffuserSameAsERLength = false;
    
#endif
    

    
    

}





/*
 * generate randomised delay tap times between minDelay and maxDelay
 * See (http://users.spa.aalto.fi/mak/PUB/AES_Jarvelainen_velvet.pdf)
 *
 * @param minDelay    min delay time in samples
 * @param maxDelay    max delay time in samples
 * @param tapTimes    array where output will be stored
 * @param numTaps     the length of tapTimes and number of times to generate
 */
void FDN::velvetNoiseTapTimes(float minDelay,
                         float maxDelay,
                         size_t* tapTimes,
                         size_t numTaps){
    assert(maxDelay > minDelay);
    
    float minDelayF = (float)minDelay;
    float maxDelayF = (float)maxDelay;
    float numTapsF = (float)numTaps;
    float range = maxDelayF-minDelayF;
    float tapSpacing = range / numTapsF;
    randomSeed = 1145;
    updateRand();
    // seed the random number generator for consistent results
    std::srand(11);
    
    for(size_t i=0; i<numTaps; i++){
        float evenlySpacedTap = minDelay + tapSpacing * ((float)i+0.5);
        float randomJitter = tapSpacing * (float)std::rand() / (float)RAND_MAX;
        float randomlySpacedTap = evenlySpacedTap + randomJitter;
        tapTimes[i] = randomlySpacedTap;
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










void FDN::impulseResponse(long numSamples, std::ofstream* outputFile){
    
    clock_t begin = clock();

    

    
    
    float fdnOut, erOut, diffuserOut, directRayOut, nothing, input;
    float zero = 0.0f;
    float one = 1.0f;
    

    /*
     *  Gain model on: ER and diffusion run in parallel, diffusion output to the
     *  FDN input.
     */
    #ifdef FDN_GAIN_MODEL_ON
    printf("GAIN ON\n");
    
    for (int i = 0; i < numSamples; i++){
        // set the input to one for the first sample, zero otherwise.
        input = i==0 ? one : zero;
        
        // direct rays
        directRayOut = directRays.process(input * directRayAttenuation);
        
        // early reflections
        BMMultiTapDelay_processMono(&earlyReflector, input, &erOut);
        
        #ifdef FDN_NO_DIFFUSER
            // undo the beta, multiply by energy of from RT60 calculation
            float erOutUndoBeta = erOut / GainValues.averageBeta * totalRt60Energy;
        
            // pass to fdn
            processAudio(&erOutUndoBeta, &fdnOut, &nothing);
        
        #elif FDN_NO_DIFFUSER_WITH_GAMMA
            // undo the beta, multiply by energy of from RT60 calculation
            float erOutUndoBeta = erOut / GainValues.averageBeta * GainValues.averageGamma;
        
            // pass to fdn
            processAudio(&erOutUndoBeta, &fdnOut, &nothing);
        
        #else
            // diffusion
            BMMultiTapDelay_processMono(&diffuser, input, &diffuserOut);
        
            // diffusion goes into the fdn
            processAudio(&diffuserOut, &fdnOut, &nothing);
        
        #endif
        
        // scale the FDN output to represent collection at the Listener
        fdnOut *= GainValues.averageUpsilon;
        
        // write output
        *outputFile << erOut+fdnOut+directRayOut << ",";
    }
    
    
    /*
     * Gain model off, early reflections go into FDN
     */
    #else
    

    for (int i = 0; i < numSamples; i++){

        // set the input to one for the first sample, zero otherwise.
        input = i==0 ? one : zero;
        
        
        //Channel processing
        float fdnTankOutsNew[CHANNELS] = {0};
        float directRaysOutput[2] = { input *directRayAttenuation, input *directRayAttenuation };
        processTankOut(fdnTankOutsNew); //convert to 8 channels from outputsBinaural
        
        float fdnTankOutLeft[CHANNELS] = {};
        float fdnTankOutRight[CHANNELS] = {};
        filterChannels(fdnTankOutsNew, directRaysOutput, fdnTankOutLeft, fdnTankOutRight); //HRTF
        
        //get the output of multitappeddelay to outputsBinaural
//        BMMultiTapDelay_processMultiChannelOut(<#BMMultiTapDelay *This#>, <#float input#>, <#float *output#>)
        
        processDirectRays(directRaysOutput, directRaysOutput); //delays the rays
        vDSP_sve(fdnTankOutLeft, 1, &reverbOut[0], CHANNELS);
        vDSP_sve(fdnTankOutRight, 1, &reverbOut[1], CHANNELS);

        // direct rays
        directRayOut = directRays.process(input * directRayAttenuation);
        
        // early reflections
        BMMultiTapDelay_processMono(&earlyReflector, input, &erOut);
        
//        BMMultiTapDelay_processMono(&diffuser, input, &diffuserOut); //Comment this line, this is for timing
        
        // output of ER goes into fdn
        processAudio(&erOut, &fdnOut, &nothing);
        
        // write output
        *outputFile << erOut+fdnOut+directRayOut << ",";
    }
    printf("GAIN OFF \n");
    #endif
    
    clock_t end = clock();
    double elapsed_msecs = double(end - begin) / CLOCKS_PER_SEC * 1000.f;
    time_elapsed_msecs += elapsed_msecs;
    printf("Time elapsed: %f ms\n", time_elapsed_msecs);
    
    *outputFile << "\n";
    

}



//ITD for direct rays
inline void FDN::processDirectRays(float* input, float* directRaysOutput){
    directRaysOutput[0] = directRaysStereo[0].process(input[0]);
    directRaysOutput[1] = directRaysStereo[1].process(input[1]);
}

//Multiplexer
inline void FDN::processTankOut(float fdnTankOut[CHANNELS]){
    for (size_t i = 0; i < numDelays; i++){
        size_t channel = delayTimesChannel[i];
        fdnTankOut[channel] += outputsBinaural[i];
    }
}

//HRTF
inline void FDN::filterChannels(float fdnTankOut[8], float directRay[2], float fdnTankOutLeft[8], float fdnTankOutRight[8]){
    //Filter left & right ear
    for (size_t i = 0; i < CHANNELS; i++){
        fdnTankOutLeft[i] = leftEarFilter[i].process(fdnTankOut[i]);
        fdnTankOutRight[i] = rightEarFilter[i].process(fdnTankOut[i]);
    }
    //Filter direct rays
    directRay[0] = directRayFilter[0].process(directRay[0]);
    directRay[1] = directRayFilter[1].process(directRay[1]);
}



FDN::FDN(int rvType)
{
    numDelays = abs(rvType);
//    numDelays = ER_NUMTAPS;
    setGains();

//    GainValues.calculateGains(Room.segmentedSides, lLoc, ssLoc, numDelays, false); //always set to false, to make upsilon has no loss at all
    
    
    bool gainModelOn = false;
#ifdef FDN_GAIN_MODEL_ON
    gainModelOn = true;
#endif
    
    #ifdef FDN_ER_SECOND_ORDER
    
        #ifdef DIFFUSER_VELVET
            printf("DIFFUSER VELVET used, no need to calculate RT60Gains with angle\n");
    
        #elif DIFFUSER_GAMMA
            printf("DIFFUSER GAMMA used, no need to calculate RT60Gains \n");
    
    
        #else
    printf("Using allgains method here for timing \n");
    
    clock_t begin = clock();

    
    GainValues.calculateAllGainsSecondOrder(RT60GainsForMultiTapDelay, Room.segmentedSides, lLoc, ssLoc, numDelays, diffuserSameAsERLength, GAINRT60, DMIN, gainModelOn);
//            GainValues.calculateRT60GainsSecondOrder(RT60GainsForMultiTapDelay,
//                                                 ssLoc,
//                                                 Room.segmentedSides,
//                                                 GAINRT60, DMIN, lLoc, diffuserSameAsERLength);
    
    clock_t end = clock();
    double elapsed_msecs = double(end - begin) / CLOCKS_PER_SEC * 1000.f;
    time_elapsed_msecs += elapsed_msecs;

        #endif
    
    #else
        #ifdef DIFFUSER_VELVET
            printf("DIFFUSER VELVET used, no need to calculate RT60Gains with angle\n");
        #elif DIFFUSER_GAMMA
            printf("DIFFUSER GAMMA used, no need to calculate RT60Gains \n");
        #else
            GainValues.calculateRT60Gains(RT60GainsForMultiTapDelay,
                                      ssLoc,
                                      Room.segmentedSides,
                                      GAINRT60, DMIN, lLoc, diffuserSameAsERLength);
        #endif
    
    #endif
    

    totalRt60Energy = 0.0f;
    for (int i = 0; i<ER_NUMTAPS; i++){
       //printf("RT gains %d : %0.12f \n", i,  RT60GainsForMultiTapDelay[i] );
        totalRt60Energy += powf(RT60GainsForMultiTapDelay[i], 2);
        
    }
    totalRt60Energy = sqrtf(totalRt60Energy);

    directRayAttenuation =  expf(-ALPHA * ssLoc.distance(lLoc)) * DMIN/ssLoc.distance(lLoc);
    
    float delay = ssLoc.distance(lLoc) / 340.f;
    directRays.setTimeSafe(delay);
    
    #ifdef FDN_ER_SECOND_ORDER
        Room.getSecondOrderDelayValues(erDelayTapTimes, lLoc, ssLoc, FDN_SAMPLE_RATE); //FOR SECOND ORDER
    #else
        Room.getDelayValues(erDelayTapTimes, lLoc, ssLoc, FDN_SAMPLE_RATE, 0);   //FOR FIRST ORDER
    #endif
    

    
    
    /***********************************
     *    setup the early reflector    *
     ***********************************/
    // set delay tap gains to unity for early reflections
    for (int i = 0; i<ER_NUMTAPS; i++)
        delayTapGains[i] = GainValues.beta[i];
    
    // initialise the early reflections
    BMMultiTapDelay_init(&earlyReflector, false, erDelayTapTimes, 0, delayTapGains, 0, ER_NUMTAPS,0);

    
    
    /****************************
     *    setup the diffuser    *
     ****************************/
    #ifdef DIFFUSER_VELVET

        printf("Set diffuser velvet\n");
    
        float minDelay = CGFLOAT_MAX;
        float maxDelay = 0;
        float centreOfMass=0.0f;
        float totalMass = 0.0f;
        float mass;
    
        #ifdef FDN_ER_SECOND_ORDER
            for (int patchBefore = 0; patchBefore < Room.elements; patchBefore++){
                for (int patchAfter = 0; patchAfter < Room.elements; patchAfter++){
                    if (patchAfter != patchBefore && !(Room.segmentedSides[patchBefore].normal.normalize().equal(Room.segmentedSides[patchAfter].normal.normalize()))){
                        float distance = ssLoc.distance(Room.segmentedSides[patchBefore].getMidpoint()) + Room.segmentedSides[patchBefore].getMidpoint().distance(Room.segmentedSides[patchAfter].getMidpoint());
                        float timeInSamples = distance / 340.f * FDN_SAMPLE_RATE;
                        if (timeInSamples < minDelay)
                            minDelay = timeInSamples;
                        if (timeInSamples > maxDelay)
                            maxDelay = timeInSamples;
                        
                        // find the centre of mass of the output taps; an average
                        // weighted by output gain.
                        mass = 1.0f; // set gamma here to get a weighted average
                        centreOfMass += timeInSamples*mass;
                        totalMass += mass;
                    }
                }
            }
            // find the centre of mass, the weighted average delay time
            centreOfMass /= totalMass;
    
            // set the maxDelay such that the midpoint between min and max is at
            // the centreOfMass
            maxDelay = minDelay + 2*(centreOfMass-minDelay);
    
    
    
        #else

            for (int i = 0; i<ER_NUMTAPS; i++){
                float time = erDelayTapTimes[i] - (lLoc.distance(Room.segmentedSides[i].getMidpoint()) / 340.f * FDN_SAMPLE_RATE);
                if (time < minDelay){
                    minDelay = time;
                }
            if (time > maxDelay){
                maxDelay = time;
                }
            }
        #endif
    
        printf("Max length delay in distance is : %f \n", maxDelay * SOUNDSPEED / FDN_SAMPLE_RATE);
        printf("Min length delay in distance is : %f \n", minDelay * SOUNDSPEED / FDN_SAMPLE_RATE);
    
        //set diffuser delay times
        velvetNoiseTapTimes(minDelay, maxDelay, diffuserTapTimes, DIFFUSER_VELVET_NUMTAPS);
    
    
        // set delay tap gains for the diffuser based on RT60 decay time
        float* diffuserTapGains = (float*)malloc(sizeof(float)*DIFFUSER_VELVET_NUMTAPS);
        for(size_t i=0; i<DIFFUSER_VELVET_NUMTAPS; i++){
            diffuserTapGains[i] = gain(FDN_RT60, diffuserTapTimes[i], FDN_SAMPLE_RATE);

        }

    
        // normalize the diffuser tap gains to unity * rt60Decay
        float normalisation = 1.0f / sqrtf((float)DIFFUSER_VELVET_NUMTAPS) * sqrtf(4.f * M_PI * DMIN * DMIN);
        vDSP_vsmul(diffuserTapGains, 1, &normalisation, diffuserTapGains, 1, DIFFUSER_VELVET_NUMTAPS);
    
        float velvetEnergy = 0.0f;
            for(size_t i=0; i<DIFFUSER_VELVET_NUMTAPS; i++){
                velvetEnergy += diffuserTapGains[i] * diffuserTapGains[i];
            }
    
        printf("The velvet energy is: %f ", velvetEnergy);
        BMMultiTapDelay_init(&diffuser, false, diffuserTapTimes, 0, diffuserTapGains, 0, DIFFUSER_VELVET_NUMTAPS, 0);
    
        // clear the diffuser tap gains
        free(diffuserTapGains);
    
    
    #elif DIFFUSER_SAME_AS_ER
            for (int i = 0; i<DIFFUSER_NUMTAPS;i++){
                diffuserTapTimes[i] = erDelayTapTimes[i];
            }
            printf("Set diffuser same as ER\n");
    
            // set delay tap gains for the diffuser based on RT60 decay time
            float* diffuserTapGains = (float*)malloc(sizeof(float)*DIFFUSER_NUMTAPS);
            for(size_t i=0; i<DIFFUSER_NUMTAPS; i++){
                diffuserTapGains[i] = RT60GainsForMultiTapDelay[i];
            }
    
        // normalize the diffuser tap gains to unity * rt60Decay
        BMMultiTapDelay_init(&diffuser, false, diffuserTapTimes, 0, diffuserTapGains, 0, DIFFUSER_NUMTAPS, 0);
    
        // clear the diffuser tap gains
        free(diffuserTapGains);
    
    
    #elif DIFFUSER_SHORT_TAP
        printf("SHORT TAP IS ON \n");

        #ifdef FDN_ER_SECOND_ORDER
            int index = 0;
            //set as same as first order, dont wory because rt60 gains will be zero for this
            for (int i = 0; i<NUM_FIRSTORDER; i++){
                diffuserTapTimes[i] = erDelayTapTimes[i];
                    index ++;
            }
    
    
            for (int patchBefore = 0; patchBefore < Room.elements; patchBefore++){
                    for (int patchAfter = 0; patchAfter < Room.elements; patchAfter++){
                            if (patchAfter != patchBefore){
                                float distance = ssLoc.distance(Room.segmentedSides[patchBefore].getMidpoint()) + Room.segmentedSides[patchBefore].getMidpoint().distance(Room.segmentedSides[patchAfter].getMidpoint());
                                    size_t timeInSamples = distance / (size_t)SOUNDSPEED * FDN_SAMPLE_RATE;
                                    diffuserTapTimes[index] = timeInSamples;
                                    index++;
                            }
                    }
            }
    
    
    
    
            BMMultiTapDelay_init(&diffuser, false, diffuserTapTimes, 0, RT60GainsForMultiTapDelay, 0, DIFFUSER_NUMTAPS, 0);
        #else
    
            for (int i = 0; i<DIFFUSER_NUMTAPS;i++){
                    diffuserTapTimes[i] = erDelayTapTimes[i] - Room.segmentedSides[i].getMidpoint().distance(lLoc) / 340.f * 44100.f;
                    //            printf("initialTap length %d, diffuser length %d RT60gains %d %f beta %f\n", erDelayTapTimes[i], diffuserTapTimes[i], i, RT60GainsForMultiTapDelay[i], GainValues.beta[i]);
            }
    
            BMMultiTapDelay_init(&diffuser, false, diffuserTapTimes, 0, RT60GainsForMultiTapDelay, 0, DIFFUSER_NUMTAPS, 0);
        #endif
    
    
    #elif FDN_NO_DIFFUSER
        printf("No diffuser used\n");
    
    
    #elif DIFFUSER_GAMMA
        printf("Setting up diffuser with gain as gamma\n");
        #ifdef FDN_ER_SECOND_ORDER
            int index = 0;
            //set as same as first order, dont wory because rt60 gains will be zero for this
            for (int i = 0; i<NUM_FIRSTORDER; i++){
                diffuserTapTimes[i] = erDelayTapTimes[i];
                index ++;
            }
            
            
            for (int patchBefore = 0; patchBefore < Room.elements; patchBefore++){
                for (int patchAfter = 0; patchAfter < Room.elements; patchAfter++){
                    if (patchAfter != patchBefore){
                        float distance = ssLoc.distance(Room.segmentedSides[patchBefore].getMidpoint()) + Room.segmentedSides[patchBefore].getMidpoint().distance(Room.segmentedSides[patchAfter].getMidpoint());
                        size_t timeInSamples = distance / (size_t)SOUNDSPEED * FDN_SAMPLE_RATE;
                        diffuserTapTimes[index] = timeInSamples;
                        index++;
                    }
                }
            }
        #else
            for (int i = 0; i<DIFFUSER_NUMTAPS;i++){
                 diffuserTapTimes[i] = erDelayTapTimes[i] - Room.segmentedSides[i].getMidpoint().distance(lLoc) / 340.f * 44100.f;
            }
        #endif
    
//        for (int i = 0; i<DIFFUSER_NUMTAPS;i++){
//            diffuserTapTimes[i] = erDelayTapTimes[i];
//        }
    
        // set delay tap gains for the diffuser based on RT60 decay time
        float* diffuserTapGains = (float*)malloc(sizeof(float)*DIFFUSER_NUMTAPS);
        for(size_t i=0; i<DIFFUSER_NUMTAPS; i++){
            diffuserTapGains[i] = GainValues.gamma[i];
        }
    
        //set the diffuser
        BMMultiTapDelay_init(&diffuser, false, diffuserTapTimes, 0, diffuserTapGains, 0, DIFFUSER_NUMTAPS, 0);
    
        // clear the diffuser tap gains
        free(diffuserTapGains);
    
    #endif

         this->rvType = rvType;
    if (rvType > 0)
        // rvType must be a power of 2 larger than 2^4==16 for hadamard mixing
        assert(isPowerOfTwo(rvType) && rvType >= 16);
    else
        // rvType must be a multiple of 4 for sparse block circulant mixing
        assert(rvType % 4 == 0);
    
    // generate randomised delay tap outputs.
    float minDelayTime = minDim/340.f*44100.f;
    float maxDelayTime = maxDim/340.f*44100.f;

    velvetNoiseTapTimes(minDelayTime, maxDelayTime, delayTimes, numDelays);
    
    // count total delay time
    totalDelayTime = 0;
    for (size_t i = 0; i < numDelays; i++){
        totalDelayTime += delayTimes[i];
        std::cout << delayTimes[i] << ", ";
    }


    
    // randomly shuffle the order of delay times in the array
    randomPermutation1Channel(delayTimes, numDelays, 1);
    
    // declare and initialise the delay buffers
    delayBuffers = new float[totalDelayTime];
    for (int i = 0; i < totalDelayTime; i++) delayBuffers[i] = 0;
    
    // calculate start and end indices for each delay line
    startIndices[0] = 0;
    endIndices[0] = delayTimes[0];
    for (int i = 1; i < numDelays; i++){
        startIndices[i] = startIndices[i - 1] + delayTimes[i - 1];
        endIndices[i] = startIndices[i] + delayTimes[i];
    }
    
    // since we copy one input channel into many delay lines, we need to attenuate the input to avoid clipping.
    fdnMultiplexInputAttenuation = 1.0f / sqrt((float)numDelays);
    
    wetPct = .2;
    dryPct = sqrt(1.0- (wetPct*wetPct));
    
    // we specify our feedback matrix with elements of 1, and -1.  In order to keep it unitary,
    // we need to scale it down.
    // here we calculate the dividend needed to scale the feedback matrix in order to keep it unitary.

    if (rvType > 0) // fast hadamard transform
        matrixAttenuation = 1.0f / sqrt(numDelays);
    else  // 4x4 circular network
        matrixAttenuation = 1.0f / sqrt(4.0f);
    
    un_matrixAttenuation = 1.0 / matrixAttenuation;
    
    // randomise output tap signs
    for (int i = 0; i < numDelays; i++){
        tapSigns[i] = pow(-1.0, i * 31 % 71);
        //printf("tapsigns %d %f \n", i, tapSigns[i]);
    }
    
    hadamardTransform(tapSigns, temp3, numDelays);
    memcpy(tapSigns, temp3, sizeof(float) * numDelays);
    
    resetDelay();

}









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
    for (size_t i = 0; i < numDelays; i++){
        rwIndices[i] = startIndices[i];
    }
}










// larger values decay more slowly
void FDN::setDecayTime(double rt60){
    if (rt60 < 0) {
        for (int i = 0; i < numDelays; i++){
            tapGains[i] = 1.0f;
        }
    }
    else if (rt60 == 0) {
        for (int i = 0; i < numDelays; i++){
            tapGains[i] = 0.0f;
        }
    }
    else {
        for (int i = 0; i < numDelays; i++) {
            tapGains[i] = gain(rt60, delayTimes[i],FDN_SAMPLE_RATE);
//            printf("TapGains: %f delayTime %d \n", tapGains[i], delayTimes[i]);
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
    
    setDecayTime(FDN_RT60);
    std::srand(0);
}






inline void FDN::incrementIndices(){
    for (int i = 0; i < numDelays; i++){
        rwIndices[i]++;
        if (rwIndices[i] >= endIndices[i])
            rwIndices[i] = startIndices[i];
    }
}






void FDN::testProcess(long numSamples){
    float in, left, right;
    in = left = right = 0.0f;
    for (long i = 0; i < numSamples; i++){
        processAudio(&in, &left, &right);
    }
}








void FDN::processAudio(float* pInput, float* pOutputL, float* pOutputR)
{

    // Read the Input
    *pOutputL = 0;
    *pOutputR = 0;
    int i;
    float xn = *pInput * fdnMultiplexInputAttenuation;
    
    
    // mono output copied to both channels
    for (i = 0; i < numDelays; i++) {
        // copy delay line signals to output buffer and attenuate
        outputs[i] = delayBuffers[rwIndices[i]] * tapGains[i];
    }

    
    // process the feedback matrix
    if (rvType > 0) {
        
        for (i =0; i<numDelays; i++)
        {
            *pOutputL += outputs[i] ;
        }
        
        // Filling up temp3 with copies of the input sample
        vDSP_vfill(&xn, temp3, 1, numDelays);
        // Randomizing the sign of the input sample vector
        vDSP_vmul(temp3, 1, tapSigns, 1, temp3, 1, numDelays);
        // Mixing the input sample vector with the feedback vector
        vDSP_vadd(temp3, 1, outputs, 1, outputs, 1, numDelays);
        
        // Matrix mixing operation
        hadamardTransform(outputs, temp3, numDelays);
        
        // Write into the delay line inputs
        for (int i = 0; i < numDelays; i++) delayBuffers[rwIndices[i]] = temp3[i];
        
    }
    else {
        for (int i = 0; i < numDelays; i += 4){
            ppxx_xn = outputs[i] + outputs[i + 1] + xn;
            xxpp = outputs[i + 2] + outputs[i + 3];
            pnxx_xn = outputs[i] - outputs[i + 1] + xn;
            xxpn = outputs[i + 2] - outputs[i + 3];
            if (i + 7 < numDelays){
                delayBuffers[rwIndices[(i + 4)]] = ppxx_xn + xxpp; //++++
                delayBuffers[rwIndices[(i + 5)]] = ppxx_xn - xxpp; //++--
                delayBuffers[rwIndices[(i + 6)]] = pnxx_xn + xxpn; //+-+-
                delayBuffers[rwIndices[(i + 7)]] = pnxx_xn - xxpn; //+--+
            }
        }
        delayBuffers[rwIndices[0]] = ppxx_xn + xxpp; //++++
        delayBuffers[rwIndices[1]] = ppxx_xn - xxpp; //++--
        delayBuffers[rwIndices[2]] = pnxx_xn + xxpn; //+-+-
        delayBuffers[rwIndices[3]] = pnxx_xn - xxpn; //+--+
    }
    
    incrementIndices();


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
        swap(a, i, j + (rand % (length / audioChannels))*audioChannels);
    }
}






void FDN::fileResponse(long numSamples, std::ifstream* inputFile, std::ofstream* outputFile){
    float left, right;
    float zero = 0.0f;
    uint32_t samplesRead = 0;
    
    // process reverb on the file input
    wav_hdr wavHeader;
    int headerSize = sizeof(wav_hdr);
    
    //Read the header
    inputFile->read((char*)&wavHeader,headerSize);
    samplesRead += (headerSize / 2);
    
    //Read and process the audio data in the file
    byteTo16 sampleb;
    float samplef;
    while (inputFile->read((char*)&sampleb.bytes,2))
    {
        samplesRead++;
        
        // convert from int16 to float
        samplef = (float)sampleb.sixteenBit / (float)INT16_MAX;

        // process reverb on the sample
        processAudio(&samplef, &left, &right);
        *outputFile << left << "," << right << "\n";
    }
    
    // process the reverb tail, input is zero
    for (long i = samplesRead; i < numSamples; i++){
        processAudio(&zero, &left, &right);

        
        *outputFile << left << "," << right << "\n";
    }
    
    *outputFile << "\n";
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

void FDN::densityResponse(long numSamples, std::ofstream* outputFile){
    float left, right;
    float zero = 0.0f;
    float one = 1.0f;
    
    float window[DENSITY_WINDOW_SIZE];
    
    // insert the impulse
    processAudio(&one, &left, &right);
    
    // process the reverb tail
    for (int i = 0; i < numSamples; i += DENSITY_WINDOW_SIZE){
        for (int j = 0; j < DENSITY_WINDOW_SIZE; j++){
            processAudio(&zero, window + j, &right);
        }
        float threshold = max(window, DENSITY_WINDOW_SIZE) / 100.0;
        float echoDensity = (float)thresholdCount(window, DENSITY_WINDOW_SIZE, threshold) / ((float)DENSITY_WINDOW_SIZE / 44100.0);
        *outputFile << echoDensity << ",";
    }
    *outputFile << 0.0;
}

//
//// see: http://software.intel.com/en-us/articles/fast-random-number-generator-on-the-intel-pentiumr-4-processor/
//// rather than returning an output, this function updates a class variable so that we only have to generate 1 random number for each sample.
inline void FDN::updateRand()
{
    randomSeed = (214013 * randomSeed + 2531011);
    rand = (randomSeed >> 16) & 0x7FFF;
}


//Helper functions
float FDN::channeltoangleNormal(size_t channel){
    float channelWidth = 360.0f / (float)CHANNELS;
    return (float)channel * channelWidth + (channelWidth / 2.0f);
}

float FDN::channelToAngle(size_t channel){
    float channelWidth = 360.0f / (float)CHANNELS;
    size_t midChannel = CHANNELS / 2 - 1;
    // return (float)channel * channelWidth + (channelWidth / 2.0f);
    if (channel <= midChannel){
        return (float)channel * channelWidth + (channelWidth / 2.0f);
    }
    else{
        return (-1.f  * (float) (CHANNELS - channel) * channelWidth )+ (channelWidth / 2.0f);
    }
}

size_t FDN::angleToChannel(float angleInDegrees){
    
    float channelWidth = 360.0f / CHANNELS;
    
    //Channel 0 to 7 from 12 o'clock, clockwise according from DFX
    return (size_t)floorf(angleInDegrees / channelWidth);
}

//azimuth is clockwise 0 to 360 degree, p. 149 DAFX
size_t FDN::determineChannel(float x,float y){
    float xL = lLoc.x;
    float yL = lLoc.y;
    
    float xDistance = x - xL;
    float yDistance = y - yL;
    
    float angle2 = atan2f(xDistance, yDistance) * 180.0f / M_PI;
    if (angle2 < 0.0f ){
        angle2 += 360.0f;
    }
    return angleToChannel(angle2);
    
}

void FDN::setHRTFFilters(){
    //Right is negative, p.151 DAFX
    for (size_t i = 0; i < CHANNELS; i++){
        leftEarFilter[i].setAngle(-channelToAngle(i), FDN_SAMPLE_RATE, false);
        rightEarFilter[i].setAngle(channelToAngle(i), FDN_SAMPLE_RATE, true);
    }
    
    directRayFilter[0].setAngle(-directRayAngles[0], FDN_SAMPLE_RATE,false);
    directRayFilter[1].setAngle(directRayAngles[1], FDN_SAMPLE_RATE, true);
    
}

void FDN::setDirectSingleTapDelay(){
    directRaysStereo[0].setTimeSafe(directDelayTimes[0]);
    directRaysStereo[1].setTimeSafe(directDelayTimes[1]);
}


inline void FDN::setDirectDelayTimes()
{
    
    //Calculate delay from source to receiver
    float directDelayLeft = ssLoc.distance(lLocLE) / SOUNDSPEED;
    float directDelayRight = ssLoc.distance(lLocRE) / SOUNDSPEED;
    directDelayTimes[0] = directDelayLeft;
    
    directDelayTimes[1] = directDelayRight;
}

void FDN::setDelayChannels(){
    
    for (size_t i = 0; i < ER_NUMTAPS; i++){
        delayTimesChannel[i] = determineChannel(Room.segmentedSides[i].getMidpoint().x, Room.segmentedSides[i].getMidpoint().y);
    }
    delayTimes[ER_NUMTAPS] = 0;
    
    
    
    
}
