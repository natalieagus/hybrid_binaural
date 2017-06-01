//
//  Gains.hpp
//  surfaceIntegration
//
//  Created by Natalie Agus on 31/8/16.
//  Copyright © 2016 Hans. All rights reserved.
//

#ifndef Gains_hpp
#define Gains_hpp

#include "Cuboid.hpp"
#include "Vector3D.hpp"
#include "Plane3D.hpp"
#include "Cuboid.hpp"
#include <math.h>
#include <stdio.h>
//#import "energyTerm.hpp"
#include "MonteCarlo.hpp"
#include "TransposeGain.hpp"

//#define DIFFUSER_NUMTAPS 576//6*4*4 //128 //6*4*4//576 //CHANGE THIS TO 36 FOR SECOND ORDER
//#define FDN_NUMDELAYS DIFFUSER_NUMTAPS
//#define ER_NUMTAPS 36//2916 //576
//#define DIFFUSER_NUMTAPS ER_NUMTAPS
//#define DELAYSFDN 20
#define NUM_MONTECARLO 50

//#define NUM_FIRSTORDER 6//54 //24
#define DIFFUSER_VELVET_NUMTAPS ER_NUMTAPS


//change them to make the first order reflections lose energy if you have to. 
#define KD 1.0f
#define KS 0.0f
#define ENERGYINITIAL 1.0f
#define ALPHA 0.005f //0.0 is absorb nothing, Nepers per meter, a/20log10e, a = decibels per meter

typedef struct Gains{
    
    Gains(){
    
    }
    //Constructor with arguments
    Gains(float dmin, int numberDelays, float totalSurfaceArea, float totalVolume, float RT60, bool secondOrder){
        this->dmin = dmin ;
        this->numberDelays = numberDelays;
        this->totalSurfaceArea = totalSurfaceArea;
        this->totalInputEnergy = 0.f;
        this->correctInputEnergy = 4.f * M_PI * powf(dmin, 2) * ENERGYINITIAL;
        this->secondOrder = secondOrder;
        this->rt60 = RT60;
        this->totalVolume = totalVolume;
        
        //sabine's equation
        absorptionCoefficient = 0.1611f * this->totalVolume / (this->rt60 * this->totalSurfaceArea);
        printf("Energy left is : %f \n", 1.f-absorptionCoefficient);
        //uncomment this if we wanna deactivate absorption coefficient
//        absorptionCoefficient = 0.f;
        
        printf("Abs coeff : %f \n", absorptionCoefficient);
    };
    
    
    ////****METHODS****////
    void calculateUpsilon(Plane3D* surfaces, Vector3D L, int segments, float surfaceArea);
    void calculateAllGainsSecondOrder(float*gains, Plane3D* surfaces, Vector3D L, Vector3D S,float rt60, float dminimum, size_t numberOfFirstOrder, float* gainCoefficients, Vector3D* ARESources);
    float calculateBetaFirstOrder(Plane3D* surfaces, Vector3D L, Vector3D S);
    void monteCarloBeta(Vector3D* points, Vector3D L, Vector3D S, Vector3D N, int numPoints, float* beta, float area, float dminimum);
    float pointCollectionFunction(Vector3D x, Vector3D L, Vector3D N, float visibility, float absorption);
    float reflectionKernel(Vector3D x, Vector3D L, Vector3D S, Vector3D N, float visibility);
    Vector3D Lambda(Vector3D u, Vector3D x);


    ////****VARIABLES****//
    bool secondOrder;
    float absorptionCoefficient;
    float dmin;
    float totalSurfaceArea;
    int numberDelays;
    float* upsilon;
    float* beta;
    float* gamma;
    float rt60;
    float totalVolume;

    float totalInputEnergy;
    float correctInputEnergy;
    
    
}Gains;

#endif /* Gains_hpp */
