
//  Gains.cpp
//  surfaceIntegration
//
//  Created by Natalie Agus on 31/8/16.
//  Copyright © 2016 Hans. All rights reserved.
//

#include "Gains.hpp"
#include <iostream>
#include <algorithm>


//#define ENERGYLOSS 1.0f  //Decreasing this increases the D50 for originalDesign

#define ENERGYDETECTEDBYLISTENER 1.0f

#import <Accelerate/Accelerate.h>
using namespace std;


float Gains::pointCollectionFunction(Vector3D x, Vector3D L, Vector3D N, float visibility, float absorptionRate){
//    Vector3D xLnormalized = Lambda(x, L);
    Vector3D xL = L.subtract(x);
    float g = xL.normalize().dotProduct(N) / powf(xL.magnitude(), 2);
//    printf("xL  %f %f %f , N %f %f %f g %f x %f %f %f L %f %f %f \n", xL.x, xL.y, xL.z, N.x, N.y, N.z , g, x.x, x.y, x.z, L.x, L.y, L.z);
    return visibility * g;
}


float Gains::reflectionKernel(Vector3D x, Vector3D L, Vector3D S, Vector3D N, float visibility){
    Vector3D xSnormalized = Lambda(x, S);
    Vector3D xS = S.subtract(x);
    float g = xSnormalized.dotProduct(N) / powf(xS.magnitude(), 2);
    return visibility *  (1.f-absorptionCoefficient) / M_PI * g;
    
}

Vector3D Gains::Lambda(Vector3D u, Vector3D x){
    Vector3D ux = x.subtract(u);
    ux = ux.normalize();
    return ux;
}

void Gains::monteCarloBeta(Vector3D *points, Vector3D L, Vector3D S, Vector3D N, int numPoints, float* beta, float area, float dminimum){
    
    //Integrate h * R over the surface patch
    float hRInt = 0.0f;
    for (int i = 0; i<numPoints; i++){
//        printf("points  %f %f %f \n", points[i].x, points[i].y, points[i].z);
        float h = pointCollectionFunction(points[i], L, N, 1.0f, ALPHA);
        float R = reflectionKernel(points[i], L, S, N, 1.0f);
        hRInt += (area * h * R)/(float)numPoints * ENERGYDETECTEDBYLISTENER;
    }
     *beta = (dminimum*dminimum * hRInt);
}


void Gains:: calculateUpsilon(Plane3D* surfaces, Vector3D L, int segments, float surfaceArea){
    upsilon =  static_cast<float*>(malloc(sizeof(float)*segments));
    
    
    Vector3D **points = static_cast<Vector3D**>(malloc(segments * sizeof(Vector3D *)));
    for (int i=0; i<segments; i++)
        points[i] = static_cast<Vector3D*>(malloc(NUM_MONTECARLO * sizeof(Vector3D)));
    
    //Create random Points across the room
    for (int i = 0; i< segments; i++){
        Vector3D c2 = surfaces[i].corner;
        Vector3D s12 = surfaces[i].S1;
        Vector3D s22 = surfaces[i].S2;
        randomPointsOnRectangle(c2, s12, s22, points[i], NUM_MONTECARLO); //points on second surface

    }
    
    float averageUpsilon = 0.0f;
    
    for (int i = 0; i<segments; i++){
        float hInt = 0.0f;
        for(int j = 0; j<NUM_MONTECARLO; j++){
            //Integrate h over the surface patch
            float airAbsorption = 0.0f;
            
            airAbsorption = expf(-0.0f * points[i][j].distance(L));
            
            
            hInt += (surfaces[i].getArea() * pointCollectionFunction(points[i][j], L, surfaces[i].normal, 1.0f, 0.0f)  * airAbsorption)/(float)(NUM_MONTECARLO) * ENERGYDETECTEDBYLISTENER;
        
        }
        upsilon[i] = sqrtf(((float) segments / (M_PI * surfaceArea)) * hInt);
        averageUpsilon += upsilon[i] * upsilon[i];
    }
    
    averageUpsilon = sqrtf(averageUpsilon/segments);
    
    for (int i = 0; i<segments; i++){
        upsilon[i] /= averageUpsilon;
    }
}

//float Gains::calculateBetaFirstOrder(Plane3D *surfaces, Vector3D L, Vector3D S){
//    
//    beta = new float[ER_NUMTAPS];
//
//    Vector3D points [ER_NUMTAPS][NUM_MONTECARLO];
//
//
//    for (int i = 0; i<numberDelays; i++){
//            beta[i] = 0.0f;
//    }
//    
//        
//    for (int i = 0; i < numberDelays; i++){
//        Vector3D c = surfaces[i].corner;
//        Vector3D s1 = surfaces[i].S1;
//        Vector3D s2 = surfaces[i].S2;
// 
//        //get monte carlo points
//        randomPointsOnRectangle(c, s1, s2, points[i], NUM_MONTECARLO);
//        
//
//        //Calculate Beta
//        //Integrate h * R over the surface patch
//        float hRInt = 0.0f;
//        for (int j = 0; j<NUM_MONTECARLO; j++){
//            float airAbsorptionSx = expf(-ALPHA * S.distance(points[i][j]));
//            float airAbsorptionxL = expf(-ALPHA * L.distance(points[i][j]));
//            float h = pointCollectionFunction(points[i][j], L, surfaces[i].normal, 1.0f, ALPHA);
//            float R = reflectionKernel(points[i][j], L, S, surfaces[i].normal, 1.0f);
//            hRInt += (surfaces[i].getArea() * h * R)/(float)NUM_MONTECARLO * ENERGYDETECTEDBYLISTENER * airAbsorptionSx * airAbsorptionxL;
//        }
//        beta[i] = (dmin*dmin * hRInt);
//
//
//    }
//    
//    for (int i = 0; i<numberDelays; i++){
//        beta[i] = sqrtf(beta[i]);
//
//    }
//
//    
//
//    return 1;
//    
//    
//}


// computes the appropriate feedback gain attenuation
// to get a decay envelope with the specified RT60 time (in seconds)
// from a delay line of the specified length.
//
// This formula comes from solving EQ 11.33 in DESIGNING AUDIO EFFECT PLUG-INS IN C++ by Will Pirkle
// which is attributed to Jot, originally.
double get_gain(double rt60, double delayLengthInSamples) {
    //    printf("Delay length in samples: %f\n", delayLengthInSamples);
    return pow(10.f, (-3.0 * delayLengthInSamples) / (rt60)); //no division by 44100 here, because input delay length in samples is not multiplied by 44100
}


void Gains::calculateAllGainsSecondOrder(float*gains, Plane3D* surfaces, Vector3D L, Vector3D S,float rt60, float dminimum, size_t numberOfFirstOrder, float* gainCoefficients, Vector3D* ARESources)
{
   
    numberDelays = (int)numberOfFirstOrder;
    size_t numberSecondOrder = numberOfFirstOrder * numberOfFirstOrder;
    
//    beta = new float[numberSecondOrder];
//    upsilon = new float[numberSecondOrder];
//    gamma = new float[numberSecondOrder];
    
    
    for(int i = 0; i<numberSecondOrder; i++){
//        gains[i] = 0.0f;
//        beta[i]  = 0.0f;
//        upsilon[i] = 0.0f;
        gainCoefficients[i] = 0.0f;
    }
    
    Vector3D **points = static_cast<Vector3D**>(malloc(numberOfFirstOrder * sizeof(Vector3D *)));
    for (int i=0; i<numberOfFirstOrder; i++)
        points[i] = static_cast<Vector3D*>(malloc(NUM_MONTECARLO * sizeof(Vector3D)));
    
    
//    float gxTempValues[numberDelays*NUM_MONTECARLO];
//    float gx1x2gxTempValues[numberDelays*NUM_MONTECARLO];
//    float gx1x2gxhTempValues[numberDelays*NUM_MONTECARLO];
//    
//    // clean tempvalues
//    for (int t = 0; t<numberDelays*NUM_MONTECARLO; t++){
//        gx1x2gxTempValues[t] = 0.0f;
//        gxTempValues[t] = 0.0f;
//        gx1x2gxhTempValues[t] = 0.f;
//    }
//    
    
    //get total area of the room
    float totalArea = 0.0f;
    for (int i = 0; i<numberDelays; i++){
        totalArea += surfaces[i].getArea();
    }
    
    
    //Create random Points across the room
    for (int i = 0; i< numberDelays; i++){
        Vector3D c2 = surfaces[i].corner;
        Vector3D s12 = surfaces[i].S1;
        Vector3D s22 = surfaces[i].S2;
        randomPointsOnRectangle(c2, s12, s22, points[i], NUM_MONTECARLO); //points on second surface
    }
    
    //Handle Beta for first order
    for (int i = 0; i < numberDelays; i++){
//        printf("Surfaces %f %f %f, %f %f %f, %f %f %f \n", surfaces[i].corner.x, surfaces[i].corner.y, surfaces[i].corner.z, surfaces[i].S1.x, surfaces[i].S1.y, surfaces[i].S1.z, surfaces[i].S2.x,surfaces[i].S2.y, surfaces[i].S2.z  );
        monteCarloBeta(points[i], L, S, surfaces[i].normal, NUM_MONTECARLO, &gainCoefficients[i], surfaces[i].getArea(), dminimum);
//        printf("First order Beta: %f \n", gainCoefficients[i]);
        ARESources[i] = Vector3D(surfaces[i].getMidpoint().x, surfaces[i].getMidpoint().y, surfaces[i].getMidpoint().z);
    }
    
    
    printf("Total pure energy is: %f \n", 4 * M_PI * dmin * dmin);
    
    float totalEnergy = 0.0f;
    
    int idx = numberDelays;
    //Calculating Gamma and Beta
    for (int patch = 0; patch < numberDelays; patch++){ //to
        for (int patchBefore = 0 ; patchBefore < numberDelays; patchBefore++){ //from
            
            double patchBeta = 0.f;
            if (patch != patchBefore){
                float MonteCarloCoeff = (surfaces[patch].getArea() * surfaces[patchBefore].getArea()) /(NUM_MONTECARLO * NUM_MONTECARLO);

     
                
                for (int pointsPatch = 0 ; pointsPatch < NUM_MONTECARLO; pointsPatch ++){ //to
                    for (int pointsPatchBefore = 0 ;pointsPatchBefore < NUM_MONTECARLO; pointsPatchBefore ++){ //from
                        
                        Vector3D xS = S.subtract(points[patchBefore][pointsPatchBefore]);
                        Vector3D xPatchBefore = points[patchBefore][pointsPatchBefore];
                        Vector3D xSnormalized = Lambda(xPatchBefore, S);

                        float gx = xSnormalized.dotProduct(surfaces[patchBefore].normal) / powf(xS.magnitude(), 2);
                        float sxAirAbsorption = expf(-ALPHA * S.distance(xPatchBefore));
                        
                        Vector3D xPatch = points[patch][pointsPatch];
                        Vector3D x1x2normalized = Lambda(xPatchBefore, xPatch);
                        float dotx1x2 = x1x2normalized.dotProduct(surfaces[patchBefore].normal);
                        
                        Vector3D x2x1normalized = Lambda(xPatch, xPatchBefore);
                        float dotx2x1 = x2x1normalized.dotProduct(surfaces[patch].normal);
                        
                        Vector3D x1x2 = xPatch.subtract(xPatchBefore);
                        float gx1x2 = dotx1x2 * dotx2x1 / powf(x1x2.magnitude(), 2.f);
                        float x1x2AirAbsorption = expf(-ALPHA * xPatchBefore.distance(xPatch));
                        
                        float h = pointCollectionFunction(points[patch][pointsPatch], L, surfaces[patch].normal, 1.0f, ALPHA);
                        float x2LAirAbsorption = expf(-ALPHA * xPatch.distance(L));
                        
                        patchBeta += x2LAirAbsorption * h *  (1.f-absorptionCoefficient)/M_PI * x1x2AirAbsorption * gx1x2 *  (1.f-absorptionCoefficient)/M_PI* sxAirAbsorption * gx * dmin * dmin * MonteCarloCoeff * ENERGYDETECTEDBYLISTENER;
                        

                        
                    }
                    
                    
                }

                ARESources[idx] = Vector3D(surfaces[patch].getMidpoint().x, surfaces[patch].getMidpoint().y, surfaces[patch].getMidpoint().z);
                gainCoefficients[idx] = patchBeta;
                idx++;
                
            }
        }
        
    }
    
    for (int i = 0; i<numberSecondOrder; i++){
        if (i >=  numberOfFirstOrder)
            gainCoefficients[i] = -1*sqrtf(gainCoefficients[i]);
        else{
              gainCoefficients[i] = sqrtf(gainCoefficients[i]);
        }
//        printf("i: %d Beta %f \n", i, gainCoefficients[i]);
    }

    printf("Final idx : %d \n", idx-1);
    printf("Total energy after 2nd order : %f \n", totalEnergy);
    
}


