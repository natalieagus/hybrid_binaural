
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
    Vector3D xLnormalized = Lambda(x, L);
    Vector3D xL = L.subtract(x);
//    printf("x %f %f %f, L %f %f %f, xL %f %f %f norm %f\t", x.x, x.y, x.z, L.x, L.y, L.z, xL.x, xL.y, xL.z);
    float g = xLnormalized.dotProduct(N) / powf(xL.magnitude(), 2);
    return visibility * g;
}

void Gains::monteCarloGammaFirstOrder(Vector3D* points, Vector3D S, Vector3D N, float surfaceArea){
    float gInt = 0;
    
    float monteCarloCoeff = surfaceArea / (float)(NUM_MONTECARLO);
    

    //    printf("area/numpoints %f \n", area/numPoints);
    for (int i = 0; i< NUM_MONTECARLO; i++){
        
        Vector3D xSnormalized = Lambda(points[i], S);
        Vector3D xS = S.subtract(points[i]);
        float g = xSnormalized.dotProduct(N) / powf(xS.magnitude(), 2);
        float airAbsorption = expf(-ALPHA * S.distance(points[i]));
        gInt += g * dmin * dmin * (1.f-absorptionCoefficient) * airAbsorption * monteCarloCoeff;
    }
    
//    *gamma = (dmin * dmin * gInt);
    averageGamma += gInt;
    
    
}



//float Gains::getDBRDF(){
//    return 1.0f/M_PI * ENERGYLOSS;
//}


float Gains::reflectionKernel(Vector3D x, Vector3D L, Vector3D S, Vector3D N, float visibility){
    Vector3D xSnormalized = Lambda(x, S);
    Vector3D xS = S.subtract(x);
    float g = xSnormalized.dotProduct(N) / powf(xS.magnitude(), 2);
    return visibility *  (1.f-absorptionCoefficient) / M_PI * g;
    
}

Vector3D Gains::Lambda(Vector3D u, Vector3D x){
    Vector3D ux = x.subtract(u);
//               printf("ux %f %f %f  x %f %f %f u %f %f %f \t", ux.x, ux.y, ux.z, x.x, x.y, x.z, u.x, u.y, u.z);
    ux = ux.normalize();
    return ux;
}




void Gains::monteCarloBeta(Vector3D *points, Vector3D L, Vector3D S, Vector3D N, int numPoints, float* beta, float area){
    
    //Integrate h * R over the surface patch
    float hRInt = 0.0f;
    for (int i = 0; i<numPoints; i++){
        float h = pointCollectionFunction(points[i], L, N, 1.0f, ALPHA);
        float R = reflectionKernel(points[i], L, S, N, 1.0f);
        hRInt += (area * h * R)/(float)numPoints * ENERGYDETECTEDBYLISTENER;
    }
    

     *beta = (dmin*dmin * hRInt);
    
}





float Gains::calculateGains(Plane3D *surfaces, Vector3D L, Vector3D S, int numFDNdelays, bool useERLengthGamma){
    
//    printf("SourceLoc : %f %f %f lLoc %f %f %f \n", S.x, S.y, S.z, L.x, L.y, L.z);
//        printf("dmin %f \n", dmin);
    averageUpsilon = 0.0f;
    averageGamma = 0.0f;
    averageBeta = 0.0f;
    
//    printf("Calculating gains...");

    
    beta = new float[ER_NUMTAPS];
    upsilon = new float[ER_NUMTAPS];
    gamma = new float[ER_NUMTAPS];
    

    
//    printf("Number delays first order: %d \n", numberDelays);
    
    
    if (secondOrder){
        
        for (int i = 0; i<ER_NUMTAPS; i++){
            beta[i] = 0.0f;
            upsilon[i] = 0.0f;
            gamma[i] = 0.0f;
        }
        
    //Handle Beta for second order + totalGamma + totalUpsilon
        monteCarloSecondOrder(surfaces, L, S, useERLengthGamma);
//square root the beta, and upsilon
        averageUpsilon = sqrtf(averageUpsilon);

        
        averageGamma = 0.0f;
        averageBeta = 0.0f;
        
        for (int i = 0; i<ER_NUMTAPS; i++){
            averageGamma += gamma[i];
            averageBeta += beta[i];
            beta[i] = sqrtf(beta[i]);
//            printf("beta %d is %f , ", i,beta[i]);
            gamma[i] = sqrtf(gamma[i]);
            
        }
        printf("average gamma in energy %f\n", averageGamma );
        printf("average beta %f \n", averageBeta);
        averageGamma = sqrtf(averageGamma);
        averageBeta = sqrtf(averageBeta);
//        printf("average gamma in energy\n", averageGamma );
    }
    

    
    else{
        
        Vector3D points [ER_NUMTAPS][NUM_MONTECARLO];
 //=======FOR FIRST ORDER =======

        for (int i = 0; i<numberDelays; i++){
            beta[i] = 0.0f;
            upsilon[i] = 0.0f;
            gamma[i] = 0.0f;
        }
    
        
    for (int i = 0; i < numberDelays; i++){
        Vector3D c = surfaces[i].corner;
        Vector3D s1 = surfaces[i].S1;
        Vector3D s2 = surfaces[i].S2;
//        printf("{{%f, %f, %f}, {%f, %f ,%f}, {%f, %f ,%f} } normal {%f %f %f},\n", c.x, c.y, c.z, s1.x, s1.y, s1.z, s2.x, s2.y, s2.z, surfaces[i].normal.x, surfaces[i].normal.y, surfaces[i].normal.z);
        
        //get monte carlo points
        randomPointsOnRectangle(c, s1, s2, points[i], NUM_MONTECARLO);
        
//        for (int j =0; j<NUM_MONTECARLO; j++){
//      printf("{%f, %f, %f}, ", points[i][j].x, points[i][j].y, points[i][j].z);
//        }

        //Calculate Beta
        //Integrate h * R over the surface patch
        float hRInt = 0.0f;
        for (int j = 0; j<NUM_MONTECARLO; j++){
            float airAbsorptionSx = expf(-ALPHA * S.distance(points[i][j]));
            float airAbsorptionxL = expf(-ALPHA * L.distance(points[i][j]));
            float h = pointCollectionFunction(points[i][j], L, surfaces[i].normal, 1.0f, ALPHA);
            float R = reflectionKernel(points[i][j], L, S, surfaces[i].normal, 1.0f);
//            printf("beta R %f  point %f %f %f \n ", h, points[i][j].x, points[i][j].y, points[i][j].z);
            hRInt += (surfaces[i].getArea() * h * R)/(float)NUM_MONTECARLO * ENERGYDETECTEDBYLISTENER * airAbsorptionSx * airAbsorptionxL;
        }
        beta[i] = (dmin*dmin * hRInt);
//        printf("Beta i %f \n",i, sqrtf(beta[i]));
 
        averageBeta += (dmin*dmin * hRInt);
//        printf("beta ...");
        
        //Calculate Gamma
        float gInt = 0.f;
        for (int j = 0; j< NUM_MONTECARLO; j++){
            Vector3D xSnormalized = Lambda(points[i][j], S);
            Vector3D xS = S.subtract(points[i][j]);
            float g = xSnormalized.dotProduct(surfaces[i].normal) / powf(xS.magnitude(), 2);
            float airAbsorption = expf(-ALPHA * S.distance(points[i][j]));
            gInt += g * dmin * dmin * 1.0f * airAbsorption * surfaces[i].getArea() / (float) NUM_MONTECARLO *  (1.f-absorptionCoefficient); //have one energy loss here
        }
        gamma[i] = gInt;
        averageGamma += gInt;
//        printf("gamma ...");
        
        //Calculate Upsilon
        float hInt = 0.f;
        for(int j = 0; j<NUM_MONTECARLO; j++){
//                    printf("Point : %f %f %f \t",points[i][j].x, points[i][j].y, points[i][j].z );
            //Integrate h over the surface patch
            
            float airAbsorption = 0.0f;
            if (!useERLengthGamma)
                airAbsorption = expf(-0.0f * points[i][j].distance(L));
            else
                airAbsorption = expf(ALPHA * points[i][j].distance(L));
            hInt += (surfaces[i].getArea() * pointCollectionFunction(points[i][j], L, surfaces[i].normal, 1.0f, 0.0f) * airAbsorption)/((float)(NUM_MONTECARLO))  ;
//            printf("PCF %f  \t", pointCollectionFunction(points[i][j], L, surfaces[i].normal, 1.0f, 0.0f));

//            hInt += (surfaces[i].getArea() * pointCollectionFunction(points[i][j], L, surfaces[i].normal, 1.0f, 0.0f) * airAbsorption)/(float)(NUM_MONTECARLO) * ENERGYDETECTEDBYLISTENER;
            
        }
        upsilon[i] =  numFDNdelays * hInt  * ( 1.0f / ( totalSurfaceArea) )* 1.0f /M_PI;

        averageUpsilon += upsilon[i];
        
//        printf("Surface %d done\n", i);

    }
    
    totalInputEnergy = dmin * dmin * 4*M_PI;
//    printf("Total Energy %f \n", totalInputEnergy);
    

        
    for (int i = 0; i<numberDelays; i++){
        beta[i] = sqrtf(beta[i]);
        upsilon[i] = sqrtf(upsilon[i]);
        gamma[i] = sqrtf(gamma[i]);
        //        beta[i] = 1.0f;
//        upsilon[i] =sqrtf(upsilon[i]);
        //        printf("upsilon %d : %f \n",i, upsilon[i]);
    }
        printf("average gamma in energy %f\n", averageGamma );
//    printf("Avg upsilon %f Avg beta %f Avg gamma%f \n", averageUpsilon/numberDelays, averageBeta, averageGamma);
    averageBeta = sqrtf(averageBeta);
    averageGamma = sqrtf(averageGamma);
        

//    float ratio = transposeEnergyRatio(upsilon, (size_t) ER_NUMTAPS);
//    averageUpsilon = sqrtf(averageUpsilon/numFDNdelays * ratio); // to be applied at total output of FDN : k y^2, hence must be divided by K

    averageUpsilon = sqrtf(averageUpsilon/numFDNdelays);

    
    //    printf("averageGamma * averageUpsilon * 1/sqrtf(MUL)  %f \n", averageGamma * averageUpsilon * 1/sqrtf(MULTITAPDELAYS));

    
//    ====END OF FOR FIRST ORDER

    

    
    }

    return 1;
    
    
}


//Calculates gamma, beta, and upsilon for up to second order reflection
void Gains::monteCarloSecondOrder(Plane3D *surfaces, Vector3D L, Vector3D S, bool useERLengthGamma){

    averageGamma = 0.0f;
    
    Vector3D points [NUM_FIRSTORDER][NUM_MONTECARLO];

    
    float gxTempValues[numberDelays*NUM_MONTECARLO];
    float gx1x2gxTempValues[numberDelays*NUM_MONTECARLO];
    float gx1x2gxhTempValues[numberDelays*NUM_MONTECARLO];
    
    // clean tempvalues
    for (int t = 0; t<numberDelays*NUM_MONTECARLO; t++){
        gx1x2gxTempValues[t] = 0.0f;
        gxTempValues[t] = 0.0f;
        gx1x2gxhTempValues[t] = 0.f;
    }
    
    
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
        monteCarloBeta(points[i], L, S, surfaces[i].normal, NUM_MONTECARLO, &beta[i], surfaces[i].getArea());
    }

    
    //Calculate Upsilon

    for (int i = 0; i<numberDelays; i++){
        float hInt = 0.0f;
        for(int j = 0; j<NUM_MONTECARLO; j++){
            //Integrate h over the surface patch
            float airAbsorption = 0.0f;
            if (!useERLengthGamma)
                airAbsorption = expf(-0.0f * points[i][j].distance(L));
            else
                airAbsorption = expf(-ALPHA * points[i][j].distance(L));
            
            hInt += (surfaces[i].getArea() * pointCollectionFunction(points[i][j], L, surfaces[i].normal, 1.0f, 0.0f)  * airAbsorption)/(float)(NUM_MONTECARLO) * ENERGYDETECTEDBYLISTENER;
            
//            hInt += (surfaces[i].getArea() * pointCollectionFunction(points[i][j], L, surfaces[i].normal, 1.0f, 0.0f))/(float)(NUM_MONTECARLO) * ENERGYDETECTEDBYLISTENER;
            //                   printf("HInt %f ", hInt);
        }
        averageUpsilon += (( 1.0f/ ( totalSurfaceArea)) * 1.0/M_PI * hInt);
 
        
    }
//
//    printf("Average Upsilon: %f \n", averageUpsilon);
    

    int idx = numberDelays;
    //Calculating Gamma and Beta
    for (int patch = 0; patch < numberDelays; patch++){ //to
        for (int patchBefore = 0 ; patchBefore < numberDelays; patchBefore++){ //from
            
            double patchBeta = 0.f;
            double patchGamma = 0.f;
            if (patch != patchBefore){

                float MonteCarloCoeff = (surfaces[patch].getArea() * surfaces[patchBefore].getArea()) /(NUM_MONTECARLO * NUM_MONTECARLO);
//                printf("From patch : %d to patch %d, monte carlo coeff is %f \n", patchBefore, patch, MonteCarloCoeff);
                
                for (int pointsPatch = 0 ; pointsPatch < NUM_MONTECARLO; pointsPatch ++){ //to
                    for (int pointsPatchBefore = 0 ;pointsPatchBefore < NUM_MONTECARLO; pointsPatchBefore ++){ //from
                    

                        Vector3D xPatchBefore = points[patchBefore][pointsPatchBefore];
                        Vector3D xSnormalized = Lambda(xPatchBefore, S);
                        Vector3D xS = S.subtract(xPatchBefore);
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
                        

                        patchGamma += gx1x2 * x1x2AirAbsorption *  (1.f-absorptionCoefficient)/M_PI * sxAirAbsorption * gx * dmin * dmin  * MonteCarloCoeff *  (1.f-absorptionCoefficient);

                        float h = pointCollectionFunction(points[patch][pointsPatch], L, surfaces[patch].normal, 1.0f, ALPHA);
                        float x2LAirAbsorption = expf(-ALPHA * xPatch.distance(L));
                        
                        patchBeta += x2LAirAbsorption * h *  (1.f-absorptionCoefficient)/M_PI * x1x2AirAbsorption * gx1x2 *  (1.f-absorptionCoefficient)/M_PI* sxAirAbsorption * gx * dmin * dmin * MonteCarloCoeff * ENERGYDETECTEDBYLISTENER;
                        
//                        printf("delta gamma %0.12f, delta beta %0.12f  \n", gx1x2 * x1x2AirAbsorption * ENERGYLOSS/M_PI * sxAirAbsorption * gx * dmin * dmin  * MonteCarloCoeff, x2LAirAbsorption * h * ENERGYLOSS/M_PI * x1x2AirAbsorption * gx1x2 * ENERGYLOSS/M_PI* sxAirAbsorption * gx * dmin * dmin * MonteCarloCoeff);
                    
                }
         
            
                }
                gamma[idx] = patchGamma;
                beta[idx] = patchBeta;
                idx++;
    
            }
        }
    
    }

    
//    printf("total gamma : %0.10f\n", averageGamma );
//
//    for (int i = 0; i<MULTITAPDELAYS; i++)
//    {
//        printf("i %d beta: %0.10f \n", i, beta[i]);
//    }
////
//    printf("total upsilon: %0.10f \n", averageUpsilon);
    
}






float Gains::checkEnergy(){
//    float totalUp = 0.0f;
    averageBeta = 0.0f;
//    float totalGamma = 0.0f;
    for (int i = 0; i<numberDelays; i++){
//        totalUp += upsilon[i] ;
        averageBeta += beta[i] ;

    }
    
    float lateEnergy = averageGamma * averageUpsilon;

    float difference = lateEnergy - averageBeta;
    
    averageBeta = sqrtf(averageBeta);
    
    printf("Average gamma, %f, average Beta, %f, average Upsilon: %f \n", averageGamma, averageBeta, averageUpsilon);
    return abs(difference);
    
}


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


void Gains::calculateRT60Gains(float* gains, Vector3D S, Plane3D* surfaces, float rt60, float dmin, Vector3D L, bool useERLength){
    
    Vector3D points[NUM_MONTECARLO];
    
    float totalEnergy = 0.0f;
    
    for (int i = 0; i<numberDelays; i++){
        //get the surface information
        Vector3D c = surfaces[i].corner;
        Vector3D s1 = surfaces[i].S1;
        Vector3D s2 = surfaces[i].S2;
        //        printf("{{%f, %f, %f}, {%f, %f ,%f}, {%f, %f ,%f} } normal {%f %f %f},\n", c.x, c.y, c.z, s1.x, s1.y, s1.z, s2.x, s2.y, s2.z, surfaces[i].normal.x, surfaces[i].normal.y, surfaces[i].normal.z);
        
        //get monte carlo points
        randomPointsOnRectangle(c, s1, s2, points, NUM_MONTECARLO);
        
        //calculate through the points
        float gainInt = 0.0f;
        float normalizationTerm = (surfaces[i].getArea() / (float)NUM_MONTECARLO);
        
        for(int j = 0; j<NUM_MONTECARLO; j++){
            float r = S.distance(points[j]);
            float rSquared = powf(r,2);
          
            float rSL = 0.0f;
        if(useERLength)
             rSL = r + L.distance(points[j]);
        else
             rSL = r;
            
            
            //get the dot product
            Vector3D xS = S.subtract(points[j]);
            xS = xS.normalize();
            float dotProduct = xS.dotProduct(surfaces[i].normal);
            
            float gainReduction = powf(get_gain(rt60, (r+rSL)/343.f),2);

            gainInt += (gainReduction / rSquared) * dotProduct * normalizationTerm * dmin * dmin;
        }
        totalEnergy += gainInt;
        gains[i] = sqrtf(gainInt);
        
    }
    
    
    printf("Total energy by RT60decay is : %f \n", totalEnergy);
    
}


void Gains::calculateRT60GainsSecondOrder(float* gains, Vector3D S, Plane3D* surfaces, float rt60, float dmin, Vector3D L, bool useERLength){
    
    Vector3D points[NUM_FIRSTORDER][NUM_MONTECARLO];
    
    //Create random Points across the room
    for (int i = 0; i< numberDelays; i++){
        Vector3D c2 = surfaces[i].corner;
        Vector3D s12 = surfaces[i].S1;
        Vector3D s22 = surfaces[i].S2;
        randomPointsOnRectangle(c2, s12, s22, points[i], NUM_MONTECARLO); //points on second surface
    }
    
    for(int i = 0; i<NUM_FIRSTORDER; i++){
        gains[i] = 0.0f;
    }
    
//    float tempEnergy[NUM_FIRSTORDER];
//    float totalEnergy = 0.0f;
//    //calculate first bounce energy on each patch
//    for (int i = 0; i<numberDelays; i++){
//        //get the surface information
//        //calculate through the points
//        float gainInt = 0.0f;
//        float normalizationTerm = (surfaces[i].getArea() / (float)NUM_MONTECARLO);
//        
//        for(int j = 0; j<NUM_MONTECARLO; j++){
//            float r = S.distance(points[i][j]);
//            float rSquared = powf(r,2);
//            
//            
//            float rSL = r ;//+ L.distance(points[j]);
//            
//            //get the dot product
//            Vector3D xS = S.subtract(points[i][j]);
//            xS = xS.normalize();
//            float dotProduct = xS.dotProduct(surfaces[i].normal);
//            
//            float gainReduction = powf(get_gain(rt60, (r+rSL)/343.f * 44100.f),2);
//            
//            gainInt += (gainReduction / rSquared) * dotProduct * normalizationTerm * dmin * dmin;
//        }
//        
//        tempEnergy[i] = gainInt;
//        totalEnergy += gainInt;
////        printf("TempEnergy : %f \n", tempEnergy[i]);
//    }
//    printf("TotalEnergy : %f \n", totalEnergy);
//    //calculate energy at the second patch
//    
    
    printf("Total pure energy is: %f \n", 4 * M_PI * dmin * dmin);
    
    int index = NUM_FIRSTORDER;
    float totalEnergy = 0.0f;
//    totalEnergy = 0.0f;
    //iterate through the points
    for (int patchbefore = 0; patchbefore < numberDelays; patchbefore++){
        for (int patchafter = 0; patchafter < numberDelays; patchafter++){
            if (patchafter != patchbefore){
                double gainInt = 0.0f;
                
//                float monteCarloNormalizationTerm = surfaces[patchafter].getArea() / (NUM_MONTECARLO);
                float monteCarloNormalizationTerm = surfaces[patchafter].getArea() * surfaces[patchbefore].getArea() / (NUM_MONTECARLO * NUM_MONTECARLO);
                for(int i = 0; i < NUM_MONTECARLO; i++){
                    for (int j = 0 ; j<NUM_MONTECARLO ;j++){
                        
                        float distanceSx = S.distance(points[patchbefore][i]);
                        float distanceSquaredSx = powf(distanceSx, 2);
                        float gainReductionSx = powf(get_gain(rt60, distanceSx / 343.f), 2);

                        Vector3D xS = S.subtract(points[patchbefore][i]);
                        xS = xS.normalize();
                        float dotProductSx1 = xS.dotProduct(surfaces[patchbefore].normal);
                        
                        
                        
                        float dotProductx1x2 = points[patchbefore][i].subtract(points[patchafter][j]).normalize().dotProduct(surfaces[patchafter].normal);
                        
                        float dotProductx2x1 = points[patchafter][j].subtract(points[patchbefore][i]).normalize().dotProduct(surfaces[patchbefore].normal);
                        
                        float distancexx = 0.0f;
                        
                        if (useERLength)
                            distancexx = points[patchbefore][i].distance(points[patchafter][j]) + points[patchafter][j].distance(L);
                            //this is like saying the spread amount from that point x on the wall to any location thats distance L away from that point on the wall,
                        //so its all the energy in the enclosed hemisphere (which is energy at that point x) but decayed by distance L according to the RT60 formula
                        //it is integration of dot product from x to dL, with energy Ex/pi * gainRt60, and the integration is pi, so we just have Ex*gainRT60
                        else
                            distancexx = points[patchbefore][i].distance(points[patchafter][j]);
                        
                        
                        float distanceSquaredxx = powf( points[patchbefore][i].distance(points[patchafter][j]), 2);
                        float gainReductionxx = powf(get_gain(rt60, distancexx / 343.f), 2);
                        
//                        gainInt  +=   (gainReductionxx / distanceSquaredxx) * dotProductx1x2 * (1.f/(2*M_PI)) * (gainReductionSx / distanceSquaredSx) * dotProductSx1 * dmin * dmin * monteCarloNormalizationTerm;
                        
                        gainInt  +=   (gainReductionxx / distanceSquaredxx) * dotProductx2x1 * dotProductx1x2 * (1.f/(M_PI)) * (gainReductionSx / distanceSquaredSx) * dotProductSx1 * dmin * dmin * monteCarloNormalizationTerm;
                        //use the above because on a surface with a normal, you technically do not spread energy to the direction perpendicular to the normal, so it is division by pi and having the dot product all over.
                        
                        
                        
//                        printf("gainreduction xx %f distanceSquared xx %f dotProductx1x2 %f gainReductionSx %f distanceSquaredSx %f dotProductSx1 %f montecarlonormterm %f  dmin %f Total : %f \n", gainReductionxx, distanceSquaredxx, dotProductx1x2, gainReductionSx, distanceSquaredSx, dotProductSx1, monteCarloNormalizationTerm, dmin,  (gainReductionxx / distanceSquaredxx) * dotProductx1x2 * (1.f/(2*M_PI)) * (gainReductionSx / distanceSquaredSx) * dotProductSx1 * dmin * dmin * monteCarloNormalizationTerm);
                        
//                        float distancexx = 0.0f;
//                        if (useERLength)
//                            distancexx = surfaces[patchbefore].getMidpoint().distance(points[patchafter][j]) + points[patchafter][j].distance(L);
//                        else
//                            distancexx = surfaces[patchbefore].getMidpoint().distance(points[patchafter][j]);
//                        
//                        
//                        float distanceSquaredxx = powf(surfaces[patchbefore].getMidpoint().distance(points[patchafter][j]), 2);
//                        float gainReductionxx = powf(get_gain(rt60, distancexx / 343.f * 44100.f), 2);
//
//                        
//                        
//                        
//                        float distance = S.distance(points[patchbefore][i]) + points[patchbefore][i].distance(points[patchafter][j]);// + points[patchafter][j].distance(L);
//                        float distanceSquared = powf(distance, 2);
                        //                             float dotProductx1x2 = surfaces[patchbefore].getMidpoint().subtract(points[patchafter][j]).normalize().dotProduct(surfaces[patchafter].normal);
                        
//

//                        printf("%f %f %f\n", points[patchbefore][i].subtract(points[patchafter][j]).normalize().x, points[patchbefore][i].subtract(points[patchafter][j]).normalize().y, points[patchbefore][i].subtract(points[patchafter][j]).normalize().z);
//
                    
//                        gainInt += tempEnergy[patchbefore] / (2 * M_PI) * (gainReductionxx / distanceSquaredxx * dotProductx1x2) * monteCarloNormalizationTerm ;
  
//                        printf("gainInt %f \n", gainInt);

                    }
                }
                
                
                totalEnergy += gainInt;
                gains[index] = sqrtf(gainInt);
                index ++;
            }

            
      
        }
    }
    
    printf("Total energy after 2nd order : %f \n", totalEnergy);
    
    
    
}

void Gains::calculateAllGainsSecondOrder(float*gains,Plane3D* surfaces, Vector3D L, Vector3D S, int numFDNdelays, bool useERLength,float rt60, float dmin, bool gainOn)
{
    
    //    printf("SourceLoc : %f %f %f lLoc %f %f %f \n", S.x, S.y, S.z, L.x, L.y, L.z);
    //        printf("dmin %f \n", dmin);
    averageUpsilon = 0.0f;
    averageGamma = 0.0f;
    averageBeta = 0.0f;
    
    //    printf("Calculating gains...");
    
//    Vector3D points [ER_NUMTAPS][NUM_MONTECARLO];
    
    beta = new float[ER_NUMTAPS];
    upsilon = new float[ER_NUMTAPS];
    gamma = new float[ER_NUMTAPS];
    
    
    for(int i = 0; i<ER_NUMTAPS; i++){
        gains[i] = 0.0f;
        beta[i]  = 0.0f;
        upsilon[i] = 0.0f;
    }
    
    averageGamma = 0.0f;
    
    Vector3D points [NUM_FIRSTORDER][NUM_MONTECARLO];
    
    
    float gxTempValues[numberDelays*NUM_MONTECARLO];
    float gx1x2gxTempValues[numberDelays*NUM_MONTECARLO];
    float gx1x2gxhTempValues[numberDelays*NUM_MONTECARLO];
    
    // clean tempvalues
    for (int t = 0; t<numberDelays*NUM_MONTECARLO; t++){
        gx1x2gxTempValues[t] = 0.0f;
        gxTempValues[t] = 0.0f;
        gx1x2gxhTempValues[t] = 0.f;
    }
    
    
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
        monteCarloBeta(points[i], L, S, surfaces[i].normal, NUM_MONTECARLO, &beta[i], surfaces[i].getArea());
    }
    
    
    //Calculate Upsilon
    if (gainOn){
    for (int i = 0; i<numberDelays; i++){
        float hInt = 0.0f;
        for(int j = 0; j<NUM_MONTECARLO; j++){
            //Integrate h over the surface patch
            float airAbsorption = 0.0f;

                airAbsorption = expf(-0.0f * points[i][j].distance(L));

            
            hInt += (surfaces[i].getArea() * pointCollectionFunction(points[i][j], L, surfaces[i].normal, 1.0f, 0.0f)  * airAbsorption)/(float)(NUM_MONTECARLO) * ENERGYDETECTEDBYLISTENER;
            
            //            hInt += (surfaces[i].getArea() * pointCollectionFunction(points[i][j], L, surfaces[i].normal, 1.0f, 0.0f))/(float)(NUM_MONTECARLO) * ENERGYDETECTEDBYLISTENER;
            //                   printf("HInt %f ", hInt);
        }
        averageUpsilon += (( 1.0f/ ( totalSurfaceArea)) * 1.0/M_PI * hInt);
        
        
    }
    }
    //
    //    printf("Average Upsilon: %f \n", averageUpsilon);
    
    
    printf("Total pure energy is: %f \n", 4 * M_PI * dmin * dmin);
    
    float totalEnergy = 0.0f;
    
    int idx = numberDelays;
    //Calculating Gamma and Beta
    for (int patch = 0; patch < numberDelays; patch++){ //to
        for (int patchBefore = 0 ; patchBefore < numberDelays; patchBefore++){ //from
            
            double patchBeta = 0.f;
            double patchGamma = 0.f;
            if (patch != patchBefore){
                double gainInt = 0.0f;
                float MonteCarloCoeff = (surfaces[patch].getArea() * surfaces[patchBefore].getArea()) /(NUM_MONTECARLO * NUM_MONTECARLO);
                //                printf("From patch : %d to patch %d, monte carlo coeff is %f \n", patchBefore, patch, MonteCarloCoeff);
         
                
                for (int pointsPatch = 0 ; pointsPatch < NUM_MONTECARLO; pointsPatch ++){ //to
                    for (int pointsPatchBefore = 0 ;pointsPatchBefore < NUM_MONTECARLO; pointsPatchBefore ++){ //from
                        
                        
                        
                        Vector3D xS = S.subtract(points[patchBefore][pointsPatchBefore]);
                        Vector3D xPatchBefore = points[patchBefore][pointsPatchBefore];
                        Vector3D xSnormalized = Lambda(xPatchBefore, S);
//                        Vector3D xS = S.subtract(xPatchBefore);
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
                        
//                        if (gainOn){
//                            patchGamma += gx1x2 * x1x2AirAbsorption *  (1.f-absorptionCoefficient)/M_PI * sxAirAbsorption * gx * dmin * dmin  * MonteCarloCoeff*  (1.f-absorptionCoefficient); //got twice absorption loss here, since upsilon has the division by pi
//                        }
                        
                        float h = pointCollectionFunction(points[patch][pointsPatch], L, surfaces[patch].normal, 1.0f, ALPHA);
                        float x2LAirAbsorption = expf(-ALPHA * xPatch.distance(L));
                        
                        patchBeta += x2LAirAbsorption * h *  (1.f-absorptionCoefficient)/M_PI * x1x2AirAbsorption * gx1x2 *  (1.f-absorptionCoefficient)/M_PI* sxAirAbsorption * gx * dmin * dmin * MonteCarloCoeff * ENERGYDETECTEDBYLISTENER;
                        
                        //                        printf("delta gamma %0.12f, delta beta %0.12f  \n", gx1x2 * x1x2AirAbsorption * ENERGYLOSS/M_PI * sxAirAbsorption * gx * dmin * dmin  * MonteCarloCoeff, x2LAirAbsorption * h * ENERGYLOSS/M_PI * x1x2AirAbsorption * gx1x2 * ENERGYLOSS/M_PI* sxAirAbsorption * gx * dmin * dmin * MonteCarloCoeff);
                     
                        float distanceSx = S.distance(points[patch][pointsPatch]);
                        float distanceSquaredSx = powf(distanceSx, 2);
//                        float gainReductionSx = powf(get_gain(rt60, distanceSx / 343.f * 44100.f), 2);
                        
                        xS = S.subtract(points[patch][pointsPatch]);
                        xS = xS.normalize();
                        float dotProductSx1 = xS.dotProduct(surfaces[patch].normal);
                        
                        
                        
                        float dotProductx1x2 = points[patch][pointsPatch].subtract(points[patchBefore][pointsPatchBefore]).normalize().dotProduct(surfaces[patchBefore].normal);
                        
                        float dotProductx2x1 = points[patchBefore][pointsPatchBefore].subtract(points[patch][pointsPatch]).normalize().dotProduct(surfaces[patch].normal);
                        
                        float distancexx = 0.0f;
                        
                        if (useERLength)
                            distancexx = points[patch][pointsPatch].distance(points[patchBefore][pointsPatchBefore]) + points[patchBefore][pointsPatchBefore].distance(L);
                        //this is like saying the spread amount from that point x on the wall to any location thats distance L away from that point on the wall,
                        //so its all the energy in the enclosed hemisphere (which is energy at that point x) but decayed by distance L according to the RT60 formula
                        //it is integration of dot product from x to dL, with energy Ex/pi * gainRt60, and the integration is pi, so we just have Ex*gainRT60
                        else
                            distancexx = points[patch][pointsPatch].distance(points[patchBefore][pointsPatchBefore]);
                        
                        
                        float distanceSquaredxx = powf( points[patch][pointsPatch].distance(points[patchBefore][pointsPatchBefore]), 2);
//                        float gainReductionxx = powf(get_gain(rt60, distancexx / 343.f * 44100.f), 2);
                        
                        //                        gainInt  +=   (gainReductionxx / distanceSquaredxx) * dotProductx1x2 * (1.f/(2*M_PI)) * (gainReductionSx / distanceSquaredSx) * dotProductSx1 * dmin * dmin * monteCarloNormalizationTerm;
                        
                        if (gainOn){
                            float gainReductionSx = powf(get_gain(rt60, distanceSx / 343.f ), 2.f);
                            float gainReductionxx = powf(get_gain(rt60, distancexx / 343.f ), 2.f);
                            gainInt  +=   (gainReductionxx / distanceSquaredxx) * dotProductx2x1 * dotProductx1x2 * (1.f/(M_PI)) * (gainReductionSx / distanceSquaredSx) * dotProductSx1 * dmin * dmin * MonteCarloCoeff;
                        
                        }
                        
                    }
                    
                    
                }
                
                if(gainOn){
                    totalEnergy += gainInt;
                    gains[idx] = sqrtf(gainInt);
                
//                    gamma[idx] = patchGamma;
                }
                
                beta[idx] = patchBeta;
                idx++;
                
            }
        }
        
    }
    
    averageUpsilon = sqrtf(averageUpsilon);
    averageGamma = 0.0f;
    averageBeta = 0.0f;
    
    for (int i = 0; i<ER_NUMTAPS; i++){
        averageGamma += gamma[i];
        averageBeta += beta[i];
        beta[i] = sqrtf(beta[i]);
//        printf("beta %d is %f , ", i,beta[i]);
//        gamma[i] = sqrtf(gamma[i]);
        
    }
    printf("average upsilon %f \n", averageUpsilon);
    printf("average beta %f \n", averageBeta);
//    averageGamma = sqrtf(averageGamma);
    averageBeta = sqrtf(averageBeta);
//    printf("average gamma in energy %f\n", averageGamma * averageGamma );
    printf("Total energy after 2nd order : %f \n", totalEnergy);
    
}


