//
//  HRTFFilter.c
//  Head Related Transfer Function, first order model
//
//  Created by Hans on 22/5/17.
//  Everyone may use this code freely without restrictions.
//

#include <math.h>
#include "HRTFFilter.h"




#ifdef __cplusplus
extern "C" {
#endif
    
   
    
/*
 * process one sample of input through the filter
 */
float HRTF_process(HRTFFilter* This, float sample){
    float out = sample*This->b0 + This->zb*This->b1 - This->za*This->a1;
    This->zb = sample;
    This->za = out;
    return out;
}






/*
 * Sets the angle of incidence for a first order
 * Head Related Transfer Function filter model.
 *
 * @param theta - angle from the source to the ear with 0 pointing in
 *                  the direction out from the ear
 * @param fs    - sample rate
 */
void HRTF_setAngle(HRTFFilter* This, float theta, float fs){
    
    float theta0 = 150.0f ;
    float alpha_min = 0.1f ;
    float c = 334.0f; // speed of sound
    float a = 0.08f; // radius of head
    float w0 = c/a;
    float alpha = (1.0f + alpha_min/2.0f) +
    (1.0f - alpha_min/2.0f) * cos((theta / theta0) * M_PI);
    //    printf("Angle: %f alpha %f \n", theta, alpha);
    
    This->b0 =  (alpha+w0/fs)/(1.0f+w0/fs);
    This->b1 = (-alpha+w0/fs)/(1.0f+w0/fs);
    This->a1 =   -(1.0f-w0/fs)/(1.0f+w0/fs);
    
    This->za = This->zb = 0.0f;
}


//Input: 0 to 360 degrees, output 0 to 180 or -180 to 0
float convertAzimuthToLeftEar(float azimuth){
    
    //change azimuth to be between 0 to 180 and 0 to -180
    if (azimuth >= 360.f){
        azimuth -= 360.f;
    }
    
    if (azimuth > 180.f){
        azimuth = azimuth - 360.f;
    }
    
    if(azimuth > 90.f && azimuth <= 180.f){
        float leftEarAngle = azimuth + 90.f - 360.f;
        return leftEarAngle;
    }
    
    return azimuth + 90.f;
    
    
}


//Input: 0 to 360 degrees, output 0 to 180 or -180 to 0
float convertAzimuthToRightEar(float azimuth){
    
    //change azimuth to be between 0 to 180 and 0 to -180
    if (azimuth >= 360.f){
        azimuth -= 360.f;
    }
    
    
    if (azimuth > 180.f){
        azimuth = azimuth - 360.f;
    }
    
    if(azimuth >= -180.f && azimuth < -90.f){
        float rightEarAngle = azimuth - 90.f + 360.f;
        return rightEarAngle;
    }
    
    return azimuth - 90.f;
    
}

//Using Arctan x/y here to find how much azimuth angle is object from listener
//Returns: 0 to 360 degrees
float angleAzimuth(float xL, float yL, float xObject, float yObject){
    float angle = 0.0f;
    
    float xDistance = xObject - xL;
    float yDistance = yObject - yL;
    
    //Special cases
    if (fabs(yDistance - 0.00001) <= 0.0001f){
        // zero y, can be at 90 or -90
        if (xDistance >= 0){
            angle = 90.f;
        }
        else{
            angle = -90.f;
        }
    }
    
    if (fabs(xDistance - 0.00001) <= 0.0001f){
        //zero x, can be at 0 or 180;
        if (yDistance >= 0){
            angle = 0.0f;
        }
        else{
            angle = 180.f;
        }
    }
    
    //Normal case
    //x over y because we want azimuth
    else{
        angle = atan2f(xDistance, yDistance) * 180.0f / M_PI;
    }
    
    //shift between 0 to 360
    if (angle < 0.f){
        angle += 360.f;
    }
    
    return angle;
}

    
#ifdef __cplusplus
}
#endif
