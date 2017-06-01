//
//  Cuboid.cpp
//  surfaceIntegration
//
//  Created by Natalie Agus on 31/8/16.
//  Copyright © 2016 Hans. All rights reserved.
//

#include "Cuboid.hpp"
#include <algorithm>    // std::sort
using namespace std;

void Cuboid::segmentCube(int tilesPerSide){
    
    if (tilesPerSide == 1){
        segmentedSides = sides;
        elements = 6;
    }
    //Make sure tilesPerSide can be exactly square-rooted
    else{
    try
    {
        if(sqrtf(tilesPerSide) != floorf(sqrtf(tilesPerSide))){
            throw 1;
    }
    }
    catch (int e){
        printf("Please enter some number that can be square rooted for tilesPerSide \n");
        exit(1);
    }
    

    
    segmentedSides = new Plane3D[tilesPerSide*6];
    
    //now divide each side into tiles per side
    dimensions = static_cast<int>(sqrtf(tilesPerSide));

    
    int index = 0;
    for (int k = 0; k<6; k++){
        
        //get the plane
        Plane3D plane = sides[k];
        
        Vector3D s1 = plane.S1.scalarMult((float)1/dimensions);
        Vector3D s2 = plane.S2.scalarMult((float)1/dimensions);
        
        for(int i = 0; i<dimensions; i++){
            for (int j = 0; j<dimensions; j++){
                //Get the corner
                Vector3D c = plane.corner.add(plane.S1.scalarMult((float)j/dimensions));
                c = c.add(plane.S2.scalarMult((float)i/dimensions));
                
                segmentedSides[index] = Plane3D(c,s1,s2);
                segmentedSides[index].normal = plane.normal;
                index = index + 1;
            }
        }
    }
    
    elements = 6*tilesPerSide;
    }
    
    //Printing the result
//    for (int i = 0; i< 6*tilesPerSide ; i++){
//        printf("c : {%f %f %f}, s1: {%f %f %f}, s2: {%f %f %f}, n: {%f %f %f} \n", segmentedSides[i].corner.x, segmentedSides[i].corner.y, segmentedSides[i].corner.z, segmentedSides[i].S1.x, segmentedSides[i].S1.y, segmentedSides[i].S1.z, segmentedSides[i].S2.x, segmentedSides[i].S2.y, segmentedSides[i].S2.z, segmentedSides[i].normal.x, segmentedSides[i].normal.y, segmentedSides[i].normal.z);
//    }
//    
//    for (int i = 0; i<(index); i++){
//        printf("{{%f, %f, %f}, {%f, %f, %f}, {%f, %f, %f}},", segmentedSides[i].corner.x, segmentedSides[i].corner.y, segmentedSides[i].corner.z, segmentedSides[i].S1.x , segmentedSides[i].S1.y, segmentedSides[i].S1.z, segmentedSides[i].S2.x, segmentedSides[i].S2.y, segmentedSides[i].S2.z );
//    }
    
}



void Cuboid::sliceCube(int slices){
    int toSlice = 0;
    
    while(true){
        printf("toSlice %d, slices %d \n", toSlice, slices);
        if (6*toSlice * toSlice < slices){
            toSlice += 1;
        }
        if (6*toSlice * toSlice > slices){
            toSlice -= 1;
            break;
        }
        if (6*toSlice * toSlice == slices){
            break;
        }
        
    }
    
    printf("toSlice %d \n", toSlice);
    segmentCube(toSlice*toSlice);
    
    int cubicTiles = 6*toSlice*toSlice;
    Plane3D* newSegmentedSides = new Plane3D[slices];
    
    //copy into newSegmentedSides
    for (int i = 0; i<cubicTiles; i++){
        newSegmentedSides[i] = Plane3D(segmentedSides[i].corner, segmentedSides[i].S1, segmentedSides[i].S2);
    }
    
    elements = slices;
    int remainder = slices - cubicTiles;
    int index = cubicTiles;
    
//    srand(11);
    int randNum = rand()%(cubicTiles-1);
    //divide some plane into two
    for (int i = 0; i<remainder; i++){
        Plane3D refPlane = segmentedSides[randNum];
        //divide into two along S1
        if (randNum %2 == 0){
            Vector3D newCorner = refPlane.corner.add(refPlane.S1.scalarMult(0.5f));
            Vector3D newS1 = refPlane.S1.scalarMult(0.5f);
            newSegmentedSides[randNum] = Plane3D(refPlane.corner, newS1, refPlane.S2);
            newSegmentedSides[index] = Plane3D(newCorner, newS1, refPlane.S2);
            index++;
            randNum = rand()%(cubicTiles-1);
        }
        else{
            Vector3D newCorner = refPlane.corner.add(refPlane.S2.scalarMult(0.5f));
            Vector3D newS2 = refPlane.S2.scalarMult(0.5f);
            newSegmentedSides[randNum] = Plane3D(refPlane.corner,  refPlane.S1, newS2);
            newSegmentedSides[index] = Plane3D(newCorner,refPlane.S1, newS2);
            index++;
            randNum = rand()%(cubicTiles-1);
        }
    }
    
    segmentedSides = new Plane3D[slices];
    //copy back into segmentedsides
    for (int i = 0; i<elements; i++){
        segmentedSides[i] = Plane3D(newSegmentedSides[i].corner, newSegmentedSides[i].S1, newSegmentedSides[i].S2);
    }
    
    free(newSegmentedSides);
    
//    ////
//        for (int i = 0; i<(slices); i++){
//            printf("{{%f, %f, %f}, {%f, %f, %f}, {%f, %f, %f}},", segmentedSides[i].corner.x, segmentedSides[i].corner.y, segmentedSides[i].corner.z, segmentedSides[i].S1.x , segmentedSides[i].S1.y, segmentedSides[i].S1.z, segmentedSides[i].S2.x, segmentedSides[i].S2.y, segmentedSides[i].S2.z );
//        }
    
}


float Cuboid::segmentCubeOnce(){
    //segment floor and ceiling
    segmentedSides = new Plane3D[6];
    int index = 0;

    for (int k = 0; k<6; k++){
        segmentedSides[index] = sides[k];
        index ++;
    }
    
//    for (int i = 0; i< 6; i++){
//        printf("segmentCubeOnce () %f %f %f, %f %f %f, %f %f %f \n ",segmentedSides[i].corner.x, segmentedSides[i].corner.y, segmentedSides[i].corner.z, segmentedSides[i].S1.x,  segmentedSides[i].S1.y, segmentedSides[i].S1.z, segmentedSides[i].S2.x, segmentedSides[i].S2.y, segmentedSides[i].S2.z);
//    }
//    
//    for (int i= 0; i< 6; i++){
//        printf("Area of side %d is %f \n", i, segmentedSides[i].getArea());
//    }
     return index;
    
}


void Cuboid::segmentCubeBasedOnProjectedArea(int numDelays, Vector3D S, Vector3D L){
    segmentedSides = new Plane3D[numDelays];
    elements = numDelays;
    
  //  int index = segmentCubeOnce(3);
    int index = segmentCubeOnce();
    
//    printf("Index is finally: %d \n", index);
    

//        for (int i = 0; i< 6; i++){
//            printf("segmentCubeOnce () %f %f %f, %f %f %f, %f %f %f \n ",segmentedSides[i].corner.x, segmentedSides[i].corner.y, segmentedSides[i].corner.z, segmentedSides[i].S1.x,  segmentedSides[i].S1.y, segmentedSides[i].S1.z, segmentedSides[i].S2.x, segmentedSides[i].S2.y, segmentedSides[i].S2.z);
//        }
   // int delaysCeilingFloor = (2*3*3);
    int delaysCeilingFloor = 0;
    
    while(index<numDelays){
        
        //get the side with the largest projected area
        int maxIndex = 0;
        float maxArea = 0;
        
        for (int i = delaysCeilingFloor; i<index; i++){

            //add is better
            float comparisonValue =   fabs(ProjectedArea_rectangleSubDiv(segmentedSides[i], S, 5)) + fabs(ProjectedArea_rectangleSubDiv(segmentedSides[i], L, 5)) ;

           // printf("%d subdivided ? \n ", segmentedSides[i].subdivided);
            if (comparisonValue > maxArea && segmentedSides[i].subdivided < numDelays){
                maxArea = comparisonValue;
                maxIndex = i;
            }
            
//            else if (segmentedSides[i].subdivided >= 5){
//                printf("SUBDIVIDED TOO MUCH... \n");
//            }
        }
        
        //split that patch
        index = dividePlane(segmentedSides[maxIndex], index, maxIndex, S, L);
        
     //   index = dividePlaneAlongS1(segmentedSides[maxIndex], index, maxIndex);
        
    }
//        for (int i = 0; i<(index); i++){
//            printf("{{%f, %f, %f}, {%f, %f, %f}, {%f, %f, %f}},", segmentedSides[i].corner.x, segmentedSides[i].corner.y, segmentedSides[i].corner.z, segmentedSides[i].S1.x , segmentedSides[i].S1.y, segmentedSides[i].S1.z, segmentedSides[i].S2.x, segmentedSides[i].S2.y, segmentedSides[i].S2.z );
//        }


}

//
////NEED CHANGES
//float Cuboid::projectedAreaOfAPlane(Vector3D S, Vector3D L, Plane3D patch){
//    // (projectedArea(An,L) + projectedArea(An,S)) / delayLength^2
// //   float areaAnL = patch.getArea() * fabs(L.normalize().dotProduct(patch.normal)) / powf(L.distance(patch.getMidpoint()), 2);
//    //float areaAnS = patch.getArea() * fabs(S.normalize().dotProduct(patch.normal)) / powf(S.distance(patch.getMidpoint()), 2);
//    float delayLength = S.distance(patch.getMidpoint()) + L.distance(patch.getMidpoint());
//    return 1 / powf(delayLength, 2);
//  //  return patch.getArea() * fabs(S.normalize().dotProduct(patch.normal)) * (1.f/powf(S.distance(patch.getMidpoint()),2));
//}


int Cuboid::dividePlaneAlongS1(Plane3D divide, int index, int sourceIndex){
    
    Vector3D newCorner = divide.corner.add(divide.S2.scalarMult(0.5));
    Vector3D newS2 = divide.S2.scalarMult(0.5);
    segmentedSides[sourceIndex] = Plane3D(divide.corner, divide.S1, newS2);
    
    segmentedSides[index] = Plane3D(newCorner, divide.S1, newS2);
    index++;
    
    return index;
}

int Cuboid::dividePlane(Plane3D divide, int index, int sourceIndex, Vector3D S, Vector3D L ){
//    
    Vector3D omegaS1L= L.subtract(divide.corner.add(divide.S1.scalarMult(0.5))).normalize();
    Vector3D omegaS1S= S.subtract(divide.corner.add(divide.S1.scalarMult(0.5))).normalize();
    float projectedLengthS1 = divide.S1.magnitude() * omegaS1L.dotProduct(divide.normal) + divide.S1.magnitude() * omegaS1S.dotProduct(divide.normal);
    
    Vector3D omegaS2L = L.subtract(divide.corner.add(divide.S2.scalarMult(0.5))).normalize();
    Vector3D omegaS2S = S.subtract(divide.corner.add(divide.S2.scalarMult(0.5))).normalize();
    float projectedLengthS2 = divide.S2.magnitude() * omegaS2L.dotProduct(divide.normal) + divide.S2.magnitude() * omegaS2S.dotProduct(divide.normal);
//    
//    //project S1 into omega
//    float projectedAreaS1 = omega.normalize().scalarMult( divide.S1.dotProduct(omega)).magnitude();
//    //project S2 into omega
//    float projectedAreaS2 = omega.normalize().scalarMult( divide.S2.dotProduct(omega)).magnitude();
//    
//    
//    
//    
    Vector3D c2 = divide.corner.add(divide.S1.scalarMult(0.5));
    Vector3D s1 = divide.S1.scalarMult(0.5);
    Vector3D c22 = divide.corner.add(divide.S2.scalarMult(0.5));
    Vector3D s2 = divide.S2.scalarMult(0.5);

//    float projectedAreaS1 = fabs(ProjectedArea_rectangleSubDiv(Plane3D(divide.corner, s1, divide.S2), L, 2)) * fabs(ProjectedArea_rectangleSubDiv(Plane3D(divide.corner, s1, divide.S2), S, 2))  + fabs(ProjectedArea_rectangleSubDiv(Plane3D(c2, s1, divide.S2), L, 2)) * fabs(ProjectedArea_rectangleSubDiv(Plane3D(c2, s1, divide.S2), S, 2));
//
//
//    float projectedAreaS2 = fabs(ProjectedArea_rectangleSubDiv(Plane3D(divide.corner, divide.S1, s2), L, 2)) * fabs(ProjectedArea_rectangleSubDiv(Plane3D(divide.corner, divide.S1, s2), S, 2)) + fabs(ProjectedArea_rectangleSubDiv(Plane3D(c22, divide.S1, s2), L, 2)) * fabs(ProjectedArea_rectangleSubDiv(Plane3D(c22, divide.S1, s2), S, 2));
//    
//  

    if (projectedLengthS1 > projectedLengthS2){
        //divide along S1
        int subdivision = segmentedSides[sourceIndex].subdivided;
        segmentedSides[sourceIndex] = Plane3D(divide.corner, s1, divide.S2);
        segmentedSides[sourceIndex].subdivided = subdivision + 1;

        segmentedSides[index] = Plane3D(c2, s1, divide.S2);
        segmentedSides[index].subdivided = segmentedSides[sourceIndex].subdivided;
        
        index++;
    }
    else{
        //Divide along S2
        int subdivision = segmentedSides[sourceIndex].subdivided;
        segmentedSides[sourceIndex] = Plane3D(divide.corner, divide.S1, s2);
        segmentedSides[sourceIndex].subdivided = subdivision + 1;

        segmentedSides[index] = Plane3D(c22, divide.S1, s2);
        segmentedSides[index].subdivided = segmentedSides[sourceIndex].subdivided;
        index++;
        

    }
    return index;
}

//Either S1 or S2
bool Cuboid::longestDimension(Plane3D patch){
    if (patch.S1.magnitude() > patch.S2.magnitude()){
        return true;
    }
    else{
        return false;
    }
}


void Cuboid::getDelayValues(size_t *delayValues, Vector3D leftEar, Vector3D rightEar, Vector3D S, int Hz, int startIndex){
    for (int i = startIndex; i< elements + startIndex; i++){
        Vector3D p = segmentedSides[i].getMidpoint();

        float d1 = S.subtract(p).magnitude();
        float left = leftEar.subtract(p).magnitude();
        float right = rightEar.subtract(p).magnitude();
        
        if (right < left){

            delayValues[i] = static_cast<int>((d1+right)/SOUNDSPEED*Hz);
        }
        else if (left <= right){
            delayValues[i] = static_cast<int>((d1+left)/SOUNDSPEED*Hz);
        }
        
    }
}

void Cuboid::getSecondOrderDelayValues(float *delayValues, Vector3D L, Vector3D S, int Hz){
    
    
//    for (int i = 0; i< 6; i++){
//        printf(" i: %d , def %f %f %f, %f %f %f, %f %f %f \n ", i, segmentedSides[i].corner.x, segmentedSides[i].corner.y, segmentedSides[i].corner.z, segmentedSides[i].S1.x,  segmentedSides[i].S1.y, segmentedSides[i].S1.z, segmentedSides[i].S2.x, segmentedSides[i].S2.y, segmentedSides[i].S2.z);
//    }

    
    //First order
    int index = 0;
    for (int i =0; i< elements; i++){
        Vector3D p = segmentedSides[i].getMidpoint();
        
        float d1 = S.subtract(p).magnitude();
        float d2 = L.subtract(p).magnitude();
    
        delayValues[i] = ((d1+d2)/SOUNDSPEED);
//        printf("%d delay value : %zu \n", index, delayValues[i]);
        index ++;
    }
    
    //Second order
    for (int i =0; i< elements; i++){
        
        Vector3D wall1 = segmentedSides[i].getMidpoint();
        float d1 = wall1.subtract(S).magnitude();
        
        for (int j = 0; j <elements; j++){
            if (i!=j){
            Vector3D wall2 = segmentedSides[j].getMidpoint();
            
            float d2 = wall1.subtract(wall2).magnitude();
            
            float d3 = wall2.subtract(L).magnitude();
        
            delayValues[index] = ((d1+d2+d3)/SOUNDSPEED);
//            printf("2nd order %d delay value : %zu \n", index, delayValues[index]);
            index++;
            }
        }
    }
    
}



float Cuboid::projAreaSubSec(Plane3D r, Vector3D L){
    // centre of r
    Vector3D c = r.getMidpoint();
    
    // area of r
    float A = r.getArea();
    
    // surface normal of r
    Vector3D n =   r.S1.crossProduct(r.S2).normalize();//  normalize(crossProduct(r.s1, r.s2));
    
    // L-c
    Vector3D cToL = L.subtract(c);
    
    // distance from L to centre of r
    float d = cToL.magnitude();
    
    // direction from c to L
    Vector3D Omega = cToL.normalize();
    
    // projected area = Area (n . Omega) / d^2
    return (A * fabs(n.dotProduct(Omega))) / (d*d);
}




/*!
 * Estimate the area of the rectangle r projected onto the unit sphere
 * centred at L.
 *
 * @param r         the rectangle
 * @param L         point at centre of the unit sphere of projection
 * @param divisions number of times to subdivide the rectangle in each
 *                  direction (l x w). higher values are more accurate.
 */
float Cuboid::ProjectedArea_rectangleSubDiv(Plane3D r, Vector3D L, size_t divisions){
    // how many sub patches?
    size_t numPatches = divisions*divisions;
    
    // allocate memory
    Plane3D* subMesh = new Plane3D[numPatches];
    
    // subdivide r
    Vector3D s1s = r.S1.scalarMult(1.0f/(float)divisions);//  scalarDivide(r.s1,divisions);
    Vector3D s2s = r.S2.scalarMult(1.0f/(float)divisions);;
    for(size_t i=0; i<divisions; i++){
        for(size_t j=0; j<divisions; j++){
            Vector3D cij =  r.corner.add(s1s.scalarMult(i)).add(s2s.scalarMult(j));// add3(r.c,
                             //  scalarProduct(i,s1s),
                             //  scalarProduct(j,s2s));
            Plane3D rij = Plane3D(cij, s1s, s2s);
            subMesh[divisions*i + j] = rij;
        }
    }
    
    // get projected area of each subdivision
    float projArea = 0.0f;
    for(size_t i=0; i<numPatches; i++)
        projArea += projAreaSubSec(subMesh[i],L);
    
    //
    free(subMesh);
    
    return projArea;
}


