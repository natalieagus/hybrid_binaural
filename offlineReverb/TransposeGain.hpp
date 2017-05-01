//
//  TransposeGain.h
//  MatrixGainTransposition
//
//  Created by Hans on 15/3/17.
//  Copyright © 2017 Hans. All rights reserved.
//



#include <stdio.h>

/*
 * returns (energy with upsilon at input) / (energy with upsilon at output)
 *
 *      Norm[M.U.ones] / Norm[U.M.ones]
 *
 * to get the gain coefficient to compensate for this, take the square root
 * of the output.
 *
 * @param upsilon - an array representing the output tap coefficients
 *                  for an imagingary FDN
 * @length        - the length of upsilon
 */
float transposeEnergyRatio(float* upsilon, size_t matrixSize);


