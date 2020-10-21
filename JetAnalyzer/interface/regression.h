#ifndef REGRESSION_H
#define REGRESSION_H



#include <list>
#include <string>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>


class regression {
	public:
		regression(std::string fileName);
		void predict(float * inputVal, float * outputVal, int dim);

	private:
		float * weights[6]; //store addresses of weight arrays
		float * biases[6]; //store addresses of bias arrays

		float normInfo1[29][2];
		float normInfo2[4][2];

		float weight0[50][29];
		float bias0[50][1];
		float weight1[50][50];
		float bias1[50][1];
		float weight2[50][50];
		float bias2[50][1];
		float weight3[50][50];
		float bias3[50][1];
		float weight4[50][50];
		float bias4[50][1];
		float weight5[4][50];
		float bias5[4][1];

		/*const static int dim = 1;
		const static int inputDim = 29;
        const static int outputDim = 1;
        const static int hiddenLayers = 5;
        const static int layerWdith = 50;

        //std::cout << "A" << std::endl;
        float innerVal1[layerWdith][dim];
        float innerVal2[layerWdith][dim];
        float inputCopy[dim][inputDim];*/

		void extractMatrix(std::string fileName, float * matrix, int dim1);
		void matrixMult(float * a, float * b, float * c, int dim1, int dim2, int dim3);
		void matrixAdd(float * a, float * b, int dim1, int dim2);
		void leakyReLU(float * a, int dim1, int dim2);
		void roundPredictions(float * a, int dim1, int dim2);
		void copyMatrix(float * a, float * b, int dim1, int dim2);
		void normalize(float * norm, float * a, int dim1, int dim2);
		void unnormalize(float * a, int dim1, int dim2);
};

#endif // REGRESSION_H
