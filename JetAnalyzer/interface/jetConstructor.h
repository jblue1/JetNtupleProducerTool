#ifndef JETCONSTRUCTOR_H
#define JETCONSTRUCTOR_H

#include "classifier.h"
#include "regression.h"

class jetConstructor {
    public:
        jetConstructor();
        virtual ~jetConstructor();
        int run(float ** particles, int particleCount, float ** outputJets);

    protected:

    private:
        float dR2(float eta1, float eta2, float phi1, float phi2);
        float dPhi(float phi1, float phi2);

        regression regressionNet = regression("/home/brkronheim/Honors/algorithm/jetAlgorithm/data/regression");
		classifier classifierNet = classifier("/home/brkronheim/Honors/algorithm/jetAlgorithm/data/classification");

};

#endif // JETCONSTRUCTOR_H
