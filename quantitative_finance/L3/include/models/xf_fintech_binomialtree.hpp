/*
 * Copyright 2019 Xilinx, Inc.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef _XF_FINTECH_BINOMIAL_TREE_H_
#define _XF_FINTECH_BINOMIAL_TREE_H_

#include <chrono>
#include <vector>

#include <inaccel/coral>

namespace xf {
namespace fintech {

const static int BinomialTreeEuropeanPut = 1;
const static int BinomialTreeEuropeanCall = 2;
const static int BinomialTreeAmericanPut = 3;
const static int BinomialTreeAmericanCall = 4;

struct BinomialTreeInputDataType {
    double S;
    double K;
    double T;
    double rf;
    double V;
    double q;
    int N;
    int packed[3]; // Bitwidth of (packed) data on axi master must be power of 2
                   // (set to support float)
};

/**
 * @class BinomialTree
 *
 * @brief This class implements the Binomial Tree Model.
 *
 * It is intended that the user will populate the inputData structure with
 * appropriate asset data prior to calling run() method. When the run completes
 * the calculated output data (one or more options) will be available to the user.
 */

class BinomialTree {
   public:
    BinomialTree();
    virtual ~BinomialTree() { }

    /**
     * Calculate one or more options based on input data and option type
     *
     * @param inputData structure to be populated with the asset data
     * @param outputData one or more calculated option values returned
     * @param optionType option type is American/European Call or Put
     * @param numOptions number of options to be calculate
     */
    void run(xf::fintech::BinomialTreeInputDataType* inputData,
            double* outputData,
            int optionType,
            int numOptions);

    /**
     * This method returns the time the execution of the last call to run() took.
     */
    long long int getLastRunTime(void);

   private:
    static const int MAX_OPTION_CALCULATIONS = 1024;

    inaccel::vector<xf::fintech::BinomialTreeInputDataType> m_hostInputBuffer;
    inaccel::vector<double> m_hostOutputBuffer;

   private:
    std::chrono::time_point<std::chrono::high_resolution_clock> m_runStartTime;
    std::chrono::time_point<std::chrono::high_resolution_clock> m_runEndTime;
};

} // end namespace fintech
} // end namespace xf

#endif /* _XF_FINTECH_BINOMIAL_TREE_H_ */
