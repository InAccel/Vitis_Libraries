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

#include "models/xf_fintech_binomialtree.hpp"
#include "xf_fintech_binomialtree_kernel_constants.hpp"

using namespace xf::fintech;

BinomialTree::BinomialTree() {
    m_hostInputBuffer.resize(MAX_OPTION_CALCULATIONS);
    m_hostOutputBuffer.resize(MAX_OPTION_CALCULATIONS);
}

void BinomialTree::run(xf::fintech::BinomialTreeInputDataType* inputData,
                      double* outputData,
                      int optionType,
                      int numOptions) {

    if (numOptions % 8 != 0) {
		throw std::runtime_error("[XLNX] BinomialTree::run - number of options to calculate should be a multiple of 8");
    }

    int num_options = numOptions;
    int start_index = 0;
    int option_type = optionType;

    m_runStartTime = std::chrono::high_resolution_clock::now();

    // prepare the data
    for (int i = 0; i < numOptions; i++) {
        m_hostInputBuffer[i].S = inputData[i].S;
        m_hostInputBuffer[i].K = inputData[i].K;
        m_hostInputBuffer[i].T = inputData[i].T;
        m_hostInputBuffer[i].rf = inputData[i].rf;
        m_hostInputBuffer[i].V = inputData[i].V;
        m_hostInputBuffer[i].q = inputData[i].q;
        m_hostInputBuffer[i].N = inputData[i].N;
    }

	inaccel::request request("com.xilinx.vitis.quantitativeFinance.binomialTree.engine");
    // Set the arguments
	request.arg_array(m_hostInputBuffer.data(), m_hostInputBuffer.data() + numOptions)
		.arg_array(m_hostOutputBuffer.data(), m_hostOutputBuffer.data() + numOptions)
		.arg(option_type)
		.arg(num_options)
		.arg(start_index);

	inaccel::submit(request).get();

    // --------------------------------
    // Give the caller back the results
    // --------------------------------
    for (int i = 0; i < numOptions; i++) {
        outputData[i] = m_hostOutputBuffer[i];
    }

    m_runEndTime = std::chrono::high_resolution_clock::now();
}

long long int BinomialTree::getLastRunTime(void) {
    long long int duration = 0;
    duration =
        (long long int)std::chrono::duration_cast<std::chrono::microseconds>(m_runEndTime - m_runStartTime).count();
    return duration;
}
