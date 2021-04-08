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

#ifndef _XF_FINTECH_CF_BLACK_SCHOLES_H_
#define _XF_FINTECH_CF_BLACK_SCHOLES_H_

#include <chrono>

#include <inaccel/coral>

namespace xf {
namespace fintech {

typedef enum {
    Call = 0,
    Put = 1

} OptionType;

/**
 * @class CFBlackScholes
 *
 * @brief This class implements the Closed Form Black Scholes model.
 *
 * @details The parameter passed to the constructor controls the size of the
 * underlying buffers that will be allocated.
 * This prameter therefore controls the maximum number of assets that can be
 * processed per call to run()
 *
 * It is intended that the user will populate the input buffers with appropriate
 * asset data prior to calling run()
 * When run completes, the calculated output data will be available in the
 * relevant output buffers.
 */
class CFBlackScholes {
   public:
    CFBlackScholes(unsigned int maxAssetsPerRun);
    virtual ~CFBlackScholes();

   public:
    /**
     * @param KDataType This is the data type that the underlying HW kernel has
     * been built with.
     *
     */
    typedef float KDataType;

   public: // INPUT BUFFERS
    KDataType* stockPrice;
    KDataType* strikePrice;
    KDataType* volatility;
    KDataType* riskFreeRate;
    KDataType* timeToMaturity;

   public: // OUTPUT BUFFERS
    KDataType* optionPrice;
    KDataType* delta;
    KDataType* gamma;
    KDataType* vega;
    KDataType* theta;
    KDataType* rho;

   public:
    /**
     * This method is used to begin processing the asset data that is in the input
     * buffers.
     * If this function returns successfully, calculated results are available in
     * the output buffers.
     *
     * @param optionType The option type of ALL the assets data
     * @param numAssets The number of assets to process.
     */
    void run(OptionType optionType, unsigned int numAssets);

	/**
     * This method is used to begin processing the asset data that is in the input
     * buffers.
     * This functions is asynchonous so a call to wait must be made to wait the
     * calculation of the results
     *
     * @param optionType The option type of ALL the assets data
     */
    void runAsync(OptionType optionType, unsigned int numAssets);

	/**
	 * This method is used to wait the calculation of the results
	 * If this function returns successfully, calculated results are available in
	 * the output buffers.
	 */
	void wait();

   public:
    /**
     * This method returns the time the execution of the last call to run() took
     *
     * @returns Execution time in microseconds
     */
    long long int getLastRunTime(void); // in microseconds

   protected:
    void allocateBuffers(unsigned int numRequestedElements);
    void deallocateBuffers(void);

   protected:
    unsigned int calculatePaddedNumElements(unsigned int numRequestedElements);

   protected:
    unsigned int m_numPaddedBufferElements;

   private:
    static const unsigned int KERNEL_PARAMETER_BITWIDTH = 512;
    static const unsigned int NUM_ELEMENTS_PER_BUFFER_CHUNK;

	std::future<void> response;

   protected:
    std::chrono::time_point<std::chrono::high_resolution_clock> m_runStartTime;
    std::chrono::time_point<std::chrono::high_resolution_clock> m_runEndTime;
};

} // end namespace fintech
} // end namespace xf

#endif
