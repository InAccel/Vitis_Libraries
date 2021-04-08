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

#ifndef _XF_FINTECH_MC_AMERICAN_H_
#define _XF_FINTECH_MC_AMERICAN_H_

#include <chrono>

namespace xf {
namespace fintech {

typedef enum {
    Call = 0,
    Put = 1

} OptionType;

/**
 * @class MCAmerican
 *
 * @brief This class implements the Monte-Carlo American model.
 *
 */
class MCAmerican {
   public:
    MCAmerican();
    virtual ~MCAmerican();

   public:
    /**
     * Run a single asset until the REQUIRED TOLERANCE is met...
     *
     * @param optionType either American/European Call or Put
     * @param stockPrice the stock price
     * @param strikePrice the strike price
     * @param riskFreeRate the risk free interest rate
     * @param dividendYield the dividend yield
     * @param volatility the volatility
     * @param timeToMaturity the time to maturity
     * @param requiredTolerance the tolerance
     * @param pOptionPrice the returned option price
     *
     */
    void run(OptionType optionType,
            double stockPrice,
            double strikePrice,
            double riskFreeRate,
            double dividendYield,
            double volatility,
            double timeToMaturity,
            double requiredTolerance,
            double* pOptionPrice);

    /**
     * Run a single asset for the REQUIRED NUMBER OF SAMPLES...
     *
     * @param optionType either American/European Call or Put
     * @param stockPrice the stock price
     * @param strikePrice the strike price
     * @param riskFreeRate the risk free interest rate
     * @param dividendYield the dividend yield
     * @param volatility the volatility
     * @param timeToMaturity the time to maturity
     * @param requiredSamples the number of samples
     * @param pOptionPrice the returned option price
     *
     */
    void run(OptionType optionType,
            double stockPrice,
            double strikePrice,
            double riskFreeRate,
            double dividendYield,
            double volatility,
            double timeToMaturity,
            unsigned int requiredSamples,
            double* pOptionPrice);

   public:
    /**
     * This method returns the time the execution of the last call to run() took
     *
     * @returns Execution time in microseconds
     */
    long long int getLastRunTime(void);

   private:
    void runInternal(OptionType optionType,
                    double stockPrice,
                    double strikePrice,
                    double riskFreeRate,
                    double dividendYield,
                    double volatility,
                    double timeToMaturity,
                    double requiredTolerance,
                    unsigned int requiredSamples,
                    double* pOptionPrice);

   private:
    uint8_t* m_hostOutputPricesBuffer;
    uint8_t* m_hostOutputMatrixBuffer;
    uint8_t* m_hostCoeffBuffer;
    void* m_hostOutputBuffer1;
    void* m_hostOutputBuffer2;

    int flag = 1;

   private:
    std::chrono::time_point<std::chrono::high_resolution_clock> m_runStartTime;
    std::chrono::time_point<std::chrono::high_resolution_clock> m_runEndTime;
};

} // end namespace fintech
} // end namespace xf

#endif
