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

#include <inaccel/coral>

#include "models/xf_fintech_mc_american.hpp"

#include "xf_fintech_mc_american_kernel_constants.hpp"

using namespace xf::fintech;

static const size_t PRICE_ELEMENT_SIZE = sizeof(KDataType) * UN_K1;
static const size_t MATRIX_ELEMENT_SIZE = sizeof(KDataType);
static const size_t COEFF_ELEMENT_SIZE = sizeof(KDataType) * COEF;

static const size_t PRICE_NUM_ELEMENTS = DEPTH_P;
static const size_t MATRIX_NUM_ELEMENTS = DEPTH_M;
static const size_t COEFF_NUM_ELEMENTS = TIMESTEPS - 1;

MCAmerican::MCAmerican() {
    inaccel::allocator<uint8_t> u8_allocator;
    inaccel::allocator<KDataType> kdatatype_allocator;

    m_hostOutputPricesBuffer = u8_allocator.allocate(PRICE_ELEMENT_SIZE * PRICE_NUM_ELEMENTS);
    m_hostOutputMatrixBuffer = u8_allocator.allocate(MATRIX_ELEMENT_SIZE * MATRIX_NUM_ELEMENTS);
    m_hostCoeffBuffer = u8_allocator.allocate(COEFF_ELEMENT_SIZE * COEFF_NUM_ELEMENTS);
    m_hostOutputBuffer1 = kdatatype_allocator.allocate(1);
    m_hostOutputBuffer2 = kdatatype_allocator.allocate(1);
}

MCAmerican::~MCAmerican() {
    inaccel::allocator<uint8_t> u8_allocator;
    inaccel::allocator<KDataType> kdatatype_allocator;

    u8_allocator.deallocate(m_hostOutputPricesBuffer, PRICE_ELEMENT_SIZE * PRICE_NUM_ELEMENTS);
    u8_allocator.deallocate(m_hostOutputMatrixBuffer, MATRIX_ELEMENT_SIZE * MATRIX_NUM_ELEMENTS);
    u8_allocator.deallocate(m_hostCoeffBuffer, COEFF_ELEMENT_SIZE * COEFF_NUM_ELEMENTS);

    kdatatype_allocator.deallocate((KDataType*)m_hostOutputBuffer1, 1);
    kdatatype_allocator.deallocate((KDataType*)m_hostOutputBuffer2, 1);
}

void MCAmerican::run(OptionType optionType,
                    double stockPrice,
                    double strikePrice,
                    double riskFreeRate,
                    double dividendYield,
                    double volatility,
                    double timeToMaturity,
                    double requiredTolerance,
                    double* pOptionPrice) {
    unsigned int requiredSamples;

    // The kernels take in BOTH requiredTolerance AND requiredSamples.
    // However only ONE is used during processing...
    // If requiredSamples > 0, the model will run for that number of samples
    // If requiredSamples == 0, the model will run for as long as necessary to
    // meet requiredTolerance

    // since this method only exposes requiredTolerance, we must set
    // requiredSamples = 0
    requiredSamples = 0;

    runInternal(optionType, stockPrice, strikePrice, riskFreeRate, dividendYield, volatility, timeToMaturity,
                         requiredTolerance, requiredSamples, pOptionPrice);
}

void MCAmerican::run(OptionType optionType,
                    double stockPrice,
                    double strikePrice,
                    double riskFreeRate,
                    double dividendYield,
                    double volatility,
                    double timeToMaturity,
                    unsigned int requiredSamples,
                    double* pOptionPrice) {
    double requiredTolerance;

    // The kernels take in BOTH requiredTolerance AND requiredSamples.
    // However only ONE is used during processing...
    // If requiredSamples > 0, the model will run for that number of samples
    // If requiredSamples == 0, the model will run for as long as necessary to
    // meet requiredTolerance

    // since this method only exposes requiredSamples, we must set
    // requiredTolerance = 0.0
    requiredTolerance = 0.0;

    runInternal(optionType, stockPrice, strikePrice, riskFreeRate, dividendYield, volatility, timeToMaturity,
                         requiredTolerance, requiredSamples, pOptionPrice);
}

void MCAmerican::runInternal(OptionType optionType,
                            double stockPrice,
                            double strikePrice,
                            double riskFreeRate,
                            double dividendYield,
                            double volatility,
                            double timeToMaturity,
                            double requiredTolerance,
                            unsigned int requiredSamples,
                            double* pOptionPrice) {
    unsigned int timeSteps = TIMESTEPS;

    unsigned int calibrateSamples = 4096;

    m_runStartTime = std::chrono::high_resolution_clock::now();

    // --------------------
    // Run PRESAMPLE kernel
    // --------------------
	inaccel::request preSample("com.xilinx.vitis.quantitativeFinance.monteCarlo.PreSample");
	preSample.arg((KDataType)stockPrice)
		.arg((KDataType)volatility)
		.arg((KDataType)riskFreeRate)
		.arg((KDataType)dividendYield)
		.arg((KDataType)timeToMaturity)
		.arg((KDataType)strikePrice)
		.arg(optionType)
		.arg_array(m_hostOutputPricesBuffer, m_hostOutputPricesBuffer + PRICE_ELEMENT_SIZE * PRICE_NUM_ELEMENTS)
		.arg_array(m_hostOutputMatrixBuffer, m_hostOutputMatrixBuffer + MATRIX_ELEMENT_SIZE * MATRIX_NUM_ELEMENTS)
		.arg(calibrateSamples)
		.arg(timeSteps);

    // ----------------------
    // Run CALIBRATION kernel
    // ----------------------
	inaccel::request calibration("com.xilinx.vitis.quantitativeFinance.monteCarlo.Calibration");
	calibration.arg((KDataType)timeToMaturity)
		.arg((KDataType)riskFreeRate)
		.arg((KDataType)strikePrice)
		.arg(optionType)
		.arg_array(m_hostOutputPricesBuffer, m_hostOutputPricesBuffer + PRICE_ELEMENT_SIZE * PRICE_NUM_ELEMENTS)
		.arg_array(m_hostOutputMatrixBuffer, m_hostOutputMatrixBuffer + MATRIX_ELEMENT_SIZE * MATRIX_NUM_ELEMENTS)
		.arg_array(m_hostCoeffBuffer, m_hostCoeffBuffer + COEFF_ELEMENT_SIZE * COEFF_NUM_ELEMENTS)
		.arg(calibrateSamples)
		.arg(timeSteps);

    // -------------------
    // Run PRICING kernels
    // -------------------
	inaccel::request pricing1("com.xilinx.vitis.quantitativeFinance.monteCarlo.Pricing1");
	pricing1.arg((KDataType)stockPrice)
		.arg((KDataType)volatility)
		.arg((KDataType)dividendYield)
		.arg((KDataType)riskFreeRate)
		.arg((KDataType)timeToMaturity)
		.arg((KDataType)strikePrice)
		.arg(optionType)
		.arg_array(m_hostCoeffBuffer, m_hostCoeffBuffer + COEFF_ELEMENT_SIZE * COEFF_NUM_ELEMENTS)
		.arg_array((KDataType *) m_hostOutputBuffer1, (KDataType *) m_hostOutputBuffer1 + 1)
		.arg((KDataType)requiredTolerance)
		.arg(requiredSamples)
		.arg(timeSteps);

	inaccel::submit(preSample).get();
	inaccel::submit(calibration).get();
	std::future<void> prcng1 = inaccel::submit(pricing1);

	if (this->flag) {
		try {
			inaccel::request pricing2("com.xilinx.vitis.quantitativeFinance.monteCarlo.Pricing2");
			pricing2.arg((KDataType)stockPrice)
			.arg((KDataType)volatility)
			.arg((KDataType)dividendYield)
			.arg((KDataType)riskFreeRate)
			.arg((KDataType)timeToMaturity)
			.arg((KDataType)strikePrice)
			.arg(optionType)
			.arg_array(m_hostCoeffBuffer , m_hostCoeffBuffer + COEFF_ELEMENT_SIZE * COEFF_NUM_ELEMENTS)
			.arg_array((KDataType *) m_hostOutputBuffer2, (KDataType *) m_hostOutputBuffer2 + 1)
			.arg((KDataType)requiredTolerance)
			.arg(requiredSamples)
			.arg(timeSteps);

			inaccel::submit(pricing2).get();
		} catch (std::exception e) {
			this->flag = 0;
		}
	}

	prcng1.get();

    // ----------------------------------------------------------------------------------------
    // Average the outputs from the two pricing kernels, and give the result
    // back to the caller
    // ----------------------------------------------------------------------------------------
	if (this->flag)
		*pOptionPrice = (((KDataType*)m_hostOutputBuffer1)[0] + ((KDataType*)m_hostOutputBuffer2)[0]) / 2.0;
	else
		*pOptionPrice = ((KDataType*)m_hostOutputBuffer1)[0];

    m_runEndTime = std::chrono::high_resolution_clock::now();
}

long long int MCAmerican::getLastRunTime(void) {
    long long int duration = 0;

    duration =
        (long long int)std::chrono::duration_cast<std::chrono::microseconds>(m_runEndTime - m_runStartTime).count();

    return duration;
}
