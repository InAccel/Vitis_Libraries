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

#include "models/xf_fintech_cf_black_scholes.hpp"

using namespace xf::fintech;

// For efficient memory bandwidth utilization, the HW kernel operates on data
// elements
// that are KERNEL_PARAMETER_WIDTH (512 at time of writing) bits wide.
//
// This is used for both input and output buffers.
//
// This means that each call to the kernel needs to pass pointers to buffers
// that at least 512 bits wide.
// i.e. we need our buffers to be allocated in chunks of KERNEL_PARAMETER_WIDTH
// / sizeof(float) elements.

const unsigned int CFBlackScholes::NUM_ELEMENTS_PER_BUFFER_CHUNK =
    CFBlackScholes::KERNEL_PARAMETER_BITWIDTH / sizeof(CFBlackScholes::KDataType);

CFBlackScholes::CFBlackScholes(unsigned int maxNumAssets) {
    this->allocateBuffers(maxNumAssets);
}

CFBlackScholes::~CFBlackScholes() {
    this->deallocateBuffers();
}

void CFBlackScholes::allocateBuffers(unsigned int numRequestedElements) {
    inaccel::allocator<KDataType> allocator;

    m_numPaddedBufferElements = calculatePaddedNumElements(numRequestedElements);

    this->stockPrice = allocator.allocate(m_numPaddedBufferElements);
    this->strikePrice = allocator.allocate(m_numPaddedBufferElements);
    this->volatility = allocator.allocate(m_numPaddedBufferElements);
    this->riskFreeRate = allocator.allocate(m_numPaddedBufferElements);
    this->timeToMaturity = allocator.allocate(m_numPaddedBufferElements);

    this->optionPrice = allocator.allocate(m_numPaddedBufferElements);

    this->delta = allocator.allocate(m_numPaddedBufferElements);
    this->gamma = allocator.allocate(m_numPaddedBufferElements);
    this->vega = allocator.allocate(m_numPaddedBufferElements);
    this->theta = allocator.allocate(m_numPaddedBufferElements);
    this->rho = allocator.allocate(m_numPaddedBufferElements);
}

void CFBlackScholes::deallocateBuffers(void) {
    inaccel::allocator<KDataType> allocator;

    if (this->stockPrice != nullptr) {
        allocator.deallocate(this->stockPrice, m_numPaddedBufferElements);
        this->stockPrice = nullptr;
    }

    if (this->strikePrice != nullptr) {
        allocator.deallocate(this->strikePrice, m_numPaddedBufferElements);
        this->strikePrice = nullptr;
    }

    if (this->volatility != nullptr) {
        allocator.deallocate(this->volatility, m_numPaddedBufferElements);
        this->volatility = nullptr;
    }

    if (this->riskFreeRate != nullptr) {
        allocator.deallocate(this->riskFreeRate, m_numPaddedBufferElements);
        this->riskFreeRate = nullptr;
    }

    if (this->timeToMaturity != nullptr) {
        allocator.deallocate(this->timeToMaturity, m_numPaddedBufferElements);
        this->timeToMaturity = nullptr;
    }

    if (this->optionPrice != nullptr) {
        allocator.deallocate(this->optionPrice, m_numPaddedBufferElements);
        this->optionPrice = nullptr;
    }

    if (this->delta != nullptr) {
        allocator.deallocate(this->delta, m_numPaddedBufferElements);
        this->delta = nullptr;
    }

    if (this->gamma != nullptr) {
        allocator.deallocate(this->gamma, m_numPaddedBufferElements);
        this->gamma = nullptr;
    }

    if (this->vega != nullptr) {
        allocator.deallocate(this->vega, m_numPaddedBufferElements);
        this->vega = nullptr;
    }

    if (this->theta != nullptr) {
        allocator.deallocate(this->theta, m_numPaddedBufferElements);
        this->theta = nullptr;
    }

    if (this->rho != nullptr) {
        allocator.deallocate(this->rho, m_numPaddedBufferElements);
        this->rho = nullptr;
    }

    m_numPaddedBufferElements = 0;
}

unsigned int CFBlackScholes::calculatePaddedNumElements(unsigned int numRequestedElements) {
    unsigned int numChunks;
    unsigned int numPaddedElements;

    // due to the way the HW processes data, the number of elements in a buffer
    // needs to be multiples of NUM_ELEMENTS_PER_BUFFER_CHUNK.
    // so we need to round up the amount to the next nearest whole number of
    // chunks

    numChunks = numRequestedElements + (NUM_ELEMENTS_PER_BUFFER_CHUNK - 1) / NUM_ELEMENTS_PER_BUFFER_CHUNK;

    numPaddedElements = numChunks * NUM_ELEMENTS_PER_BUFFER_CHUNK;

    return numPaddedElements;
}

void CFBlackScholes::run(OptionType optionType, unsigned int numAssets) {
	this->runAsync(optionType, numAssets);
	this->wait();
}

void CFBlackScholes::runAsync(OptionType optionType, unsigned int numAssets) {
	unsigned int optionFlag;

	unsigned int numPaddedAssets;

    m_runStartTime = std::chrono::high_resolution_clock::now();

    if (optionType == OptionType::Call) {
        optionFlag = 1;
    } else {
        optionFlag = 0;
    }

    numPaddedAssets = calculatePaddedNumElements(numAssets);

	inaccel::request req("com.xilinx.vitis.quantitativeFinance.blackScholes.calculator");

	req.arg_array(stockPrice, stockPrice + m_numPaddedBufferElements)
		.arg_array(volatility, volatility + m_numPaddedBufferElements)
		.arg_array(riskFreeRate, riskFreeRate + m_numPaddedBufferElements)
		.arg_array(timeToMaturity, timeToMaturity + m_numPaddedBufferElements)
		.arg_array(strikePrice, strikePrice + m_numPaddedBufferElements)
		.arg(optionFlag)
		.arg(numPaddedAssets)
		.arg_array(optionPrice, optionPrice + m_numPaddedBufferElements)
		.arg_array(delta, delta + m_numPaddedBufferElements)
		.arg_array(gamma, gamma + m_numPaddedBufferElements)
		.arg_array(vega, vega + m_numPaddedBufferElements)
		.arg_array(theta, theta + m_numPaddedBufferElements)
		.arg_array(rho, rho + m_numPaddedBufferElements);

	response = inaccel::submit(req);
}

void CFBlackScholes::wait() {
	response.get();

	m_runEndTime = std::chrono::high_resolution_clock::now();
}

long long int CFBlackScholes::getLastRunTime(void) {
    long long int duration = 0;

    duration =
        (long long int)std::chrono::duration_cast<std::chrono::microseconds>(m_runEndTime - m_runStartTime).count();

    return duration;
}
