#
# Copyright 2019 Xilinx, Inc.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

import inaccel.coral as inaccel
import numpy as np
import time

class MCAmerican:
	UN_K1 = 2
	COEF = 4
	TIMESTEPS = 100
	ITERATION = 4
	SZ = 8 * np.double().nbytes

	PRICE_ELEMENT_SIZE = np.double().nbytes * UN_K1
	MATRIX_ELEMENT_SIZE = np.double().nbytes
	COEFF_ELEMENT_SIZE = np.double().nbytes * COEF

	PRICE_NUM_ELEMENTS = 1024 * TIMESTEPS * ITERATION
	MATRIX_NUM_ELEMENTS = 9 * TIMESTEPS
	COEFF_NUM_ELEMENTS = TIMESTEPS - 1

	def __init__(self):
		self.flag = True
		with inaccel.allocator:
			self.m_hostOutputPricesBuffer = np.ndarray(self.PRICE_ELEMENT_SIZE * self.PRICE_NUM_ELEMENTS, dtype = np.uint8)
			self.m_hostOutputMatrixBuffer = np.ndarray(self.MATRIX_ELEMENT_SIZE * self.MATRIX_NUM_ELEMENTS, dtype = np.uint8)
			self.m_hostCoeffBuffer = np.ndarray(self.COEFF_ELEMENT_SIZE * self.COEFF_NUM_ELEMENTS, dtype = np.uint8)
			self.m_hostOutputBuffer1 = np.ndarray(1, dtype = np.double)
			self.m_hostOutputBuffer2 = np.ndarray(1, dtype = np.double)

	def run(self, optionType, _stockPrice, _strikePrice, _riskFreeRate, _dividendYield, _volatility, _timeToMaturity, samplesOrTolerance):
		self.m_runStartTime = int(round(time.time() * 1000000))

		stockPrice = np.double(_stockPrice)
		strikePrice = np.double(_strikePrice)
		riskFreeRate = np.double(_riskFreeRate)
		dividendYield = np.double(_dividendYield)
		volatility = np.double(_volatility)
		timeToMaturity = np.double(_timeToMaturity)

		if isinstance(samplesOrTolerance, float):
			requiredTolerance = np.double(samplesOrTolerance)
			requiredSamples = np.int32(0)
		else:
			requiredTolerance = np.double(0.0)
			requiredSamples = np.int32(samplesOrTolerance)

		calibrateSamples = np.int32(4096)

		timeSteps = np.int32(self.TIMESTEPS)

		# --------------------
		# Run PRESAMPLE kernel
		# --------------------
		preSample = inaccel.request("com.xilinx.vitis.quantitativeFinance.monteCarlo.PreSample")

		preSample.arg(stockPrice).arg(volatility).arg(riskFreeRate).arg(dividendYield).arg(timeToMaturity).arg(strikePrice).arg(optionType.value).arg(self.m_hostOutputPricesBuffer).arg(self.m_hostOutputMatrixBuffer).arg(calibrateSamples).arg(timeSteps)

		# ----------------------
		# Run CALIBRATION kernel
		# ----------------------
		calibration = inaccel.request("com.xilinx.vitis.quantitativeFinance.monteCarlo.Calibration");
		calibration.arg(timeToMaturity).arg(riskFreeRate).arg(strikePrice).arg(optionType.value).arg(self.m_hostOutputPricesBuffer).arg(self.m_hostOutputMatrixBuffer).arg(self.m_hostCoeffBuffer).arg(calibrateSamples).arg(timeSteps);

		# -------------------
		# Run PRICING kernels
		# -------------------
		pricing1 = inaccel.request("com.xilinx.vitis.quantitativeFinance.monteCarlo.Pricing1");
		pricing1.arg(stockPrice).arg(volatility).arg(dividendYield).arg(riskFreeRate).arg(timeToMaturity).arg(strikePrice).arg(optionType.value).arg(self.m_hostCoeffBuffer).arg(self.m_hostOutputBuffer1).arg(requiredTolerance).arg(requiredSamples).arg(timeSteps);

		inaccel.submit(preSample).result()
		inaccel.submit(calibration).result()
		prcng1 = inaccel.submit(pricing1)

		if (self.flag):
			try:
				pricing2 = inaccel.request("com.xilinx.vitis.quantitativeFinance.monteCarlo.Pricing2")
				pricing2.arg(stockPrice).arg(volatility).arg(dividendYield).arg(riskFreeRate).arg(timeToMaturity).arg(strikePrice).arg(optionType.value).arg(self.m_hostCoeffBuffer).arg(self.m_hostOutputBuffer2).arg(requiredTolerance).arg(requiredSamples).arg(timeSteps)

				inaccel.submit(pricing2).result()
			except Exception:
				self.flag = False

		prcng1.result()

		self.m_runEndTime = int(round(time.time() * 1000000))

		# ----------------------------------------------------------------------------------------
		# Average the outputs from the two pricing kernels, and give the result
		# back to the caller
		# ----------------------------------------------------------------------------------------
		if (self.flag):
			return (self.m_hostOutputBuffer1[0] + self.m_hostOutputBuffer2[0]) / 2.0
		else:
			return self.m_hostOutputBuffer1[0]

	def lastruntime(self):
		duration = self.m_runEndTime - self.m_runStartTime

		return duration
