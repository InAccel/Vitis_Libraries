#
# Copyright 2019 Xilinx, Inc.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
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

class CFBlackScholes:
	NUM_ELEMENTS_PER_BUFFER_CHUNK = 16;

	def __init__(self, maxNumAssets):
		self.numAssets = maxNumAssets
		self.allocateBuffers(maxNumAssets)

	def allocateBuffers(self, numRequestedElements):
		m_numPaddedBufferElements = self.calculatePaddedNumElements(numRequestedElements)

		with inaccel.allocator:
			self.stockPrice = np.ndarray(m_numPaddedBufferElements, dtype = np.float32)
			self.strikePrice = np.ndarray(m_numPaddedBufferElements, dtype = np.float32)
			self.volatility = np.ndarray(m_numPaddedBufferElements, dtype = np.float32)
			self.riskFreeRate = np.ndarray(m_numPaddedBufferElements, dtype = np.float32)
			self.timeToMaturity = np.ndarray(m_numPaddedBufferElements, dtype = np.float32)

			self.optionPrice = np.ndarray(m_numPaddedBufferElements, dtype = np.float32)

			self.delta = np.ndarray(m_numPaddedBufferElements, dtype = np.float32)
			self.gamma = np.ndarray(m_numPaddedBufferElements, dtype = np.float32)
			self.vega = np.ndarray(m_numPaddedBufferElements, dtype = np.float32)
			self.theta = np.ndarray(m_numPaddedBufferElements, dtype = np.float32)
			self.rho = np.ndarray(m_numPaddedBufferElements, dtype = np.float32)

	def calculatePaddedNumElements(self, numRequestedElements):
		numChunks = int((numRequestedElements + (self.NUM_ELEMENTS_PER_BUFFER_CHUNK - 1)) / self.NUM_ELEMENTS_PER_BUFFER_CHUNK)

		numPaddedElements = numChunks * self.NUM_ELEMENTS_PER_BUFFER_CHUNK

		return numPaddedElements

	def run(self, optionType, numAssets):
		self.runAsync(optionType, numAssets)
		self.wait()

	def runAsync(self, optionType, numAssets):
		self.m_runStartTime = int(round(time.time() * 1000000))

		numPaddedAssets = np.int32(self.calculatePaddedNumElements(numAssets))

		req = inaccel.request("com.xilinx.vitis.quantitativeFinance.blackScholes.calculator")

		req.arg(self.stockPrice).arg(self.volatility).arg(self.riskFreeRate).arg(self.timeToMaturity).arg(self.strikePrice).arg(optionType.value).arg(numPaddedAssets).arg(self.optionPrice).arg(self.delta).arg(self.gamma).arg(self.vega).arg(self.theta).arg(self.rho)

		self.response = inaccel.submit(req)

	def wait(self):
		self.response.result()

		self.m_runEndTime = int(round(time.time() * 1000000))

	def lastruntime(self):
		duration = self.m_runEndTime - self.m_runStartTime

		return duration
