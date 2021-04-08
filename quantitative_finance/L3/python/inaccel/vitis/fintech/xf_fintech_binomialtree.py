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

BinomialTreeInputDataTypeDouble = np.dtype([('S', np.double), ('K', np.double), ('T', np.double), ('rf', np.double), ('V', np.double), ('q', np.double), ('N', np.int32), ('packed', np.int32, 3)])

class BinomialTree:
	MAX_OPTION_CALCULATIONS = 1024;

	BinomialTreeEuropeanPut = np.int32(1);
	BinomialTreeEuropeanCall = np.int32(2);
	BinomialTreeAmericanPut = np.int32(3);
	BinomialTreeAmericanCall = np.int32(4);

	def __init__(self):
		with inaccel.allocator:
			self.inputBuffer = np.ndarray(self.MAX_OPTION_CALCULATIONS, dtype = BinomialTreeInputDataTypeDouble)
			self.outputBuffer = np.ndarray(self.MAX_OPTION_CALCULATIONS, dtype = np.double)

	def run(self, optionType):
		self.m_runStartTime = int(round(time.time() * 1000000))

		numOptions = np.int32(self.inputBuffer.size)
		startIndex = np.int32(0)

		if ((numOptions % 8) != 0):
			raise RuntimeError("[XLNX] BinomialTree::run - number of options to calculate should be a multiple of 8")

		req = inaccel.request("com.xilinx.vitis.quantitativeFinance.binomialTree.engine")

		req.arg(self.inputBuffer).arg(self.outputBuffer).arg(optionType).arg(numOptions).arg(startIndex)

		inaccel.submit(req).result()

		self.m_runEndTime = int(round(time.time() * 1000000))

	def lastruntime(self):
		duration = self.m_runEndTime - self.m_runStartTime

		return duration
