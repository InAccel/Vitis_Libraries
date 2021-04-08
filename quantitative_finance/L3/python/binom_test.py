#!/usr/bin/env python3

# Ensure environmental variables i.e. paths are set to used the modules
from inaccel.vitis.fintech import BinomialTree, OptionType

# State test financial model
print("\nThe Binomial financial model\n========================================\n")

bto = BinomialTree()

for i in range(0, bto.inputBuffer.size):
	bto.inputBuffer[i] = (110.0, 100.0 + i, 1.0, 0.05, 0.2, 0, 1024, (0, 0, 0))

bto.run(BinomialTree.BinomialTreeEuropeanPut)

print("Output Result: num ", bto.outputBuffer.size)
for i in range(bto.outputBuffer.size):
	print(bto.outputBuffer[i])
