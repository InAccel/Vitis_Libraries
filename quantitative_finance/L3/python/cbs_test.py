#!/usr/bin/env python3

# Ensure environmental variables i.e. paths are set to the named the modules
from inaccel.vitis.fintech import CFBlackScholes, OptionType

# State test financial model
print("\nThe CFBlack Scholes financial model\n==================================================\n")

# Declaring Variables
# Example financial data to test the module as used in the C++ example script
numAssets = 100  # reduced from 100000 to 100 for clarity of script output - tested at 100000 samples

# Selecting and loading into FPGA on chosen card the financial model to be used
CFBlackScholes = CFBlackScholes(numAssets)   # warning the lower levels to accomodate at least this figure

# Inputs
CFBlackScholes.stockPrice.fill(100.0)
CFBlackScholes.strikePrice.fill(100.0)
CFBlackScholes.volatility.fill(0.1)
CFBlackScholes.riskFreeRate.fill(0.025)
CFBlackScholes.timeToMaturity.fill(1.0)

#Feed in the data and request the result using tolerance method
print("\nRunning...")
CFBlackScholes.run(OptionType.Put, numAssets)
print("Done")
runtime = CFBlackScholes.lastruntime()

#Format output to match the example in C++, simply to aid comparison of results
print("+-------+-----------+----------------+--------------+---------------+---------------+---------------+")
print("| Index | Price     |     Delta      |     Gamma    |     Vega      |     Theta     |     Rho       |")
print("+-------+-----------+----------------+--------------+---------------+---------------+---------------+")
for loop in range(0, numAssets) :
    print(loop,"\t%9.5f"%CFBlackScholes.optionPrice[loop],"\t%9.5f"%CFBlackScholes.delta[loop],"\t%9.5f"%CFBlackScholes.gamma[loop],"\t%9.5f"%CFBlackScholes.vega[loop],"\t%9.5f"%CFBlackScholes.theta[loop],"\t%9.5f"%CFBlackScholes.rho[loop])



print("\nThis run took", str(runtime), "microseconds")
