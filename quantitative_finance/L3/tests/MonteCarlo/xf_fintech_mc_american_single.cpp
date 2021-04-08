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

#include <stdio.h>
#include <string.h>
#include <string>

#include <chrono>
#include <vector>

#include "models/xf_fintech_mc_american.hpp"

using namespace xf::fintech;

MCAmerican mcAmerican;

double varianceMultiplier;

static OptionType optionType = Put;

static const double initialStockPrice = 36.0;
static const double initialStrikePrice = 40.0;
static const double initialRiskFreeRate = 0.06;
static const double initialDividendYield = 0.0;
static const double initialVolatility = 0.20;
static const double initialTimeToMaturity = 1.0; /* in years */
static const double initialRequiredTolerance = 0.02;

/* The following variable is used to vary our input data for each run....*/
static const double varianceFactor = 0.01;

static double stockPrice;
static double strikePrice;
static double riskFreeRate;
static double dividendYield;
static double volatility;
static double timeToMaturity;
static double requiredTolerance;

/* The following variable will hold our calculated option price... */
static double optionPrice;



int main() {
    int i;
    std::chrono::time_point<std::chrono::high_resolution_clock> start;
    std::chrono::time_point<std::chrono::high_resolution_clock> end;

    printf("\n\n\n");

    printf(
        "[XLNX] "
        "***************************************************************\n");
    printf("[XLNX] Running MC AMERICAN SINGLE ASSET...\n");
    printf(
        "[XLNX] "
        "***************************************************************\n");

    //
    // Run the model a few times...
    //
    printf(
        "[XLNX] "
        "+-----------+-------------+--------------+-----------+------------+---"
        "---------+--------------+----------------+\n");
    printf(
        "[XLNX] | Iteration | Stock Price | Strike Price | Risk Free | Div. "
        "Yield | Volatility | Option Price | Execution Time |\n");
    printf(
        "[XLNX] "
        "+-----------+-------------+--------------+-----------+------------+---"
        "---------+--------------+----------------+\n");

    for (i = 0; i < 100; i++) {
        /* We will apply some variance to our data here so we are not cacheing any
         * values... */
        double variance = (1.0 + (varianceFactor * i));

        stockPrice = initialStockPrice * variance;
        strikePrice = initialStrikePrice * variance;
        riskFreeRate = initialRiskFreeRate * variance;
        dividendYield = initialDividendYield * variance;
        volatility = initialVolatility * variance;

        timeToMaturity = initialTimeToMaturity;
        requiredTolerance = initialRequiredTolerance;

        mcAmerican.run(optionType, stockPrice, strikePrice, riskFreeRate, dividendYield, volatility, timeToMaturity, requiredTolerance, &optionPrice);

        printf(
            "[XLNX] | %9d | %11.4f | %12.4f | %9.4f | %10.4f | %10.4f | %12.4f "
            "| %11lld us |\n",
            i, stockPrice, strikePrice, riskFreeRate, dividendYield, volatility, optionPrice, mcAmerican.getLastRunTime());
    }

        printf(
            "[XLNX] "
            "+-----------+-------------+--------------+-----------+------------+---"
            "---------+--------------+----------------+\n");

    return 0;
}
