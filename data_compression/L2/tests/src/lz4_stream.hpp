/*
 * (c) Copyright 2019 Xilinx, Inc. All rights reserved.
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
 *
 */
/**
 * @file xil_lz4.hpp
 * @brief Header for LZ4 host functionality
 *
 * This file is part of Vitis Data Compression Library host code for lz4 compression.
 */

#ifndef _XFCOMPRESSION_LZ4_STREAM_HPP_
#define _XFCOMPRESSION_LZ4_STREAM_HPP_

#include <assert.h>
#include <iomanip>
#include <iostream>
#include <stdint.h>
#include <vector>
#include <math.h>
#include <time.h>
#include <string>
#include <fstream>
#include <chrono>

/**
 * Default block size
 */
#define BLOCK_SIZE_IN_KB 64

/**
 * Below are the codes as per LZ4 standard for
 * various maximum block sizes supported.
 */
#define BSIZE_STD_64KB 0x40
#define BSIZE_STD_256KB 0x50
#define BSIZE_STD_1024KB 0x60
#define BSIZE_STD_4096KB 0x70

/**
 * Maximum block sizes supported by LZ4
 */
#define MAX_BSIZE_64KB (64 * 1024)
#define MAX_BSIZE_256KB (256 * 1024)
#define MAX_BSIZE_1024KB (1024 * 1024)
#define MAX_BSIZE_4096KB (4096 * 1024)

/**
 * This value is used to set
 * uncompressed block size value.
 * 4th byte is always set to below
 * and placed as uncompressed byte
 */
#define NO_COMPRESS_BIT 128

/**
 * In case of uncompressed block
 * Values below are used to set
 * 3rd byte to following values
 * w.r.t various maximum block sizes
 * supported by standard
 */
#define BSIZE_NCOMP_64 1
#define BSIZE_NCOMP_256 4
#define BSIZE_NCOMP_1024 16
#define BSIZE_NCOMP_4096 64

/**
 * @brief Validate the compressed file.
 *
 * @param inFile_name input file name
 * @param outFile_name output file name
 */
int validate(std::string& inFile_name, std::string& outFile_name);

static uint64_t getFileSize(std::ifstream& file) {
    file.seekg(0, file.end);
    uint64_t file_size = file.tellg();
    file.seekg(0, file.beg);
    return file_size;
}

/**
 *  xfLz4Streaming class. Class containing methods for LZ4
 * compression and decompression to be executed on host side.
 */
class xfLz4Streaming {
   public:
    /**
     * @brief Compress streaming.
     *
     * @param in input byte sequence
     * @param out output byte sequence
     * @param actual_size input size
     * @param host_buffer_size host buffer size
     */
    uint64_t compressStream(uint8_t* in, uint8_t* out, uint64_t actual_size);

    /**
     * @brief Compress the input file.
     *
     * @param inFile_name input file name
     * @param outFile_name output file name
     * @param actual_size input size
     */
    uint64_t compressFile(std::string& inFile_name, std::string& outFile_name, uint64_t actual_size);

    /**
     * @brief Decompress the input file.
     *
     * @param inFile_name input file name
     * @param outFile_name output file name
     * @param actual_size input size
     */
    uint64_t decompressFile(std::string& inFile_name, std::string& outFile_name, uint64_t actual_size);

    /**
     * @brief Decompress streaming.
     *
     * @param in input byte sequence
     * @param out output byte sequence
     * @param actual_size input size
     * @param original_size original size
     * @param host_buffer_size host buffer size
     */
    uint64_t decompressStream(uint8_t* in, uint8_t* out, uint64_t actual_size, uint64_t original_size);

    /**
     * Block Size
     */
    uint32_t m_BlockSizeInKb;

    /**
     * Switch between FPGA/Standard flows
     */
    bool m_switch_flow;

    uint32_t m_req_num;
    /**
     * @brief Class constructor
     *
     * @param binaryFile file to be read
     * @param flow 0 for compression, 1 for decompression
     */
    xfLz4Streaming(uint32_t req_num);

    /**
     * @brief Class destructor.
     */
    ~xfLz4Streaming();
};

#endif // _XFCOMPRESSION_LZ4_STREAM_HPP_
