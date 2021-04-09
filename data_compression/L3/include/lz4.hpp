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

#ifndef _XFCOMPRESSION_LZ4_HPP_
#define _XFCOMPRESSION_LZ4_HPP_

#include <cstdint>
#include <string>

/**
 * Maximum compute units supported
 */
#if (C_COMPUTE_UNIT > D_COMPUTE_UNIT)
#define MAX_COMPUTE_UNITS C_COMPUTE_UNIT
#else
#define MAX_COMPUTE_UNITS D_COMPUTE_UNIT
#endif

/**
 * Maximum host buffer used to operate per kernel invocation
 */
#define HOST_BUFFER_SIZE (2 * 1024 * 1024)

/**
 * Default block size
 */
#ifndef BLOCK_SIZE_IN_KB
#define BLOCK_SIZE_IN_KB 64
#endif
/**
 * Value below is used to associate with
 * Overlapped buffers, ideally overlapped
 * execution requires 2 resources per invocation
 */
#define OVERLAP_BUF_COUNT 2

namespace xf {
namespace compression {
/**
 *  xfLz4 class. Class containing methods for LZ4
 * compression and decompression to be executed on host side.
 */
class xfLz4 {
   public:

    /**
     * @brief This module does the sequential execution of compression
     * where all the I/O operations and kernel execution are done one
     * after another in sequential order.
     *
     * @param in input byte sequence
     * @param out output byte sequence
     * @param actual_size input size
     * @param host_buffer_size host buffer size
     */
    uint64_t compress(uint8_t* in, uint8_t* out, uint64_t actual_size, bool file_list_flag);

    /**
    * @brief Decompress.
    *
    * @param in input byte sequence
    * @param out output byte sequence
    * @param actual_size input size
    * @param original_size original size
    * @param host_buffer_size host buffer size
    */

    uint64_t decompress(uint8_t* in,
                        uint8_t* out,
                        uint64_t actual_size,
                        uint64_t original_size,
                        bool file_list_flag);

    /**
     * @brief This module does the memory mapped execution of decompression
     * where the I/O operations and kernel execution is done in sequential order
     *
     * @param in input byte sequence
     * @param out output byte sequence
     * @param actual_size input size
     * @param original_size original size
     * @param host_buffer_size host buffer size
     */

    uint64_t decompressFile(std::string& inFile_name,
                            std::string& outFile_name,
                            uint64_t actual_size,
                            bool file_list_flag);

    /**
     * @brief This module is provided to support compress API and
     * it's not recommended to use for high throughput.
     *
     * @param inFile_name input file name
     * @param outFile_name output file name
     * @param actual_size input size
     */

    uint64_t compressFile(std::string& inFile_name,
                          std::string& outFile_name,
                          uint64_t actual_size,
                          bool file_list_flag);

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
     */
    xfLz4(uint32_t req_num);
    xfLz4();

    /**
     * @brief Class destructor.
     */
    ~xfLz4();
};

} // end namespace compression
} // end namespace xf
#endif // _XFCOMPRESSION_LZ4_HPP_
