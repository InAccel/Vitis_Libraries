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
#include <inaccel/coral>
#include "lz4_stream.hpp"
#include "xxhash.h"

#define BLOCK_SIZE 64
#define KB 1024
#define MAGIC_HEADER_SIZE 4
#define MAGIC_BYTE_1 4
#define MAGIC_BYTE_2 34
#define MAGIC_BYTE_3 77
#define MAGIC_BYTE_4 24
#define FLG_BYTE 104

int validate(std::string& inFile_name, std::string& outFile_name) {
    std::string command = "cmp " + inFile_name + " " + outFile_name;
    int ret = system(command.c_str());
    return ret;
}

// Constructor
xfLz4Streaming::xfLz4Streaming(uint32_t req_num) { m_req_num = req_num; }

// Destructor
xfLz4Streaming::~xfLz4Streaming() {}

uint64_t xfLz4Streaming::compressFile(std::string& inFile_name, std::string& outFile_name, uint64_t input_size) {
    if (m_switch_flow == 0) { // Xilinx FPGA compression flow
        std::ifstream inFile(inFile_name.c_str(), std::ifstream::binary);
        std::ofstream outFile(outFile_name.c_str(), std::ofstream::binary);

        if (!inFile) {
            std::cout << "Unable to open file";
            exit(1);
        }

        std::vector<uint8_t> in(input_size);
        std::vector<uint8_t> out(input_size);

        inFile.read((char*)in.data(), input_size);
        // LZ4 header
        outFile.put(MAGIC_BYTE_1);
        outFile.put(MAGIC_BYTE_2);
        outFile.put(MAGIC_BYTE_3);
        outFile.put(MAGIC_BYTE_4);

        // FLG & BD bytes
        // --no-frame-crc flow
        // --content-size
        outFile.put(FLG_BYTE);

        // Default value 64K
        uint8_t block_size_header = 0;
        switch (m_BlockSizeInKb) {
            case 64:
                outFile.put(BSIZE_STD_64KB);
                block_size_header = BSIZE_STD_64KB;
                break;
            case 256:
                outFile.put(BSIZE_STD_256KB);
                block_size_header = BSIZE_STD_256KB;
                break;
            case 1024:
                outFile.put(BSIZE_STD_1024KB);
                block_size_header = BSIZE_STD_1024KB;
                break;
            case 4096:
                outFile.put(BSIZE_STD_4096KB);
                block_size_header = BSIZE_STD_4096KB;
                break;
            default:
                std::cout << "Invalid Block Size" << std::endl;
                break;
        }

        // uint32_t host_buffer_size = (m_BlockSizeInKb * 1024) * 32;
        // if ((m_BlockSizeInKb * 1024) > input_size) host_buffer_size = m_BlockSizeInKb * 1024;

        uint8_t temp_buff[10] = {FLG_BYTE,         block_size_header, input_size,       input_size >> 8,
                                 input_size >> 16, input_size >> 24,  input_size >> 32, input_size >> 40,
                                 input_size >> 48, input_size >> 56};

        // xxhash is used to calculate hash value
        uint32_t xxh = XXH32(temp_buff, 10, 0);
        uint64_t enbytes;
        outFile.write((char*)&temp_buff[2], 8);

        // Header CRC
        outFile.put((uint8_t)(xxh >> 8));
        // LZ4 streaming compression
        enbytes = compressStream(in.data(), out.data(), input_size);
        // Writing compressed data
        outFile.write((char*)out.data(), enbytes);

        outFile.put(0);
        outFile.put(0);
        outFile.put(0);
        outFile.put(0);

        // Close file
        inFile.close();
        outFile.close();
        return enbytes;
    } else { // Standard LZ4 flow
        std::string command = "../../../common/lz4/lz4 --content-size -f -q " + inFile_name;
        system(command.c_str());
        std::string output = inFile_name + ".lz4";
        std::string rout = inFile_name + ".std.lz4";
        std::string rename = "mv " + output + " " + rout;
        system(rename.c_str());
        return 0;
    }
}

// Note: Various block sizes supported by LZ4 standard are not applicable to
// this function. It just supports Block Size 64KB
uint64_t xfLz4Streaming::compressStream(uint8_t* in, uint8_t* out, uint64_t input_size) {
    uint32_t host_buffer_size = m_BlockSizeInKb * 1024;
    uint32_t total_block_count = (input_size - 1) / host_buffer_size + 1;

    // Index calculation
    std::vector<inaccel::vector<uint8_t>> blocks_in(total_block_count);
    std::vector<inaccel::vector<uint8_t>> blocks_out(total_block_count);
    std::vector<inaccel::vector<uint32_t>> compressed_block_sizes_out(total_block_count);

    std::chrono::duration<double, std::nano> kernel_time_ns_1(0);

    // copy input to input buffer
    for (uint64_t i = 0; i < total_block_count; i++) {
        uint32_t c_input_size = host_buffer_size;
        if (i == total_block_count - 1) c_input_size = input_size - (host_buffer_size * i);
        blocks_in[i].resize(c_input_size);
        blocks_out[i].resize(c_input_size);
        compressed_block_sizes_out[i].resize(1);

        std::memcpy(blocks_in[i].data(), in + (host_buffer_size * i), c_input_size);
    }

    std::vector<std::future<void>> responses(m_req_num);

    // send requests to the Coral manager, in blocks of m_req_num async requests
    for (uint64_t i = 0; i < total_block_count; i+=m_req_num) {
        auto kernel_start = std::chrono::high_resolution_clock::now();

        for (unsigned j = 0; j < m_req_num; j++) {
            if(i+j < total_block_count) {
                uint32_t c_input_size = host_buffer_size;
                if (i+j == total_block_count - 1) c_input_size = input_size - (host_buffer_size * (i+j));
                inaccel::request req("com.xilinx.vitis.dataCompression.lz4.compressStream");
                req.arg(blocks_in[i+j]);
                req.arg(blocks_out[i+j]);
                req.arg(compressed_block_sizes_out[i+j]);
                req.arg(c_input_size);
                responses[j] = inaccel::submit(req);
            }
        }

        for (unsigned j = 0; j < m_req_num; j++)
            if(i+j < total_block_count)
                responses[j].get();

        auto kernel_end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration<double, std::nano>(kernel_end - kernel_start);
        kernel_time_ns_1 += duration;
    }

    // read the data to output buffer
    uint64_t outIdx = 0;
    for (uint64_t i = 0; i < total_block_count; ++i) {
        // copy the compressed data to out pointer
        uint32_t compressedSize = compressed_block_sizes_out[i][0];
        // std::cout << "Compressed size: " << compressedSize << std::endl;
        // current block input size
        uint32_t c_input_size = host_buffer_size;
        if (i == total_block_count - 1) c_input_size = input_size - (host_buffer_size * i);

        if (c_input_size > compressedSize) {
            // copy the compressed data
            std::memcpy(out + outIdx, &compressedSize, 4);
            outIdx += 4;
            std::memcpy(out + outIdx, blocks_out[i].data(), compressedSize);
            outIdx += compressedSize;
        } else {
            // copy the original data, since no compression
            if (c_input_size == host_buffer_size) {
                out[outIdx++] = 0;
                out[outIdx++] = 0;

                switch (c_input_size) {
                    case MAX_BSIZE_64KB:
                        out[outIdx++] = BSIZE_NCOMP_64;
                        break;
                    case MAX_BSIZE_256KB:
                        out[outIdx++] = BSIZE_NCOMP_256;
                        break;
                    case MAX_BSIZE_1024KB:
                        out[outIdx++] = BSIZE_NCOMP_1024;
                        break;
                    case MAX_BSIZE_4096KB:
                        out[outIdx++] = BSIZE_NCOMP_4096;
                        break;
                    default:
                        out[outIdx++] = BSIZE_NCOMP_64;
                        break;
                }

                out[outIdx++] = NO_COMPRESS_BIT;
            } else {
                uint8_t tmp = c_input_size;
                out[outIdx++] = tmp;
                tmp = c_input_size >> 8;
                out[outIdx++] = tmp;
                tmp = c_input_size >> 16;
                out[outIdx++] = tmp;
                out[outIdx++] = NO_COMPRESS_BIT;
            }
            std::memcpy(out + outIdx, blocks_in[i].data(), c_input_size);
            outIdx += c_input_size;
        }
        blocks_in[i].resize(0);
        blocks_out[i].resize(0);
        compressed_block_sizes_out[i].resize(0);
        blocks_in[i].shrink_to_fit();
        blocks_out[i].shrink_to_fit();
        compressed_block_sizes_out[i].shrink_to_fit();
    }

    float throughput_in_mbps_1 = (float)input_size * 1000 / kernel_time_ns_1.count();
    std::cout << std::fixed << std::setprecision(2) << "KT(MBps)\t\t:" << throughput_in_mbps_1 << std::endl;

    return outIdx;
}

uint64_t xfLz4Streaming::decompressFile(std::string& inFile_name, std::string& outFile_name, uint64_t input_size) {
    if (m_switch_flow == 0) {
        std::ifstream inFile(inFile_name.c_str(), std::ifstream::binary);
        std::ofstream outFile(outFile_name.c_str(), std::ofstream::binary);

        if (!inFile) {
            std::cout << "Unable to open file";
            exit(1);
        }

        std::vector<uint8_t> in(input_size);

        // Read magic header 4 bytes
        char c = 0;
        char magic_hdr[] = {MAGIC_BYTE_1, MAGIC_BYTE_2, MAGIC_BYTE_3, MAGIC_BYTE_4};
        for (uint32_t i = 0; i < MAGIC_HEADER_SIZE; i++) {
            inFile.get(c);
            if (c == magic_hdr[i])
                continue;
            else {
                std::cout << "Problem with magic header " << c << " " << i << std::endl;
                exit(1);
            }
        }

        // Header Checksum
        inFile.get(c);

        // Check if block size is 64 KB
        inFile.get(c);
        // printf("block_size %d \n", c);

        switch (c) {
            case BSIZE_STD_64KB:
                m_BlockSizeInKb = 64;
                break;
            case BSIZE_STD_256KB:
                m_BlockSizeInKb = 256;
                break;
            case BSIZE_STD_1024KB:
                m_BlockSizeInKb = 1024;
                break;
            case BSIZE_STD_4096KB:
                m_BlockSizeInKb = 4096;
                break;
            default:
                std::cout << "Invalid Block Size" << std::endl;
                break;
        }
        // printf("m_BlockSizeInKb %d \n", m_BlockSizeInKb);

        // Original size
        uint64_t original_size = 0;
        inFile.read((char*)&original_size, 8);
        inFile.get(c);
        // printf("original_size %d \n", original_size);

        // Allocat output size
        std::vector<uint8_t > out(original_size);

        // Read block data from compressed stream .lz4
        inFile.read((char*)in.data(), (input_size - 15));

        // uint32_t host_buffer_size = (m_BlockSizeInKb * 1024) * 32;

        // if ((m_BlockSizeInKb * 1024) > original_size) host_buffer_size = m_BlockSizeInKb * 1024;

        uint64_t debytes;
        // Decompression Streaming
        debytes = decompressStream(in.data(), out.data(), (input_size - 15), original_size);
        outFile.write((char*)out.data(), debytes);
        // Close file
        inFile.close();
        outFile.close();
        return debytes;
    } else {
        std::string command = "../../../common/lz4/lz4 --content-size -f -q -d " + inFile_name;
        system(command.c_str());
        return 0;
    }
}

// Note: Various block sizes supported by LZ4 standard are not applicable to
// this function. It just supports Block Size 64KB
uint64_t xfLz4Streaming::decompressStream(uint8_t* in, uint8_t* out, uint64_t input_size, uint64_t original_size) {
    uint32_t host_buffer_size = m_BlockSizeInKb * 1024;
    uint32_t total_block_count = (original_size - 1) / host_buffer_size + 1;

    inaccel::allocator<uint8_t> uint8_alloc;

    std::vector<uint8_t*> blocks_in(total_block_count);
    std::vector<uint8_t*> blocks_out(total_block_count);
    std::vector<uint32_t> compressedSizes(total_block_count);

    uint64_t cIdx = 0;
    uint64_t total_decompressed_size = 0;
    std::chrono::duration<double, std::nano> kernel_time_ns_1(0);

    // copy input to input buffer
    for (uint64_t i = 0; i < total_block_count; i++) {
        // get the size of the compressed block
        compressedSizes[i] = 0;

        std::memcpy(&compressedSizes[i], in + cIdx, 4);
        cIdx += 4;
        uint32_t tmp = compressedSizes[i];
        tmp >>= 24;

        if (tmp == 128) {
            uint8_t b1 = compressedSizes[i];
            uint8_t b2 = compressedSizes[i] >> 8;
            uint8_t b3 = compressedSizes[i] >> 16;
            // uint8_t b4 = compressedSizes[i] >> 24;

            if (b3 == 1) {
                compressedSizes[i] = host_buffer_size;
            } else {
                uint16_t size = 0;
                size = b2;
                size <<= 8;
                uint16_t temp = b1;
                size |= temp;
                compressedSizes[i] = size;
            }
        }
        // decompressed block input size
        uint32_t dBlockSize = host_buffer_size;
        if (i == total_block_count - 1) dBlockSize = original_size - (host_buffer_size * i);

        // If compressed size is less than original block size
        if (compressedSizes[i] < dBlockSize) {
            blocks_in[i] = uint8_alloc.allocate(compressedSizes[i]);
            blocks_out[i] = uint8_alloc.allocate(dBlockSize);
            // copy data to buffer
            std::memcpy(blocks_in[i], in + cIdx, compressedSizes[i]);
            cIdx += compressedSizes[i];
        } else if (compressedSizes[i] == dBlockSize) {
            blocks_in[i] = NULL;
            blocks_out[i] = NULL;
            // no compression, copy as it is to output
            std::memcpy(out + (host_buffer_size * i), in + cIdx, dBlockSize);
            cIdx += dBlockSize;
            total_decompressed_size += dBlockSize;
        } else {
            assert(0);
        }
    }

    std::vector<std::future<void>> responses(m_req_num);

    // send requests to the Coral manager, in blocks of m_req_num async requests
    for (uint64_t i = 0; i < total_block_count; i+=m_req_num) {
        auto kernel_start = std::chrono::high_resolution_clock::now();

        for (unsigned j = 0; j < m_req_num; j++) {
            if(i+j < total_block_count && blocks_in[i+j] != NULL) {
                uint32_t dBlockSize = host_buffer_size;
                if (i+j == total_block_count - 1) dBlockSize = original_size - (host_buffer_size * (i+j));

                inaccel::request req("com.xilinx.vitis.dataCompression.lz4.decompressStream");
                req.arg_array(&blocks_in[i+j], &blocks_in[i+j] + compressedSizes[i+j]);
                req.arg_array(&blocks_out[i+j], &blocks_out[i+j] + dBlockSize);
                req.arg(dBlockSize);
                req.arg(compressedSizes[i+j]);
                responses[j] = inaccel::submit(req);
            }
        }

        for (unsigned j = 0; j < m_req_num; j++)
            if(i+j < total_block_count && blocks_in[i+j] != NULL)
                responses[j].get();

        auto kernel_end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration<double, std::nano>(kernel_end - kernel_start);
        kernel_time_ns_1 += duration;
    }


    for (uint64_t i = 0; i < total_block_count; i++) {
        if(blocks_out[i] != NULL) {
            uint32_t dBlockSize = host_buffer_size;
            if (i == total_block_count - 1) dBlockSize = original_size - (host_buffer_size * i);
            // copy data to output
            std::memcpy(out + (host_buffer_size * i), blocks_out[i], dBlockSize);

            total_decompressed_size += dBlockSize;
            uint8_alloc.deallocate(blocks_in[i], dBlockSize);
            uint8_alloc.deallocate(blocks_out[i], dBlockSize);
            // std::cout << "Compressed size: " << compressedSize << std::endl;
        }
    }


    // for (uint64_t buf_indx = 0, blk_idx = 0; buf_indx < original_size; buf_indx += host_buffer_size, ++blk_idx) {
    //     // get the size of the compressed block
    //     uint32_t compressedSize = 0;

    //     std::memcpy(&compressedSize, in + cIdx, 4);
    //     cIdx += 4;

    //     uint32_t tmp = compressedSize;
    //     tmp >>= 24;

    //     if (tmp == 128) {
    //         uint8_t b1 = compressedSize;
    //         uint8_t b2 = compressedSize >> 8;
    //         uint8_t b3 = compressedSize >> 16;
    //         // uint8_t b4 = compressedSize >> 24;

    //         if (b3 == 1) {
    //             compressedSize = host_buffer_size;
    //         } else {
    //             uint16_t size = 0;
    //             size = b2;
    //             size <<= 8;
    //             uint16_t temp = b1;
    //             size |= temp;
    //             compressedSize = size;
    //         }
    //     }
    //     // decompressed block input size
    //     uint32_t dBlockSize = host_buffer_size;
    //     if (blk_idx == total_block_count - 1) dBlockSize = original_size - (host_buffer_size * blk_idx);

    //     // If compressed size is less than original block size
    //     if (compressedSize < dBlockSize) {
       //      blocks_in[blk_idx].resize(compressedSize);
       //      blocks_out[blk_idx].resize(dBlockSize);
    //         // copy data to buffer
    //         std::memcpy(blocks_in[blk_idx].data(), in + cIdx, compressedSize);
    //         // set kernel arguments
    //         inaccel::request req("com.xilinx.vitis.dataCompression.lz4.decompressStream");

    //         req.arg(blocks_in[blk_idx]);
    //         req.arg(blocks_out[blk_idx]);
    //         req.arg(dBlockSize);
    //         req.arg(compressedSize);

    //         auto kernel_start = std::chrono::high_resolution_clock::now();

    //         inaccel::wait(inaccel::submit(req));

    //         auto kernel_end = std::chrono::high_resolution_clock::now();
    //         auto duration = std::chrono::duration<double, std::nano>(kernel_end - kernel_start);
    //         kernel_time_ns_1 += duration;

    //         // copy data to output
    //         std::memcpy(out + buf_indx, blocks_out[blk_idx].data(), dBlockSize);
    //         cIdx += compressedSize;

    //         // std::cout << "Compressed size: " << compressedSize << std::endl;
    //     } else if (compressedSize == dBlockSize) {
    //         // no compression, copy as it is to output
    //         std::memcpy(out + buf_indx, in + cIdx, dBlockSize);
    //         cIdx += dBlockSize;
    //     } else {
    //         assert(0);
    //     }
    //     total_decompressed_size += dBlockSize;
    // }
    float throughput_in_mbps_1 = (float)total_decompressed_size * 1000 / kernel_time_ns_1.count();
    std::cout << std::fixed << std::setprecision(2) << "KT(MBps)\t\t:" << throughput_in_mbps_1 << std::endl;

    return original_size;
} // End of decompress
