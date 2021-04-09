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
#include "lz4.hpp"
#include "xxhash.h"
#include <inaccel/coral>

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
xfLz4::xfLz4(uint32_t req_num): m_req_num(req_num) {}

// Destructor
xfLz4::~xfLz4() {}

uint64_t xfLz4::compressFile(std::string& inFile_name, std::string& outFile_name, uint64_t input_size) {
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

        uint8_t temp_buff[10] = {FLG_BYTE,         block_size_header, input_size,       input_size >> 8,
                                 input_size >> 16, input_size >> 24,  input_size >> 32, input_size >> 40,
                                 input_size >> 48, input_size >> 56};

        // xxhash is used to calculate hash value
        uint32_t xxh = XXH32(temp_buff, 10, 0);
        uint64_t enbytes;
        outFile.write((char*)&temp_buff[2], 8);

        // Header CRC
        outFile.put((uint8_t)(xxh >> 8));
        // LZ4 multiple/single cu sequential version
        enbytes = compressSequential(in.data(), out.data(), input_size);
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
uint64_t xfLz4::compressSequential(uint8_t* in, uint8_t* out, uint64_t input_size) {
    uint32_t block_size_in_bytes = m_BlockSizeInKb * 1024;

    uint32_t blocks_per_chunk = 1;
    if(block_size_in_bytes < input_size) blocks_per_chunk = 32;

    uint32_t chunk_size = block_size_in_bytes * blocks_per_chunk;

    uint32_t chunk_num = (input_size - 1) / chunk_size + 1;

    std::vector<inaccel::vector<uint8_t>> buffers_in(chunk_num);
    std::vector<inaccel::vector<uint8_t>> buffers_out(chunk_num);
    std::vector<inaccel::vector<uint32_t>> compressed_block_sizes_out(chunk_num);
    std::vector<inaccel::vector<uint32_t>> block_sizes_in(chunk_num);

    std::chrono::duration<double, std::nano> kernel_time_ns_1(0);

    std::vector<std::future<void>> responses(m_req_num);

    for (uint32_t chunkIdx = 0; chunkIdx < chunk_num; ++chunkIdx) {
        uint32_t curr_chunk_size = chunk_size;
        if(chunkIdx*chunk_size + curr_chunk_size > input_size)
            curr_chunk_size = input_size - chunkIdx*chunk_size;

        uint32_t curr_blocks = (curr_chunk_size - 1) / block_size_in_bytes + 1;

        buffers_in[chunkIdx].resize(curr_chunk_size);
        buffers_out[chunkIdx].resize(curr_chunk_size);
        std::memcpy(buffers_in[chunkIdx].data(), &in[chunkIdx*chunk_size], curr_chunk_size);

        compressed_block_sizes_out[chunkIdx].resize(curr_blocks);
        block_sizes_in[chunkIdx].resize(curr_blocks);

        for (uint32_t blockIdx = 0; blockIdx < curr_blocks; ++blockIdx) {
            uint32_t block_size = block_size_in_bytes;
            if (blockIdx*block_size_in_bytes + block_size > curr_chunk_size) {
                block_size = curr_chunk_size - blockIdx*block_size_in_bytes;
            }
            block_sizes_in[chunkIdx][blockIdx] = block_size;
        }
    }

    for (uint32_t chunkIdx = 0; chunkIdx < chunk_num; chunkIdx+=m_req_num) {

        auto kernel_start = std::chrono::high_resolution_clock::now();

        for (unsigned req = 0; req < m_req_num; req++) {
            if(chunkIdx+req < chunk_num) {
                inaccel::request req("com.xilinx.vitis.dataCompression.lz4.compress");
                req.arg(buffers_in[chunkIdx+req]);
                req.arg(buffers_out[chunkIdx+req]);
                req.arg(compressed_block_sizes_out[chunkIdx+req]);
                req.arg(block_sizes_in[chunkIdx+req]);
                req.arg(m_BlockSizeInKb);
                req.arg((uint32_t)buffers_in[chunkIdx+req].size());

                responses[req] = inaccel::submit(req);
            }
        }

        for (unsigned req = 0; req < m_req_num; req++) {
            if(chunkIdx+req < chunk_num)
                responses[req].get();
        }

        auto kernel_end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration<double, std::nano>(kernel_end - kernel_start);
        kernel_time_ns_1 += duration;
    }

    // Keeps track of output buffer index
    uint64_t outIdx = 0;
    for (uint32_t chunkIdx = 0; chunkIdx < chunk_num; ++chunkIdx) {
        uint32_t curr_chunk_size = chunk_size;
        if(chunkIdx*chunk_size + curr_chunk_size > input_size)
            curr_chunk_size = input_size - chunkIdx*chunk_size;

        uint32_t curr_blocks = (curr_chunk_size - 1) / block_size_in_bytes + 1;

        for (uint32_t blockIdx = 0; blockIdx < curr_blocks; ++blockIdx) {
            uint32_t block_size = block_size_in_bytes;
            if (blockIdx*block_size_in_bytes + block_size > curr_chunk_size) {
                block_size = curr_chunk_size - blockIdx*block_size_in_bytes;
            }

            uint32_t compressed_size = compressed_block_sizes_out[chunkIdx][blockIdx];
            assert(compressed_size != 0);

            int perc_cal = curr_chunk_size * 10;
            perc_cal = perc_cal / block_size;

            if (compressed_size < block_size && perc_cal >= 10) {
                std::memcpy(out + outIdx, &compressed_size, 4);
                outIdx += 4;
                std::memcpy(out + outIdx, &(buffers_out[chunkIdx][blockIdx*block_size_in_bytes]), compressed_size);
                outIdx += compressed_size;
            } else {
                if (block_size == block_size_in_bytes) {
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
                std::memcpy(out + outIdx, &(buffers_in[chunkIdx][blockIdx*block_size_in_bytes]), block_size);
                outIdx += block_size;
            }
        } // End of chunk (block by block) copy to output buffer


        buffers_in[chunkIdx].resize(0);
        buffers_out[chunkIdx].resize(0);
        compressed_block_sizes_out[chunkIdx].resize(0);
        block_sizes_in[chunkIdx].resize(0);
        blocks_in[chunkIdx].shrink_to_fit();
        buffers_out[chunkIdx].shrink_to_fit();
        compressed_block_sizes_out[i].shrink_to_fit();
        block_sizes_in[chunkIdx].shrink_to_fit();
    }
    float throughput_in_mbps_1 = (float)input_size * 1000 / kernel_time_ns_1.count();
    std::cout << std::fixed << std::setprecision(2) << "KT(MBps)\t\t:" << throughput_in_mbps_1 << std::endl;
    return outIdx;
}

uint64_t xfLz4::decompressFile(std::string& inFile_name, std::string& outFile_name, uint64_t input_size) {
    if (m_switch_flow == 0) {
        std::ifstream inFile(inFile_name.c_str(), std::ifstream::binary);
        std::ofstream outFile(outFile_name.c_str(), std::ofstream::binary);

        if (!inFile) {
            std::cout << "Unable to open file";
            exit(1);
        }

        std::vector<uint8_t > in(input_size);

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
        std::vector<uint8_t> out(original_size);

        // Read block data from compressed stream .lz4
        inFile.read((char*)in.data(), (input_size - 15));

        uint64_t debytes;
        // Decompression Sequential multiple cus.
        debytes = decompressSequential(in.data(), out.data(), (input_size - 15), original_size);
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
uint64_t xfLz4::decompressSequential(
    uint8_t* in, uint8_t* out, uint64_t input_size, uint64_t original_size) {
    uint32_t block_size_in_bytes = m_BlockSizeInKb * 1024;

    uint32_t blocks_per_chunk = 1;
    if(block_size_in_bytes < input_size) blocks_per_chunk = 32;

    uint32_t chunk_size = block_size_in_bytes * blocks_per_chunk;

    uint32_t chunk_num = (input_size - 1) / chunk_size + 1;

    uint32_t total_block_count = (original_size - 1) / block_size_in_bytes + 1;
    uint32_t block_cntr = 0;
    uint32_t done_block_cntr = 0;

    std::vector<inaccel::vector<uint8_t>> buffers_in(chunk_num);
    std::vector<inaccel::vector<uint8_t>> buffers_out(chunk_num);
    std::vector<inaccel::vector<uint32_t>> compressed_block_sizes_in(chunk_num);
    std::vector<inaccel::vector<uint32_t>> block_sizes_in(chunk_num);

    uint32_t no_compress_case = 0;
    std::chrono::duration<double, std::nano> kernel_time_ns_1(0);

    std::vector<std::future<void>> responses(m_req_num);

    uint64_t inIdx = 0;
    uint64_t total_decompressed_size = 0;

    for (uint32_t chunkIdx = 0; chunkIdx < chunk_num; ++chunkIdx) {
        buffers_in[chunkIdx].resize(blocks_per_chunk*block_size_in_bytes);
        buffers_out[chunkIdx].resize(blocks_per_chunk*block_size_in_bytes);
        compressed_block_sizes_in[chunkIdx].resize(blocks_per_chunk);
        block_sizes_in[chunkIdx].resize(blocks_per_chunk);

        uint64_t curr_chunk_size = 0;
        uint64_t cIdx = 0;

        //loop to find chunk size and compressed block sizes
        for (uint32_t blockIdx = 0; blockIdx < blocks_per_chunk; ++blockIdx) {
            std::memcpy(&compressed_block_sizes_in[chunkIdx][blockIdx], &in[inIdx], 4);
            inIdx += 4;
            uint32_t tmp = compressed_block_sizes_in[chunkIdx][blockIdx];
            tmp >>= 24;

            if (tmp == 128) {
                uint8_t b1 = compressed_block_sizes_in[chunkIdx][blockIdx];
                uint8_t b2 = compressed_block_sizes_in[chunkIdx][blockIdx] >> 8;
                uint8_t b3 = compressed_block_sizes_in[chunkIdx][blockIdx] >> 16;
                // uint8_t b4 = compressed_block_sizes_in[chunkIdx][blockIdx] >> 24;

                if (b3 == 1) {
                    compressed_block_sizes_in[chunkIdx][blockIdx] = block_size_in_bytes;
                } else {
                    uint16_t size = 0;
                    size = b2;
                    size <<= 8;
                    uint16_t temp = b1;
                    size |= temp;
                    compressed_block_sizes_in[chunkIdx][blockIdx] = size;
                }
            }

            curr_chunk_size += compressed_block_sizes_in[chunkIdx][blockIdx];

            block_sizes_in[chunkIdx][blockIdx] = block_size_in_bytes;
            if(chunkIdx*blocks_per_chunk + blockIdx == total_block_count - 1) {
                block_sizes_in[chunkIdx][blockIdx] = original_size -
                    (chunkIdx*blocks_per_chunk + blockIdx)*block_size_in_bytes;
            }

            // If compressed size is less than original block size
            if (compressed_block_sizes_in[chunkIdx][blockIdx] < block_sizes_in[chunkIdx][blockIdx]) {
                std::memcpy(&(buffers_in[chunkIdx][blockIdx*block_size_in_bytes]),
                            &in[inIdx], compressed_block_sizes_in[chunkIdx][blockIdx]);
                inIdx += compressed_block_sizes_in[chunkIdx][blockIdx];
            } else if (compressed_block_sizes_in[chunkIdx][blockIdx] == block_sizes_in[chunkIdx][blockIdx]) {
                no_compress_case++;
                // No compression block
                std::memcpy(out + (chunkIdx*blocks_per_chunk + blockIdx)*block_size_in_bytes,
                            in + inIdx, block_sizes_in[chunkIdx][blockIdx]);
                total_decompressed_size += block_sizes_in[chunkIdx][blockIdx];
            } else {
                assert(0);
            }
            inIdx += compressed_block_sizes_in[chunkIdx][blockIdx];

            if(inIdx >= input_size || (chunkIdx*blocks_per_chunk + blockIdx == total_block_count - 1))
                break;
        }
    }


    for (uint32_t chunkIdx = 0; chunkIdx < chunk_num; chunkIdx+=m_req_num) {

        auto kernel_start = std::chrono::high_resolution_clock::now();
        for (unsigned req = 0; req < m_req_num; req++) {
            if(chunkIdx+req < chunk_num) {
                uint32_t bufblocks = 0;
                for (uint32_t blockIdx = 0; blockIdx < blocks_per_chunk; ++blockIdx)
                    if (compressed_block_sizes_in[chunkIdx+req][blockIdx] <
                            block_sizes_in[chunkIdx+req][blockIdx])
                        bufblocks++;

                if(bufblocks) {
                    inaccel::request req("com.xilinx.vitis.dataCompression.lz4.decompress");
                    req.arg(buffers_in[chunkIdx+req]);
                    req.arg(buffers_out[chunkIdx+req]);
                    req.arg(block_sizes_in[chunkIdx+req]);
                    req.arg(compressed_block_sizes_in[chunkIdx+req]);
                    req.arg(m_BlockSizeInKb);
                    req.arg(bufblocks);

                    responses[req] = inaccel::submit(req);
                }
                else responses[req] = NULL;
            }
        }

        for (unsigned req = 0; req < m_req_num; req++) {
            if(chunkIdx+req < chunk_num && responses[req])
                responses[req].get();
        }

        auto kernel_end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration<double, std::nano>(kernel_end - kernel_start);
        kernel_time_ns_1 += duration;
    }

    for (uint32_t chunkIdx = 0; chunkIdx < chunk_num; ++chunkIdx) {
        for (uint32_t blockIdx = 0; blockIdx < blocks_per_chunk; ++blockIdx) {
            if (compressed_block_sizes_in[chunkIdx][blockIdx] <
                    block_sizes_in[chunkIdx][blockIdx]) {
                std::memcpy(out + (chunkIdx*blocks_per_chunk + blockIdx)*block_size_in_bytes,
                            buffers_out[chunkIdx] + blockIdx*block_size_in_bytes,
                            block_sizes_in[chunkIdx][blockIdx]);
                total_decompressed_size += block_sizes_in[chunkIdx][blockIdx];
            }
        }

        buffers_in[chunkIdx].resize(0);
        buffers_out[chunkIdx].resize(0);
        block_sizes_in[chunkIdx].resize(0);
        compressed_block_sizes_in[chunkIdx].resize(0);
        blocks_in[chunkIdx].shrink_to_fit();
        buffers_out[chunkIdx].shrink_to_fit();
        block_sizes_in[chunkIdx].shrink_to_fit();
        compressed_block_sizes_in[i].shrink_to_fit();
    }

    float throughput_in_mbps_1 = (float)total_decompressed_size * 1000 / kernel_time_ns_1.count();
    std::cout << std::fixed << std::setprecision(2) << "KT(MBps)\t\t:" << throughput_in_mbps_1 << std::endl;

    return original_size;
} // End of decompress
