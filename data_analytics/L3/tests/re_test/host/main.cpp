/*
 * Copyright 2020 Xilinx, Inc.
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
#include <cstdlib>
#include <fstream>

#include "xf_data_analytics/text/regex_engine.hpp"
#include "general_config.hpp"

// for validating result.
extern "C" {
#include "oniguruma.h"
}

enum {
    MAX_MSG_DEPTH = 250000000,   // Max number of messages in a section
    MAX_MSG_LEN = 65536,         // Max length of message in byte
    MAX_LNM = 6000000,           // Max number of lines in a single section
    MAX_OUT_DEPTH = MAX_LNM * 20 // 20 for 19 capturing groups at most
};
int check_result(std::string pattern,
                 uint64_t* msg_buff,
                 uint32_t* offt_buff,
                 uint16_t* len_buff,
                 uint32_t* out_buff,
                 uint32_t lnm,
                 uint32_t cpgp_nm) {
    int r;
    unsigned char *start, *range, *end;
    regex_t* reg;
    OnigErrorInfo einfo;
    OnigRegion* region = onig_region_new();
    OnigEncoding use_encs[1];

    use_encs[0] = ONIG_ENCODING_ASCII;
    onig_initialize(use_encs, sizeof(use_encs) / sizeof(use_encs[0]));

    UChar* pattern_c = (UChar*)pattern.c_str();

    r = onig_new(&reg, pattern_c, pattern_c + strlen((char*)pattern_c), ONIG_OPTION_DEFAULT, ONIG_ENCODING_ASCII,
                 ONIG_SYNTAX_DEFAULT, &einfo);
    if (r != ONIG_NORMAL) {
        char s[ONIG_MAX_ERROR_MESSAGE_LEN];
        onig_error_code_to_str((UChar*)s, r, &einfo);
        fprintf(stderr, "ERROR: %s\n", s);
        return -1;
    }
    unsigned char* max_str = (unsigned char*)malloc(MAX_MSG_LEN);
    for (unsigned i = 0; i < lnm; ++i) {
        // generate referecne
        int offt = offt_buff[i];
        memcpy(max_str, &msg_buff[offt], len_buff[i]);
        max_str[len_buff[i]] = '\0';
        UChar* str = (UChar*)max_str;
        end = str + strlen((char*)str);
        start = str;
        range = end;
        r = onig_search(reg, str, end, start, range, region, ONIG_OPTION_NONE);
        // compare with actual result
        // match
        if (r == 0) {
            if (out_buff[i * (cpgp_nm + 1)] == 1) {
                for (unsigned j = 0; j < cpgp_nm; ++j) {
                    if ((unsigned) region->beg[j] != out_buff[i * (cpgp_nm + 1) + j + 1] % 65536 ||
                        (unsigned) region->end[j] != out_buff[i * (cpgp_nm + 1) + j + 1] / 65536) {
                        fprintf(stderr, "ERROR: msg: %d, capture group: %d, ref:[%d, %d], act:[%d, %d]\n", i, j,
                                region->beg[j], region->end[j], out_buff[i * (cpgp_nm + 1) + j + 1] % 65536,
                                out_buff[i * (cpgp_nm + 1) + j + 1] / 65536);
                        return -1;
                    }
                }
            } else {
                fprintf(stderr, "ERROR: msg: %d, ref: %d, act: %d\n", i, 1, out_buff[i * (cpgp_nm + 1)]);
                return -1;
            }
            // mismatch
        } else if (r == ONIG_MISMATCH) {
            if (out_buff[i * (cpgp_nm + 1)] != 0) {
                fprintf(stderr, "ERROR: msg: %d, ref: %d, act: %d\n", i, 0, out_buff[i * (cpgp_nm + 1)]);
                return -1;
            }
        } else {
            char s[ONIG_MAX_ERROR_MESSAGE_LEN];
            onig_error_code_to_str((UChar*)s, r);
            fprintf(stderr, "ERROR: %s\n", s);
            onig_region_free(region, 1 /* 1:free self, 0:free contents only */);
            onig_free(reg);
            onig_end();
            return -1;
        }
    }
    onig_region_free(region, 1 /* 1:free self, 0:free contents only */);
    onig_free(reg);
    onig_end();
    return 0;
}
void store_dat(std::ofstream& out_file, uint32_t* out_buff, uint32_t lnm, uint32_t cpgp_nm) {
    typedef union {
        int16_t int_a[2];
        uint32_t d;
    } uint32_un;
    for (unsigned int i = 0; i < lnm; ++i) {
        for (unsigned int j = 0; j < cpgp_nm + 1; ++j) {
            if (i * (cpgp_nm + 1) + j < MAX_OUT_DEPTH) {
                uint32_t out = out_buff[i * (cpgp_nm + 1) + j];
                // match result
                if (j == 0) {
                    if (out == 0)
                        out_file << "Mismatch\n";
                    else
                        out_file << "Match\n";
                } else {
                    // offset of capture group
                    uint32_un tmp;
                    tmp.d = out;
                    out_file << j - 1 << ": [" << tmp.int_a[0] << ", " << tmp.int_a[1] << "]\n";
                }
            }
        }
    }
}
int load_dat(
    std::ifstream& log_file, uint64_t* msg_buff, uint32_t* offt_buff, uint16_t* len_buff, uint32_t& lnm, int limit_ln) {
    typedef union {
        char c_a[8];
        uint64_t d;
    } uint64_un;

    lnm = 0;
    std::string line;
    uint32_t offt = 0;
    while (!log_file.eof() && (offt < (MAX_MSG_DEPTH - MAX_MSG_LEN / 8)) && lnm < MAX_LNM &&
           (limit_ln == -1 || (int) lnm < limit_ln)) {
        getline(log_file, line);
        size_t sz = line.size();
        // max line
        if (sz >= MAX_MSG_LEN) {
            std::cerr << "ERROR: length of line exceeds " << MAX_MSG_LEN << ".\n";
            return -1;
            // ignore empty line
        } else if (sz > 0) {
            offt_buff[lnm] = offt;
            len_buff[lnm] = sz;
            for (unsigned i = 0; i < (sz + 7) / 8; ++i) {
                uint64_un out;
                for (unsigned int j = 0; j < 8; ++j) {
                    if (i * 8 + j < sz) {
                        out.c_a[j] = line[i * 8 + j];
                    } else {
                        out.c_a[j] = ' ';
                    }
                }
                msg_buff[offt++] = out.d;
            }
            lnm++;
        }
    }
    // one more
    offt_buff[lnm] = offt;
    return 0;
}
int main(int argc, const char* argv[]) {
    std::cout << "----------------------------------------log analytics with regex----------------------------------------" << std::endl;
    // command argument parser
    // TODO use new argument parser from Utility library.
    xf::data_analytics::text::details::ArgParser parser(argc, argv);

    std::string log_path;
    if (!parser.getCmdOption("-in", log_path)) {
        std::cout << "ERROR:  input log path is not specified.\n";
        return 1;
    }
    std::string out_path;
    if (!parser.getCmdOption("-out", out_path)) {
        std::cout << "ERROR:  output path is not specified.\n";
        return 1;
    }
    std::string ln_nm;
    int limit_ln = -1;
    if (parser.getCmdOption("-lnm", ln_nm)) {
        try {
            limit_ln = std::stoi(ln_nm);
        } catch (...) {
            limit_ln = -1;
        }
    }
    std::string pattern =
        "^(?<remote>[^ ]*) (?<host>[^ ]*) (?<user>[^ ]*) \\[(?<time>[^\\]]*)\\] \"(?<method>\\S+)(?: "
        "+(?<path>[^\\\"]*?)(?: +\\S*)?)?\" (?<code>[^ ]*) (?<size>[^ ]*)(?: \"(?<referer>[^\\\"]*)\" "
        "\"(?<agent>[^\\\"]*)\"(?:\\s+(?<http_x_forwarded_for>[^ ]+))?)?$";

    // allocate the in-memory buffer
    uint64_t* msg_buff = (uint64_t*) malloc(MAX_MSG_DEPTH * sizeof(uint64_t));
    uint32_t* offt_buff = (uint32_t*) malloc(MAX_LNM * sizeof(uint32_t));
    uint16_t* len_buff = (uint16_t*) malloc(MAX_LNM * sizeof(uint16_t));
    uint32_t* out_buff = (uint32_t*) malloc(MAX_OUT_DEPTH * sizeof(uint32_t));
    // constructor of reEngine
    xf::data_analytics::text::re::RegexEngine reInst(INSTR_DEPTH, CCLASS_NM, CPGP_NM, MSG_SZ, // re limits
                                                     SLICE_SZ, SLICE_NM);                     // processing
    xf::data_analytics::text::re::ErrCode err_code;
    // compile pattern
    err_code = reInst.compile(pattern);
    if (err_code != 0) return -1;

    // get capture group number
    uint16_t cpgp_nm = reInst.getCpgpNm();

    // load data from disk to in-memory buffer
    std::ifstream log_file(log_path);
    std::ofstream out_file(out_path);
    uint32_t lnm = 0;
    if (!log_file.is_open()) {
        std::cerr << "ERROR: " << log_path << " cannot be opened for read.\n";
        return -1;
    }
    if (!out_file.is_open()) {
        std::cerr << "ERROR: " << out_path << " cannot be opened for write.\n";
        return -1;
    }
    while (!log_file.eof() && (limit_ln == -1 || (int) lnm < limit_ln)) {
        // load data
        if (load_dat(log_file, msg_buff, offt_buff, len_buff, lnm, limit_ln) == -1) {
            return -1;
        }
        if (lnm > 0) {
            struct timeval re_start, re_end;
            // call reInst to do regex
            gettimeofday(&re_start, 0);
            err_code = reInst.match(lnm, msg_buff, offt_buff, len_buff, out_buff);
            gettimeofday(&re_end, 0);
            double re_tvtime = xf::data_analytics::text::details::tvdiff(re_start, re_end);
            double total_log_size = (double)offt_buff[lnm - 1] * 8 / 1024 / 1024;
            std::cout << "HW regex pipelined, time: " << (double)re_tvtime / 1000 << " ms, size: " << total_log_size
                      << " MB, throughput: " << total_log_size / 1024 / ((double)re_tvtime / 1000000) << " GB/s" << std::endl;
            if (err_code) {
                std::cerr << "ERROR: match failed.\n";
                return -1;
            }
            // write data to disk
            store_dat(out_file, out_buff, lnm, cpgp_nm);
            // if check is open, check the result with golden
            gettimeofday(&re_start, 0);
            int r = check_result(pattern, msg_buff, offt_buff, len_buff, out_buff, lnm, cpgp_nm);
            gettimeofday(&re_end, 0);
            double _re_tvtime = xf::data_analytics::text::details::tvdiff(re_start, re_end);
            std::cout << "SW regex, time: " << (double)_re_tvtime / 1000 << " ms, size: " << total_log_size
                      << " MB, throughput: " << total_log_size / 1024 / ((double)_re_tvtime / 1000000) << " GB/s" << std::endl;
            if (r == -1) {
                fprintf(stderr, "ERROR: result mismatch, speedup: %.2f\n", _re_tvtime/re_tvtime);
            } else {
                fprintf(stdout, "SUCCESS: result match, speedup: %.2f\n", _re_tvtime/re_tvtime);
            }
        }
    }
    std::cout << "--------------------------------------Finished regex pipelined test-------------------------------------" << std::endl;
    log_file.close();
    out_file.close();
    std::quick_exit(0);
    return 0;
}
