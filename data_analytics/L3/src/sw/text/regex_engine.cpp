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
#include "xf_data_analytics/text/regex_engine.hpp"

extern "C" {
#include "oniguruma.h"
}

#include <inaccel/coral>
#include <omp.h>

namespace xf {
namespace data_analytics {
namespace text {
namespace re {

// constructor
RegexEngine::RegexEngine(const int instr_depth,
                         const int char_class_num,
                         const int capture_grp_num,
                         const int msg_size,
                         const int max_slice_size,
                         const int max_slice_num)
    : reCfg(instr_depth, char_class_num, capture_grp_num),
      kInstrDepth(instr_depth),
      kCharClassNum(char_class_num),
      kCaptureGrpNum(capture_grp_num),
      kMsgSize(msg_size),
      kMaxSliceSize(max_slice_size),
      kMaxSliceNum(max_slice_num) {}

ErrCode RegexEngine::match_all(const uint64_t* msg_buff,
		uint32_t* offt_buff,
		uint16_t* len_buff,
		uint32_t* out_buff,
		const uint64_t* re_cfg,
		uint32_t total_lnm) {
	uint32_t cpgp_nm = reCfg.getCpgpNm();

	// calculate the section number
	timeval tv_start, tv_end, re_start, re_end;

	gettimeofday(&tv_start, 0);
	uint32_t max_slice_lnm = 0;

	uint16_t* lnm_per_sec = (uint16_t*) malloc(kMaxSliceNum * sizeof(uint16_t));
	uint32_t* pos_per_sec = (uint32_t*) malloc(kMaxSliceNum * sizeof(uint32_t));


	uint32_t sec_num = findSecNum(len_buff, total_lnm, &max_slice_lnm, lnm_per_sec, pos_per_sec);

	gettimeofday(&tv_end, 0);
	double tvtime = details::tvdiff(tv_start, tv_end);
	fprintf(stdout, "The log file is partitioned into %d section(s) with max_slice_lnm %d. Partitioning took %f ms.\n", sec_num, max_slice_lnm, tvtime / 1000);

	unsigned cu_num = 24;
	std::vector<uint64_t*> cfg_in_slice(cu_num);
	std::vector<uint64_t*> msg_in_slice(cu_num);
	std::vector<uint16_t*> len_in_slice(cu_num);
	std::vector<uint32_t*> out_slice(cu_num);

	for (unsigned i = 0; i < cu_num; i++) {
		cfg_in_slice[i] = (uint64_t*) inaccel_alloc((kInstrDepth + kCharClassNum * 4 + 2) * sizeof(uint64_t));
		memcpy(cfg_in_slice[i], re_cfg, kInstrDepth + kCharClassNum * 4 + 2);
		msg_in_slice[i] = (uint64_t*) inaccel_alloc((kMaxSliceSize / 8 + 1) * sizeof(uint64_t));
		len_in_slice[i] = (uint16_t*) inaccel_alloc((max_slice_lnm + 2) * sizeof(uint16_t));
		out_slice[i] = (uint32_t*) inaccel_alloc(((cpgp_nm + 1) * max_slice_lnm + 1) * sizeof(uint32_t));
	}

	// initialize oniguruma regex
	OnigRegion* region = onig_region_new();
	regex_t* reg;
	OnigEncoding use_encs[1];
	use_encs[0] = ONIG_ENCODING_ASCII;
	onig_initialize(use_encs, sizeof(use_encs) / sizeof(use_encs[0]));
	UChar* pattern_c = (UChar*)reCfg.pattern.c_str();
	OnigErrorInfo einfo;
	onig_new(&reg, pattern_c, pattern_c + strlen((char*)pattern_c), ONIG_OPTION_DEFAULT,
									ONIG_ENCODING_ASCII, ONIG_SYNTAX_DEFAULT, &einfo);

	gettimeofday(&re_start, 0);

	omp_set_num_threads(cu_num);
	#pragma omp parallel for
	for (unsigned sec = 0; sec < sec_num; sec ++) {
		int tid = omp_get_thread_num();
		size_t slice_sz = 0;
		for (int i = 0; i < lnm_per_sec[sec]; ++i) {
			size_t tmp = slice_sz;
			if (len_buff[pos_per_sec[sec] + i] < kMsgSize) {
				tmp += (len_buff[pos_per_sec[sec] + i] + 7) / 8;
			}
			// reach the end of file or stop when reach the limitation of slice size
			if (pos_per_sec[sec] + i >= total_lnm || tmp > kMaxSliceSize / 8) {
				break;
			} else {
				// if the message length exceed the maxinum, set len = 0;
				if (len_buff[pos_per_sec[sec] + i] > kMsgSize) {
					len_in_slice[tid][i + 2] = 0;
				} else {
					len_in_slice[tid][i + 2] = len_buff[pos_per_sec[sec] + i];
					memcpy(msg_in_slice[tid] + slice_sz + 1, msg_buff + offt_buff[pos_per_sec[sec] + i], len_buff[pos_per_sec[sec] + i]);
				}
				slice_sz = tmp;
			}
		}

		msg_in_slice[tid][0] = (uint64_t)(slice_sz + 1);
		len_in_slice[tid][0] = (lnm_per_sec[sec] + 2) / 65536;
		len_in_slice[tid][1] = (lnm_per_sec[sec] + 2) % 65536;

		inaccel::request reEngine("com.xilinx.dataAnalytics.reEngine");
		reEngine.arg_array(cfg_in_slice[tid], cfg_in_slice[tid] + kInstrDepth + kCharClassNum * 4 + 2)
				.arg_array(msg_in_slice[tid], msg_in_slice[tid] + kMaxSliceSize / 8 + 1)
				.arg_array(len_in_slice[tid], len_in_slice[tid] + max_slice_lnm + 2)
				.arg_array(out_slice[tid], out_slice[tid] + (cpgp_nm + 1) * max_slice_lnm + 1);

		inaccel::submit(reEngine).get();

		unsigned char* max_str = (unsigned char*)malloc(65536);
		for (unsigned int i = 0; i < lnm_per_sec[sec]; ++i) {
			uint8_t result = out_slice[tid][i * (cpgp_nm + 1) + 1];
			// step 1: stack overflow or large message
			if (result == 2 || result == 3) {
				// step 2: find the position and length message
				size_t msg_pos = offt_buff[pos_per_sec[sec] + i];
				const uint64_t* msg = &msg_buff[msg_pos];
				uint16_t msg_len = len_buff[pos_per_sec[sec] + i];
				memcpy(max_str, msg, msg_len);
				max_str[msg_len] = '\0';
				UChar* str = (UChar*)max_str;
				unsigned char* end = str + strlen((char*)str);
				unsigned char* start = str;
				unsigned char* range = end;
				int r = onig_search(reg, str, end, start, range, region, ONIG_OPTION_NONE);
				// printf("[DEBUG], post_proc: %d, r: %d\n", *(q.start_pos) + i, r);
				// step 4: insert the result back to out_buff
				if (r == 0) {
					out_slice[tid][i * (cpgp_nm + 1) + 1] = 1;
					for (unsigned j = 0; j < cpgp_nm; ++j) {
						uint32_t out = region->end[j] * 65536 + region->beg[j];
						out_slice[tid][i * (cpgp_nm + 1) + 2 + j] = out;
					}
				} else if (r == ONIG_MISMATCH) {
					out_slice[tid][i * (cpgp_nm + 1) + 1] = 0;
				}
			}
		}

		size_t sz = lnm_per_sec[sec] * (cpgp_nm + 1) * sizeof(out_slice[tid][0]);
		// copy result to output buffer
		memcpy(out_buff + pos_per_sec[sec] * (cpgp_nm + 1), out_slice[tid] + 1, sz);
	}
	gettimeofday(&re_end, 0);

	onig_region_free(region, 1);
	onig_free(reg);
	onig_end();

	double re_tvtime = details::tvdiff(re_start, re_end);
	double total_log_size = (double)offt_buff[total_lnm - 1] * 8 / 1024 / 1024;
	std::cout << "regex pipelined, time: " << (double)re_tvtime / 1000 << " ms, size: " << total_log_size
			  << " MB, throughput: " << total_log_size / 1024 / ((double)re_tvtime / 1000000) << " GB/s" << std::endl;
	std::cout
			<< "-----------------------------Finished regex pipelined test----------------------------------------------"
			<< std::endl
			<< std::endl;

	return SUCCESS;
}
ErrCode RegexEngine::compile(std::string pattern) {
    return reCfg.compile(pattern);
}
ErrCode RegexEngine::match(
    uint32_t total_lnm, const uint64_t* msg_buff, uint32_t* offt_buff, uint16_t* len_buff, uint32_t* out_buff) {
    const uint64_t* cfg_buff = reCfg.getConfigBits();
    // match
    return match_all(msg_buff, offt_buff, len_buff, out_buff, cfg_buff, total_lnm);
}
uint32_t RegexEngine::getCpgpNm() const {
    return reCfg.getCpgpNm();
}
uint32_t RegexEngine::findSecNum(
    uint16_t* len_buff, uint32_t lnm, uint32_t* slice_lnm, uint16_t* lnm_per_sec, uint32_t* pos_per_sec) {
    uint32_t sec_sz = 0;
    uint32_t sec_nm = 0;
    uint32_t start_lnm = 0;
    uint32_t end_lnm = 0;
    uint32_t tmp_slice_nm = 0;
    for (unsigned int i = 0; i < lnm; ++i) {
        if (len_buff[i] < kMsgSize) {
            sec_sz += (len_buff[i] + 7) / 8;
            if (sec_sz > kMaxSliceSize / 8) {
                start_lnm = end_lnm;
                end_lnm = i;
                if (end_lnm - start_lnm > tmp_slice_nm) tmp_slice_nm = end_lnm - start_lnm;
                lnm_per_sec[sec_nm] = end_lnm - start_lnm;
                pos_per_sec[sec_nm] = start_lnm;
                sec_nm++;
                sec_sz = (len_buff[i] + 7) / 8;
            } else if (i == lnm - 1) {
                start_lnm = end_lnm;
                end_lnm = lnm;
                lnm_per_sec[sec_nm] = end_lnm - start_lnm;
                pos_per_sec[sec_nm] = start_lnm;
                sec_nm++;
                if (end_lnm - start_lnm > tmp_slice_nm) tmp_slice_nm = end_lnm - start_lnm;
            }
        }
    }
    *slice_lnm = tmp_slice_nm + 2;
    return sec_nm;
}

} // namesapce re
} // namespace text
} // namespace data_analytics
} // namespace xf
