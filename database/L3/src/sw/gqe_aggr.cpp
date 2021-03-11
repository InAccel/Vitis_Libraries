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

#include <unordered_map>
// L2
#include "xf_database/meta_table.hpp"
#include "xf_database/gqe_utils.hpp"
// L3
#include "xf_database/gqe_aggr.hpp"

#include <inaccel/coral.h>
#include <sys/time.h>
#include <ctime>

namespace xf {
namespace database {
namespace gqe {
Aggregator::Aggregator() {}

Aggregator::~Aggregator() {};

ErrCode Aggregator::aggregate(Table& tab_in,
                              std::vector<EvaluationInfo> evals_info,
                              std::string filter_str,
                              std::string group_keys_str,
                              std::string output_str,
                              Table& tab_out,
                              AggrStrategyBase* strategyimp) {
    // strategy
    bool new_s = false;
    if (strategyimp == nullptr) {
        strategyimp = new AggrStrategyBase();
        new_s = true;
    }
    auto params = strategyimp->getSolutionParams(tab_in);
    // for deug
    // cfg
    AggrConfig aggr_config(tab_in, evals_info, filter_str, group_keys_str, output_str, (params[0] == 1));
    // join
    ErrCode err = aggr_all(tab_in, tab_out, aggr_config, params);
    if (new_s) delete strategyimp;
    return err;
}

ErrCode Aggregator::aggr_all(Table& tab_in, Table& tab_out, AggrConfig& aggr_cfg, std::vector<size_t> params) {
    ErrCode err;
    size_t _solution = params[0];
    if (_solution == 0) {
        std::cout << "direct aggregate" << std::endl;
        err = aggr_sol0(tab_in, tab_out, aggr_cfg, params);
    } else if (_solution == 1) {
        std::cout << "pipelined aggregate" << std::endl;
        err = aggr_sol1(tab_in, tab_out, aggr_cfg, params);
    } else if (_solution == 2) {
        std::cout << "partition + pipelined aggregate" << std::endl;
        err = aggr_sol2(tab_in, tab_out, aggr_cfg, params);
    } else {
        return PARAM_ERROR;
    }
    return err;
}

ErrCode Aggregator::aggr_sol0(Table& tab_in, Table& tab_out, AggrConfig& aggr_cfg, std::vector<size_t> params) {
    const int size_of_apu512 = sizeof(ap_uint<512>);
    const int VEC_LEN = size_of_apu512 / sizeof(int);
    int tab_in_col_num = tab_in.getColNum();
    int tab_in_row_num = tab_in.getRowNum();
    int tab_out_col_num = tab_out.getColNum();
    int result_nrow = tab_out.getRowNum();
    int tb_in_col_depth = (tab_in_row_num + VEC_LEN - 1) / VEC_LEN;
    int tb_in_col_size = size_of_apu512 * tb_in_col_depth;
    int* tab_in_col_type = new int[tab_in_col_num];
    int* tab_out_col_type = new int[tab_out_col_num];
    for (int i = 0; i < tab_in_col_num; i++) {
        tab_in_col_type[i] = tab_in.getColTypeSize(i);
    }
    for (int i = 0; i < tab_out_col_num; i++) {
        tab_out_col_type[i] = tab_out.getColTypeSize(i);
    }

    //--------------- get sw scan input host bufer -----------------
    std::vector<int8_t> scan_list = aggr_cfg.getScanList();
#ifdef USER_DEBUG
    std::cout << "table in info: " << tab_in_row_num << " rows, " << tab_in_col_num << " cols." << std::endl;
    std::cout << "table out info: " << result_nrow << " rows, " << tab_out_col_num << " cols." << std::endl;
    std::cout << "scan_list:" << std::endl;
    for (size_t i = 0; i < scan_list.size(); i++) {
        std::cout << (int)scan_list[i] << " ";
    }
    std::cout << std::endl;
#endif
    char* table_in_col[tab_in_col_num];
    for (int i = 0; i < tab_in_col_num; i++) {
        table_in_col[i] = tab_in.getColPointer(scan_list[i]);
    }

    //--------------- get output host bufer -----------------
    size_t table_result_depth = (result_nrow + VEC_LEN - 1) / VEC_LEN;
    size_t table_result_size = table_result_depth * size_of_apu512;
	char* table_out_col[16];
    for (int i = 0; i < 16; i++) {
        table_out_col[i] = (char*) cube_alloc(table_result_size);
	if(!table_out_col[i]){
		std::cout<< "Error at alloc" << std::endl;
		exit(1);
	}
    }

    MetaTable meta_aggr_in;
    meta_aggr_in.setColNum(tab_in_col_num);
    for (int i = 0; i < tab_in_col_num; i++) {
        meta_aggr_in.setCol(i, i, tab_in_row_num);
    }
    MetaTable meta_aggr_out;
    meta_aggr_out.setColNum(8);
    for (int i = 0; i < 8; i++) {
        meta_aggr_out.setCol(i, i, result_nrow);
    }
    ap_uint<32>* table_cfg = aggr_cfg.getAggrConfigBits();

    ap_uint<32>* table_cfg_out = aggr_cfg.getAggrConfigOutBits();

#ifdef USER_DEBUG
    std::cout << "debug 0" << std::endl;
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 4; j++) {
            std::cout << table_cfg[i].range(j * 8 + 7, j * 8) << " ";
        }
    }
    std::cout << std::endl;
#endif

    //--------------- get output host bufer -----------------
    char* mext_table_in_col[8];

    for (int i = 0; i < 8; i++) {
        mext_table_in_col[i] = (char*) cube_alloc(tb_in_col_size);
	if(!mext_table_in_col[i]){
                std::cout<< "Error at alloc" << std::endl;
                exit(1);
        }
        if( i < tab_in_col_num){
                memcpy(mext_table_in_col[i], table_in_col[i], tb_in_col_size);
        }
    }

    ap_uint<512>* meta_in = meta_aggr_in.meta();
    ap_uint<512>* meta_out = meta_aggr_out.meta();

    struct timeval time_now{};
	gettimeofday(&time_now, nullptr);
    time_t start = (time_now.tv_sec * 1000) + (time_now.tv_usec / 1000);

    request request = request_create("com.xilinx.vitis.database.gqeAggr");
    int j = 0;
    for(int i = 0; i < 8; i++){
            request_arg(request, j++, 0, mext_table_in_col[i]);
    }
    request_arg(request, j++, 0, meta_in);
    request_arg(request, j++, 0, meta_out);
    for( int i = 0; i < 16; i++){
            request_arg(request, j++, 0, table_out_col[i]);
    }
    request_arg(request, j++, 0, table_cfg);
    request_arg(request, j++, 0, table_cfg_out);

    session response = inaccel_submit(request);
	if(!response){
		printf("Error at submit\n");
		exit(1);
	}
    printf("%d\n", inaccel_wait(response));
    request_free(request);


    gettimeofday(&time_now, nullptr);
    time_t end = (time_now.tv_sec * 1000) + (time_now.tv_usec / 1000);

    std::cout << "finished data transfer d2h" << std::endl;
    int nrow = meta_aggr_out.getColLen();
    std::cout << "After Aggr Row Num:" << nrow << std::endl;

    std::vector<std::vector<int> > merge_info;
    int output_col_num = aggr_cfg.getOutputColNum();
#ifdef USER_DEBUG
    std::cout << "Merging into " << output_col_num << " cols, info:" << std::endl;
#endif
    for (int i = 0; i < output_col_num; i++) {
        merge_info.push_back(aggr_cfg.getResults(i));
    }
    double l_input_memcpy_size = 0;
    double l_output_memcpy_size = 0;
    for (int i = 0; i < tab_in_col_num; i++) {
        l_input_memcpy_size += (double)tab_in_row_num * tab_in_col_type[i];
    }
    for (int i = 0; i < tab_out_col_num; i++) {
        l_output_memcpy_size += (double)nrow * tab_out_col_type[i];
    }
    l_input_memcpy_size = (double)l_input_memcpy_size / 1024 / 1024;
    l_output_memcpy_size = (double)l_output_memcpy_size / 1024 / 1024;

    std::cout << "-----------------------Data Transfer Info-----------------------" << std::endl;
    std::cout << "H2D size = " << l_input_memcpy_size << " MB" << std::endl;
    std::cout << "D2H size = " << l_output_memcpy_size << " MB" << std::endl;

    std::cout << "------------------------Performance Info------------------------" << std::endl;

    std::cout << std::dec << "Total execution time " << end - start << " ms" << std::endl;

    char** cos_of_table_out = new char*[output_col_num];
    std::cout << "output_col_num:" << output_col_num << std::endl;
    for (int i = 0; i < output_col_num; i++) {
        cos_of_table_out[i] = tab_out.getColPointer(i);
    }
    for (int j = 0; j < nrow; j++) {
        for (int i = 0; i < output_col_num; i++) {
            std::vector<int> index = merge_info[i];
            int col_size = sizeof(int);
            if (index.size() == 1) {
                int tmp = 0;
                memcpy(&tmp, table_out_col[index[0]] + j * col_size, col_size);
                int64_t tmp_64b = tmp;
                memcpy(cos_of_table_out[i] + j * tab_out_col_type[i], &tmp_64b, tab_out_col_type[j]);
            } else if (index.size() == 2) {
                ap_uint<32> low_bits = 0;
                ap_uint<32> high_bits = 0;
                memcpy(&low_bits, table_out_col[index[0]] + j * col_size, col_size);
                memcpy(&high_bits, table_out_col[index[1]] + j * col_size, col_size);
                uint64_t merge_result = (ap_uint<64>)(high_bits, low_bits);
                memcpy(cos_of_table_out[i] + j * tab_out_col_type[i], &merge_result, tab_out_col_type[j]);
            } else {
                std::cout << "Error:Invalid Result index" << std::endl;
                exit(1);
            }
        }
    }
    tab_out.setRowNum(nrow);
    std::cout << "Aggr done, table saved: " << tab_out.getRowNum() << " rows," << tab_out.getColNum() << " cols"
              << std::endl;

    return SUCCESS;
}
/*
struct queue_struct {
    // the sec index
    int sec;
    // the partition index
    int p;
    // the nrow setup of MetaTable, only the last round nrow is different to per_slice_nrow in probe
    int meta_nrow;
    // updating meta info (nrow) for each partition&slice, due to async, this change is done in threads
    MetaTable* meta;
    // dependency event num
    int num_event_wait_list;
    // dependency events
    cl_event* event_wait_list;
    // user event to trace current memcpy operation
    cl_event* event;
    // memcpy src locations
    char* ptr_src[16];
    // ----- part o memcpy in used -----
    // data size of memcpy in
    int size;
    // memcpy dst locations
    char* ptr_dst[16];
    // ----- part o memcpy out used -----
    int partition_num;
    // the allocated size (nrow) of each partititon out buffer
    int part_max_nrow_512;
    // memcpy dst locations, used in part memcpy out
    char*** part_ptr_dst;
    // ----- probe memcpy used -----
    int slice;
    // the nrow of first (slice_num - 1) rounds, only valid in probe memcpy in
    int per_slice_nrow;
    int* tab_col_type_size;
    int* tab_part_sec_nrow;
    // for contiguous cols, input
    //
    //------------------------add by changg----------------------//
    int valid_col_num;
    int key_num;
    int pld_num;
    std::vector<bool> write_flag;
    std::vector<std::vector<int> > merge_info;
};

class threading_pool_for_aggr_pip {
   public:
    std::thread part_l_in_ping_t;
    std::thread part_l_in_pong_t;
    std::thread part_l_out_ping_t;
    std::thread part_l_out_pong_t;

    std::thread aggr_in_ping_t;
    std::thread aggr_in_pong_t;
    std::thread aggr_out_ping_t;
    std::thread aggr_out_pong_t;

    std::queue<queue_struct> q1_ping; // aggr memcpy in used queue
    std::queue<queue_struct> q1_pong; // aggr memcpy in used queue

    std::queue<queue_struct> q2_ping; // aggr memcpy out queue
    std::queue<queue_struct> q2_pong; // aggr memcpy out queue

    // the flag indicate each thread is running
    std::atomic<bool> q1_ping_run;
    std::atomic<bool> q1_pong_run;
    std::atomic<bool> q2_ping_run;
    std::atomic<bool> q2_pong_run;
    //------------------------add by changg----------------------//
    std::unordered_map<Key, Payloads, KeyHasher> ping_merge_map;
    std::unordered_map<Key, Payloads, KeyHasher> pong_merge_map;

    // the total aggr num
    std::atomic<int64_t> aggr_sum_nrow;

    // constructor
    threading_pool_for_aggr_pip() { aggr_sum_nrow = 0; };

    void aggr_memcpy_in_ping_t() {
        while (q1_ping_run) {
#ifdef Valgrind_debug
            sleep(1);
#endif
            while (!q1_ping.empty()) {
                queue_struct q = q1_ping.front();
                clWaitForEvents(q.num_event_wait_list, q.event_wait_list);

                for (int i = 0; i < q.valid_col_num; i++) {
                    memcpy(q.ptr_dst[i], q.ptr_src[i], q.size);
                }

                q.meta->setColNum(q.valid_col_num);
                for (int i = 0; i < q.valid_col_num; i++) {
                    q.meta->setCol(i, i, q.meta_nrow);
                }
                q.meta->meta();

                clSetUserEventStatus(q.event[0], CL_COMPLETE);
                // remove the first element after processing it.
                q1_ping.pop();
            }
        }
    };

    void aggr_memcpy_in_pong_t() {
        while (q1_pong_run) {
#ifdef Valgrind_debug
            sleep(1);
#endif
            while (!q1_pong.empty()) {
                queue_struct q = q1_pong.front();
                clWaitForEvents(q.num_event_wait_list, q.event_wait_list);

                for (int i = 0; i < q.valid_col_num; i++) {
                    memcpy(q.ptr_dst[i], q.ptr_src[i], q.size);
                }

                q.meta->setColNum(q.valid_col_num);
                for (int i = 0; i < q.valid_col_num; i++) {
                    q.meta->setCol(i, i, q.meta_nrow);
                }
                q.meta->meta();

                clSetUserEventStatus(q.event[0], CL_COMPLETE);
                // remove the first element after processing it.
                q1_pong.pop();
            }
        }
    };

    // aggr memcpy out thread
    // int nrow = q.meta->getColLen();
    // for (int c = 0; c < 16; c++) {
    //     if (c != 7)
    //         memcpy(tab_aggr_res_col[c] + aggr_sum_nrow, reinterpret_cast<int*>(tab_aggr_out_col[c][kid]),
    //                nrow * sizeof(int));
    // }
    // aggr_sum_nrow += nrow;
    // std::cout << "output nrow[" << p << "] = " << nrow << std::endl;

    void aggr_memcpy_out_ping_t() {
        while (q2_ping_run) {
#ifdef Valgrind_debug
            sleep(1);
#endif
            while (!q2_ping.empty()) {
                queue_struct q = q2_ping.front();
                clWaitForEvents(q.num_event_wait_list, q.event_wait_list);

                int nrow = q.meta->getColLen();
                std::cout << "output nrow[" << q.p << "] = " << nrow << std::endl;
                char* tab_aggr_res_col[16];
                for (int c = 0; c < 16; c++) {
                    if (q.write_flag[c]) tab_aggr_res_col[c] = q.ptr_src[c];
                }
                for (int i = 0; i < nrow; i++) {
                    Key key;
                    key.key_num = q.key_num;
                    for (int k = 0; k < key.key_num; k++) {
                        std::vector<int> index = q.merge_info[k];
                        if (index.size() > 1) {
                            std::cout << "please ensure grp keys in first cols" << std::endl;
                            exit(1);
                        }
                        // memcpy(&key.keys[k], tab_aggr_res_col[index[0]] + i * q.tab_col_type_size[k],
                        // q.tab_col_type_size[k]);
                        uint32_t tmp = 0;
                        memcpy(&(tmp), tab_aggr_res_col[index[0]] + i * sizeof(int), sizeof(int));
                        key.keys[k] = tmp;
                    }
                    if (ping_merge_map.find(key) != ping_merge_map.end()) {
                        for (int p = 0; p < q.pld_num; p++) {
                            // int col_size = q.tab_col_type_size[p - key.key_num];
                            int col_size = sizeof(int);
                            std::vector<int> index = q.merge_info[p + q.key_num];
                            if (index.size() == 1) {
                                uint32_t tmp = 0;
                                memcpy(&tmp, tab_aggr_res_col[index[0]] + i * col_size, col_size);
                                ping_merge_map[key].values[p] += tmp;
                            } else if (index.size() == 2 || index.size() == 3) {
                                ap_uint<32> low_bits = 0;
                                ap_uint<32> high_bits = 0;
                                memcpy(&high_bits, tab_aggr_res_col[index[1]] + i * col_size, col_size);
                                memcpy(&low_bits, tab_aggr_res_col[index[0]] + i * col_size, col_size);
                                ping_merge_map[key].values[p] += (ap_uint<64>)(high_bits, low_bits);
                            }
                        }
                    } else {
                        Payloads pld;
                        for (int p = 0; p < q.pld_num; p++) {
                            // int col_size = q.tab_col_type_size[p - key.key_num];
                            int col_size = sizeof(int);
                            std::vector<int> index = q.merge_info[p + q.key_num];
                            if (index.size() == 1) {
                                uint32_t tmp = 0;
                                memcpy(&tmp, tab_aggr_res_col[index[0]] + i * col_size, col_size);
                                pld.values[p] = tmp;
                            } else if (index.size() == 2 || index.size() == 3) {
                                ap_uint<32> low_bits = 0;
                                ap_uint<32> high_bits = 0;
                                memcpy(&high_bits, tab_aggr_res_col[index[1]] + i * col_size, col_size);
                                memcpy(&low_bits, tab_aggr_res_col[index[0]] + i * col_size, col_size);
                                pld.values[p] = (ap_uint<64>)(high_bits, low_bits);
                            }
                        }
                        ping_merge_map.insert(std::make_pair(key, pld));
                    }
                }
                aggr_sum_nrow += nrow;
                std::cout << "output aggr_sum_nrow[" << q.p << "] = " << aggr_sum_nrow << std::endl;

                clSetUserEventStatus(q.event[0], CL_COMPLETE);
                // remove the first element after processing it.
                q2_ping.pop();
            }
        }
    }

    void aggr_memcpy_out_pong_t() {
        while (q2_pong_run) {
#ifdef Valgrind_debug
            sleep(1);
#endif
            while (!q2_pong.empty()) {
                queue_struct q = q2_pong.front();
                clWaitForEvents(q.num_event_wait_list, q.event_wait_list);

                int nrow = q.meta->getColLen();
                std::cout << "output nrow[" << q.p << "] = " << nrow << std::endl;
                char* tab_aggr_res_col[16];
                for (int c = 0; c < 16; c++) {
                    if (q.write_flag[c]) tab_aggr_res_col[c] = q.ptr_src[c];
                }
                // std::cout << "&&&&&&&&&&" << q.key_num << std::endl;
                for (int i = 0; i < nrow; i++) {
                    Key key;
                    key.key_num = q.key_num;
                    for (int k = 0; k < key.key_num; k++) {
                        std::vector<int> index = q.merge_info[k];
                        if (index.size() > 1) {
                            std::cout << "please ensure grp keys in first cols" << std::endl;
                            exit(1);
                        }
                        // memcpy(&key.keys[k], tab_aggr_res_col[index[0]] + i * q.tab_col_type_size[k],
                        //       q.tab_col_type_size[k]);
                        uint32_t tmp = 0;
                        memcpy(&tmp, tab_aggr_res_col[index[0]] + i * sizeof(int), sizeof(int));
                        key.keys[k] = tmp;
                    }
                    if (pong_merge_map.find(key) != pong_merge_map.end()) {
                        for (int p = 0; p < q.pld_num; p++) {
                            std::vector<int> index = q.merge_info[p + q.key_num];
                            // int col_size = q.tab_col_type_size[p - key.key_num];
                            int col_size = sizeof(int);
                            if (index.size() == 1) {
                                uint32_t tmp = 0;
                                memcpy(&tmp, tab_aggr_res_col[index[0]] + i * col_size, col_size);
                                pong_merge_map[key].values[p] += tmp;
                            } else if (index.size() == 2 || index.size() == 3) {
                                ap_uint<32> low_bits;
                                ap_uint<32> high_bits;
                                memcpy(&high_bits, tab_aggr_res_col[index[1]] + i * col_size, col_size);
                                memcpy(&low_bits, tab_aggr_res_col[index[0]] + i * col_size, col_size);
                                pong_merge_map[key].values[p] += (ap_uint<64>)(high_bits, low_bits);
                            }
                        }
                    } else {
                        Payloads pld;
                        for (int p = 0; p < q.pld_num; p++) {
                            std::vector<int> index = q.merge_info[p + q.key_num];
                            // int col_size = q.tab_col_type_size[p - key.key_num];
                            int col_size = sizeof(int);
                            if (index.size() == 1) {
                                uint32_t tmp = 0;
                                memcpy(&tmp, tab_aggr_res_col[index[0]] + i * col_size, col_size);
                                pld.values[p] = tmp;
                            } else if (index.size() == 2 || index.size() == 3) {
                                ap_uint<32> low_bits = 0;
                                ap_uint<32> high_bits = 0;
                                memcpy(&high_bits, tab_aggr_res_col[index[1]] + i * col_size, col_size);
                                memcpy(&low_bits, tab_aggr_res_col[index[0]] + i * col_size, col_size);
                                pld.values[p] = (ap_uint<64>)(high_bits, low_bits);
                            }
                        }
                        pong_merge_map.insert(std::make_pair(key, pld));
                    }
                }
                aggr_sum_nrow += nrow;
                std::cout << "output aggr_sum_nrow[" << q.p << "] = " << aggr_sum_nrow << std::endl;

                clSetUserEventStatus(q.event[0], CL_COMPLETE);
                // remove the first element after processing it.
                q2_pong.pop();
            }
        }
    }

    // initialize the table L aggr threads
    void aggr_init() {
        aggr_sum_nrow = 0;
        // start the part o memcpy in thread and non-stop running
        q1_ping_run = 1;
        aggr_in_ping_t = std::thread(&threading_pool_for_aggr_pip::aggr_memcpy_in_ping_t, this);

        // start the part o memcpy in thread and non-stop running
        q1_pong_run = 1;
        aggr_in_pong_t = std::thread(&threading_pool_for_aggr_pip::aggr_memcpy_in_pong_t, this);

        // start the part o memcpy in thread and non-stop running
        q2_ping_run = 1;
        aggr_out_ping_t = std::thread(&threading_pool_for_aggr_pip::aggr_memcpy_out_ping_t, this);

        // start the part o memcpy in thread and non-stop running
        q2_pong_run = 1;
        aggr_out_pong_t = std::thread(&threading_pool_for_aggr_pip::aggr_memcpy_out_pong_t, this);
    }
};
*/
ErrCode Aggregator::aggr_sol1(Table& tab_in, Table& tab_out, AggrConfig& aggr_cfg, std::vector<size_t> params) {
/*
    const int VEC_LEN = 16;
    const int size_of_apu512 = sizeof(ap_uint<512>);

    const int sec_num = params[2];
    int l_nrow = tab_in.getRowNum();
    int l_ncol = tab_in.getColNum();
#ifdef USER_DEBUG
    std::cout << "tab_in info: " << l_nrow << " rows, " << l_ncol << " cols." << std::endl;
#endif

    // start threading pool threads
    gqe::utils::MM mm;
    threading_pool_for_aggr_pip pool;
    pool.aggr_init();

    // aggr kernel
    cl_kernel aggrkernel[2];
    aggrkernel[0] = clCreateKernel(prg, "gqeAggr", &err);
    aggrkernel[1] = clCreateKernel(prg, "gqeAggr", &err);
    // divide table L into many sections
    // the col nrow of each section
    int tab_sec_nrow_each = (l_nrow + sec_num - 1) / sec_num;
    if (tab_sec_nrow_each < sec_num) {
        std::cerr << "Error: Table L section size is smaller than section number!!!";
        std::cerr << "sec size of Table L: " << tab_sec_nrow_each << ", ";
        std::cerr << "sec number of Table L: " << sec_num << std::endl;
        exit(1);
    }
    int* tab_in_col_size = new int[l_ncol];
    int* tab_sec_size = new int[l_ncol];
    for (int i = 0; i < l_ncol; i++) {
        tab_in_col_size[i] = tab_in.getColTypeSize(i);
        tab_sec_size[i] = tab_sec_nrow_each * tab_in_col_size[i];
    }

    int tab_sec_nrow[sec_num];
    for (int sec = 0; sec < sec_num - 1; sec++) {
        tab_sec_nrow[sec] = tab_sec_nrow_each;
    }
    tab_sec_nrow[sec_num - 1] = l_nrow - tab_sec_nrow_each * (sec_num - 1);

    // define and load lineitem table
    // data load from disk. due to table size, data read into several sections
    std::vector<int8_t> scan_list = aggr_cfg.getScanList();
    char* tab_in_user_col_sec[8][sec_num];
    for (int i = 0; i < 8; i++) {
        if (i < l_ncol) {
            for (int j = 0; j < sec_num; j++) {
                tab_in_user_col_sec[i][j] = tab_in.getColPointer(scan_list[i], sec_num, j);
            }
        } else {
            for (int j = 0; j < sec_num; j++) {
                tab_in_user_col_sec[i][j] = mm.aligned_alloc<char>(VEC_LEN);
            }
        }
    }

    // L host side pinned buffers for aggr kernel
    char* tab_in_col[8][2];
    for (int i = 0; i < 8; i++) {
        if (i < l_ncol) {
            tab_in_col[i][0] = mm.aligned_alloc<char>(tab_sec_size[i]);
            tab_in_col[i][1] = mm.aligned_alloc<char>(tab_sec_size[i]);
            memset(tab_in_col[i][0], 0, tab_sec_size[i]);
            memset(tab_in_col[i][1], 0, tab_sec_size[i]);
        } else {
            tab_in_col[i][0] = mm.aligned_alloc<char>(VEC_LEN);
            tab_in_col[i][1] = mm.aligned_alloc<char>(VEC_LEN);
        }
    }

    // define the nrow of aggr result
    // int aggr_result_nrow = tab_sec_nrow_each;
    int aggr_result_nrow = tab_out.getRowNum();
    int out_ncol = aggr_cfg.getOutputColNum();
    int* tab_out_col_type = new int[out_ncol];
    for (int i = 0; i < out_ncol; i++) {
        tab_out_col_type[i] = tab_out.getColTypeSize(i);
    }
    int key_num = aggr_cfg.getGrpKeyNum();
    std::vector<std::vector<int> > merge_info;
    for (int i = 0; i < out_ncol; i++) {
        merge_info.push_back(aggr_cfg.getResults(i));
    }
    std::vector<bool> write_flag = aggr_cfg.getWriteFlag();

#ifdef USER_DEBUG
    for (auto l : write_flag) {
        std::cout << l << " ";
    }
    for (int kk = 0; kk < write_flag.size(); kk++) {
        std::cout << write_flag[kk] << " ";
    }
    std::cout << std::endl;
#endif
    size_t aggr_result_nrow_512 = (aggr_result_nrow + VEC_LEN - 1) / VEC_LEN;
    size_t aggr_result_nrow_512_size = aggr_result_nrow_512 * size_of_apu512;

    char* tab_aggr_out_col[16][2];
    for (int k = 0; k < 2; k++) {
        for (int i = 0; i < 16; i++) {
            if (write_flag[i]) {
                tab_aggr_out_col[i][k] = mm.aligned_alloc<char>(aggr_result_nrow_512_size);
                memset(tab_aggr_out_col[i][k], 0, aggr_result_nrow_512_size);
            } else {
                tab_aggr_out_col[i][k] = mm.aligned_alloc<char>(VEC_LEN);
            }
        }
    }
    ap_uint<32>* cfg_aggr = aggr_cfg.getAggrConfigBits();
    ap_uint<32>* cfg_aggr_out = aggr_cfg.getAggrConfigOutBits();

    //--------------- meta setup -----------------
    // setup meta input and output
    MetaTable meta_aggr_in[2];
    for (int i = 0; i < 2; i++) {
        meta_aggr_in[i].setColNum(l_ncol);
        for (int j = 0; j < l_ncol; j++) {
            meta_aggr_in[i].setCol(j, j, tab_sec_nrow_each);
        }
        meta_aggr_in[i].meta();
    }
    MetaTable meta_aggr_out[2];
    meta_aggr_out[0].setColNum(16);
    meta_aggr_out[1].setColNum(16);
    for (int c = 0; c < 16; c++) {
        meta_aggr_out[0].setCol(c, c, aggr_result_nrow_512);
        meta_aggr_out[1].setCol(c, c, aggr_result_nrow_512);
    }
    meta_aggr_out[0].meta();
    meta_aggr_out[1].meta();

    cl_mem_ext_ptr_t mext_tab_aggr_in_col[8][2];
    cl_mem_ext_ptr_t mext_meta_aggr_in[2], mext_meta_aggr_out[2];
    cl_mem_ext_ptr_t mext_cfg_aggr, mext_cfg_aggr_out;
    cl_mem_ext_ptr_t mext_tab_aggr_out[16][2];
    cl_mem_ext_ptr_t mext_aggr_tmp[8];

    int agg_i = 0;
    for (int k = 0; k < 2; k++) {
        agg_i = 0;
        for (int c = 0; c < 8; c++) {
            mext_tab_aggr_in_col[c][k] = {agg_i++, tab_in_col[c][k], aggrkernel[k]};
        }
        mext_meta_aggr_in[k] = {agg_i++, meta_aggr_in[k].meta(), aggrkernel[k]};
        mext_meta_aggr_out[k] = {agg_i++, meta_aggr_out[k].meta(), aggrkernel[k]};
    }

    for (int k = 0; k < 2; k++) {
        agg_i = 10;
        for (int i = 0; i < 16; ++i) {
            mext_tab_aggr_out[i][k] = {agg_i++, tab_aggr_out_col[i][k], aggrkernel[k]};
        }
    }

    mext_cfg_aggr = {agg_i++, cfg_aggr, aggrkernel[0]};
    mext_cfg_aggr_out = {agg_i++, cfg_aggr_out, aggrkernel[0]};

    for (int c = 0; c < 8; c++) {
        mext_aggr_tmp[c] = {agg_i++, nullptr, aggrkernel[0]};
    }

    cl_mem buf_tab_aggr_in_col[8][2];
    for (int k = 0; k < 2; k++) {
        for (int i = 0; i < 8; ++i) {
            if (i < l_ncol) {
                buf_tab_aggr_in_col[i][k] =
                    clCreateBuffer(ctx, CL_MEM_EXT_PTR_XILINX | CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY, tab_sec_size[i],
                                   &mext_tab_aggr_in_col[i][k], &err);
            } else {
                buf_tab_aggr_in_col[i][k] =
                    clCreateBuffer(ctx, CL_MEM_EXT_PTR_XILINX | CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY, VEC_LEN,
                                   &mext_tab_aggr_in_col[i][k], &err);
            }
        }
    }

    cl_mem buf_aggr_tmp[8];
    for (int i = 0; i < 8; i++) {
        buf_aggr_tmp[i] = clCreateBuffer(ctx, CL_MEM_READ_WRITE | CL_MEM_HOST_NO_ACCESS | CL_MEM_EXT_PTR_XILINX,
                                         (size_t)(8 * S_BUFF_DEPTH), &mext_aggr_tmp[i], &err);
    }

    cl_mem buf_meta_aggr_in[2];
    cl_mem buf_meta_aggr_out[2];
    for (int i = 0; i < 2; i++) {
        buf_meta_aggr_in[i] = clCreateBuffer(ctx, CL_MEM_EXT_PTR_XILINX | CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY,
                                             (sizeof(ap_uint<512>) * 24), &mext_meta_aggr_in[i], &err);
        buf_meta_aggr_out[i] = clCreateBuffer(ctx, CL_MEM_EXT_PTR_XILINX | CL_MEM_USE_HOST_PTR | CL_MEM_READ_WRITE,
                                              (sizeof(ap_uint<512>) * 24), &mext_meta_aggr_out[i], &err);
    }

    cl_mem buf_cfg_aggr = clCreateBuffer(ctx, CL_MEM_EXT_PTR_XILINX | CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY,
                                         size_t(4 * 128), &mext_cfg_aggr, &err);
    cl_mem buf_cfg_aggr_out = clCreateBuffer(ctx, CL_MEM_EXT_PTR_XILINX | CL_MEM_USE_HOST_PTR | CL_MEM_READ_WRITE,
                                             size_t(4 * 128), &mext_cfg_aggr_out, &err);
    cl_mem buf_tab_aggr_out_col[16][2];
    for (int k = 0; k < 2; k++) {
        for (int i = 0; i < 16; ++i) {
            if (write_flag[i])
                buf_tab_aggr_out_col[i][k] =
                    clCreateBuffer(ctx, CL_MEM_EXT_PTR_XILINX | CL_MEM_USE_HOST_PTR | CL_MEM_READ_WRITE,
                                   aggr_result_nrow_512_size, &mext_tab_aggr_out[i][k], &err);
            else
                buf_tab_aggr_out_col[i][k] =
                    clCreateBuffer(ctx, CL_MEM_EXT_PTR_XILINX | CL_MEM_USE_HOST_PTR | CL_MEM_READ_WRITE, VEC_LEN,
                                   &mext_tab_aggr_out[i][k], &err);
        }
    }

    // set args and enqueue kernel
    for (int k = 0; k < 2; k++) {
        int j = 0;
        clSetKernelArg(aggrkernel[k], j++, sizeof(cl_mem), &buf_tab_aggr_in_col[0][k]);
        clSetKernelArg(aggrkernel[k], j++, sizeof(cl_mem), &buf_tab_aggr_in_col[1][k]);
        clSetKernelArg(aggrkernel[k], j++, sizeof(cl_mem), &buf_tab_aggr_in_col[2][k]);
        clSetKernelArg(aggrkernel[k], j++, sizeof(cl_mem), &buf_tab_aggr_in_col[3][k]);
        clSetKernelArg(aggrkernel[k], j++, sizeof(cl_mem), &buf_tab_aggr_in_col[4][k]);
        clSetKernelArg(aggrkernel[k], j++, sizeof(cl_mem), &buf_tab_aggr_in_col[5][k]);
        clSetKernelArg(aggrkernel[k], j++, sizeof(cl_mem), &buf_tab_aggr_in_col[6][k]);
        clSetKernelArg(aggrkernel[k], j++, sizeof(cl_mem), &buf_tab_aggr_in_col[7][k]);
        clSetKernelArg(aggrkernel[k], j++, sizeof(cl_mem), &buf_meta_aggr_in[k]);
        clSetKernelArg(aggrkernel[k], j++, sizeof(cl_mem), &buf_meta_aggr_out[k]);
        for (int c = 0; c < 16; c++) {
            clSetKernelArg(aggrkernel[k], j++, sizeof(cl_mem), &buf_tab_aggr_out_col[c][k]);
        }
        clSetKernelArg(aggrkernel[k], j++, sizeof(cl_mem), &buf_cfg_aggr);
        clSetKernelArg(aggrkernel[k], j++, sizeof(cl_mem), &buf_cfg_aggr_out);
        for (int c = 0; c < 8; c++) {
            clSetKernelArg(aggrkernel[k], j++, sizeof(cl_mem), &buf_aggr_tmp[c]);
        }
    }

    clEnqueueMigrateMemObjects(cq, 1, &buf_meta_aggr_out[0], 0, 0, nullptr, nullptr);
    clEnqueueMigrateMemObjects(cq, 1, &buf_meta_aggr_out[1], 0, 0, nullptr, nullptr);

    std::vector<cl_mem> aggr_in_vec[2];
    for (int k = 0; k < 2; k++) {
        for (int i = 0; i < l_ncol; i++) {
            aggr_in_vec[k].push_back(buf_tab_aggr_in_col[i][k]);
        }
        aggr_in_vec[k].push_back(buf_meta_aggr_in[k]);
        aggr_in_vec[k].push_back(buf_cfg_aggr);
    }

    std::vector<cl_mem> aggr_out_vec[2];
    for (int k = 0; k < 2; k++) {
        for (int i = 0; i < 16; ++i) {
            aggr_out_vec[k].push_back(buf_tab_aggr_out_col[i][k]);
        }
        aggr_out_vec[k].push_back(buf_meta_aggr_out[k]);
    }
    clEnqueueMigrateMemObjects(cq, aggr_in_vec[0].size(), aggr_in_vec[0].data(),
                               CL_MIGRATE_MEM_OBJECT_CONTENT_UNDEFINED, 0, nullptr, nullptr);
    clEnqueueMigrateMemObjects(cq, aggr_in_vec[1].size(), aggr_in_vec[1].data(),
                               CL_MIGRATE_MEM_OBJECT_CONTENT_UNDEFINED, 0, nullptr, nullptr);
    clEnqueueMigrateMemObjects(cq, aggr_out_vec[0].size(), aggr_out_vec[0].data(),
                               CL_MIGRATE_MEM_OBJECT_CONTENT_UNDEFINED, 0, nullptr, nullptr);
    clEnqueueMigrateMemObjects(cq, aggr_out_vec[1].size(), aggr_out_vec[1].data(),
                               CL_MIGRATE_MEM_OBJECT_CONTENT_UNDEFINED, 0, nullptr, nullptr);

    std::vector<std::vector<cl_event> > evt_aggr_memcpy_in(sec_num);
    std::vector<std::vector<cl_event> > evt_aggr_h2d(sec_num);
    std::vector<std::vector<cl_event> > evt_aggr_krn(sec_num);
    std::vector<std::vector<cl_event> > evt_aggr_d2h(sec_num);
    std::vector<std::vector<cl_event> > evt_aggr_memcpy_out(sec_num);
    for (int p = 0; p < sec_num; p++) {
        evt_aggr_memcpy_in[p].resize(1);
        evt_aggr_h2d[p].resize(1);
        evt_aggr_krn[p].resize(1);
        evt_aggr_d2h[p].resize(1);
        evt_aggr_memcpy_out[p].resize(1);

        evt_aggr_memcpy_in[p][0] = clCreateUserEvent(ctx, &err);
        evt_aggr_memcpy_out[p][0] = clCreateUserEvent(ctx, &err);
    }

    std::vector<std::vector<cl_event> > evt_aggr_h2d_dep(sec_num);
    evt_aggr_h2d_dep[0].resize(1);
    for (int i = 1; i < sec_num; ++i) {
        if (i == 1)
            evt_aggr_h2d_dep[i].resize(1);
        else
            evt_aggr_h2d_dep[i].resize(2);
    }
    std::vector<std::vector<cl_event> > evt_aggr_krn_dep(sec_num);
    evt_aggr_krn_dep[0].resize(1);
    for (int i = 1; i < sec_num; ++i) {
        if (i == 1)
            evt_aggr_krn_dep[i].resize(2);
        else
            evt_aggr_krn_dep[i].resize(3);
    }
    std::vector<std::vector<cl_event> > evt_aggr_d2h_dep(sec_num);
    evt_aggr_d2h_dep[0].resize(1);
    for (int i = 1; i < sec_num; ++i) {
        if (i == 1)
            evt_aggr_d2h_dep[i].resize(1);
        else
            evt_aggr_d2h_dep[i].resize(2);
    }

    queue_struct aggr_min[sec_num];
    queue_struct aggr_mout[sec_num];

    // because pure device buf is used in aggr kernel, the kernel needs to be run invalid for 1 round
    {
        std::cout << "xxxxxxxxxxxxxxxxxxxxxx invalid below xxxxxxxxxxxxxxxxxxxxxxxx" << std::endl;
        meta_aggr_in[0].setColNum(1);
        meta_aggr_in[0].setCol(0, 0, 15);
        meta_aggr_in[0].meta();
        clEnqueueMigrateMemObjects(cq, aggr_in_vec[0].size(), aggr_in_vec[0].data(), 0, 0, nullptr,
                                   &evt_aggr_h2d[0][0]);

        clEnqueueTask(cq, aggrkernel[0], evt_aggr_h2d[0].size(), evt_aggr_h2d[0].data(), nullptr);
        clFinish(cq);
        std::cout << "xxxxxxxxxxxxxxxxxxxxxx invalid above xxxxxxxxxxxxxxxxxxxxxxxx" << std::endl;
    }
    gqe::utils::Timer timer_aggr;
    timer_aggr.add(); // 0

    // aggr loop run
    for (int p = 0; p < sec_num; p++) {
        int kid = p % 2;
        // 1)memcpy in, copy the data from tab_part_new_col[partition_num][8][l_new_part_nrow_512] to
        // tab_aggr_in_col0-7
        aggr_min[p].p = p;
        aggr_min[p].valid_col_num = l_ncol;
        aggr_min[p].size = tab_sec_size[p];
        aggr_min[p].event = &evt_aggr_memcpy_in[p][0];
        aggr_min[p].meta = &meta_aggr_in[kid];
        aggr_min[p].meta_nrow = tab_sec_nrow[p];
        for (int i = 0; i < l_ncol; i++) {
            aggr_min[p].ptr_src[i] = tab_in_user_col_sec[i][p];
            aggr_min[p].ptr_dst[i] = tab_in_col[i][kid];
        }
        if (p > 1) {
            aggr_min[p].num_event_wait_list = evt_aggr_h2d[p - 2].size();
            aggr_min[p].event_wait_list = evt_aggr_h2d[p - 2].data();
        } else {
            aggr_min[p].num_event_wait_list = 0;
            aggr_min[p].event_wait_list = nullptr;
        }
        if (kid == 0) pool.q1_ping.push(aggr_min[p]);
        if (kid == 1) pool.q1_pong.push(aggr_min[p]);

        // 2)h2d
        evt_aggr_h2d_dep[p][0] = evt_aggr_memcpy_in[p][0];
        if (p > 1) {
            evt_aggr_h2d_dep[p][1] = evt_aggr_krn[p - 2][0];
        }
        clEnqueueMigrateMemObjects(cq, aggr_in_vec[kid].size(), aggr_in_vec[kid].data(), 0, evt_aggr_h2d_dep[p].size(),
                                   evt_aggr_h2d_dep[p].data(), &evt_aggr_h2d[p][0]);

        // 3)aggr kernel
        evt_aggr_krn_dep[p][0] = evt_aggr_h2d[p][0];
        if (p > 0) {
            evt_aggr_krn_dep[p][1] = evt_aggr_krn[p - 1][0];
        }
        if (p > 1) {
            evt_aggr_krn_dep[p][2] = evt_aggr_d2h[p - 2][0];
        }
        clEnqueueTask(cq, aggrkernel[kid], evt_aggr_krn_dep[p].size(), evt_aggr_krn_dep[p].data(), &evt_aggr_krn[p][0]);
        // 4)d2h
        evt_aggr_d2h_dep[p][0] = evt_aggr_krn[p][0];
        if (p > 1) {
            evt_aggr_d2h_dep[p][1] = evt_aggr_memcpy_out[p - 2][0];
        }
        clEnqueueMigrateMemObjects(cq, aggr_out_vec[kid].size(), aggr_out_vec[kid].data(), CL_MIGRATE_MEM_OBJECT_HOST,
                                   evt_aggr_d2h_dep[p].size(), evt_aggr_d2h_dep[p].data(), &evt_aggr_d2h[p][0]);
        // 5)memcpy out
        aggr_mout[p].p = p;
        aggr_mout[p].write_flag = write_flag;
        aggr_mout[p].merge_info = merge_info;
        aggr_mout[p].tab_col_type_size = tab_out_col_type;
        aggr_mout[p].key_num = key_num;
        aggr_mout[p].pld_num = out_ncol - key_num;
        aggr_mout[p].event = &evt_aggr_memcpy_out[p][0];
        aggr_mout[p].meta = &meta_aggr_out[kid];
        for (int c = 0; c < 16; c++) {
            if (write_flag[c]) {
                aggr_mout[p].ptr_src[c] = tab_aggr_out_col[c][kid];
            }
        }
        aggr_mout[p].num_event_wait_list = evt_aggr_d2h[p].size();
        aggr_mout[p].event_wait_list = evt_aggr_d2h[p].data();
        if (kid == 0) pool.q2_ping.push(aggr_mout[p]);
        if (kid == 1) pool.q2_pong.push(aggr_mout[p]);
    }
    clWaitForEvents(evt_aggr_memcpy_out[sec_num - 1].size(), evt_aggr_memcpy_out[sec_num - 1].data());
    if (sec_num > 1) {
        clWaitForEvents(evt_aggr_memcpy_out[sec_num - 2].size(), evt_aggr_memcpy_out[sec_num - 2].data());
    }
    timer_aggr.add(); // 1

#ifdef AGGR_PERF_PROFILE
    cl_ulong start, end;
    long evt_ns;
    for (int p = 0; p < sec_num; p++) {
        double input_memcpy_size = 0;
        for (int ii = 0; ii < l_ncol; ii++) {
            input_memcpy_size += (double)tab_sec_nrow[p] * tab_in_col_type[ii];
        }
        input_memcpy_size = input_memcpy_size / 1024 / 1024;

        clGetEventProfilingInfo(evt_aggr_h2d[p][0], CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, NULL);
        clGetEventProfilingInfo(evt_aggr_h2d[p][0], CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &end, NULL);
        evt_ns = end - start;
        std::cout << "Sec: " << p << ", h2d, size: " << input_memcpy_size << " MB, time: " << double(evt_ns) / 1000000
                  << " ms, throughput: " << input_memcpy_size / 1024 / ((double)evt_ns / 1000000000) << " GB/s "
                  << std::endl;

        clGetEventProfilingInfo(evt_aggr_krn[p][0], CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, NULL);
        clGetEventProfilingInfo(evt_aggr_krn[p][0], CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &end, NULL);
        evt_ns = end - start;
        std::cout << "Part: " << p << ", krn, size: " << input_memcpy_size << " MB, time: " << double(evt_ns) / 1000000
                  << " ms, throughput: " << input_memcpy_size / 1024 / ((double)evt_ns / 1000000000) << " GB/s "
                  << std::endl;

        clGetEventProfilingInfo(evt_aggr_d2h[p][0], CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, NULL);
        clGetEventProfilingInfo(evt_aggr_d2h[p][0], CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &end, NULL);
        evt_ns = end - start;
        double output_memcpy_size = (double)aggr_result_nrow_512_size * out_ncol / 1024 / 1024;
        std::cout << "Part: " << p << ", d2h, size: " << output_memcpy_size << " MB, time: " << double(evt_ns) / 1000000
                  << " ms, throughput: " << output_memcpy_size / 1024 / ((double)evt_ns / 1000000000) << " GB/s "
                  << std::endl;
    }

    std::cout << "output aggr_sum_nrow = " << pool.aggr_sum_nrow << std::endl;
    std::cout << "ping merge map size = " << pool.ping_merge_map.size() << std::endl;
    std::cout << "pong merge map size = " << pool.pong_merge_map.size() << std::endl;
    std::cout << "-----------------------------------------------------" << std::endl;
#endif
    double tvtime_aggr = timer_aggr.getMilliSec();

    pool.q1_ping_run = 0;
    pool.q1_pong_run = 0;
    pool.q2_ping_run = 0;
    pool.q2_pong_run = 0;

    pool.aggr_in_ping_t.join();
    pool.aggr_in_pong_t.join();
    pool.aggr_out_ping_t.join();
    pool.aggr_out_pong_t.join();
    // when merging the application (0: low bit, 1: high bit, 0~1 sum, 2: count)
    std::unordered_map<Key, int, KeyHasher> key_map;
    std::vector<bool> route_list(merge_info.size());
    std::vector<std::vector<int> > route_info(merge_info.size());
    // merge_info size = output size,
    // route_list get which key needs combination
    // key_map save all the sum locations, using (high,low)->index, using it to find the sum location in
    for (size_t i = 0; i < merge_info.size(); i++) {
        std::vector<int> index = merge_info[i];
        if (index.size() <= 2) {
            route_list[i] = false;
            Key key;
            key.key_num = index.size();
            for (size_t ii = 0; ii < index.size(); ii++) {
                key.keys[ii] = index[ii];
            }
            key_map.insert(std::make_pair(key, i));
        } else {
            route_list[i] = true;
        }
    }
    for (size_t i = 0; i < merge_info.size(); i++) {
        std::vector<int> index = merge_info[i];
        if (index.size() == 3) {
            Key key_sum;
            key_sum.key_num = 2;
            for (int ii = 0; ii < 2; ii++) {
                key_sum.keys[ii] = index[ii];
            }
            Key key_counter;
            key_counter.key_num = 1;
            key_counter.keys[0] = index[2];
            int sum_ind = i;
            int sum_counter = i;
            if (key_map.find(key_sum) != key_map.end()) {
                sum_ind = key_map[key_sum];
            }
            if (key_map.find(key_counter) != key_map.end()) {
                sum_counter = key_map[key_counter];
            }
            route_info[i] = {sum_ind, sum_counter};
        }
    }
    gqe::utils::Timer timer_merge;
    timer_merge.add(); // 0
    for (auto it = pool.pong_merge_map.begin(); it != pool.pong_merge_map.end(); it++) {
        Key key = it->first;
        Payloads pld = it->second;
        if (pool.ping_merge_map.find(key) != pool.ping_merge_map.end()) {
            for (int i = 0; i < out_ncol - key_num; i++) {
                pool.ping_merge_map[key].values[i] += pld.values[i];
            }
        } else {
            for (int i = 0; i < out_ncol - key_num; i++) {
                pool.ping_merge_map[key].values[i] = pld.values[i];
            }
            pool.ping_merge_map.insert(std::make_pair(key, pld));
        }
    }

    char** output_ptr = new char*[out_ncol];
    for (int i = 0; i < out_ncol; i++) {
        output_ptr[i] = tab_out.getColPointer(i);
    }

    tab_out.setRowNum(pool.ping_merge_map.size());
    int32_t row_counter = 0;
    for (auto it = pool.ping_merge_map.begin(); it != pool.ping_merge_map.end(); it++) {
        Key key = it->first;
        for (size_t i = 0; i < merge_info.size(); i++) {
            if (route_list[i] == true) {
                int sum_ind = route_info[i][0];
                int sum_counter = route_info[i][1];
                pool.ping_merge_map[key].values[i - key.key_num] =
                    pool.ping_merge_map[key].values[sum_ind - key.key_num] /
                    pool.ping_merge_map[key].values[sum_counter - key.key_num];
            }
            if (i < (size_t)key.key_num) {
                // std::cout << key.keys[i] << std::endl;
                int64_t key_v = key.keys[i];
                memcpy(output_ptr[i] + row_counter * tab_out_col_type[i], &key_v, tab_out_col_type[i]);

            } else {
                memcpy(output_ptr[i] + row_counter * tab_out_col_type[i - key.key_num],
                       &(pool.ping_merge_map[key].values[i - key.key_num]), tab_out_col_type[i - key.key_num]);
            }
        }
        row_counter++;
    }
    timer_merge.add(); // 1
    std::cout << merge_info.size() << " cols," << row_counter << " rows" << std::endl;
    double in1_bytes = (double)l_nrow * sizeof(int) * l_ncol / 1024 / 1024;
    double out_bytes = (double)pool.ping_merge_map.size() * sizeof(int) * out_ncol / 1024 / 1024;
    std::cout << "-----------------------Data Transfer Info-----------------------" << std::endl;
    std::cout << "H2D size = " << in1_bytes << " MB" << std::endl;
    std::cout << "D2H size = " << out_bytes << " MB" << std::endl;
    std::cout << "-----------------------Performance Info-----------------------" << std::endl;
    double tvtime_merge = timer_merge.getMilliSec();
    std::cout << "All time: " << (double)(tvtime_aggr + tvtime_merge) / 1000
              << " ms, throughput: " << in1_bytes / 1024 / ((double)(tvtime_aggr + tvtime_merge) / 1000000) << " GB/s"
              << std::endl;
    std::cout << "Output number = " << pool.ping_merge_map.size() << std::endl;
    std::cout << "Aggr done, table saved: " << tab_out.getRowNum() << " rows," << tab_out.getColNum() << " cols"
              << std::endl;
*/
    return SUCCESS;
}

inline void zipInt64(void* dest, const void* low, const void* high, size_t n) {
    int64_t* dptr = (int64_t*)dest;
    if (high == nullptr) {
        for (size_t i = 0; i < n; ++i) {
            int64_t tmp = *((const int32_t*)low + i); // sign-ext 32b to 64b
            *(dptr + i) = tmp;
        }
    } else {
        // zip two 32b to one 64b
        for (size_t i = 0; i < n; ++i) {
            uint64_t tmp = *((const uint32_t*)low + i); // zero-ext 32b to 64b
            tmp |= (uint64_t)(*((const uint32_t*)high + i)) << 32;
            *(dptr + i) = (int64_t)tmp; // cast unsigned to signed
        }
    }
}

ErrCode Aggregator::aggr_sol2(Table& tab_in, Table& tab_out, AggrConfig& aggr_cfg, std::vector<size_t> params) {
    const int VEC_LEN = 16;
    const int size_of_apu512 = sizeof(ap_uint<512>);
    gqe::utils::MM mm;
    int l_nrow = tab_in.getRowNum();
    int l_ncol = tab_in.getColNum();

    const int tab_part_sec_num = params[1];
    int log_partition_num = params[3];

#ifdef USER_DEBUG
    std::cout << "tab_in info: " << l_nrow << " rows, " << l_ncol << " cols." << std::endl;
#endif
    // ------------------------------------------
    // --------- partitioning Table L ----------
    // partition setups
    int k_depth = 512;
    const int partition_num = 1 << log_partition_num;

    // divide table L into many sections
    // the col nrow of each section
    int tab_part_sec_nrow_each = (l_nrow + tab_part_sec_num - 1) / tab_part_sec_num;
    if (tab_part_sec_nrow_each < tab_part_sec_num) {
        std::cerr << "Error: Table L section size is smaller than section number!!!";
        std::cerr << "sec size of Table L: " << tab_part_sec_nrow_each << ", ";
        std::cerr << "sec number of Table L: " << tab_part_sec_num << std::endl;
        exit(1);
    }

    int tab_part_sec_nrow[tab_part_sec_num];
    for (int sec = 0; sec < tab_part_sec_num - 1; sec++) {
        tab_part_sec_nrow[sec] = tab_part_sec_nrow_each;
    }
    tab_part_sec_nrow[tab_part_sec_num - 1] = l_nrow - tab_part_sec_nrow_each * (tab_part_sec_num - 1);

    int* tab_in_col_size = new int[l_ncol];
    int* tab_part_sec_size = new int[l_ncol];
    for (int i = 0; i < l_ncol; i++) {
        tab_in_col_size[i] = tab_in.getColTypeSize(i);
        tab_part_sec_size[i] = tab_part_sec_nrow_each * tab_in_col_size[i];
    }

    // todo: change config
    std::vector<int8_t> scan_list = aggr_cfg.getScanList();
    std::vector<int8_t> part_list = aggr_cfg.getPartList();

    char* tab_part_in_col[8][tab_part_sec_num];
    for (int i = 0; i < 8; i++) {
		for (int j = 0; j < tab_part_sec_num; j++) {
			if (i < l_ncol) {
				tab_part_in_col[i][j] = (char*) cube_alloc(tab_part_sec_size[i]);
				if(!tab_part_in_col[i][j]){
					std::cout<< "Error at alloc" << std::endl;
					exit(1);
				}

				memcpy(tab_part_in_col[i][j], tab_in.getColPointer(i, tab_part_sec_num, j), tab_part_sec_nrow[j] * tab_in_col_size[i]);
	        } else {
				tab_part_in_col[i][j] = (char*) cube_alloc(VEC_LEN);
				if(!tab_part_in_col[i][j]){
					std::cout<< "Error at alloc" << std::endl;
					exit(1);
				}

	        }
		}
    }

    // partition output data
    int tab_part_out_col_nrow_512_init = tab_part_sec_nrow_each * 4 / VEC_LEN;
    assert(tab_part_out_col_nrow_512_init > 0 && "Error: table output col size must be greater than 0");
    // the depth of each partition in each col
    int tab_part_out_col_eachpart_nrow_512 = (tab_part_out_col_nrow_512_init + partition_num - 1) / partition_num;
    int l_partition_out_col_part_nrow_max = tab_part_out_col_eachpart_nrow_512 * 16;
    // update depth to make sure the buffer size is aligned by parttion_num * tab_part_out_col_eachpart_nrow_512
    int tab_part_out_col_nrow_512 = partition_num * tab_part_out_col_eachpart_nrow_512;
    int tab_part_out_col_size = tab_part_out_col_nrow_512 * size_of_apu512;

	char* tab_part_out_col[8][tab_part_sec_num];
	for (int i = 0; i < 8; i++) {
		for (int j = 0; j < tab_part_sec_num; j++) {
			if (i < l_ncol) {
				tab_part_out_col[i][j] = (char*) cube_alloc(tab_part_out_col_size);
				if(!tab_part_out_col[i][j]){
					std::cout<< "Error at alloc" << std::endl;
					exit(1);
				}

			} else {
				tab_part_out_col[i][j] = (char*) cube_alloc(VEC_LEN);
				if(!tab_part_out_col[i][j]){
                                        std::cout<< "Error at alloc" << std::endl;
                                        exit(1);
                                }
			}
		}
	}

    ap_uint<512>* cfg_part = aggr_cfg.getPartConfigBits();

    //--------------- metabuffer setup L -----------------
    MetaTable meta_part_in[tab_part_sec_num];

    // setup partition kernel used meta output
    MetaTable meta_part_out[tab_part_sec_num];
    for (int k = 0; k < tab_part_sec_num; k++) {
        meta_part_out[k].setColNum(l_ncol);
        meta_part_out[k].setPartition(partition_num, tab_part_out_col_eachpart_nrow_512);
        for (int j = 0; j < l_ncol; j++) {
            meta_part_out[k].setCol(j, j, tab_part_out_col_nrow_512);
        }
        meta_part_out[k].meta();
    }


	// combine sec0_partition0, sec1_parttion0, ...secN_partition0 in 1 buffer. the depth is
    // int l_new_part_nrow_512 = tab_part_out_col_eachpart_nrow_512 * tab_part_sec_num;
    int l_new_part_nrow_512_size = tab_part_out_col_eachpart_nrow_512 * tab_part_sec_num * size_of_apu512;

    // create user partition res cols
    // all sections partition 0 output to same 8-col buffers
    // ap_uint<512> tab_part_new_col[partition_num][8][l_new_part_nrow_512]
    char* tab_part_new_col[partition_num][8];
    for (int p = 0; p < partition_num; ++p) {
        for (int i = 0; i < 8; ++i) {
            tab_part_new_col[p][i] = (char*) cube_alloc(l_new_part_nrow_512_size);
	    if(!tab_part_new_col[p][i]){
                std::cout<< "Error at alloc" << std::endl;
                exit(1);
    	    }
        }
    }

    std::cout << "=================" << std::endl;

        std::atomic<int> l_new_part_offset[256];
        memset(l_new_part_offset, 0, sizeof(int) * 256);

	session part_response[tab_part_sec_num];
	request part_request[tab_part_sec_num];
	int idx = 0;
    for (int sec = 0; sec < tab_part_sec_num; sec++) {
		meta_part_in[sec].setColNum(l_ncol);
        for (int i = 0; i < l_ncol; i++) {
            meta_part_in[sec].setCol(i, i, tab_part_sec_nrow[sec]);
        }

        // run kernel
		part_request[sec] = request_create("com.xilinx.vitis.database.gqePart");
		int j = 0;
		request_arg(part_request[sec], j++, sizeof(int), &k_depth);
		request_arg(part_request[sec], j++, sizeof(int), &idx);
		request_arg(part_request[sec], j++, sizeof(int), &log_partition_num);

        for (int vv = 0; vv < l_ncol; vv++) {
            request_arg(part_request[sec], j++, 0, tab_part_in_col[part_list[vv]][sec]);
        }
        for (int vv = l_ncol; vv < 8; vv++) {
            request_arg(part_request[sec], j++, 0, tab_part_in_col[vv][sec]);
        }
        request_arg(part_request[sec], j++, 0, meta_part_in[sec].meta()); // size 24*sizeof(ap_uint<512>) with cube_alloc
		request_arg(part_request[sec], j++, 0, meta_part_out[sec].meta()); // size 24*sizeof(ap_uint<512>) with cube_alloc
        for (int vv = 0; vv < l_ncol; vv++) {
            request_arg(part_request[sec], j++, 0, tab_part_out_col[part_list[vv]][sec]);
        }
        for (int vv = l_ncol; vv < 8; vv++) {
            request_arg(part_request[sec], j++, 0, tab_part_out_col[vv][sec]);
        }
        request_arg(part_request[sec], j++, 0, cfg_part);

		part_response[sec] = inaccel_submit(part_request[sec]);
		if(!part_response[sec]){
			printf("Error at submit\n");
			exit(1);
		}
	}

	for (int sec = 0; sec < tab_part_sec_num; sec++) {
		if(inaccel_wait(part_response[sec])){
			std::cout << "Error at Part wait " << sec << std::endl;
			exit(1);
		}
        // memcpy out
		int l_partition_out_col_part_depth = tab_part_out_col_eachpart_nrow_512;
		int* nrow_per_part_l = meta_part_out[sec].getPartLen();

		for (int p = 0; p < partition_num; ++p) {
			int sec_partitioned_res_part_nrow = nrow_per_part_l[p];
			if (sec_partitioned_res_part_nrow > l_partition_out_col_part_nrow_max) {
				std::cerr << "partition out nrow: " << sec_partitioned_res_part_nrow
						  << ", buffer size: " << l_partition_out_col_part_nrow_max << std::endl;
				std::cerr << "ERROR: Table L output partition size is smaller than required!" << std::endl;
				exit(1);
			}

			int offset = l_new_part_offset[p];
			l_new_part_offset[p] += sec_partitioned_res_part_nrow;

			for (int i = 0; i < l_ncol; i++) {
				memcpy(tab_part_new_col[p][i] + offset * tab_in_col_size[i],
					   tab_part_out_col[i][sec] + p * l_partition_out_col_part_depth * sizeof(ap_uint<512>),
					   tab_in_col_size[i] * sec_partitioned_res_part_nrow);
			}
		}

		request_free(part_request[sec]);

    }


#ifdef XDEBUG
    // print new_part column data
    for (int p = 0; p < partition_num; p++) {
        int part_nrow = l_new_part_offset[p];
        std::cout << "----------------------------------" << std::endl;
        std::cout << "Tab L, p: " << p << ", nrow = " << part_nrow << std::endl;
        for (int c = 0; c < 7; c++) {
            std::cout << "col: " << c << std::endl;
            int nrow_16 = std::min((part_nrow + 15) / 16, 10);
            for (int ss = 0; ss < nrow_16; ss++) {
                for (int m = 0; m < 16; m++) {
                    std::cout << tab_part_new_col[p][c][ss](32 * m + 31, 32 * m) << ", ";
                }
                std::cout << std::endl;
            }
        }
    }
    std::cout << "-------------------------------------------------------------" << std::endl;
#endif

#ifdef AGGR_PERF_PROFILE
    for (int sec = 0; sec < tab_part_sec_num; sec++) {
        double input_memcpy_size = 0;
        for (int i = 0; i < l_ncol; i++) {
            input_memcpy_size += tab_part_sec_size[i];
        }
        input_memcpy_size = input_memcpy_size / 1024 / 1024;

        std::cout << "Tab L sec: " << sec << ", h2d, size: " << input_memcpy_size
                  << " MB"
                  << std::endl;

        double output_memcpy_size = (double)tab_part_out_col_size * l_ncol / 1024 / 1024;
        std::cout << "Tab L sec: " << sec << ", d2h, size: " << output_memcpy_size
                  << " MB"
                  << std::endl;
    }

    for (int p = 0; p < partition_num; ++p) {
        std::cout << "Tab L, p: " << p << ", nrow = " << l_new_part_offset[p] << std::endl;
    }
#endif


    //=====================================================
    //==================== group aggr =====================
    //=====================================================

    ap_uint<32>* cfg_aggr = aggr_cfg.getAggrConfigBits();
    ap_uint<32>* cfg_aggr_out = aggr_cfg.getAggrConfigOutBits();
    int out_ncol = aggr_cfg.getOutputColNum();
    if (out_ncol < (int)tab_out.getColNum()) {
        std::cout << "Coumn number of Out table should >= " << out_ncol << std::endl;
        exit(1);
    }
    int* tab_out_col_size = new int[out_ncol];
    for (int i = 0; i < out_ncol; i++) {
        tab_out_col_size[i] = tab_out.getColTypeSize(i);
    }
    int key_num = aggr_cfg.getGrpKeyNum();
    std::vector<std::vector<int> > merge_info;
#ifdef USER_DEBUG
    std::cout << "Merging into " << out_ncol << " cols, info:" << std::endl;
#endif
    for (int i = 0; i < out_ncol; i++) {
        merge_info.push_back(aggr_cfg.getResults(i));
    }
    std::vector<bool> write_flag = aggr_cfg.getWriteFlag();

    // aggr input hbuf
    int aggr_in_nrow_max = 0;
    int aggr_in_nrow[partition_num];
    for (int p = 0; p < partition_num; p++) {
        aggr_in_nrow[p] = l_new_part_offset[p];
        aggr_in_nrow_max = std::max(aggr_in_nrow_max, aggr_in_nrow[p]);
    }
    std::cout << "aggr_in_nrow_max = " << aggr_in_nrow_max << std::endl;
    int* aggr_in_nrow_max_size = new int[l_ncol];
    for (int i = 0; i < l_ncol; i++) {
        aggr_in_nrow_max_size[i] = aggr_in_nrow_max * tab_in_col_size[i];
    }

    char* tab_aggr_in_col[8][partition_num];
    for (int k = 0; k < partition_num; k++) {
        for (int c = 0; c < 8; c++) {
            if (c < l_ncol){
				tab_aggr_in_col[c][k] = (char*) tab_part_new_col[k][scan_list[c]];
			}
            else{
				tab_aggr_in_col[c][k] = (char*) cube_alloc(VEC_LEN);
				if(!tab_aggr_in_col[c][k]){
                                        std::cout<< "Error at alloc" << std::endl;
                                        exit(1);
                                }

			}
        }
    }

    // define the nrow of aggr result
    int aggr_result_nrow = aggr_in_nrow_max;
    size_t aggr_result_nrow_512 = (aggr_result_nrow + VEC_LEN - 1) / VEC_LEN;
    size_t aggr_result_nrow_512_size = aggr_result_nrow_512 * sizeof(ap_uint<512>);
    // in int32 imp, suppose all the data bit is 32 bit
    char* tab_aggr_out_col[16][partition_num];
    for (int k = 0; k < partition_num; k++) {
        for (int i = 0; i < 16; i++) {
            if (write_flag[i]) {
                tab_aggr_out_col[i][k] = (char*) cube_alloc(aggr_result_nrow_512_size);
		if(!tab_aggr_out_col[i][k]){
                	std::cout<< "Error at alloc" << std::endl;
        	        exit(1);
		}

            } else {
                tab_aggr_out_col[i][k] = (char*) cube_alloc(VEC_LEN);
		if(!tab_aggr_out_col[i][k]){
                        std::cout<< "Error at alloc" << std::endl;
                        exit(1);
                }
            }
        }
    }

	//--------------- meta setup -----------------
    // setup meta input and output
	MetaTable meta_aggr_in[partition_num+1];
	MetaTable meta_aggr_out[partition_num];
	for (int i = 0; i < partition_num; i++) {
		meta_aggr_out[i].setColNum(16);
		for (int c = 0; c < 16; c++) {
			meta_aggr_out[i].setCol(c, c, aggr_result_nrow_512);
		}
		meta_aggr_out[i].meta();
	}



    // the collection of each partition aggr results
    // in int32 imp, suppose all the data bit is 32 bit
    char* tab_aggr_res_col[16];
    for (int i = 0; i < 16; i++) {
		void* ptr = nullptr;
		if (posix_memalign(&ptr, 4096, aggr_result_nrow_512_size)) throw std::bad_alloc();
		tab_aggr_res_col[i] = reinterpret_cast<char*>(ptr);
    }

    // because pure device buf is used in aggr kernel, the kernel needs to be run invalid for 1 round
    {
        std::cout << "xxxxxxxxxxxxxxxxxxxxxx invalid below xxxxxxxxxxxxxxxxxxxxxxxx" << std::endl;
		char* tmp_in_col[8];
		for(int i = 0; i < l_ncol; i++){
			tmp_in_col[i] = (char*) cube_alloc(aggr_in_nrow_max_size[i]);
			if(!tmp_in_col[i]){
				std::cout<< "Error at alloc" << std::endl;
				exit(1);
			}

		}
		for( int i = l_ncol; i < 8; i++){
			tmp_in_col[i] = (char*) cube_alloc(VEC_LEN);
			if(!tmp_in_col[i]){
                                std::cout<< "Error at alloc" << std::endl;
                                exit(1);
                        }

		}
        meta_aggr_in[partition_num].setColNum(1);
        meta_aggr_in[partition_num].setCol(0, 0, 15);
	meta_aggr_in[partition_num].meta();
		request request = request_create("com.xilinx.vitis.database.gqeAggr");
		int j = 0;
		for(int i = 0; i < l_ncol; i++){
				request_arg(request, j++, 0, tmp_in_col[i]);
		}
		for( int i = l_ncol; i < 8; i++){
			request_arg(request, j++, 0, tmp_in_col[i]);
		}
		request_arg(request, j++, 0, meta_aggr_in[partition_num].meta());
		request_arg(request, j++, 0, meta_aggr_out[0].meta());
		for( int i = 0; i < 16; i++){
			request_arg(request, j++, 0, tab_aggr_out_col[i][0]);
		}
		request_arg(request, j++, 0, cfg_aggr);
		request_arg(request, j++, 0, cfg_aggr_out);

		session response = inaccel_submit(request);
		if(!response){
			printf("Error at submit\n");
			exit(1);
		}
		if(inaccel_wait(response)){
			std::cout << "Error at Aggr wait" << std::endl;
			exit(1);
		}
		request_free(request);
        std::cout << "xxxxxxxxxxxxxxxxxxxxxx invalid above xxxxxxxxxxxxxxxxxxxxxxxx" << std::endl;
    }

	// the total aggr num
	std::atomic<int64_t> aggr_sum_nrow;
	aggr_sum_nrow = 0;
	session aggr_response[partition_num];
	request aggr_request[partition_num];
    // aggr loop run
    for (int p = 0; p < partition_num; p++) {
	meta_aggr_in[p].setColNum(l_ncol);
        for (int i = 0; i < l_ncol; i++) {
            meta_aggr_in[p].setCol(i, i, aggr_in_nrow[p]);
        }
        // run kernel
		aggr_request[p] = request_create("com.xilinx.vitis.database.gqeAggr");
		int j = 0;
		for(int i = 0; i < l_ncol; i++){
			request_arg(aggr_request[p], j++, 0, tab_aggr_in_col[i][p]);
		}
		for(int i = l_ncol; i < 8; i++){
			request_arg(aggr_request[p], j++, 0, tab_aggr_in_col[i][p]);
		}
		request_arg(aggr_request[p], j++, 0, meta_aggr_in[p].meta());
		request_arg(aggr_request[p], j++, 0, meta_aggr_out[p].meta());
		for( int i = 0; i < 16; i++){
			request_arg(aggr_request[p], j++, 0, tab_aggr_out_col[i][p]);
		}
		request_arg(aggr_request[p], j++, 0, cfg_aggr);
		request_arg(aggr_request[p], j++, 0, cfg_aggr_out);

		aggr_response[p] = inaccel_submit(aggr_request[p]);
		if(!aggr_response[p]){
			printf("Error at submit\n");
			exit(1);
		}
	}
	for (int p = 0; p < partition_num; p++) {
		if( inaccel_wait(aggr_response[p])){
			std::cout << "Error at Aggr wait " << p  << std::endl;
			exit(1);
		}
		request_free(aggr_request[p]);

        // 5)memcpy out
		int nrow = meta_aggr_out[p].getColLen();
		std::cout << "output nrow[" << p << "] = " << nrow << std::endl;

		// in int32 imp, suppose all the col are int32, after aggr, when int64, use the
		// tab_out.tab_col_type_size
		for (int c = 0; c < 16; c++) {
			if (write_flag[c])
				memcpy(tab_aggr_res_col[c] + aggr_sum_nrow * sizeof(int), tab_aggr_out_col[c][p], nrow * sizeof(int));
		}
		aggr_sum_nrow += nrow;
		std::cout << "output aggr_sum_nrow[" << p << "] = " << aggr_sum_nrow << std::endl;
    }


#ifdef AGGR_PERF_PROFILE
    for (int p = 0; p < partition_num; p++) {
        double input_memcpy_size = 0;
        for (int pp = 0; pp < l_ncol; pp++) {
            input_memcpy_size += (double)aggr_in_nrow[p] * tab_in_col_size[pp];
        }
        input_memcpy_size = input_memcpy_size / 1024 / 1024;

        std::cout << "Part: " << p << ", h2d, size: " << input_memcpy_size << " MB"
                  << std::endl;

        double output_memcpy_size = 0;
        for (int pp = 0; pp < 16; pp++) {
            if (write_flag[pp]) output_memcpy_size += (double)aggr_result_nrow_512_size;
        }
        output_memcpy_size = output_memcpy_size / 1024 / 1024;
        std::cout << "Part: " << p << ", d2h, size: " << output_memcpy_size << " MB"
                  << std::endl;
    }

#endif

    double l_input_memcpy_size = 0;
    for (int i = 0; i < l_ncol; i++) {
        l_input_memcpy_size += (double)l_nrow * tab_in_col_size[i];
    }
    l_input_memcpy_size = l_input_memcpy_size / 1024 / 1024;
    double l_output_memcpy_size = (double)aggr_sum_nrow * out_ncol * sizeof(int64_t) / 1024 / 1024;
    std::cout << "-----------------------Data Transfer Info-----------------------" << std::endl;
    std::cout << "H2D size = " << l_input_memcpy_size << " MB" << std::endl;
    std::cout << "D2H size = " << l_output_memcpy_size << " MB" << std::endl;

    tab_out.setRowNum(aggr_sum_nrow);
    // tab_out.setWColNum(out_ncol);
    char** cos_of_table_out = new char*[out_ncol];
    for (int i = 0; i < out_ncol; i++) {
        cos_of_table_out[i] = tab_out.getColPointer(i);
    }

    // XXX For gqeAggr v2, the high/low bits are written to different buffers for 64b output.
    // This should be addressed by 64b native kernels.
    for (int i = 0; i < out_ncol; ++i) {
        auto index = merge_info[i];
        if (index.size() == 1) {
            if (tab_out_col_size[i] == 4) {
                memcpy(cos_of_table_out[i], tab_aggr_res_col[index[0]], 4 * aggr_sum_nrow);
            } else if (tab_out_col_size[i] == 8) {
                // sign extend the low
                zipInt64(cos_of_table_out[i], tab_aggr_res_col[index[0]], nullptr, aggr_sum_nrow);
            } else {
                assert(0 && "Illegal column size");
            }
        } else if (index.size() == 2) {
            if (tab_out_col_size[i] == 4) {
                // truncate and just keep the low bits.
                memcpy(cos_of_table_out[i], tab_aggr_res_col[index[0]], 4 * aggr_sum_nrow);
            } else if (tab_out_col_size[i] == 8) {
                // zip
                zipInt64(cos_of_table_out[i], tab_aggr_res_col[index[0]], tab_aggr_res_col[index[1]],
                         aggr_sum_nrow);
            } else {
                assert(0 && "Illegal column size");
            }
        } else {
            assert(0 && "Unsupported output content handling");
        }
    }

    std::cout << "Aggr done, table saved: " << tab_out.getRowNum() << " rows," << tab_out.getColNum() << " cols"
              << std::endl;
    return SUCCESS;
}

} // database
} // gqe
} // xf
