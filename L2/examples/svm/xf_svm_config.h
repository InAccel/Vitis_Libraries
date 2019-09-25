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

#ifndef _XF_SVM_CONFIG_H_
#define _XF_SVM_CONFIG_H_

/////  To set the parameters in the test-bench /////

#include "ap_int.h"
#include "hls_stream.h"
#include "common/xf_common.hpp"
#include "common/xf_utility.hpp"
#include "xf_config_params.h"
#include "imgproc/xf_svm.hpp"

// Set input and output pixel depth:
#define IN_TYPE XF_16SC1
#define PTR_IN_WIDTH 16
#define OUT_TYPE XF_32SC1
#define PTR_OUT_WIDTH 32

// Set the optimization type:
#define NPC1 XF_NPPC1

#endif // end of _XF_SVM_CONFIG_H_
