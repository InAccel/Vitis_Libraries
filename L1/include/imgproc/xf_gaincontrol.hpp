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

#ifndef _XF_GC_HPP_
#define _XF_GC_HPP_

#include "hls_stream.h"
#include "common/xf_common.hpp"

#ifndef XF_IN_STEP
#define XF_IN_STEP 8
#endif
#ifndef XF_OUT_STEP
#define XF_OUT_STEP 8
#endif
/* calculates the weighted sum of 2 inut images */

namespace xf {

namespace cv {
template <int BFORMAT,
          int SRC_T,
          int ROWS,
          int COLS,
          int NPC,
          int PLANES,
          int DEPTH_SRC,
          int DEPTH_DST,
          int WORDWIDTH_SRC,
          int WORDWIDTH_DST,
          int TC>
void gaincontrolkernel(xf::cv::Mat<SRC_T, ROWS, COLS, NPC>& src1,
                       xf::cv::Mat<SRC_T, ROWS, COLS, NPC>& dst,
                       uint16_t height,
                       uint16_t width) {
    ap_uint<13> i, j, k, l;

    int STEP = XF_PIXELWIDTH(SRC_T, NPC) / PLANES;

    XF_SNAME(WORDWIDTH_DST) pxl_pack_out;
    XF_SNAME(WORDWIDTH_SRC) pxl_pack1, pxl_pack2;
RowLoop:
    for (i = 0; i < height; i++) {
#pragma HLS LOOP_TRIPCOUNT min = ROWS max = ROWS
#pragma HLS LOOP_FLATTEN OFF
    ColLoop:
        for (j = 0; j < width; j++) {
#pragma HLS LOOP_TRIPCOUNT min = TC max = TC
#pragma HLS pipeline

            pxl_pack1 = (XF_SNAME(WORDWIDTH_SRC))(src1.read(i * width + j)); // reading from 1st input stream

        ProcLoop:
            for (k = 0, l = 0; k < ((8 << XF_BITSHIFT(NPC)) * PLANES); k += XF_IN_STEP, l += XF_OUT_STEP) {
                XF_PTNAME(DEPTH_SRC) pxl1 = pxl_pack1.range(k + 7, k); // extracting each pixel in case of 8-pixel mode
                XF_PTNAME(DEPTH_SRC) t;
                if ((BFORMAT == XF_BAYER_RG) && (NPC == XF_NPPC2)) {
                    if (i % 2 != 0 && k == 8) {
                        XF_PTNAME(DEPTH_SRC) v1 = pxl1;
                        short v2 = (short)(v1 * 1.34375);
                        if (v2 > 255)
                            t = 255;
                        else
                            t = v2;
                    } else if (i % 2 == 0 && k == 0) {
                        XF_PTNAME(DEPTH_SRC) v1 = pxl1;
                        short v2 = (short)(v1 * 1.94921875);
                        if (v2 > 255)
                            t = 255;
                        else
                            t = v2;
                    } else {
                        t = pxl1;
                    }
                }

                if ((BFORMAT == XF_BAYER_RG) && (NPC == XF_NPPC1)) {
                    if (i % 2 == 0 && j == 0) {
                        XF_PTNAME(DEPTH_SRC) v1 = pxl1;
                        short v2 = (short)(v1 * 1.94921875); //(v1*1.34375);
                        if (v2 > 255)
                            t = 255;
                        else
                            t = v2;
                    } else if (i % 2 != 0 && j != 0) {
                        XF_PTNAME(DEPTH_SRC) v1 = pxl1;
                        short v2 = (short)(v1 * 1.34375);
                        if (v2 > 255)
                            t = 255;
                        else
                            t = v2;
                    } else {
                        t = pxl1;
                    }
                }

                pxl_pack_out.range(l + XF_OUT_STEP - 1, l) = t;
            }

            dst.write(i * width + j, (XF_SNAME(WORDWIDTH_DST))pxl_pack_out); // writing into ouput stream
        }
    }
}

template <int BFORMAT, int SRC_T, int ROWS, int COLS, int NPC = 1>
void gaincontrol(xf::cv::Mat<SRC_T, ROWS, COLS, NPC>& src1, xf::cv::Mat<SRC_T, ROWS, COLS, NPC>& dst) {
#pragma HLS INLINE OFF
#ifndef __SYNTHESIS__
    assert(((src1.rows == dst.rows) && (src1.cols == dst.cols)) && "Input and output image should be of same size");
    assert(((src1.rows <= ROWS) && (src1.cols <= COLS)) && "ROWS and COLS should be greater than input image");
    // assert(((NPC == XF_NPPC1) || (NPC == XF_NPPC8) ) && "NPC must be XF_NPPC1, XF_NPPC8 ");
#endif
    short width = src1.cols >> XF_BITSHIFT(NPC);

    gaincontrolkernel<BFORMAT, SRC_T, ROWS, COLS, NPC, XF_CHANNELS(SRC_T, NPC), XF_DEPTH(SRC_T, NPC),
                      XF_DEPTH(SRC_T, NPC), XF_WORDWIDTH(SRC_T, NPC), XF_WORDWIDTH(SRC_T, NPC),
                      (COLS >> XF_BITSHIFT(NPC))>(src1, dst, src1.rows, width);
}
} // namespace cv
} // namespace xf
#endif //_XF_GC_HPP_
