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

#include "common/xf_headers.hpp"
#include "xf_integral_image_config.h"

int main(int argc, char** argv) {
    cv::Mat in_img, in_img1, out_img, ocv_ref, ocv_ref1;
    cv::Mat in_gray, in_gray1, diff;

    if (argc != 2) {
        fprintf(stderr, "Usage: <executable> <input image>\n");
        return -1;
    }

    // Read input image
    in_img = cv::imread(argv[1], 0);
    if (in_img.data == NULL) {
        // cout << "Can't open image !!" << endl;
        return -1;
    }

    // create memory for output images
    ocv_ref.create(in_img.rows, in_img.cols, CV_32S);
    ocv_ref1.create(in_img.rows, in_img.cols, CV_32S);
    cv::integral(in_img, ocv_ref, -1);

    for (int i = 0; i < in_img.rows; i++) {
        for (int j = 0; j < in_img.cols; j++) {
            ocv_ref1.at<unsigned int>(i, j) = ocv_ref.at<unsigned int>(i + 1, j + 1);
        }
    }

    imwrite("out_ocv.png", ocv_ref1);

    // create memory for output image
    diff.create(in_img.rows, in_img.cols, CV_32S);
    out_img.create(in_img.rows, in_img.cols, CV_32S);

    // image height and width
    uint16_t height = in_img.rows;
    uint16_t width = in_img.cols;

    static xf::cv::Mat<XF_8UC1, HEIGHT, WIDTH, XF_NPPC1> imgInput(in_img.rows, in_img.cols);
    static xf::cv::Mat<XF_32UC1, HEIGHT, WIDTH, XF_NPPC1> imgOutput(in_img.rows, in_img.cols);

    imgInput.copyTo(in_img.data);

    integral_accel(imgInput, imgOutput);

    // Write output image
    xf::cv::imwrite("hls_out.jpg", imgOutput);

    // Compute absolute difference image
    xf::cv::absDiff(ocv_ref1, imgOutput, diff);

    // Save the difference image
    imwrite("diff.png", diff);

    float err_per;
    double minval = 256, maxval = 0;
    int cnt = 0;
    for (int i = 0; i < in_img.rows; i++) {
        for (int j = 0; j < in_img.cols; j++) {
            unsigned int v = diff.at<unsigned int>(i, j);

            if (v > 0) cnt++;
            if (minval > v) minval = v;
            if (maxval < v) maxval = v;
        }
    }

    err_per = 100.0 * (float)cnt / (in_img.rows * in_img.cols);
    fprintf(stderr,
            "Minimum error in intensity = %f\nMaximum error in intensity = %f\nPercentage of pixels above error "
            "threshold = %f\n",
            minval, maxval, err_per);

    if (err_per > 0.0f) {
        return 1;
    }

    return 0;
}
