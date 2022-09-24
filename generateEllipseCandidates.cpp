//#include "mex.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <climits>
#include <cfloat>
#include <iostream>
//#include "lapack.h"  //matlab 
//#include "include/lapacke_config.h"  //lapack�ֶ���δ�ɹ� (//lapack manual, unsuccessful)
//#include "include/lapacke.h"
#include "/usr/local/Cellar/opencv@2/2.4.13.7_12/include/opencv2/core/core.hpp"
#include "/usr/local/Cellar/opencv@2/2.4.13.7_12/include/opencv2/features2d/features2d.hpp"
#include "/usr/local/Cellar/opencv@2/2.4.13.7_12/include/opencv2/nonfree/features2d.hpp"
#include "/usr/local/Cellar/opencv@2/2.4.13.7_12/include/opencv2/highgui/highgui.hpp"
#include </usr/local/Cellar/opencv@2/2.4.13.7_12/include/opencv2/opencv.hpp>
#include "lapack.h"
///#include "/Applications/MATLAB_R2018a.app/extern/include/mex.hpp"
//#include "/Applications/MATLAB_R2018a.app/extern/include/mexAdapter.hpp"

#include "HelperFunctions.h"
#include "Structs.h"
#include "MagicNumbers.h"
#include "Gradient.h"
#include "Region.h"
#include "LineSegment.h"
#include "NFA.h"
#include "Cluster.h"


using namespace cv;

//=============================================================================
//��Ҫ��������ͷ�ļ� (Need to include the following header files)
//#include <opencv2\opencv.hpp>
//using namespace cv;
void cvCanny3(const void *srcarr, void *dstarr,
              void *dxarr, void *dyarr,
              int aperture_size) {
    //cv::Ptr<CvMat> dx, dy;
    cv::AutoBuffer<char> buffer;
    std::vector<uchar *> stack;
    uchar **stack_top = nullptr, **stack_bottom = nullptr;

    CvMat srcstub, *src = cvGetMat(srcarr, &srcstub);
    CvMat dststub, *dst = cvGetMat(dstarr, &dststub);

    CvMat dxstub, *dx = cvGetMat(dxarr, &dxstub);
    CvMat dystub, *dy = cvGetMat(dyarr, &dystub);


    CvSize size;
    int flags = aperture_size;
    int low, high;
    int *mag_buf[3];
    uchar *map;
    ptrdiff_t mapstep;
    int maxsize;
    int i, j;
    CvMat mag_row;

    if (CV_MAT_TYPE(src->type) != CV_8UC1 ||
        CV_MAT_TYPE(dst->type) != CV_8UC1 ||
        CV_MAT_TYPE(dx->type) != CV_16SC1 ||
        CV_MAT_TYPE(dy->type) != CV_16SC1)
        CV_Error(CV_StsUnsupportedFormat, "");

    if (!CV_ARE_SIZES_EQ(src, dst))
        CV_Error(CV_StsUnmatchedSizes, "");

    aperture_size &= INT_MAX;
    if ((aperture_size & 1) == 0 || aperture_size < 3 || aperture_size > 7)
        CV_Error(CV_StsBadFlag, "");


    size.width = src->cols;
    size.height = src->rows;

    //aperture_size = -1; //SCHARR
    cvSobel(src, dx, 1, 0, aperture_size);
    cvSobel(src, dy, 0, 1, aperture_size);

    Mat1f magGrad(size.height, size.width, 0.f);
    float maxGrad(0);
    float val(0);
    for (i = 0; i < size.height; ++i) {
        float *_pmag = magGrad.ptr<float>(i);
        const short *_dx = (short *) (dx->data.ptr + dx->step * i);
        const short *_dy = (short *) (dy->data.ptr + dy->step * i);
        for (j = 0; j < size.width; ++j) {
            val = float(abs(_dx[j]) + abs(_dy[j]));
            _pmag[j] = val;
            maxGrad = (val > maxGrad) ? val : maxGrad;
        }
    }

    //% Normalize for threshold selection
    //normalize(magGrad, magGrad, 0.0, 1.0, NORM_MINMAX);

    //% Determine Hysteresis Thresholds

    //set magic numbers
    const int NUM_BINS = 64;
    const double percent_of_pixels_not_edges = 0.9;
    const double threshold_ratio = 0.3;

    //compute histogram
    int bin_size = cvFloor(maxGrad / float(NUM_BINS) + 0.5f) + 1;
    if (bin_size < 1) bin_size = 1;
    int bins[NUM_BINS] = {0};
    for (i = 0; i < size.height; ++i) {
        float *_pmag = magGrad.ptr<float>(i);
        for (j = 0; j < size.width; ++j) {
            int hgf = int(_pmag[j]);
            bins[int(_pmag[j]) / bin_size]++;
        }
    }




    //% Select the thresholds
    float total(0.f);
    float target = float(size.height * size.width * percent_of_pixels_not_edges);
    int low_thresh, high_thresh(0);

    while (total < target) {
        total += bins[high_thresh];
        high_thresh++;
    }
    high_thresh *= bin_size;
    low_thresh = cvFloor(threshold_ratio * float(high_thresh));

    if (flags & CV_CANNY_L2_GRADIENT) {
        Cv32suf ul, uh;
        ul.f = (float) low_thresh;
        uh.f = (float) high_thresh;

        low = ul.i;
        high = uh.i;
    } else {
        low = cvFloor(low_thresh);
        high = cvFloor(high_thresh);
    }


    buffer.allocate((size.width + 2) * (size.height + 2) + (size.width + 2) * 3 * sizeof(int));
    mag_buf[0] = (int *) (char *) buffer;
    mag_buf[1] = mag_buf[0] + size.width + 2;
    mag_buf[2] = mag_buf[1] + size.width + 2;
    map = (uchar * )(mag_buf[2] + size.width + 2);
    mapstep = size.width + 2;

    maxsize = MAX(1 << 10, size.width * size.height / 10);
    stack.resize(maxsize);
    stack_top = stack_bottom = &stack[0];

    memset(mag_buf[0], 0, (size.width + 2) * sizeof(int));
    memset(map, 1, mapstep);
    memset(map + mapstep * (size.height + 1), 1, mapstep);

    /* sector numbers
       (Top-Left Origin)

        1   2   3
         *  *  *
          * * *
        0*******0
          * * *
         *  *  *
        3   2   1
    */

#define CANNY_PUSH(d)    *(d) = (uchar)2, *stack_top++ = (d)
#define CANNY_POP(d)     ((d) = *--stack_top)

    mag_row = cvMat(1, size.width, CV_32F);

    // calculate magnitude and angle of gradient, perform non-maxima supression.
    // fill the map with one of the following values:
    //   0 - the pixel might belong to an edge
    //   1 - the pixel can not belong to an edge
    //   2 - the pixel does belong to an edge
    for (i = 0; i <= size.height; i++) {
        int *_mag = mag_buf[(i > 0) + 1] + 1;
        float *_magf = (float *) _mag;
        const short *_dx = (short *) (dx->data.ptr + dx->step * i);
        const short *_dy = (short *) (dy->data.ptr + dy->step * i);
        uchar *_map;
        int x, y;
        ptrdiff_t magstep1, magstep2;
        int prev_flag = 0;

        if (i < size.height) {
            _mag[-1] = _mag[size.width] = 0;

            if (!(flags & CV_CANNY_L2_GRADIENT))
                for (j = 0; j < size.width; j++)
                    _mag[j] = abs(_dx[j]) + abs(_dy[j]);

            else {
                for (j = 0; j < size.width; j++) {
                    x = _dx[j];
                    y = _dy[j];
                    _magf[j] = (float) std::sqrt((double) x * x + (double) y * y);
                }
            }
        } else
            memset(_mag - 1, 0, (size.width + 2) * sizeof(int));

        // at the very beginning we do not have a complete ring
        // buffer of 3 magnitude rows for non-maxima suppression
        if (i == 0)
            continue;

        _map = map + mapstep * i + 1;
        _map[-1] = _map[size.width] = 1;

        _mag = mag_buf[1] + 1; // take the central row
        _dx = (short *) (dx->data.ptr + dx->step * (i - 1));
        _dy = (short *) (dy->data.ptr + dy->step * (i - 1));

        magstep1 = mag_buf[2] - mag_buf[1];
        magstep2 = mag_buf[0] - mag_buf[1];

        if ((stack_top - stack_bottom) + size.width > maxsize) {
            int sz = (int) (stack_top - stack_bottom);
            maxsize = MAX(maxsize * 3 / 2, maxsize + 8);
            stack.resize(maxsize);
            stack_bottom = &stack[0];
            stack_top = stack_bottom + sz;
        }

        for (j = 0; j < size.width; j++) {
#define CANNY_SHIFT 15
#define TG22  (int)(0.4142135623730950488016887242097*(1<<CANNY_SHIFT) + 0.5)

            x = _dx[j];
            y = _dy[j];
            int s = x ^y;
            int m = _mag[j];

            x = abs(x);
            y = abs(y);
            if (m > low) {
                int tg22x = x * TG22;
                int tg67x = tg22x + ((x + x) << CANNY_SHIFT);

                y <<= CANNY_SHIFT;

                if (y < tg22x) {
                    if (m > _mag[j - 1] && m >= _mag[j + 1]) {
                        if (m > high && !prev_flag && _map[j - mapstep] != 2) {
                            CANNY_PUSH(_map + j);
                            prev_flag = 1;
                        } else
                            _map[j] = (uchar) 0;
                        continue;
                    }
                } else if (y > tg67x) {
                    if (m > _mag[j + magstep2] && m >= _mag[j + magstep1]) {
                        if (m > high && !prev_flag && _map[j - mapstep] != 2) {
                            CANNY_PUSH(_map + j);
                            prev_flag = 1;
                        } else
                            _map[j] = (uchar) 0;
                        continue;
                    }
                } else {
                    s = s < 0 ? -1 : 1;
                    if (m > _mag[j + magstep2 - s] && m > _mag[j + magstep1 + s]) {
                        if (m > high && !prev_flag && _map[j - mapstep] != 2) {
                            CANNY_PUSH(_map + j);
                            prev_flag = 1;
                        } else
                            _map[j] = (uchar) 0;
                        continue;
                    }
                }
            }
            prev_flag = 0;
            _map[j] = (uchar) 1;
        }

        // scroll the ring buffer
        _mag = mag_buf[0];
        mag_buf[0] = mag_buf[1];
        mag_buf[1] = mag_buf[2];
        mag_buf[2] = _mag;
    }

    // now track the edges (hysteresis thresholding)
    while (stack_top > stack_bottom) {
        uchar *m;
        if ((stack_top - stack_bottom) + 8 > maxsize) {
            int sz = (int) (stack_top - stack_bottom);
            maxsize = MAX(maxsize * 3 / 2, maxsize + 8);
            stack.resize(maxsize);
            stack_bottom = &stack[0];
            stack_top = stack_bottom + sz;
        }

        CANNY_POP(m);

        if (!m[-1])
            CANNY_PUSH(m - 1);
        if (!m[1])
            CANNY_PUSH(m + 1);
        if (!m[-mapstep - 1])
            CANNY_PUSH(m - mapstep - 1);
        if (!m[-mapstep])
            CANNY_PUSH(m - mapstep);
        if (!m[-mapstep + 1])
            CANNY_PUSH(m - mapstep + 1);
        if (!m[mapstep - 1])
            CANNY_PUSH(m + mapstep - 1);
        if (!m[mapstep])
            CANNY_PUSH(m + mapstep);
        if (!m[mapstep + 1])
            CANNY_PUSH(m + mapstep + 1);
    }

    // the final pass, form the final image
    for (i = 0; i < size.height; i++) {
        const uchar *_map = map + mapstep * (i + 1) + 1;
        uchar *_dst = dst->data.ptr + dst->step * i;

        for (j = 0; j < size.width; j++) {
            _dst[j] = (uchar) - (_map[j] >> 1);
        }
    }
}

void Canny3(InputArray image, OutputArray _edges,
            OutputArray _sobel_x, OutputArray _sobel_y,
            int apertureSize, bool L2gradient) {
    Mat src = image.getMat();
    _edges.create(src.size(), CV_8U);
    _sobel_x.create(src.size(), CV_16S);
    _sobel_y.create(src.size(), CV_16S);


    CvMat c_src = src, c_dst = _edges.getMat();
    CvMat c_dx = _sobel_x.getMat();
    CvMat c_dy = _sobel_y.getMat();


    cvCanny3(&c_src, &c_dst,
             &c_dx, &c_dy,
             apertureSize + (L2gradient ? CV_CANNY_L2_GRADIENT : 0));
}


//=============================================================================
/** Convert ellipse from matrix form to common form:
    ellipse = (centrex,centrey,ax,ay,orientation).
 */
int ellipse2Param(const double *p, double param[]) {
    // ax^2 + bxy + cy^2 + dx + ey + f = 0
    double a, b, c, d, e, f;
    double thetarad, cost, sint, cos_squared, sin_squared, cos_sin, Ao, Au, Av, Auu, Avv, tuCentre, tvCentre, wCentre, uCentre, vCentre, Ru, Rv;
    a = p[0];
    b = p[1];
    c = p[2];
    d = p[3];
    e = p[4];
    f = p[5];

    thetarad = 0.5 * atan2(b, a - c);
    cost = cos(thetarad);
    sint = sin(thetarad);
    sin_squared = sint * sint;
    cos_squared = cost * cost;
    cos_sin = sint * cost;
    Ao = f;
    Au = d * cost + e * sint;
    Av = -d * sint + e * cost;
    Auu = a * cos_squared + c * sin_squared + b * cos_sin;
    Avv = a * sin_squared + c * cos_squared - b * cos_sin;

    if (Auu == 0 || Avv == 0) {
        param[0] = 0;
        param[1] = 0;
        param[2] = 0;
        param[3] = 0;
        param[4] = 0;
        return 0;
    } else {
        tuCentre = -Au / (2. * Auu);
        tvCentre = -Av / (2. * Avv);
        wCentre = Ao - Auu * tuCentre * tuCentre - Avv * tvCentre * tvCentre;
        uCentre = tuCentre * cost - tvCentre * sint;
        vCentre = tuCentre * sint + tvCentre * cost;
        Ru = -wCentre / Auu;
        Rv = -wCentre / Avv;
        //     if (Ru>0) Ru=pow(Ru,0.5);
        //     else Ru=-pow(-Ru,0.5);
        //     if (Rv>0) Rv=pow(Rv,0.5);
        //     else Rv=-pow(-Rv,0.5);
        if (Ru <= 0 || Rv <= 0)//������С��0�����������
            return 0;
        Ru = sqrt(Ru);
        Rv = sqrt(Rv);
        param[0] = uCentre;
        param[1] = vCentre;
        param[2] = Ru;
        param[3] = Rv;
        param[4] = thetarad;
        //�����Ru < Rv������Ե�һ�� (There will be Ru < Rv situation, swap it)
        if (Ru < Rv) {
            param[2] = Rv;
            param[3] = Ru;
            if (thetarad <
                0)//���������ᣬʹ�õ���������Ϊ���ᣬ���ĸ�Ϊ���� (Swap the long and short axes so that the third parameter is the long axis and the fourth is the short axis)
                param[4] += M_1_2_PI;
            else
                param[4] -= M_1_2_PI;
            if (thetarad <
                -M_1_2_PI)//��������޶���-pi/2 ~ pi/2���߱�Ψһ�� (The long-axis inclination is limited to -pi/2 ~ pi/2, which is unique)
                param[4] += M_PI;
            if (thetarad > M_1_2_PI)
                param[4] -= M_PI;
        }
    }
    return 1;
}

//input : (xi,yi)
//output: x0,y0,a,b,phi,ellipara��Ҫ���������ڴ� (x0,y0,a,b,phi,ellipara need to apply for memory in advance)
//successfull, return 1; else return 0
int fitEllipse(point2d *dataxy, int datanum, double *ellipara) {
    auto *D = (double *) malloc(datanum * 6 * sizeof(double));
    double S[36];
    double C[36];
    memset(D, 0, sizeof(double) * datanum);
    memset(S, 0, sizeof(double) * 36);
    memset(C, 0, sizeof(double) * 36);
    for (int i = 0; i < datanum; i++) {
        D[i * 6] = dataxy[i].x * dataxy[i].x;
        D[i * 6 + 1] = dataxy[i].x * dataxy[i].y;
        D[i * 6 + 2] = dataxy[i].y * dataxy[i].y;
        D[i * 6 + 3] = dataxy[i].x;
        D[i * 6 + 4] = dataxy[i].y;
        D[i * 6 + 5] = 1;
    }
    for (int i = 0; i < 6; i++)
        for (int j = i; j < 6; j++) {
            //S[i*6+j]
            for (int k = 0; k < datanum; k++)
                S[i * 6 + j] += D[k * 6 + i] * D[k * 6 + j];
        }
    free(D);//�ͷ��ڴ� (free memory)
    //�Գƾ���ֵ (Symmetric Matrix Assignment)
    for (int i = 0; i < 6; i++)
        for (int j = 0; j < i; j++)
            S[i * 6 + j] = S[j * 6 + i];
    C[0 * 6 + 2] = 2;
    C[1 * 6 + 1] = -1;
    C[2 * 6 + 0] = 2;
    // eig(S,C) eig(inv(S)*C)
    double alphar[6], alphai[6], beta[6];
    double vl[36] = {0};//�˴����� (not used here)
    double vr[36] = {0};
    char JOBVL = 'N';
    char JOBVR = 'V';
    ptrdiff_t fitN = 6;
    double fitWork[64];
    ptrdiff_t workLen = 64;
    ptrdiff_t info;
    //info = LAPACKE_dggev(LAPACK_ROW_MAJOR,'N','V',6,S,6,C,6,alphar,alphai,beta,vl,6,vr,6);
    //ע��SΪ�Գƾ��󣬹�ת�ú���ڱ�����������ȣ�S���Բ���
    //(Note that S is a symmetric matrix, so it is equal to itself after transposition, becoming column priority, and S can be unchanged)
    dggev(&JOBVL, &JOBVR, &fitN, S, &fitN, C, &fitN, alphar, alphai, beta, vl, &fitN, vr, &fitN, fitWork, &workLen,
          &info);
    if (info == 0) {
        int index = -1;
        for (int i = 0; i < 6; i++)
            if ((alphar[i] >= -(2.2204460492503131e-014)) && (alphai[i] == 0) &&
                (beta[i] != 0)) // 100*DBL_EPSILON, eigenvalue = (alphar + i*alphai)/beta
                index = i;//vr[:,i],vr��i�ж�Ӧ������������Ϊ��ϲ��� (vr[:,i], the eigenvector corresponding to the i-th column of vr is the fitting parameter)
        if (index ==
            -1)//����һ�Σ��ſ��ʵ��>0��Լ�����ſ�>-0.005 (try again, relax the constraint on the real part > 0 to >-0.005)
        {
            double temp = -0.005;//��������ܹؼ� (//This parameter is critical)
            for (int i = 0; i < 6; i++)
                if ((alphar[i] >= temp) && (alphai[i] == 0) &&
                    (beta[i] != 0)) // 100*DBL_EPSILON, eigenvalue = (alphar + i*alphai)/beta
                {
                    temp = alphar[i];
                    index = i;//vr[:,i],vr��i�ж�Ӧ������������Ϊ��ϲ��� (vr[:,i], the eigenvector corresponding to the i-th column of vr is the fitting parameter)
                }
        }
        if (index != -1) {
            //�˴�����beta�����ݲ��� (Borrow beta to pass parameters here)
            //beta[0] = vr[6*0+index];
            //beta[1] = vr[6*1+index];
            //beta[2] = vr[6*2+index];
            //beta[3] = vr[6*3+index];
            //beta[4] = vr[6*4+index];
            //beta[5] = vr[6*5+index];
            beta[0] = vr[6 * index + 0];
            beta[1] = vr[6 * index + 1];
            beta[2] = vr[6 * index + 2];
            beta[3] = vr[6 * index + 3];
            beta[4] = vr[6 * index + 4];
            beta[5] = vr[6 * index + 5];
            ellipse2Param(beta, ellipara);//ax^2 + bxy + cy^2 + dx + ey + f = 0, transform to (x0,y0,a,b,phi)
            return 1;
        }
    }
    return 0;
}

//input: dataxyΪ���ݵ�(xi,yi),�ܹ���datanum�� (dataxy is the data point (xi,yi), there are a total of datanum)
//output: ��Ͼ���S. ע�⣺S��Ҫ���������ڴ棬double S[36]. (Fitting matrix S. Note: S needs to apply for memory in advance, double S[36].)
inline void calcuFitMatrix(point2d *dataxy, int datanum, double *S) {
    auto *D = (double *) malloc(datanum * 6 * sizeof(double));
    memset(D, 0, sizeof(double) * datanum);
    for (int i = 0; i < datanum; i++) {
        D[i * 6] = dataxy[i].x * dataxy[i].x;
        D[i * 6 + 1] = dataxy[i].x * dataxy[i].y;
        D[i * 6 + 2] = dataxy[i].y * dataxy[i].y;
        D[i * 6 + 3] = dataxy[i].x;
        D[i * 6 + 4] = dataxy[i].y;
        D[i * 6 + 5] = 1;
    }
    for (int i = 0; i < 6; i++) {
        for (int j = i; j < 6; j++) {
            //S[i*6+j]
            for (int k = 0; k < datanum; k++)
                S[i * 6 + j] += D[k * 6 + i] * D[k * 6 + j];
        }
    }
    free(D);//�ͷ��ڴ� (free memory)
    //�Գƾ���ֵ (Symmetric Matrix Assignment)
    for (int i = 0; i < 6; i++)
        for (int j = 0; j < i; j++)
            S[i * 6 + j] = S[j * 6 + i];
}

//input: fit matrixes S1,S2. length is 36.
//output: fit matrix S_out. S_out = S1 + S2.
//S_out������Ҫ�����ڴ� (S_out needs to apply for memory in advance)
inline void addFitMatrix(const double *S1, const double *S2, double *S_out) {
    int ind;
    for (int i = 0; i < 6; i++)
        for (int j = i; j < 6; j++) {
            ind = i * 6 + j;
            S_out[ind] = S1[ind] + S2[ind];
        }
    //�Գƾ���ֵ (Symmetric Matrix Assignment)
    for (int i = 0; i < 6; i++)
        for (int j = 0; j < i; j++)
            S_out[i * 6 + j] = S_out[j * 6 + i];
}

//input : S����6 x 6 = 36
//output: (A,B,C,D,E,F)��A>0, ellicoeff��Ҫ���������ڴ�. ��Ҫת����(x0,y0,a,b,phi)
//		  ((A,B,C,D,E,F) and A>0, ellicoeff needs to apply for memory in advance. When converting to (x0,y0,a,b,phi))
//ellipse2Param(ellicoeff,ellipara); ax^2 + bxy + cy^2 + dx + ey + f = 0, transform to (x0,y0,a,b,phi)
//successfull, return 1; else return 0
int fitEllipse2(const double *S, double *ellicoeff) {
    double C[36];
    memset(C, 0, sizeof(double) * 36);

    C[0 * 6 + 2] = 2;
    C[1 * 6 + 1] = -1;
    C[2 * 6 + 0] = 2;
    // eig(S,C) eig(inv(S)*C)
    double alphar[6], alphai[6], beta[6];
    double vl[36] = {0};//�˴����� (not used here)
    double vr[36] = {0};
    char JOBVL = 'N';
    char JOBVR = 'V';
    ptrdiff_t fitN = 6;
    double fitWork[64];
    ptrdiff_t workLen = 64;
    ptrdiff_t info;
    //info = LAPACKE_dggev(LAPACK_ROW_MAJOR,'N','V',6,S,6,C,6,alphar,alphai,beta,vl,6,vr,6);
    dggev(&JOBVL, &JOBVR, &fitN, S, &fitN, C, &fitN, alphar, alphai, beta, vl, &fitN, vr, &fitN, fitWork, &workLen,
          &info);
    if (info == 0) {
        int index = -1;
        for (int i = 0; i < 6; i++)
            if ((alphar[i] >= -(2.2204460492503131e-014)) && (alphai[i] == 0) &&
                (beta[i] != 0)) // 100*DBL_EPSILON, eigenvalue = (alphar + i*alphai)/beta
                index = i;//vr[:,i],vr��i�ж�Ӧ������������Ϊ��ϲ��� (vr[:,i], the eigenvector corresponding to the i-th column of vr is the fitting parameter)
        if (index ==
            -1)//����һ�Σ��ſ��ʵ��>0��Լ�����ſ�>-0.005 (try again, relax the constraint on the real part > 0 to >-0.005)
        {
            double temp = -0.005;//��������ܹؼ�
            for (int i = 0; i < 6; i++)
                if ((alphar[i] >= temp) && (alphai[i] == 0) &&
                    (beta[i] != 0)) // 100*DBL_EPSILON, eigenvalue = (alphar + i*alphai)/beta
                {
                    temp = alphar[i];
                    index = i;//vr[:,i],vr��i�ж�Ӧ������������Ϊ��ϲ��� (vr[:,i], the eigenvector corresponding to the i-th column of vr is the fitting parameter)
                }
        }
        if (index != -1) {
            //�˴�����beta�����ݲ��� (Borrow beta to pass parameters here)
            if (vr[6 * index + 0] < 0)//ע�������� (Pay attention to column first)
            {
                ellicoeff[0] = -vr[6 * index + 0]; //-vr[6*0+index];
                ellicoeff[1] = -vr[6 * index + 1]; //-vr[6*1+index];
                ellicoeff[2] = -vr[6 * index + 2]; //-vr[6*2+index];
                ellicoeff[3] = -vr[6 * index + 3]; //-vr[6*3+index];
                ellicoeff[4] = -vr[6 * index + 4]; //-vr[6*4+index];
                ellicoeff[5] = -vr[6 * index + 5]; //-vr[6*5+index];
            } else {
                ellicoeff[0] = vr[6 * index + 0];//vr[6*0+index];
                ellicoeff[1] = vr[6 * index + 1];//vr[6*1+index];
                ellicoeff[2] = vr[6 * index + 2];//vr[6*2+index];
                ellicoeff[3] = vr[6 * index + 3];//vr[6*3+index];
                ellicoeff[4] = vr[6 * index + 4];//vr[6*4+index];
                ellicoeff[5] = vr[6 * index + 5];//vr[6*5+index];
            }
            return 1;
        }
    }
    return 0;
}

//��Σ�e1 = (x1,y1,a1,b1,phi1), e2 = (x2,y2,a2,b2,phi2) (Input parameters: e1 = (x1,y1,a1,b1,phi1), e2 = (x2,y2,a2,b2,phi2))
//��������Ϊ1������Ϊ0 (Output: 1 if equal, 0 otherwise)
inline bool
isEllipseEqual(double *ellipse1, double *ellipse2, double centers_distance_threshold, double semimajor_errorratio,
               double semiminor_errorratio, double angle_errorratio, double iscircle_ratio) {
    bool con1 = (abs(ellipse1[0] - ellipse2[0]) < centers_distance_threshold &&
                 abs(ellipse1[1] - ellipse2[1]) < centers_distance_threshold &&
                 abs(ellipse1[2] - ellipse2[2]) / max(ellipse1[2], ellipse2[2]) < semimajor_errorratio &&
                 abs(ellipse1[3] - ellipse2[3]) / min(ellipse1[3], ellipse2[3]) < semiminor_errorratio);
    bool con2 = (ellipse1[3] / ellipse1[2] >= iscircle_ratio);//0.9 0.85
    bool con3 = (ellipse2[3] / ellipse2[2] >= iscircle_ratio);
    bool con4 = ((con2 && con3) ||
                 (!con2 && !con3 && abs(ellipse1[4] - ellipse2[4]) <= angle_errorratio * M_PI));
    return (con1 && con4);
}

inline bool
regionLimitation(point2d point_g1s, point2d g1s_ls_dir, point2d point_g1e, point2d g1e_ls_dir, point2d point_g2s,
                 point2d g2s_ls_dir, point2d point_g2e, point2d g2e_ls_dir, double polarity,
                 double region_limitation_dis_tolerance) {
    point2d g1m_ls_dir{}, g2m_ls_dir{};
    point2d g1s_arc_dir{}, g1e_arc_dir{}, g1m_arc_dir{}, g2s_arc_dir{}, g2e_arc_dir{}, g2m_arc_dir{};
    point2d test_vec1{}, test_vec2{}, test_vec3{}; //��ָ��Բ�ĵ������Ͳ������� (The vector of the arc pointing to the center of the circle and the test vector)
    //���pend<-pstart���ɵ�����Ϊgim_arc_dir (The vector composed of pend<-pstart of the group is gim_arc_dir)
    double xdelta, ydelta, theta;
    xdelta = point_g1e.x - point_g1s.x;
    ydelta = point_g1e.y - point_g1s.y;
    theta = atan2(ydelta, xdelta);
    g1m_ls_dir.x = cos(theta);
    g1m_ls_dir.y = sin(theta);
    xdelta = point_g2e.x - point_g2s.x;
    ydelta = point_g2e.y - point_g2s.y;
    theta = atan2(ydelta, xdelta);
    g2m_ls_dir.x = cos(theta);
    g2m_ls_dir.y = sin(theta);

    if (polarity == 1)// polarity is equal 1, arc vector = (dy,-dx)
    {
        g1s_arc_dir.x = g1s_ls_dir.y;
        g1s_arc_dir.y = -g1s_ls_dir.x;
        g1e_arc_dir.x = g1e_ls_dir.y;
        g1e_arc_dir.y = -g1e_ls_dir.x;
        g1m_arc_dir.x = g1m_ls_dir.y;
        g1m_arc_dir.y = -g1m_ls_dir.x;
        g2s_arc_dir.x = g2s_ls_dir.y;
        g2s_arc_dir.y = -g2s_ls_dir.x;
        g2e_arc_dir.x = g2e_ls_dir.y;
        g2e_arc_dir.y = -g2e_ls_dir.x;
        g2m_arc_dir.x = g2m_ls_dir.y;
        g2m_arc_dir.y = -g2m_ls_dir.x;
    } else// polarity is equal -1, arc vector = (-dy,dx)
    {
        g1s_arc_dir.x = -g1s_ls_dir.y;
        g1s_arc_dir.y = g1s_ls_dir.x;
        g1e_arc_dir.x = -g1e_ls_dir.y;
        g1e_arc_dir.y = g1e_ls_dir.x;
        g1m_arc_dir.x = -g1m_ls_dir.y;
        g1m_arc_dir.y = g1m_ls_dir.x;
        g2s_arc_dir.x = -g2s_ls_dir.y;
        g2s_arc_dir.y = g2s_ls_dir.x;
        g2e_arc_dir.x = -g2e_ls_dir.y;
        g2e_arc_dir.y = g2e_ls_dir.x;
        g2m_arc_dir.x = -g2m_ls_dir.y;
        g2m_arc_dir.y = g2m_ls_dir.x;
    }
    test_vec1.x = (point_g2e.x - point_g1s.x);
    test_vec1.y = (point_g2e.y - point_g1s.y);
    test_vec2.x = (point_g2s.x - point_g1e.x);
    test_vec2.y = (point_g2s.y - point_g1e.y);
    test_vec3.x = (test_vec1.x + test_vec2.x) / 2;
    test_vec3.y = (test_vec1.y + test_vec2.y) / 2;
    double t1, t2, t3, t4, t5, t6;
    t1 = dotProduct(g1s_arc_dir, test_vec1);
    t2 = dotProduct(g1e_arc_dir, test_vec2);
    t3 = dotProduct(g1m_arc_dir, test_vec3);
    t4 = -dotProduct(g2e_arc_dir, test_vec1);
    t5 = -dotProduct(g2s_arc_dir, test_vec2);
    t6 = -dotProduct(g2m_arc_dir, test_vec3);

    if (dotProduct(g1s_arc_dir, test_vec1) >= region_limitation_dis_tolerance && \
         dotProduct(g1e_arc_dir, test_vec2) >= region_limitation_dis_tolerance && \
         dotProduct(g1m_arc_dir, test_vec3) >= region_limitation_dis_tolerance && \
        -dotProduct(g2e_arc_dir, test_vec1) >= region_limitation_dis_tolerance && \
        -dotProduct(g2s_arc_dir, test_vec2) >= region_limitation_dis_tolerance && \
        -dotProduct(g2m_arc_dir, test_vec3) >= region_limitation_dis_tolerance
            )
        return TRUE;
    return FALSE;
}

/*
void drawEllipse(Mat img, double * ellipara)
{
  Point peliicenter(ellipara[0],ellipara[1]);
  Size  saxis(ellipara[2],ellipara[3]);
  //Mat ellimat = Mat::zeros(img.rows,img.cols,CV_8UC3);
  //ellimat.setTo(255);
  static int ccc = 0;
  static unsigned int cnt = 0;
  if(cnt % 2 == 0 )
	  ccc = 0;
  else
  {
	  ccc = 255;
	  cout<<cnt/2<<'\t'<<ellipara[0]<<'\t'<<ellipara[1]<<"\t"<<ellipara[2]<<'\t'<<ellipara[3]<<'\t'<<ellipara[4]<<endl;
  }
  cnt++;

  Mat imgtemp = img.clone();
  ellipse(imgtemp,peliicenter,saxis,ellipara[4]*180/M_PI,0,360,(Scalar(0,255,ccc)),2);
  namedWindow("w1");
  imshow("w1",imgtemp);
  //waitKey(0);
}
void drawEdge(Mat img, point2d * dataxy, int num)
{
	 static int ccc = 0;
     static int cnt = 0;
     cnt++;
     if(cnt % 2 == 0 )
	     ccc = 0;
     else
	  ccc = 255;
	Mat imgtemp = img.clone();
	for (int i = 0; i<num; i++)
	{
		imgtemp.at<Vec3b>(dataxy[i].y,dataxy[i].x) = (Vec3b(ccc,255,0));
	}
	namedWindow("w2");
    imshow("w2",imgtemp);
}
*/

/*----------------------------------------------------------------------------*/
/** Approximate the distance between a point and an ellipse using Rosin distance.
 */
inline double d_rosin(double *param, double x, double y) {
    double ae2 = param[2] * param[2];
    double be2 = param[3] * param[3];
    x = x - param[0];
    y = y - param[1];
    double xp = x * cos(-param[4]) - y * sin(-param[4]);
    double yp = x * sin(-param[4]) + y * cos(-param[4]);
    double fe2;
    fe2 = ae2 - be2;
    double X = xp * xp;
    double Y = yp * yp;
    double delta = (X + Y + fe2) * (X + Y + fe2) - 4 * X * fe2;
    double A = (X + Y + fe2 - sqrt(delta)) / 2.0;
    double ah = sqrt(A);
    double bh2 = fe2 - A;
    double term = (A * be2 + ae2 * bh2);
    double xi = ah * sqrt(ae2 * (be2 + bh2) / term);
    double yi = param[3] * sqrt(bh2 * (ae2 - A) / term);
    double d[4], dmin;


    d[0] = dist(xp, yp, xi, yi);
    d[1] = dist(xp, yp, xi, -yi);
    d[2] = dist(xp, yp, -xi, yi);
    d[3] = dist(xp, yp, -xi, -yi);
    dmin = DBL_MAX;
    for (double i : d) {
        if (i <= dmin)
            dmin = i;
    }
//  if (X+Y>xi*xi+yi*yi)
//    return dmin;
//  else return -dmin;
    return dmin;
}
/*----------------------------------------------------------------------------*/

//���� (enter)
//lsd�㷨���õ����߶μ���lines������line_num��return�ķ���ֵ��line_nums���߶Σ�Ϊһάdouble������lines������Ϊ8*n��ÿ8��Ϊһ��
// (The number line_num of the line segment set lines detected by the lsd algorithm, the return value of return is line_nums
// line segments, which is a one-dimensional double-type array lines, with a length of 8*n, and each 8 is a group)
//����x1,y1,x2,y2,dx,dy,length,polarity (Store x1,y1,x2,y2,dx,dy,length,polarity)
//groups: �߶η��飬ÿ����水�ռ��ηֲ�˳��˳ʱ�������ʱ��洢���߶��������߶�������Χ��0~line_num-1. ����������ָ�룬ʹ��ʱҪע��(*group)
// (groups: Line segments are grouped. Each group stores the line segment index clockwise or counterclockwise according to
// the geometric distribution order. The range of the line segment index is 0~line_num-1. Since it is a pointer, pay attention when using it (*group))
//first_group_ind��second_group_ind��ƥ����ӵ�����������ȡsalient hypothesisʱ��second_group_ind = -1, fit_matrix2 = NULL.
// (first_group_ind, second_group_ind are the indices of matching teams, when extracting the salient hypothesis, second_group_ind = -1, fit_matrix2 = NULL.)
//fit_matrix1, fit_matrix2, �ֱ�����ӵĶ�Ӧ����Ͼ��� (are the corresponding fitting matrices of the team, respectively)
//angles, �Ǳ�Ե��ͼ+�ݶȷ��� �ޱ�Ե��ʱ��NODEF (is the edge point map + gradient direction. NODEF when there are no edge points)
//distance_tolerance:
//group_inliers_num:��¼�Ÿ������֧���ڵ����������飬ʵʱ���£���ʼʱΪ0
// (group_inliers_num: An array that records the number of inliers supported by each group, updated in real time, initially 0)
//��� (output)
//ellipara
bool calcEllipseParametersAndValidate(double *lines, int line_num, vector <vector<int>> *groups, int first_group_ind,
                                      int second_group_ind, double *fit_matrix1, double *fit_matrix2,
                                      image_double angles, double distance_tolerance, unsigned int *group_inliers_num,
                                      point5d *ellipara) {
    double S[36]; //��Ͼ���S (fit matrix S)
    double Coefficients[6] = {0, 0, 0, 0, 0, 0};// ax^2 + bxy + cy^2 + dx + ey + f = 0
    double param[5], param2[5];
    int info, addr;
    rect rec{};
    rect_iter *iter;
    int rec_support_cnt, rec_inliers_cnt;
    bool flag1 = TRUE, flag2 = TRUE;
    double point_normalx, point_normaly, point_normal, temp;
    std::vector<point2i> first_group_inliers, second_group_inliers;
    point2i pixel_temp{};
    double semimajor_errorratio, semiminor_errorratio, iscircle_ratio;
    if (fit_matrix2 == nullptr ||
        second_group_ind == -1)//ֻ��һ�����ǶȽϴ���������� (Fit only one group with greater coverage)
    {
        for (int i = 0; i < 36; i++)
            S[i] = fit_matrix1[i];
    } else {
        addFitMatrix(fit_matrix1, fit_matrix2, S);//����Խ�����ϣ� S = fit_matrix1 + fit_matrix2 (Fitting pairs of groups)
    }
    info = fitEllipse2(S, Coefficients);// ax^2 + bxy + cy^2 + dx + ey + f = 0, a > 0
    if (info == 0)//���ʧ�� (Fit failed)
    {
        ellipara = nullptr;
        return FALSE;
    }
    ellipse2Param(Coefficients, param);// (x0,y0,a,b,phi)
    if (min(param[2], param[3]) < 3 * distance_tolerance ||
        max(param[2], param[3]) > min(angles->xsize, angles->ysize) || param[0] < 0 || param[0] > angles->xsize ||
        param[1] < 0 || param[1] > angles->ysize) {
        ellipara = nullptr;
        return FALSE;
    }
    //if ( first_group_ind == 2 && second_group_ind == 8)
    //drawEllipse(img,param);
    //����е� first group�Ƚ����ڵ�׼����֤�����Ҹ������֧���ڵ�����
    // (The first group in the team is verified by the interior point criterion first, and the number of supported interior points of the group is updated)
    for (unsigned int i = 0; i < (*groups)[first_group_ind].size(); i++) {
        addr = (*groups)[first_group_ind][i] *
               8; //��first_group_ind����ĵ�i���߶�����*8 (The i-th line segment index of the first_group_ind grouping*8)
        rec.x1 = lines[addr];
        rec.y1 = lines[addr + 1];
        rec.x2 = lines[addr + 2];
        rec.y2 = lines[addr + 3];
        rec.x = (rec.x1 + rec.x2) / 2;
        rec.y = (rec.y1 + rec.y2) / 2;
        rec.dx = lines[addr + 4];
        rec.dy = lines[addr + 5];
        rec.width = 3 * distance_tolerance;
        //line_length[i] = (int)lines[addr+6];//��¼�߶γ��ȵ�����line_length[i] (Record the length of the line segment to the array line_length[i])
        rec_support_cnt = rec_inliers_cnt = 0;//�������Ҫ (clearing is important)
        if (lines[addr + 7] == 1) //����һ�� (same polarity)
        {
            for (iter = ri_ini(&rec); !ri_end(iter); ri_inc(iter))//�߶�1 (line segment 1)
            {
                //��Ӿ��ο��ܻ�Խ�� (The bounding rectangle may be out of bounds)
                if (iter->x >= 0 && iter->y >= 0 && iter->x < angles->xsize && iter->y < angles->ysize) {
                    temp = angles->data[iter->y * angles->xsize +
                                        iter->x];//�ڵ���ݶȷ��� (Gradient direction of interior points)
                    if (temp != NOTDEF) {
                        //test point's normal is (ax0+by0/2+d/2, cy0+bx0/2+e/2)
                        point_normalx = Coefficients[0] * iter->x + (Coefficients[1] * iter->y + Coefficients[3]) / 2;
                        point_normaly = Coefficients[2] * iter->y + (Coefficients[1] * iter->x + Coefficients[4]) / 2;
                        point_normal = atan2(-point_normaly,
                                             -point_normalx); //��Ե��ķ��߷���,ָ����Բ�ڲ� (The normal direction of the edge point, pointing to the inside of the ellipse)
                        rec_inliers_cnt++;
                        if (angle_diff(point_normal, temp) <= M_1_8_PI) //+- 22.5���� �� || d - r || < 3 dis_t
                        {
                            rec_support_cnt++;
                            pixel_temp.x = iter->x;
                            pixel_temp.y = iter->y;
                            first_group_inliers.push_back(
                                    pixel_temp);//��Ӹ��߶ζ�Ӧ���ڵ� (Add the interior point corresponding to the line segment)
                        }
                    }
                }
            }
        } else//�����෴ (opposite polarity)
        {
            for (iter = ri_ini(&rec); !ri_end(iter); ri_inc(iter))//�߶�1 (line segment 1)
            {
                //��Ӿ��ο��ܻ�Խ�� (The bounding rectangle may be out of bounds)
                if (iter->x >= 0 && iter->y >= 0 && iter->x < angles->xsize && iter->y < angles->ysize) {
                    temp = angles->data[iter->y * angles->xsize +
                                        iter->x];//�ڵ���ݶȷ��� (Gradient direction of interior points)
                    if (temp != NOTDEF) {
                        //test point's normal is (ax0+by0/2+d/2, cy0+bx0/2+e/2)
                        point_normalx = Coefficients[0] * iter->x + (Coefficients[1] * iter->y + Coefficients[3]) / 2;
                        point_normaly = Coefficients[2] * iter->y + (Coefficients[1] * iter->x + Coefficients[4]) / 2;
                        point_normal = atan2(point_normaly,
                                             point_normalx); //��Ե��ķ��߷���,ָ����Բ��� (The normal direction of the edge point, pointing to the outside of the ellipse)
                        rec_inliers_cnt++;
                        if (angle_diff(point_normal, temp) <=
                            M_1_8_PI) //+- 22.5���� �� || d - r || < 3 dis_t (+- within 22.5�� and || d - r || < 3 dis_t)
                        {
                            rec_support_cnt++;
                            pixel_temp.x = iter->x;
                            pixel_temp.y = iter->y;
                            first_group_inliers.push_back(
                                    pixel_temp);//��Ӹ��߶ζ�Ӧ���ڵ� (Add the interior point corresponding to the line segment)
                        }
                    }
                }
            }
        }
        if (!(rec_support_cnt > 0 &&
              (rec_support_cnt >= 0.8 * lines[addr + 6] || rec_support_cnt * 1.0 / rec_inliers_cnt >= 0.6))) {
            flag1 = FALSE; //flag1 ��ʼ��ʱΪTRUE, һ��������һ���߶β�����Ҫ��ֱ��false, �ڵ�׼����֤��ͨ��
            //(flag1 is TRUE when initialized, once there is a line segment in the group that does not meet the requirements, it is directly false, and the interior point criterion verification fails)
            break;
        }
    }
    if (flag1 == TRUE && first_group_inliers.size() >= 0.8 *
                                                       group_inliers_num[first_group_ind])//�������ͳ�ƹ����ڵ�,ͨ����֤ (Close to the largest statistical interior point, pass the verification)
    {
        if (first_group_inliers.size() >=
            group_inliers_num[first_group_ind])//��������ֹ�������ڵ��� (The maximum number of inliers the update group has ever seen)
            group_inliers_num[first_group_ind] = first_group_inliers.size();
    } else
        flag1 = FALSE;
    //��һ���������֤ (The first group completes verification)
    if (second_group_ind == -1 ||
        fit_matrix2 == nullptr)//ֻ��һ�����ǶȽϴ���������� (Fit only one group with greater coverage)
    {
        ellipara->x = param[0];//��Ϊ������Σ�����Ҫ����������ǿ����Բ (Because in any case, it is necessary to return a highly significant ellipse)
        ellipara->y = param[1];
        ellipara->a = param[2];
        ellipara->b = param[3];
        ellipara->phi = param[4];
        if (flag1 == TRUE)//ͨ���ڵ��ٴ���ϣ�������� (Refit with interior points to improve quality)
        {
            auto *dataxy = (point2d *) malloc(sizeof(point2d) * first_group_inliers.size());
            for (unsigned int i = 0; i < first_group_inliers.size(); i++) {
                dataxy[i].x = first_group_inliers[i].x;
                dataxy[i].y = first_group_inliers[i].y;
            }
            info = fitEllipse(dataxy, first_group_inliers.size(), param2);
            free(dataxy); //�ͷ��ڴ� (free memory)
            if (info == 1 && isEllipseEqual(param2, param, 3 * distance_tolerance, 0.1, 0.1, 0.1, 0.9)) {
                ellipara->x = param2[0];//������Բ�����Ʒ�� (Update ellipse, improve quality)
                ellipara->y = param2[1];
                ellipara->a = param2[2];
                ellipara->b = param2[3];
                ellipara->phi = param2[4];
                //drawEllipse(img,param2);
            }
        }
        return TRUE;//����ֻ��һ�������ȡ��Բ����ʱֱ�ӷ��� (For the extraction ellipse with only one group, return directly at this time)
    }
    //��������������е� second group�����ڵ�׼����֤�����Ҹ������֧���ڵ�����
    // (Next, verify the interior point criterion for the second group in the team, and update the number of supported interior points for the group)
    if (flag1 ==
        FALSE)//����������У������һ���鶼�޷������ڵ�Ҫ��ֱ�ӷ���false (In the team operation, if the first group cannot meet the interior point requirements, return false directly)
        return FALSE;
    for (unsigned int i = 0; i < (*groups)[second_group_ind].size(); i++) {
        addr = (*groups)[second_group_ind][i] *
               8; //��first_group_ind����ĵ�i���߶�����*8 (The i-th line segment index of the first_group_ind grouping*8)
        rec.x1 = lines[addr];
        rec.y1 = lines[addr + 1];
        rec.x2 = lines[addr + 2];
        rec.y2 = lines[addr + 3];
        rec.x = (rec.x1 + rec.x2) / 2;
        rec.y = (rec.y1 + rec.y2) / 2;
        rec.dx = lines[addr + 4];
        rec.dy = lines[addr + 5];
        rec.width = 3 * distance_tolerance;
        //line_length[i] = (int)lines[addr+6];//��¼�߶γ��ȵ�����line_length[i] (Record the length of the line segment to the array line_length[i])
        rec_support_cnt = rec_inliers_cnt = 0;//�������Ҫ (clearing is important)
        if (lines[addr + 7] == 1) //����һ�� (same polarity)
        {
            for (iter = ri_ini(&rec); !ri_end(iter); ri_inc(iter))//�߶�1 (�߶�1)
            {
                //��Ӿ��ο��ܻ�Խ�� (The bounding rectangle may be out of bounds)
                if (iter->x >= 0 && iter->y >= 0 && iter->x < angles->xsize && iter->y < angles->ysize) {
                    temp = angles->data[iter->y * angles->xsize +
                                        iter->x];//�ڵ���ݶȷ��� (Gradient direction of interior points)
                    if (temp != NOTDEF) {
                        //test point's normal is (ax0+by0/2+d/2, cy0+bx0/2+e/2)
                        point_normalx = Coefficients[0] * iter->x + (Coefficients[1] * iter->y + Coefficients[3]) / 2;
                        point_normaly = Coefficients[2] * iter->y + (Coefficients[1] * iter->x + Coefficients[4]) / 2;
                        point_normal = atan2(-point_normaly,
                                             -point_normalx); //��Ե��ķ��߷���,ָ����Բ�ڲ� (The normal direction of the edge point, pointing to the inside of the ellipse)
                        rec_inliers_cnt++;
                        if (angle_diff(point_normal, temp) <=
                            M_1_8_PI) //+- 22.5���� �� || d - r || < 3 dis_t (+- within 22.5�� and || d - r || < 3 dis_t)
                        {
                            rec_support_cnt++;
                            pixel_temp.x = iter->x;
                            pixel_temp.y = iter->y;
                            second_group_inliers.push_back(
                                    pixel_temp);//��Ӹ��߶ζ�Ӧ���ڵ� (Add the interior point corresponding to the line segment)
                        }
                    }
                }
            }
        } else//�����෴ (opposite polarity)
        {
            for (iter = ri_ini(&rec); !ri_end(iter); ri_inc(iter))//�߶�1 (line segment 1)
            {
                //��Ӿ��ο��ܻ�Խ�� (The bounding rectangle may be out of bounds)
                if (iter->x >= 0 && iter->y >= 0 && iter->x < angles->xsize && iter->y < angles->ysize) {
                    temp = angles->data[iter->y * angles->xsize +
                                        iter->x];//�ڵ���ݶȷ��� (Gradient direction of interior points)
                    if (temp != NOTDEF) {
                        //test point's normal is (ax0+by0/2+d/2, cy0+bx0/2+e/2)
                        point_normalx = Coefficients[0] * iter->x + (Coefficients[1] * iter->y + Coefficients[3]) / 2;
                        point_normaly = Coefficients[2] * iter->y + (Coefficients[1] * iter->x + Coefficients[4]) / 2;
                        point_normal = atan2(point_normaly,
                                             point_normalx); //��Ե��ķ��߷���,ָ����Բ��� (The normal direction of the edge point, pointing to the outside of the ellipse)
                        rec_inliers_cnt++;
                        if (angle_diff(point_normal, temp) <=
                            M_1_8_PI) //+- 22.5���� �� || d - r || < 3 dis_t (+- within 22.5�� and || d - r || < 3 dis_t)
                        {
                            rec_support_cnt++;
                            pixel_temp.x = iter->x;
                            pixel_temp.y = iter->y;
                            second_group_inliers.push_back(
                                    pixel_temp);//��Ӹ��߶ζ�Ӧ���ڵ� (Add the interior point corresponding to the line segment)
                        }
                    }
                }
            }
        }
        if (!(rec_support_cnt > 0 &&
              (rec_support_cnt >= 0.8 * lines[addr + 6] || rec_support_cnt * 1.0 / rec_inliers_cnt >= 0.6))) {
            flag2 = FALSE; //flag1 ��ʼ��ʱΪTRUE, һ��������һ���߶β�����Ҫ��ֱ��false, �ڵ�׼����֤��ͨ��
            // (flag1 is TRUE when initialized, once there is a line segment in the group that does not meet the requirements, it is directly false, and the interior point criterion verification fails)
            break;
        }
    }
    if (flag2 == TRUE && second_group_inliers.size() >= 0.8 *
                                                        group_inliers_num[second_group_ind])//�������ͳ�ƹ����ڵ�,ͨ����֤ (Close to the largest statistical interior point, pass the verification)
    {
        if (second_group_inliers.size() >=
            group_inliers_num[second_group_ind])//��������ֹ�������ڵ��� (The maximum number of inliers the update group has ever seen)
            group_inliers_num[second_group_ind] = second_group_inliers.size();
    } else
        flag2 = FALSE;
    //�ڶ����������֤ (The second group completes the verification)
    if (flag1 == TRUE && flag2 == TRUE) {
        auto *dataxy = (point2d *) malloc(
                sizeof(point2d) * (first_group_inliers.size() + second_group_inliers.size()));
        for (unsigned int i = 0; i < first_group_inliers.size(); i++) {
            dataxy[i].x = first_group_inliers[i].x;
            dataxy[i].y = first_group_inliers[i].y;
        }
        addr = first_group_inliers.size();
        for (unsigned int i = 0; i <
                                 second_group_inliers.size(); i++)//������������ʱһ��Ҫע��������Χ (Be sure to pay attention to the index range when concatenating two arrays)
        {
            dataxy[addr + i].x = second_group_inliers[i].x;
            dataxy[addr + i].y = second_group_inliers[i].y;
        }
//		drawEdge(img,dataxy,(first_group_inliers.size() + second_group_inliers.size()));
        info = fitEllipse(dataxy, (first_group_inliers.size() + second_group_inliers.size()), param2);
        free(dataxy); //�ͷ��ڴ� (free memory)
        //С���������Բ��Ҫ�ſ���� (Ellipses with small major and minor axes need to relax the parameters)
        if (param[2] <= 50)
            semimajor_errorratio = 0.25;
        else if (param[2] <= 100)
            semimajor_errorratio = 0.15;
        else
            semimajor_errorratio = 0.1;
        if (param[3] <= 50)
            semiminor_errorratio = 0.25;
        else if (param[3] <= 100)
            semiminor_errorratio = 0.15;
        else
            semiminor_errorratio = 0.1;
        if (param[2] <= 50 && param[3] <= 50)
            iscircle_ratio = 0.75;
        else if (param[2] >= 50 && param[2] <= 100 && param[3] >= 50 && param[3] <= 100)
            iscircle_ratio = 0.85;
        else
            iscircle_ratio = 0.9;
        if (info == 1 &&
            isEllipseEqual(param2, param, 3 * distance_tolerance, semimajor_errorratio, semiminor_errorratio, 0.1,
                           iscircle_ratio)) {
            ellipara->x = param2[0];//������Բ�����Ʒ�� (������Բ�����Ʒ��)
            ellipara->y = param2[1];
            ellipara->a = param2[2];
            ellipara->b = param2[3];
            ellipara->phi = param2[4];
            //drawEllipse(img,param2);
            return TRUE;
        }
    }
    return FALSE;
}


//���� (enter)
//lsd�㷨���õ����߶μ���lines������line_num��return�ķ���ֵ��line_nums���߶Σ�Ϊһάdouble������lines������Ϊ8*n��ÿ8��Ϊһ��
// (The number line_num of the line segment set lines detected by the lsd algorithm, the return value of return is line_nums line segments,
// which is a one-dimensional double-type array lines, with a length of 8*n, and each 8 is a group)
//����x1,y1,x2,y2,dx,dy,length,polarity (Store x1,y1,x2,y2,dx,dy,length,polarity)
//groups: �߶η��飬ÿ����水�ռ��ηֲ�˳��˳ʱ�������ʱ��洢���߶��������߶�������Χ��0~line_num-1
// (groups: line segment grouping, each group stores the line segment index clockwise or counterclockwise according
// 			to the geometric distribution order, the line segment index range is 0~line_num-1)
//coverages: ÿ������ĽǶȸ��Ƿ�Χ0~2pi���������ֻ��1���߶Σ����ǽǶ�Ϊ0�����鳤�ȵ��ڷ����������
// (coverages: The angle coverage of each group is 0~2pi. If there is only one line segment in the
// 			   group, the coverage angle is 0. The array length is equal to the number of groups.)
//angles ���Ե����ݶȷ���gradient direction, �ޱ�Ե��λNOTDEF
// (angles store the gradient direction of the edge point gradient direction, no edge point NOTDEF)
//����ֵ PairedGroupList* list ���ص��ǳ�ʼ��Բ���ϵ����飬����list->length.
// (Return value PairedGroupList* list Returns an array of initial ellipse sets, with length list->length.)
//�мǣ����ڴ��ں��������룬����ú����ǵ��ͷ��ڴ棬���ú���freePairedSegmentList()�����ͷ�
// (Remember, the memory is applied in the function, remember to release the memory after using the function, call the function freePairedSegmentList() to release)

PairGroupList *getValidInitialEllipseSet(const double *lines, int line_num, vector <vector<int>> *groups, const double *coverages,
                                         image_double angles, double distance_tolerance, int specified_polarity) {
    //���ټ��� (Accelerated Computing)
    //int* lineInliersIndex = (int*)malloc(sizeof(int)*line_num);//�����i���߶��ҵ����ڵ㣬���¼������Ϊj = length(supportInliers),��supportInliers.at(j)���Ÿ��߶ε�֧���ڵ�,û�ҵ��ڵ���߶ζ�Ӧ����Ϊ��ʼֵ-1.
    // (If the i-th line segment finds an interior point, record its index as j = length(supportInliers), that is, supportInliers.at(j) stores the support interior point of the line segment,
    // 		and the corresponding index of the line segment where the interior point is not found is the initial value - 1.)
    //vector<vector<point2d>> supportInliers;//������Ӧ�߶ε�֧���ڵ� (Save the supported interior points of the corresponding line segment)
    //memset(lineInliersIndex,-1,sizeof(int)*line_num);//�˴�Ҫʵ��ȷʵ���У������������Գ�ʼ��Ϊ0��-1.���ڸ�������ֻ����Ϊ0.
    // (It is indeed feasible to practice here. For integers, it can be initialized to 0, -1. For floating-point numbers, it can only be 0.)

    PairGroupList *pairGroupList = nullptr;
    PairGroupNode *head, *tail;
    int pairlength = 0;
    point2d pointG1s{}, pointG1e{}, pointG2s{}, pointG2e{}, g1s_ls_dir{}, g1e_ls_dir{}, g2s_ls_dir{}, g2e_ls_dir{};
    double polarity;
    point5d ellipara{};
    int groupsNum = (*groups).size();//������� (number of groups)
    auto *fitMatrixes = (double *) malloc(sizeof(double) * groupsNum *
                                            36);//������Ͼ���S_{6 x 6}. ÿ���鶼��һ����Ͼ��� (Define the fit matrix S_{6 x 6}. Each group has a fit matrix)
    auto *supportInliersNum = (unsigned int *) malloc(sizeof(int) *
                                                              groupsNum);//���ڴ洢ÿ�������������ֵ�֧���ڵ����� (Used to store the maximum number of supported inliers that have ever occurred for each group)
    memset(fitMatrixes, 0, sizeof(double) * groupsNum * 36);
    memset(supportInliersNum, 0, sizeof(unsigned int) * groupsNum);//��ʼ��Ϊ0. (Initialized to 0.)
    //double distance_tolerance = max( 2.0, 0.005*min(angles->xsize,angles->ysize) ); // 0.005%*min(xsize,ysize)
    int i, j;
    int cnt_temp, ind_start, ind_end;
    bool info;

    //ʵ������Ͼ���Si (Instantiate the fit matrix Si)
    auto *dataxy = (point2d *) malloc(sizeof(point2d) * line_num *
                                         2);//�����㹻���ڴ�, line_num���߶Σ�����2line_num���˵� (Apply for a large enough memory, line_num line segments, a total of 2line_num endpoints)
    for (i = 0; i < groupsNum; i++) {
        cnt_temp = 0;//ǧ��ע��Ҫ��0 (Be careful to clear 0)
        for (j = 0; j < (*groups)[i].size(); j++) {
            //ÿһ���߶���2���˵� (Each line segment has 2 endpoints)
            dataxy[cnt_temp].x = lines[(*groups)[i][j] * 8];
            dataxy[cnt_temp++].y = lines[(*groups)[i][j] * 8 + 1];
            dataxy[cnt_temp].x = lines[(*groups)[i][j] * 8 + 2];
            dataxy[cnt_temp++].y = lines[(*groups)[i][j] * 8 + 3];
        }
        calcuFitMatrix(dataxy, cnt_temp, fitMatrixes + i * 36);
    }
    free(dataxy);//�ͷ��ڴ� (free memory)

    head = tail = nullptr;//����ʼ��Բ���ϴ洢�������� (Store the initial set of ellipses into a linked list)
    //selection of salient elliptic hypothesis
    for (i = 0; i < groupsNum; i++) {
        if (coverages[i] >= M_4_9_PI)//����ĸ��ǽǶ�>= 4pi/9 = 80��, ������Ϊ���кܴ�������ԣ���ֱ�������ȡ
            // (When the coverage angle of the group is >= 4pi/9 = 80��, we think it has great significance and can be directly fitted and extracted)
        {
            //���뼫���ж�,ֻ��ȡָ�����Ե���Բ (Add polarity judgment, only extract ellipse with specified polarity)
            if (specified_polarity == 0 || (lines[(*groups)[i][0] * 8 + 7] == specified_polarity)) {
                //�����Դ�ĳ�ʼ��Բ��ȡ��һ���᷵��TRUE�����û��Ҫ���ж� (The initial ellipse extraction with large significance will definitely return TRUE, so there is no need to judge again)
                info = calcEllipseParametersAndValidate(lines, line_num, groups, i, -1, (fitMatrixes + i * 36), NULL,
                                                        angles, distance_tolerance, supportInliersNum, &ellipara);
                if (info == FALSE) {
                    continue;
                    error("getValidInitialEllipseSet, selection of salient ellipses failed!");//�����������֣���,��54.jpg���ָ����� (Will this happen? ? , the problem occurs when running 54.jpg)
                }
                auto *node = (PairGroupNode *) malloc(sizeof(PairGroupNode));
                node->center.x = ellipara.x;
                node->center.y = ellipara.y;
                node->axis.x = ellipara.a;
                node->axis.y = ellipara.b;
                node->phi = ellipara.phi;
                node->pairGroupInd.x = i;
                node->pairGroupInd.y = -1;//�� (without)
                if (head != nullptr) {
                    tail->next = node;
                    tail = node;
                } else {
                    head = tail = node;
                }
                pairlength++;
            }
        }
    }
    //selection of pair group hypothesis
    for (i = 0; i < groupsNum - 1; i++)
        for (j = i + 1; j < groupsNum; j++) {
            //���뼫���ж�,ֻ��ȡָ�����Ե���Բ (Add polarity judgment, only extract ellipse with specified polarity)
            if (specified_polarity == 0 || (lines[(*groups)[i][0] * 8 + 7] == specified_polarity)) {
                //group i 's polarity is the same as group j; and the number of two paired groups should be >= 3.
                if (lines[(*groups)[i][0] * 8 + 7] == lines[(*groups)[j][0] * 8 + 7] &&
                    ((*groups)[i].size() + (*groups)[j].size()) >= 3) {
                    ind_start = (*groups)[i][0];//��i����ʼһ���߶����� (The first line segment index of the i-th group)
                    ind_end = (*groups)[i][(*groups)[i].size() -
                                           1];//��i������һ���߶����� (last line segment index of group i)
                    pointG1s.x = lines[ind_start * 8];
                    pointG1s.y = lines[ind_start * 8 + 1];
                    g1s_ls_dir.x = lines[ind_start * 8 + 4];
                    g1s_ls_dir.y = lines[ind_start * 8 + 5];
                    pointG1e.x = lines[ind_end * 8 + 2];
                    pointG1e.y = lines[ind_end * 8 + 3];
                    g1e_ls_dir.x = lines[ind_end * 8 + 4];
                    g1e_ls_dir.y = lines[ind_end * 8 + 5];
                    ind_start = (*groups)[j][0];//��j����ʼһ���߶����� (The first line segment index of the jth group)
                    ind_end = (*groups)[j][(*groups)[j].size() -
                                           1];//��j������һ���߶����� (last line segment index of the jth group)
                    pointG2s.x = lines[ind_start * 8];
                    pointG2s.y = lines[ind_start * 8 + 1];
                    g2s_ls_dir.x = lines[ind_start * 8 + 4];
                    g2s_ls_dir.y = lines[ind_start * 8 + 5];
                    pointG2e.x = lines[ind_end * 8 + 2];
                    pointG2e.y = lines[ind_end * 8 + 3];
                    g2e_ls_dir.x = lines[ind_end * 8 + 4];
                    g2e_ls_dir.y = lines[ind_end * 8 + 5];
                    polarity = lines[ind_start * 8 + 7]; //i,j����ļ��� (The polarities of the two groups i, j)
                    if (regionLimitation(pointG1s, g1s_ls_dir, pointG1e, g1e_ls_dir, pointG2s, g2s_ls_dir, pointG2e,
                                         g2e_ls_dir, polarity, -3 * distance_tolerance))//���ڱ˴˵�����������
                    {
                        //if ( i == 2)
                        //	drawPairGroup(img,lines,(*groups),i,j);

                        if (calcEllipseParametersAndValidate(lines, line_num, groups, i, j, (fitMatrixes + i * 36),
                                                             (fitMatrixes + j * 36), angles, distance_tolerance,
                                                             supportInliersNum,
                                                             &ellipara))//����һ�㷽��������⣬�߶ε��ڵ�֧�ֱ���
                        {
                            auto *node = (PairGroupNode *) malloc(sizeof(PairGroupNode));
                            node->center.x = ellipara.x;
                            node->center.y = ellipara.y;
                            node->axis.x = ellipara.a;
                            node->axis.y = ellipara.b;
                            node->phi = ellipara.phi;
                            node->pairGroupInd.x = i;
                            node->pairGroupInd.y = -1;//�� (without)
                            if (head != nullptr) {
                                tail->next = node;
                                tail = node;
                            } else {
                                head = tail = node;
                            }
                            pairlength++;
                        }
                    }

                }
            }
        }
    if (pairlength > 0) {
        PairGroupNode *p;
        p = head;
        pairGroupList = pairGroupListInit(pairlength);
        for (int i = 0; i < pairGroupList->length; i++) {
            pairGroupList->pairGroup[i].center.x = p->center.x;
            pairGroupList->pairGroup[i].center.y = p->center.y;
            pairGroupList->pairGroup[i].axis.x = p->axis.x;
            pairGroupList->pairGroup[i].axis.y = p->axis.y;
            pairGroupList->pairGroup[i].phi = p->phi;
            pairGroupList->pairGroup[i].pairGroupInd.x = p->pairGroupInd.x;//��¼���(i,j),��groups�еĵ�i����͵�j���鹹�ɵ�ƥ�����������Ч��Բ���� (Record the group pair (i, j), the matching group consisting of the i-th group and the j-th group in groups produces the valid ellipse parameter)
            pairGroupList->pairGroup[i].pairGroupInd.y = p->pairGroupInd.y;
            p = p->next;
        }
        tail->next = nullptr;
        while (head != nullptr) {
            p = head;
            head = head->next;
            free(p);
        }
    }
    //supportInliers.resize(0);
    //free(lineInliersIndex);//�ͷ��߶��ڵ������ (release the index of the point within the segment)
    free(supportInliersNum);//�ͷŴ洢�������֧���ڵ����������� (Free the array that stores the number of supported interior points for each group)
    free(fitMatrixes);//�ͷŴ洢���������Ͼ��� (release the fit matrix that stores the individual groups)
    return pairGroupList;
}


void generateEllipseCandidates(PairGroupList *pairGroupList, double distance_tolerance, double *&ellipse_candidates,
                               int *candidates_num) {
    if (pairGroupList->length <=
        0)//��⣬����Ҫ��1����������������ѡ (detection, at least 1 sample must be used to generate candidates)
    {
        ellipse_candidates = nullptr;
        (*candidates_num) = 0;
        return;
    }
    double *centers;
    int center_num; //��Բ����(xi,yi)�ľ������� (Number of clusters at the center of the ellipse (xi,yi))
    double *phis;
    int phi_num;    //���ÿһ����Բ����(xi,yi)����б�Ƕ�phi�ľ������� (For each ellipse center (xi,yi), the number of clusters at the tilt angle phi)
    double *axises;
    int axis_num;   //���ÿһ����Բ���ĺ����(xi,yi,phi),���̰���(a,b)�ľ������� (For each ellipse center and inclination (xi, yi, phi), the number of clusters of the major and minor semi-axes (a, b))
    auto *bufferXY = (double *) calloc(pairGroupList->length * 2, sizeof(double));
    auto *bufferPhi = (double *) calloc(pairGroupList->length, sizeof(double));
    auto *bufferAB = (double *) calloc(pairGroupList->length * 2, sizeof(double));
    auto *bufferIndexes = (point2i *) calloc(pairGroupList->length,
                                                sizeof(point2i));//point[i].x��¼��i��������bufferXX�е���ʼ����λ�ã�point[i].y��¼��i��������bufferXX�еĳ���
    // (point[i].x records the starting index position of the i-th category in bufferXX, point[i].y records the length of the i-th category in bufferXX)
    auto *buffer2AB = (double *) calloc(pairGroupList->length * 2, sizeof(double));
    auto *buffer2Indexes = (point2i *) calloc(pairGroupList->length,
                                                 sizeof(point2i));//point[i].x��¼��i��������bufferXX�е���ʼ����λ�ã�point[i].y��¼��i��������bufferXX�еĳ���
    // (point[i].x records the starting index position of the i-th category in bufferXX, point[i].y records the length of the i-th category in bufferXX)
    int *buffer_temp = (int *) calloc(pairGroupList->length, sizeof(int));
    int addr, addr2, info, ind;
    double dis_min, dis_temp;
    if (bufferXY == nullptr || bufferPhi == nullptr || bufferAB == nullptr || bufferIndexes == nullptr ||
        buffer2AB == nullptr || buffer2Indexes == nullptr || buffer_temp == nullptr
            ) {
        ellipse_candidates = nullptr;
        (*candidates_num) = 0;
        error("generateEllipseCandidates, not enough memory");
    }
    (*candidates_num) = 0; //��ѡ��Բ��������ʼ��Ϊ0,�ǳ���Ҫ (The number of candidate ellipses, initialized to 0, very important)
    //copy
    for (int i = 0; i < pairGroupList->length; i++) {
        addr = 2 * i;
        bufferXY[addr] = pairGroupList->pairGroup[i].center.x;
        bufferXY[addr + 1] = pairGroupList->pairGroup[i].center.y;
    }
    //cluster the ellipses' centers
    info = cluster2DPoints(bufferXY, pairGroupList->length, distance_tolerance, centers, &center_num);
    if (info == 0) {
        ellipse_candidates = nullptr;
        (*candidates_num) = 0;
        error("generateEllipseCandidates, cluster2DPoints, error in clustering elliptic centers");
    }
    //classification,Ѱ��ÿ��������ľ������� (classification, find the cluster center to which each point belongs)
    for (int i = 0; i < pairGroupList->length; i++) {
        dis_min = DBL_MAX;
        ind = -1;
        for (int j = 0; j < center_num; j++) {
            addr = 2 * j;
            dis_temp = (pairGroupList->pairGroup[i].center.x - centers[addr]) *
                       (pairGroupList->pairGroup[i].center.x - centers[addr]) +
                       (pairGroupList->pairGroup[i].center.y - centers[addr + 1]) *
                       (pairGroupList->pairGroup[i].center.y - centers[addr + 1]);
            if (dis_temp < dis_min) {
                dis_min = dis_temp;
                ind = j; //record the nearest center's index
            }
        }
        buffer_temp[i] = ind; //�˴�����buffer2�����µ�i����ʼ��Բ��Ӧ��ind����Բ�������� (Here, buffer2 is used to record the i-th initial ellipse corresponding to the ind-th ellipse cluster center)
    }
    //����������˳��浽bufferXY,bufferPhi,bufferAB�У���bufferIndexes[i]���ŵ�i���������ĵ���ʼ����λ�úͳ���
    // (Store the classification results in bufferXY, bufferPhi, bufferAB in order, and bufferIndexes[i]
    //  stores the starting index position and length of the i-th cluster center)
    memset(bufferIndexes, 0, sizeof(point2i) * pairGroupList->length);
    ind = 0;//���㣬��������ʼλ�ã�����λ����ind*2,�����Ļ�ַ (Cleared, the starting position of the sample point, the index position is ind*2, the base address of the partition)
    for (int i = 0; i < center_num; i++) {
        bufferIndexes[i].x = ind;
        for (int j = 0; j < pairGroupList->length; j++) {
            if (buffer_temp[j] == i) {
                addr = ind *
                       2;//�мǳ��̰�����һ��һ��索�ģ���Ҫ x 2 (Remember that the long and short semi-axes are stored in groups of inches, and you need x 2)
                addr2 = bufferIndexes[i].y * 2;
                bufferPhi[ind + bufferIndexes[i].y] = pairGroupList->pairGroup[j].phi;
                bufferAB[addr + addr2] = pairGroupList->pairGroup[j].axis.x;
                bufferAB[addr + addr2 + 1] = pairGroupList->pairGroup[j].axis.y;
                bufferIndexes[i].y++;//��i������������Χ�ĵ�������1 (Add 1 to the number of points around the i-th cluster center)
            }
        }
        if (bufferIndexes[i].y == 0)//����������Χû�п����ĵ� (There are no close points around the cluster center)
        {
            error("generateEllipseCandidates, no XY points near to the clustering center");
        }
        ind += bufferIndexes[i].y;
    }
    //cout<<"2D cluster centers over"<<endl;
    //��ÿһ����Բ���ĵ���Χ�ĵ������Ǿ��� (Dip clustering of points around the center of each ellipse)
    //��i����Բ�������ģ����ڽ����������Χ�ǣ�bufferIndexs[i].x ~ (bufferIndex[i].x + bufferIndex[i].y-1)
    // (The i-th ellipse cluster center, the index range of its neighboring points is: bufferIndexs[i].x ~ (bufferIndex[i].x + bufferIndex[i].y-1))
    for (int i = 0; i < center_num; i++) {


        double *phi_pointer_temp = bufferPhi + bufferIndexes[i].x;//���ָ�� (Tilt pointer)
        double *ab_pointer_temp = bufferAB + bufferIndexes[i].x *
                                             2;//���̰����ָ��,��ס x 2 (The pointer of the long and short semi-axis, remember x 2)
        info = cluster1DDatas(phi_pointer_temp, bufferIndexes[i].y, 0.0873, phis,
                              &phi_num);//��phi����, pi/180*5 = 0.0873, 5����� (For phi clustering, pi/180*5 = 0.0873, 5�� error)
        if (info ==
            0) //����Ϊʲô����������centers[i]����Χ����û��������ĵ�,����bufferIndexes[i].y = 0 (I don't know why, the cluster center centers[i] may not have the closest point around it, the number bufferIndexes[i].y = 0)
        {
            //cout<<"generateEllipseCandidates, cluster2DPoints, error in clustering elliptic phis"<<endl;
            continue;
            //error("generateEllipseCandidates, cluster2DPoints, error in clustering elliptic phis");
        }
        //classification,Ѱ��ÿ��������ľ������� (classification, find the cluster center to which each point belongs)
        for (int j = 0; j < bufferIndexes[i].y; j++) {
            dis_min = DBL_MAX;
            ind = -1;
            for (int k = 0; k < phi_num; k++) {
                dis_temp = (*(phi_pointer_temp + j) - phis[k]) * (*(phi_pointer_temp + j) - phis[k]);
                if (dis_temp < dis_min) {
                    dis_min = dis_temp;
                    ind = k;//record the nearest phi's index
                }
            }
            buffer_temp[j] = ind;
        }
        //����������˳��洢��buffer2AB�У���buffer2Indexes[j].x��Ӧ��i��phi�ľ���������ʼ�㣬buffer2Indexes[j].y��Ӧ����(����)
        // (Store the classification results in buffer2AB in order, and buffer2Indexes[j].x corresponds to the starting
        //  point of the cluster center of the i-th phi, and buffer2Indexes[j].y corresponds to the number (length))
        memset(buffer2Indexes, 0, sizeof(point2i) * bufferIndexes[i].y);
        ind = 0;
        for (int j = 0; j < phi_num; j++) {
            buffer2Indexes[j].x = ind;//��ʼ�� (starting point)
            for (int k = 0; k < bufferIndexes[i].y; k++) {
                if (buffer_temp[k] == j) {
                    addr = ind * 2;
                    addr2 = buffer2Indexes[j].y * 2;
                    buffer2AB[addr + addr2] = *(ab_pointer_temp + k * 2);
                    buffer2AB[addr + addr2 + 1] = *(ab_pointer_temp + k * 2 + 1);
                    buffer2Indexes[j].y++;//���ȼ�1 (length plus 1)
                }
            }
            ind += buffer2Indexes[j].y;
        }
        for (int j = 0; j < phi_num; j++) {
            double *ab_pointer_temp2 = buffer2AB + buffer2Indexes[j].x *
                                                   2; //���̰����ָ��,��ס x 2 (The pointer of the long and short semi-axis, remember x 2)
            info = cluster2DPoints(ab_pointer_temp2, buffer2Indexes[j].y, distance_tolerance, axises, &axis_num);
            if (info ==
                0) //����Ϊʲô����������phi_j����Χ����û��������ĵ�,����buffer2Indexes[j].y = 0 (I don't know why, the cluster center phi_j may not have the closest point around it, the number buffer2Indexes[j].y = 0)
            {
                //cout<<"generateEllipseCandidates, cluster2DPoints, error in clustering elliptic axises"<<endl;
                continue;
                //error("generateEllipseCandidates, cluster2DPoints, error in clustering elliptic axises");
            }
            //����ѡ��Բ��д��bufferXY,bufferPhi,bufferAB����, ��ѡ��Բ����(*candidates_num)++
            // (Rewrite the candidate ellipse into bufferXY, bufferPhi, bufferAB, the number of candidate ellipses (*candidates_num)++)
            for (int k = 0; k < axis_num; k++) {
                addr = (*candidates_num) * 2;
                bufferXY[addr] = centers[i * 2];
                bufferXY[addr + 1] = centers[i * 2 + 1];
                bufferPhi[(*candidates_num)] = phis[j];
                bufferAB[addr] = axises[k * 2];
                bufferAB[addr + 1] = axises[k * 2 + 1];
                (*candidates_num)++;
            }
            free(axises);//cluster2DPoints�ϸ�Ҫ������axises����Ҫ�ͷź����ڲ�������ڴ�
            // (Cluster2DPoints is strictly required. After axises are used up, the memory allocated inside the function needs to be released)
        }
        free(phis);//cluster1DDatas�ϸ�Ҫ������phis����Ҫ�ͷź����ڲ�������ڴ�
        // (Cluster1DDatas is strictly required. After using phis, you need to release the memory requested by the function.)
    }
    free(centers);//cluster2DPoints�ϸ�Ҫ������centers����Ҫ�ͷź����ڲ�������ڴ�
    // (Cluster2DPoints is strictly required. After running out of centers, you need to release the memory requested by the function.)
    //�ͷ��ں�����ͷ����Ĳ����ڴ� (Free part of the memory allocated at the beginning of the function)
    free(buffer_temp); //�˴��ͷų����� (issue here)
    free(buffer2Indexes);
    free(buffer2AB);
    free(bufferIndexes);
    ellipse_candidates = (double *) malloc(sizeof(double) * (*candidates_num) * 5);
    for (int i = 0; i < (*candidates_num); i++) {
        addr = 2 * i;
        ellipse_candidates[i * 5] = bufferXY[addr];
        ellipse_candidates[i * 5 + 1] = bufferXY[addr + 1];
        ellipse_candidates[i * 5 + 2] = bufferAB[addr];
        ellipse_candidates[i * 5 + 3] = bufferAB[addr + 1];
        ellipse_candidates[i * 5 + 4] = bufferPhi[i];
    }
    //�ͷ��ں�����ͷ������ڴ� (Free the memory allocated at the beginning of the function)
    free(bufferAB);
    free(bufferPhi);
    free(bufferXY);
    if ((*candidates_num) <= 0) {
        *candidates_num = 0;
        ellipse_candidates = nullptr;
        //cout<<"no any candidates generated!"<<endl;
    }
}







//==========================================END=======================================================================
/**
���룺(enter)
prhs[0]: ����ĻҶ�ͼ�񣬵�ͨ������С��imgy x imgx (input grayscale image, single channel, size imgy x imgx)
prhs[1]: ��Ե��ȡѡ��1 canny; 2 sobel (edge extraction selection, 1 canny; 2 sobel)
prhs[2]: ���ָ������Բ���� (Detects the specified ellipse polarity)
�����(output)
plhs[0]: ��ѡ��Բ���(xi,yi,ai,bi,phi_i)', 5 x m (Candidate ellipse combination (xi,yi,ai,bi,phi_i)', 5 x m)
plhs[1]: ��Եͼ����С��imgy x imgx�����Ե������Ϊ edgepix_n. ��ֵ����0 ���� 255
		 (Edge map, the size is imgy x imgx, let the total number of edge points be edgepix_n. Binarization, 0 or 255)
plhs[2]: ��Ե����ݶ��������󣬴�С�� 2 x edgepix_n, (cos(theta_rad),sin(theta_rad))'...
		 (Matrix of gradient vectors of edge points of size 2 x edgepix_n, (cos(theta_rad),sin(theta_rad))'...)
plhs[3]: �߶�ͼ����С��imgy x imgx (Line segment plot, size is imgy x imgx)
*/
/*
compile��
mex generateEllipseCandidates.cpp -IF:\OpenCV\opencv2.4.9\build\include -IF:\OpenCV\opencv2.4.9\build\include\opencv -IF:\OpenCV\opencv2.4.9\build\include\opencv2 -LF:\OpenCV\opencv2.4.9\build\x64\vc11\lib -IF:\Matlab\settlein\extern\include -LF:\Matlab\settlein\extern\lib\win64\microsoft -lopencv_core249 -lopencv_highgui249 -lopencv_imgproc249 -llibmwlapack.lib
*/
//======================================MEX function==================================================================

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if (nrhs != 3)
        mexErrMsgIdAndTxt("MATLAB:revord:invalidNumInputs", "One input required.");
    else if (nlhs > 4)
        mexErrMsgIdAndTxt("MATLAB:revord:maxlhs", "Too many output arguments.");
    uchar *inputimg = (uchar *) mxGetData(prhs[0]);
    int imgy, imgx;
    int edge_process_select = (int) mxGetScalar(
            prhs[1]);//��Ե��ȡѡ��1 canny; 2 sobel (edge extraction selection, 1 canny; 2 sobel)
    int specified_polarity = (int) mxGetScalar(prhs[2]);//1,ָ��������Բ����ҪΪ��; -1ָ������Ϊ��; 0��ʾ���ּ�����Բ�����
    // (1, the polarity of the specified ellipse to be detected must be positive; -1 specifies that the polarity is negative; 0 means that both polar ellipses are detected)
    imgy = (int) mxGetM(prhs[0]);
    imgx = (int) mxGetN(prhs[0]);
    auto *data = (double *) malloc(imgy * imgx *
                                     sizeof(double));//����������е�ͼ������ת�浽һά������ (Dump the image data in the input matrix into a 1D array)
    for (int c = 0; c < imgx; c++) {
        for (int r = 0; r < imgy; r++) {
            data[c + r * imgx] = inputimg[r + c * imgy];
        }
    }
    int n;//�߶����� (Number of line segments)
    //int new_n;
    vector <vector<int>> groups;
    double *coverages;
    int *reg;
    int reg_x;
    int reg_y;
    double *out = mylsd(&n, data, imgx, imgy, &reg, &reg_x, &reg_y);
    groupLSs(out, n, reg, reg_x, reg_y, &groups);//���� (grouping)
    free(reg); //�ͷ��ڴ� (free memory)
    calcuGroupCoverage(out, n, groups, coverages);//����ÿ����ĸ��ǽǶ� (Calculate the coverage angle for each group)

    printf("The number of output arc-support line segments: %i\n", n);
    printf("The number of arc-support groups:%i\n", groups.size());
    /*int groups_t = 0;
	for (int i = 0; i<groups.size(); i++)
	{
		groups_t+= groups[i].size();
	}
	printf("Groups' total ls num:%i\n",groups_t);*/

    image_double angles;
    if (edge_process_select == 1)
        calculateGradient2(data, imgx, imgy, &angles); //version2, sobel; version 3 canny
    else
        calculateGradient3(data, imgx, imgy, &angles); //version2, sobel; version 3 canny
    PairGroupList *pairGroupList;
    double distance_tolerance = 2;//max( 2.0, 0.005*min(angles->xsize,angles->ysize) ); // 0.005%*min(xsize,ysize)
    double *candidates; //��ѡ��Բ (candidate ellipse)
    double *candidates_out;//�����ѡ��Բָ�� (output candidate ellipse pointer)
    int candidates_num = 0;//��ѡ��Բ���� (Number of candidate ellipses)
    //rejectShortLines(out,n,&new_n);
    pairGroupList = getValidInitialEllipseSet(out, n, &groups, coverages, angles, distance_tolerance,
                                              specified_polarity);
    if (pairGroupList != nullptr) {
        printf("The number of initial ellipses��%i \n", pairGroupList->length);
        generateEllipseCandidates(pairGroupList, distance_tolerance, candidates, &candidates_num);
        printf("The number of ellipse candidates: %i \n", candidates_num);

        plhs[0] = mxCreateDoubleMatrix(5, candidates_num, mxREAL);
        candidates_out = (double *) mxGetPr(plhs[0]);
        //��ѡԲ���(xi,yi,ai,bi,phi_i)', 5 x candidates_num, ���Ƶ�����candidates_out��
        // (Candidate circle combination (xi,yi,ai,bi,phi_i)', 5 x candidates_num, copied to matrix candidates_out)
        memcpy(candidates_out, candidates, sizeof(double) * 5 * candidates_num);

        freePairGroupList(pairGroupList);
        free(candidates);
    } else {
        printf("The number of initial ellipses��%i \n", 0);
        double *candidates_out;
        plhs[0] = mxCreateDoubleMatrix(5, 1, mxREAL);
        candidates_out = (double *) mxGetPr(plhs[0]);
        candidates_out[0] = candidates_out[1] = candidates_out[2] = candidates_out[3] = candidates_out[4] = 0;
    }
    uchar *edgeimg_out;
    unsigned long edge_pixels_total_num = 0;//��Ե������ (total edge pixels)
    double *gradient_vec_out;
    plhs[1] = mxCreateNumericMatrix(imgy, imgx, mxUINT8_CLASS, mxREAL);
    edgeimg_out = (uchar *) mxGetData(plhs[1]);
    //����Եͼ���Ƶ�����edgeimg_out�� (Copy the edge map into the matrix edgeimg_out)
    //���ݶ������浽����gradient_vec_out�� (Store the gradient vector in the matrix gradient_vec_out)
    unsigned long addr, g_cnt = 0;
    for (int c = 0; c < imgx; c++)
        for (int r = 0; r < imgy; r++) {
            addr = r * imgx + c;
            if (angles->data[addr] == NOTDEF)
                edgeimg_out[c * imgy + r] = 0;
            else {
                edgeimg_out[c * imgy + r] = 255;//Ϊ��Ե�㣬��ֵΪ��ɫ (is the edge point, assigned as white)
                //------------------------------------------------
                edge_pixels_total_num++;
            }
        }
    printf("edge pixel number: %i\n", edge_pixels_total_num);
    //����edge_pixels_total_num x 2 ������ÿһ����Ե����ݶ�����������Ϊ���ȣ�����matlab��ϰ��
    // (Apply edge_pixels_total_num x 2 to save the gradient vector of each edge point, listed first, in line with the habit of matlab)
    plhs[2] = mxCreateDoubleMatrix(2, edge_pixels_total_num, mxREAL);
    gradient_vec_out = (double *) mxGetPr(plhs[2]);
    for (int c = 0; c < imgx; c++)
        for (int r = 0; r < imgy; r++) {
            addr = r * imgx + c;
            if (angles->data[addr] != NOTDEF) {
                gradient_vec_out[g_cnt++] = cos(angles->data[addr]);
                gradient_vec_out[g_cnt++] = sin(angles->data[addr]);
            }
        }
    //---------------------------------------------------------------------
    //����߶μ���ͼ�� (Output image for line segment detection)
    if (nlhs == 4) {
        Mat ls_mat = Mat::zeros(imgy, imgx, CV_8UC1);
        for (int i = 0; i < n; i++)//draw lines
        {
            Point2d p1(out[8 * i], out[8 * i + 1]), p2(out[8 * i + 2], out[8 * i + 3]);
            line(ls_mat, p1, p2, Scalar(255, 0, 0));
        }
        if (candidates_num > 0)//draw ellipses
        {
            for (int i = 0; i < candidates_num; i++)
                ellipse(ls_mat, cv::Point((int) candidates_out[i * 5], (int) candidates_out[i * 5 + 1]),
                        cv::Size(candidates_out[i * 5 + 2], candidates_out[i * 5 + 3]),
                        candidates_out[i * 5 + 4] * 180 / M_PI, 0, 360, (Scalar(255, 0, 0)), 1);
        }
        plhs[3] = mxCreateDoubleMatrix(imgy, imgx, mxREAL);
        double *ls_img_out = (double *) mxGetPr(plhs[3]);
        //memcpy(ls_out_mat,ls_mat.data ,sizeof(unsigned char)*M*N);
        for (int i = 0; i < imgx; i++)
            for (int j = 0; j < imgy; j++)
                ls_img_out[i * imgy + j] = ls_mat.data[j * imgx + i];
    }
    //---------------------------------------------------------------------
    //�����free���ͷų��������ڲ�����ѡԲ���õ���һϵ���ڴ�
    // (The free here is to release a series of memory used in the program to generate candidate circles)
    free(data);
    free(coverages);
    free(out);
    free_image_double(angles);

}










/*
void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
	int m = mxGetM(prhs[0]);
	int n = mxGetN(prhs[0]);
	double * p = (double*)mxGetData(prhs[0]);
	int sum = 0;
	for (int c = 0; c<n; c++)
		for ( int r = 0; r<m; r++)
			sum += p[c*m+r];
	plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
	double *pout = mxGetPr(plhs[0]);
	*pout = sum;

}
*/
