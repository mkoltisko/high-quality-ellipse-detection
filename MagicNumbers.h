

#ifndef MAGICNUMBERS_H
#define MAGICNUMBERS_H

#ifndef FALSE
#define FALSE 0
#endif /* !FALSE */

#ifndef TRUE
#define TRUE 1
#endif /* !TRUE */

/** Label for pixels with undefined gradient. */
#define NOTDEF (-1024.0)
/** PI */
#ifndef M_PI
#define M_PI   3.14159265358979323846
#endif /* !M_PI */
#define M_1_2_PI 1.57079632679489661923
#define M_1_4_PI 0.785398163

#define M_3_4_PI 2.35619449

#define M_1_8_PI 0.392699081
#define M_3_8_PI 1.178097245
#define M_5_8_PI 1.963495408
#define M_7_8_PI 2.748893572
#define M_4_9_PI 1.396263401595464  //80°
#define M_1_9_PI  0.34906585  //20°
#define M_1_10_PI 0.314159265358979323846   //18°
#define M_1_12_PI 0.261799387   //15°
#define M_1_15_PI 0.20943951    //12°
#define M_1_18_PI 0.174532925   //10°
/** 3/2 pi */
#define M_3_2_PI 4.71238898038
/** 2 pi */
#define M_2__PI  6.28318530718
/** Doubles relative error factor
 */
#define RELATIVE_ERROR_FACTOR 100.0


//=================================================================================================================
//===========================================LSD functions=========================================================
/** ln(10) */
#ifndef M_LN10
#define M_LN10 2.30258509299404568402    //ln10
#endif /* !M_LN10 */

/** Label for pixels not used in yet. */
#define NOTUSED 0

/** Label for pixels already used in detection. */
#define USED    1

//对于构成圆弧的像素标记极性，如果梯度的方向和弧的方向指向一致，则为SAME_POLE,否则为OPP_POLE,该标记初始是为0
// (For the polarity of the pixel marker that forms the arc, if the direction of the gradient is the same as
//  the direction of the arc, it is SAME_POLE, otherwise it is OPP_POLE, and the marker is initially 0)
#define NOTDEF_POL 0
#define SAME_POL 1
#define OPP_POL  (-1)

#endif //_MAGICNUMBERS_H
