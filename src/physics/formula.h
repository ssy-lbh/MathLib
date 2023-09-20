#ifndef FORMULA_H
#define FORMULA_H

#include "math/math_base.h"
#include "math/matrix.h"
#include "math/complex.h"

#include "units.h"

// g_uv = diag(1, 1, 1, -1) 平直时空，时空间隔为 s^2 = (x, y, z, ct)^\tilde \cdot g \cdot (x, y, z, ct)
// F_uv = ((E_x, B_xy, B_xz, iB_xt), (B_xy, E_y, B_yz, iB_yt), (B_xz, B_yz, E_z, iB_zt), (-iB_xt, -iB_yt, -iB_zt, -E_t))
// 黎曼曲率张量 R_uv = \partial_u \Gamma_v - \partial_v \Gamma_u + \Gamma_u \Gamma_v - \Gamma_v \Gamma_u

//TODO 复矩阵计算库

#endif