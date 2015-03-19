#ifndef __HILBERT_CURVE_H__
#define __HILBERT_CURVE_H__

/*
 * This converts 2d coordinates x and y into the location in the hilbert curve.
 * n indicates the size of the square (2^n x 2^n).
 */
int hilbert_xy2d (int n, int x, int y);

#endif
