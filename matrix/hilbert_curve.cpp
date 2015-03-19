/*
 * The code in this file is copied from wikipedia.
 * http://en.wikipedia.org/wiki/Hilbert_curve
 */
#include "hilbert_curve.h"

//rotate/flip a quadrant appropriately
void hilbert_rot(int n, int *x, int *y, int rx, int ry)
{
	if (ry == 0) {
		if (rx == 1) {
			*x = n-1 - *x;
			*y = n-1 - *y;
		}

		//Swap x and y
		int t  = *x;
		*x = *y;
		*y = t;
	}
}

//convert (x,y) to d
int hilbert_xy2d (int n, int x, int y)
{
	int rx, ry, s, d=0;
	for (s=n/2; s>0; s/=2) {
		rx = (x & s) > 0;
		ry = (y & s) > 0;
		d += s * s * ((3 * rx) ^ ry);
		hilbert_rot(s, &x, &y, rx, ry);
	}
	return d;
}
