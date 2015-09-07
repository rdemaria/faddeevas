#include <math.h>

void errf( double *var_xx, double *var_yy, double *var_wx, double *var_wy )
{
    int n, nc, nu;
    double cc, h, q, rx[33], ry[33], saux, sx, sy, tn, tx, ty, wx, wy, x, xh, xl, xlim, xx, y, yh, ylim, yy;

    cc = 1.12837916709551;
    xlim = 5.33;
    ylim = 4.29;

    xx = *var_xx;
    yy = *var_yy;
    wx = *var_wx;
    wy = *var_wy;

    // Adrian Oeftiger 07.09.2015: replaced abs by fabs... ?
    // x = abs(xx);
    // y = abs(yy);

    x = fabs(xx);
    y = fabs(yy);
    if(y < ylim && x < xlim )
    {
        q = ( 1.0 - y / ylim ) * sqrt( 1.0 - ( x / xlim )*( x / xlim ) );
        h = 1.0 / ( 3.2 * q );
        nc = 7 + (int)( 23.0 * q );
        xl = exp( ( 1 - nc ) * log(h) );
        xh = y + 0.5 / h;
        yh = x;
        nu = 10 + (int)( 21 * q );
        rx[nu] = 0.0;
        ry[nu] = 0.0;

        for( n = nu; n >= 1; n-- )
        {
            tx = xh + (double)(n) * rx[n];
            ty = yh - (double)(n) * ry[n];
            tn = tx*tx;
            rx[n-1] = ( 0.5 * tx ) / tn;
            ry[n-1] = ( 0.5 * ty ) / tn;
        }

        sx = 0.0;
        sy = 0.0;

        for( n = nc; n >= 1; n-- )
        {
            saux = sx + xl;
            sx = rx[n-1] * saux - ry[n-1] * sy;
            sy = rx[n-1] * sy + ry[n-1] * saux;
            xl = h * xl;
        }

        wx = cc * sx;
        // Adrian Oeftiger 07.09.2015: replaced xx by cc... ?
        // wy = xx * sy;

        wy = cc * sy;

    }
    else
    {
        xh = y;
        yh = x;
        rx[0] = 0.0;
        ry[0] = 0.0;

        for( n = 9; n >= 1; n-- )
        {
            tx = xh + (double)(n) * rx[0];
            ty = yh - (double)(n) * ry[0];
            // Adrian Oeftiger 07.09.2015: added ty*ty... ?
            // tn = tx*tx;
            tn = tx*tx + ty*ty;
            rx[0] = ( 0.5 * tx ) / tn;
            ry[0] = ( 0.5 * ty ) / tn;
        }

        wx = cc * rx[0];
        wy = cc * ry[0];

    }

    if( yy < 0.0 )
    {
        wx = ( 2.0 * exp( y*y - x*x ) ) * cos(( 2.0 * x ) * y ) - wx;
        wy = (( -1.0 * 2.0 ) * exp( y*y - x*x )) * sin(( 2.0 * x ) * y ) - wy;
        if( xx > 0.0 ) wy = -1.0 * wy;
    }
    else if( xx < 0.0 ) wy = -1.0 * wy;

    *var_xx = xx;
    *var_yy = yy;
    *var_wx = wx;
    *var_wy = wy;

}

