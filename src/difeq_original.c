void difeq(int k, int k1, int k2, int jsf, int is1, int isf, int ne, matrix_t D, matrix_t Y)
{
    real_t xt, y0t, y1t, y2t, kp, km, alpha, beta, gamma;
    kp = fabsl(m + s) / 2.0;
    km = fabsl(m - s) / 2.0;

    if (k == k1)
    {
        // x = -1 boundary
        MATRIX_SET(D, 2, 3, (powl(c, 2.0) + 2.0 * c * s - (km + kp - s) * (1.0 + km + kp + s) + MATRIX_GET(Y, 2, 0)) / (2.0 + 4.0 * km));
        MATRIX_SET(D, 2, 4, 1.0);
        MATRIX_SET(D, 2, 5, MATRIX_GET(Y, 0, 0) / (2.0 + 4.0 * km));

        MATRIX_SET(D, 2, jsf - 1, MATRIX_GET(Y, 1, 0) + MATRIX_GET(Y, 0, 0) * MATRIX_SET(D, 2, 3));
    }
    else if (k > k2)
    {
        // x = +1 boundary
        MATRIX_SET(D, 0, 3, (2.0 * km * (1 + km + kp) - powl(c - s, 2.0) - s - MATRIX_GET(Y, 2, N_PTS - 1)) / (2.0 + 4.0 * kp));
        MATRIX_SET(D, 0, 4, 1.0);
        MATRIX_SET(D, 0, 5, MATRIX_GET(Y, 0, N_PTS - 1) / (2.0 + 4.0 * km));

        MATRIX_SET(D, 0, jsf - 1,
                   MATRIX_GET(Y, 1, N_PTS - 1) - MATRIX_GET(Y, 0, N_PTS - 1) * ((4.0 * MATRIX_GET(Y, 2, N_PTS - 1) + 4.0 * c * c - 4.0 * kp - 4.0 * kp * kp - m * m + 4.0 * s - 8.0 * c * s + 2.0 * m * s + 3.0 * s * s) / (16.0 * kp + 8.0 * kp * kp - 2.0 * (-4.0 + m * m + 2.0 * m * s + s * s)) + km / 2.0));

        MATRIX_SET(D, 1, 3, 1.0);
        MATRIX_SET(D, 1, 4, 0.0);
        MATRIX_SET(D, 1, 5, 0.0);

        MATRIX_SET(D, 1, jsf - 1, MATRIX_GET(Y, 0, N_PTS - 1) - anormr);
    }
    else
    {
        // Internal points
        MATRIX_SET(D, 0, 0, -1.0);
        MATRIX_SET(D, 0, 1, -0.5 * STEP);
        MATRIX_SET(D, 0, 2, 0.0);
        MATRIX_SET(D, 0, 3, 1.0);
        MATRIX_SET(D, 0, 4, -0.5 * STEP);
        MATRIX_SET(D, 0, 5, 0.0);

        xt = (VECTOR_GET(x, k) + VECTOR_GET(x, k - 1)) / 2.0;
        y0t = (MATRIX_GET(Y, 0, k) + MATRIX_GET(Y, 0, k - 1)) / 2.0;
        y1t = (MATRIX_GET(Y, 1, k) + MATRIX_GET(Y, 1, k - 1)) / 2.0;
        y2t = (MATRIX_GET(Y, 2, k) + MATRIX_GET(Y, 2, k - 1)) / 2.0;

        MATRIX_SET(D, 1, 0, 0.5 * STEP * (powl(c * xt, 2.0) - 2.0 * c * s * xt + y2t - (km + kp - s) * (1.0 + km + kp + s)) / (1 - powl(xt, 2.0)));
        MATRIX_SET(D, 1, 1, -1.0 - STEP * (-km + kp + (1.0 + km + kp) * xt) / (1.0 - powl(xt, 2.0)));
        MATRIX_SET(D, 1, 2, 0.5 * STEP * y0t / (1.0 - powl(xt, 2.0)));
        MATRIX_SET(D, 1, 3, MATRIX_GET(D, 1, 0));
        MATRIX_SET(D, 1, 4, MATRIX_GET(D, 1, 1) + 2.0);
        MATRIX_SET(D, 1, 5, MATRIX_GET(D, 1, 2));

        MATRIX_SET(D, 2, 0, 0.0);
        MATRIX_SET(D, 2, 1, 0.0);
        MATRIX_SET(D, 2, 2, -1.0);
        MATRIX_SET(D, 2, 3, 0.0);
        MATRIX_SET(D, 2, 4, 0.0);
        MATRIX_SET(D, 2, 5, 1.0);

        MATRIX_SET(D, 0, jsf - 1,
                   MATRIX_GET(Y, 0, k - 1) - MATRIX_GET(Y, 0, k - 2) - 0.5 * STEP * (MATRIX_GET(Y, 1, k - 1) + MATRIX_GET(Y, 1, k - 2)));
        MATRIX_SET(D, 1, jsf - 1,
                   MATRIX_GET(Y, 1, k - 1) - MATRIX_GET(Y, 1, k - 2) - STEP * ((1.0 / powl(1.0 - xk * xk, 2.0)) * (2.0 * (kp - km * (1.0 - xk) + xk + kp * xk) * (1.0 - xk * xk) * y1t + y0t * (kp + m * m - s - km * km * powl(1.0 - xk, 2.0) + 2.0 * c * s * xk + 2.0 * m * s * xk + xk * xk * (-c * c - kp + s * (1.0 + s) - 2.0 * c * s * xk + c * c * xk * xk) - kp * kp * powl(1.0 + xk, 2.0) + km * (1.0 + 2.0 * kp) * (1.0 - xk * xk) - y2t + xk * xk * y2t))));
        MATRIX_SET(D, 2, jsf - 1,
                   MATRIX_GET(Y, 2, k - 1) - MATRIX_GET(Y, 2, k - 2));
    }
}