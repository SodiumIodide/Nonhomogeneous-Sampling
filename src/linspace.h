#ifndef _LINSPACE_H
#define _LINSPACE_H

void linspace(double *arr_x, double x_start, double x_fin, int x_len) {
    double delta_x = (x_fin - x_start) / (x_len - 1);

    for (int i = 0; i < x_len; i++) {
        *(arr_x + i) = x_start + i * delta_x;
    }
}

#endif
