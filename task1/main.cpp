#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>

float formJacobi(float **alpha, float *x, float *x1, float *beta, int n, int core);

int methodJacobi(float **A, float *b, float *x, int n, float epsilon, int core);

int main() {
//    omp_set_dynamic(0);
    float **A, *b, *x;
    float epsilon = 1e-6;

    int cores[] = {1, 2, 3, 4, 5, 6, 7, 8};
    int size_matrix[] = {500, 1000, 5000,  10000,15000};
    int i, j;

    for (auto size : size_matrix) {
        std::cout << "N = " << size << ":" << std::endl;
        for (auto core : cores) {
            omp_set_num_threads(core);
            A = new float *[size];
            for (i = 0; i < size; i++)
                A[i] = new float[size];
            b = new float[size];
            x = new float[size];
            for (i = 0; i < size; A[i][i] = 1, i++)
                for (j = 0; j < size; j++)
                    if (i != j)
                        A[i][j] = 0.1 / (i + j);
            for (i = 0; i < size; i++)
                b[i] = sin(i);

            auto t1 = omp_get_wtime();
            methodJacobi(A, b, x, size, epsilon, core);
            auto t2 = omp_get_wtime();
            std::cout << t2 - t1 <<
            ' ';
            std::cout.flush();
            delete[] A;
            delete[] b;
            delete[] x;
        }
        std::cout << std::endl;
    }
    return 0;
}

float formJacobi(float **alpha, float *x, float *x1, float *beta, int n, int core) {
    int i, j;
    float s, max;
#pragma omp parallel for schedule(guided) private(i, j, s)
    for (i = 0; i < n; ++i) {
        s = 0;
        for (j = 0; j < n; ++j) {
            s += alpha[i][j] * x[j];
        }
        s += beta[i];
        if (i == 0) {
            max = fabs(x[i] - s);
        } else if (fabs(x[i] - s) > max) {
            max = fabs(x[i] - s);
        }
        x1[i] = s;
    }
    return max;
}

int methodJacobi(float **A, float *b, float *x, int n, float epsilon, int core) {
    float **alpha, *beta, *x1, max;
    int i, j, kvo;
    alpha = new float *[n];
    for (i = 0; i < n; ++i)
        alpha[i] = new float[n];
    beta = new float[n];
    x1 = new float[n];
#pragma omp parallel for schedule(guided)
    for (i = 0; i < n; ++i) {
#pragma omp parallel for schedule(guided)
        for (j = 0; j < n; ++j) {
            if (i == j)
                alpha[i][j] = 0;
            else
                alpha[i][j] = -A[i][j] / A[i][i];
        }
        beta[i] = b[i] / A[i][i];
    }
    for (i = 0; i < n; i++)
        x1[i] = beta[i];
    kvo = 0;
    max = 5 * epsilon;
    while (max > epsilon) {
        for (i = 0; i < n; ++i)
            x[i] = x1[i];
        max = formJacobi(alpha, x, x1, beta, n, core);
        kvo++;
    }
    delete[] alpha;
    delete[] beta;
    delete[] x1;
    return kvo;
}