#include <iostream>
#include <mpi.h>
#include <vector>
#include <math.h>

float formJacobi(std::vector<std::vector<double> >& alpha, std::vector<double>&x,
                 std::vector<double>&x1, std::vector<double>&beta, int n, int rank, int num_procs);

std::vector<double> methodJacobi(std::vector<std::vector<double> >& A, std::vector<double>&b,
                                 std::vector<double>&x, int n, float epsilon, int rank, int num_procs);

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    int rank, num_procs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    float epsilon = 1e-6;

    int size = 10000;
    MPI_Bcast(&size, 1, MPI_INT, 0, MPI_COMM_WORLD);

    std::vector<std::vector<double> > A(size, std::vector<double>(size));
    std::vector<double> b(size);
    std::vector<double> x(size);
    for (int i = 0; i < size; A[i][i] = 1, i++)
        for (int j = 0; j < size; j++)
            if (i != j)
                A[i][j] = 0.1 / (i + j);
    for (int i = 0; i < size; i++)
        b[i] = sin(i);

    for (int i = 0; i < size; i++) {
        MPI_Bcast(&A[i][0], size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
    MPI_Bcast(&b[0], size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double start = MPI_Wtime();
    x = methodJacobi(A, b, x, size, epsilon, rank, num_procs);
    double finish = MPI_Wtime();
    std::cout << finish - start << std::endl;
    std::cout.flush();

    std::cout << std::endl;

    MPI_Finalize();
    return 0;
}

float formJacobi(std::vector<std::vector<double> >& alpha, std::vector<double>&x,
                 std::vector<double>&x1, std::vector<double>&beta, int n, int rank, int num_procs) {

    std::vector<int> prefix(num_procs + 1);
    prefix[0] = int(n / num_procs);
    for (int i = 1; i <= num_procs; i++)
        prefix[i] = prefix[i - 1] + int(n / num_procs);
    prefix[num_procs] += int(n % num_procs);

    int rank_start = prefix[rank - 1];
    int rank_stop = prefix[rank];
    float s, max;
    for (int i = rank_start; i < rank_stop; ++i) {
        s = 0;
        for (int j = 0; j < n; ++j) {
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

//    MPI_Allreduce(MPI_IN_PLACE, &max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    return max;
}

std::vector<double> methodJacobi(std::vector<std::vector<double> >& A, std::vector<double>&b,
                                 std::vector<double>&x, int n, float epsilon, int rank, int num_procs) {
    float max, norm;

    std::vector<int> prefix(num_procs + 1);
    prefix[0] = int(n / num_procs);
    for (int i = 1; i <= num_procs; i++)
        prefix[i] = prefix[i - 1] + int(n / num_procs);
    prefix[num_procs] += int(n % num_procs);

    int rank_start = prefix[rank - 1];
    int rank_stop = prefix[rank];

    std::vector<std::vector<double> > alpha(n, std::vector<double>(n));
    std::vector<double> beta(n);
    std::vector<double> x1(n);
    std::vector<double> result_x(n);

    for (int i = rank_start; i < rank_stop; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i == j)
                alpha[i][j] = 0;
            else
                alpha[i][j] = -A[i][j] / A[i][i];
        }
        beta[i] = b[i] / A[i][i];
    }
    for (int i = rank_start; i < rank_stop; i++)
        x1[i] = beta[i];

    max = 5 * epsilon;

    while (max > epsilon) {

        for (int i = rank_start; i < rank_stop; ++i)
            x[i] = x1[i];
//        norm = formJacobi(alpha, x, x1, beta, n, rank, num_procs);
        max = formJacobi(alpha, x, x1, beta, n, rank, num_procs);
//        MPI_Allreduce(&norm, &max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
//        MPI_Bcast(&max, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Allgather(&x[rank_start], prefix[rank], MPI_DOUBLE,
                      &result_x[rank_start], n,MPI_DOUBLE, MPI_COMM_WORLD);
    }

    return result_x;
}