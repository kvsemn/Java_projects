package com.company;

public class Main {

    public static void show_matrix(double[][] matrix) {
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix.length; j++) {
                if (matrix[i][j] < 0)
                    System.out.print(matrix[i][j] + " ");
                else
                    System.out.print(" " + matrix[i][j] + " ");
            }
            System.out.println();
        }
    }
    public static void show_row(double[] row) {
        for (int i=0; i < row.length; i++) {
            if (row[i] < 0) {
                System.out.print(row[i] + " ");
                System.out.println();
            }
            else {
                System.out.print(" " + row[i] + " ");
                System.out.println();
            }
        }
    }

    public static void LU(double[][] A, int n, double [] B) {
        int l2 = 0;
        double [][] mu = new double[n][n];
        double [][] E = new double[n][n];
        double [][] LU = new double[A.length][A.length];

        for (int i=0; i<E.length; i++) {
            for (int j=0; j<E.length; j++) {
                if(j==i) E[i][j] = 1;
                else E[i][j] =0;
            }
        }
        for (int l1=0; l1 < A.length; l1++) {
            for (int i = l1; i < A.length; i++) {
                if (Math.abs(A[i][l2]) > Math.abs(A[l1][l2])) {
                    //выделение главного элемента
                    double [] temp = A[i];
                    A[i] = A[l1];
                    A[l1] = temp;

                    double t = B[i];
                    B[i] = B[l1];
                    B[l1] = t;

                    double[] e = E[i];
                    E[i] = E[l1];
                    E[l1] = e;

                }

            }

            for (int i=l1+1; i < A.length; i++) {
                mu[i][l2] = A[i][l2] / A[l1][l2]; //коэффициент
                B[i] = B[i] - mu[i][l2]*B[l1];
                for (int j=l2; j < A.length; j++) {
                    A[i][j] = A[i][j] - mu[i][l2] * A[l1][j];
                }
            }
            l2++;
        }

        for (int i=0; i<LU.length; i++) {
            for (int j=0; j<LU.length; j++) {
                if(i>j) LU[i][j] = mu[i][j];
                else LU[i][j] = A[i][j];
            }
        }
        System.out.println("Матрица LU-разложения: ");
        show_matrix(LU);
        System.out.println();
        /*System.out.println("Матрица U: ");
        show_matrix(A);

        double [][] L = new double[n][n];
        for (int i=0; i<L.length; i++) {
            for (int j=0; j<L.length; j++) {
                if(i<j) L[i][j] = 0;
                if(i==j) L[i][j] = 1;
                if(i>j) L[i][j] = mu[i][j];

            }
        }
        System.out.println();
        System.out.println("Матрица L: ");
        show_matrix(L);
        System.out.println();
        System.out.println("Преобразованная матрица B: ");
        show_row(B);*/

        //решение СЛАУ методом LU-разложения
        //Ux=z
        double [] X = new double[n];
        for (int i=n-1; i >= 0; i--) {
            double tmp = B[i];
            for (int j = n-1; j >= 0; j--) {
                if (i != j)
                    tmp -= A[i][j] * X[j];
            }
            X[i] = 1.0/A[i][i]*tmp;
        }
        System.out.println();
        System.out.println("Решение СЛАУ:");
        show_row(X);
        System.out.println();

    }

    public static double determinant(double[][] matrix, int n) {
        double det = 0;
        if (n == 1) {
            det = matrix[0][0];
        } else if (n == 2) {
            det = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
        } else
         {
            for (int i = 0; i < n; i++) {
                double[][] subMatrix = new double[n - 1][n - 1];
                for (int j = 1; j < n; j++) {
                    int k = 0;
                    for (int l = 0; l < n; l++) {
                        if (l == i) {
                            continue;
                        }
                        subMatrix[j - 1][k] = matrix[j][l];
                        k++;
                    }
                }
                det += Math.pow(-1, i) * matrix[0][i] * determinant(subMatrix, n - 1);
            }
        }
        return det;
    }

    public static double[][] minors(double[][] matrix) {
        int n = matrix.length;
        double[][] minors = new double[n][n];

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                // Get cofactor of matrix[i][j]
                double[][] temp = new double[n - 1][n - 1];
                getCofactor(matrix, temp, i, j, n);

                // Calculate the determinant of the cofactor
                minors[i][j] = Math.pow(-1, i+j) * determinant(temp, n - 1);
            }
        }
        return minors;
    }

    public static void getCofactor(double[][] matrix, double[][] temp, int p, int q, int n) {
        int i = 0, j = 0;

        // Looping for each element of the matrix
        for (int row = 0; row < n; row++) {
            for (int col = 0; col < n; col++) {
                // Copying into temporary matrix only those element
                // which are not in given row and column
                if (row != p && col != q) {
                    temp[i][j++] = matrix[row][col];

                    // Row is filled, so increase row index and
                    // reset col index
                    if (j == n - 1) {
                        j = 0;
                        i++;
                    }
                }
            }
        }
    }

    public static double[][] transpose(double[][] matrix) {
        int n = matrix.length;
        double[][] transpose = new double[n][n];

        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                transpose[i][j] = matrix[j][i];

        return transpose;

    }

    public static double norma(double[][] matrix) {
        double temp = 0.0;
        double norm = 0.0;
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix.length; j++) {
                temp += Math.abs(matrix[i][j]);
            }
            if (norm < temp) norm = temp;
            temp = 0.0;
        }
        return norm;
    }

    public static double[] Gaus_Method(double [][] A, double [] B) {
        double[] X = new double[B.length];
        double [][] mu = new double[A.length][A.length];
        int l2=0;
        for (int l1=0; l1 < A.length; l1++) {
            for (int i = l1; i < A.length; i++) {
                if (A[i][l2] == 0) {
                    //выделение ненулевого элемента (если он есть)
                    double [] temp = A[i];
                    A[i] = A[l1];
                    A[l1] = temp;
                    double t = B[i];
                    B[i] = B[l1];
                    B[l1] = t;
                }
            }
            for (int i=l1+1; i < A.length; i++) {
                mu[i][l2] = A[i][l2] / A[l1][l2]; //коэффициент
                B[i] = B[i] - mu[i][l2]*B[l1];
                for (int j=l2; j < A.length; j++) {
                    A[i][j] = A[i][j] - mu[i][l2] * A[l1][j];
                }
            }
            l2++;
        }
        for (int i= A.length-1; i >= 0; i--) {
            double tmp = B[i];
            for (int j = A.length-1; j >= 0; j--) {
                if (i != j)
                    tmp -= A[i][j] * X[j];
            }
            X[i] = 1.0/A[i][i]*tmp;
        }
        System.out.println("Решение СЛАУ методом Гаусса:");
        return X;
    }

        public static void main (String[]args) {

            double[][] A = {{3, -8, 1, -7}, {6, 4, 8, 5}, {-1, 1, -9, -3}, {-6, 6, 9, -4}};
            //double[][] A = {{10, 1, 1}, {2, 10, 1}, {2, 2, 10}};
            double[] B = {96, -13, -54, 82};
            //double[] B = {12, 13, 14};
            int n = 4;
            System.out.println("Исходная матрица А:");
            show_matrix(A);

            double detA = determinant(A,n);
            System.out.println("Определитель матрицы A: " + detA);
            System.out.println("Обратная матрица A: ");
            double[][] tr_A = transpose(A);
            double[][] minors = minors(tr_A);
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    minors[i][j] = 1/detA * minors[i][j];
                }
            }
            show_matrix(minors);
            System.out.println();

            //Проверка обусловленности матрицы A
            double cond = norma(A) * norma(minors);
            System.out.println("Обусловленность матрицы A: " + cond);
            System.out.println();

            //решение слау методом Гаусса
            double [] X = Gaus_Method(A, B);
            show_row(X);

            //LU-разложение
            double[][] A1 = {{3, -8, 1, -7}, {6, 4, 8, 5}, {-1, 1, -9, -3}, {-6, 6, 9, -4}};
            double[] B1 = {96, -13, -54, 82};
            LU(A1, n, B1);
        }

}
