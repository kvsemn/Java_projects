package com.company;

public class Main {
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

    public static double norma_row(double[] matrix) {
        double temp = 0.0;
        double norm = 0.0;
        for (int i = 0; i < matrix.length; i++) {
            temp += Math.abs(matrix[i]);
            if (norm < temp) norm = temp;
            temp = 0.0;
        }
        return norm;
    }

    public static void show_row(double[] row) {
        for (int i = 0; i < row.length; i++) {
            if (row[i] < 0) {
                System.out.print(row[i] + " ");
                System.out.println();
            } else {
                System.out.print(" " + row[i] + " ");
                System.out.println();
            }
        }
    }

    public static double[] multiply(double[][] A, double[] B) {
        int n = B.length;
        double[] C = new double[n];
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) {
                C[j] += A[j][k] * B[k];
            }
        }
        return C;
    }

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

    public static void iterations(double[][] alfa, double[] betta, double eps, int maxIterations) {
        int iterations=0;
        double eps_k;
        double[][] X_k = new double[(int) maxIterations][betta.length];
        //опрежеделие количества действительных итераций
        for (int i = 1; i < maxIterations; i++) {
            for (int j = 0; j < alfa.length; j++) {
                if (i == 1) {
                    double[] C = multiply(alfa, betta);
                    X_k[i][j] = betta[j] + C[j];
                } else {
                    double[] C = multiply(alfa, X_k[i - 1]);
                    X_k[i][j] = betta[j] + C[j];
                }
            }

            double temp = 0.0;
            double normX = 0.0;
            for (int r = 0; r < betta.length; r++) {
                if (i == 1) temp = Math.abs(X_k[i][r] - betta[r]);
                else
                    temp = Math.abs(X_k[i][r] - X_k[i - 1][r]);
                if (normX < temp) normX = temp;
            }

            eps_k = norma(alfa) / (1 - norma(alfa)) * normX;
            iterations++;
            if (eps_k < eps) {
                System.out.println("Решение методом простых итераций:");
                show_row(X_k[iterations]);
                System.out.println("Вычислительный процесс завершен за " + iterations + " итерации.");
                break;
            }

        }

    }


    public static void solve(double[][] A, double[] B, double[] X, double eps, int maxIterations) {
        int n = A.length;
        int iterations = 0;
        double[] prevX = new double[n];
        do {
            iterations++;
            for (int i = 0; i < n; i++) {
                prevX[i] = X[i];
                X[i] = B[i];
                for (int j = 0; j < n; j++) {
                    if (j != i) {
                        X[i] -= A[i][j] * X[j];
                    }
                }
                X[i] /= A[i][i];
            }
            if (norm(X, prevX) < eps) {
                break;
            }
        } while (iterations < maxIterations);
        // Выводим результат
        System.out.println("Решение системы уравнений методом Зейделя:");
        for (int i = 0; i < n; i++) {
            System.out.println("x" + (i + 1) + " = " + X[i]);
        }
        System.out.println("Вычислительный процесс завершен за " + (iterations-1) + " итерации.");
    }

    public static double norm(double[] X, double[] prevX) {
        double norm = 0;
        for (int i = 0; i < X.length; i++) {
            norm += Math.abs(X[i] - prevX[i]);
        }
        return norm;
    }

    public static void main(String[] args) {
	// write your code here
        double[][] A = {{20, 5, 7, 1}, {-1, 13, 0, -7}, {4, -6, 17, 5}, {-9, 8, 4, -25}};
        double[] B = {-117, -1, 49, -21};

        double[][] alfa = new double[A.length][A.length];
        double[] betta = new double[B.length];
        //Эквивалентный вид
        for (int i = 0; i < A.length; i++) {
            betta[i] = B[i] / A[i][i];
            for (int j = 0; j < A.length; j++) {
                if (i == j) alfa[i][j] = 0;
                else alfa[i][j] = -A[i][j] / A[i][i];
            }
        }
        System.out.println("Эквивалентный вид (матрица alfa):");
        show_matrix(alfa);
        System.out.println();

        //достаточное условие сходимости метода
        System.out.println("Достаточное условие сходимости метода итераций выполнено: " + norma(alfa) + "<1");
        System.out.println();
        double eps = 0.01;
        //априорная оценка необходимого количества итераций
        double maxIterations = Math.round(((Math.log(eps) - Math.log(norma_row(betta)) + Math.log(1 - norma(alfa))) / Math.log(norma(alfa))) - 1);
        // Вызываем функцию метода итераций
        iterations(alfa, betta, eps, (int)maxIterations);

        double[] X = new double[B.length];
        // Вызываем функцию решения системы уравнений методом Зейделя
        solve(A, B, X, eps, (int)maxIterations);
    }
}
