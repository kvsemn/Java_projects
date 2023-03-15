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

    private static double sign(double x) {
        if (x > 0)
            return 1;
        else if (x < 0)
            return -1;
        return 0;
    }

    public static double[][] multiplyV(double [] matrix) {
        double C[][] = new double[matrix.length][matrix.length];
        for (int i=0; i< matrix.length; i++) {
            for (int j=0; j< matrix.length; j++) {
                C[i][j] += Math.round(matrix[i] * matrix[j]*100.0)/100.0;
            }
        }
        return C;
    }
    public static double multiply_rowV(double [] matrix) {
        double C = 0;
        for (int i=0; i< matrix.length; i++) {
            C += Math.round(Math.pow(matrix[i], 2) *100.0)/100.0;
        }
        return C;
    }
    public static double[][] multiply(double [][] A, double [][] B) {
        int n = A.length;
        double [][] C = new double[A.length][A.length];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                for (int k = 0; k < n; k++) {
                    C[i][j] += Math.round(A[i][k] * B[k][j] * 100.0)/100.0;
                }
            }
        }
        return C;
    }

    public static void QR(double[][] A, double eps){
        //единичная матрица
        double [][] E = new double[A.length][A.length];
        for (int i=0; i<A.length; i++) {
            for (int j=0; j<A.length; j++) {
                if(i==j) E[i][j] = 1;
            }
        }
        int iterations = 0;
        double [][] A_k;
        double [] lymbda = new double[A.length];
        double [] lymbda_prev = new double[A.length];
        double eps_k = 1;
        boolean run = true;

        while(run) {

            //ортогональная матрица
            double [][] Q = E;

            for (int n = 0; n < A.length - 1; n++) {
                double[] v = new double[A.length];
                double temp = 0;
                //находим вектор v
                for (int j = n; j < A.length; j++)
                    temp += Math.pow(A[j][n], 2);
                for (int i = 0; i < A.length; i++) {
                    if (n > 0 && i == 0) v[n - 1] = 0;
                    else if (i == n) v[i] = A[n][n] + sign(A[n][n]) * Math.sqrt(temp);
                    else v[i] = A[i][n];
                }

                //матрица Хаусхолдера
                double[][] H = new double[A.length][A.length];
                double[][] v_k = multiplyV(v);
                double v_t = multiply_rowV(v);
                for (int i = 0; i < H.length; i++) {
                    for (int j = 0; j < H.length; j++) {
                        H[i][j] = E[i][j] - 2 * (v_k[i][j] / v_t);
                    }
                }

                //получение матрицы A_k (по завершении цикла - матрица R)
                A = multiply(H, A);

                //получение ортогональной матрицы
                Q = multiply(Q, H);
            }

            //Проверка завершения предыдущего цикла
            if (iterations == 0) {
                A_k = multiply(Q,A);
                System.out.println("QR-разложение для матрицы A_0");
                System.out.println("Матрица Q_0:");
                show_matrix(Q);
                System.out.println("\nМатрица R_0:");
                show_matrix(A);
                System.out.println("\nМатрица Q_0 x R_0:");
                show_matrix(A_k);
                System.out.println();
            }

            //матрица A_k (R_k x Q_k)
            A_k = multiply(A, Q);
            //округление элементов матрицы до сотых
            for (int i=0; i< A_k.length; i++) {
                for (int j=0; j< A_k.length; j++) {
                    A_k[i][j] = Math.round(A_k[i][j] * 100.0) / 100.0;
                }
            }

            //нахожу решения уравнения
            int n = A_k.length-2;

            for (int j=0; j< A_k.length; j++) {

                double D = Math.pow((A_k[n-1][n-1] + A_k[n + 1][n + 1]), 2)
                        - 4 * (A_k[n][n + 1] * A_k[n + 1][n]);

                if(j==0) lymbda[j] = A_k[0][0];
                else{
                    lymbda[j] = (A_k[n][n] + A_k[n + 1][n + 1] + Math.sqrt(D))/2;
                    lymbda[j+1] = (A_k[n][n] + A_k[n + 1][n + 1] - Math.sqrt(D))/2;
                    j+=1;
                }
            }

            for (int j=0; j<lymbda.length; j++) {

                if (Math.abs(lymbda[j] - lymbda_prev[j]) < eps || Math.sqrt(eps_k) < eps)
                    run = false;
            }

            for (int j=0; j<lymbda.length; j++) {
                lymbda_prev[j] = lymbda[j];
                for (int i=0; i <lymbda.length; i++){
                    eps_k += Math.pow(A_k[i][i], 2);
                }
            }
            //для продолжения итерационного процесса
            A = A_k;
            iterations++;
        }
        //Вывод данных
        System.out.println("Количество проделанных итераций: " + iterations);
        System.out.println("Решение задачи:");
        for (int j=0; j<lymbda.length; j++) {
            System.out.println(lymbda[j]);
        }

    }

    public static void main(String[] args) {
	// write your code here
        //double [][] A = {{1,3,1}, {1,1,4}, {4,3,1}};
        double [][] A = {{5,8,-2}, {7,-2,-4}, {5,8,-1}};
        double eps = 0.01;
        QR(A,eps);

    }
}
