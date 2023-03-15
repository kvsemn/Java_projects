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
    public static double[][] multiply(double [][] A, double[][] B) {
        int n = A.length;
        double [][] C = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                for (int k = 0; k < n; k++) {
                    C[i][j] += A[i][k] * B[k][j];
                }
            }
        }
        return C;
    }

    public static double[][] transpose(double[][] matrix) {
        int n = matrix.length;
        double[][] transpose = new double[n][n];

        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                transpose[i][j] = matrix[j][i];

        return transpose;

    }

    public static void Yakobi(double [][] A, double eps) {
        double t_A=1;
        int iterations = 0;

        double [][] U_k = new double[A.length][A.length];
        for (int l=0; l< U_k.length; l++) {
            for (int m=0; m< U_k.length; m++) {
                if (l==m) U_k[l][m] = 1; //единичная матрица
            }
        }

        while (t_A > eps){
            //максимальный по модулю внедиагональный элемент матрицы
            double maxA=0;
            int i=0, j=0;
            for (int l=0; l< A.length; l++) {
                for (int m=0; m< A.length; m++) {
                    if(l<m) {
                        double temp = Math.abs(A[l][m]);
                        if (maxA < temp) {
                            maxA = temp;
                            i = l; j = m;
                        }
                    }
                }
            }

            //находим соответствующую матрицу вращения
            double phi = 0;
            if (A[i][i] == A[j][j]) phi = Math.PI/4;
            else phi = Math.atan(2*maxA/(A[i][i]-A[j][j]))/2;
            double [][] U = new double[A.length][A.length];
            U[i][j] = -Math.sin(phi);
            U[j][i] = Math.sin(phi);
            U[i][i] = Math.cos(phi);
            U[j][j] = Math.cos(phi);
            for (int m=0; m<U.length; m++)
                if(m!= i && m!=j) U[m][m] =1;

            // A_k
            double[][] U_t = transpose(U);
            A = multiply(multiply(U_t, A), U);

            //Собственные векторы (*=U_k) - произведение матриц на k-ой итерации
            U_k = multiply(U_k, U);

            //t(A_k)
            /*for (int l=0; l<A.length; l++) {
                for(int m=0; m<A.length; m++) {
                    if(l<m) {
                        t_A += Math.pow(A[l][m],2);

                    }
                }
            }*/
            t_A = Math.pow(A[0][1], 2) + Math.pow(A[0][2], 2) + Math.pow(A[1][2], 2);
            //t_A = Math.round(Math.sqrt(t_A));
            t_A = Math.sqrt(t_A);
            //System.out.println(t_A);
            iterations++;
        }
        //искомые собственные значения (диагональные элементы матрицы A_k)
        double [] lymbda = new double[A.length];
        System.out.println("Искомые собственные значения:");
        for (int l=0; l< A.length; l++) {
            lymbda[l] = A[l][l];
            System.out.println(lymbda[l]);
        }
        System.out.println();
        //выводим собственные векторы
        System.out.println("Cобственные векторы:");
        double[] X = new double[U_k.length];
        for (int l=0; l< U_k.length; l++) {
            System.out.println("X"+(l+1)+": ");
            for (int m=0; m< U_k.length; m++) {
                System.out.println(U_k[m][l]);
            }
        }
        System.out.println("Задача была решена за " + iterations + " итерации.");
    }

    public static void main(String[] args) {
	// write your code here
        double [][] A = {{4,2,1}, {2,5,3}, {1,3,6}}; //матрица из методички
        //double [][] A = {{0,-7,7}, {-7,-9,-5}, {7,-5,-1}};
        double eps = 0.01;
        Yakobi(A, eps);
    }
}
