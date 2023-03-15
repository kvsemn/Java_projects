package com.company;

public class Main {

    public static double[] run_through(double[][] A, double[] B) {
        int n = A.length;
        double[] x = new double[n];

        double [] a = new double[n];
        double [] b = new double[n];
        double [] c = new double[n];

        double [] p = new double[n];
        double [] q = new double[n];

        a[0] = 0;
        b[0] = A[0][0];
        c[0] = A[0][1];
        c[n-1] = 0;
        b[n-2] = A[n-2][n-2];
        b[n-1] = A[n-1][n-1];
        for (int i = 1; i < A.length-1; i++) {
            for (int j = 0; j < A.length; j++) {
                if (A[i][j]!=0) {
                    a[i] = A[i][j];
                    b[i] = A[i][j+1];
                    c[i] = A[i][j+2];
                    j=j+2;
                }
            }
        }

        /*for (int i=0; i< A.length; i++) {
            if (Math.abs(b[i]) > Math.abs(a[i])+Math.abs(c[i])) {
                System.out.println("Достаточное условие выполнено");
            }
        }*/

        p[0]=-c[0]/b[0];
        q[0]=B[0]/b[0];
        for (int i = 1; i < A.length; i++) {
            p[i] = -c[i]/(b[i]+a[i]*p[i-1]);
            q[i] = (B[i]-a[i]*q[i-1])/(b[i]+a[i]*p[i-1]);
        }
        for (int i=(n-1); i>=0; i--) {
            if(i == n-1) x[i] = q[i];
            else x[i] = p[i] * x[i+1] + q[i];
        }
        return x;
    }

    public static void main(String[] args) {
	// write your code here
        double [][] A = {{8,4,0,0,0}, {-5,22,8,0,0}, {0,-5,-11,1,0}, {0,0,-9,-15,1}, {0,0,0,1,7}};
        double [] B = {48, 125, -43, 18, -23};

        System.out.println("Решение СЛАУ:");
        double [] X = run_through(A,B);
        for (int i=0; i<X.length; i++) {
            System.out.println("X"+(i+1)+": "+X[i]);
        }
    }
}
