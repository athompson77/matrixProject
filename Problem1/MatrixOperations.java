import java.util.*;
import java.text.*;

public class MatrixOperations {
    
    public static void printMatrix(double[][] A) {
        int n = A.length;
        System.out.println();
        DecimalFormat fmt = new DecimalFormat("0.####");
        fmt.setMinimumFractionDigits(4);
        for (int i=0; i< n; i++) {          
            System.out.print("[ ");
            for (int j=0; j < n; j++) {
                System.out.print("{");
                System.out.print(fmt.format(A[i][j]));
                System.out.print("} ");
            }           
            System.out.print("]");
            System.out.println();
            System.out.println();
        }
    }

    public static double norm(double[] vector1) {
        double result = 0;
        for (int i=0; i<vector1.length; i++) {
            result = result + (vector1[i] * vector1[i]);      
        }
        result = Math.sqrt(result);
        return result;
    }

    public static double norm2D(double[][] vector1) {
        double result = 0;
        int n = vector1.length;
        for (int i = 0; i < n; i++) {
            result = result + (vector1[i][0] * vector1[i][0]);
        }
        result = Math.sqrt(result);
        return result;
    }

    public static double[][] transposeMatrix(double[][] A) {
        int n = A.length;
        double[][] temp = new double[n][n];
        for (int i=0; i<n; i++) {
            for (int j=0; j<n; j++) {
                temp[i][j] = A[i][j];              
            }
        }
        for (int i=0; i<n; i++) {
            for (int j=0; j<n; j++) {
                if (i <= j) {
                    temp[i][j] = A[j][i];
                }
                if (i < j) {
                    temp[j][i] = A[i][j];
                }
            }
        }
        return temp;
    }

    public static double[][] multiplyTranspose(double[][] vector1, double[][] vector2) {
        int n = vector1.length;
        double[][] Matrix = new double[(int)n][(int)n];
        for(int i=0; i<n; i++){
            for(int j=0; j<n; j++){
                for(int k=0; k<1; k++){
                    Matrix[i][j] += vector1[i][k]*vector2[k][j];
                }
            }
        }
        return Matrix;
    }

    public static double[] subtractVector(double[] vector1, double[] vector2) {
        double n = vector1.length;
        double[] result = new double[(int)n];
        for (int i = 0; i<n; i++) {
            result[i] = vector1[i] - vector2[i];
        }
        return result;  
    }

    public static double[][] subtractMatrix(double[][] A, double[][] B) {
        int n = A.length;
        double[][] result = new double[(int)n][(int)n];
        for (int i = 0; i<n; i++) {
            for (int j=0; j<n; j++) {
                result[i][j] = A[i][j] - B[i][j];
            }
        }
        return result;  
    }

    public static double normOfInfinity(double[][] A) {
        int n = A.length;
        double row = 0;
        for (int i = 0; i < n; i++) {
            double row2 = 0;
            for (int j = 0; j < n; j++) {
                row2 = row2 + Math.abs(A[i][j]);
            }
            row = Math.max(row,row2);
        }
        return row;
    }

    public static double[][] matrixMultiplication(double[][] A, double[][] B) {
        int m = A.length;
        int n = A[0].length;
        int o = B.length;
        int p = B[0].length;

        double[][] product = new double[m][p];
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < p; j++) {
                product[i][j] = 0.00000;
            }
        }

        for (int i = 0; i<m; i++) {
            for (int j = 0; j < p; j++) {
                for (int k = 0; k < n; k++) {
                product[i][j] += A[i][k] * B[k][j];
                }
            }
        }
        return product;
    }

    public static double[][] scalarMultiplication(double[][] Matrix1, double scalar) {
        int n = Matrix1.length;
        double result[][] = new double[(int)n][(int)n];   
            for(int i = 0; i < Matrix1.length; i++) {
                for(int j = 0; j < Matrix1.length; j++) {         
                    result[i][j] += Matrix1[i][j]*scalar;               
                }  
            }  
            return result;
    }

    public static void print3DMatrix(double[][][] A) {
        int n = A.length;
        System.out.println();
        DecimalFormat fmt = new DecimalFormat("0.####");
        fmt.setMinimumFractionDigits(4);
        for (int i=0; i< n-1; i++) {
            int h = i+1;
            System.out.println("H" + h + ":");
            for (int j=0; j < n; j++) {
                System.out.print("[");
                for (int k=0; k<n; k++) {
                    System.out.print("{");
                    System.out.print(fmt.format(A[j][k][i]));
                    System.out.print("} ");
                }
                System.out.println("]");
            }           
            System.out.println();       
        }
    } 
}