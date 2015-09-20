import java.util.Scanner;
import java.text.*;

public class lu_fact extends MatrixOperations {
 
    public static double[][] calculateLU(double[][] A, String returned) {
 
        DecimalFormat fmt = new DecimalFormat("0.####");
 
        int n = A.length;
        System.out.println();
        double[][] hil = new double[n][n];
        double[][] HCopy = new double[n][n];
 

        // VECTOR B:
        double[] bVector = new double[n];
        for (int i=0; i<n; i++) {         
                double Bi = (Math.pow(.1, ((double)n)/3));
                bVector[i] = Bi;
        }

 
        //THE IDENTITY MATRIX:
        System.out.println();
        double[][] I = new double[n][n];
        for (int i=0; i<n; i++){
            for (int j=0; j<n; j++){
                if (i==j){
                    I[i][j] = 1;
                }
                else{
                    I[i][j] = 0;
                }
            }
        }
     
     
        // ROW-REDUCING THE HILBERT MATRIX:
        // FINDING U:
         
         
        double[][] L = findL(HCopy);
        double[][] U = findU(hil);
        double[][] LU = matrixMultiplication(L, U);
         
         
        /**System.out.println("L=");
        *printMatrix(L);
        *System.out.println("U=");
        *printMatrix(U);*/

        double[][] errorMatrix = subtractMatrix(LU, hil);
        double errorLU = normOfInfinity(errorMatrix);
        /**System.out.println("The error of [LU - H]:");
        *System.out.println(errorLU); */

        double[][] errorValue = new double[1][1];
        errorValue[0][0] = errorLU;       
    
        if (returned.equals("L")) {
            return L;
        } else if (returned.equals("U")) {
            return U;
        } else if (returned.equals("E")) {
            return errorValue;
        } else if (returned.equals("H")) {
            return hil;
        } else {
            return null;
        }
    } 

     
    public static double[][] findL( double[][] Matrix) {
        int n = Matrix.length;
        double[][] L = new double[n][n];
        for (int i=0; i<n; i++) {
            L[i][0] = (Matrix[i][0]/Matrix[0][0]);
        }
         
        double[][] tempL = new double[n][n];
        double[][] I = new double[n][n];
        for (int id=0; id<n; id++){
            for (int jd=0; jd<n; jd++){
                if (id==jd){
                    I[id][jd] = 1;
                }
                else{
                    I[id][jd] = 0;
                }
            }
        }
        for (int copy1=0; copy1<n; copy1++) {
            for (int copy2=0; copy2<n; copy2++) {   
                tempL[copy1][copy2] = I[copy1][copy2];
            }
        }
         
        for (int j=0; j < (n-1); j++) { 
            for (int i=j+1; i < n; i++) {
                double a = Matrix[i][j]/Matrix[j][j];
                Matrix[i][j] = 0;
                for (int k=j+1; k < n; k++) {
                    Matrix[i][k] = Matrix[i][k] - a*Matrix[j][k];
                }
            }
            for (int l=j+1; l<n-1; l++) {
                tempL[l+1][j+1]=(double)Matrix[l+1][j+1]/(double)Matrix[j+1][j+1];          
            }
            L = matrixMultiplication(L, tempL);
        }   
        return L;
    }
 
 
    public static double[][] findU( double[][] Matrix) {
        int n = Matrix.length;
        for (int j=0; j < (n-1); j++) { //<==looking at the jth column of the matrix
            for (int i=j+1; i < n; i++) {//subtracting row j
                double a = Matrix[i][j]/Matrix[j][j];
                Matrix[i][j] = 0;
                for (int k=j+1; k < n; k++) {
                    Matrix[i][k] = Matrix[i][k] - a*Matrix[j][k];
                }
            }
        }
        return Matrix;
    }
     
    public static double[][] calculateLUHIL(int dimension, String returned) {
 
        DecimalFormat fmt = new DecimalFormat("0.####");
 
        /**Scanner scan = new Scanner (System.in);
        System.out.println("Dimensions of nxn Hilbert matrix n?");*/
        int n = dimension;
        System.out.println();
        double[][] hil = new double[n][n];
        double[][] HCopy = new double[n][n];
 
        // Construction of the Hilbert Matrix, our starting point:
        //System.out.println();
        for (int i=0; i<n; i++) {           
            for (int j=0; j<n; j++) {
                double Hij = (1/(((double)i+1)+((double)j+1)-1));
                hil[i][j] = Hij;
                HCopy[i][j] = Hij;
                fmt.setMinimumFractionDigits(4);
            }           

        }
        // VECTOR B:
        double[] bVector = new double[n];
        for (int i=0; i<n; i++) {         
                double Bi = (Math.pow(.1, ((double)n)/3));
                bVector[i] = Bi;
        }

 
        //THE IDENTITY MATRIX:
        //System.out.println();
        double[][] I = new double[n][n];
        for (int i=0; i<n; i++){
            for (int j=0; j<n; j++){
                if (i==j){
                    I[i][j] = 1;
                }
                else{
                    I[i][j] = 0;
                }
            }
        }
     
     
        // ROW-REDUCING THE HILBERT MATRIX:
        // FINDING U:
         
         
        double[][] L = findL(HCopy);
        double[][] U = findU(hil);
        double[][] LU = matrixMultiplication(L, U);
         
         
        /**System.out.println("L=");
        printMatrix(L);
        System.out.println("U=");
        printMatrix(U);*/

        double[][] errorMatrix = subtractMatrix(LU, hil);
        double errorLU = normOfInfinity(errorMatrix);
        /**System.out.println("The error of [LU - H]:");
        System.out.println(errorLU); */

        double[][] errorValue = new double[1][1];
        errorValue[0][0] = errorLU;       
    
        if (returned.equals("lower")) {
            return L;
        } else if (returned.equals("upper")) {
            return U;
        } else if (returned.equals("error")) {
            return errorValue;
        } else if (returned.equals("Hilbert")) {
            return hil;
        } else {
            return null;
        }
        
    }
} //end of class