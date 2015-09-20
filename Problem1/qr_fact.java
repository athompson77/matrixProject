import java.text.*;
import java.util.*;
 
public class qr_fact extends MatrixOperations {
 
    public static double[][] hhCalculateQR(double[][] aMatrix, String returned) {
 
        DecimalFormat fmt = new DecimalFormat("0.####");

        int n = aMatrix.length;
        System.out.println();
        double[][] hil = aMatrix;
        double[][] hil2 = aMatrix;   
        double[][] Q = new double[n][n];
        double[][] R = new double[n][n];
         
        //Starting Q out as the identity:
        for (int i=0; i<n; i++){
            for (int j=0; j<n; j++){
                if (i==j){
                    Q[i][j] = 1;
                }
                else{
                    Q[i][j] = 0;
                }
            }
        }
         
        //Creating the vector b, given in the problem statement.
        double[] b = new double[n];
        for (int i=0; i<n; i++) {           
 
            double Bi = (Math.pow(.1, ((double)n)/3));
            b[i] = Bi;             
        }     
         
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
         
     
        double[][] hCopy = hil;
        for (int i=0; i<n; i++) {
            for (int j=0; j<n; j++) {
                hCopy[i][j] = hil[i][j];
            }
        }
         
        /** This code creates a vector of matrices so that Hn does not get lost while Householder
        * is calculating the reflection matrices. For example, as we proceed through the algorithm,
        * Hn will be modified in each step along the way, finally resulting in Q. Then if we would 
        * like to multiply Hn...H3*H2*H1 to find R, we won't be able to if Hn is not stored somewhere,
        * so I used this vector of matrices to store each Hn along the way in order.
        */
     
        double[][][] VectorofHs = new double[n][n][n];
        double[][] A = new double[n][n];
        for (int copy1=0; copy1<n; copy1++) {
            for (int copy2=0; copy2<n; copy2++) {   
                    A[copy1][copy2] = hil[copy1][copy2];
            }
        }
         
        // Beginning the Householder reflections:
         
        //For matrices ai:
        for (int i=0; i<n; i++) {       
            double[] X = new double[n];
            //For vector X:
            for (int x=0; x<n; x++) {
                X[x] = hCopy[x][i];
                if (x<i) {
                    X[x] = 0;
                }
            }
             
            int y = i+1;
             
            double[] C = new double[n];
            for (int c=0; c<n; c++) {
                if (c==i) {
                    C[i] = norm(X);
                }
                else {
                    C[c] = 0;
                }
            }
             
         
            if (i!= n-1) {
                double[][] U = new double[n][1];
                double[] U1 = subtractVector(X, C);
         
                //Turning U into a 2D array vector:
                for (int l=0; l<n; l++) {
                    U[l][0] = U1[l];
                }
                double[][] Ut = new double[1][n];
                double[] u1 = subtractVector(X, C);
                //Turning Ut into a 2D array vector:
                for (int m=0; m<n; m++) {
                    Ut[0][m] = u1[m];
                }
                double[][] UUt = multiplyTranspose(U, Ut);
                double TwoDivNorm2 = 2/(norm2D(U)*norm2D(U));
                double[][] TwoDivNormTimesUUt = scalarMultiplication(UUt, TwoDivNorm2);
                double[][]Han = subtractMatrix(I, TwoDivNormTimesUUt);                          
         
                 
                for (int row=0; row<n; row++) {
                    for (int column=0; column<n; column++) {
                        if (row==column && row<=(i-1)) {
                            Han[row][column] = 1;
                        } 
                    }
                }
                 
                double[][] Hn = matrixMultiplication(Han, A);
                R = matrixMultiplication(Han, A);
                 
                for (int copy1=0; copy1<n; copy1++) {
                    for (int copy2=0; copy2<n; copy2++) {   
                        A[copy1][copy2] = Hn[copy1][copy2];
                    }
                }
                 
                //Setting Hn back to the HCopy to start over:
                for (int copy1=0; copy1<n; copy1++) {
                    for (int copy2=0; copy2<n; copy2++) {   
                        hCopy[copy1][copy2] = Hn[copy1][copy2];
                    }
                }                               
                //Inserting Hn into the matrix Vector.
                for (int copy1=0; copy1<n; copy1++) {
                    for (int copy2=0; copy2<n; copy2++) {
                        VectorofHs[copy1][copy2][i] = Hn[copy1][copy2];
                    }
                }   
                 
                for (int copy1=0; copy1<n; copy1++) {
                    for (int copy2=0; copy2<n; copy2++) {   
                        R[copy1][copy2] = Hn[copy1][copy2];
                    }
                }                           
                Q = matrixMultiplication(Q, Han);
            }           
        }
         
        /** System.out.print("Q:");
        printMatrix(Q);
        System.out.println();
        System.out.print("R:");
        printMatrix(R); */
         
        double[][] Z = matrixMultiplication(Q, R);
         
        // THE HILBERT MATRIX AGAIN:
        System.out.println();
        for (int i=0; i<n; i++) {           
            for (int j=0; j<n; j++) {
                double Hij2 = (1/(((double)i+1)+((double)j+1)-1));
                hil2[i][j] = Hij2;
                fmt.setMinimumFractionDigits(4);
            }
        }  

        double[][] errorMatrix = subtractMatrix(Z, hil2);
        double errorQR= normOfInfinity(errorMatrix);
        /**System.out.println("The error of [QR - H]:");
        System.out.println(errorQR); */

        double[][] errorValue = new double[1][1];
        errorValue[0][0] = errorQR; 

        if (returned.equals("Q")) {
            return Q;
        } else if (returned.equals("R")) {
            return R;
        } else if (returned.equals("error")) {
            return errorValue;
        } else if (returned.equals("Hilbert")) {
            return hil;
        } else {
            return null;
        }
    } 
    public static double[][] hhCalculateQRHIL(int dimension, String returned) {
 
 
        DecimalFormat fmt = new DecimalFormat("0.####");
        
        int n = dimension;
        System.out.println();
        double[][] hil = new double[n][n];  
        double[][] hil2 = new double[n][n];     
        double[][] Q = new double[n][n];
        double[][] R = new double[n][n];
         
        //Starting Q out as the identity:
        for (int i=0; i<n; i++){
            for (int j=0; j<n; j++){
                if (i==j){
                    Q[i][j] = 1;
                }
                else{
                    Q[i][j] = 0;
                }
            }
        }
         
        //Creating the vector b, given in the problem statement.
        double[] b = new double[n];
        for (int i=0; i<n; i++) {           
 
            double Bi = (Math.pow(.1, ((double)n)/3));
            b[i] = Bi;             
        }       
    
 
        // THE HILBERT MATRIX:
        for (int i=0; i<n; i++) {           

            for (int j=0; j<n; j++) {
                double Hij = (1/(((double)i+1)+((double)j+1)-1));
                hil[i][j] = Hij;
                fmt.setMinimumFractionDigits(4);
            }
        }

         
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
         
     
        double[][] hCopy = hil;
        for (int i=0; i<n; i++) {
            for (int j=0; j<n; j++) {
                hCopy[i][j] = hil[i][j];
            }
        }
         
        /** This code creates a vector of matrices so that Hn does not get lost while Householder
        * is calculating the reflection matrices. For example, as we proceed through the algorithm,
        * Hn will be modified in each step along the way, finally resulting in Q. Then if we would 
        * like to multiply Hn...H3*H2*H1 to find R, we won't be able to if Hn is not stored somewhere,
        * so I used this vector of matrices to store each Hn along the way in order.
        */
     
        double[][][] VectorofHs = new double[n][n][n];
        double[][] A = new double[n][n];
        for (int copy1=0; copy1<n; copy1++) {
            for (int copy2=0; copy2<n; copy2++) {   
                    A[copy1][copy2] = hil[copy1][copy2];
            }
        }
         
        // Beginning the Householder reflections:
         
        //For matrices ai:
        for (int i=0; i<n; i++) {       
            double[] X = new double[n];
            //For vector X:
            for (int x=0; x<n; x++) {
                X[x] = hCopy[x][i];
                if (x<i) {
                    X[x] = 0;
                }
            }
             
            int y = i+1;
             
            double[] C = new double[n];
            for (int c=0; c<n; c++) {
                if (c==i) {
                    C[i] = norm(X);
                }
                else {
                    C[c] = 0;
                }
            }
             
         
            if (i!= n-1) {
                double[][] U = new double[n][1];
                double[] U1 = subtractVector(X, C);
         
                //Turning U into a 2D array vector:
                for (int l=0; l<n; l++) {
                    U[l][0] = U1[l];
                }
                double[][] Ut = new double[1][n];
                double[] u1 = subtractVector(X, C);
                //Turning Ut into a 2D array vector:
                for (int m=0; m<n; m++) {
                    Ut[0][m] = u1[m];
                }
                double[][] UUt = multiplyTranspose(U, Ut);
                double TwoDivNorm2 = 2/(norm2D(U)*norm2D(U));
                double[][] TwoDivNormTimesUUt = scalarMultiplication(UUt, TwoDivNorm2);
                double[][]Han = subtractMatrix(I, TwoDivNormTimesUUt);                          
         
                 
                for (int row=0; row<n; row++) {
                    for (int column=0; column<n; column++) {
                        if (row==column && row<=(i-1)) {
                            Han[row][column] = 1;
                        } 
                    }
                }
                 
                double[][] Hn = matrixMultiplication(Han, A);
                R = matrixMultiplication(Han, A);
                 
                for (int copy1=0; copy1<n; copy1++) {
                    for (int copy2=0; copy2<n; copy2++) {   
                        A[copy1][copy2] = Hn[copy1][copy2];
                    }
                }
                 
                //Setting Hn back to the HCopy to start over:
                for (int copy1=0; copy1<n; copy1++) {
                    for (int copy2=0; copy2<n; copy2++) {   
                        hCopy[copy1][copy2] = Hn[copy1][copy2];
                    }
                }                               
                //Inserting Hn into the matrix Vector.
                for (int copy1=0; copy1<n; copy1++) {
                    for (int copy2=0; copy2<n; copy2++) {
                        VectorofHs[copy1][copy2][i] = Hn[copy1][copy2];
                    }
                }   
                 
                for (int copy1=0; copy1<n; copy1++) {
                    for (int copy2=0; copy2<n; copy2++) {   
                        R[copy1][copy2] = Hn[copy1][copy2];
                    }
                }                           
                Q = matrixMultiplication(Q, Han);
            }           
        }
         
        /**System.out.print("Q:");
        printMatrix(Q);
        System.out.println();
        System.out.print("R:");
        printMatrix(R);*/
         
        double[][] Z = matrixMultiplication(Q, R);
         
        // THE HILBERT MATRIX AGAIN:
        System.out.println();
        for (int i=0; i<n; i++) {           
            for (int j=0; j<n; j++) {
                double Hij2 = (1/(((double)i+1)+((double)j+1)-1));
                hil2[i][j] = Hij2;
                fmt.setMinimumFractionDigits(4);
            }
        }


        double[][] errorMatrix = subtractMatrix(Z, hil2);
        double errorQR= normOfInfinity(errorMatrix);
        /** System.out.println("The error of [QR - H]:");
        System.out.println(errorQR); */

        double[][] errorValue = new double[1][1];
        errorValue[0][0] = errorQR; 

        if (returned.equals("Q")) {
            return Q;
        } else if (returned.equals("R")) {
            return R;
        } else if (returned.equals("error")) {
            return errorValue;
        } else if (returned.equals("Hilbert")) {
            return hil;
        } else {
            return null;
        } 
    }

    public static double[][] grCalculateQR(double[][] A, String returned) {
    
 
        DecimalFormat fmt = new DecimalFormat("0.####");
        
        
        int n = A.length;

        double[][] hil = A;
        double[][] Gn = new double[n][n];
        double[][] An = new double[n][n];
        double[][] Q = new double[n][n];
         

 
        // THE HILBERT MATRIX:
        System.out.println();
        for (int i=0; i<n; i++) {
 
            for (int j=0; j<n; j++) {
                double Hij = (1/(((double)i+1)+((double)j+1)-1));
                hil[i][j] = Hij;
                An[i][j] = Hij;
                fmt.setMinimumFractionDigits(4);
            }
        } 
         
 
        for (int i=0; i<n; i++){
            for (int j=0; j<n; j++){
                if (i==j){
                    Gn[i][j] = 1;
                    Q[i][j] = 1;
                }
                else{
                    Gn[i][j] = 0;
                    Q[i][j] = 0;
                }
            }
        }
         
        int iteration = 1;
        double a = An[0][n-2];
        double b = An[0][n-1];
        double cosX;
        double sinX;
         
        // The for loops that begin the Givens rotation matrices.
 
        for (int i=0; i<n; i++) {
                for (int j=(n-1); j>i; j--) {                                       
                     
                    a = An[j-1][i];
                    b = An[j][i];   
                    cosX = a/(Math.sqrt(a*a+b*b));
                    sinX = -b/(Math.sqrt(a*a+b*b));
                     
                    Gn[j][j] = cosX;
                    Gn[j][j-1] = sinX;
                    Gn[j-1][j] = -sinX;
                    Gn[j-1][j-1] = cosX;
                              
 
                    An = matrixMultiplication(Gn, An);
                           
                     
                    Q = matrixMultiplication(Gn, Q);
 
                    // Turning the Gn matrix back into the identity.
                     
                    for (int ident=0; ident<n; ident++){
                        for (int ident2=0; ident2<n; ident2++){
                            if (ident==ident2)
                                Gn[ident][ident2] = 1;
                            else
                                Gn[ident][ident2] = 0;
                        }
                    }           
                    iteration += 1;     
                }//end j    
        }//end i
         
        
         
        Q = transposeMatrix(Q);
        double[][] answer = matrixMultiplication(Q, An);

        // Calculating the maximum norm of QR-H.
         
        answer = subtractMatrix(answer, hil);
        double maxNorm = normOfInfinity(answer);
        
        double[][] errorValue = new double[1][1];
        errorValue[0][0] = maxNorm; 

        if (returned.equals("Q")) {
            return Q;
        } else if (returned.equals("R")) {
            return An;
        } else if (returned.equals("error")) {
            return errorValue;
        } else if (returned.equals("Hilbert")) {
            return hil;
        } else {
            return null;
        }    
    } //end of main

    public static double[][] grCalculateQRHIL(int dimension, String returned) {
     
        /** Asking for the number of dimensions to work with the computing factorization of the Hilbert matrix,
        * formatting the numbers and matrices, and computing/printing the Hilbert matrix:
        */
 
        DecimalFormat fmt = new DecimalFormat("0.####");

        int n = dimension;

        double[][] hil = new double[n][n];
        double[][] Gn = new double[n][n];
        double[][] An = new double[n][n];
        double[][] Q = new double[n][n];
         

 
        // THE HILBERT MATRIX:
        System.out.println();
        for (int i=0; i<n; i++) {
            for (int j=0; j<n; j++) {
                double Hij = (1/(((double)i+1)+((double)j+1)-1));
                hil[i][j] = Hij;
                An[i][j] = Hij;
                fmt.setMinimumFractionDigits(4);
            }       
        } 
         
 
        for (int i=0; i<n; i++){
            for (int j=0; j<n; j++){
                if (i==j){
                    Gn[i][j] = 1;
                    Q[i][j] = 1;
                }
                else{
                    Gn[i][j] = 0;
                    Q[i][j] = 0;
                }
            }
        }
         
        int iteration = 1;
        double a = An[0][n-2];
        double b = An[0][n-1];
        double cosX;
        double sinX;
         
        // The for loops that begin the Givens rotation matrices.
 
        for (int i=0; i<n; i++) {
                for (int j=(n-1); j>i; j--) {                                       
                     
                    a = An[j-1][i];
                    b = An[j][i];   
                    cosX = a/(Math.sqrt(a*a+b*b));
                    sinX = -b/(Math.sqrt(a*a+b*b));
                     
                    Gn[j][j] = cosX;
                    Gn[j][j-1] = sinX;
                    Gn[j-1][j] = -sinX;
                    Gn[j-1][j-1] = cosX;
                              
 
                    An = matrixMultiplication(Gn, An);
                           
                     
                    Q = matrixMultiplication(Gn, Q);
 
                    // Turning the Gn matrix back into the identity.
                     
                    for (int ident=0; ident<n; ident++){
                        for (int ident2=0; ident2<n; ident2++){
                            if (ident==ident2)
                                Gn[ident][ident2] = 1;
                            else
                                Gn[ident][ident2] = 0;
                        }
                    }           
                    iteration += 1;     
                }//end j    
        }//end i
         
         
        Q = transposeMatrix(Q);
        double[][] answer = matrixMultiplication(Q, An);

        // Calculating the maximum norm of QR-H.
         
        answer = subtractMatrix(answer, hil);
        double maxNorm = normOfInfinity(answer);
        
        double[][] errorValue = new double[1][1];
        errorValue[0][0] = maxNorm; 

        if (returned.equals("Q")) {
            return Q;
        } else if (returned.equals("R")) {
            return An;
        } else if (returned.equals("error")) {
            return errorValue;
        } else if (returned.equals("Hilbert")) {
            return hil;
        } else {
            return null;
        } 
    } //end of main
} 