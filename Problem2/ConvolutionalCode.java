import Jama.Matrix;
import java.util.Scanner;
import java.util.Random;
import java.util.ArrayList;

public class ConvolutionalCode {
    
    public static void main(String[] args) {
        Scanner keyboard = new Scanner(System.in);

        boolean b = true;
        while(b) {
            System.out.println("Type 'e' to encode a random x stream, 'd' to decode your own y stream");
            String ans = keyboard.nextLine();
            if (ans.equals("e")) {
                Matrix y = encode();
                System.out.println("Encoded result:");
                y.print(0,0);

                System.out.println("Would you like to re-decode this using iterative methods? (y/n)");
                String choice = keyboard.nextLine();            
                if (choice.equals("y")) {
                    decode(y);
                }
                b = false;
            }

            else if (ans.equals("d")) {
                decode();
                b = false;
            }

            else {System.out.println("Invalid input. Try again.");}
        }
    }

//--------------------------------------------------------------------------
    private static Matrix encode() {
        Scanner keyboard = new Scanner(System.in);
        
        System.out.println();
        System.out.println("+------------------+");
        System.out.println("| Encoding Problem |");
        System.out.println("+------------------+");
        System.out.println();
        System.out.println("Input length (4 or above) of random binary stream x:");
        int n = keyboard.nextInt();
        keyboard.nextLine();
        Matrix x = new Matrix(n, 1);

        //CREATES RANDOM BINARY VECTOR OF LENGTH N
        for (int i = 0; i < n; i = i + 1) {
            x.set(i, 0, new Random().nextInt(2));
        }
        System.out.println();
        System.out.println("Created x stream:");
        x.print(0,0);

        //ENCODES X INTO Y0 AND Y1 BASED ON FORMULA
        int ylength = n+3;
        Matrix y0 = new Matrix(ylength, 1); 
        for (int i = 0; i < ylength; i = i + 1) {
            double j00, j02, j03;
            if (i < n) {
                j00 = x.get(i,0);
            } else {
                j00 = 0;
            }
            if (i<=(n+1) && i>1) {
                j02 = x.get((i-2),0);
            } else {
                j02 = 0;
            }
            if (i<=(n+2) && i>2) {
                j03 = x.get((i-3),0);
            } else {
                j03 = 0;
            }
            double value0 = (j00 + j02 + j03) % 2;
            y0.set(i, 0, value0);
        }
        Matrix y1 = new Matrix(ylength, 1); 
        for (int i = 0; i < ylength; i = i + 1) {
            double j10, j11, j13;
            if (i < n) {
                j10 = x.get(i,0);
            } else {
                j10 = 0;
            }
            if (i<=(n) && i>0) {
                j11 = x.get((i-1),0);
            } else {
                j11 = 0;
            }
            if (i<=(n+2) && i>2) {
                j13 = x.get((i-3),0);
            } else {
                j13 = 0;
            }
            double value1 = (j10 + j11 + j13) % 2;
            y1.set(i, 0, value1);
        }

        //COMBINES Y0 AND Y1 INTO A SINGLE VECTOR
        Matrix y = new Matrix(ylength, 2);
        for(int i = 0; i < ylength; i = i + 1) {
            y.set(i, 0, y0.get(i,0));
            y.set(i, 1, y1.get(i,0));
        }
        return y;
    }

//-----------------------------------------------------------------------
    private static void decode() {
        Scanner keyboard = new Scanner(System.in);

        System.out.println();
        System.out.println("+------------------+");
        System.out.println("| Decoding Problem |");
        System.out.println("+------------------+");
        System.out.println();

        System.out.println("You must input y stream (length 7 or above, 1/0s only) in following format:");
        System.out.println("10 11 01 01 10 11 10");
        System.out.println("Please input below:");
        String ans = keyboard.nextLine();
        String[] elements = ans.split(" ");
        int n = elements.length;

        //CREATES Y0 AND Y1 VECTORS FROM USER INPUTTED VECTOR IF FORMATTED CORRECTLY
        Matrix y0 = new Matrix(n, 1);
        Matrix y1 = new Matrix(n, 1);
        for (int i = 0; i < n; i = i + 1) {
            String first = elements[i].substring(0,1);
            y0.set(i, 0, Double.parseDouble(first));
            String second = elements[i].substring(1,2);
            y1.set(i, 0, Double.parseDouble(second));
        }

        //MAKES INITIAL GUESS VECTOR FROM USER INPUT IF FORMATTED CORRECTLY
        System.out.println("You must input initial guess x0 (same length as y, 1/0s only) in following format:");
        System.out.println("1 1 0 1 1 1 0");
        System.out.println("Please input below:");
        String guess = keyboard.nextLine();
        String[] digits = guess.split(" ");
        int nx = digits.length;

        Matrix x0 = new Matrix(nx+3, 1);
        for (int i = 0; i < nx; i = i + 1) {
            x0.set(i, 0, Double.parseDouble(digits[i]));
        }

        //GETS ERROR BOUND
        System.out.println("Type error tolerance number in decimal form:");
        String tolerance = keyboard.nextLine();
        double tol = Double.parseDouble(tolerance);

        //MAKES A0 AND A1 MATRICES
        Matrix a0 = new Matrix(n+3, n);
        Matrix a1 = new Matrix(n+3, n);
        for (int i = 0; i < n; i = i + 1) {
            a0.set(i, i, 1);
            a0.set(i+2, i, 1);
            a0.set(i+3, i, 1);
        }
        for (int i = 0; i < n; i = i + 1) {
            a1.set(i, i, 1);
            a1.set(i+1, i, 1);
            a1.set(i+3, i, 1);
        }
        a0 = a0.getMatrix(0, n-1, 0, n-1);
        a1 = a1.getMatrix(0, n-1, 0, n-1);

        //DECODES VIA JACOBI AND GAUSS_SEIDEL
        System.out.println("x estimation using Gauss-Seidel and A0(x)=y0:");
        gauss_seidel(a0, y0, tol);
        System.out.println("x estimation using Gauss-Seidel and A1(x)=y1:");
        gauss_seidel(a1, y1, tol);
        System.out.println("x estimation using Jacobi and A0(x)=y0:");
        jacobi(a0, y0, tol);
        System.out.println("x estimation using Jacobi and A1(x)=y1:");
        jacobi(a1, y1, tol);
    }

//----------------------------------------------------------------------------------
    private static void decode(Matrix y) { //IF USING VECTOR CREATED BY ENCODE()
        Scanner keyboard = new Scanner(System.in);

        System.out.println();
        System.out.println("+------------------+");
        System.out.println("| Decoding Problem |");
        System.out.println("+------------------+");
        System.out.println();

        //RE-CREATES Y0 AND Y1 VECTORS FROM THE VECTOR CREATED BY ENCODE()
        int n = y.getRowDimension();
        Matrix y0 = new Matrix(n,1);
        Matrix y1 = new Matrix(n,1);

        for (int i = 0; i < n; i = i + 1) {
            y0.set(i,0, y.get(i,0));
            y1.set(i,0, y.get(i,1));
        }

        //GETS ERROR BOUND
        System.out.println("Type error tolerance number in decimal form:");
        String tolerance = keyboard.nextLine();
        double tol = Double.parseDouble(tolerance);

        //MAKES A0 AND A1 MATRICES
        Matrix a0 = new Matrix(n+3, n);
        Matrix a1 = new Matrix(n+3, n);
        for (int i = 0; i < n; i = i + 1) {
            a0.set(i, i, 1);
            a0.set(i+2, i, 1);
            a0.set(i+3, i, 1);
        }
        for (int i = 0; i < n; i = i + 1) {
            a1.set(i, i, 1);
            a1.set(i+1, i, 1);
            a1.set(i+3, i, 1);
        }
        a0 = a0.getMatrix(0, n-1, 0, n-1);
        a1 = a1.getMatrix(0, n-1, 0, n-1);

        //DECODES VIA JACOBI AND GAUSS_SEIDEL
        System.out.println("x estimation using Jacobi and A0(x)=y0:");
        System.out.println();        
        jacobi(a0, y0, tol);
        System.out.println("x estimation using Jacobi and A1(x)=y1:");
        System.out.println();        
        jacobi(a1, y1, tol);
        System.out.println("x estimation using Gauss-Seidel and A0(x)=y0:");
        System.out.println();
        gauss_seidel(a0, y0, tol);
        System.out.println("x estimation using Gauss-Seidel and A1(x)=y1:");
        gauss_seidel(a1, y1, tol);
    }

    public static void jacobi(Matrix a, Matrix y, double tol) {

        int n = a.getRowDimension();
        Matrix x = new Matrix(n,1,0);

        //SEPARATE A INTO D AND L+U
        Matrix d = new Matrix(n,n);
        Matrix lu = a.copy();
        for (int i = 0; i < n; i = i + 1) {
            d.set(i,i, a.get(i,i));
            lu.set(i,i, 0);
        }

        int times = 0;

        double iterations = getIterations(d,lu,tol);

        if (iterations == -1) {
            System.out.println("This method doesn't converge within a fixed amount of iterations.%n");
            return;
        }

        Matrix b = new Matrix(n,1);
        while (times <= iterations) {

            b = mult(lu.uminus(), x);
            b.plusEquals(y);

            //x = mult(d, b);
            x = forwardSub(d, b);

            times = times + 1;
        }
        x = mod(x);
        x.print(0,0);
        System.out.printf("Number of iterations taken: %.3f%n", iterations);

    }

    public static void gauss_seidel(Matrix a, Matrix y, double tol) {

        int n = a.getRowDimension();
        Matrix x = new Matrix(n,1,0);

        //SEPARATE A INTO D+L and U
        Matrix dl = a.copy();
        Matrix u = a.copy();
        for (int i = 0; i < n; i = i + 1) {
            for (int j = 0; j < n; j = j + 1) {
                if ((j-i) > 0) {
                    dl.set(i,j, 0);
                } else {
                    u.set(i,j, 0);
                }
            }
        }

        //FIND INVERSE OF LEFT SIDE
        Matrix dlI = invertLower(dl);
        dlI = mod(dlI);
        double iterations = getIterations(dlI, u, tol);

        int times = 0;
        Matrix b = new Matrix(n,1);
        while(times <= iterations) {

            b = mult(u.uminus(),x);
            b.plusEquals(y);
            x = forwardSub(dl, b);

            times = times + 1;
        }

        x = mod(x);
        x.print(0,0);

        System.out.printf("Number of iterations taken: %.3f%n", iterations);
    }

//--------HELPER FUNCTIONS------------------------------------------

    //TO USE AFTER ADDING/SUBTRACTION
    public static Matrix mod(Matrix m) {
        int r = m.getRowDimension();
        int c = m.getColumnDimension();
        for (int i = 0; i < r; i = i + 1) {
            for (int j = 0; j < c; j = j + 1) {
                m.set(i,j, Math.abs(m.get(i,j))%2);
            }
        }
        return m;
    }

    public static double getIterations(Matrix s, Matrix t, double tol) {

        Matrix a = mult(s,t);

        int rows = a.getRowDimension();
        for (int i = 0; i < rows; i = i + 1){
            for (int j = 0; j < rows; j = j + 1) {
                a.set(i,j, Math.abs(a.get(i,j)));
            }
        }
        Matrix u = new Matrix(rows, 1, 1);

        //GET LARGEST EIGENVALUE (power method)
        double iterations = 0;
        Matrix previous = new Matrix(1,1);
        double topvalue = 1;
        double previousvalue = 1;
        double currentvalue = 1;
        double oldapprox = 1;
        double approximation = 1;
        double subtract = 1;
        double absolute = Math.abs(subtract);
        double eigenvalue = 1;
        Matrix uvector = u.copy();
        while(absolute > tol) {
            topvalue = uvector.get(0,0);
            //System.out.println("Part 1");
            previous.set(0,0,topvalue);
            previousvalue = previous.get(0,0);
            //System.out.println("Part 2");
            uvector = mult(a, uvector);
            //System.out.println("Part 3");
            currentvalue = uvector.get(0,0);
            //System.out.println("Part 4");
            approximation = currentvalue / previousvalue;
            //System.out.println("Part 5");
            subtract = (approximation - oldapprox);
            oldapprox = approximation;
            //System.out.println("Part 6");
            absolute = Math.abs(subtract);
            //System.out.println("Part 7");
            iterations++;
        }
        eigenvalue = Math.abs(approximation);

        if (Double.isNaN(eigenvalue) || eigenvalue == 0) {
            System.out.println("Norm = 0, method converges, cannot solve for number iterations");
            System.out.println("Running default 100 iterations...");
            return 100;
        }
        System.out.println("Norm: " + eigenvalue);

        if (eigenvalue >=1) {
            return -1;
        }

        return ((double) Math.log(tol) / (double) Math.log(eigenvalue));
    }

    public static Matrix mult(Matrix a, Matrix b) {
        int r = a.getRowDimension();
        int c = a.getColumnDimension();
        int e = b.getColumnDimension();
        Matrix fin = new Matrix(r,e);
        for (int x = 0; x < e; x = x + 1) {
            for (int i = 0; i < r; i = i + 1) {
                double ans = 0;
                for (int j = 0; j < c; j = j + 1) {
                    ans = ans + (a.get(i,j) * b.get(j,x));
                }
                fin.set(i,x, ans);
            }
        }
        return fin;
    }

    //SOLVES A SOLUTION THAT CAN USE FORWARD SUBSTITUTION
    public static Matrix forwardSub(Matrix m, Matrix b) {
        int n = m.getRowDimension();
        Matrix x = new Matrix(n,1);
        x.set(0,0, b.get(0,0)/m.get(0,0));
        for (int i = 1; i < n; i = i + 1) {
            double add = 0;
            for(int j = i - 1; j >= 0; j = j - 1) {
                add = add + m.get(i,j)*x.get(j,0);
            }
            double value = b.get(i,0) - add;
            value = value / m.get(i,i);
            x.set(i,0, value);
        }
        return x;
    }

    //USE FORWARD SUB WITH IDENTITY MATRIX TO GET INVERSE
    public static Matrix invertLower(Matrix m) {
        int n = m.getRowDimension();
        Matrix x = new Matrix(n,n);
        for(int col = 0; col < n; col = col + 1) {
            x.set(col,col, 1/m.get(col,col));
            for (int i = col+1; i < n; i = i + 1) {
                double add = 0;
                for(int j = i - 1; j >= 0; j = j - 1) {
                    add = add + m.get(i,j)*x.get(j,col);
                }
                double value = 0 - add;
                value = value / m.get(i,i);

                x.set(i,col, value);
            }
        }
        return x;
    }

    //takes the diagonal elements and makes them inverse
    public static Matrix invertDiag(Matrix m) {
        int n = m.getRowDimension();
        Matrix k = m.copy();
        for (int i = 0; i < n; i = i + 1) {
            k.set(i,i, (1/m.get(i,i)));
        }
        return k;
    }

}