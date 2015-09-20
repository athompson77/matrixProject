import Jama.Matrix;
import java.io.File;
import java.util.Scanner;
import java.util.ArrayList;
import java.lang.Math;

public class jacobi {
    
    public static void main(String[] args) throws Exception {
        System.out.println();
        System.out.println("+-------------------+");
        System.out.println("| Jacobi Method Sim |");
        System.out.println("+-------------------+");
        System.out.println();

        Scanner keyboard = new Scanner(System.in);
        System.out.println("Type filename of the augmented [A|b] (nx[n+1]) matrix (including full path and extension):");
        String filename = keyboard.nextLine();
        Matrix ay = getMatrix(filename);
        System.out.println("Type error tolerance number in decimal form:");
        String tolerance = keyboard.nextLine();
        double tol = Double.parseDouble(tolerance);

        //SEPARATE AY INTO A AND Y
        int n = ay.getRowDimension();
        Matrix x0 = new Matrix(n,1,0);
        Matrix a = ay.getMatrix(0,n-1, 0, n-1);
        Matrix y = ay.getMatrix(0,n-1, n, n);

        solve(a, y, x0, tol);
    }

    public static void solve(Matrix a, Matrix y, Matrix x, double tol) {

        int n = a.getRowDimension();

        //SEPARATE A INTO D AND L+U
        Matrix d = new Matrix(n,n);
        Matrix lu = a.copy();;
        for (int i = 0; i < n; i = i + 1) {
            d.set(i,i, a.get(i,i));
            lu.set(i,i, 0);
        }

        int times = 0;
        Matrix dI = invertDiag(d);

        double iterations = getIterations(dI,lu,tol);
        if (iterations == -1) {
            System.out.println("This method doesn't converge within a fixed amount of iterations.");
            return;
        }

        Matrix b = new Matrix(n,1);
        while (times <= iterations){

            b = mult(lu.uminus(), x);
            b.plusEquals(y);
            x = mult(dI, b);

            times = times + 1;
        }
        System.out.println("Estimated x:");
        x.print(1, 5);
        System.out.printf("Estimated number of iterations to solve: %.3f", iterations);

    }

//--------HELPER FUNCTIONS------------------------------------------

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

    private static Matrix getMatrix(String filename) throws Exception {

        Scanner fileReader = new Scanner(new File(filename));

        ArrayList<String> augmented = new ArrayList<>();
        while (fileReader.hasNext()) {
            augmented.add(fileReader.nextLine());
        }

        int rows = augmented.size();
        String[] element = augmented.get(0).split(" ");
        int columns = element.length;

        Matrix m = new Matrix(rows, columns);
        for (int i = 0; i < rows; i = i + 1) {
            String[] elements = augmented.get(i).split(" ");
            for (int j = 0; j < columns; j = j + 1) {
                double value = Double.parseDouble(elements[j]);
                m.set(i, j, value);
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

        //GET EIGENVALUE (power method)
        double iterations = 0;
        Matrix previous = new Matrix(1,1);
        double topvalue = 1;
        double previousvalue = 1;
        double currentvalue = 1;
        double oldapprox = 1;
        double approximation = 1;
        double subtract = 1;
        double absolute = Math.abs(subtract);
        double eigen = 1;
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
        double eigenvalue = Math.abs(approximation);
        if (Double.isNaN(eigenvalue) || eigenvalue == 0) {
            System.out.println("Norm = 0, method converges, cannot solve for number iterations");
            System.out.println("Running default 100 iterations");
            return 100;
        }

        if (eigenvalue >=1) {
            return -1;
        }

        return ((double) Math.log(tol) / (double) Math.log(eigenvalue));
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