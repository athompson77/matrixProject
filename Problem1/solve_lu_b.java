import java.util.Scanner;
import java.io.File;
import java.util.ArrayList;
import java.lang.Math;
import java.io.PrintWriter;

public class solve_lu_b extends lu_fact {
    
    public static void main(String[] args) throws Exception {
        
        
        double[][] A = getMatrix(args[0]);
        double[][] b = getMatrix(args[1]);
        String dimensionString = args[0];
        int dimension = Integer.parseInt(dimensionString);
        double[][] L = calculateLU(A, "L");
        double[][] U = calculateLU(A, "U");

        if (args[2].equals("hil")) {
            double[][] hilL = calculateLUHIL(dimension, "L");
            double[][] hilU = calculateLUHIL(dimension, "U");
            double[][] y = forwardSub(hilL, b);
            double[][] x = backwardSub(hilU, y);
            printMatrix(x);
        } else {
            double[][] y = forwardSub(L, b);
            double[][] x = backwardSub(U, y);
            printMatrix(x);
        }
        
    }

    public static double[][] getMatrix(String filename) throws Exception {

        Scanner fileReader = new Scanner(new File(filename));

        ArrayList<String> augmented = new ArrayList<>();
        while (fileReader.hasNext()) {
            augmented.add(fileReader.nextLine());
        }

        int rows = augmented.size();
        String[] element = augmented.get(0).split(" ");
        int columns = element.length;

        double[][] m = new double[rows][columns];
        for (int i = 0; i < rows; i = i + 1) {
            String[] elements = augmented.get(i).split(" ");
            for (int j = 0; j < columns; j = j + 1) {
                double value = Double.parseDouble(elements[j]);
                m[i][j] = value;
            }
        }
        return m;
    }
    public static double[][] forwardSub(double[][] m, double[][] b) {
        int n = m.length;
        double[][] x = new double[n][1];
        x[0][0] = b[0][0]/m[0][0];
        for (int i = 1; i < n; i = i + 1) {
            double add = 0;
            for(int j = i - 1; j >= 0; j = j - 1) {
                add = add + m[i][j]*x[i][j];
            }
            double value = b[i][0] - add;
            value = value / m[i][i];
            x[i][0] =  value;
        }
        return x;
    }
    public static double[][] backwardSub(double[][] m, double[][] b) {
        int n = m.length;
        double[][] x = new double[n][1];
        x[0][0] = b[0][0]/m[0][0];
        for (int i = n; i > -1; i--) {
            double add = 0;
            for(int j = i - 1; j >= 0; j = j - 1) {
                add = add + m[j][i]*x[i][j];
            }
            double value = b[i][0] - add;
            value = value / m[i][i];
            x[i][0] =  value;
        }
        return x;
    }
}