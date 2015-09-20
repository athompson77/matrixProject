import java.util.Scanner;
import java.io.File;
import java.util.ArrayList;
import java.lang.Math;
import java.io.PrintWriter;

public class solve_qr_b extends qr_fact {
    
    public static void main(String[] args) throws Exception {
        
        double[][] A = getMatrix(args[0]);
        double[][] b = getMatrix(args[1]);
        String dimensionString = args[2];
        int dimension = Integer.parseInt(dimensionString);
        
        double[][] Q;
        double[][] R; 

        if (args[3].equals("hh")) {
            Q = hhCalculateQR(A,"Q");
            R = hhCalculateQR(A,"R");
        } else {
            Q = grCalculateQR(A,"Q");
            R = grCalculateQR(A,"R");
        }

        if (args[4].equals("hil") && args[1].equals("hh")) {
            double[][] hilQ = hhCalculateQRHIL(dimension, "Q");
            double[][] hilR = hhCalculateQRHIL(dimension, "R");
            double[][] qTranspose = transposeMatrix(hilQ);
            double[][] x = backwardSub(R, qTranspose);
            printMatrix(x);
        } else {
            double[][] qTranspose = transposeMatrix(Q);
            double[][] x = backwardSub(R, qTranspose);
            printMatrix(x);
        }
    }

    private static double[][] getMatrix(String filename) throws Exception {

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