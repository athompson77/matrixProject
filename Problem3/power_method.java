import Jama.Matrix;
import java.io.File;
import java.util.Scanner;
import java.util.ArrayList;
import java.lang.Math;

public class power_method {
	public static void main(String[] args) throws Exception{
		System.out.println();
        System.out.println("+-----------------------+");
        System.out.println("|  Power Iteration Sim  |");
        System.out.println("+-----------------------+");
        System.out.println();

        Scanner input = new Scanner(System.in);
        
		System.out.println("Type filename of the SQUARE matrix (including full path and extension):");
		String filename = input.nextLine();
        Matrix powermatrix = getMatrix(filename);
		System.out.println("Verified. Note that, for simplicity, the simulator will set the vector u-0 to a vector of ones.");
		System.out.println("...");
		System.out.println("");
		
		int rows = powermatrix.getRowDimension();
		Matrix u = new Matrix(rows, 1, 1);
		Matrix zero = new Matrix(rows, 1, 0);
		
        System.out.println("Set an error tolerance in decimal form.");
		String tolerance = input.nextLine();
		double errortol = Double.parseDouble(tolerance);
		
		System.out.println("Verified. Note that the program will only iterate 75 times.");
		System.out.println("...");
		System.out.println("");
		
		int iterations = 0;
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
		
		while(absolute > errortol) {
			topvalue = uvector.get(0,0);
			previous.set(0,0,topvalue);
			previousvalue = previous.get(0,0);
			uvector = mult(powermatrix, uvector);
			currentvalue = uvector.get(0,0);
			approximation = currentvalue / previousvalue;
			subtract = (approximation - oldapprox);
			oldapprox = approximation;
			absolute = Math.abs(subtract);
			iterations++;
		}
		eigen = approximation;
		
		if(iterations < 75) {
			System.out.printf("The approximated eigenvalue is %.3f ." , eigen);	
			System.out.println(" You went through " + iterations + " iterations.");
			
		}
		else {
			System.out.println("The program could not approximate an eigenvalue within 75 iterations.");
		}
	}

//-------HELPER FUNCTIONS-----------------------------------

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

        //GET EIGENVALUE (power method)


        double eigenvalue = (1.0/3);

        return (Math.log10(tol) / Math.log10(eigenvalue));
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

}