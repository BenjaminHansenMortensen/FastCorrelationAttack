import java.util.*;
import java.util.stream.IntStream;

public class FastCorrelationAttack {
    private int[] keyStream;
    private int[] c_n;
    private int[] x_m;
    // Remember to change boolean function as well in FilterGeneratorF and FilterGeneratorG


    private int[] f;  // Boolean function output for binary inputs X_1, ... , X_n from 0 to 2^m - 1, where m is the number of inputs to F
    private int termsInCharPoly;  // Refactor into calS and calAPosteriorProb

    public FastCorrelationAttack(int[] keyStream, int[] c_n, int[] x_m) {
        this.keyStream = keyStream;
        this.c_n = c_n;
        this.x_m = x_m;

        analysis();
    }

    private void analysis() {
        FilterGeneratorF generatorF = new FilterGeneratorF();
        double p = calPriorProbability();
        ArrayList<int[]> relations = getEquationsOfSpareRelations();
        double[] posteriorProbabilities = calPostPriorProbabilities(relations, p);
        composeAndSolveRelations(posteriorProbabilities, generatorF);
    }

    // Composes the relations for V_n then finds solutions for S_0, ... , S_n and checks if it's the correct initial state of the LFSR
    private void composeAndSolveRelations(double[] posteriorProbabilities, FilterGeneratorF generatorF) {
        int[] indicesOfGreatestValues = getIndicesOfGreatestValues(posteriorProbabilities, c_n.length + 1);
        double p = getMin(indicesOfGreatestValues, posteriorProbabilities);
        double q = 1 - p;
        System.out.println();
        System.out.println("q* of value: " + q);

        FilterGeneratorG g = new FilterGeneratorG();
        int[][] u_n = new int[indicesOfGreatestValues.length][c_n.length];
        int[] x_n = new int[c_n.length + 1];
        int n = 0;
        for (int i = 0; i < keyStream.length; i++) {
            int[] u_i = g.clockRegisters();
            if (contains(indicesOfGreatestValues, i)) {
                u_n[n] = u_i;
                x_n[n] = keyStream[i];
                n++;
            }
        }

        MatrixOverF2 systemOfEquations = new MatrixOverF2(u_n);
        systemOfEquations.printMatrix();
        ArrayList<int[]> vectors = getVectorsWithWeightQ(q);
        for (int[] vector : vectors){
            int[] b = g.xor(x_n, vector);
            ArrayList<int[]> solutions = systemOfEquations.augmentSol(b);
            for (int[] sol : solutions ) {
                if (systemOfEquations.isConsistent()) {
                    if (generatorF.testInitialState(sol)) {
                        System.out.println();
                        System.out.println("Tested vector: " + Arrays.toString(vector) + ", with solution: " + Arrays.toString(sol));
                        System.out.println("Found Initial state of the Filter Generator: " + Arrays.toString(sol));
                        return;
                    }
                }
            }
        }
    }



    // Finds all vectors of weight Q
    private ArrayList<int[]> getVectorsWithWeightQ(double q) {
        ArrayList<Integer> weights = new ArrayList<>();
        int n = 1;
        double weightOfN = n / (c_n.length + 1);  // n + delta, where delta = 1
        while (weightOfN < q) {
            n++;
            weightOfN = (double) n / (c_n.length + 1);
        }
        weights.add(n - 1);
        weights.add(n);


        ArrayList<int[]> vectorsWithWeightQ = new ArrayList<>();
        vectorsWithWeightQ.add(new int[c_n.length + 1]);
        for (int k : weights) {
            for (int[] vector : kbits(c_n.length + 1, k))
                vectorsWithWeightQ.add(vector);
        }

        weights.add(0); // Mannually adds it here since combinations can't does with n = 0
        System.out.println("Therefore trying vectors of weight " + weights);
        return vectorsWithWeightQ;
    }


    // Generates all binary vectors of length n with weight k, source for python code translated into java below
    // Source: https://stackoverflow.com/questions/1851134/generate-all-binary-strings-of-length-n-with-k-bits-set
    private  ArrayList<int[]> kbits(int n, int k) {
        ArrayList<int[]> result = new ArrayList<>();
        for (int[] bits : generate(n, k)) {
            int[] s = new int[n];
            for (int bit : bits)
                s[bit] = 1;
            result.add(s);
        }
        return result;
    }

    // Generates combinations of r elements from the n elements, source for method below
    // Source: https://www.baeldung.com/java-combinations-algorithm
    private ArrayList<int[]> generate(int n, int r) {
        ArrayList<int[]> combinations = new ArrayList<>();
        int[] combination = new int[r];

        // initialize with lowest lexicographic combination
        for (int i = 0; i < r; i++) {
            combination[i] = i;
        }

        while (combination[r - 1] < n) {
            combinations.add(combination.clone());

            // generate next combination in lexicographic order
            int t = r - 1;
            while (t != 0 && combination[t] == n - r + t) {
                t--;
            }
            combination[t]++;
            for (int i = t + 1; i < r; i++) {
                combination[i] = combination[i - 1] + 1;
            }
        }

        return combinations;
    }

    private boolean contains(int[] arr, int key) {
        return Arrays.stream(arr).anyMatch( x -> x == key);
    }

    // Gets minimum value of indices from array
    private double getMin(int[] indices, double[] arr) {
        double min = Arrays.stream(arr).max().getAsDouble();
        for (int index : indices)
            if (arr[index] < min)
                min = arr[index];

        return min;
    }

    //  Gets the indices of the greatest n values in array, source for method below
    //  Source: https://stackoverflow.com/questions/17623468/get-indices-of-n-maximums-in-java-array
    private int[] getIndicesOfGreatestValues(double[] array, int n) {  // From stackoverflow
        //create sort able array with index and value pair
        IndexValuePair[] pairs = new IndexValuePair[array.length];
        for (int i = 0; i < array.length; i++) {
            pairs[i] = new IndexValuePair(i, array[i]);
        }

        //sort
        Arrays.sort(pairs, new Comparator<IndexValuePair>() {
            public int compare(IndexValuePair o1, IndexValuePair o2) {
                return Double.compare(o2.value, o1.value);
            }
        });

        //extract the indices
        int[] result = new int[n];
        for (int i = 0; i < n; i++) {
            result[i] = pairs[i].index;
        }
        return result;
    }

    private class IndexValuePair {
        private int index;
        private double value;

        public IndexValuePair(int index, double value) {
            this.index = index;
            this.value = value;
        }
    }

    private ArrayList<int[]> getEquationsOfSpareRelations() {
        ArrayList<Integer> charPolyExponents = getExponentsOfTermsInTheCharacteristicPolynomial();

        int last = termsInCharPoly - 1;
        ArrayList<int[]> relations = new ArrayList<>();
        while (charPolyExponents.get(last) < keyStream.length) {
            ArrayList<Integer> exponentsCopy = new ArrayList<>(charPolyExponents);
            while (exponentsCopy.get(last) < keyStream.length) {
                int[] relation = new int[keyStream.length + 1]; // + 1 for RHS
                int rightHandSide = 0;
                for (int exponent : exponentsCopy) {
                    relation[exponent] = 1;
                    rightHandSide += keyStream[exponent];
                }
                relation[keyStream.length] = rightHandSide % 2;
                relations.add(relation);

                for (int x = 0; x < exponentsCopy.size(); x++) {
                    int xPower = exponentsCopy.get(x);
                    exponentsCopy.set(x, xPower + 1);
                }
            }
            for (int x = 0; x < charPolyExponents.size(); x++) {
                int xPower = charPolyExponents.get(x);
                charPolyExponents.set(x, xPower * 2);
            }
        }
        return relations;
    }

    private ArrayList<Integer> getExponentsOfTermsInTheCharacteristicPolynomial() {
        ArrayList<Integer> charPolyExponents = new ArrayList<>();
        for (int i = c_n.length - 1; i >= 0 ; i--) {
            if (c_n[i] != 0)
                charPolyExponents.add(c_n.length - 1 - i);
        }
        charPolyExponents.add(c_n.length);
        termsInCharPoly = charPolyExponents.size();
        return charPolyExponents;
    }

    private double[] calPostPriorProbabilities(ArrayList<int[]> relations, double p) {
        int[] k = IntStream.range(0, keyStream.length).toArray();
        int[] mk = new int[keyStream.length];
        int[] hk = new int[keyStream.length];

        for (int[] relation : relations)
            for (int i = 0; i < keyStream.length; i++)
                if (relation[i] != 0) {
                    mk[i] += 1;
                    if (relation[keyStream.length] == 0)  // RHS of relation
                        hk[i] += 1;
                }

        double s = calS(p);
        double[] pk = new double[keyStream.length];
        for (int i = 0; i < pk.length; i++) {
            pk[i] = (p * Math.pow(s, hk[i]) * Math.pow(1 - s, mk[i] - hk[i])) /
                    (p * Math.pow(s, hk[i]) * Math.pow(1 - s, mk[i] - hk[i]) +
                     (1 - p) * Math.pow(1- s, hk[i]) * Math.pow(s, mk[i] - hk[i]));
        }

        System.out.println();
        System.out.println("A posterior probabilities: with S(p,t) = " + s);
        System.out.println("Maximum p: " + Arrays.stream(pk).max().getAsDouble());
        System.out.println("Pk" + Arrays.toString(pk));
        return pk;
    }

    // Calculates the value S for S(p,t)
    private double calS(double p) {
        double s = p;  // Recursion 1
        for (int i = 1; i < termsInCharPoly - 1; i++) {
            s = p * s + (1 - p) * (1 - s);
        }
        return s;
    }

    // Calculates the prior probability and finds the best approximation of F
    private double calPriorProbability() {  // OBS: Does not handle affine approximations correctly
        double[] walshHadamardSpectrum = calWalshHadamardSpecter();

        double max = 0.0;
        int maxIndex = 0;
        boolean affine = false;
        boolean negative;
        for (int i = 0; i < walshHadamardSpectrum.length; i++) {
            if (walshHadamardSpectrum[i] < 0)
                negative = true;
            else
                negative = false;

            walshHadamardSpectrum[i] = ((1 + Math.abs(walshHadamardSpectrum[i])) / 2);
            if (walshHadamardSpectrum[i] > max) {
                max = walshHadamardSpectrum[i];
                maxIndex = i;
                if (negative)
                    affine = true;
                else
                    affine = false;
            }
        }
        findApproximation(maxIndex, affine, max);

        return max;
    }

    // Binary expansion to find affine/linear approximation of F
    private void findApproximation(int maxIndex, boolean affine, double max) {  // OBS: Need to manually update boolean function in G
        String approx = "";
        String[] binaryExpansion = String.format("%" + x_m.length +"s", Integer.toBinaryString(maxIndex)).replace(' ', '0').split("");
        for (int i = 0; i < binaryExpansion.length; i++) {
            if (binaryExpansion[i].equals("1")) {
                if (!approx.equals(""))
                    approx += " + X_" + (i + 1);
                else
                    approx += "X_" + (i + 1);
            }
        }

        if (affine) {
            for (int i = 0; i < keyStream.length; i++) {  // Lucky, affine not handled in boolean function which it should have been...
                keyStream[i] = keyStream[i] ^ 1;
            }
            approx += " + 1";
        }
        System.out.println();
        System.out.println("Best affine/linear approximation of F: " + approx + ", with a prior probability: " + max);
    }

    // Calculates the Walsh Hadamard Spectrum for F
    private double[] calWalshHadamardSpecter() {
        double[] walshHadamardSpectrum = new double[f.length];
        for (int i = 0; i < f.length; i++)  // Initialization
            walshHadamardSpectrum[i] = Math.pow(-1, f[i]);

        int lastSteps = 1;
        for (int steps = 2; steps <= f.length; steps *= 2) {
            for (int i = 0; i < f.length - lastSteps; i += steps) {
                for (int j = 0; j < lastSteps; j++) {
                    double buff = walshHadamardSpectrum[i + j];
                    walshHadamardSpectrum[i + j] = walshHadamardSpectrum[i + j] + walshHadamardSpectrum[i + j + lastSteps];
                    walshHadamardSpectrum[i + j + lastSteps] = buff - walshHadamardSpectrum[i + j + lastSteps];
                }
            }
            lastSteps = steps;
        }

        for (int i = 0; i < f.length; i++)
            walshHadamardSpectrum[i] = walshHadamardSpectrum[i] / f.length;

        return walshHadamardSpectrum;
    }

    private class FilterGeneratorG {
        int numOfInputVariablesToG;  // The number of input variables to the boolean function F
        int degreeOfCharPoly;  // Degree of the characteristic polynomial for the LFSR (number of registers )
        LinkedList<int[]> registers = new LinkedList<>();  // Registers in the LFSR
        ArrayList<Integer> tapsIndices = new ArrayList<>();  // Taps for recurrence for the LFSR

        public FilterGeneratorG() {
            this.numOfInputVariablesToG = x_m.length;
            this.degreeOfCharPoly = c_n.length;
            setTaps();  // Sets the taps of the LFSR
            setRegisters();  // Initialises the values of the LFSR
        }

        // Sets the taps of the LFSR
        private void setTaps() {
            int[] reversed = IntStream.range(0, degreeOfCharPoly).map(i -> c_n[degreeOfCharPoly -i-1]).toArray();
            for (int i = 0; i < reversed.length; i++) {
                if (reversed[i] != 0)
                    tapsIndices.add(i);
            }
        }

        // Sets the registers index equal to their binary representation
        private void setRegisters() {
            for (int i = 0; i < c_n.length; i++) {
                int[] s_n = new int[c_n.length];
                s_n[i] = 1;
                registers.add(s_n);
            }
            return;
        }

        // Clocks the register, XORing according to taps to find each registers composition at time t
        private int[] clockRegisters() {  // [ S_0 + S_1, ... , S_39 + S_40 ] at t = 41
            ArrayList<int[]> x_n = new ArrayList<>();
            for (Integer x_i : x_m) {
                x_n.add(registers.get(x_i));
            }
            int[] s_n = booleanFunction(x_n);

            int[] nextS_n = new int[c_n.length];
            for (int tap : tapsIndices)
                nextS_n = xor(nextS_n, registers.get(tap));

            registers.add(nextS_n);
            registers.removeFirst();
            return s_n;
        }


        // ANF of the boolean functions over F_2
        private int[] booleanFunction(ArrayList<int[]> x) {  // Hard coded values
            return xor(x.get(2), xor(x.get(3), x.get(4)));  // All X_n are n - 1 since index starts at 0 not 1
        }

        // XORs together two vectors
        private int[] xor(int[] a1, int[] a2) {
            if (a1.length != a2.length)
                throw new ArrayIndexOutOfBoundsException();

            int[] product = new int[a1.length];

            for (int i = 0; i < a1.length; i++) {
                if (a1[i] != -1 && a2[i] != -1)
                    product[i] = a1[i] ^ a2[i];
                else
                    product[i] = -1;
            }

            return product;
        }

    }

    private class FilterGeneratorF {
        int numOfInputVariablesToF;  // The number of input variables to the boolean function F
        int degreeOfCharPoly;  // Degree of the characteristic polynomial for the LFSR (number of registers )
        LinkedList<Integer> registers = new LinkedList<>();  // Registers in the LFSR
        ArrayList<Integer> tapsIndices = new ArrayList<>();  // Taps for recurrence for the LFSR
        ArrayList<Integer> sequence = new ArrayList<>();  // Output sequence for the filter generator (key stream)


        public FilterGeneratorF() {
            this.numOfInputVariablesToF = x_m.length;
            this.degreeOfCharPoly = c_n.length;
            findTruthTableValues();
            setTaps();
        }

        // Tests if an initial state for the filter generator produces the key stream
        private boolean testInitialState(int[] initialValues) {
            registers.clear();
            sequence.clear();

            setRegisters(initialValues);
            for (int i = 0; i < keyStream.length; i++) {
                clockRegisters();
            }

            return equals(keyStream, sequence.stream().mapToInt(i -> i).toArray());
        }

        // Sets the taps of the LFSR
        private void setTaps() {
            int[] reversed = IntStream.range(0, degreeOfCharPoly).map(i -> c_n[degreeOfCharPoly -i-1]).toArray();
            for (int i = 0; i < reversed.length; i++) {
                if (reversed[i] != 0)
                    tapsIndices.add(i);
            }
        }

        // Calculates the boolean function values for F from 0 to 2^n - 1, where n is the numbers of input variables to F
        private void findTruthTableValues() {
            f = new int[(int) Math.pow(2, x_m.length)];
            for (int i = 0; i < (int) Math.pow(2, x_m.length); i++) {
                f[i]  = booleanFunction(Arrays.stream(String.format("%" + x_m.length +"s", Integer.toBinaryString(i)).replace(' ', '0').split("")).mapToInt(x -> Integer.parseInt(x)).toArray());
            }

            return;
        }

        // Sets the initial values of the registers in the LFSR
        private void setRegisters(int[] initialValues) {  // Input in form: {S_n, ... , S_1, S_0}
            for (int i = initialValues.length - 1; i >= 0; i--) {
                registers.add(initialValues[i]);
            }
            return;
        }

        // Clocks the register, adding X_i to the sequence of outputs
        private void clockRegisters() {
            int[] x_n = new int[numOfInputVariablesToF];
            for (int i = 0; i < numOfInputVariablesToF; i++) {
                x_n[i] = registers.get(x_m[i]);
            }
            sequence.add(booleanFunction(x_n));

            int nextS_n = 0;
            for (int tap : tapsIndices)
                nextS_n += registers.get(tap);

            registers.add(nextS_n % 2);
            registers.removeFirst();
            return;
        }

        // ANF of the boolean functions over F_2
        private int booleanFunction(int[] x) {  // Hard coded values
            return (x[0] * x[1] + x[2] + x[3] + x[4]) % 2; // All X_n are n - 1 since index starts at 0 not 1
        }

        // Checks if all values with same index are equal
        private boolean equals(int[] group, int[] equals) {
            if (group.length != equals.length)
                return false;

            for (int bit = 0; bit < group.length; bit++)
                if (group[bit] != equals[bit])
                    return false;

            return true;
        }
    }

    private class MatrixOverF2 {  // Rows and Columns start at 0, not 1
        private int m;  // Number of rows
        private int n;  // Number of columns
        private int rank = 0;  // The rank of the matrix
        private int[][] matrix;  // The matrix
        private int[][] a;  // Copy of the matrix where row reduction os performed
        private int [] b;  // A * x = b, b is augmented and included when row reductions is performed

        public MatrixOverF2(int[][] matrix) {
            this.matrix = matrix;
            this.m = matrix.length;
            this.n = matrix[0].length;
        }

        // Takes A and augments it with b. Then row finds all solutions
        public ArrayList<int[]> augmentSol(int[] b) {
            this.b = b;
            this.a = this.matrix.clone();
            toRowEchelonForm();
            return findAllSolutions();
        }

        // Reduces A to row echelon form
        private void toRowEchelonForm(){
            int pivot = 0;

            for (int j = 0; j < n; j++) {
                for (int i = pivot; i < m; i++) {
                    if (a[i][j] == 0)
                        continue;
                    else {
                        swapRows(pivot, i);
                        for (int u = pivot + 1; u < m; u++) {
                            if (a[u][j] == 1)
                                addRows(u, pivot);
                        }

                        pivot++;
                        break;
                    }
                }
            }

            this.rank = pivot;
            return;
        }

        // Finds all solutions for the system when in echelon form
        private ArrayList<int[]> findAllSolutions() {
            ArrayList<int[]> solutions = new ArrayList<>();
            solutions.add(new int[n]);
            int i = rank - 1;

            for (int j = n - 1; j >= 0; j--)  {
                ArrayList<int[]> buff = new ArrayList<>();
                for (int[] sol : solutions) {
                    int sum = 0;
                    for (int u = j + 1; u < n; u++) {
                        sum += a[i][u] * sol[u];
                    }
                    if (isPivot(i,j)) {
                        sol[j] = (b[i] + sum) % 2;
                        buff.add(sol);
                    }
                    else {
                        int[] copy = sol.clone();
                        sol[j] = 0;
                        copy[j] = 1;
                        buff.add(sol);
                        buff.add(copy);
                    }
                }
                if (isPivot(i,j))
                    i--;
                solutions = buff;
            }
            // Reverse state to match LFSR (S_n-1, ... , S_1, S_0)
            for (int x = 0; x < solutions.size(); x++) {
                int[] sol = solutions.get(x);
                solutions.set(x, IntStream.range(0, sol.length).map(j -> sol[sol.length-j-1]).toArray());
            }

            return solutions;
        }

        // Checks if an element at location i,j in A is a pivot element
        private boolean isPivot(int i, int j) {
            if (a[i][j] != 1)
                return false;

            for (int x = 0; x < j; x++) {
                if (a[i][x] != 0)
                    return false;
            }

            return true;
        }


        // Checks is system is consistent
        private boolean isConsistent() {
            for (int i = rank; i < b.length; i++) {
                if (b[i] != 0)
                    return false;
            }
            return true;
        }

        // Adds two rows together in a
        private void addRows(int row, int addedWith) {
            for (int i = 0; i < n; i++) {
                a[row][i] = (a[row][i] + a[addedWith][i]) % 2;
            }
            b[row] = (b[row] + b[addedWith]) % 2;
        }

        // Swaps two rows in a
        private void swapRows(int row, int swappedWith) {
            int[] buffer = a[row];
            a[row] = a[swappedWith];
            a[swappedWith] = buffer;

            int bufferB = b[row];
            b[row] = b[swappedWith];
            b[swappedWith] = bufferB;
        }

        // Swaps two columns in a
        private void swapColumns(int column, int swappedWith) {
            for (int i = 0; i < m; i++) {
                int buffer = a[i][column];
                a[i][column] = a[i][swappedWith];
                a[i][swappedWith] = buffer;
            }
        }

        // Prints the matrix, not A
        public void printMatrix() {
            for (int[] row : this.matrix) {
                System.out.println(Arrays.toString(row));
            }
            System.out.println("\n");
            return;
        }
    }
}

