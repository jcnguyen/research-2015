/**
 * Arrays2.java
 *
 * Contains various methods for manipulating arrays.
 *
 * @author Ludo Waltman
 * @author Nees Jan van Eck
 * @version 1.3.1, 11/17/14
 *
 * @author Complex Network Research
 *         http://research.pomona.edu/complexnetworks/
 * @version 11/18/15
 */

import java.util.Arrays;
import java.util.Random;

public class Arrays2 {

    /**
     * Calculate the sum of the values of an array.
     *
     * @param value  the array
     * @return the sum of the array
     **/
    public static double calcSum(double[] value) {
        double sum;
        int i;

        sum = 0;
        for (i = 0; i < value.length; i++)
            sum += value[i];
        return sum;
    }

    /**
     * Calculate the sum of the values between two indices in an array.
     *
     * @param value       the array
     * @param beginIndex  the beginning index to start the sum (included)
     * @param endIndex    the final index to end the sum (not included)
     * @return the sum between two indices of the array
     **/
    public static double calcSum(double[] value, int beginIndex, int endIndex) {
        double sum;
        int i;

        sum = 0;
        for (i = beginIndex; i < endIndex; i++)
            sum += value[i];
        return sum;
    }

    /**
     * Calculate the average of the values in an array.
     *
     * @param value  the array
     * @return the average of the array
     **/
    public static double calcAverage(double[] value) {
        double average;
        int i;

        average = 0;
        for (i = 0; i < value.length; i++)
            average += value[i];
        average /= value.length;
        return average;
    }

    /**
     * Determine the median value of an array.
     *
     * @param value  the array
     * @return the median value of the array
     **/
    public static double calcMedian(double[] value) {
        double median;
        double[] sortedValue;

        sortedValue = (double[])value.clone();
        Arrays.sort(sortedValue);
        if (sortedValue.length % 2 == 1)
            median = sortedValue[(sortedValue.length - 1) / 2];
        else
            median = (sortedValue[sortedValue.length / 2 - 1] + 
                sortedValue[sortedValue.length / 2]) / 2;
        return median;
    }

    /**
     * Determine the minimum value of an array.
     *
     * @param value  the array
     * @return the minimum vallue of the array
     **/
    public static double calcMinimum(double[] value) {
        double minimum;
        int i;

        minimum = value[0];
        for (i = 1; i < value.length; i++)
            minimum = Math.min(minimum, value[i]);
        return minimum;
    }

    /**
     * Determine the minimum value of an array.
     * Accounts for NaN such that NaN < someDouble returns false.
     *
     * @param value  the array
     * @return the minimum vallue of the array
     **/
    public static double calcMinimum2(double[] value) {
        double minimum;
        int i;

        for (i = 0; i < value.length; i++) {
            if (((Double) value[i]).isNaN()) {
                value[i] = Double.POSITIVE_INFINITY;
            }
        }

        minimum = value[0];
        for (i = 1; i < value.length; i++) {
            minimum = Math.min(minimum, value[i]);
        }

        return minimum;
    }


    /**
     * Determines the maximum value of a double array.
     *
     * @param value  the array
     * @return the maximum value of the array
     **/
    public static double calcMaximum(double[] value) {
        double maximum;
        int i;

        maximum = value[0];
        for (i = 1; i < value.length; i++)
            maximum = Math.max(maximum, value[i]);
        return maximum;
    }

    /**
     * Determines the maximum value of an int array.
     *
     * @param value  the array
     * @return the maximum value of the array
     **/
    public static int calcMaximum(int[] value) {
        int i, maximum;

        maximum = value[0];
        for (i = 1; i < value.length; i++)
            maximum = Math.max(maximum, value[i]);
        return maximum;
    }

    /**
     * Generates a new array of permutations with 
     * integer values [0, size of new array).
     *
     * @param nElements  the size of the new array
     * @return an array of permutations
     **/
    public static int[] generateRandomPermutation(int nElements) {
        return generateRandomPermutation(nElements, new Random());
    }

    /**
     * Generates a new array of permutations with 
     * integer values [0, size of new array).
     *
     * @param nElements  the size of the new array
     * @param random     the random generator
     * @return an array of permutations
     **/
    public static int[] generateRandomPermutation(int nElements, Random random) {
        int i, j, k;
        int[] permutation;

        permutation = new int[nElements];
        for (i = 0; i < nElements; i++)
            permutation[i] = i;
        for (i = 0; i < nElements; i++) {
            j = random.nextInt(nElements);
            k = permutation[i];
            permutation[i] = permutation[j];
            permutation[j] = k;
        }
        return permutation;
    }

}
