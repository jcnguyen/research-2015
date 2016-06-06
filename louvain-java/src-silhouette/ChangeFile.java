import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.Console;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.Random;
import java.util.HashSet;
import java.util.Set;


public class ChangeFile {


    public static void main(String[] args) throws IOException {
        makeLeftLessThanRight(args[0]);
    }


    private static void makeLeftLessThanRight(String fileName) throws IOException {

        String[] splittedLine;
        int start;
        int destination;
        double weight = 1;

        BufferedReader bufferedReader = new BufferedReader(new FileReader(fileName));
        while (bufferedReader.readLine() != null) {
            splittedLine = bufferedReader.readLine().split(" ");
            start = Integer.parseInt(splittedLine[0]);
            destination = Integer.parseInt(splittedLine[1]);

            if (start > destination) {
                int temp = start;
                start = destination;
                destination = temp;
            }

            if (splittedLine.length > 2) {
                weight = Double.parseDouble(splittedLine[2]);
                System.out.println(start + " " + destination + " " + weight);
            } else {
                System.out.println(start + " " + destination);
            }
        }
        bufferedReader.close();

    }

    public static void stripDuplicatesFromFile(String fileName) throws IOException {

        // get the number of lines in the file
        BufferedReader reader = new BufferedReader(new FileReader(fileName));
        int nLines = 0;
        while (reader.readLine() != null)
            nLines++;
        reader.close();

        reader = new BufferedReader(new FileReader(fileName));
        Set<String> lines = new HashSet<String>(nLines); // maybe should be bigger
        String line;
        while ((line = reader.readLine()) != null) {
            lines.add(line);
        }
        reader.close();

        for (String unique : lines) {
            System.out.println(unique);
        }

    }

}
