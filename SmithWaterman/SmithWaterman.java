import java.io.BufferedReader;
import java.io.FileReader; // para ler o ficheiro que contem o Blossum50

import java.util.List;
import java.util.ArrayList;
import java.lang.Math;
import java.lang.StringBuilder;

public class SmithWaterman
{
    private static List<Integer> s1 = new ArrayList<Integer>();
    private static List<Integer> s2 = new ArrayList<Integer>();

    private static List<List<Integer>> scoreMatrix;

    private static Integer gapPenalty;
    private static List<String> aminoAcids = new ArrayList<String>();

    private List< Integer> maxScore = new ArrayList<Integer>();
    private List < List <Integer> > parents;
    private List < List <Integer> > smithMatrix;

    private int maxGlobal = 0;

    public SmithWaterman() {

        parents = new ArrayList < List < Integer > >(s1.size() + 1);
        smithMatrix = new ArrayList < List < Integer > >(s1.size() + 1);
        for(int i = 0; i < s1.size() + 1; i++) {
            List< Integer> construtorMatrices = new ArrayList<Integer>(s2.size() + 1);

            for(int j = 0; j < s2.size() +1; j++) {
                construtorMatrices.add(-1); // mesma coisa que não estar inicializado
            }
            parents.add(construtorMatrices);
        }
        for(int i = 0; i < s1.size() + 1; i++) {
            List< Integer> construtorMatrices = new ArrayList<Integer>(s2.size() + 1);

            for(int j = 0; j < s2.size() +1; j++) {
                construtorMatrices.add(-1);
            }
            smithMatrix.add(construtorMatrices);
        }
    }

    public static void main (String[] args)
    {
        scoreMatrix = new ArrayList < List <Integer>  >();
        try(BufferedReader br = new BufferedReader(new FileReader(args[3]))) { // input s1, s1, gapPenalty, Blossum50.txt
            String line;

            String first = br.readLine();
            String[] aminoAux = first.split(" ");
            for (String amino : aminoAux) {
                aminoAcids.add(amino);
            }
            while((line = br.readLine()) != null) {
                String[] row = line.split(" "); // remover os epaços e colocar o resto numa array

                List<Integer> list = new ArrayList<Integer>();

                for (String value : row) {
                    list.add(new Integer(value));
                }

                scoreMatrix.add(list);
            }
        } catch(Exception e) {
            System.out.println(e);
        }

        for(String i : args[0].split("")) {
            s1.add(new Integer(aminoAcids.indexOf(i)));
        }

        for(String i : args[1].split("")) {
            s2.add(new Integer(aminoAcids.indexOf(i)));
        }

        gapPenalty = new Integer(args[2]);

        SmithWaterman sm = new SmithWaterman();

        sm.smithWaterman(s1.size(), s2.size());

        List < List <Integer>> positions = sm.maxScorePositions();

        System.out.print("The score of the optimal alignment is: ");
        System.out.println(sm.maxGlobal);

        System.out.println("\nThe possible optimal alignments between s1 and s2 are:\n");
        sm.traceback(positions);
    }

    public int smithWaterman (int i, int j){
        int parent = 0;
        Integer up, left, diagonal;
        if(i == 0 || j == 0) {
            smithMatrix.get(i).set(j, new Integer(0));
            return 0;
        }

        up = smithWaterman(i - 1, j) + gapPenalty;
        left = smithWaterman(i, j -1) + gapPenalty;
        diagonal = smithWaterman(i-1, j-1) + scoreMatrix.get(s1.get(i-1)).get(s2.get(j-1)).intValue();

        int max = Math.max(up, Math.max(0, Math.max(left, diagonal)));

        smithMatrix.get(i).set(j, new Integer(max));

        if(left == max) {
            parent += 1;
        }
        if(up == max) {
            parent += 3;
        }
        if(diagonal == max) {
            parent += 5;
        }

        parents.get(i).set(j, new Integer(parent));

        if(max > maxGlobal) maxGlobal = max;

        return max;
    }

    public List < List < Integer >> maxScorePositions() {
        List < List < Integer >> positions = new ArrayList<List<Integer>> ();
        for(int i = 1; i < s1.size() + 1; i++) {
            for(int j = 1; j < s2.size() + 1; j++) {
                if(smithMatrix.get(i).get(j) == maxGlobal) {
                    List<Integer> aux = new ArrayList<Integer>();

                    aux.add(i); aux.add(j);

                    positions.add(aux);
                }
            }
        }
        return positions;
    }

    public void traceback(List<List<Integer>> positions) {
        System.out.println("Optimal alignments:");
        for(int i = 0; i < positions.size(); i++) {
            dynamicPrintAlignments(positions.get(i).get(0).intValue(), positions.get(i).get(1).intValue(), new StringBuilder(), new StringBuilder());
        }
    }

    public void dynamicPrintAlignments(int i, int j, StringBuilder sb1, StringBuilder sb2) {
        if(smithMatrix.get(i).get(j) != 0) {
            if(parents.get(i).get(j) == 1) {
                dynamicPrintAlignments(i, j-1,
                    sb1.append("-"),
                    sb2.append(aminoAcids.get(s2.get(j-1)))
                );
            } else if(parents.get(i).get(j) == 3) {
                dynamicPrintAlignments(i-1, j,
                    sb1.append(aminoAcids.get(s1.get(i-1))),
                    sb2.append("-")
                );
            } else if(parents.get(i).get(j) == 5) {
                dynamicPrintAlignments(i-1, j-1,
                    sb1.append(aminoAcids.get(s1.get(i-1))),
                    sb2.append(aminoAcids.get(s2.get(j-1)))
                );
            } else if(parents.get(i).get(j) == 4) {
                dynamicPrintAlignments(i, j-1,
                    (new StringBuilder(sb1.toString())).append("-"),
                    (new StringBuilder(sb2.toString())).append(aminoAcids.get(s2.get(j-1)))
                );
                dynamicPrintAlignments(i-1, j,
                    sb1.append(aminoAcids.get(s1.get(i-1))),
                    sb2.append("-")
                );
            } else if(parents.get(i).get(j) == 6) {
                dynamicPrintAlignments(i-1, j-1,
                    (new StringBuilder(sb1.toString())).append(aminoAcids.get(s1.get(i-1))),
                    (new StringBuilder(sb2.toString())).append(aminoAcids.get(s2.get(j-1)))
                );
                dynamicPrintAlignments(i, j-1,
                    sb1.append("-"),
                    sb2.append(aminoAcids.get(s2.get(j-1)))
                );
            } else if(parents.get(i).get(j) == 8) {
                dynamicPrintAlignments(i-1, j-1,
                    (new StringBuilder(sb1.toString())).append(aminoAcids.get(s1.get(i-1))),
                    (new StringBuilder(sb2.toString())).append(aminoAcids.get(s2.get(j-1)))
                );
                dynamicPrintAlignments(i-1, j,
                    sb1.append(aminoAcids.get(s1.get(i-1))),
                    sb2.append("-")
                );
            } else if(parents.get(i).get(j) == 9) {
                dynamicPrintAlignments(i, j-1,
                    (new StringBuilder(sb1.toString())).append("-"),
                    (new StringBuilder(sb2.toString())).append(aminoAcids.get(s2.get(j-1)))
                );
                dynamicPrintAlignments(i-1, j-1,
                    (new StringBuilder(sb1.toString())).append(aminoAcids.get(s1.get(i-1))),
                    (new StringBuilder(sb2.toString())).append(aminoAcids.get(s2.get(j-1)))
                );
                dynamicPrintAlignments(i-1, j,
                    sb1.append(aminoAcids.get(s1.get(i-1))),
                    sb2.append("-")
                );
            }

        } else {
            System.out.println();
            System.out.println(sb1.reverse());
            System.out.println(sb2.reverse());
            System.out.println();
        }
    }
}
