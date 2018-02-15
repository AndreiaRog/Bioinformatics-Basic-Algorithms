import java.io.BufferedReader;
import java.io.FileReader; // para ler o ficheiro que contem o Blossum50

import java.util.List;
import java.util.ArrayList;
import java.lang.Math;
import java.lang.StringBuilder;

public class Viterbi
{
    //first and second chain of characters
    private static List<Integer> seqInput = new ArrayList<Integer>();

    private static List<List< Integer>> backTrace;
    private double[][] transitions;

    private double[][] emissions;

    public Viterbi() {

        backTrace = new ArrayList<List<Integer>>();
        for(int i = 0; i < 3; i++) {
            List<Integer> aux = new ArrayList<Integer>();

            backTrace.add(aux);
        }

        transitions = new double[3][3];

        //entry 0 is the state start
        transitions[0][0] = 0.6;
        transitions[0][1] = 0.4;
        transitions[0][2] = 0;

        //entry 1 us state exinb
        transitions[1][0] = 0.25;
        transitions[1][1] = 0.5;
        transitions[1][2] = 0.25;


        //entry 2 is t«state intron»
        transitions[2][0] = 0.25;
        transitions[2][1] = 0.25;
        transitions[2][2] = 0.5;

        emissions = new double[3][4];

        emissions[0][0] = 0.4;
        emissions[0][1] = 0.3;
        emissions[0][2] = 0;
        emissions[0][3] = 0.3;

        emissions[1][0] = 0.1;
        emissions[1][1] = 0.1;
        emissions[1][2] = 0.4;
        emissions[1][3] = 0.4;

        emissions[2][0] = 0.4;
        emissions[2][1] = 0.3;
        emissions[2][2] = 0.3;
        emissions[2][3] = 0;
    }


    public static void main (String[] args)
    {
        // receiving a sequence we transform the characters in their index
        //reads the first chain of amino acids
        String s = new String(args[0]);

        //cycle to analyse each character in the string and convert it to an integer to the list seqInput
    	for(int i=0; i< s.length(); i++){
    		if(s.charAt(i) == 'A') { 
                seqInput.add(0);
            }

            else if(s.charAt(i) == 'T'){
                seqInput.add(1);
            }

            else if(s.charAt(i) == 'C'){
                seqInput.add(2);
            }

            else {
                seqInput.add(3);
            }
    	}

    	System.out.println(s.length());
        Viterbi v = new Viterbi();

        //-2 is the arbitrary value we chose to represent the final state of the Viterbi algorithm
        double value = v.viterbi(-2, s.length());

        System.out.println(backTrace.get(0).size());
        System.out.println(backTrace.get(1).size());
        System.out.println(backTrace.get(2).size());
    }

    public double viterbi (int state, int index){

        if(index == 0) {
            return 0.333333333333333333333333; //this is just the best way to represent 1 / 3 more threes could be added :)
        } else if(index == 1) {                //we make a different approach to when the index is 1
            return viterbi(-1, index-1) * emissions[state][seqInput.get(index-1)]; 
        }

        // o estado só é -1 por imposição neste if, que corresponde a ter o index 1 portanto considerar o index 0 na proxima.

        // e quando o state é -2 e o index 1 entra no segundo e nao no terceiro

        //-2 final state
        if(state == -2) {
            double prob0 = viterbi(0, index);
            double prob1 = viterbi(1, index);
            double prob2 = viterbi(2, index);

            return Math.max(prob0, Math.max(prob1, prob2));


        }

        double prob0 = viterbi(0, index-1) * transitions[0][state] * emissions[state][seqInput.get(index-1)];
        double prob1 = viterbi(1, index-1) * transitions[1][state] * emissions[state][seqInput.get(index-1)];
        double prob2 = viterbi(2, index-1) * transitions[2][state] * emissions[state][seqInput.get(index-1)];

	    if(state>-1){
	        if (prob0>prob1){
	        	if (prob0>prob2){
	        		backTrace.get(state).add(0); 
	        	}
	        }
	        else if (prob1>prob0){
	        	if (prob1>prob2){
	        		backTrace.get(state).add(1);
	        	}
	        }
	        else if (prob2>prob0){
	        	if (prob2>prob1){
	        		backTrace.get(state).add(2);
	        	}
	        }
	    }
	    
        return Math.max(prob0, Math.max(prob1, prob2));

        
        // add to end of List of BackTrack, with #number equal to state we are in, the state of the maximum probability resulting from this calculation
    }
}
