/**
 * Name: Seonjun Mun
 */

package com.company;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

// Class that implements Smith-Waterman local alignment, scoring, traceback, and empirical p-value calculation.
// We mainly use 10 proteins with their sequences and the BLOSUM62 scoring matrix from NCBI for this project
// and compute those values with the Smith-Waterman algorithm, and placed the result in the output.txt file.
public class SmithWaterman {

    private int[][] score_matrix;
    private static int[][] blossom_matrix;
    private static Map<String,Integer> blossom_map;
    private String str1;
    private String str2;
    private String id1;
    private String id2;
    private int m;
    private int n;
    private static Map<String,String> sequences_map;
    private static List<String> identifiers_list;

    // Constructor for SmithWaterman with two sequences and their identifiers
    public SmithWaterman(String str1, String str2, String id1, String id2) {
        this.str1 = str1.toUpperCase();
        this.str2 = str2.toUpperCase();
        this.id1 = id1;
        this.id2 = id2;
        m = str1.length();
        n = str2.length();
        score_matrix = new int[m+1][n+1];
        fillScoreMatrix();
    }

    // Main program wrapper to call the routine to generate cases
    public static void main(String[] args) {
        // Fill the blossom map and the blossom matrix, and fill the sequences map
        fillBlossomMap();
        fillBlossomMatrix();
        fillSequences();

        // First test case
        String w1 = "deadly", w2 = "ddgearlyk";
        SmithWaterman smithWaterman1 = new SmithWaterman(w1, w2, "","");

        int[] maxIJ = smithWaterman1.getMaxPosition();
        String[] opt_alignments1 = smithWaterman1.getOptAlignment(maxIJ[0], maxIJ[1], "", "");

        System.out.println(opt_alignments1[0]);
        System.out.println(opt_alignments1[1]);

        System.out.println("Max score: " + smithWaterman1.getMaxScore());
        int orig_score = smithWaterman1.getMaxScore();
        smithWaterman1.printScoreMatrix();
        System.out.println("p-value: " + smithWaterman1.getEmpiricalPValue(orig_score, 999));
        System.out.println("\n");


        // Generate each pair of two proteins out of 10 proteins that we have
        for (int i = 0; i < identifiers_list.size(); i++) {
            for (int j = i+1; j < identifiers_list.size(); j++) {
                String input1 = sequences_map.get(identifiers_list.get(i));
                String input2 = sequences_map.get(identifiers_list.get(j));

                SmithWaterman sw = new SmithWaterman(input1, input2, identifiers_list.get(i), identifiers_list.get(j));
                int[] maxPos = sw.getMaxPosition();
                String[] opt_alignments = sw.getOptAlignment(maxPos[0], maxPos[1], "", "");

                int counter1 = 0, counter2 = 0;
                int idx1 = 0, idx2 = 0;

                String s1 = opt_alignments[0];
                String s2 = opt_alignments[1];

                // To make printed outputs 60 characters per line with the correct indices
                String[] list1 = new String[s1.length() / 61 + 1];
                String[] list2 = new String[s2.length() / 61 + 1];
                int[] idxList1 = new int[list1.length];
                int[] idxList2 = new int[list2.length];

                for (int i1 = 0; i1 < list1.length; i1++) {
                    idxList1[i1] = idx1+1;
                    int temp = 0;
                    StringBuilder sb = new StringBuilder();
                    while (temp <= 60 && counter1 < s1.length()) {
                        if (s1.charAt(counter1) != '-') idx1++;
                        sb.append(s1.charAt(counter1++));
                        temp++;
                    }
                    list1[i1] = sb.toString();
                }
                for (int i2 = 0; i2 < list2.length; i2++) {
                    idxList2[i2] = idx2+1;
                    int temp = 0;
                    StringBuilder sb = new StringBuilder();
                    while (temp <= 60 && counter2 < s2.length()) {
                        if (s2.charAt(counter2) != '-') idx2++;
                        sb.append(s2.charAt(counter2++));
                        temp++;
                    }
                    list2[i2] = sb.toString();
                }

                for (int k = 0; k < list1.length; k++) {
                    System.out.print(sw.id1 + ":  ");
                    System.out.println(idxList1[k] + "  " + list1[k]);
                    System.out.println();
                    System.out.print(sw.id2 + ":  ");
                    System.out.println(idxList2[k] + "  " + list2[k]);
                    System.out.println();
                }

                System.out.println("Optimal score: " + sw.getMaxScore());

                int original_score = sw.getMaxScore();
                sw.printScoreMatrix();
                System.out.println("p-value: " + sw.getEmpiricalPValue(original_score, 999));
                System.out.println("\n");
            }
        }

        // Calculate empirical p-values for the significance of the alignments of
        // P15172 versus Q10574 and for P15172 versus O95363.
        String P15172 = sequences_map.get("P15172");
        String Q10574 = sequences_map.get("Q10574");
        String O95363 = sequences_map.get("O95363");
        SmithWaterman sw1 = new SmithWaterman(P15172, Q10574, "P15172", "Q10574");
        SmithWaterman sw2 = new SmithWaterman(P15172, O95363, "P15172", "O95363");

        int score1 = sw1.getMaxScore();
        int score2 = sw2.getMaxScore();
        System.out.println("p-value for P15172 versus Q10574 with 999 random trials: "
                + sw1.getEmpiricalPValue(score1, 999));
        System.out.println("p-value for P15172 versus O95363 with 999 random trials: "
                + sw2.getEmpiricalPValue(score2, 999));
    }

    // Construct the sequence map with identifiers as keys and sequences as values
    private static void fillSequences() {
        sequences_map = new HashMap<>();
        sequences_map.put("P10085","MELLSPPLRDIDLTGPDGSLCSFETADDFYDDPCFDSPDLRFFEDLDPRLVHMGALLKPE" +
                "EHAHFPTAVHPGPGAREDEHVRAPSGHHQAGRCLLWACKACKRKTTNADRRKAATMRERR" +
                "RLSKVNEAFETLKRCTSSNPNQRLPKVEILRNAIRYIEGLQALLRDQDAAPPGAAAFYAP" +
                "GPLPPGRGSEHYSGDSDASSPRSNCSDGMMDYSGPPSGPRRQNGYDTAYYSEAARESRPG" +
                "KSAAVSSLDCLSSIVERISTDSPAAPALLLADAPPESPPGPPEGASLSDTEQGTQTPSPD" +
                "AAPQCPAGSNPNAIYQVL");
        sequences_map.put("P13904","MELLPPPLRDMEVTEGSLCAFPTPDDFYDDPCFNTSDMSFFEDLDPRLVHVTLLKPEEPH" +
                "HNEDEHVRAPSGHHQAGRCLLWACKACKRKTTNADRRKAATMRERRRLSKVNEAFETLKR" +
                "YTSTNPNQRLPKVEILRNAIRYIESLQALLHDQDEAFYPVLEHYSGDSDASSPRSNCSDG" +
                "MMDYNSPPCGSRRRNSYDSSFYSDSPNDSRLGKSSVISSLDCLSSIVERISTQSPSCPVP" +
                "TAVDSGSEGSPCSPLQGETLSERVITIPSPSNTCTQLSQDPSSTIYHVL");
        sequences_map.put("P15172","MELLSPPLRDVDLTAPDGSLCSFATTDDFYDDPCFDSPDLRFFEDLDPRLMHVGALLKPE" +
                "EHSHFPAAVHPAPGAREDEHVRAPSGHHQAGRCLLWACKACKRKTTNADRRKAATMRERR" +
                "RLSKVNEAFETLKRCTSSNPNQRLPKVEILRNAIRYIEGLQALLRDQDAAPPGAAAAFYA" +
                "PGPLPPGRGGEHYSGDSDASSPRSNCSDGMMDYSGPPSGARRRNCYEGAYYNEAPSEPRP" +
                "GKSAAVSSLDCLSSIVERISTESPAAPALLLADVPSESPPRRQEAAAPSEGESSGDPTQS" +
                "PDAAPQCPAGANPNPIYQVL");
        sequences_map.put("P16075","MDLLGPMEMTEGSLCSFTAADDFYDDPCFNTSDMHFFEDLDPRLVHVGGLLKPEEHPHTR" +
                "APPREPTEEEHVRAPSGHHQAGRCLLWACKACKRKTTNADRRKAATMRERRRLSKVNEAF" +
                "ETLKRCTSTNPNQRLPKVEILRNAIRYIESLQALLREQEDAYYPVLEHYSGESDASSPRS" +
                "NCSDGMMEYSGPPCSSRRRNSYDSSYYTESPNDPKHGKSSVVSSLDCLSSIVERISTDNS" +
                "TCPILPPAEAVAEGSPCSPQEGGNLSDSGAQIPSPTNCTPLPQESSSSSSSNPIYQVL");
        sequences_map.put("P17542","MTERPPSEAARSDPQLEGRDAAEASMAPPHLVLLNGVAKETSRAAAAEPPVIELGARGGP" +
                "GGGPAGGGGAARDLKGRDAATAEARHRVPTTELCRPPGPAPAPAPASVTAELPGDGRMVQ" +
                "LSPPALAAPAAPGRALLYSLSQPLASLGSGFFGEPDAFPMFTTNNRVKRRPSPYEMEITD" +
                "GPHTKVVRRIFTNSRERWRQQNVNGAFAELRKLIPTHPPDKKLSKNEILRLAMKYINFLA" +
                "KLLNDQEEEGTQRAKTGKDPVVGAGGGGGGGGGGAPPDDLLQDVLSPNSSCGSSLDGAAS" +
                "PDSYTEEPAPKHTARSLHPAMLPAADGAGPR");
        sequences_map.put("P22816","MTKYNSGSSEMPAAQTIKQEYHNGYGQPTHPGYGFSAYSQQNPIAHPGQNPHQTLQNFFS" +
                "RFNAVGDASAGNGGAASISANGSGSSCNYSHANHHPAELDKPLGMNMTPSPIYTTDYDDE" +
                "NSSLSSEEHVLAPLVCSSAQSSRPCLTWACKACKKKSVTVDRRKAATMRERRRLRKVNEA" +
                "FEILKRRTSSNPNQRLPKVEILRNAIEYIESLEDLLQESSTTRDGDNLAPSLSGKSCQSD" +
                "YLSSYAGAYLEDKLSFYNKHMEKYGQFTDFDGNANGSSLDCLNLIVQSINKSTTSPIQNK" +
                "ATPSASDTQSPPSSGATAPTSLHVNFKRKCST");
        sequences_map.put("Q8IU24","MEFVELSSCRFDATPTFCDRPAAPNATVLPGEHFPVPNGSYEDQGDGHVLAPGPSFHGPG" +
                "RCLLWACKACKKKTVPIDRRKAATMRERRRLVKVNEAFDILKKKSCANPNQRLPKVEILR" +
                "NAISYIEQLHKLLRDSKENSSGEVSDTSAPSPGSCSDGMAAHSPHSFCTDTSGNSSWEQG" +
                "DGQPGNGYENQSCGNTVSSLDCLSLIVQSISTIEGEENNNASNTPR");
        sequences_map.put("Q10574","MSWEQYQMYVPQCHPSFMYQGSIQSTMTTPLQSPNFSLDSPNYPDSLSNGGGKDDKKKCR" +
                "RYKTPSPQLLRMRRSAANERERRRMNTLNVAYDELREVLPEIDSGKKLSKFETLQMAQKY" +
                "IECLSQILKQDSKNENLKSKSG");
        sequences_map.put("Q90477","MELSDIPFPIPSADDFYDDPCFNTNDMHFFEDLDPRLVHVSLLKPDEHHHIEDEHVRAPS" +
                "GHHQAGRCLLWACKACKRKTTNADRRKAATMRERRRLSKVNDAFETLKRCTSTNPNQRLP" +
                "KVEILRNAISYIESLQALLRSQEDNYYPVLEHYSGDSDASSPRSNCSDGMMDFMGPTCQT" +
                "RRRNSYDSSYFNDTPNADARNNKNSVVSSLDCLSSIVERISTETPACPVLSVPEGHEESP" +
                "CSPHEGSVLSDTGTTAPSPTSCPQQQAQETIYQVL");
        sequences_map.put("O95363","MVGSALRRGAHAYVYLVSKASHISRGHQHQAWGSRPPAAECATQRAPGSVVELLGKSYPQ" +
                "DDHSNLTRKVLTRVGRNLHNQQHHPLWLIKERVKEHFYKQYVGRFGTPLFSVYDNLSPVV" +
                "TTWQNFDSLLIPADHPSRKKGDNYYLNRTHMLRAHTSAHQWDLLHAGLDAFLVVGDVYRR" +
                "DQIDSQHYPIFHQLEAVRLFSKHELFAGIKDGESLQLFEQSSRSAHKQETHTMEAVKLVE" +
                "FDLKQTLTRLMAHLFGDELEIRWVDCYFPFTHPSFEMEINFHGEWLEVLGCGVMEQQLVN" +
                "SAGAQDRIGWAFGLGLERLAMILYDIPDIRLFWCEDERFLKQFCVSNINQKVKFQPLSKY" +
                "PAVINDISFWLPSENYAENDFYDLVRTIGGDLVEKVDLIDKFVHPKTHKTSHCYRITYRH" +
                "MERTLSQREVRHIHQALQEAAVQLLGVEGRF");

        identifiers_list = new ArrayList<>();
        for (String key: sequences_map.keySet()) {
            identifiers_list.add(key);
        }
    }

    // Construct the blossom map with letters as keys and indices of the letters as values
    private static void fillBlossomMap() {
        blossom_matrix = new int[20][20];
        blossom_map = new HashMap<>();
        blossom_map.put("A", 0);
        blossom_map.put("R", 1);
        blossom_map.put("N", 2);
        blossom_map.put("D", 3);
        blossom_map.put("C", 4);
        blossom_map.put("Q", 5);
        blossom_map.put("E", 6);
        blossom_map.put("G", 7);
        blossom_map.put("H", 8);
        blossom_map.put("I", 9);
        blossom_map.put("L", 10);
        blossom_map.put("K", 11);
        blossom_map.put("M", 12);
        blossom_map.put("F", 13);
        blossom_map.put("P", 14);
        blossom_map.put("S", 15);
        blossom_map.put("T", 16);
        blossom_map.put("W", 17);
        blossom_map.put("Y", 18);
        blossom_map.put("V", 19);
    }

    // Fill the blossom matrix that we will use to calculate optimal scores
    private static void fillBlossomMatrix() {
        BufferedReader reader;
        List<String> lines = new ArrayList<>();
        try {
            reader = new BufferedReader(new FileReader("BLOSUM62.txt"));
            String line = reader.readLine();
            while (line != null) {
                lines.add(line);
                line = reader.readLine();
            }
            reader.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
        int r = 0, c = 0;
        for (String line: lines) {
            c = 0;
            String[] arr = line.split(" ");
            for (String s: arr) {
                if (s.length() == 0) continue;
                int num = Integer.valueOf(s);
                if (r >= 20 || c >= 20) continue;
                blossom_matrix[r][c++] = num;
            }
            r++;
        }
    }


    // Computer the sigma value for str1(i) and str2(j)
    private int sigma(int i, int j) {
        return blossom_matrix[blossom_map.get(String.valueOf(str1.charAt(i-1)))]
                [blossom_map.get(String.valueOf(str2.charAt(j-1)))];
    }

    // Fill the score matrix that we will use to compute the scores for each position
    private void fillScoreMatrix() {
        for (int i = 1; i <= m; i++) {
            score_matrix[i][0] = 0;
        }
        for (int j = 1; j <= n; j++) {
            score_matrix[0][j] = 0;
        }
        for (int i = 1; i <= m; i++) {
            for (int j = 1; j <= n; j++) {
                int diagonal = score_matrix[i-1][j-1] + sigma(i, j);
                int left = score_matrix[i-1][j] -4;
                int up = score_matrix[i][j-1] -4;
                score_matrix[i][j] = Math.max(diagonal, Math.max(left, Math.max(up, 0)));
            }
        }
    }

    // Get the max score of the alignment in the matrix
    private int getMaxScore() {
        int max_score = 0;
        for (int i = 1; i <= m; i++) {
            for (int j = 1; j <= n; j++) {
                if (score_matrix[i][j] > max_score) {
                    max_score = score_matrix[i][j];
                }
            }
        }
        return max_score;
    }

    // Get the position of the max score of the alignment in the matrix
    private int[] getMaxPosition() {
        int max_score = 0;
        int maxI = 0, maxJ = 0;
        for (int i = 1; i <= m; i++) {
            for (int j = 1; j <= n; j++) {
                if (score_matrix[i][j] > max_score) {
                    max_score = score_matrix[i][j];
                    maxI = i;
                    maxJ = j;
                }
            }
        }
        return new int[]{maxI, maxJ};
    }

    // Get the optimal alignments of two sequences by appending to suffix1 and suffix2
    // using the traceback technic from the max score position to the top-left starting point
    private String[] getOptAlignment(int i, int j, String suffix1, String suffix2) {
        if (score_matrix[i][j] == 0) {  // need to fix
            String[] res = new String[2];
            res[0] = suffix1;
            res[1] = suffix2;
            return res;
        }
        if (score_matrix[i][j] + 4 == score_matrix[i-1][j]) {
            return getOptAlignment(i-1, j, str1.charAt(i-1) + suffix1, "-" + suffix2);
        }
        else if (score_matrix[i][j] + 4 == score_matrix[i][j-1]) {
            return getOptAlignment(i, j-1, "-" + suffix1, str2.charAt(j-1) + suffix2);
        }
        else if (score_matrix[i][j] - sigma(i, j) == score_matrix[i-1][j-1]) {
            return getOptAlignment(i - 1, j - 1, str1.charAt(i - 1) + suffix1, str2.charAt(j - 1) + suffix2);
        } else {
            return new String[2];
        }
    }

    // Print the score matrix if the lengths of the both strings are less than 15
    private void printScoreMatrix() {
        if (str1.length() >= 15 || str2.length() >= 15) {
            System.out.println("Strings are longer than 15; Matrix not printed");
            return;
        }

        for (int j = 0; j <= n; j++) {
            if (j == 0) {
                System.out.print("    ");
            } else {
                System.out.print(str2.charAt(j-1) + " ");
            }
        }
        System.out.println();
        for (int i = 0; i <= m; i++) {
            if (i == 0) {
                System.out.print("  ");
            } else {
                System.out.print(str1.charAt(i-1) + " ");
            }
            for (int j = 0; j <= n; j++) {
                System.out.print(score_matrix[i][j] + " ");
            }
            System.out.println();
        }
    }

    // Get the empirical p-value with the original score and the number of permutations (N)
    private double getEmpiricalPValue(int original_score, int num_of_permutations) {
        if (num_of_permutations <= 0) {
            System.out.println("N must be greater than 0.");
        }
        double k = 0;
        for (int i = 1; i <= num_of_permutations; i++) {
            SmithWaterman sw = new SmithWaterman(str1, randomSwitchS2(str2), id1, id2);
            if (sw.getMaxScore() >= original_score) {
                k++;
            }
        }
        double p_value = (k+1) / (num_of_permutations+1);
        return p_value;
    }

    // Switch the characters of the string randomly
    private String randomSwitchS2(String str2) {
        char[] arr = str2.toCharArray();
        Random rand = new Random();
        for (int i = n-1; i > 0; i--){
            int j = rand.nextInt(i + 1);
            swap(arr, i, j);
        }
        return String.valueOf(arr);
    }

    // Swap arr[i] and arr[j]
    private void swap(char[] arr, int i, int j) {
        char temp = arr[i];
        arr[i] = arr[j];
        arr[j] = temp;
    }
}
