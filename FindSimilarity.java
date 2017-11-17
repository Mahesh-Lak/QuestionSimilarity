import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;

import javax.lang.model.type.IntersectionType;

import org.apache.lucene.analysis.Analyzer;
import org.apache.lucene.analysis.TokenStream;
import org.apache.lucene.analysis.core.StopAnalyzer;
import org.apache.lucene.analysis.core.StopFilter;
import org.apache.lucene.analysis.en.PorterStemFilter;
import org.apache.lucene.analysis.standard.StandardAnalyzer;
import org.apache.lucene.analysis.tokenattributes.CharTermAttribute;

import com.opencsv.CSVReader;
import com.opencsv.CSVWriter;

class QuestionSet{
	String questionId;
	String questionOne;
	String questionTwo;
	
	void questionSetter(String questionId, String questionOne, String questionTwo) {
		this.questionId = questionId;
		this.questionOne = questionOne;
		this.questionTwo = questionTwo;
	}
	String getQuestionId() {
		return questionId;
	}
	String getQuestionOne() {
		return questionOne;
	}
	String getQuestionTwo() {
		return questionTwo;
	}
}

public class CosineSimilarity {
	public static void main(String args[]) throws Exception {
		CSVReader reader = new CSVReader(new FileReader("/Users/test.csv"));
		List<String[]> csvEntries = new ArrayList<String[]>();
		csvEntries.add(new String[] {"id", "is_duplicate"});
		String [] nextLine;
		int i=0;
	      while ((nextLine = reader.readNext()) != null) {
	    	if(!nextLine[0].equals("id")) {
		        String questionId = nextLine[0];
		        String questionOne = stem(nextLine[3].toLowerCase());
		        String questionTwo = stem(nextLine[4].toLowerCase());
		        Double cosValue = cosineSimilarity(questionOne, questionTwo);
		        Integer isDuplicate;
		        
//		        double NGRAMValue = NGramDistance(questionOne, questionTwo);
//		        double LevenValue = LevensteinDistance(questionOne, questionTwo);
		        double jaccardValue = jaccardDistance(questionOne, questionTwo);
		        System.out.println(++i+"  "+jaccardValue);
		        if(jaccardValue < 0.40) {
		        	isDuplicate = 1;
		        }
		        else {
		        	isDuplicate = 0;
		        }
//		        if(cosValue> 0.60) {
//		        	
//		        	if(NGRAMValue < 0.50) {
//		        		isDuplicate = 1;
//		        	}
//		        }
//		        else {
//		        	isDuplicate= 0;
//		        }
		        csvEntries.add(new String[] {questionId, isDuplicate.toString()});
	    	}
	      }
	      writeToCSV(csvEntries);
	      
	}
	/**
	 * To eliminate stop words and perform stemming
	 *
	 */
	@SuppressWarnings("resource")
	public static String stem(String term) throws Exception {
	    Analyzer analyzer = new StandardAnalyzer();
	    String resultStemStop = "";
	    TokenStream result = analyzer.tokenStream(null, term);
	    result = new PorterStemFilter(result);
	    result = new StopFilter(result, StopAnalyzer.ENGLISH_STOP_WORDS_SET);
	    CharTermAttribute resultAttr = result.addAttribute(CharTermAttribute.class);
	    result.reset();
	    while (result.incrementToken()) {
	    	resultStemStop = resultStemStop + " " + resultAttr.toString();
	    }
	    return resultStemStop;
	}
	public static Map<String, Integer> getTermFrequencyMap(String[] terms) {
	        Map<String, Integer> termFrequencyMap = new HashMap<>();
	        for (String term : terms) {
	            Integer n = termFrequencyMap.get(term);
	            n = (n == null) ? 1 : ++n;
	            termFrequencyMap.put(term, n);
	        }
	        return termFrequencyMap;
	 }

	    /**
	     * @param text1 
	     * @param text2 
	     * @return cosine similarity of text1 and text2
	     */
	    public static double cosineSimilarity(String text1, String text2) {
	        //Get vectors
	        Map<String, Integer> a = getTermFrequencyMap(text1.split("\\W+"));
	        Map<String, Integer> b = getTermFrequencyMap(text2.split("\\W+"));

	        //Get unique words from both sequences
	        HashSet<String> intersection = new HashSet<>(a.keySet());
	        intersection.retainAll(b.keySet());

	        double dotProduct = 0, magnitudeA = 0, magnitudeB = 0;

	        //Calculate dot product
	        for (String item : intersection) {
	            dotProduct += a.get(item) * b.get(item);
	        }

	        //Calculate magnitude a
	        for (String k : a.keySet()) {
	            magnitudeA += Math.pow(a.get(k), 2);
	        }

	        //Calculate magnitude b
	        for (String k : b.keySet()) {
	            magnitudeB += Math.pow(b.get(k), 2);
	        }

	        //return cosine similarity
	        return dotProduct / Math.sqrt(magnitudeA * magnitudeB);
	    }
	    static void writeToCSV(List<String[]> data) throws Exception{
	    	String csv = "C:/result.csv";
	    	CSVWriter writer = new CSVWriter(new FileWriter(csv));
	    	writer.writeAll(data); 
	    	writer.close();
	    }
	    static public final double NGramDistance(final String s0, final String s1) {
	    	int n=2;
	        if (s0 == null) {
	            throw new NullPointerException("s0 must not be null");
	        }

	        if (s1 == null) {
	            throw new NullPointerException("s1 must not be null");
	        }

	        if (s0.equals(s1)) {
	            return 0;
	        }

	        final char special = '\n';
	        final int sl = s0.length();
	        final int tl = s1.length();

	        if (sl == 0 || tl == 0) {
	            return 1;
	        }

	        int cost = 0;
	        if (sl < n || tl < n) {
	            for (int i = 0, ni = Math.min(sl, tl); i < ni; i++) {
	                if (s0.charAt(i) == s1.charAt(i)) {
	                    cost++;
	                }
	            }
	            return (float) cost / Math.max(sl, tl);
	        }

	        char[] sa = new char[sl + n - 1];
	        float[] p; //'previous' cost array, horizontally
	        float[] d; // cost array, horizontally
	        float[] d2; //placeholder to assist in swapping p and d

	        //construct sa with prefix
	        for (int i = 0; i < sa.length; i++) {
	            if (i < n - 1) {
	                sa[i] = special; //add prefix
	            } else {
	                sa[i] = s0.charAt(i - n + 1);
	            }
	        }
	        p = new float[sl + 1];
	        d = new float[sl + 1];

	        // indexes into strings s and t
	        int i; // iterates through source
	        int j; // iterates through target

	        char[] t_j = new char[n]; // jth n-gram of t

	        for (i = 0; i <= sl; i++) {
	            p[i] = i;
	        }

	        for (j = 1; j <= tl; j++) {
	            //construct t_j n-gram
	            if (j < n) {
	                for (int ti = 0; ti < n - j; ti++) {
	                    t_j[ti] = special; //add prefix
	                }
	                for (int ti = n - j; ti < n; ti++) {
	                    t_j[ti] = s1.charAt(ti - (n - j));
	                }
	            } else {
	                t_j = s1.substring(j - n, j).toCharArray();
	            }
	            d[0] = j;
	            for (i = 1; i <= sl; i++) {
	                cost = 0;
	                int tn = n;
	                //compare sa to t_j
	                for (int ni = 0; ni < n; ni++) {
	                    if (sa[i - 1 + ni] != t_j[ni]) {
	                        cost++;
	                    } else if (sa[i - 1 + ni] == special) {
	                        //discount matches on prefix
	                        tn--;
	                    }
	                }
	                float ec = (float) cost / tn;
	                // minimum of cell to the left+1, to the top+1,
	                // diagonally left and up +cost
	                d[i] = Math.min(
	                        Math.min(d[i - 1] + 1, p[i] + 1), p[i - 1] + ec);
	            }
	            // copy current distance counts to 'previous row' distance counts
	            d2 = p;
	            p = d;
	            d = d2;
	        }

	        // our last action in the above loop was to switch d and p, so p now
	        // actually has the most recent cost counts
	        return p[sl] / Math.max(tl, sl);
	    }
	    public static double LevensteinDistance(final String s1, final String s2) {
	        if (s1 == null) {
	            throw new NullPointerException("s1 must not be null");
	        }

	        if (s2 == null) {
	            throw new NullPointerException("s2 must not be null");
	        }

	        if (s1.equals(s2)) {
	            return 0;
	        }

	        if (s1.length() == 0) {
	            return s2.length();
	        }

	        if (s2.length() == 0) {
	            return s1.length();
	        }

	        // create two work vectors of integer distances
	        int[] v0 = new int[s2.length() + 1];
	        int[] v1 = new int[s2.length() + 1];
	        int[] vtemp;

	        // initialize v0 (the previous row of distances)
	        // this row is A[0][i]: edit distance for an empty s
	        // the distance is just the number of characters to delete from t
	        for (int i = 0; i < v0.length; i++) {
	            v0[i] = i;
	        }

	        for (int i = 0; i < s1.length(); i++) {
	            // calculate v1 (current row distances) from the previous row v0
	            // first element of v1 is A[i+1][0]
	            //   edit distance is delete (i+1) chars from s to match empty t
	            v1[0] = i + 1;

	            // use formula to fill in the rest of the row
	            for (int j = 0; j < s2.length(); j++) {
	                int cost = 1;
	                if (s1.charAt(i) == s2.charAt(j)) {
	                    cost = 0;
	                }
	                v1[j + 1] = Math.min(
	                        v1[j] + 1,              // Cost of insertion
	                        Math.min(
	                                v0[j + 1] + 1,  // Cost of remove
	                                v0[j] + cost)); // Cost of substitution
	            }

	            // copy v1 (current row) to v0 (previous row) for next iteration
	            //System.arraycopy(v1, 0, v0, 0, v0.length);

	            // Flip references to current and previous row
	            vtemp = v0;
	            v0 = v1;
	            v1 = vtemp;

	        }

	        return v0[s2.length()];
	    }
	    public static final double jaccardSimilarity(final String s1, final String s2) {

	        if (s1 == null) {
	            throw new NullPointerException("s1 must not be null");
	        }

	        if (s2 == null) {
	            throw new NullPointerException("s2 must not be null");
	        }

	        if (s1.equals(s2)) {
	            return 1;
	        }

	        Map<String, Integer> profile1 = getProfile(s1);
	        Map<String, Integer> profile2 = getProfile(s2);

	        Set<String> union = new HashSet<String>();
	        union.addAll(profile1.keySet());
	        union.addAll(profile2.keySet());

	        int inter = 0;

	        for (String key : union) {
	            if (profile1.containsKey(key) && profile2.containsKey(key)) {
	                inter++;
	            }
	        }

	        return 2.0 * inter / (profile1.size() + profile2.size());
	    }
	    public static final Map<String, Integer> getProfile(final String string) {
	    	int k = 3;
	    	final Pattern SPACE_REG = Pattern.compile("\\s+");
	        HashMap<String, Integer> shingles = new HashMap<String, Integer>();

	        String string_no_space = SPACE_REG.matcher(string).replaceAll(" ");
	        for (int i = 0; i < (string_no_space.length() - k + 1); i++) {
	            String shingle = string_no_space.substring(i, i + k);
	            Integer old = shingles.get(shingle);
	            if (old != null) {
	                shingles.put(shingle, old + 1);
	            } else {
	                shingles.put(shingle, 1);
	            }
	        }

	        return Collections.unmodifiableMap(shingles);
	    }
	    public static final double jaccardDistance(final String s1, final String s2) {
	        return 1.0 - jaccardSimilarity(s1, s2);
	    }
	    


	 }
