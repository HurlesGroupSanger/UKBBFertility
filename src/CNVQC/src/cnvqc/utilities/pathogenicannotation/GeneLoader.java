package cnvqc.utilities.pathogenicannotation;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Map;

import htsjdk.samtools.util.IntervalTree;

public class GeneLoader {

	Map<String, IntervalTree<String>> geneIntervals;
	
	public GeneLoader() throws IOException {
		
		BufferedReader geneReader = new BufferedReader(new InputStreamReader(this.getClass().getResourceAsStream("/cnvqc/utilities/pathogenicannotation/geneList_hg19.bed"), "UTF-8"));
		geneIntervals = BuildIntervalTree(geneReader);
		geneReader.close();
				
	}
	
	public Map<String, IntervalTree<String>> getGeneIntervals() {
		return geneIntervals;
	}
	
	private Map<String, IntervalTree<String>> BuildIntervalTree(BufferedReader geneReader) throws IOException {
		Map<String, IntervalTree<String>> genes = new HashMap<String, IntervalTree<String>>();
		
		String line;
		String data[];
		
		while ((line = geneReader.readLine()) != null) {
			data = line.split("\t");
			String chr = data[0];
			int start = Integer.parseInt(data[1]);
			int stop = Integer.parseInt(data[2]);
			String geneName = data[3];
			
			if (genes.containsKey(chr)) {
				
				genes.get(chr).put(start, stop, geneName);
				
			} else {
				
				IntervalTree<String> currentChr = new IntervalTree<String>();
				currentChr.put(start, stop, geneName);
				genes.put(chr, currentChr);
				
			}
		}
				
		return genes;
	}
	
}
