package cnvqc.merger;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.exec.ExecuteException;

import cnvqc.utilities.Bedtools;
import cnvqc.utilities.CNV;
import htsjdk.samtools.util.IntervalTree;

public class InitialMergerEngine {
	
	public static Map<String,IntervalTree<List<CNV>>> mergeCNVsbyOverlap(List<CNV> filteredCNVs, File tmpDir) throws IOException {
		
		int totalAdded = filteredCNVs.size();
		Map<String, CNV> mergedCNVs = Bedtools.printBED(filteredCNVs, tmpDir);
		Map<String,IntervalTree<List<CNV>>> mergedClusters = mergeOverlap(mergedCNVs, tmpDir);
		System.err.println("bedtools merge complete. Total CNVs added : " + totalAdded);
		 
		return mergedClusters;
		
	}
		
	private static Map<String,IntervalTree<List<CNV>>> mergeOverlap(Map<String, CNV> mergedCNVs, File tmpDir) throws ExecuteException, IOException {
		
		Bedtools.runBedTools("bedtools merge -i " + tmpDir.getAbsolutePath() + "/merge_cnv.sorted.bed -c 4 -o collapse", new File(tmpDir.getAbsolutePath() + "/cnv.out.bed"));		
		
		BufferedReader bedReader = new BufferedReader(new FileReader(new File(tmpDir.getAbsolutePath() + "/cnv.out.bed")));
		
		String line;
		
		Map<String,IntervalTree<List<CNV>>> cnvIntervals = new HashMap<String,IntervalTree<List<CNV>>>();
		
		while ((line = bedReader.readLine()) != null) {
			
			String data[] = line.split("\t");
			String samples[] = data[3].split(",");
			
			List<CNV> finalCNVs = new ArrayList<CNV>();
			
			for (String sample : samples) {
				finalCNVs.add(mergedCNVs.get(sample));
			}
			
			if (cnvIntervals.containsKey(data[0])) {
				cnvIntervals.get(data[0]).put(Integer.parseInt(data[1]), Integer.parseInt(data[2]), finalCNVs);
			} else {
				IntervalTree<List<CNV>> currentTree = new IntervalTree<List<CNV>>();
				currentTree.put(Integer.parseInt(data[1]), Integer.parseInt(data[2]), finalCNVs);
				cnvIntervals.put(data[0], currentTree);
			}
			
		}
		
		bedReader.close();
		
//		new File(tmpDir.getAbsolutePath() + "/merge_cnv.in.bed").delete();
//		new File(tmpDir.getAbsolutePath() + "/merge_cnv.sorted.bed").delete();
//		new File(tmpDir.getAbsolutePath() + "/merge_cnv.out.bed").delete();
		
		return cnvIntervals;
	}
	
}
