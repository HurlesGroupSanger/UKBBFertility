package cnvqc.merger;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

import org.apache.commons.exec.ExecuteException;

import cnvqc.utilities.Bedtools;
import cnvqc.utilities.CNV;
import cnvqc.utilities.CNVInterval;

public class ReciprocalClustering {

	private File tmpDir;
		
	public ReciprocalClustering (File tmpDir) throws IOException {
		
		this.tmpDir = tmpDir;
				
	}
	
	public Map<Integer, CNVInterval> getMergedCNVs(List<CNV> currentCNVList) throws IOException {
		Map<String, CNV> toMerge = Bedtools.printBED(currentCNVList, tmpDir);
		Map<Integer, CNVInterval> finalCNVs = mergeReciprocal(toMerge, tmpDir);
		return finalCNVs;
	}
	
	private Map<Integer, CNVInterval> mergeReciprocal(Map<String, CNV> mergedCNVs, File tmpDir) throws ExecuteException, IOException {
		
		BufferedReader bedReader = new BufferedReader(new FileReader(new File(tmpDir.getAbsolutePath() + "/merge_cnv.sorted.bed")));
		
		String line;
		String data[];

		Map<String,CNV> recovered = new ConcurrentHashMap<String,CNV>();
		
		while ((line = bedReader.readLine()) != null) {
			
			data = line.split("\t");
			
			CNV cnv = mergedCNVs.get(data[3]);
			
			recovered.put(data[3], cnv);
			
		}
		
		int mergeNum = 1;
		Map<Integer, CNVInterval> finalCNVs = new HashMap<Integer, CNVInterval>();
				
		for (Map.Entry<String, CNV> itrOneEntry : recovered.entrySet()) {
				
			String currentTopLevel = itrOneEntry.getKey();
			CNV toAdd = itrOneEntry.getValue();
			
			if (recovered.remove(currentTopLevel) != null) {

				finalCNVs.put(mergeNum, new CNVInterval(toAdd));
			
				for (Map.Entry<String, CNV> itrTwoEntry : recovered.entrySet()) {
					
					boolean allMatch = true;
					String currentItr = itrTwoEntry.getKey();
					CNV cnvOne = itrTwoEntry.getValue();
					
					for (CNV cnvTwo : finalCNVs.get(mergeNum).getAttachedCNVs()) {
					
						if (new CNVInterval(cnvOne).overlaps(new CNVInterval(cnvTwo))) {
						
							if (!calcRecip(cnvOne, cnvTwo)) {
								allMatch = false;
								break;
							}
							
						} else {
							
							allMatch = false;
							
						}
						
					}
					
					if (allMatch == true) {
						
						CNVInterval intOne = finalCNVs.get(mergeNum);
						CNVInterval intTwo = new CNVInterval(cnvOne);
						finalCNVs.put(mergeNum, intOne.union(intTwo));
						recovered.remove(currentItr);
						
					}
					
				}
				mergeNum++;
			}

		}
		
		bedReader.close();
		return finalCNVs;		
		
	}
	
	private boolean calcRecip(CNV cnvOne, CNV cnvTwo) {
		
		Double maxStart = Math.max((double)cnvOne.getStart(), (double)cnvTwo.getStart());
		Double minEnd = Math.min((double)cnvOne.getEnd(), (double)cnvTwo.getEnd());
		Double overlapLen = minEnd - maxStart;
		
		double recipOne = overlapLen / cnvOne.getLength();
		double recipTwo = overlapLen / cnvTwo.getLength();
//		System.out.println(recipOne + "\t" + recipTwo);
		if (recipOne >= 0.75 && recipTwo >= 0.75) {
			return true;
		} else {
			return false;
		}
		
	}
	
}
