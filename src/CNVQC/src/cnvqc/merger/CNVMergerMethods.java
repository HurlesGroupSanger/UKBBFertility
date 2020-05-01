package cnvqc.merger;

import java.io.BufferedWriter;
import java.io.Closeable;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import cnvqc.utilities.CNV;
import cnvqc.utilities.CNV.CopyType;
import cnvqc.utilities.CNVInterval;
import cnvqc.utilities.pathogenicannotation.PathogenicAnnotator;
import cnvqc.utilities.pathogenicannotation.PathogenicAnnotator.OverlapError;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.IntervalTree.Node;

public class CNVMergerMethods implements Closeable {

	private double epsilon;
	private VCFEngine vcfEngine;
	private File output;
	private File tmpDir;
	private List<CNV> cnvs;
	private ReciprocalClustering clusterer;
	private BufferedWriter perIndividualWriter;
	private PathogenicAnnotator rawCNVAnnotator;
	
	public CNVMergerMethods(File output, File fastaRef, File fastaIndex, Set<String> samples, Set<String> chromosomes, File tmpDir, List<CNV> cnvs, ReciprocalClustering clusterer) throws IOException {
		
		this.output = output;
		this.tmpDir = tmpDir;
		this.cnvs = cnvs;
		this.clusterer = clusterer;
		vcfEngine = new VCFEngine(new File(output.getAbsolutePath() + ".reciprocal"),
				fastaRef,
				fastaIndex,
				samples, 
				chromosomes);
		rawCNVAnnotator = new PathogenicAnnotator();
		perIndividualWriter = new BufferedWriter(new FileWriter(new File(output.getAbsolutePath() + ".qcdMerged.txt")));
				
	}
	
	@Override
	public void close() throws IOException {
		vcfEngine.close();
		perIndividualWriter.close();
	}
	
	public List<CNV> MergeCNVs() throws IOException, OverlapError {
		
		//Merge CNVs by simple overlap (just required > 0 overlap)
		BufferedWriter mergedWriter = new BufferedWriter(new FileWriter(new File(output.getAbsolutePath() + ".reciprocal.merged.bed")));
		BufferedWriter rawWriter = new BufferedWriter(new FileWriter(new File(output.getAbsolutePath() + ".reciprocal.raw_cnv.bed")));
		
		Map<String,IntervalTree<List<CNV>>> mergeIntervals = InitialMergerEngine.mergeCNVsbyOverlap(cnvs, tmpDir);
		
		mergedWriter.write("track name=Merged" + " description=\"Merged CNVs INTERVAL\" itemRgb=\"On\"\n");
		rawWriter.write("track name=Original" + " description=\"Unmerged CNVs INTERVAL\" itemRgb=\"On\"\n");
				
		List<CNV> finalCNVs = new ArrayList<CNV>();
		
		int mergeNum = 1;
		
		for (Map.Entry<String, IntervalTree<List<CNV>>> currentChrEntry : mergeIntervals.entrySet()) {
			
			String chr = currentChrEntry.getKey();
			Iterator<Node<List<CNV>>> treeItr = currentChrEntry.getValue().iterator();
			
			int totalMerged = 0;
			while (treeItr.hasNext()) {

				Node<List<CNV>> currentNode = treeItr.next();
				
				Map<Integer, CNVInterval> mergedCNVs = clusterer.getMergedCNVs(currentNode.getValue());
				
				for (Map.Entry<Integer, CNVInterval> cnvEntry : mergedCNVs.entrySet()) {
					
					CNVInterval currentInterval = cnvEntry.getValue();

					String color = checkColor(mergeNum);
										
					mergedWriter.write(chr + "\t" + currentInterval.getStart() + "\t" + currentInterval.getEnd() + "\tCNV_" + mergeNum + "\t1000\t+\t" + currentInterval.getStart() + "\t" + currentInterval.getEnd() + "\t" + color + "\n");

					Map<String, IndividualRecord> individualCopyNumbers = new HashMap<String, IndividualRecord>();
					List<CNV> currentRawCNVs = currentInterval.getAttachedCNVs();
					
					Set<String> pathogenic = rawCNVAnnotator.parsePathogenic(chr, currentInterval.getStart(), currentInterval.getEnd());
							
					String pathString = deconvolutePathogenic(pathogenic); //We should realistically only have one term here!
					
					boolean caughtDEL = false;
					boolean caughtDUP = false;
					
					for (CNV cnv : currentRawCNVs) {
					
						cnv.attachMergeGroup(mergeNum);
						finalCNVs.add(cnv);
						rawWriter.write(chr + "\t" + cnv.getStart() + "\t" + cnv.getEnd() + "\tCNV_" + mergeNum + "-" + cnv.getSampleInformation().getSplitFile().getName() + "\t1000\t+\t" + cnv.getStart() + "\t" + cnv.getEnd() + "\t" + color + "\n");
						totalMerged++;
						individualCopyNumbers.put(cnv.getSampleInformation().getSplitFile().getName(), new IndividualRecord(cnv.getCopyNumber(),cnv.getConfidence()));
						perIndividualWriter.write(cnv.getPrintable() + "\t" + chr + ":" + currentInterval.getStart() + "-" + currentInterval.getEnd() + "\t" + pathString);
						perIndividualWriter.newLine();
						
						if (cnv.getCopyType() == CopyType.DEL) {
							caughtDEL = true;
						} else {
							caughtDUP = true;
						}
						
					}
					
					perIndividualWriter.flush();
					vcfEngine.addRecord(chr, currentInterval.getStart(), currentInterval.getEnd(), individualCopyNumbers, pathString, caughtDEL && caughtDUP);
					mergeNum++;
				}
				
				mergedWriter.flush();
				rawWriter.flush();
				
			}
			
			System.err.println("stats\t" + epsilon + "\t" + chr + "\t" + totalMerged);
			
		}
		
		mergedWriter.close();
		rawWriter.flush();
		rawWriter.close();
		
		return finalCNVs;
		
	}

	public class IndividualRecord {
		
		private int copyNumber;
		private double confidence;
				
		public IndividualRecord(int copyNumber, double confidence) {
			this.copyNumber = copyNumber;
			this.confidence = confidence;
		}

		public int getCopyNumber() {
			return copyNumber;
		}
		public double getConfidence() {
			return confidence;
		}
		
	}
	
	private String checkColor(int iteration) {
		
		Map<Integer, String> color = new HashMap<Integer, String>();
		color.put(0, "255,0,0");
		color.put(1, "255,127,0");
		color.put(2, "255,255,0");
		color.put(3, "0,255,0");
		color.put(4, "0,0,255");
		color.put(5, "75,0,130");
		color.put(6, "143,0,255");
		
		int mod = iteration % 7;
		return color.get(mod);
				
	}

	private String deconvolutePathogenic(Set<String> pathogenic) {
		
		//This just sets up a heirarchy for commonly overlapped calls
		if (pathogenic.size() == 1) {
			return(pathogenic.iterator().next());
		} else {
			if (pathogenic.contains("Large")) {
				return("Large");
			} else if (pathogenic.contains("15q13.3")) {
				return("15q13.3");
			} else if (pathogenic.contains("2q13")) {
				return("2q13");
			} else {
				return(null);
			}
		}
		
		
	}
	
}
