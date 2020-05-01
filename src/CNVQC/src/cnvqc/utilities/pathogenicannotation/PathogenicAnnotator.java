package cnvqc.utilities.pathogenicannotation;

import java.io.IOException;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

import cnvqc.utilities.pathogenicannotation.PathogenicLoader.PathogenicLocus;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.IntervalTree.Node;

public class PathogenicAnnotator {

	private Map<String, IntervalTree<String>> geneIntervals;
	private PathogenicLoader pathogenic;

	public PathogenicAnnotator () throws IOException {
		GeneLoader gene = new GeneLoader();
		geneIntervals = gene.getGeneIntervals();
		pathogenic = new PathogenicLoader();
	}

	public Set<String> parsePathogenic(String chr, int start, int end) throws OverlapError {

		Set<String> genes = parseGenes(chr, start, end);
		
		Iterator<Node<PathogenicLocus>> pathItr;
		pathItr = pathogenic.getPathogenicIntervals().get(chr).overlappers(start, end);
		Set<String> loci = new HashSet<String>();
		
		while (pathItr.hasNext()) {
			Node<PathogenicLocus> currentNode = pathItr.next();
			int locusStart = currentNode.getStart();
			int locusEnd = currentNode.getEnd();
			PathogenicLocus currentLocus = currentNode.getValue();
			int length = currentLocus.getLength();
			int cnvLength = end - start;
			
			double overlap;
			double bpOverlap;
			if (start > locusStart && end < locusEnd) {
				bpOverlap = cnvLength;
				overlap = bpOverlap / length;
			} else if (start < locusStart && end > locusEnd) {
				bpOverlap = length;
				overlap = 1.0;
			} else if (end > locusEnd) {
				bpOverlap = locusEnd - start;
				overlap = bpOverlap / length;
			} else if (start < locusStart) {
				bpOverlap = end - locusStart;
				overlap = bpOverlap / length;
			} else if (start == locusStart && end == locusEnd) {
				bpOverlap = length;
				overlap = 1.0;
			} else if (start == locusStart && end != locusEnd) {
				if (end > locusEnd) {
					bpOverlap = length;
					overlap = 1.0;
				} else {
					bpOverlap = cnvLength;
					overlap = bpOverlap / length;
				}
			} else if (start != locusStart && end == locusEnd) {
				if (start < locusStart) {
					bpOverlap = length;
					overlap = 1.0;
				} else {
					bpOverlap = cnvLength;
					overlap = bpOverlap / length;
				}
			} else {
				bpOverlap = 0;
				overlap = 0;
				throw new OverlapError();
			}

			boolean foundGenes;
			boolean foundExons;
			String locusName = currentLocus.getName();
			switch (currentLocus.getPathogenicType()) {
			case FIFTY:
				if (overlap > 0.50) {
					loci.add(locusName);
				}
				break;
			case ONEHUNDRED:
				if (overlap > 0.95) {
					loci.add(locusName);
				}
				break;
			case GENEREQUIRED:
				foundGenes = CheckGenes(currentLocus.getRequiredGenes(), genes);
				if (overlap > 0.50 && foundGenes) {
					loci.add(locusName);
				}
				break;
			case MBGENE:
				foundGenes = CheckGenes(currentLocus.getRequiredGenes(), genes);
				if (cnvLength >= 1000000 && foundGenes) {
					loci.add(locusName);
				}
				break;
			case EXONS:
				foundExons = CheckExons(currentLocus.getExons(), start, end);
				StringBuilder locus = new StringBuilder();
				if (foundExons) {
					locus.append(locusName + "DEL");
				}
				foundGenes = CheckGenes(currentLocus.getRequiredGenes(), genes);
				if (overlap > 0.95) {
					if (locus.length() > 0) {
						locus.append("|" + locusName + "DUP");
					} else {
						locus.append(locusName + "DUP");
					}
				}
				if (locus.length() > 0) {
					loci.add(locus.toString());
				}
				break;
			case MB:
				if (bpOverlap >= 1000000) {
					loci.add(locusName);
				}
				break;
			case SEGDUP:
				if (bpOverlap >= 1000000) {
					loci.add(locusName);
				}
				break;
			default:
				break;
			}
		}
		if ((end - start) > 25000000) {
			loci.add("Large");
		}
		return loci;
		
	}
 	public Set<String> parseGenes(String chr, int start, int end) {
 		Set<String> genes = new HashSet<String>();
		IntervalTree<String> currentTree = geneIntervals.get(chr);
		Iterator<Node<String>> geneHits = currentTree.overlappers(start, end);
		while (geneHits.hasNext()) {
			Node<String> currentNode = geneHits.next();
			genes.add(currentNode.getValue());
		}
		return genes;
 	}
	
 	private boolean CheckGenes (Set<String> requiredGenes, Set<String> foundGenes) {
 		
 		int toFind = requiredGenes.size();
 		int found = 0;
 		for (String g : requiredGenes) {
			if (foundGenes.contains(g)) {
				found++;
			}
		}
 		return toFind == found;
 		
 	}
 	private boolean CheckExons (IntervalTree<Integer> exons, int start, int end) {
	
 		Iterator<Node<Integer>> foundExons = exons.overlappers(start, end);
 		int totalFound = 0;
 		while (foundExons.hasNext()) {
 			foundExons.next();
 			totalFound++;
 		}
 		return totalFound > 0;
 		
 	}
 	
 	public class OverlapError extends Exception {

		private static final long serialVersionUID = 1L;
 		
		public OverlapError() {
			super("Overlap possibility not found!");
		}
 		
 		
 	}
 	
}
