package snvcounter.parsers;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;

import htsjdk.tribble.readers.TabixReader;
import htsjdk.tribble.readers.TabixReader.Iterator;
import snvcounter.annotator.Gene;
import snvcounter.annotator.PrintableVariant;

public class AnnotationParser implements Closeable {
	
	protected TabixReader annotationReader;
	private AnnotationType anoteType;
	protected File annotationFile;
	protected Gene gene;
	
	public AnnotationParser(File annotationFile, AnnotationType anoteType, Gene gene) throws IOException {
				
		File annotationIndex = new File(annotationFile.getAbsolutePath() + ".tbi");
		annotationReader = new TabixReader(annotationFile.getAbsolutePath(), annotationIndex.getAbsolutePath());
		this.anoteType = anoteType;
		this.annotationFile = annotationFile;
		this.gene = gene;
		
	}
	
	public void close() {
		annotationReader.close();
	}
	
	public double getAnnotation(PrintableVariant variant) throws IOException {
						
		String ref = variant.getRefBaseString();
		String alt = variant.getAltBaseString();
		int pos = variant.getStart();
		String chr = variant.getParsedChr();
		
		String line;
		String data[];
		
		double caughtScore = Double.NaN;
		
		Iterator annotationItr = annotationReader.query(chr, pos - 1, pos + 1);
			
		while((line = annotationItr.next()) != null) {
			
			data = line.split("\t");
			int tabixPos = Integer.parseInt(data[1]);
			String tabixRef = data[2];
			String tabixAlt = data[3];
			
			try {
				// Some MPC scores are NA, change to NaN so they can be parsed by my code...
				if (data[anoteType.getColumnNumber()].equals("NA")) {
					data[anoteType.getColumnNumber()] = "NaN";
				}
				double rawScore = Double.parseDouble(data[anoteType.getColumnNumber()]);
				
				int sizeDiff = Math.abs(alt.length() - ref.length());
				
				//PEXT will not match InDels. Fix here by just taking the PEXT score of the InDel's start position. Should be reasonably accurate.
				if (pos == tabixPos && ((ref.equals(tabixRef) && alt.equals(tabixAlt)) || (anoteType == AnnotationType.PEXT && sizeDiff > 0))) {
					
					// PEXT annotates on the per-gene level so match here
					if (anoteType.doGeneMatch()) {
					
						if (data[6].equals(gene.getEnsID())) {
							caughtScore = rawScore;
							break;
						}
						
					} else {
						caughtScore = rawScore;
						break;
					}
					
				}
			} catch (NumberFormatException e) {
				if (data[anoteType.getColumnNumber()].equals(".")) {
					caughtScore = 0;
				} else {
					System.err.println(anoteType.toString() + "\t" + data[anoteType.getColumnNumber()]);
					throw new NumberFormatException();
				}
			}

		}

		return caughtScore;
		
	}
		
	public enum AnnotationType {
		CADD(5, false),GNOMAD(5, false),MPC(18, false),VQSR(5, false),VEP(0, false),PEXT(5, true);
		
		private int columnNumber;
		private boolean geneMatch;
		
		private AnnotationType(int columnNumber, boolean geneMatch) {
			this.columnNumber = columnNumber;
			this.geneMatch = geneMatch;
		}
		
		public int getColumnNumber() {
			return columnNumber;
		}
		public boolean doGeneMatch() {
			return geneMatch;
		}
	}
	
}
