package snvcounter.parsers;

import java.io.File;
import java.io.IOException;

import org.apache.commons.math3.fraction.Fraction;

import htsjdk.tribble.readers.TabixReader.Iterator;
import snvcounter.annotator.Gene;
import snvcounter.annotator.PrintableVariant;

public class VEPAnnotator extends AnnotationParser {
		
	// This is a special extension of my default annotation parser specifically for VEP annotation
	public VEPAnnotator(File annotationFile, AnnotationType anoteType, Gene gene) throws IOException {
		super(annotationFile, anoteType, gene);
	}
	
	public VEPAnnotation getVEPAnnotation (PrintableVariant variant) throws IOException, InterruptedException {
				
		String ref = variant.getRefBaseString();
		String alt = variant.getAltBaseString();
		int pos = variant.getStart();
		String chr = variant.getParsedChr();
		
		
		VEPVariant originalVar = new VEPVariant(ref, alt, pos);
		VEPVariant adjustedVar = NormalizeIndel(ref, alt, pos);
					
		VEPAnnotation orig = IterateOverVariants(originalVar, chr);
				
		if (orig.getCsq() == Consequence.FAIL) {
			// Try the alternate annotation b/c it is likely to be an InDel
			VEPAnnotation rescue = IterateOverVariants(adjustedVar, chr);
			return rescue;
		} else {
			return orig;
		}
		
	}

	private VEPAnnotation IterateOverVariants (VEPVariant var, String chr) throws NumberFormatException, IOException {
		
		Iterator annotationItr = annotationReader.query(chr, var.getPosition() - 1, var.getPosition() + 1);
		
		String line;
		String data[];
						
		Consequence canonicalCsq = Consequence.FAIL;
		boolean isLastExon = false;
		boolean isLastIntron = false;
		boolean foundVariant = false;
		
		while((line = annotationItr.next()) != null) {
			
			data = line.split("\t");
			
			// Deidentify all information:
			// These columns MUST be correct in the VEP annotation that we create
			int tabixPos = Integer.parseInt(data[1]);
			String tabixRef = data[2];
			String tabixAlt = data[3];
			String vepGeneID = data[4];
			boolean isCanonical = data[6].matches("YES");
			String consequence = data[7];
			String loftee = data[8];
			String exonContext = data[9];
			String intronContext = data[10];
			
			// Don't want non-coding transcripts...
			if (consequence.contains("NMD_transcript_variant") || consequence.contains("non_coding_transcript")) {continue;}
			
			if (var.getPosition() == tabixPos && var.getRef().equals(tabixRef) && var.getAlt().equals(tabixAlt)) {			

				if (vepGeneID.equals(gene.getEnsID())) {
					
					// This is to check and see if we did find the variant but did not find a canonical transcript.
					foundVariant = true;
					
					Consequence currentCsq;
					
					if (consequence.contains("synonymous_variant")) {
						currentCsq = Consequence.SYN;
					} else if (consequence.contains("splice_acceptor_variant") || 
							consequence.contains("splice_donor_variant") ||
							consequence.contains("stop_gained") || 
							consequence.contains("frameshift_variant")) {
						if (loftee.equals("HC")) {
							currentCsq = Consequence.LOF_HC;
						} else if (loftee.equals("LC")) {
							currentCsq = Consequence.LOF_LC;
						} else {
							currentCsq = Consequence.LOF_NONCODING; // Happens when a non-coding transcript. Also happens in very rare cases for seemingly normal transcripts
						}
					} else if (consequence.contains("missense_variant")) {
						currentCsq = Consequence.MIS;
					} else if (consequence.contains("inframe")) {
						currentCsq = Consequence.INFRAME;
					} else {
						currentCsq = Consequence.NONE;
					}
					
					// Only care about exon positioning if canonical transcript
					if (isCanonical) {
						
						// Also set the canonical consequence flag
						canonicalCsq = currentCsq;
						
						if (exonContext.equals(".")) {
							isLastExon = false;
						} else {
							String exonParsed[] = exonContext.split("/");
							String multiExon[] = exonParsed[0].split("-"); //This is necessary as occasionally large InDel variants can span multiple introns/exons
							Fraction exonFrac = Fraction.getReducedFraction(Integer.parseInt(multiExon[0]), Integer.parseInt(exonParsed[1]));
							
							if (exonFrac.doubleValue() == 1) {
								isLastExon = true;
							}
						}
						if (intronContext.equals(".")) {
							isLastIntron = false;
						} else {
							String intronParsed[] = intronContext.split("/");
							String multiIntron[] = intronParsed[0].split("-");
							Fraction intronFrac = Fraction.getReducedFraction(Integer.parseInt(multiIntron[0]), Integer.parseInt(intronParsed[1]));
							
							if (intronFrac.doubleValue() == 1) {
								isLastIntron = true;
							}
						}
					}
										
				}
			}
		}

		// Set a flag here if we actually found the canonical transcript for this particular variant/gene pair.
		if (foundVariant == true && canonicalCsq == Consequence.FAIL) {
			canonicalCsq = Consequence.NOCAN;
		}
		
		return new VEPAnnotation(canonicalCsq, isLastExon, isLastIntron);
		
	}
		
	// We need to left normalize InDels so that they can be annotated properly
	private VEPVariant NormalizeIndel (String ref, String alt, int position) {
		
		// Deletion - take 1 base off the ref
		if (ref.length() > alt.length()) {

			if (alt.length() > 1) {
			
				ref = ref.substring(1, ref.length());
				alt = alt.substring(1, alt.length());
				position = position+=1;
				
			} else {
				ref = ref.substring(1, ref.length());
				alt = "-";
				position = position+=1;
			}
				
		// Insertion - take 1 base off the alt
		} else if (ref.length() < alt.length()) {
			
			if (ref.length() > 1) {
				
				ref = ref.substring(1, ref.length());
				alt = alt.substring(1, alt.length());
				position = position+=1;
				
			} else {
				ref = "-";
				alt = alt.substring(1, alt.length());
			}
			
		}
		// Do nothing if they are equal, is a SNV and do not need to normalize
		
		return(new VEPVariant(ref, alt, position));
		
	}
	
	// Class holds my relevant variant-level information
	private class VEPVariant {
		
		private String ref;
		private String alt;
		private int position;
		
		private VEPVariant(String ref, String alt, int position) {
			
			this.ref = ref;
			this.alt = alt;
			this.position = position;
			
		}

		public String getRef() {
			return ref;
		}

		public String getAlt() {
			return alt;
		}

		public int getPosition() {
			return position;
		}
		
	}
	
	// Class holds actual annotation that we care about:
	public class VEPAnnotation {
		
		Consequence csq;
		boolean isLastExon;
		boolean isLastIntron;
		
		public VEPAnnotation(Consequence csq, boolean isLastExon, boolean isLastIntron) {
			
			this.csq = csq;
			this.isLastExon = isLastExon;
			this.isLastIntron = isLastIntron;
			
		}

		public Consequence getCsq() {
			return csq;
		}
		
		public boolean isRelevant() {
			if (csq == Consequence.LOF_NONCODING || csq == Consequence.LOF_LC || csq == Consequence.NOCAN || csq == Consequence.NONE || csq == Consequence.FAIL) {
				return false;
			} else {
				return true;
			}
		}

		public boolean isLastExon() {
			return isLastExon;
		}
		
		public boolean isLastIntron() {
			return isLastIntron;
		}
	
	}
	
	// Possible VEP consequences that are relevant to the project
	public enum Consequence {
		
		SYN,LOF_HC,LOF_LC,LOF_NONCODING,MIS,INFRAME,NONE,NOCAN,FAIL;
				
	}
	
}
