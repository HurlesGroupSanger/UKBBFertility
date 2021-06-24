package snvcounter.annotator;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Map.Entry;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import snvcounter.parsers.AnnotationParser.AnnotationType;
import snvcounter.parsers.VEPAnnotator.VEPAnnotation;
import snvcounter.utilities.Combine;
import snvcounter.utilities.SNVCounterOptions.GenomeVersion;

public class PrintableVariant extends VariantContext {

	private static final long serialVersionUID = 1L;
	
	private VEPAnnotation vepAnnotation;
	private String parsedChr;
	private String refString;
	private String altString;
	private long adjustedStart;
	private double cadd;
	private double mpc;
	private double gnomadAF;
	private double VQSR;
	private double pextScore;
	private int altNumber;
	private static DecimalFormat df;
		
	// Is a class for holding a single variant (i.e. 1 row in a VCF file) and all relevant annotation.
	public PrintableVariant(VariantContext variantContext, GenomeVersion genomeVers, int altNumber) throws Exception {
		
		super(variantContext);
		parsedChr = checkChrAnnotation(genomeVers);
		this.altNumber = altNumber;
		this.refString = this.getReference().getBaseString();
		this.altString = this.getAlternateAllele(altNumber).getBaseString();
		this.adjustedStart = this.start;
		df = new DecimalFormat("#");
			
	}
	public void leftCorrectVariant() {
		
		if (this.refString.length() == this.altString.length()) {
			if (this.refString.length() > 1) {
				fixEqualLengthSNV();
			}
		} else {
			if (this.refString.length() != 1 && this.altString.length() != 1) {
				leftCorrectInDel();
			}
		}
		
	}
	// Is called after annotating VEP and VQSR. 
	// Fixes SNVs that were deemed multiallelic with an InDel and thus have are 2+ bases long with only one actual change
	private void fixEqualLengthSNV() {
		
		String ref = this.refString;
		String alt = this.altString;
		long currentPos = this.getStart();
		
		char[] refVec = ref.toCharArray();
		char[] altVec = alt.toCharArray();
		
		for (int x = 0; x < refVec.length; x++) {
			
			if (refVec[x] != altVec[x]) {
				ref = String.valueOf(refVec[x]);
				alt = String.valueOf(altVec[x]);
				currentPos += x;
				break;
			}				
		}
		
		this.refString = ref;
		this.altString = alt;	
		this.adjustedStart = currentPos;
		
	}
	private void leftCorrectInDel() {
		
		char[] refVec = this.refString.toCharArray();
		char[] altVec = this.altString.toCharArray();
		long currentPos = this.getStart();
		
		StringBuilder refBuild = new StringBuilder();
		StringBuilder altBuild = new StringBuilder();
		
		char lastMatch = Character.MIN_VALUE;
		
		// Deletions
		if (refVec.length > altVec.length) {
			
			for (int refLoc = 0, altLoc = 0; refLoc < refVec.length; refLoc++) {
				
				if (refLoc >= altVec.length) {
					if (refBuild.length() == 0) {
						refBuild.append(lastMatch);
						altBuild.append(lastMatch);
						currentPos += (refLoc - 1);
					}
					refBuild.append(refVec[refLoc]);	

				} else if (refVec[refLoc] == altVec[altLoc]) {
					lastMatch = refVec[refLoc];
					altLoc++;
				} else {
					
					if (refBuild.length() == 0) {
						refBuild.append(lastMatch);
						altBuild.append(lastMatch);
						currentPos += (refLoc - 1);
					}
					refBuild.append(refVec[refLoc]);				
				}
			}
						
		// Insertions	
		} else {
			for (int refLoc = 0, altLoc = 0; altLoc < altVec.length; altLoc++) {
				
				if (altLoc >= refVec.length) {
					if (altBuild.length() == 0) {
						refBuild.append(lastMatch);
						altBuild.append(lastMatch);
						currentPos += (altLoc - 1);
					}
					altBuild.append(altVec[altLoc]);			

				} else if (refVec[refLoc] == altVec[altLoc]) {
					lastMatch = altVec[altLoc];
					refLoc++;
				} else {
					
					if (altBuild.length() == 0) {
						refBuild.append(lastMatch);
						altBuild.append(lastMatch);
						currentPos += (altLoc - 1);
					}
					altBuild.append(altVec[altLoc]);			
				}
			}
		}
		
		this.refString = refBuild.toString();
		this.altString = altBuild.toString();
		this.adjustedStart = currentPos;
				
	}
	
	// Just check if passes VQSR
	public boolean isPrintable() {
		
		//First term is SNVs, apply SNV VQSR; Second term is InDels, apply InDel VQSR			
		if ((this.getRefBaseString().length() == this.getAltBaseString().length() && VQSR >= -3.6037) || 
				(this.getRefBaseString().length() != this.getAltBaseString().length() && VQSR >= -5.7733)) { 
			
			return(true);
		} else {
			return(false);
		}
		
	}
	
	// Just format the variant for proper print format.
	public String printableVariant(Entry<String, Integer> sampleEntry, double AC, double AN, int AP, double AF, String geneID) {
		
		List<String> printableString = new ArrayList<String>();
		
		printableString.add(sampleEntry.getKey());
		printableString.add(sampleEntry.getValue().toString());
		printableString.add(parsedChr);
		printableString.add(Integer.toString(this.getStart()));
		printableString.add(this.getRefBaseString());
		printableString.add(this.getAltBaseString());
		
		printableString.add(printDoubleScore(cadd, AnnotationType.CADD));
		printableString.add(printDoubleScore(mpc, AnnotationType.MPC));
		printableString.add(printDoubleScore(gnomadAF, AnnotationType.GNOMAD));
		printableString.add(printDoubleScore(VQSR, AnnotationType.VQSR));
		printableString.add(printDoubleScore(pextScore, AnnotationType.PEXT));
		
		printableString.add(vepAnnotation.getCsq().toString());
		printableString.add(df.format(AC));
		printableString.add(df.format(AN));
		printableString.add(String.format("%.3e", AF));
		printableString.add(Integer.toString(AP));
		printableString.add(geneID);
		printableString.add(Boolean.toString(vepAnnotation.isLastExon()));
		printableString.add(Boolean.toString(vepAnnotation.isLastIntron()));
		
		return(Combine.combineList(printableString, "\t"));
		
	}
	
	// Just formats scores for printing and checks if we need to query CADD files to recover
	private String printDoubleScore(double score, AnnotationType annoteType) {
		
		String printableString;
		
		if (Double.isNaN(score)) {
			
			// I think allow NaN? Can adjust for each annotation type. 
			// CADD should ALWAYS be present
			// MPC is going to be unrecoverable I think
			// PEXT should ALWAYS be present
			// GnomadAF can definitely be NaN
			// VQSR should ALWAYS be present
			printableString = Double.toString(score);
				
		} else {
			printableString = Double.toString(score);
		}
		
		return(printableString);
		
	}

	// Holders for various basic information
 	public String getParsedChr() {
		return parsedChr;
	}
	public String getRefBaseString() {
		return(refString);	
	}
	public String getAltBaseString() {
		return(altString);
	}
	
	@Override
	public int getStart() {
		if (adjustedStart == this.start) {
			return super.getStart();
		} else {
			return (int)adjustedStart;
		}
		
	}
	public int getAltNumber() {
		return altNumber;
	}
	public Allele getAltAllele() {
		return this.getAlternateAllele(altNumber);
	}
	public boolean isRelevant() {
		return vepAnnotation.isRelevant();
	}
	
	// For setting various annotations
	public void setVepAnnotation(VEPAnnotation vepAnnotation) {
		this.vepAnnotation = vepAnnotation;
	}
	public void setCadd(double cadd) {
		this.cadd = cadd;
	}
	public void setMpc(double mpc) {
		this.mpc = mpc;
	}
	public void setGnomadAF(double gnomadAF) {
		this.gnomadAF = gnomadAF;
	}
	public void setVQSR(double vQSR) {
		VQSR = vQSR;
	}
	public void setPextScore(double pextScore) {
		this.pextScore = pextScore;
	}
	
	// Private helpers to the constructor
	// to check to make sure 'chr' annotation in chromosome names matches
	private String checkChrAnnotation(GenomeVersion genomeVers) throws Exception {
		String variantChr;
		if (genomeVers == GenomeVersion.HG38) {
			variantChr = parseChr(this.getContig());
		} else {
			variantChr = this.getContig();
		}
		return(variantChr);
	}
	private String parseChr(String variantChr) throws Exception {
		Matcher chrMatcher = Pattern.compile("chr(\\S+)").matcher(variantChr);
		if (chrMatcher.matches()) {
			variantChr = chrMatcher.group(1);
		} else {
			throw new Exception("Chromosome " + variantChr + " not parseable");
		}
		
		return variantChr;
	}
	
}
