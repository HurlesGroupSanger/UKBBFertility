package snvcounter.annotator;

import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Map.Entry;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import snvcounter.utilities.SNVCounterOptions.GenomeVersion;
import snvcounter.parsers.AnnotationParser;
import snvcounter.parsers.AnnotationParser.AnnotationType;
import snvcounter.parsers.VEPAnnotator.VEPAnnotation;
import snvcounter.utilities.Combine;

public class PrintableVariant extends VariantContext {

	private static final long serialVersionUID = 1L;
	
	private VEPAnnotation vepAnnotation;
	private String parsedChr;
	private String refString;
	private String altString;
	private double cadd;
	private double mpc;
	private double gnomadAF;
	private double VQSR;
	private double pextScore;
	private int altNumber;
	private AnnotationParser caddRecovery;
	private static DecimalFormat df;
		
	public PrintableVariant(VariantContext variantContext, GenomeVersion genomeVers, int altNumber, AnnotationParser caddRecovery) throws Exception {
		
		super(variantContext);
		parsedChr = checkChrAnnotation(genomeVers);
		this.altNumber = altNumber;
		this.refString = this.getReference().getBaseString();
		this.altString = this.getAlternateAllele(altNumber).getBaseString();
		this.caddRecovery = caddRecovery;
		df = new DecimalFormat("#");
			
	}
	
	// Is called after annotating VEP and VQSR
	public void fixEqualLengthSNV() {
		
		String ref = this.refString;
		String alt = this.altString;
		
		if (ref.length() > 1) {
			if (ref.length() == alt.length()) {
				
				char[] refVec = ref.toCharArray();
				char[] altVec = alt.toCharArray();
				
				for (int x = 0; x < refVec.length; x++) {
					
					if (refVec[x] != altVec[x]) {
						
						ref = String.valueOf(refVec[x]);
						alt = String.valueOf(altVec[x]);
						break;
					}				
				}
			}
		}
		
		this.refString = ref;
		this.altString = alt;		
		
	}
	
	public boolean isPrintable() {
		
		//First term is SNVs, apply SNV VQSR; Second term is InDels, apply InDel VQSR			
		if ((this.getRefBaseString().length() == this.getAltBaseString().length() && VQSR >= -3.6037) || 
				(this.getRefBaseString().length() != this.getAltBaseString().length() && VQSR >= -5.7733)) { 
			
			return(true);
		} else {
			return(false);
		}
		
	}
	public String printableVariant(Entry<String, Integer> sampleEntry, double AC, double AN, int AP, double AF, String geneID) throws ScoreNotFoundException {
		
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
	private String printDoubleScore(double score, AnnotationType annoteType) throws ScoreNotFoundException {
		
		String printableString;
		
		if (Double.isNaN(score)) {
			
			if (annoteType != AnnotationType.CADD) {
				// I think allow NaN? Can adjust for each annotation type. 
				// MPC is going to be unrecoverable I think
				// PEXT should ALWAYS be present
				// GnomadAF can definitely be NaN
				// VQSR should ALWAYS be present
				printableString = Double.toString(score);
				
			} else {
				//Try and recover CADD buy running the CADD API I've built...
				try {
					printableString = Double.toString(caddRecovery.getAnnotation(this));
				} catch (IOException e) {
					e.printStackTrace();
					throw new ScoreNotFoundException(annoteType, this);			
				}
			}
		} else {
			printableString = Double.toString(score);
		}
		
		return(printableString);
		
	}

 	public String getParsedChr() {
		return parsedChr;
	}
	
	public String getRefBaseString() {
		return(refString);	
	}
	public String getAltBaseString() {
		return(altString);
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
	// For setting different annotations
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
	
	private class ScoreNotFoundException extends Exception {

		private static final long serialVersionUID = 1L;
		
		public ScoreNotFoundException(AnnotationType annoteType, PrintableVariant printableVar) {
			super(annoteType.toString() + " score for variant " + printableVar.getParsedChr() + ":" + printableVar.getStart() + " " + printableVar.getRefBaseString() + "-" + printableVar.getAltBaseString() +  " not found!\n");
			// Extra for tabix printable if I ever need...
			// tabix " + annotationFile.getAbsolutePath() + " \"" + chr + ":" + (pos - 1) + "-" + (pos + 1) + "\"\n"
		}
		
	}
	
}
