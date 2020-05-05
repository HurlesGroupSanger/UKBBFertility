package snvcounter.annotator;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math3.stat.inference.AlternativeHypothesis;
import org.apache.commons.math3.stat.inference.BinomialTest;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;

public class Genotyper {

	private PrintableVariant printableVar;
	private Set<String> sampleIDs;
	private double AN;
	private double AC;
	private int AP;
	
	private Map<String, Integer> samplesWithVariant;
		
	// This class goes through all individual genotypes ONLY for samples in the 'sampleIDs' list
	public Genotyper (PrintableVariant printableVar, Set<String> sampleIDs) {
		
		this.printableVar = printableVar;
		this.sampleIDs = sampleIDs;
		AN = 0;
		AC = 0;
		AP = 0;
		samplesWithVariant = new HashMap<String, Integer>();
		genotype();
		
	}
		
	public double getAN() {
		return AN;
	}
	public double getAC() {
		return AC;
	}
	public int getAP() {
		return AP;
	}
	public double getAF() {
		return AC / AN;
	}
	public Map<String, Integer> getSamplesWithVariant() {
		return samplesWithVariant;
	}

	private void genotype() {
		
		GenotypesContext currentGenotypes = printableVar.getGenotypes(sampleIDs);
		Iterator<Genotype> genoItr = currentGenotypes.iterator();
		
		// Now iterate through genotypes if we found a relevant consequence:		
		while (genoItr.hasNext()) {
			Genotype geno = genoItr.next();
			if (geno.isNoCall()) {
				AP+=2;
			} else {
				int DP = geno.getDP();
				int GQ = geno.getGQ();
				int AD[] = geno.getAD();
				
				// This is necessary to be able to run UKBB data since it does not have annotated DP/GQ/AD values, 
				// just converts them to values that pass our defaults
				DP = DP == -1 ? 10 : DP;
				GQ = GQ == -1 ? 40 : GQ;
				int altDepth = AD == null ? 5 : AD[printableVar.getAltNumber() + 1];
						
				int alleleCount = geno.countAllele(printableVar.getAltAllele());
				double adPval = alleleCount == 1 ? calculateADpval(DP, altDepth) : 1.0;
				
				AP += 2;
				AN += 2;
				
				// Apply Genotype Level Filters
				if (DP >= 7 && adPval > 0.001 && GQ >= 20) {
					AC += alleleCount;
									
					if (alleleCount > 0) {
						samplesWithVariant.put(geno.getSampleName(), alleleCount);
					}
				} else {
					AN -= 2;
				}
			}
		}		
		
	}
	
	private static double calculateADpval(int DP, int AD) {
		
		double adPval;
		try {
			adPval = new BinomialTest().binomialTest(DP, AD, 0.5, AlternativeHypothesis.TWO_SIDED);
		} catch (Exception e) {
			adPval = 0.0000001;
		}
		return adPval;
	}
	
}
