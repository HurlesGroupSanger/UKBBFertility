package cnvqc.merger;

import java.io.BufferedWriter;
import java.io.Closeable;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DateFormat;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.GregorianCalendar;
import java.util.List;
import java.util.Map;
import java.util.Set;

import cnvqc.merger.CNVMergerMethods.IndividualRecord;
import cnvqc.utilities.Combine;
import cnvqc.utilities.pathogenicannotation.PathogenicAnnotator.OverlapError;
import htsjdk.samtools.reference.FastaSequenceIndex;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class VCFEngine implements Closeable {

	private BufferedWriter vcfWriter;
	private IndexedFastaSequenceFile fastaRef;
	private Set<String> samples;
	private Set<String> chrs;
	private DecimalFormat df;
	
	public VCFEngine(File output, File fastaRef, File fastaIndex, Set<String> samples, Set<String> chrs) throws IOException {
		
		vcfWriter = new BufferedWriter(new FileWriter(new File(output.getAbsolutePath() + ".vcf")));
		this.fastaRef = new IndexedFastaSequenceFile(fastaRef, new FastaSequenceIndex(fastaIndex));
		this.samples = samples;
		this.chrs = chrs;
		df = new DecimalFormat("#.##");
		writeHeader();
		
	}
	
	@Override
	public void close() throws IOException {
		vcfWriter.close();
	}
	
	public void writeHeader() throws IOException {
		
		Date currentDate = GregorianCalendar.getInstance().getTime();
		DateFormat df = DateFormat.getDateTimeInstance();
		vcfWriter.write("##fileformat=VCFv4.2\n");
		vcfWriter.write("##fileDate=" + df.format(currentDate) + "\n");
		vcfWriter.write("##source=CNVPolisherv1.0\n");
		for (String chr : chrs) {
			int len = fastaRef.getSequence(chr).length();
			vcfWriter.write("##contig=<ID=" + chr + ",length=" + len + ">\n");
		}
		vcfWriter.write("##ALT=<ID=CN0,Description=\"Deletion\">\n");
		vcfWriter.write("##ALT=<ID=DUP,Description=\"Duplication\">\n");
		vcfWriter.write("##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles; If unknown, will be -1\">\n");
		vcfWriter.write("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">\n");
		vcfWriter.write("##INFO=<ID=AC,Number=1,Type=Integer,Description=\"Allele count for this record\">\n");
		vcfWriter.write("##INFO=<ID=DELAC,Number=1,Type=Integer,Description=\"DEL Allele count for this record\">\n");
		vcfWriter.write("##INFO=<ID=DUPAC,Number=1,Type=Integer,Description=\"DUP Allele count for this record\">\n");
		vcfWriter.write("##INFO=<ID=GENE,Number=.,Type=String,Description=\"Gene annotation for this variant\">\n");
		vcfWriter.write("##INFO=<ID=PATH,Number=.,Type=String,Description=\"Pathogenic annotation for this variant\">\n");
		vcfWriter.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
		vcfWriter.write("##FORMAT=<ID=WC,Number=1,Type=Float,Description=\"WES Overlap Confidence by Random Forest\">\n");
		vcfWriter.write("##FORMAT=<ID=DS,Number=1,Type=Integer,Description=\"Genotype Dosage\">\n");
		vcfWriter.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + Combine.combineSet(samples, "\t") + "\n");
		vcfWriter.flush();
		
	}
	
	public void addRecord(String chr, int start, int end, Map<String, IndividualRecord> cnvs, String pathString, boolean isMultiallelic) throws IOException, OverlapError {
		
		List<String> genotypeList = new ArrayList<String>();
				
		int dupAC = 0;
		int delAC = 0;
		
		for (String sample : samples) {

			String genotypeString;
			if (cnvs.containsKey(sample)) {
				
				String genotype;
				String confidence = df.format(cnvs.get(sample).getConfidence());
				double dosage;
				switch (cnvs.get(sample).getCopyNumber()) {
					case 0:
						genotype = "1/1";
						dosage = 0;
						delAC+=2;
						break;
					case 1:
						genotype = "0/1";
						dosage = 1;
						delAC+=1;
						break;
					case 3:
						if (isMultiallelic) {
							genotype = "0/2";
						} else {
							genotype = "0/1";
						}
						dosage = 3;
						dupAC+=1;
						break;
					case 4:
						if (isMultiallelic) {
							genotype = "2/2";
						} else {
							genotype = "1/1";
						}
						dosage = 4;
						dupAC+=2;
						break;
					default:
						genotype = "./.";
						dosage = Double.NaN;
						break;
						
				}
				genotypeString = genotype + ":" + confidence + ":" + df.format(dosage);
			} else {
				genotypeString = "0/0:.:2";
			}
		
			genotypeList.add(genotypeString);
			
		}
		
		writeRecord(genotypeList, chr, start, end, delAC, dupAC, pathString);
				
	}
	
	private void writeRecord(List<String> genotypes, String chr, int start, int stop, int delAC, int dupAC, String pathString) throws IOException, OverlapError {

		String ref = new String(fastaRef.getSubsequenceAt(chr, (long) start, (long) start).getBases());
		
		int len = stop - start;
		
		String altField;
		String acField = ";AC=" + (dupAC + delAC) + ";DELAC=" + delAC + ";DUPAC=" + dupAC;
		if (delAC == 0 && dupAC > 0) {
			altField = "<DUP>";
		} else if (delAC > 0 && dupAC == 0) {
			altField = "<CN0>";
		} else {
			altField = "<CN0>,<DUP>";
		}
		
		vcfWriter.write(chr + "\t" + start + "\t.\t" + ref + "\t" + altField + "\t.\t.\tEND=" + stop + ";SVLEN=" + len + acField + ";PATH=" + pathString + "\tGT:WC:DS\t" + Combine.combineList(genotypes, "\t") + "\n");
		vcfWriter.flush();
			
	}
	
}
