package snvcounter.utilities;

import java.io.File;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

// This class take a bam file and checks to see if a few different nomenclatures for indices are present
public class EvaluateIndex {
	
	public static File returnIndex(File vcfFile) {
		
		File vcfIndex;
		
		File vcfIndexTBI = new File(vcfFile.getAbsolutePath() + ".tbi");
		File vcfIndexCSI = new File(vcfFile.getAbsolutePath() + ".csi");
		if (vcfIndexTBI.exists()) {
			vcfIndex = vcfIndexTBI;
		} else if (vcfIndexCSI.exists()) {
			vcfIndex = vcfIndexCSI;
		} else {
			vcfIndex = null;
		}
				
		return vcfIndex;
				
	}
	
	public static String returnAnnotationPrefix(File vcfFile) {
				
		Matcher nameMatcher = Pattern.compile("(\\S+)\\.(vcf|bed)\\.gz").matcher(vcfFile.getAbsolutePath());
		String annotationPrefix;
		
		if (nameMatcher.matches()) {
			
			File vepAnnotation = new File(nameMatcher.group(1) + ".vep/vep.tsv.gz");
			File vqsrAnnotation = new File(nameMatcher.group(1) + ".vqsr/VQSR.tsv.gz");
			File CADDAnnotation = new File(nameMatcher.group(1) + ".cadd/cadd.tsv.gz");;
			if (!vepAnnotation.exists()) {
				annotationPrefix = null;
				System.err.println("VEP annotation is expected to be at " + nameMatcher.group(1) + ".vep/vep.tsv.gz but is not there! Exiting...");
				System.exit(1);
			} else if (!vqsrAnnotation.exists()) {
				annotationPrefix = null;
				System.err.println("VQSR annotation is expected to be at " + nameMatcher.group(1) + ".vep/VQSR.tsv.gz but is not there! Exiting...");
				System.exit(1);
			} else if (!CADDAnnotation.exists()) {
				annotationPrefix = null;
				System.err.println("Supplementary CADD annotation is expected to be at " + nameMatcher.group(1) + ".cadd/cadd.tsv.gz but is not there! Exiting...");
				System.exit(1);
			} else {
				annotationPrefix = nameMatcher.group(1);
			}
			
		} else {
			
			annotationPrefix = null;
			System.err.println("VCF|BED file is not named properly as *.vcf|bed.gz! Exiting...");
			System.exit(1);
			
		}
		
		return annotationPrefix;
		
	}
	
}
