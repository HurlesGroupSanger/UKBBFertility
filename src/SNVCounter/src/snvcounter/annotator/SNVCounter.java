package snvcounter.annotator;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.List;
import java.util.Map;
import java.util.Set;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import snvcounter.parsers.AnnotationParser;
import snvcounter.parsers.AnnotationParser.AnnotationType;
import snvcounter.parsers.VEPAnnotator;
import snvcounter.utilities.EvaluateIndex;
import snvcounter.utilities.SNVCounterOptions.GenomeVersion;
import snvcounter.utilities.SNVCounterOptions.OptionHolder;

public class SNVCounter {

	OptionHolder options;
	
	public SNVCounter(OptionHolder options) {
		
		this.options = options;
		
	}
	
	public void doWork() throws Exception {

		// First get the gene that we are interested in:
		// Just to note, gene list provided here is all genes with a pLI/sHET score
		File geneList = options.getGeneList();
		int jobNumber = options.getGeneNumber();
		GenomeVersion genomeVers = options.getGenomeVersion();
		
		Gene gene = new Gene(geneList, jobNumber, genomeVers);
		
		// Open the VCF file we care about
		File vcfFile = options.getVcfFile();
		File vcfFileIndex = EvaluateIndex.returnIndex(vcfFile);
		VCFFileReader vcfReader = new VCFFileReader(vcfFile, vcfFileIndex);
				
		// Sample IDs are here.
		Set<String> sampleIDs = options.getSampleIDs();
		
		double afFilter = options.getMafFilter();
		
		File outDir = options.getOutputDir();
		BufferedWriter outWriter = new BufferedWriter(new FileWriter(new File(outDir.getAbsolutePath() + "/" + gene.getEnsID() + ".rare_variants.tsv")));
		
		File rootDir = options.getRootDir();
		
		/* We are annotating our VCFs with (file location after):
		 * 
		 * 1. MPC = annotations/<genomeVers>/mpc/fordist_constraint_official_mpc_values_v2.GRCh38.tsv.gz
		 * 2. CADD = annotations/<genomeVers>/cadd/whole_genome_SNVs.tsv.gz
		 * 3. GNOMAD AF = annotations/<genomeVers>/gnomad/gnomad.tsv.gz
		 * 4. PEXT = annotations/<genomeVers>/pext/pext.tsv.gz
		 * 5. VQSR = located with the VCF file
		 * 6. VEP = located with the VCF file
		 * 7. Supp CADD = located with the VCF file 
		 * 
		 * Build those parsers here:
		 * Actual instructions on building the files is in my RMarkdown for this project.
		 */

		// This gets the prefix necessary for fetching VEP, VQSR, and CADD annotations
		String annotationPrefix = EvaluateIndex.returnAnnotationPrefix(vcfFile);
		
		// MPC
		File mpcDict = new File(rootDir.getAbsolutePath() + "/" + genomeVers.getStringRep() + "/mpc/mpc.tsv.gz");
		AnnotationParser mpc = new AnnotationParser(mpcDict, AnnotationType.MPC, gene);
		
		// CADD -- Also includes the class for running CADD to annotate InDels
		File caddDict = new File(rootDir.getAbsolutePath() + "/" + genomeVers.getStringRep() + "/cadd/whole_genome_SNVs.tsv.gz");
		AnnotationParser cadd = new AnnotationParser(caddDict, AnnotationType.CADD, gene);
		File supplementalCADD = new File(annotationPrefix + ".cadd/cadd.tsv.gz");
		AnnotationParser caddRecovery = new AnnotationParser(supplementalCADD, AnnotationType.CADD, gene);

		// Gnomad
		File gnomadDict = new File(rootDir.getAbsolutePath() + "/" + genomeVers.getStringRep() + "/gnomad/gnomad.tsv.gz");
		AnnotationParser gnomad = new AnnotationParser(gnomadDict, AnnotationType.GNOMAD, gene);
		
		// VQSR
		File vqsrDict = new File(annotationPrefix + ".vqsr/VQSR.tsv.gz");
		AnnotationParser vqsr = new AnnotationParser(vqsrDict, AnnotationType.VQSR, gene);
		
		// PEXT
		File pextDict = new File(rootDir.getAbsolutePath() + "/" + genomeVers.getStringRep() + "/pext/pext.tsv.gz");
		AnnotationParser pext = new AnnotationParser(pextDict, AnnotationType.PEXT, gene);
		
		// VEP
		File vepDict = new File(annotationPrefix + ".vep/vep.tsv.gz");
		VEPAnnotator vep = new VEPAnnotator(vepDict, AnnotationType.VEP, gene);		
		
		// Now go through the VCF file by querying and attach information to each variant!
		CloseableIterator<VariantContext> vcfIterator = vcfReader.query(gene.getChr(), gene.getStart(), gene.getEnd());
				
		while (vcfIterator.hasNext()) {
			
			VariantContext currentVar = vcfIterator.next();
			List<Allele> altAlleles = currentVar.getAlternateAlleles();
			
			// Go through each allele properly:
			for (int alleleNum = 0; alleleNum < altAlleles.size(); alleleNum++) {
				
				PrintableVariant printableVar = new PrintableVariant(currentVar, genomeVers, alleleNum, caddRecovery);

				if (!printableVar.getAltBaseString().equals("*")) { // Remove any star alleles from analysis...
	
					// Annotate VEP
					printableVar.setVepAnnotation(vep.getVEPAnnotation(printableVar));
								
					// Have to do VQSR and VEP first since they use the VCF annotation and we have to left correct SNVs below...
					// Annotate VQSR AF
					printableVar.setVQSR(vqsr.getAnnotation(printableVar));
											
					// Method left corrects SNVs that are multiallelic with an InDel as required by the following annotations:
					printableVar.fixEqualLengthSNV();
					
					// Annotate CADD score
					printableVar.setCadd(cadd.getAnnotation(printableVar));
					
					// Annotate MPC score
					printableVar.setMpc(mpc.getAnnotation(printableVar));
											
					// Annotate Gnomad AF
					printableVar.setGnomadAF(gnomad.getAnnotation(printableVar));
					
					// Annotate PEXT
					printableVar.setPextScore(pext.getAnnotation(printableVar));
										
					// Apply Very Lenient VQSR Filter Here
					if(printableVar.isPrintable()) {
												
						Genotyper genotyper = new Genotyper(printableVar, sampleIDs);
						
						if (genotyper.getAF() < afFilter) {
						
							for (Map.Entry<String, Integer> sampleEntry : genotyper.getSamplesWithVariant().entrySet()) {
								
								outWriter.write(printableVar.printableVariant(sampleEntry, genotyper.getAC(), genotyper.getAN(), genotyper.getAP(), genotyper.getAF(), gene.getEnsID()) + "\n");
								outWriter.flush();
							}
						}
					}
				}
			}
		}
		
		mpc.close();
		cadd.close();
		gnomad.close();
		vqsr.close();
		vep.close();
		pext.close();
		vcfIterator.close();
		vcfReader.close();
		outWriter.close();
		
	}
	
}
