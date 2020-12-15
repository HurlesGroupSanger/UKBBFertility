package snvcounter.utilities;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;

import snvcounter.annotator.Gene;
import snvcounter.annotator.Gene.GeneNotFoundException;

public class SNVCounterOptions {

	private OptionHolder holder;
	
	// Is just a holder/parser for command-line options
	public SNVCounterOptions(String args[]) throws IOException, GeneNotFoundException {
		
		Options options = buildOptions();
		holder = setOptions(args, options);
		
	}
	
	public OptionHolder getOptions() {
		
		return holder;
		
	}

	private Options buildOptions() {
		
		List<Option> reqOptions = new ArrayList<Option>();
		
		reqOptions.add(new Option("genelist", true, "List of genes to annotate for rare variants."));
		reqOptions.add(new Option("gene", true, "Gene number to annotate. Must be between 1 and length(genelist)."));
		reqOptions.add(new Option("genomeversion", true, "Genome version annotating from. Must be either hg38/hg19."));
		reqOptions.add(new Option("vcf", true, "VCF file to pull annotations from."));
		reqOptions.add(new Option("samples", true, "Name of samples within the VCF file to use."));
		reqOptions.add(new Option("maf", true, "Filter all variants with MAF > this value"));
		reqOptions.add(new Option("outdir", true, "Root directory for output. Will create output for gene in gene like outdir/<GENEID>.rare_variants.tsv"));
		reqOptions.add(new Option("rootdir", true, "Root dir where annotations are found."));
		
		Options options = new Options();
		
		for (Option opt : reqOptions) {
			opt.setRequired(true);
			options.addOption(opt);
		}
		
		return options;
		
	}
	
	private OptionHolder setOptions (String args[], Options options) throws IOException, GeneNotFoundException {
		
		CommandLineParser parser = new BasicParser();
		CommandLine cmd = null;
		
		try {
			 cmd = parser.parse(options, args);
		} catch (org.apache.commons.cli.ParseException e) {
			System.out.println();
			e.printStackTrace();
			System.exit(1);
		}
	
		OptionHolder holder = new OptionHolder(cmd.getOptionValue("gene"), 
				cmd.getOptionValue("genelist"),
				cmd.getOptionValue("genomeversion"),
				cmd.getOptionValue("vcf")
				);
		
		;
		holder.setSampleIDs(cmd.getOptionValue("samples"));
		holder.setMafFilter(cmd.getOptionValue("maf"));
		holder.setOutputDir(cmd.getOptionValue("outdir"));
		holder.setRootDir(cmd.getOptionValue("rootdir"));
		
		return holder;
		
	}
		
	public class OptionHolder {
		
		private File geneList;
		private int geneNumber;
		private GenomeVersion genomeVersion;
		private Gene gene;
		private File vcfFile;
		private boolean isVcfListFile;
		
		private MultiVCFIterator vcfIterator;
		
		private Set<String> sampleIDs;
		private double mafFilter;
		private File outputDir;
		private File rootDir;
		
		public OptionHolder(String geneList, String geneString, String genomeVersionString, String vcfFileString) throws IOException, GeneNotFoundException {
			
			setGeneList(geneList);
			setGeneNumber(geneString);
			setGenomeVersion(genomeVersionString);
			setVcfFile(vcfFileString);
			setGene();
			buildVCFIterator();
			
		}

		private void setGeneList(String geneList) {
			
			File geneFile = new File(geneList);
			if (geneFile.exists()) {
				this.geneList = geneFile;
			} else {
				System.out.println("Gene file given by -genelist does not exist... Exiting!");
				System.exit(1);
			}
			
		}
		private void setGeneNumber(String geneString) {
			
			try {
				this.geneNumber = Integer.parseInt(geneString);
			} catch (NumberFormatException e) {
				System.out.println("Value provided for -gene is not an integer.");
				System.exit(1);
			}
		}

		public GenomeVersion getGenomeVersion() {
			return genomeVersion;
		}
		private void setGenomeVersion(String genomeVersionString) {
			GenomeVersion genomeVersion = searchEnum(GenomeVersion.class, genomeVersionString);
			if (genomeVersion != null) {
				this.genomeVersion = genomeVersion;
			} else {
				System.out.println("Incorrect genome version format for -genomeversion. Must be either hg19/hg38... Exiting!");
				System.exit(1);
			}
			
		}
		
		public Gene getGene() {
			return gene;
		}
		private void setGene() throws IOException, GeneNotFoundException {
			gene = new Gene(geneList, geneNumber, genomeVersion);
		}
		
		public File getVcfFile() throws IOException {
			
			return vcfFile;
			
		}	
		private void setVcfFile(String vcfFileString) {
			
			File vcfFile = new File(vcfFileString);
			if (vcfFile.exists()) {
				if (vcfFile.getName().contains("vcf.gz")) {
					isVcfListFile = false;
				} else {
					isVcfListFile = true;
				}
				this.vcfFile = vcfFile;
			} else {
				System.out.println("VCF file given by -vcf does not exist... Exiting!");
				System.exit(1);
			}
			
		}

		public Set<String> getSampleIDs() {
			return sampleIDs;
		}
		public void setSampleIDs(String sampleIDFileString) throws IOException {
			
			File sampleIDFile = new File(sampleIDFileString);
			if (sampleIDFile.exists()) {
				this.sampleIDs = parseSampleIDs(sampleIDFile);
			} else {
				System.out.println("Sample ID File file given by -sammples does not exist... Exiting!");
				System.exit(1);
			}
			
		}

		public double getMafFilter() {
			return mafFilter;
		}
		public void setMafFilter(String mafFilterString) {
			
			try {
				this.mafFilter = Double.parseDouble(mafFilterString);
			} catch (NumberFormatException e) {
				System.out.println("Value provided for -maf is not a decimal number.");
				System.exit(1);
			}
			
		}

		public File getOutputDir() {
			return outputDir;
		}
		public void setOutputDir(String outputDirString) {
			
			File outputDir = new File(outputDirString);
			if (outputDir.exists() && outputDir.isDirectory()) {
				this.outputDir = outputDir;
			} else {
				System.out.println("Directory file given by -outputdir does not exist or is not a directory... Exiting!");
				System.exit(1);
			}
			
		}

		public File getRootDir() {
			return rootDir;
		}
		public void setRootDir(String rootDirString) {
			File rootDir = new File(rootDirString);
			if (rootDir.exists() && rootDir.isDirectory()) {
				this.rootDir = rootDir;
			} else {
				System.out.println("Directory file given by -rootdir does not exist or is not a directory... Exiting!");
				System.exit(1);
			}
			this.rootDir = rootDir;
		}
			
		public void buildVCFIterator() throws IOException {
			vcfIterator = new MultiVCFIterator(vcfFile, isVcfListFile, gene);
		}
		public MultiVCFIterator getVCFIterator() {
			return vcfIterator;
		}
		
	}
	
	private static Set<String> parseSampleIDs (File sampleNameFile) throws IOException {
		
		BufferedReader sampleNameFileReader = new BufferedReader(new FileReader(sampleNameFile));
		
		String line;
		
		Set<String> sampleIDs = new HashSet<String>();
		
		while ((line = sampleNameFileReader.readLine()) != null) {
			
			sampleIDs.add(line);
			
		}
		
		sampleNameFileReader.close();
		return sampleIDs;
				
	}
	
	public enum GenomeVersion {
		
		HG19("hg19"),HG38("hg38");
		
		private String stringRep;
		
		private GenomeVersion(String stringRep) {
			this.stringRep = stringRep;
		}
			
		public String getStringRep() {
			return stringRep;
		}
		
	}
	public static <T extends Enum<?>> T searchEnum(Class<T> enumeration, String search) {
		for (T each : enumeration.getEnumConstants()) {
			if (each.name().compareToIgnoreCase(search) == 0) {
				return each;
			}
	    }
		return null;
	}
	
}
