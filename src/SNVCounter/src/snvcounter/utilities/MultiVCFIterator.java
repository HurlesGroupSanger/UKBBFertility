package snvcounter.utilities;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.tribble.readers.TabixReader;
import htsjdk.tribble.readers.TabixReader.Iterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import snvcounter.annotator.Gene;

public class MultiVCFIterator {
	
	private Gene gene;
	private CloseableIterator<VariantContext> currentItr;
	private VCFFileReader currentVCFReader;
	private int fileNum;
	private List<VCFFile> vcfFiles;
	
	public MultiVCFIterator(File vcfFile, boolean isVcfListFile, Gene gene) throws IOException, InterruptedException {
	
		super();
		
		this.gene = gene;
		vcfFiles = new ArrayList<VCFFile>();
		fileNum = 0;
				
		if (isVcfListFile == false) {
			vcfFiles.add(new VCFFile(vcfFile));
		} else {
			
			TabixReader vcfFileTabix = new TabixReader(vcfFile.getAbsolutePath(),vcfFile.getAbsolutePath() + ".tbi");
			Iterator annotationItr = vcfFileTabix.query(gene.getChr(), gene.getStart() - 10, gene.getEnd() + 10);
			String coordinate;
			String data[];
			while ((coordinate = annotationItr.next()) != null) {
				data = coordinate.split("\t");
				vcfFiles.add(new VCFFile(new File(data[3]), Integer.valueOf(data[1])));
			}
		}		
		
		// Initialise the first iterator:
		initialiseIterator();		
		
	}

	public boolean hasNext() {
		
		if (currentItr.hasNext()) {
			return true;
		} else {
			if (fileNum < vcfFiles.size()) {
				return true;
			} else {
				return false;
			}
		}
		
	}

	public VariantContext next() throws IOException, InterruptedException {
		
		if (currentItr.hasNext()) {
			
			return currentItr.next();
			
		//Check if there is another file to iterate:
		} else if (fileNum < vcfFiles.size()) {
			
			this.close();
			initialiseIterator();
			return currentItr.next();
			
		} else {
			this.close();
			return null;
		}
		
	}

	public void close() {
		currentItr.close();
		currentVCFReader.close();
	}
	
	private void initialiseIterator() throws IOException, InterruptedException {
		
		VCFFile currentVCFPack = vcfFiles.get(fileNum);
		File currentVCF = currentVCFPack.getPath();
		File currentVCFIndex = EvaluateIndex.returnIndex(currentVCF);
		
		System.out.println(currentVCF.getAbsolutePath() + "\n" + currentVCFIndex.getAbsolutePath());
		
		if (currentVCFPack.hasLength()) {
			
			currentVCFReader = new VCFFileReader(currentVCF, currentVCFIndex);
			
			if (gene.getStart() > currentVCFPack.getStartCoord()) {
				currentItr = currentVCFReader.query(gene.getChr(), (gene.getStart() - 10), (gene.getEnd()+10));
			} else {
				currentItr = currentVCFReader.iterator();
			}
		} else {
			
			// Doing a new thing where I generate a tmp slice so that I don't overburden htsjdk (because that seems to be a problem...?)
			File tmpVCF = makeTempVCF(currentVCF);
			File tmpVCFIndex = EvaluateIndex.returnIndex(tmpVCF);
			
			currentVCFReader = new VCFFileReader(tmpVCF, tmpVCFIndex);
			
			currentItr = currentVCFReader.query(gene.getChr(), (gene.getStart() - 10), (gene.getEnd()+10));
		}
		
		fileNum++;
		
	}
	
	private File makeTempVCF(File currentVCF) throws IOException, InterruptedException {
		
		File tmpVCF = new File("/tmp/snvtmp." + gene.getChr() + "_" + (gene.getStart() - 10) + "_" + (gene.getEnd() + 10) + ".vcf.gz");
//		tmpVCF.deleteOnExit();
		String locus = gene.getChr() + ":" + (gene.getStart() - 10) + "-" + (gene.getEnd() + 10);
		
		String cmd;
		cmd = "bcftools view -S /lustre/scratch115/teams/hurles/users/eg15/INTERVAL/snv_cognition/ukbb_ids.sanger.txt -O z -o " + tmpVCF.getAbsolutePath() + " -r " + "" + locus + " " + currentVCF.getAbsolutePath();
		runJob(cmd);
		cmd = "bcftools index -t " + tmpVCF.getAbsolutePath();
		runJob(cmd);
		
		return(tmpVCF);
		
	}
	
	private void runJob(String cmd) throws IOException, InterruptedException {
		
		Process pr = Runtime.getRuntime().exec(cmd);
		BufferedReader br = new BufferedReader(new InputStreamReader(pr.getErrorStream()));
		int ext = pr.waitFor();
		if (ext != 0) {
			String line;
			System.err.println("Command:");
			System.err.println("\'" + cmd + "\'");
			System.err.println("failed to run bcftools properly... Exiting!");
			System.err.println("––––––––––––– BCFTOOLS ERROR –-–––––––––––");
			while ((line = br.readLine()) != null) {
				System.err.println(line);
			}
			System.err.println("––––––––––––– BCFTOOLS ERROR –-–––––––––––");
//			System.exit(1);
		}
		
	}
	
	private class VCFFile {
		
		private File path;
		private int startCoord;
		private boolean hasLength;
		
		private VCFFile(File path) {
			this.path = path;
			this.hasLength = false;
		}
		private VCFFile(File path, int startCoord) {
			this.path = path;
			this.startCoord = startCoord;
			this.hasLength = true;
		}

		private File getPath() {
			return path;
		}
		private int getStartCoord() {
			return startCoord;
		}
		private boolean hasLength() {
			return hasLength;
		}
		
	}
	
	
	
}
