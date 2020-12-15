package snvcounter.utilities;

import java.io.File;
import java.io.IOException;
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
	private List<File> vcfFiles;
	
	public MultiVCFIterator(File vcfFile, boolean isVcfListFile, Gene gene) throws IOException {
	
		super();
		
		this.gene = gene;
		vcfFiles = new ArrayList<File>();
		fileNum = 0;
		
		if (isVcfListFile == false) {
			vcfFiles.add(vcfFile);
		} else {
			TabixReader vcfFileTabix = new TabixReader(vcfFile.getAbsolutePath(),vcfFile.getAbsolutePath() + ".tbi");
			Iterator annotationItr = vcfFileTabix.query(gene.getChr(), gene.getStart() - 10, gene.getEnd() + 10);
			
			String coordinate;
			String data[];
			while ((coordinate = annotationItr.next()) != null) {
				data = coordinate.split("\t");
				vcfFiles.add(new File(data[3]));
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

	public VariantContext next() {
		
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
	
	private void initialiseIterator() {
		
		
		File currentVCF = vcfFiles.get(fileNum);
		File currentVCFIndex = EvaluateIndex.returnIndex(currentVCF);
		currentVCFReader = new VCFFileReader(currentVCF, currentVCFIndex);	
		currentItr = currentVCFReader.query(gene.getChr(), gene.getStart(), gene.getEnd());
		
		fileNum++;
		
	}
	
	
	
	
}
