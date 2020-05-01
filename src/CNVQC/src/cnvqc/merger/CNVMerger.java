package cnvqc.merger;

import java.io.File;
import java.io.IOException;
import java.util.List;

import cnvqc.utilities.CNV;
import cnvqc.utilities.cnvreaders.ProcessedCNVReader;
import cnvqc.utilities.pathogenicannotation.PathogenicAnnotator.OverlapError;
import cnvqc.utilities.sampleprocessor.GenericSampleLoader;

public class CNVMerger {
	
	private static File tmpDir;
	
	public CNVMerger(String args[]) throws IOException, OverlapError {
		
		CNVMergerOptions options = new CNVMergerOptions(args);
		File toMerge = options.getRawCNVs();
		tmpDir = options.getTmpDirectory();
		
		GenericSampleLoader sampleLoader = new GenericSampleLoader(options.getSampleFile(), null);
				
		//This will read raw CNVS and print all raw information necessary for filtering
		ProcessedCNVReader reader = new ProcessedCNVReader(toMerge);
		List<CNV> CNVs = reader.getCNVs(sampleLoader.getSampleInformation()); //TRUE flag for only WES samples
		reader.close();
		
		//Build constructor for merging CNVs (sets up VCF writer and does initial CNV merge)
		
		ReciprocalClustering clusterer = new ReciprocalClustering(tmpDir);
		
		CNVMergerMethods mergerMethods = new CNVMergerMethods(options.getOutput(),
				options.getFastaRef(),
				options.getFastaIndex(),
				sampleLoader.getSamples(),
				reader.getChrs(),
				tmpDir,
				CNVs,
				clusterer);
		
		mergerMethods.MergeCNVs();
		
		mergerMethods.close();	
		
	}
	
}
