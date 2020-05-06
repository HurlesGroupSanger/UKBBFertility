package cnvqc.utilities.sampleprocessor;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import cnvqc.utilities.sampleprocessor.SampleInformation.Gender;

public class GenericSampleLoader {

	private Map<String, SampleInformation> sampleInformation;
	private Set<String> samples;
	
	public GenericSampleLoader(File sampleIDs, File qcsumFile) throws IOException {
		
		samples = new HashSet<String>();
		
		Map<File, String[]> qcInfo;
		if (qcsumFile != null) {
			qcInfo = buildQCInformation(qcsumFile);
		} else {
			qcInfo = null;
		}
		
		BufferedReader sampleReader = new BufferedReader(new FileReader(sampleIDs));
		
		sampleInformation = buildSampleInformation(sampleReader, qcInfo);
		
		sampleReader.close();
		
		System.err.println("Total samples parsed for sample information: " + sampleInformation.size());
		
	}
	
	public Map<String, SampleInformation> getSampleInformation() {
		return sampleInformation;
	}
	public Set<String> getSamples() {
		return samples;
	}
	
	private Map<File, String[]> buildQCInformation(File qcsumFile) throws IOException {
		
		Map<File, String[]> qcData = new HashMap<File, String[]>();
		
		BufferedReader qcReader = new BufferedReader(new FileReader(qcsumFile));
			
		String line;
		String data[];
			
		while ((line = qcReader.readLine()) != null) {
			data = line.split("\t");
			if (data[0].equals("File")) {
				continue;
			} else {
				File splitFile = new File(data[0]);
				qcData.put(splitFile, Arrays.copyOfRange(data, 1, data.length));
			}
		}
		
		qcReader.close();
		
		return qcData;
		
	}
	private Map<String, SampleInformation> buildSampleInformation(BufferedReader sampleReader, Map<File, String[]> qcInfo) throws IOException {
		
		sampleInformation = new HashMap<String, SampleInformation>();
		
		String line;
		
		while ((line = sampleReader.readLine()) != null) {
			
			String data[] = line.split("\t");
			
			String CELID = data[0];
			File splitFile = new File(data[1]);
			String wesID = data[2];
			Gender gender = Gender.valueOf(data[3]);
			Double quality = Double.parseDouble(data[4]);
			
			SampleInformation currentSamp = new SampleInformation(CELID,splitFile,wesID,gender,quality);
			
			// Add individual level QC metrics if needed:
			if (qcInfo != null) {
				
				String rawQCStats[] = qcInfo.get(splitFile);
				
				double wf = Double.parseDouble(rawQCStats[7]);
				int numCNV = Integer.parseInt(rawQCStats[8]);
				
				currentSamp.addFilterInformation(wf, numCNV);
				
				currentSamp.setLrr_mean(Double.parseDouble(rawQCStats[0]));
				currentSamp.setLrr_median(Double.parseDouble(rawQCStats[1]));
				currentSamp.setLrr_sd(Double.parseDouble(rawQCStats[2]));
				
				currentSamp.setBaf_mean(Double.parseDouble(rawQCStats[3]));
				currentSamp.setBaf_median(Double.parseDouble(rawQCStats[4]));
				currentSamp.setBaf_sd(Double.parseDouble(rawQCStats[5]));
				
				currentSamp.setBaf_drift(Double.parseDouble(rawQCStats[6]));
				
			}
			
			sampleInformation.put(splitFile.getName(), currentSamp);
			
			// Add all samples that pass QC into sample list
			if (!samples.contains(splitFile.getName())) {
				samples.add(splitFile.getName());
			}
					
		}
		
		return sampleInformation;
		
	}
	
}
