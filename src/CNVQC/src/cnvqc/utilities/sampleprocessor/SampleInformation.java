package cnvqc.utilities.sampleprocessor;

import java.io.File;

public class SampleInformation {
	
	private File splitFile;
	private String CELName;
	private String wesID;
	private boolean hasWES;
	private Gender gender;
	private Double callRate;
	
	public SampleInformation(String CELName, File splitFile, String wesID, Gender gender, double callRate) {
		
		this.CELName = CELName;
		this.splitFile = splitFile;
		this.wesID = wesID;
		this.gender = gender;
		this.callRate = callRate;
		
		if (wesID == null) {
			hasWES = false;
		} else {
			hasWES = true;
		}
		
	}

	public File getSplitFile() {
		return splitFile;
	}
	public String getCELName() {
		return CELName;
	}
	public Gender getGender() {
		return gender;
	}
	public Double getCallRate() {
		return callRate;
	}
	public String getWESID() {
		return wesID;
	}
	public boolean hasWES() {
		return hasWES;
	}
	
	private double lrr_sd;
	private double lrr_mean;
	private double lrr_median;
	private double baf_sd;
	private double baf_mean;
	private double baf_median;
	private double baf_drift;
	private double wf;
	private int numCNV;
	private boolean isIndivFiltered;
	
	public void addFilterInformation(double wf, int numCNV) {
		this.wf = wf;
		this.numCNV = numCNV;
		isIndivFiltered = checkIndivFilter(numCNV, wf);
	}

	public double getLrr_sd() {
		return lrr_sd;
	}
	public void setLrr_sd(double lrr_sd) {
		this.lrr_sd = lrr_sd;
	}
	public double getLrr_mean() {
		return lrr_mean;
	}
	public void setLrr_mean(double lrr_mean) {
		this.lrr_mean = lrr_mean;
	}
	public double getLrr_median() {
		return lrr_median;
	}
	public void setLrr_median(double lrr_median) {
		this.lrr_median = lrr_median;
	}
	public double getBaf_sd() {
		return baf_sd;
	}
	public void setBaf_sd(double baf_sd) {
		this.baf_sd = baf_sd;
	}
	public double getBaf_mean() {
		return baf_mean;
	}
	public void setBaf_mean(double baf_mean) {
		this.baf_mean = baf_mean;
	}
	public double getBaf_median() {
		return baf_median;
	}
	public void setBaf_median(double baf_median) {
		this.baf_median = baf_median;
	}
	public double getBaf_drift() {
		return baf_drift;
	}
	public void setBaf_drift(double baf_drift) {
		this.baf_drift = baf_drift;
	}

	public double getWf() {
		return wf;
	}
	public int getNumCNV() {
		return numCNV;
	}

	public boolean isIndivFiltered() {
		return isIndivFiltered;
	}
	private boolean checkIndivFilter (int numCNV, double wf) {
		
		if (numCNV <= 30 && wf > -0.03 && wf < 0.03) {
			return false;
		} else {
			return true;
		}
		
	}
	
	public enum Gender {
		MALE,
		FEMALE,
		UNKNOWN;		
	}


	
}


