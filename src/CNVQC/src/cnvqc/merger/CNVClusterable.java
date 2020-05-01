package cnvqc.merger;

import org.apache.commons.math3.ml.clustering.Clusterable;

import cnvqc.utilities.CNV;

public class CNVClusterable implements Clusterable {

	private double point[];
	private CNV cnv;
	
	public CNVClusterable(double x, double y, CNV cnv) {
		
		point = new double[]{x,y};
		this.cnv = cnv;
		
	}
	
	public double[] getPoint() {
		return point;
	}
	public CNV getCNV() {
		return cnv;
	}

}
