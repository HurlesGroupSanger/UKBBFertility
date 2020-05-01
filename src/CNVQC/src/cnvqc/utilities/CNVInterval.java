package cnvqc.utilities;

import java.util.ArrayList;
import java.util.List;

import htsjdk.samtools.util.Interval;

public class CNVInterval extends Interval {

	private List<CNV> CNVs;
	
	public CNVInterval(String chr, int start, int end) {
		super(chr, start, end);
		CNVs = new ArrayList<CNV>();
	}
	public CNVInterval(CNV cnv) {
		super(cnv.getChr(), cnv.getStart(), cnv.getEnd());
		CNVs = new ArrayList<CNV>();
		CNVs.add(cnv);
	}
	public CNVInterval(String chr, int start, int end, CNV cnv) {
		super(chr, start, end);
		CNVs = new ArrayList<CNV>();
		CNVs.add(cnv);
	}
	public CNVInterval(String chr, int start, int end, List<CNV> CNVs) {
		super(chr, start, end);
		this.CNVs = CNVs;
	}
	
	public List<CNV> getAttachedCNVs() {
		return CNVs;
	}
	
	public CNVInterval union(CNVInterval that) {
		
		CNVInterval unionInterval;
		
		if (that.getAttachedCNVs().size() != 0) {
			CNVs.addAll(that.getAttachedCNVs());
			unionInterval = new CNVInterval(this.getContig(),
					Math.min(this.getStart(), that.getStart()),
					Math.max(this.getEnd(), that.getEnd()),
					this.getAttachedCNVs());
		} else {
			unionInterval = new CNVInterval(this.getContig(),
					Math.min(this.getStart(), that.getStart()),
					Math.max(this.getEnd(), that.getEnd()));
		}
		return unionInterval;
		
	}
		
}
