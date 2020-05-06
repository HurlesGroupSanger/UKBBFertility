package cnvqc.utilities.cnvreaders;

import java.io.BufferedReader;
import java.io.Closeable;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

import cnvqc.utilities.CNV;
import cnvqc.utilities.CNV.CopyType;
import cnvqc.utilities.sampleprocessor.SampleInformation;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.IntervalTree.Node;
import htsjdk.tribble.readers.TabixReader;
import htsjdk.tribble.readers.TabixReader.Iterator;

public class RawCNVReader implements Closeable {

	private Set<String> samples;
	private Set<String> chrs;
	private Map<String, SampleInformation> sampleInformation;
	private BufferedReader cnvReader;
	private TabixReader axiomTabix;
	private TabixReader baitTabix;
	private TabixReader xhmmReader;
	private TabixReader canoesReader;
	private TabixReader clammsReader; 
	
	public RawCNVReader(File CNVs, Map<String, SampleInformation> sampleInformation, File xhmmFile, File clammsFile, File canoesFile, File axiomProbes, File exomeProbes) throws IOException {
		
		cnvReader = new BufferedReader(new FileReader(CNVs));
		samples = new HashSet<String>();
		chrs = new HashSet<String>();
		this.sampleInformation = sampleInformation;
		
		// Load WES tabix files (if applicable)
		xhmmReader = getWESTabix(xhmmFile);
		clammsReader = getWESTabix(clammsFile);
		canoesReader = getWESTabix(canoesFile);
		
		// Get Probe information for arrays and WES CNVs
		axiomTabix = new TabixReader(axiomProbes.getAbsolutePath(),axiomProbes.getAbsolutePath() + ".tbi");
		baitTabix = new TabixReader(exomeProbes.getAbsolutePath(),exomeProbes.getAbsolutePath() + ".tbi");
				
	}
	
	public Set<String> getChrs() {
		return chrs;
	}
	public Set<String> getSamples() {
		return samples;
	}
	public void close() throws IOException {
		cnvReader.close();
		baitTabix.close();
		axiomTabix.close();
	}
	
	public List<CNV> getAllCNVs(int lineStart, int lineEnd) throws NumberFormatException, IOException {
		
		List<CNV> cnvs = new ArrayList<CNV>();
		String line;
		
		int totalCNVs = 0;
		int lineNum = 1;
		
		while ((line = cnvReader.readLine()) != null) {
			if (lineNum > lineEnd) {
				break;
			} else if (lineNum >= lineStart) {
				CNV currentCNV = parseCNVLine(line);
				if (currentCNV != null) {
					totalCNVs++;
					cnvs.add(currentCNV);
				}
			}
			lineNum++;
			
		}
		
		System.err.println("Total CNVs loaded: " + totalCNVs);
		return cnvs;
		
	}
	private CNV parseCNVLine(String line) throws NumberFormatException, IOException {
		
		//PennCNV is white-space delim... whereas LRR BAF info isn't
		String lrrbaf[] = line.split("\t");
		String data[] = lrrbaf[0].split("\\s+");
		Pattern filePattern = Pattern.compile("(split\\d+\\.a\\d{6}\\S*)");
		
		//New Code to parse CNV:
		//Only want individuals that passed AFFY QC
		File splitFile = new File(data[4]);

		if (sampleInformation.containsKey(splitFile.getName())) {				
		
			SampleInformation sampInfo = sampleInformation.get(splitFile.getName());
			
			String chr = null;
			int start = -1;
			int end = -1;
			Matcher locationParser = Pattern.compile("chr([\\dXY]{1,2}):(\\d+)\\-(\\d+)").matcher(data[0]);
			if (locationParser.matches()) {
				chr = locationParser.group(1);
				start = Integer.parseInt(locationParser.group(2));
				end = Integer.parseInt(locationParser.group(3));
			}
	
			int probeCount = (int)EqualSpliter(data[1]);
			int copyNumber = (int)EqualSpliter(data[3]);
			double conf = EqualSpliter(data[7]);
			//New code to parse CNV
			
			//Captures sample names
			Matcher fileMatcher = filePattern.matcher(splitFile.getName());
			
			if (fileMatcher.matches()) {
				String fileName = fileMatcher.group(1);
				if (!samples.contains(fileName)) {
					samples.add(fileName);
				}
			}
			if (!chrs.contains(chr)) {
				chrs.add(chr);
			}

			CNV currentCNV = new CNV(
					chr,
					start,
					end,
					copyNumber,
					probeCount,
					conf,
					sampInfo,
					getLRRBAF(sampInfo.getSplitFile(), chr, start, end)
					);
			
			// Checks for intersection to WES baits -- even if the CNV doesn't have WES data (for annotation purposes)
			IntervalTree<Boolean> wesBaits = getIntersectingBaits(chr, start, end);	
			currentCNV.setTotalIntersectingBaits(wesBaits.size());
			
			//Check if there are any intersecting CNVs for this individual if individual has WES data
			if (sampInfo.hasWES()) {
				
				currentCNV.setIntersectingWESXHMMCNVs(getWESCNVs(sampInfo.getWESID(), currentCNV, xhmmReader, wesBaits)); //This checks XHMM CNVs.
				currentCNV.setIntersectingWESCLAMMSCNVs(getWESCNVs(sampInfo.getWESID(),currentCNV, clammsReader, wesBaits)); //This checks CLAMMS CNVs.
				currentCNV.setIntersectingWESCANOESCNVs(getWESCNVs(sampInfo.getWESID(),currentCNV, canoesReader, wesBaits)); //This checks CANOES CNVs.
				
			}

			return currentCNV;
			
		} else {
			
			return null;
			
		}

	}
	private IntervalTree<Boolean> getIntersectingProbes(String chr, int start, int end) throws IOException {
		
		IntervalTree<Boolean> baitTree = new IntervalTree<Boolean>();
		
		Iterator itr = axiomTabix.query(chr, start, end);

		String line;
		
		while ((line = itr.next()) != null) {

			String data[] = line.split("\t");
			
			int baitStart = Integer.parseInt(data[1]);
			int baitEnd = Integer.parseInt(data[2]);
			
			baitTree.put(baitStart, baitEnd, false);
			
		}
				
		return baitTree;
		
	}	
	private IntervalTree<Boolean> getIntersectingBaits(String chr, int start, int end) throws IOException {
		
		IntervalTree<Boolean> baitTree = new IntervalTree<Boolean>();
		
		Iterator itr = baitTabix.query(chr, start, end);

		String line;
		
		while ((line = itr.next()) != null) {

			String data[] = line.split("\t");
			
			int baitStart = Integer.parseInt(data[1]);
			int baitEnd = Integer.parseInt(data[2]);
			
			baitTree.put(baitStart, baitEnd, false);
			
		}
				
		return baitTree;
		
	}
	
	private double getWESCNVs(String ID, CNV cnv, TabixReader tabix, IntervalTree<Boolean> wesBaits) throws NumberFormatException, IOException {
		
		String line;
		String data[];
		Iterator itr = tabix.query(cnv.getChr(), cnv.getStart(), cnv.getEnd());			
		
		double totalBaits = wesBaits.size();
		double baitsIntersected = 0;
				
		double totalProbes = 0;
		double probesIntersected = 0;
								
		while ((line = itr.next()) != null) {
			
			data = line.split("\t");
				
			int startWES = Integer.parseInt(data[1]);
			int endWES = Integer.parseInt(data[2]);
			double numProbes = Double.parseDouble(data[3]);
			CopyType ctWES = CopyType.valueOf(data[5]);
			String eganIDWES = data[6];
			
			CNV wesCNV = new CNV(cnv.getChr(), startWES, endWES, 0, (int) numProbes, 0.0, null, null);
			ArrayList<CNV> toInt = new ArrayList<CNV>();
			toInt.add(cnv);
			toInt.add(wesCNV);
			Interval wesCnvInt = new Interval(cnv.getChr(), startWES, endWES);
			Interval cnvInt = new Interval(cnv.getChr(), cnv.getStart(), cnv.getEnd());
			
			if (eganIDWES.equals(ID) && cnv.getCopyType().equals(ctWES) && wesCnvInt.intersects(cnvInt)) {
								
				// This checks for the WES CNV encompassing WES baits covered by the array CNV
				java.util.Iterator<Node<Boolean>> baitOverlap = wesBaits.overlappers(startWES, endWES);
								
				while (baitOverlap.hasNext()) {
					baitOverlap.next();
					baitsIntersected++;
				}
				
				// This checks for the Array CNV encompassing array probes covered by the WES CNV
				IntervalTree<Boolean> axiomProbes = getIntersectingProbes(cnv.getChr(), startWES, endWES);
				totalProbes += axiomProbes.size();		
				
				java.util.Iterator<Node<Boolean>> probeOverlap = axiomProbes.overlappers(cnv.getStart(), cnv.getEnd() + 1);
							
				while (probeOverlap.hasNext()) {
					probeOverlap.next();
					probesIntersected++;
				}
			}
		}
		
		if (totalProbes == 0) {
			return(baitsIntersected / totalBaits);
		} else {
			return((baitsIntersected / totalBaits) * (probesIntersected / totalProbes));
		}
				
	}
	
	private LRRandBAFInformation getLRRBAF(File splitFile, String chr, int start, int end) throws IOException {
		
		TabixReader lrrbafTabixReader = new TabixReader(splitFile.getAbsolutePath() + ".sorted.bed.gz", splitFile.getAbsolutePath() + ".sorted.bed.gz.tbi");
		
		int len = end - start;
		int qStart = (start - len) < 0 ? 0 : (start - len);
		int qEnd = end + len;
	
		Iterator itr = lrrbafTabixReader.query(chr, qStart, qEnd);
		
		String line;
		String data[];
		int nLeft = 0;
		int nRight = 0;
				
		DescriptiveStatistics lrrStat = new DescriptiveStatistics();
		DescriptiveStatistics bafStat = new DescriptiveStatistics();
		
		while ((line = itr.next()) != null) {
			
			data = line.split("\t");
			int currPos = Integer.parseInt(data[1]);
			if (currPos < start) {
				nLeft++;
			} else if (currPos > (end + 1)) {
				nRight++;
			} else {
				double lrr = Double.parseDouble(data[4]);
				double baf = Double.parseDouble(data[5]);
		
				lrrStat.addValue(lrr);
				bafStat.addValue(Math.abs(0.5 - baf));
			}
			
		}
				
		LRRandBAFInformation lrrbaf = new LRRandBAFInformation(lrrStat, bafStat, nLeft, nRight);
		
		lrrbafTabixReader.close();
		return lrrbaf;
		
	}
	public class LRRandBAFInformation {
		
		private double lrrMean;
		private double lrrSD;
		private double lrrMedian;
		private double bafMean;
		private double bafSD;
		private double bafMedian;
		private int nLeft;
		private int nRight;
		private DecimalFormat df;
		
		private LRRandBAFInformation(DescriptiveStatistics lrrStat, DescriptiveStatistics bafStat, int nLeft, int nRight) {
			lrrMean = lrrStat.getMean();
			lrrSD = lrrStat.getStandardDeviation();
			lrrMedian = lrrStat.getPercentile(0.50);
			bafMean = bafStat.getMean();
			bafSD = bafStat.getStandardDeviation();
			bafMedian = bafStat.getPercentile(0.50);
			this.nLeft = nLeft;
			this.nRight = nRight;
			df = new DecimalFormat("##.#######");
		}
		
		private String format(double val) {
			if (Double.isNaN(val)) {
				return "NaN";
			} else {
				return df.format(val);
			}
		}
		public String returnPrintable() {
			return (format(bafMean) + "\t" + format(bafSD) + "\t" + format(bafMedian) + "\t" + format(lrrMean) + "\t" + format(lrrSD) + "\t" + format(lrrMedian) + "\t" + nLeft + "\t" + nRight);
		}
		
	}
	
	//Utilities to parse raw CNVs
	private double EqualSpliter (String toParse) {
		String parsed[] = toParse.split("\\=");
		String replaced = parsed[1].replaceAll(",","");
		return Double.parseDouble(replaced);
	}
	
	private TabixReader getWESTabix(File wesTabix) {
		
		TabixReader reader = null;
		
		if (wesTabix != null) {
			try {
				reader = new TabixReader(wesTabix.getAbsolutePath(),wesTabix.getAbsolutePath() + ".tbi");
			} catch (IOException e) {
				System.err.println("Could not load WES tabix file: " + wesTabix.getAbsolutePath());
				e.printStackTrace();
				System.exit(1);
			}
		}
		
		return reader;
		
	}	
}
