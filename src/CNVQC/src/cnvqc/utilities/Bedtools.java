package cnvqc.utilities;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.exec.CommandLine;
import org.apache.commons.exec.DefaultExecutor;
import org.apache.commons.exec.ExecuteException;
import org.apache.commons.exec.PumpStreamHandler;

public class Bedtools {

	public Bedtools() {
		super();
	}
	
	public static Map<String, CNV> printBED(List<CNV> cnvs, File tmpDir) throws IOException {
		
		BufferedWriter bedWriter = new BufferedWriter(new FileWriter(new File(tmpDir.getAbsolutePath() + "/merge_cnv.in.bed")));
		Map<String, CNV> mergedCNVs = new HashMap<String, CNV>();
		
		for (CNV cnv : cnvs) {
			bedWriter.write(cnv.getChr() + "\t" + cnv.getStart() + "\t" + cnv.getEnd() + "\t" + cnv.getSampleInformation().getSplitFile().getName() + "_" + cnv.getChr() + "_" + cnv.getStart() + "\n");
			mergedCNVs.put(cnv.getSampleInformation().getSplitFile().getName() + "_" + cnv.getChr() + "_" + cnv.getStart(), cnv);
		}
		
		bedWriter.flush();
		bedWriter.close();
		
		runBedTools("bedtools sort -i " + tmpDir.getAbsolutePath() + "/merge_cnv.in.bed", new File(tmpDir.getAbsolutePath() + "/merge_cnv.sorted.bed"));
		return mergedCNVs;
		
	}
	
	public static void runBedTools(String command, File outFile) throws ExecuteException, IOException {
		
		DefaultExecutor executor = new DefaultExecutor();
		FileOutputStream stdout = new FileOutputStream(outFile);
	    PumpStreamHandler psh = new PumpStreamHandler(stdout);
	    executor.setStreamHandler(psh);
	    CommandLine cmdLine = CommandLine.parse(command);
	    executor.execute(cmdLine);
		
	}
	
	
}
