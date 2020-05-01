package cnvqc;

import java.io.IOException;

import cnvqc.annotator.CNVSampleAnnotator;
import cnvqc.merger.CNVMerger;
import cnvqc.utilities.pathogenicannotation.PathogenicAnnotator.OverlapError;

public class CNVQCImplement {

	public static void main(String args[]) throws IOException, NumberFormatException, OverlapError {
		
		if (args.length == 0) {
			
			printHelp();
			
		} else {
			CNVRuntime runtime = CNVRuntime.valueOf(args[0].toUpperCase());
			String inputArgs[] = new String[args.length -1 ];
			
			for (int x = 1; x < args.length; x++) {
				inputArgs[x-1] = args[x];
			}
			
			if (runtime.equals(CNVRuntime.ANNOTATE)) {
				new CNVSampleAnnotator(inputArgs);
			}
			else if (runtime.equals(CNVRuntime.MERGE)) {
				new CNVMerger(inputArgs);
			}
			else {
				printHelp();
			}
		}
	}
	
	private static void printHelp() {
		System.err.println();
		System.err.println("Please use either \'Merge\', \'Polish\', or \'View\' as a Runtime (the first Argument)");
		System.err.println("Merge - Merge CNVs processed by Annotate and filtered by <>");
		System.err.println("Annotate - Attach sample and summary statistics to raw CNVs");
		System.err.println();
	}
	public enum CNVRuntime {
		ANNOTATE,MERGE;
	}
	
}
