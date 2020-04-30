package snvcounter.annotator;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import snvcounter.utilities.SNVCounterOptions.GenomeVersion;

public class Gene {

	private String ensID;
	private String chr;
	private int start;
	private int end;
	private int hgnc;
	private String geneName;
	
	public Gene (File geneList, int jobNumber, GenomeVersion genomeVers) throws IOException, GeneNotFoundException {
	
		boolean foundGene = getGene(geneList, jobNumber, genomeVers);
		if (foundGene == false) {
			throw new GeneNotFoundException();
		}		
	}

	public String getEnsID() {
		return ensID;
	}

	public String getChr() {
		return chr;
	}

	public int getStart() {
		return start;
	}

	public int getEnd() {
		return end;
	}

	public int getHgnc() {
		return hgnc;
	}

	public String getGeneName() {
		return geneName;
	}
	
	private boolean getGene(File geneList, int jobNumber, GenomeVersion genomeVers) throws IOException {
		
		BufferedReader geneReader = new BufferedReader(new FileReader(geneList));
		
		String line;
		String data[];
		int lineNum = 1;
		boolean foundGene = false;
		
		while ((line = geneReader.readLine()) != null) {
			
			data = line.split("\t");
			if (lineNum == jobNumber) {
				this.ensID = data[0];
				if (genomeVers == GenomeVersion.HG38) {
					this.chr = "chr" + data[1];
				} else {
					this.chr = data[1];
				}
				this.start = Integer.parseInt(data[2]);
				this.end = Integer.parseInt(data[3]);
				// Some genes do not have HGNC values for some reason...
				try {
					this.hgnc = Integer.parseInt(data[4]);
				} catch (NumberFormatException e) {
					this.hgnc = -1;
				}
				this.geneName = data[5];
				foundGene = true;
				break;
			}
			lineNum++;
			
		}
		
		geneReader.close();
		return foundGene;
		
	}
	
 	public class GeneNotFoundException extends Exception {

		/**
		 * 
		 */
		private static final long serialVersionUID = 1L;
		
		public GeneNotFoundException() {
			super("Job number extends beyond the end of the gene list file!");
		}
		
	}
	
}
