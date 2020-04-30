package snvcounter.annotator;

public class Individual {

	private int syn;
	private int lof;
	private int miss;
	
	public Individual() {
		
		this.syn = 0;
		this.lof = 0;
		this.miss = 0;
		
	}

	public int getSynonymous() {
		return syn;
	}
	public int getLoF() {
		return lof;
	}
	public int getMissense() {
		return miss;
	}
		
	public void incrementSynonymous(int toInc) {
		syn += toInc;
	}
	public void incrementLoF(int toInc) {
		lof += toInc;
	}
	public void incrementMissense(int toInc) {
		miss += toInc;
	}
	
}
