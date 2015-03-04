package ch.unil.genescore.vegas;

import java.util.LinkedHashMap;

import ch.unil.genescore.gene.Gene;
import ch.unil.genescore.main.FileExport;
import ch.unil.genescore.main.Main;
import ch.unil.genescore.main.Settings;


public class GeneResultsSnpsOutOfBounds {

	private LinkedHashMap<Gene, Integer> zeroSnpsOrAboveSnpLimit_ = null;	
	private FileExport exporter_ = null;
	
	public GeneResultsSnpsOutOfBounds(){
		exporter_ = new FileExport();
		zeroSnpsOrAboveSnpLimit_ = new LinkedHashMap<Gene, Integer>(); 
	}

	
	/** Write the list of genes that have zero or beyond max num snps */
	public void writeResultsToFile(String additionalOutputFileSuffix) {

		if (getZeroSnpsOrAboveSnpLimit().size() == 0)
			return;
		
		Main.warning("Writing genes without SNPs or exceeding the maximum number of SNPs");
		setExporter(additionalOutputFileSuffix);
		
		String header = "gene_id\tsymbol\tSNPs";
		exporter_.println(header);
		for (Gene gene : getZeroSnpsOrAboveSnpLimit().keySet()) {
			String line = gene.id_;
			line += "\t" + gene.getSymbolOrNA();
			line += "\t" + getZeroSnpsOrAboveSnpLimit().get(gene);
			exporter_.println(line);
		}
		exporter_.close();
		Main.println("");
	}
	
	public LinkedHashMap<Gene, Integer> getZeroSnpsOrAboveSnpLimit() {
		return zeroSnpsOrAboveSnpLimit_;
	}

	
	public void setZeroSnpsOrAboveSnpLimit(LinkedHashMap<Gene, Integer> zeroSnpsOrAboveSnpLimit) {
		this.zeroSnpsOrAboveSnpLimit_ = zeroSnpsOrAboveSnpLimit;
	}

	public void setExporter(FileExport exp) {
		exporter_=exp;
	}
	
	public void setExporter(String additionalOutputFileSuffix) {
		String filename = Settings.outputDirectory_ + "/" + Settings.gwasName_ + additionalOutputFileSuffix + ".numSnpError" + Settings.chromFileExtension_ + ".txt";
		exporter_.setWriter(filename);		
	}

	public void addToMap(GeneWithItsSnps geneWithItsSnps) {
		getZeroSnpsOrAboveSnpLimit().put(geneWithItsSnps.getGene(),geneWithItsSnps.getNrOfSnps());		
	}
	
	
}
