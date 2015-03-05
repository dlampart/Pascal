package ch.unil.genescore.gene;

import java.util.LinkedHashMap;

import ch.unil.genescore.main.FileParser;
import ch.unil.genescore.main.Main;
import ch.unil.genescore.main.Settings;


/**
 * Gencode gene annotation.
 * Gene types are (attribute 'gene_type'): 
 *    antisense
 *    snRNA
 *    miRNA
 *    snoRNA
 *    polymorphic_pseudogene
 *    lincRNA
 *    protein_coding
 *    3prime_overlapping_ncrna
 *    misc_RNA
 *    rRNA
 *    sense_intronic
 *    processed_transcript
 *    pseudogene
 *    sense_overlapping
 */
public class GeneAnnotationGencode extends GeneAnnotation {

	
	// ============================================================================
	// PUBLIC METHODS
	
	/** Constructor */
	public GeneAnnotationGencode() {
		
		super(Settings.gencodeAnnotationFile_);
	}
	

	// ----------------------------------------------------------------------------

	/** Load gene coordinates for the given gene set (load all genes if the set is empty) */
	public LinkedHashMap<String, Gene> loadAnnotation() {
				
		genes_ = new LinkedHashMap<String, Gene>();
		// To print gene types
		//HashSet<String> geneType = new HashSet<String>();
		
		// Open the file
		FileParser parser = new FileParser(annotationFile_);
		//GeneIdMapping mapping = GeneIdMapping.getInstance();
		
		// Skip the first 5 lines (start with #)
		String[] nextLine = parser.readLine();
		while (nextLine[0].startsWith("#"))
			nextLine = parser.readLine();
				
		while (nextLine != null) {
			// Check number of columns
			if (nextLine.length != 9)
				parser.error("Expected 9 columns");
			
			// Check that this is a gene
			if (!nextLine[2].equals("gene"))
				parser.error("Third column expected to be 'gene'");

			// Chromosome
			String chr = nextLine[0];
			// Continue if this chromosome should not be loaded based on settings
			if (skipChromosome(chr)) {
				nextLine = parser.readLine();
				continue;
			}

			if (loadOnlyProteinCoding_) {
				// Check that it's a protein coding gene
				String gene_type = getGencodeKeyValue(nextLine[8], "gene_type");
				//geneType.add(gene_type);
				
				if (!gene_type.equalsIgnoreCase("protein_coding")) {
					nextLine = parser.readLine();
					continue;
				}
				
			}
			// Gene id
			String geneId = removeEnsemblVersion(getGencodeKeyValue(nextLine[8], "gene_id").toUpperCase());
			String geneName = getGencodeKeyValue(nextLine[8], "gene_name").toUpperCase();
			
			// If a gene set to be loaded was specified and this gene is NOT in this set, continue
			if (!checkGenesToBeLoaded(geneId, geneName)) {
				nextLine = parser.readLine();
				continue;				
			}
					
			// Check that the id is unique
			if (genes_.containsKey(geneId))
				parser.error("Duplicate gene id: " + geneId);
			
			// Create the gene
			Gene nextGene = new Gene(geneId, geneName);
			genes_.put(geneId, nextGene);
				
			// Position
			int start = Integer.parseInt(nextLine[3]);
			int end = Integer.parseInt(nextLine[4]);
			boolean posStrand = GenomicElement.isPosStrand(nextLine[6]);	
			nextGene.setPosition(chr, start, end, posStrand);

			// Read next line
			nextLine = parser.readLine();
		}		
		parser.close();		

		// Print gene types
		//for (String type : geneType)
		//Ngsea.println(type);
	
		return genes_;
	}
	
	
	// ----------------------------------------------------------------------------

	/** If this is an ensembl id, remove the version number */
	public static String removeEnsemblVersion(String id) {
		
		if (id.length() > 4 && id.substring(0, 4).equals("ENSG")) {
			int dot = id.lastIndexOf(".");
			if (dot != -1)
				id = id.substring(0, dot);
		}
		return id;
	}

	
	// ============================================================================
	// PRIVATE METHODS

	/** Get the value of the given key, throw exception if not found */
	private String getGencodeKeyValue(String keyValueList, String key) {
		
		int start = keyValueList.indexOf(key + " \"");
		if (start == -1)
			Main.error("Key not found: '" + key + "\"");
		
		start = start + key.length() + 2;
		int end = keyValueList.indexOf("\"", start);
		
		String value = keyValueList.substring(start, end);
		return value;
	}
	
}
