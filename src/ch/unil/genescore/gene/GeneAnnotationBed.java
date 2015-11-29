/*******************************************************************************
 * Copyright (c) 2015 David Lamparter, Daniel Marbach
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *******************************************************************************/
package ch.unil.genescore.gene;

import java.io.File;
import java.util.HashSet;
import java.util.LinkedHashMap;

import ch.unil.genescore.main.Pascal;
import ch.unil.gpsutils.FileParser;

/**
 * bed-file annotation with genSymbol as ids
 */
public class GeneAnnotationBed extends GeneAnnotation {

	public GeneAnnotationBed(File annotationFile) {
		super(annotationFile);
		// TODO Auto-generated constructor stub
	}

	/** Constructor */
	public GeneAnnotationBed() {
		
		super(Pascal.set.bedAnnotationFile_);
		System.out.println("Load annotation: " + Pascal.set.bedAnnotationFile_);
	}

	@SuppressWarnings({ "unchecked", "rawtypes" })
	@Override
	public LinkedHashMap<String, Gene> loadAnnotation() {
		// TODO Auto-generated method stub
		genes_ = new LinkedHashMap<String, Gene>();
		HashSet<String> inconsistentEntries = new HashSet<String>();
		HashSet<String> inconsistentStrand = new HashSet<String>();
		System.out.println("Loading annotation: " + annotationFile_);
		FileParser parser = new FileParser(Pascal.log, annotationFile_);
		
		// Skip the header lines (start with #)
		String[] nextLine = parser.readLine();
		while (nextLine[0].startsWith("#"))
			nextLine = parser.readLine();

		while (nextLine != null) {
			// Check number of columns
			if (nextLine.length != 6)
				parser.error("Expected 6 columns");

			// Chromosome
			String chr = nextLine[0];
			// Continue if this chromosome should not be loaded based on settings
			if (skipChromosome(chr)) {
				nextLine = parser.readLine();
				continue;
			}
			
			// Gene id and name
			String geneId = nextLine[3];
			String geneName = nextLine[3];

			// If a gene set to be loaded was specified and this gene is NOT in this set, continue
			if (!checkGenesToBeLoaded(geneId, geneName)) {
				nextLine = parser.readLine();
				continue;				
			}

		
			// Strand, start and end
						boolean posStrand = GenomicElement.isPosStrand(nextLine[5]);
						int start = Integer.parseInt(nextLine[1]);
						int end = Integer.parseInt(nextLine[2]);

						Gene nextGene = genes_.get(geneId);
						// Create new gene
						if (nextGene == null) {
							nextGene = new Gene(geneId, geneName);
							nextGene.setPosition(chr, start, end, posStrand);
							genes_.put(geneId, nextGene);
						
						// Update existing gene
						} else {
							// Check that symbol, chr and position are consistent
							if (nextGene.chr_.equals(chr) && nextGene.symbol_.equals(geneName) &&
									(Math.abs(nextGene.start_ - start) < 1e6 || Math.abs(nextGene.end_ - end) < 1e6)) {
								
								// If only the strand is not consistent, we'll just give a warning 
								if (nextGene.posStrand_ != posStrand) {
									// We believe the longer gene
									if (end - start > nextGene.end_ - nextGene.start_)
										nextGene.posStrand_ = posStrand;
									inconsistentStrand.add(geneId);
								}
								
								// Update start and end
								if (start < nextGene.start_)
									nextGene.start_ = start;
								if (end > nextGene.end_)
									nextGene.end_ = end;
								
							} else {
								inconsistentEntries.add(geneId);
								genes_.remove(geneId);
							}
						}
						nextLine = parser.readLine();
		}
		parser.close();		

		// Sort the genes
		Genome genome = new Genome(genes_.values());
		genes_ = (LinkedHashMap) genome.asLinkedHashMap();
		
		return genes_;
							
		
	}

}
