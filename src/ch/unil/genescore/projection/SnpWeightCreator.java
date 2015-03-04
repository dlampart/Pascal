package ch.unil.genescore.projection;

import java.util.LinkedHashMap;
import java.util.TreeSet;

import ch.unil.genescore.gene.Gene;
import ch.unil.genescore.main.Settings;
import ch.unil.genescore.vegas.AllOverlappedElements;
import ch.unil.genescore.vegas.BedFileStream;
import ch.unil.genescore.vegas.GenewiseSnpWeights;
import ch.unil.genescore.vegas.OverlappedCollectionStream;
import ch.unil.genescore.vegas.OverlappedGenomicElement;
import ch.unil.genescore.vegas.OverlappedGenomicElementStream;
import ch.unil.genescore.vegas.ReferencePopulation;
import ch.unil.genescore.vegas.SnpPositionStream;

public class SnpWeightCreator {

	
	private GenewiseSnpWeights weights_ = new GenewiseSnpWeights();
	
	public SnpWeightCreator(ReferencePopulation myRefPop, LinkedHashMap<String, Gene> genes){
		
		OverlappedGenomicElementStream mySnpPosStream = null;
		if (!Settings.chromosome_.equals(""))
			mySnpPosStream = myRefPop.getSnpPositionStream(Settings.chromosome_);	
		else{
			String[] chrStrings = new String[22]; 
			for (int i=1; i<23; i++){
				chrStrings[i-1]= "chr" + i;
			}
			mySnpPosStream = myRefPop.getSnpPositionStream(chrStrings);	
			
		}
		weights_ = fillWeights(mySnpPosStream, genes);
	
	}
	private GenewiseSnpWeights fillWeights(OverlappedGenomicElementStream mySnpPosStream, LinkedHashMap<String, Gene> genes){
			// set up background weights		
			//SnpPositionStream mySnpPosStream = myRefPop.getPositionStream(chr);			
			AllOverlappedElements snpsInGenes = new AllOverlappedElements(mySnpPosStream, genes.values());
			snpsInGenes.setMods(-Settings.bedBackgroundExtension_, Settings.bedBackgroundExtension_);
			snpsInGenes.fillTreeSet();
			GenewiseSnpWeights weights = snpsInGenes.getAsGenewiseSnpWeightDataStructure(Settings.bedBackgroundWeight_, 0);
			
			if  (!Settings.bedFilePath_.equals("")) {				
				mySnpPosStream.reOpenStream();
				//TODO: add stream containing id as well here
				BedFileStream myBedStream = new BedFileStream(Settings.bedFilePath_);	
				AllOverlappedElements snpsInBed = new AllOverlappedElements(mySnpPosStream, myBedStream);
				snpsInBed.fillTreeSet();
				TreeSet<OverlappedGenomicElement> myOverlappedGenomicRegions = snpsInBed.getMyLargeList();									
				OverlappedCollectionStream myOverlappedGenomicRegionStream = new OverlappedCollectionStream(myOverlappedGenomicRegions);
				AllOverlappedElements myOverlaps = new AllOverlappedElements(myOverlappedGenomicRegionStream,genes.values());
				myOverlaps.setMods(-Settings.bedBackgroundExtension_, Settings.bedBackgroundExtension_);
				myOverlaps.fillTreeSet();
				//TODO: pruneback option here: go through each TreeSet list and drop members inside that do not conform to gene val
				if(Settings.filterOnBed_){
					myOverlaps.filterTree();
				}
				
				GenewiseSnpWeights  augStruct = myOverlaps.getAsGenewiseSnpWeightDataStructure(Settings.bedWeight_, 1);
				weights.updateWeights(augStruct);
			}
			return weights;
		}
		
		
	public GenewiseSnpWeights getWeights(){return weights_;}
	
}
