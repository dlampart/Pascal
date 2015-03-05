/*******************************************************************************
 * Copyright (c) 2015 IBM Corporation and others.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Eclipse Public License v1.0
 * which accompanies this distribution, and is available at
 * http://www.eclipse.org/legal/epl-v10.html
 *
 * Contributors:
 *     IBM Corporation - initial API and implementation
 *******************************************************************************/
package ch.unil.genescore.pathway;

import java.util.Collection;
import java.util.Deque;
import java.util.LinkedList;
import java.util.ListIterator;
import java.util.Queue;
import java.util.Stack;
import java.util.TreeSet;

import ch.unil.genescore.gene.Gene;
import ch.unil.genescore.main.Settings;

public class Damper {
	/** Damper: takes in Collection of genes and lowers chi-square scores of genes that have neighbour with higher chi-sq-scores
	 * nearby. Neighbourhood is defined as having boundaries no further away from each other than 'dist_'.
	 * deflation of values is proceeds as: newVal=oldVal*(oldVal/maxVal)^(deflationRate_).
	 * 
	 */

	TreeSet<Gene> mySet_ = null;
	Deque<Gene> myGeneDeque_=null;
	//Deque<Double> myMaxDeque_=null;
	LinkedList<Double> myMaxDeque_=null;
	int dist_;
	double deflationRate_;
	public Damper(Collection genes, int dist, double deflationRate){
		
		mySet_=new TreeSet<Gene>(genes);
		dist_ = dist; 
		deflationRate_=deflationRate;
		
	}
	
	
	public void dampSet(){
	
		myGeneDeque_ = new LinkedList<Gene>();
		myMaxDeque_ = new LinkedList<Double>();
		for (Gene el : mySet_){
			while(!dequeFits(el)){
				Gene poppedGene=null;
				Double poppedMax=null;
				poppedGene = myGeneDeque_.removeFirst();
				poppedMax = myMaxDeque_.removeFirst();
				double deflationFactor=Math.pow(poppedGene.getChi2Stat()/poppedMax,deflationRate_);
				poppedGene.setChi2Stat(poppedGene.getChi2Stat()*deflationFactor);
			}			
			if (dequeFits(el)){
				
				double dequeMax=0;
				ListIterator<Double> iter = myMaxDeque_.listIterator();
				while (iter.hasNext()){
					
					Double d = iter.next();
													
					if (d < el.getChi2Stat())
						iter.set(el.getChi2Stat());					
					if (d>dequeMax)
						dequeMax=d;
					
				}				
				
				
				Double topOf=Math.max(dequeMax,el.getChi2Stat());								
				myGeneDeque_.addLast(el);
				myMaxDeque_.addLast(topOf);
				
			}
		}
		while (!myGeneDeque_.isEmpty()){
			Gene poppedGene=null;
			Double poppedMax=null;
			poppedGene = myGeneDeque_.removeFirst();
			poppedMax = myMaxDeque_.removeFirst();
			double deflationFactor=Math.pow(poppedGene.getChi2Stat()/poppedMax,deflationRate_);
			poppedGene.setChi2Stat(poppedGene.getChi2Stat()*deflationFactor);			
		}
	}
	
	private boolean dequeFits(Gene g) {
		if (myGeneDeque_.isEmpty())
			return true;
		
		//boolean b =(myGeneDeque_.getFirst().chr_.equals(g.chr_) && myGeneDeque_.getFirst().end_ > g.start_- dist_);
		boolean b =(myGeneDeque_.getFirst().chr_.equals(g.chr_) && myGeneDeque_.getFirst().end_ > g.start_- dist_);
		return b;
		
	}

}
