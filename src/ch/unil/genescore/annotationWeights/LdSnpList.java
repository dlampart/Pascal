package  ch.unil.genescore.annotationWeights;
import java.util.ArrayDeque;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;

import ch.unil.genescore.gene.Chromosome;
import ch.unil.genescore.gene.GenomicElement;
import ch.unil.genescore.vegas.LinkageDisequilibrium;
import ch.unil.genescore.vegas.Snp;

public class LdSnpList {
HashMap<String, LinkedList<SnpLdPair>> chromosome_ = new HashMap<String,LinkedList<SnpLdPair>>();
	private int rangeCutoff_ = 0;
	private int ldCutoff_ = 0;
	public LdSnpList() {
		
		// TODO Auto-generated constructor stub
	}
	public void setupList(Chromosome snpChromosome){
		ArrayDeque<SnpLdPairs> currentLDsnps = new ArrayDeque<SnpLdPairs>();
		Snp curSnp;
		Collection<GenomicElement> iterable = snpChromosome.getElements(); 
		Iterator<GenomicElement> it = iterable.iterator();
		
		int c1=0;
		int c2=0;
		//assert(currentLDsnps.peekFirst()==null);
		//if(it.hasNext()){
		//	currentLDsnps.addFirst(new SnpLdPairs((Snp)it.next()));				
		//}
		
		while(it.hasNext()){
			System.out.println(c1++);
			curSnp=(Snp)it.next();
			currentLDsnps.addFirst(new SnpLdPairs(curSnp));
			int currentLastPos=currentLDsnps.getLast().getSnp().getStart();
			int currentEndPos=curSnp.getStart();
			//if (currentLastPos + rangeCutoff_ < currentEndPos){
				
			while(currentEndPos - currentLastPos > rangeCutoff_ ){
				SnpLdPairs SnpLDpairsToAdd = currentLDsnps.removeLast();
				chromosome_.put(SnpLDpairsToAdd.getSnpId(), SnpLDpairsToAdd.getList());				
				currentLastPos=currentLDsnps.getLast().getSnp().getStart();
			}
			assert(!currentLDsnps.isEmpty()); 	
			//SnpLdPairs newElement = new SnpLdPairs(curSnp);
			double LDval;
			Iterator<SnpLdPairs> it2 = currentLDsnps.iterator();
			SnpLdPairs iteratorPos = null;
			c2=0;
			while (it2.hasNext()){
				System.out.println(c2++);
				//System.out.println(passing through);
				iteratorPos = it2.next();					
				LDval = LinkageDisequilibrium.computeCorrelation(iteratorPos.getSnp(), curSnp);
				iteratorPos.appendToLDList(curSnp.getId(), LDval);
				if (iteratorPos.getSnp() != curSnp){//
					System.out.println("sfadasfd");
					System.out.println(iteratorPos.getSnp() != curSnp);
					currentLDsnps.getFirst().appendToLDList(iteratorPos.getSnp().getId(), LDval);					
					
				}
			}					
			
		}
		while(!currentLDsnps.isEmpty()){
			SnpLdPairs SnpLDpairsToAdd = currentLDsnps.removeLast();
			chromosome_.put(SnpLDpairsToAdd.getSnpId(), SnpLDpairsToAdd.getList());				
		}
	}
	
	//public double
		
	public void setCutoffs(int rangeCutoff,int ldCutoff){		
		rangeCutoff_ = rangeCutoff;
		ldCutoff_ = ldCutoff;		
	}
	
	

}
