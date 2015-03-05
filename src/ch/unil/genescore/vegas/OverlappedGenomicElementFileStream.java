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
package ch.unil.genescore.vegas;

public abstract class OverlappedGenomicElementFileStream {
	
	OverlappedGenomicElement currentElement_ = null;
	boolean streamInitialized_ = false;
	boolean streamOpen_= false;
	
	/** updates currentElement_ 
	 * @throws  */
	abstract protected void loadElement();
	
	public OverlappedGenomicElement getNextAsOverlappedGenomicElement(){	
		if (!streamOpen_)
			throw new RuntimeException("never use on closed stream");
			if (!streamInitialized_){//initial load
				
				loadElement();
				streamInitialized_=true;
			}
			//System.out.println(currentElement_.mainElement_.id_);
			OverlappedGenomicElement elementToBeReturned = currentElement_;
				loadElement();
			
			if (currentElement_!=null){
			//	System.out.println(currentElement_.mainElement_.id_);
				if (currentElement_.compareTo(elementToBeReturned)<0){
					System.out.println(currentElement_.mainElement_.id_);
					System.out.println(elementToBeReturned.mainElement_.id_);
					throw new RuntimeException("stream not sorted properly");
				}
			}	
			return  elementToBeReturned;
		
	}

}
