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

import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
//TODO: Dangerous class I think
public class OverlappedCollectionStream 
implements OverlappedGenomicElementStream {
	
	Collection<OverlappedGenomicElement> myCollection_ = null;
	OverlappedGenomicElement lastElement_ = null;
	OverlappedGenomicElement currentElement_ = null;
	boolean streamOpen_ = false;
	Iterator<OverlappedGenomicElement> it_ = null;
	public OverlappedCollectionStream(Collection<OverlappedGenomicElement> myCollection){
		myCollection_ = myCollection;
		streamOpen_ = true;
		it_ = myCollection_.iterator();
	}
	
	
	public OverlappedGenomicElement getNextAsOverlappedGenomicElement(){
		if (!streamOpen_)
			throw new RuntimeException("never use on closed stream");
		lastElement_ = currentElement_;
		if (it_.hasNext()){
			currentElement_=it_.next();
			
		}
		if (!it_.hasNext()) {			
			streamOpen_ = false;
		}
		if (lastElement_!=null){
			if (currentElement_.compareTo(lastElement_)<0)
				throw new RuntimeException("elements are not sorted right");
		}
		return currentElement_;
	}
	public boolean streamOpen(){
		return streamOpen_;
	}


	@Override
	public void reOpenStream() {
		if(streamOpen()){throw new RuntimeException("never call reopen on open stream");}
		OverlappedGenomicElement lastElement_ = null;
		OverlappedGenomicElement currentElement_ = null;
		streamOpen_ = true;
		it_ = myCollection_.iterator();		
	}

}
