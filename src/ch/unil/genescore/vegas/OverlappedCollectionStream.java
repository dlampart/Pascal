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
