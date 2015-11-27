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

import java.io.DataInputStream;
import java.util.ArrayList;

import ch.unil.genescore.main.Pascal;

public class WrappedSnpPositionStream 
//extends OverlappedGenomicElementFileStream
implements OverlappedGenomicElementStream {
	
	private SnpPositionStream posStream_ = null;
	private int currentIndex_ = 0;	
	private ArrayList<String> filenames_ = null;
	private String binaryFileVersionID_ = null;
	
	WrappedSnpPositionStream(ArrayList<String> filenames,String binaryFileVersionID){
		filenames_ = filenames;
		binaryFileVersionID_=binaryFileVersionID;
		currentIndex_=0;
		posStream_ = new SnpPositionStream(filenames_.get(currentIndex_),binaryFileVersionID_);
	}
	@Override
	public OverlappedGenomicElement getNextAsOverlappedGenomicElement() {
		boolean currentStreamOpen = posStream_.streamOpen();
		if (!currentStreamOpen) {
			if ((currentIndex_+1)==filenames_.size())
				throw new RuntimeException("never call on closed stream.");			
			currentIndex_++;
			posStream_ = new SnpPositionStream(filenames_.get(currentIndex_),binaryFileVersionID_);			
		}
		return posStream_.getNextAsOverlappedGenomicElement();			
	}

	@Override
	public boolean streamOpen() {
		// TODO Auto-generated method stub
		return (posStream_.streamOpen() || (currentIndex_+1)!=filenames_.size());
	}
	@Override
	public void reOpenStream() {
		if(streamOpen()){throw new RuntimeException("never call reopen on open stream");}
		currentIndex_=0;
		posStream_ = new SnpPositionStream(filenames_.get(currentIndex_),binaryFileVersionID_);
	}

	
}


