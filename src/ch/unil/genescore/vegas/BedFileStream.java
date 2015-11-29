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

import java.io.File;

import ch.unil.genescore.gene.GenomicElement;
import ch.unil.genescore.main.Pascal;
import ch.unil.gpsutils.FileParser;


public class BedFileStream 
extends OverlappedGenomicElementFileStream
implements OverlappedGenomicElementStream {
	static String fakeIds = ".";
	FileParser parser_ = null;
	File file_ = null;
	//OverlappedGenomicElement currentElement_ = null;
	
	
	public BedFileStream(File file){
		file_= file;
    	parser_ = new FileParser(Pascal.log, file);
    	streamOpen_=true;
	}

	
	public void loadElement(){
		if (!streamOpen_)
			throw new RuntimeException("never use on closed stream");
		String[] nextLine = parser_.readLine();
		if (nextLine==null){
			currentElement_ = null;
			parser_.close();
			streamOpen_=false;	
			return;
		}
		GenomicElement currentElement = null;
		if (nextLine.length>3){
			currentElement =  new GenomicElement(nextLine[3]);
		}
		else {
			currentElement =  new GenomicElement(fakeIds);
		}
		currentElement.setPosition(nextLine[0], Integer.parseInt(nextLine[1]), Integer.parseInt(nextLine[2]), true);
		currentElement_= new OverlappedGenomicElement(currentElement);			
	}
	public boolean streamOpen(){ return streamOpen_;}
	@Override
	public void reOpenStream() {
		
		if(streamOpen()){throw new RuntimeException("never call reopen on open stream");}
    	parser_ = new FileParser(Pascal.log, file_);
    	streamOpen_=true;
		// TODO Auto-generated method stub
		
	}	
}
