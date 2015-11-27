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
import java.io.IOException;

import ch.unil.genescore.main.Pascal;

public class SnpPositionStream
	extends OverlappedGenomicElementFileStream
	implements OverlappedGenomicElementStream {
	
	DataInputStream inStream_ = null;
	String binaryFileVersionID_= null;
	String filename_ = null;
	
	
	
	public SnpPositionStream(String filename,String binaryFileVersionID){
		System.out.println("reading snp positions from file:" + filename);
		DataInputStream inStream = null;
		try {						
			inStream = StreamMethods.openDataInputStream(filename, binaryFileVersionID);
			}
			catch(Exception e){
				Pascal.error(e, "Error loading SNP position (try deleting the binary reference population files)");
			}	
		inStream_ = inStream;
		streamOpen_= true;
		filename_ = filename;
		binaryFileVersionID_ = binaryFileVersionID;
	} 
	protected void loadElement()  {
		
		Snp currentSnp = null;
		if (!streamOpen_)
			throw new RuntimeException("never use on closed stream");
	
		String snpId = null;
		
		try {
			snpId = inStream_.readUTF();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		// End of file is marked with the verion ID
		if (snpId.equals(binaryFileVersionID_)){
			streamOpen_=false;
			currentElement_ = null;		
			try {
				inStream_.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			return;
		}		
		currentSnp = new Snp(snpId);		
		currentElement_ = new OverlappedGenomicElement(currentSnp);
		
		try {
			currentSnp.readPosAndAllele(inStream_);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} 
		catch (DataInconsistencyException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} 
	}
	
	public Snp loadSnp(){
		
		loadElement();		
		return ((Snp) currentElement_.getMainElement());		
	}
	
	public boolean streamOpen(){ return streamOpen_;}
	public void reOpenStream(){
		if(streamOpen()){throw new RuntimeException("never call reopen on open stream");}
		DataInputStream inStream = null;
		try {						
			inStream = StreamMethods.openDataInputStream(filename_, binaryFileVersionID_);
			}
			catch(Exception e){
				Pascal.error(e, "Error loading SNP position (try deleting the binary reference population files)");
			}	
		inStream_ = inStream;
		streamOpen_= true;
		streamInitialized_ = false;

		
	}
}
