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

import ch.unil.genescore.main.Main;

public class SnpSerializedPositionStream {

	private int lineCounter_ = 0;
	/** Next line */
	private Snp nextLine_ = null;
	/** current line */
	private Snp currentLine_ = null;
	
	DataInputStream inStream_ = null;
	String binaryFileVersionID_= null;
	String filename_ = null;
	boolean streamOpen_ = true;
	
	
	public SnpSerializedPositionStream(String filename,String binaryFileVersionID){
		System.out.println("reading snp positions from file:" + filename);
		DataInputStream inStream = null;
		try {						
			inStream = StreamMethods.openDataInputStream(filename, binaryFileVersionID);
			}
			catch(Exception e){
				Main.error(e, "Error loading SNP position (try deleting the binary reference population files)");
			}	
		inStream_ = inStream;
		streamOpen_= true;
		filename_ = filename;
		binaryFileVersionID_ = binaryFileVersionID;
		readLine();
	}
		public boolean lineAvailable(){
			return(currentLine_!=null);
		}
		
		private Snp createSnpFromStream(){
			
			if (!streamOpen_)
				throw new RuntimeException("Don't call on closed stream");
		
			String snpId = null;
			
			try {
				snpId = inStream_.readUTF();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
				System.exit(-1);
			}
			// End of file is marked with the verion ID
			if (snpId.equals(binaryFileVersionID_)){
				streamOpen_=false;
				
				try {
					inStream_.close();
				} catch (IOException e) {
					Main.error(e);
				}
				return null;
			}		
			Snp currentSnp = new Snp(snpId);		
			
			try {
				currentSnp.readPosAndMinorAllele(inStream_);
			} catch (IOException e) {
				Main.error(e, "Error loading SNP position (try deleting the binary reference population files)");				
			} 
			catch (DataInconsistencyException e) {
				Main.error(e);
			 }
			return(currentSnp);
		}
					
		public Snp readLine() {
			Snp outLine = currentLine_;
			if (lineCounter_==0){
				currentLine_ = createSnpFromStream();				
				outLine=currentLine_;
			}
			else{
				currentLine_=nextLine_;				
			}
			lineCounter_++;
			if (streamOpen_)
				nextLine_ = createSnpFromStream();				
			return outLine;	
		} 
}
