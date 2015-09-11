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
package ch.unil.genescore.vegas.test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.Reader;
import java.io.StringReader;
import java.util.ArrayList;

import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import ch.unil.genescore.main.FileParser;
import ch.unil.genescore.main.Settings;

public class FileParserTest {
	// ============================================================================
			// SETUP
			
			@BeforeClass
			public static void testSetup() {
				Settings.loadSettings();
			}

			@AfterClass
			public static void testCleanup() { }
			
			
			// ============================================================================
			// TESTS
			class BufferedReaderStub extends BufferedReader{
				
				private ArrayList<String> myOut = new ArrayList<String>();
				int counter=0;				
				public BufferedReaderStub() {
					super(new StringReader("blub.txt"));
					
				}
				public BufferedReaderStub(Reader rd) {
					super(rd);					
				}				
				public void setMyOut(ArrayList<String> out){
					myOut=out;
				}
				@Override
				public String readLine(){
					if (counter < myOut.size()){
						String outStr= myOut.get(counter);
						counter++;
						return outStr;
					}
					else{
						return null;	
						}
					}
				}
				
			
			
			
			@Test
			public void fileParserTest(){
				ArrayList<String> myStrs= new ArrayList<String>();
				String str1="first\tfirst2";
				myStrs.add(str1);
				String str2="second\tsecond2";
				myStrs.add(str2);
				String str3="third\tthird2";
				myStrs.add(str3);
				BufferedReaderStub myStub = new BufferedReaderStub();
				myStub.setMyOut(myStrs);
				FileParser fp= new FileParser(myStub);
				
				assertTrue(fp.readLineAsString().equals(str1));
				assertTrue(fp.lineAvailable());
				assertTrue(fp.readLineAsString().equals(str2));
				assertTrue(fp.lineAvailable());
				assertTrue(fp.readLineAsString().equals(str3));
				assertTrue(!fp.lineAvailable());
				assertTrue(fp.readLineAsString()==null);
				assertTrue(!fp.lineAvailable());
				
				myStub = new BufferedReaderStub();
				myStub.setMyOut(myStrs);
				fp= new FileParser(myStub);
				String[] line=fp.readLine();
				assertTrue(line[0].equals("first"));
				assertTrue(fp.lineAvailable());
				line=fp.readLine();
				assertTrue(line[0].equals("second"));
				assertTrue(fp.lineAvailable());
				line=fp.readLine();
				assertTrue(line[0].equals("third"));
				assertTrue(!fp.lineAvailable());
				assertTrue(fp.readLine()==null);
			}
}
