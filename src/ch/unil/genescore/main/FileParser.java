package ch.unil.genescore.main;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import java.util.ArrayList;
import java.util.zip.GZIPInputStream;


/**
 * Parse a tab-separated value file (other separators than tab can also be set).
 * Note: java.util.Scanner offers more sophisticated functions to parse structured text files
 */
public class FileParser {

	/** The token used to separate the columns (usually a single character) */
	String separator_ = "\t";
	/** The buffered file reader */
	BufferedReader reader_ = null;
	/** Line counter */
	private int lineCounter_ = 0;
	/** Next line */
	private String nextLine_ = null;
	/** current line */
	private String currentLine_ = null;
	
	
	// ============================================================================
	// PUBLIC METHODS
	    
	/** Constructor */
	public FileParser(){
		
	}
	public FileParser(BufferedReader reader) {
		reader_=reader;
		readFirstLine();		
	}
	public FileParser(String filename) {

		try {
			System.out.println("Reading file: " + filename);

			if (filename.endsWith(".gz")) {
				InputStream fileStream = new FileInputStream(filename);
				InputStream gzipStream = new GZIPInputStream(fileStream);
				Reader decoder = new InputStreamReader(gzipStream);
				reader_ = new BufferedReader(decoder);
			} else {
				//FileInputStream fstream = new FileInputStream(filename);
				//DataInputStream in = new DataInputStream(fstream);
				//reader_ = new BufferedReader(new InputStreamReader(in));
				reader_ = new BufferedReader(new FileReader(filename));
			}
		} catch (Exception e) {
			Main.error(e);
		}
		readFirstLine();
		
	}
	
	public FileParser(String filename, String separator){
		this(filename);
		this.setSeparator(separator);		
	}
	
    // ----------------------------------------------------------------------------
	public void readFirstLine(){
		readLineAsString();
	};
	
	/** Read and return the next line, split using the separator_. Returns null if there is no more line to read. */
	public String[] readLine() {
		if(lineAvailable())
			return readLineAsString().split(separator_, -1);
		return null;
	}
	/** Read and return the next line, split using the separator_. Returns null if there is no more line to read. */
	public String readLineAsString() {
		String outLine = currentLine_;
		try {
			if (lineCounter_==0){
				currentLine_ = reader_.readLine();				
				outLine=currentLine_;
			}
			else{
				currentLine_=nextLine_;				
			}
			lineCounter_++;
			nextLine_ = reader_.readLine();
		} catch (IOException e) {
			Main.error(e);
		}
		
		//if (currentLine_ == null)
			//return null;
		
		return outLine;
	}

    // ----------------------------------------------------------------------------

	/** Read all the lines (see readLine()) */
	public ArrayList<String[]> readAll() {
		
		ArrayList<String[]> data = new ArrayList<String[]>();
		
		try {
			while (true) {
				nextLine_ = reader_.readLine();
				if (nextLine_ == null)
					break;
				
				lineCounter_++;
				data.add(nextLine_.split(separator_));
			}
		} catch (IOException e) {
			Main.error(e);
		}
		
		return data;
	}

	
    // ----------------------------------------------------------------------------
	public boolean lineAvailable(){
		return(currentLine_!=null);
	}
	
	
	/** Be polite and close the file reader when you're done */
	public void close() {
		
		try {
			reader_.close();
		} catch (IOException e) {
			Main.error(e);
		}
	}
	  
	
    // ----------------------------------------------------------------------------

	/** Skip N lines (useful to skip headers) */
	public void skipLines(int N) {
		
		try {
			for (int i=0; i<N; i++)
				reader_.readLine();
		} catch (IOException e) {
			Main.error(e);
		}
	}

	
    // ----------------------------------------------------------------------------

	/** Throw a runtime exception indicating the current line */
	public void error(String msg) {
		
		throw new RuntimeException("Line " + lineCounter_ + ": " + msg);
	}

	
	// ============================================================================
	// STATIC METHODS

	/** Count the number of lines in the given file */
	public static int countLines(String filename) {
	
		FileParser reader = new FileParser(filename);
		int count = 0;
		
		String[] nextLine = reader.readLine();
		while (nextLine != null) {
			count++;
			nextLine = reader.readLine();  
		}
		reader.close();

		return count;
	}
	  	
	
	// ============================================================================
	// SETTERS AND GETTERS

	public void setSeparator(String separator) { separator_ = separator; }
    
	public int getLineCounter() { return lineCounter_; }
		
}
