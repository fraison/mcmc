package fileUtils;


import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.LinkedHashMap;
import java.util.logging.Logger;

import fileUtils.InOut2.fileDelimiter;

/**
 * Class to store tabular results to text file or display on screen need to add
 * table names (header), table names format and value format (
 * {@code addValFormat}) need to add values({@code addVa})
 * 
 * @author fraison
 * 
 */
public class Writer {

	/**
	 * Logger
	 */

	private static final Logger LOGGER = Logger.getLogger( "Writer");
	
	LinkedHashMap<String, String> fieldformat = new LinkedHashMap<String, String>();
	LinkedHashMap<String, Object> fieldvalues = new LinkedHashMap<String, Object>();
	LinkedHashMap<String, String> headerformat = new LinkedHashMap<String, String>();

	PrintStream outFile = null;
	StringBuilder sb;
	boolean header = true;
	String outfullfileName = null;

	/**
	 * Constructor
	 * 
	 * @param outfullfileName
	 *            in case want to save to a file can be null
	 * @throws FileNotFoundException
	 */
	public Writer(String outfullfileName) throws FileNotFoundException {
		if (outfullfileName != null) {
			this.outfullfileName = outfullfileName;
			outFile = new PrintStream(outfullfileName);
		}
	}

	/**
	 * Add a format a table name, a table name format and a format for a value
	 * 
	 * @param tableName
	 *            a table column tname
	 * @param tableValFormat
	 *            table column value format
	 * @param tableNameFormat
	 *            table column name format
	 */
	public void addValFormat(String tableName, String tableValFormat,
			String tableNameFormat) {
		fieldformat.put(tableName, tableValFormat);
		headerformat.put(tableName, tableNameFormat);
	}

	/**
	 * Add a table value
	 * 
	 * @param tabName
	 *            table column name
	 * @param sval
	 *            table column value
	 */
	public void addVal(String tabName, Object sval) {

		fieldvalues.put(tabName, sval);
	}

	/**
	 * add a line to file
	 * 
	 * @throws Exception
	 *             when file name is not set
	 */
	public void vals2File() throws Exception {
		if (outFile != null) {
			if (header) {
				String head = getHeader();
				outFile.print(head);
				header = false;
			}
			String line = getFormatedVals();
			outFile.print(line);
		} else
			throw new Exception("file name not set in writer");

	}

	/**
	 * add a line to file, with a delimiter
	 * 
	 * @param delim
	 *            file delimiter
	 * @throws Exception
	 */
	public void vals2File(fileDelimiter delim) throws Exception {
		if (outFile != null) {
			if (header) {
				String head = getHeader(delim);
				outFile.print(head);
				header = false;
			}
			String line = getFormatedVals(delim);
			outFile.print(line);
		} else
			throw new Exception("file name not set in writer");

	}
	
	/**
	 * add an empty line to file
	 * @throws Exception
	 */
	public void emptyLine2File( ) throws Exception {
		if (outFile != null) {
			outFile.print("\n");
		} else
			throw new Exception("file name not set in writer");

	}
	

	/**
	 * add a line to file without a masked column
	 * 
	 * @throws Exception
	 *             when file name is not set
	 */
	public void vals2File(String masked) throws Exception {
		if (outFile != null) {
			if (header) {
				String head = getHeader(masked);
				outFile.print(head);
				header = false;
			}
			String line = getFormatedVals(masked);
			outFile.print(line);
		} else
			throw new Exception("file name not set in writer");

	}
	

	
	
	/**
	 * add a line to file without a masked column
	 * 
	 * @throws Exception
	 *             when file name is not set
	 */
	public void vals2File(fileDelimiter delim , String masked) throws Exception {
		if (outFile != null) {
			if (header) {
				String head = getHeader(delim , masked);
				outFile.print(head);
				header = false;
			}
			String line = getFormatedVals(delim , masked);
			outFile.print(line);
		} else
			throw new Exception("file name not set in writer");

	}

	/**
	 * send values to terminal
	 * 
	 * @param printHeader
	 */
	public void vals2Term(boolean printHeader) {
		if (printHeader) {
			String head = getHeader();
			System.out.print(head);
		}
		String line = getFormatedVals();
		System.out.print(line);
	}
	
	

	/**
	 * send values to terminal with delimiter specified
	 * 
	 * @param delim
	 *            file delimiter
	 * @param printHeader
	 */
	public void vals2Term(fileDelimiter delim, boolean printHeader) {
		if (printHeader) {
			String head = getHeader(delim);
			System.out.print(head);
		}
		String line = getFormatedVals(delim);
		System.out.print(line);
	}
	
	/**
	 * send empty line to the terminal
	 * 
	 * @param delim
	 *            file delimiter
	 * @param printHeader
	 */
	public void emptyLine2Term() {
		System.out.print("\n");
	}

	/**
	 * send values to terminal without a masked column trick to avoid a column
	 * 
	 * @param printHeader
	 * @param masked
	 *            name of the key (can be partial)
	 */
	public void vals2Term(boolean printHeader, String masked) {
		if (printHeader) {
			String head = getHeader(masked);
			System.out.print(head);
		}
		String line = getFormatedVals(masked);
		System.out.print(line);
	}
	
	/**
	 * send values to terminal without a masked column trick to avoid a column
	 * 
	 * @param printHeader
	 * @param masked
	 *            name of the key (can be partial)
	 */
	public void vals2Term(fileDelimiter delim, boolean printHeader, String masked) {
		if (printHeader) {
			String head = getHeader(delim , masked);
			System.out.print(head);
		}
		String line = getFormatedVals(delim, masked);
		System.out.print(line);
	}

	/**
	 * getter for file name used to save data
	 * 
	 * @return filename
	 */
	public String getFileName() {
		return this.outfullfileName;
	}

	/**
	 * 
	 * @return String of formated values
	 */
	private String getFormatedVals() {
		sb = new StringBuilder();
		for (String key : fieldformat.keySet()) {
			sb.append(String.format(fieldformat.get(key).toString(),
					fieldvalues.get(key)));

		}
		sb.append("\n");

		return sb.toString();

	}

	/**
	 * 
	 * @return String of formated values
	 */
	private String getFormatedVals(String masked) {
		sb = new StringBuilder();
		for (String key : fieldformat.keySet()) {
			if (!key.contains(masked))
				sb.append(String.format(fieldformat.get(key).toString(),
						fieldvalues.get(key)));

		}
		sb.append("\n");

		return sb.toString();

	}

	/**
	 * 
	 * @param delimiterType
	 * @param masked
	 * @return
	 */
	private String getFormatedVals(fileDelimiter delimiterType, String masked) {

		String delim = "\\s+";
		switch (delimiterType) {
		case COMMA:
			delim = ",";
			break;
		case SPACE:
			delim = "\\s+";
			break;
		default:
			LOGGER.info("no case for this fileDelimiter");
		}

		sb = new StringBuilder();
		int count = 0; // counter
		for (String key : fieldformat.keySet()) {
			if (!key.contains(masked))
				sb.append(String.format(fieldformat.get(key).toString(),
						fieldvalues.get(key)));
			// append delimiter as long as is not last column
			if (count < fieldformat.keySet().size() - 1) {
				sb.append(delim);
			}
			count++;
		}
		sb.append("\n");

		return sb.toString();

	}

	/**
	 * 
	 * @return String of formated values
	 */
	private String getFormatedVals(fileDelimiter delimiterType) {

		String delim = "\\s+";
		switch (delimiterType) {
		case COMMA:
			delim = ",";
			break;
		case SPACE:
			delim = "\\s+";
			break;
		default:
			LOGGER.info("no case for this fileDelimiter");
		}

		sb = new StringBuilder();
		int count = 0; // counter
		for (String key : fieldformat.keySet()) {
			sb.append(String.format(fieldformat.get(key).toString(),
					fieldvalues.get(key)));
			// append delimiter as long as is not last column
			if (count < fieldformat.keySet().size() - 1) {
				sb.append(delim);
			}
			count++;
		}
		sb.append("\n");

		return sb.toString();

	}

	/**
	 * 
	 * @return formated table names
	 */
	private String getHeader() {
		sb = new StringBuilder();
		for (String key : headerformat.keySet()) {
			sb.append(String.format(headerformat.get(key).toString(), key));

		}
		sb.append("\n");

		return sb.toString();
	}

	/**
	 * 
	 * @return formated table names
	 */
	private String getHeader(fileDelimiter delimiterType) {

		String delim = "\\s+";
		switch (delimiterType) {
		case COMMA:
			delim = ",";
			break;
		case SPACE:
			delim = "\\s+";
			break;
		default:
			LOGGER.info("no case for this fileDelimiter");
		}

		sb = new StringBuilder();
		int count = 0; // counter
		for (String key : headerformat.keySet()) {
			sb.append(String.format(headerformat.get(key).toString(), key));
			// append delimiter as long as is not last column
			if (count < fieldformat.keySet().size() - 1) {
				sb.append(delim);
			}
			count++;

		}
		sb.append("\n");

		return sb.toString();
	}

	/**
	 * 
	 * @return formated table names
	 */
	private String getHeader(String masked) {
		sb = new StringBuilder();
		for (String key : headerformat.keySet()) {
			if (!key.contains(masked))
				sb.append(String.format(headerformat.get(key).toString(), key));

		}
		sb.append("\n");

		return sb.toString();
	}
	
	

	private String getHeader(fileDelimiter delimiterType , String masked) {
		
		String delim = "\\s+";
		switch (delimiterType) {
		case COMMA:
			delim = ",";
			break;
		case SPACE:
			delim = "\\s+";
			break;
		default:
			LOGGER.info("no case for this fileDelimiter");
		}
		
		sb = new StringBuilder();
		int count = 0; // counter
		for (String key : headerformat.keySet()) {
			if (!key.contains(masked))
				sb.append(String.format(headerformat.get(key).toString(), key));
			if (count < fieldformat.keySet().size() - 1) {
				sb.append(delim);
			}
		}
		sb.append("\n");

		return sb.toString();
	}
	

	/**
	 * flush file
	 */
	public void flush() {
		if (outFile != null)
			outFile.flush();
	}

	/**
	 * close file
	 */
	public void close() {
		if (outFile != null)
			outFile.close();
	}
	
	/**
	 * delete file
	 */
    public void deletefile() {
;        File file = new File(outfullfileName);
         if (file.delete())
             System.out.println("File deleted successfully"); 
         else         
        	 System.out.println("Failed to delete the file"); 	
    	
    }
    
    
	/**
	 * example showing how to use this class the writer object can be passed to
	 * any method
	 * 
	 * @param args
	 */
	public static void main(String args[]) {

		Writer w;
		try {
			w = new Writer("/home/fraison/tmp/writertry.txt");
			w.addValFormat("density", "%14.6e", "%14s");
			w.addValFormat("flux", "%12.2f", "%12s");
			w.addVal("density", 10.);
			w.addVal("density", 30.);
			w.addVal("flux", 20.);

			// FileUtilities.wtry(w);
			w.vals2Term(true);
			w.vals2File();
			w.vals2Term(true, "densit");

		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

}
