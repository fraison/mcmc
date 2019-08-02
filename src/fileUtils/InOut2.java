package fileUtils;



import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Scanner;
import java.util.logging.Logger;


/**
 * Class to help reading tabulated files
 * @author fraison
 *
 */
public class InOut2 {	

	/**
	 * TODO extend this list
	 * @author fraison
	 *
	 */
	public static enum fileDelimiter {
	    COMMA,SPACE 
	}




	/**
	 * TODO extend this list
	 * @author fraison
	 *
	 */
	public static enum columnType {
	    STRING,INTEGER,DOUBLE 
	}



	private static final Logger LOGGER = Logger.getLogger( "InOut2");
	
	
	StringBuilder sb = new StringBuilder();
	
	
	
	
	/**
	 * Parse each table of String of a list according to the format of the elements
	 * @param mylist
	 * @param type the type of each column
	 * @return a vector of lines
	 */
	public static Object[][] listParserLine(ArrayList<String[]> mylist, columnType[] type){

		
		int nline = mylist.size();
		int ncol = mylist.get(0).length;
		
		Object[][] map = new Object[nline][ncol];
		
			//Pattern.compile(pattern);
		
		for (int k=0;k<ncol;k++){	
			switch(type[k]){
			case INTEGER:for (int l=0;l<nline;l++)	map[l][k]=Integer.parseInt(mylist.get(l)[k]);break;
			case DOUBLE:for (int l=0;l<nline;l++){
				//replace e-00x or e-0xx by e-x or e-xx
				String text = mylist.get(l)[k];
				//matcher= Pattern.compile(pattern1).matcher(text);				
				//String out = matcher.replaceFirst("e-");
				
				//matcher= Pattern.compile(pattern2).matcher(out);				
				//String out2 = matcher.replaceFirst("0.0");
				
				map[l][k]=Double.parseDouble(text)  ;
			}break; 
			case STRING:for (int l=0;l<nline;l++)	map[l][k]=mylist.get(l)[k]  ;break; 
			default:LOGGER.info("no case for this columnType");
			}

		}
		return map;
	}




	/**
	 * Read any table in text file in a very generic way.
	 * @param fileName
	 * @param headerNumberOfLines number of line in the header to be skipped 
	 * @param delimiter TODO
	 * @return an a list of line whose element are store into a table of Strings
	 * @throws FileNotFoundException
	 */
	public static ArrayList<String[]> readTextFileColumns(String fileName, int headerNumberOfLines, fileDelimiter delimiter) throws FileNotFoundException{
	
	
		Scanner fscan = new Scanner(new File(fileName));//scanner for the lines
		ArrayList<String[]>listOfLine = new ArrayList<String[]>();//list of line of integers
		String myline;
	
		//skip header
		for (int i =0;i<headerNumberOfLines;i++)fscan.nextLine() ;//System.out.println(fscan.nextLine());
		
		//if (delimiter == space) 
			String delim= "\\s+";
			switch(delimiter){
			case COMMA:delim = ","; break;
			case SPACE:delim = "\\s+"; break; 			
			default:InOut2.LOGGER.info("no case for this fileDelimiter");
			}
	
		//core loop on lines
		while (fscan.hasNextLine()){
			myline  = fscan.nextLine();
			if (myline.equals("")){
				InOut2.LOGGER.info("empty line");
				break;
			}
			Scanner lscan = new Scanner(myline).useDelimiter(delim);//scanner for each string on the line
			ArrayList<String> listOfComponent = new ArrayList<String>()	;//list of strings for each line
	
			//loop on columns
			while (lscan.hasNext()){
	
				String s = lscan.next();
				listOfComponent.add(s);
			}
	
	
			String[] lineVector = new String[listOfComponent.size()];
			listOfComponent.toArray(lineVector)	; //convert  list of String into 1D array
	
			listOfLine.add(lineVector);//add the array of String to the list		
	
	
		}
		fscan.close();
		return listOfLine;
	
	}




}
