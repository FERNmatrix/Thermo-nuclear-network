// -----------------------------------------------------------------------------------
//  Class DataHolder to hold various important global data
//  in static variables
// -----------------------------------------------------------------------------------

class DataHolder {
	
	static int Znum = 110;
	static int Nnum = 200;
	
	/* 
	Dimensionality of 
            reacNum in DataHolder.java
            maxCases in PlotReactionList.java
            Color [] col = new Color [151] in RatePlotFrame.java
            imax in RatePlotFrame.java
    should be equal and determine the max number of reactions that
    can be handled by the rateplotter.
	*/
	
	static int reacNum = 151;

    // Boolean array for active reaction classes (1-8)

    static boolean [][][] includeReaction = new boolean[Znum][Nnum][9];

    // Boolean array for active reaction components

    static boolean [][][] RnotActive = new boolean[Znum][Nnum][reacNum];

    // Boolean array indicating whether isotope has been opened for
    // individual reaction selection and saved

    static boolean [][] wasOpened = new boolean[Znum][Nnum];

    // Arrays for case where which reactions to include is read in from file.  These
    // will be initialized in ChooseActiveRates.

    static int [][][] useReadRates = new int [Znum][Nnum][reacNum];
    static int [][] maxReadRates = new int [Znum][Nnum];
	
	static int [][][] activeReactionsSerialIndex = new int [Znum][Nnum][reacNum];


    // ---------------------------------------------------------------------------------
    //  Empty constructor.  This class only used to hold global data in
    //  static variables so it doesn't need to be instantiated.
    // ---------------------------------------------------------------------------------

    public DataHolder (int Z, int N) {

    }

}
