package com.vitimage.mrutils;

import java.util.ArrayList;

import com.vitimage.common.TransformUtils;
import com.vitimage.common.VitimageUtils;

import ij.macro.Interpreter;
import ij.macro.Program;
import ij.macro.Tokenizer;

public class SimplexDualCurveFitter{
	public static final boolean debugBionano=false;
	public static  int IterFactor = 500;
	protected double sigma;   
	private static final double alpha = -1.0;	  // reflection coefficient
	private static final double beta = 0.5;	  // contraction coefficient
	private static final double gamma = 2.0;	  // expansion coefficient
	private static final double root2 = 1.414214; // square root of 2
	
	protected int fit;                // Number of curve type to fit
	protected double[] trData,teData, magData;  // x,y data to fit
	protected int numPoints;          // number of data points
	protected int numParams;          // number of parametres
	protected int numVertices;        // numParams+1 (includes sumLocalResiduaalsSqrd)
	private int worst;			// worst current parametre estimates
	private int nextWorst;		// 2nd worst current parametre estimates
	private int best;			// best current parametre estimates
	protected double[][] simp; 		// the simplex (the last element of the array at each vertice is the sum of the square of the residuals)
	protected double[] next;		// new vertex to be tested
	private int numIter;		// number of iterations so far
	private int maxIter; 	// maximum number of iterations per restart
	private int restarts; 	// number of times to restart simplex after first soln.
	private static int defaultRestarts = 2;  // default number of restarts
	private int nRestarts;  // the number of restarts that occurred
	private static double maxError = 1e-10;    // maximum error tolerance
	private double[] initialParams;  // user specified initial parameters
	private long time;  //elapsed time in ms
	private String customFormula;
	private static int customParamCount;
	private double[] initialValues;
	private Interpreter macro;
	private double[] bionanoParams;
	private double[][]magData2D;
	private double[]tabTrVals;
	private double[]tabTeVals;
	private int []tabTrSeriesLength;
	private double t2;
	private double m0t2;
	private double r2;
	private double bionanoFactor=1+Math.exp(-0.25)+Math.exp(-0.5)+Math.exp(-0.75);
	private double t1;
	private double m0t1;
	
	public static double sigmaWay(double valFunk,double sigma){
		return valFunk*valFunk-2*sigma*sigma;
	}
	public static double besFunkCost(double value,double sigma) {
		double alpha=value*value/(4*sigma*sigma);
		return Math.sqrt(Math.PI*sigma*sigma/2.0)*( (1+2*alpha)*bessi0NoExp(alpha) + 2*alpha*bessi1NoExp(alpha) );
	}	
	public static double besFunkDeriv(double value,double sigma) {
		double alpha=value*value/(4*sigma*sigma);
		return Math.sqrt(Math.PI*alpha/2.0)*(bessi0NoExp(alpha) + bessi1NoExp(alpha) );
	}
	public static double bessi0NoExp( double x ){
	   double ax,ans;
	   double y;
	   ax=Math.abs(x);
	   if (ax < 3.75) {
	      y=x/3.75;
	      y=y*y;
	      ans=1/Math.exp(ax)*(1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492+y*(0.2659732+y*(0.0360768+y*0.0045813))))));
	   } else {
	      y=3.75/ax;
	      ans=(1/Math.sqrt(ax))*(0.39894228+y*(0.01328592+y*(0.00225319+y*(-0.00157565+y*(0.00916281+y*(-0.02057706+y*(0.02635537+y*(-0.01647633+y*0.00392377))))))));
	   }
	   return ans;
	}
	//Fonction de bessel modifiee de premiere espece d'ordre 1
	public static double bessi1NoExp( double x){
		double ax,ans;
		double y;
		ax=Math.abs(x);
		if (ax < 3.75) {
	      y=x/3.75;
	      y=y*y;
	      ans=ax/Math.exp(ax)*(0.5+y*(0.87890594+y*(0.51498869+y*(0.15084934+y*(0.02658733+y*(0.00301532+y*0.00032411))))));
	   } else {
	      y=3.75/ax;
	      ans=0.02282967+y*(-0.02895312+y*(0.01787654-y*0.00420059));
	      ans=0.39894228+y*(-0.03988024+y*(-0.00362018+y*(0.00163801+y*(-0.01031555+y*ans))));
	      ans *= (1/Math.sqrt(ax));
	   }
	   return x < 0.0 ? -ans : ans;
	}

    /** Construct a new SimplexCurveFitter. */
    public SimplexDualCurveFitter (double[] trData, double[]teData,double[] magData, int fitType,double sigma) {
		this.sigma=sigma;
        this.trData = trData;
        this.teData = teData;
        this.magData = magData;
        numPoints = trData.length;
        this.fit=fitType;
        initialize(fit);
    }
    
    
    public void computeBioNano(){
    	initialize2DTab();
    	inspect2DTab();
    	computeT2Bionano();
    	inspectT2EstimatedValues(); 
    	computeT1Bionano();
    }
    
	public void initialize2DTab(){
    	ArrayList<Double>listTr=new ArrayList<Double>();
    	ArrayList<Double>listTe=new ArrayList<Double>();
    	//Build a user-friendly structure to host the data
 		int nt1=0;
		double memTr=-1;double memTe=-1;
		boolean isUp=false;
		//Identify different Tr, and the number of echoes for each
		for(int c=0;c<trData.length;c++) {
			if(trData[c]>9000) {listTe.add(teData[c]);}
			if(trData[c]!=memTr) {
				memTr=trData[c];
				listTr.add(new Double(memTr));
				nt1++;
			}
		}
		tabTrVals=new double[listTr.size()]; tabTeVals=new double[listTe.size()];tabTrSeriesLength=new int[listTr.size()];
		for(int i=0;i<listTr.size();i++) {tabTrVals[i]=listTr.get(i);}
		for(int i=0;i<listTe.size();i++) {tabTeVals[i]=listTe.get(i);}
		magData2D=new double[listTr.size()][listTe.size()];for(int i=0;i<magData2D.length;i++)for(int j=0;j<magData2D[i].length;j++)magData2D[i][j]=-1;
		for(int c=0;c<trData.length;c++) {
			for(int tr=0;tr<tabTrVals.length;tr++) {
				if(trData[c]==tabTrVals[tr]) {
					for(int te=0;te<tabTeVals.length;te++) {
						if(teData[c]==tabTeVals[te]) {
							magData2D[tr][te]=magData[c];
							tabTrSeriesLength[tr]=te+1;
							if(debugBionano)System.out.println(magData2D[tr][te]);
						}
					}				
				}
			}
		}
		if(debugBionano)System.out.println("Nul ?"+(magData2D==null));
	}	
	
	public void inspect2DTab() {
		if(debugBionano)System.out.println("Nul ?"+(magData2D==null));
		for(int i=0;i<magData2D.length;i++) {
			if(debugBionano)System.out.println("\nInspecting new serie for Tr="+tabTrVals[i]+" of length "+tabTrSeriesLength[i]);
			for(int j=0;j<magData2D[i].length;j++) {
				if(debugBionano)System.out.print("  Te "+tabTeVals[j]+" : "+magData2D[i][j]);
			}
		}		
		if(debugBionano)System.out.println();
	}

	public void inspectT2EstimatedValues() {
		if(debugBionano)System.out.println("T2 vas computed. Values : T2="+this.t2+"  , R2="+this.r2+"  , M0="+this.m0t2);
	}
	
	public void computeT2Bionano() {
		int indexTr10000=magData2D.length-1;
		t2=0;
		m0t2=0;
		r2=Math.log(magData2D[indexTr10000][0]/magData2D[indexTr10000][2])/(tabTeVals[2]-tabTeVals[0]);
		if(debugBionano)System.out.println("Computed r2="+r2);
		if(r2>0.1)r2=Double.NaN;
		else if(r2<0.0006666)r2=Double.NaN;
		else {
			t2=1/r2;
			m0t2=magData2D[indexTr10000][0]*Math.exp(tabTeVals[0]*r2);
			if(debugBionano)System.out.println("Computed m0="+m0t2);
			if(debugBionano)System.out.println("Computed t2="+t2);
		}
	}	
	
	
	public void computeT1Bionano() {
		if(debugBionano)System.out.println("\n\nT1 COMPUTATION");
		//T1 COMPUTATION
    	//Sum all echoes for each Tr
		double[]magSum=new double[magData2D.length-1];
		int nForSums=200;
        for(int i=0;i<magSum.length;i++)if(tabTrSeriesLength[i]<nForSums)nForSums=tabTrSeriesLength[i];
        for(int i=0;i<magSum.length;i++)for(int j=0;(j<nForSums && j<tabTrSeriesLength[j]);j++)magSum[i]+=magData2D[i][j];
		
        double factorT2=0;
        for(int i=0;i<nForSums;i++) {
        	factorT2+=Math.exp(-tabTeVals[i]/t2);
        	if(debugBionano)System.out.println("Constitution factorT2="+factorT2+" apres iteration "+i);
        }
        
    	//Build an ArrayList of cases to compute
        ArrayList<double[]>cases=new ArrayList<double[]>();
        for(int i=0;i<magSum.length;i++) for(int j=i+1;j<magSum.length;j++) {
        	if(debugBionano)System.out.println("Checking couple "+tabTrVals[j]+" / "+tabTrVals[i]);
    		boolean newCase=false;
        	if((tabTrVals[j]/tabTrVals[i])==1.5) {cases.add(new double[] {i,j,15,0,0,0});newCase=true;}
        	if((tabTrVals[j]/tabTrVals[i])==2) {cases.add(new double[] {i,j,20,0,0,0});newCase=true;}//      0       1         2      3    4
        	if((tabTrVals[j]/tabTrVals[i])==3) {cases.add(new double[] {i,j,30,0,0,0});newCase=true;}// Indice a , indice b,   Cas,   M0 , T1,
        	if(newCase)if(debugBionano)System.out.println("i="+i+" Tr="+tabTrVals[i]+"  j="+j+" Tr="+tabTrVals[j]+"  constitued case"+TransformUtils.stringVectorN(cases.get(cases.size()-1),""));
        }
        
        //Process the cases
        for(int i=0;i<cases.size();i++) {
        	double[]tab=cases.get(i);
        	if(debugBionano)System.out.println("\nProcessing case "+i+" : "+TransformUtils.stringVectorN(cases.get(i),""));
        	double Ma=magSum[(int) tab[0]];
    		double Mb=magSum[(int) tab[1]];
    		double Ta=tabTrVals[(int) tab[0]];
    		double Tb=tabTrVals[(int) tab[1]];
    		if(debugBionano)System.out.println("Ma="+Ma+" Ta="+Ta+"  Mb="+Mb+"  Tb="+Tb);
    		try {
	    		if(tab[2]==15) {//TODO : gerer les exceptions, et envoyer une valeur négative
	        		tab[3]=-(2*Ma*Ma*Ma)/(Mb*Mb-3*Ma*Ma+Math.sqrt((-(Ma-Mb)*(Ma-Mb)*(Ma-Mb)*(3*Ma+Mb))));
	        	}
	        	if(tab[2]==20) {
	        		tab[3]=Ma*Ma/(2*Ma-Mb);
	        	}
	        	if(tab[2]==30) {
	        		tab[3]=-(2*Ma*Ma*Ma)/(Math.sqrt(-Ma*Ma*Ma*(3*Ma-4*Mb))-3*Ma*Ma);
	        	}
    		}
    		catch(Exception e) {if(debugBionano)System.out.println("Got an exception in polynom!");tab[3]=Double.NaN;}
    		try {
    			tab[4]=-Tb/(Math.log(1-Mb/tab[3]));        	
    		}
    		catch(Exception e) {if(debugBionano)System.out.println("Got an exception in log!");tab[4]=Double.NaN;}
        	//Remove all values according to conditions
            if(tab[3]<=0)tab[3]=Double.NaN;
            if(tab[4]<=0)tab[4]=Double.NaN;
            double temp=tab[4];
            if(tab[4]>3000) {temp=tab[4];tab[4]=Double.NaN;}
            if(debugBionano)System.out.println("Finally, values="+tab[3]+" , "+temp+" set to "+tab[4]);
        }
        
        
        
    	//Remove all values out from three sigma from mean
        int nM0=0;int nT1=0;
        for(int i=0;i<cases.size();i++) {
        	if(cases.get(i)[3]>0)nM0++;
        	if(cases.get(i)[4]>0)nT1++;
        }
        double[]tabM0=new double[nM0];
        double[]tabT1=new double[nT1];
        nM0=0;nT1=0;
        for(int i=0;i<cases.size();i++) {
        	if(cases.get(i)[3]>0)tabM0[nM0++]=cases.get(i)[3];
        	if(cases.get(i)[4]>0)tabT1[nT1++]=cases.get(i)[4];
        }
 
        t1=meanHampel(tabT1);
        double m0t1Temp=meanHampel(tabM0);

        m0t1=m0t1Temp/factorT2;
        if(debugBionano) System.out.println("Finally, got : t1="+t1+"  m0t1"+m0t1);
        bionanoParams=new double[] {m0t1,t1,m0t2,t2,0};
    }
    
	
	public double meanHampel(double[]tab) {
		if(tab.length==1)return tab[0];
		//copy in tab with no nan
		int nGood=0;
		if(debugBionano)System.out.println("Hampel");
		if(debugBionano)System.out.println(TransformUtils.stringVectorN(tab, "tabInit"));
		for(int i=0;i<tab.length;i++)if(!Double.isNaN(tab[i]))nGood++;
		if(debugBionano)System.out.println("Ngood="+nGood);
		double[]tab2=new double[nGood];nGood=0;
		for(int i=0;i<tab.length;i++)if(!Double.isNaN(tab[i]))tab2[nGood++]=tab[i];
		if(debugBionano)System.out.println(TransformUtils.stringVectorN(tab, "tab2"));
		
		//compute mean and std
		double[]stats=VitimageUtils.statistics1D(tab2);
		if(debugBionano)System.out.println(TransformUtils.stringVectorN(stats, "Stats tab2"));

		//copy in tab with no outliers
		nGood=0;
		for(int i=0;i<tab2.length;i++)if( (tab2[i]<(stats[0]+3*stats[1])) && (tab2[i]>(stats[0]-3*stats[1])) ) nGood++;
		double[]tabFinal=new double[nGood];nGood=0;
		if(debugBionano)System.out.println("Ngood="+nGood);
		for(int i=0;i<tab2.length;i++)if( (tab2[i]<(stats[0]+3*stats[1])) && (tab2[i]>(stats[0]-3*stats[1])) ) tabFinal[nGood++]=tab2[i];
		if(debugBionano)System.out.println(TransformUtils.stringVectorN(tabFinal, "tabFinal"));
			
		//compute mean and std
		double mean=VitimageUtils.statistics1D(tabFinal)[0];
		if(debugBionano)System.out.println("Mean ="+mean+" . Return.");
		return mean;
	}
	
	
	
    public void config(int maxIterations) {
    	IterFactor=maxIterations;
    }
    
    public void doFit() {
    	doFit(fit);
    }
    
    public void doFit(int fitType) {
    	if(fitType==MRUtils.T1T2_BIONANO) {computeBioNano();return;}
    	
    	doFit(fitType, false);
    }
    
    public void doFit(int fitType, boolean showSettings) {

        int saveFitType = fitType;
        fit = fitType;
       
 		if (initialParams!=null) {
			for (int i=0; i<numParams; i++)
				simp[0][i] = initialParams[i];
			initialParams = null;
		}
         restart(0);
        
        numIter = 0;
        boolean done = false;
        double[] center = new double[numParams];  // mean of simplex vertices
        while (!done) {
            numIter++;
            for (int i = 0; i < numParams; i++) center[i] = 0.0;
            // get mean "center" of vertices, excluding worst
            for (int i = 0; i < numVertices; i++)
                if (i != worst)
                    for (int j = 0; j < numParams; j++)
                        center[j] += simp[i][j];
            // Reflect worst vertex through centre
            for (int i = 0; i < numParams; i++) {
                center[i] /= numParams;
                next[i] = center[i] + alpha*(simp[worst][i] - center[i]);
            }
            sumResiduals(next);
            // if it's better than the best...
            if (next[numParams] <= simp[best][numParams]) {
                newVertex();
                // try expanding it
                for (int i = 0; i < numParams; i++)
                    next[i] = center[i] + gamma * (simp[worst][i] - center[i]);
                sumResiduals(next);
                // if this is even better, keep it
                if (next[numParams] <= simp[worst][numParams])
                    newVertex();
            }
            // else if better than the 2nd worst keep it...
            else if (next[numParams] <= simp[nextWorst][numParams]) {
                newVertex();
            }
            // else try to make positive contraction of the worst
            else {
                for (int i = 0; i < numParams; i++)
                    next[i] = center[i] + beta*(simp[worst][i] - center[i]);
                sumResiduals(next);
                // if this is better than the second worst, keep it.
                if (next[numParams] <= simp[nextWorst][numParams]) {
                    newVertex();
                }
                // if all else fails, contract simplex in on best
                else {
                    for (int i = 0; i < numVertices; i++) {
                        if (i != best) {
                            for (int j = 0; j < numVertices; j++)
                                simp[i][j] = beta*(simp[i][j]+simp[best][j]);
                            sumResiduals(simp[i]);
                        }
                    }
                }
            }
            order();
            
            double rtol = 2 * Math.abs(simp[best][numParams] - simp[worst][numParams]) /
            (Math.abs(simp[best][numParams]) + Math.abs(simp[worst][numParams]) + 0.0000000001);
            
            if (numIter >= maxIter)
            	done = true;
            else if (rtol < maxError) {
                restarts--;
                if (restarts < 0)
                    done = true;
                else
                    restart(best);
             }
        }
        fitType = saveFitType;
    }
        


    /** Initialise the simplex */
    protected void initialize(int fitType) {
        // Calculate some things that might be useful for predicting parametres
        numParams = getNumParams(fitType);
        numVertices = numParams + 1;      // need 1 more vertice than parametres,
        simp = new double[numVertices][numVertices];
        next = new double[numVertices];        
        maxIter = IterFactor * numParams * numParams;  // Where does this estimate come from?
        if(fitType==MRUtils.T1T2_MULTIMULTI_RICE) maxIter/=numParams;
        restarts = defaultRestarts;
        nRestarts = 0;
        switch (fit) {
           case MRUtils.T1T2_MONO_RICE:
        	   simp[0][0] = VitimageUtils.max(magData);
        	   simp[0][1] = 2000;//T1
        	   simp[0][2] = 50;//T2
               break;
           case MRUtils.T1T2_MULTI_RICE:
        	   simp[0][0] = VitimageUtils.max(magData)/2*1.1;
        	   simp[0][1] = VitimageUtils.max(magData)/2*0.9;
        	   simp[0][2] = 1000;//T1
        	   simp[0][3] = 20;//T2
        	   simp[0][4] = 100;//T2
               break;
           case MRUtils.T1T2_MULTIMULTI_RICE:
        	   simp[0][0] = VitimageUtils.max(magData)/2;
        	   simp[0][1] = VitimageUtils.max(magData)/2;
        	   simp[0][2] = 500;//T1
        	   simp[0][3] = 3000;//T1
        	   simp[0][4] = 20;//T2
        	   simp[0][5] = 100;//T2
               break;
               default:break;
        }
    }
        
    /** Restart the simplex at the nth vertex */
    void restart(int n) {
        // Copy nth vertice of simplex to first vertice
        for (int i = 0; i < numParams; i++) {
            simp[0][i] = simp[n][i];
        }
        sumResiduals(simp[0]);          // Get sum of residuals^2 for first vertex
        double[] step = new double[numParams];
        for (int i = 0; i < numParams; i++) {
            step[i] = simp[0][i] / 2.0;     // Step half the parametre value
            if (step[i] == 0.0)             // We can't have them all the same or we're going nowhere
                step[i] = 0.01;
        }
        // Some kind of factor for generating new vertices
        double[] p = new double[numParams];
        double[] q = new double[numParams];
        for (int i = 0; i < numParams; i++) {
            p[i] = step[i] * (Math.sqrt(numVertices) + numParams - 1.0)/(numParams * root2);
            q[i] = step[i] * (Math.sqrt(numVertices) - 1.0)/(numParams * root2);
        }
        // Create the other simplex vertices by modifing previous one.
        for (int i = 1; i < numVertices; i++) {
            for (int j = 0; j < numParams; j++) {
                simp[i][j] = simp[i-1][j] + q[j];
            }
            simp[i][i-1] = simp[i][i-1] + p[i-1];
            sumResiduals(simp[i]);
        }
        // Initialise current lowest/highest parametre estimates to simplex 1
        best = 0;
        worst = 0;
        nextWorst = 0;
        order();
        nRestarts++;
    }
        
    // Display simplex [Iteration: s0(p1, p2....), s1(),....] in Log window
    void showSimplex(int iter) {
        ij.IJ.log("" + iter);
        for (int i = 0; i < numVertices; i++) {
            String s = "";
            for (int j=0; j < numVertices; j++)
                s += "  "+ ij.IJ.d2s(simp[i][j], 6);
            ij.IJ.log(s);
        }
    }
        
    /** Get number of parameters for current fit formula */
    public int getNumParams(int fitType) {
        switch (fitType) {
			case MRUtils.T1T2_MONO_RICE: return 3;
			case MRUtils.T1T2_MULTI_RICE: return 5;
			case MRUtils.T1T2_MULTIMULTI_RICE: return 6;
        }
        return 0;
    }
        
	/** Returns formula value for parameters 'p' at 'x' */
	public double f(double[] p, double tr,double te) {
		return f(fit, p, tr,te);
	}

   /** Returns 'fit' formula value for parameters "p" at "x" */
    public double f(int fit, double[] p, double tr,double te) {
    	double y;
        switch (fit) {
        	case MRUtils.T1T2_MONO_RICE:return besFunkCost(p[0]*(1 - Math.exp(-(tr / p[1])))*Math.exp(-(te / p[2]) ) ,sigma);
        	case MRUtils.T1T2_MULTI_RICE:return besFunkCost(  
        			(        ( (1 - Math.exp(-(tr / p[2]))) * (  ( p[0] * Math.exp(-(te / p[3]) ) )   +  ( p[1] * Math.exp(-(te / p[4]) ) ) )    )  )  ,sigma);
        	case MRUtils.T1T2_MULTIMULTI_RICE:return besFunkCost(  
        			(        ( (1 - Math.exp(-(tr / p[2]))) * ( p[0] * Math.exp(-(te / p[4]) ) ) )  +
        					 ( (1 - Math.exp(-(tr / p[3])))*  ( p[1] * Math.exp(-(te / p[5]) ) ) )      )  ,sigma);
             default:
                return 0.0;
        }
    }
    
    /** Get the set of parameter values from the best corner of the simplex */
    public double[] getParams() {
        if(fit==MRUtils.T1T2_BIONANO)return bionanoParams;
    	order();
        return simp[best];
    }
    
	/** Returns residuals array ie. differences between data and curve. */
	public double[] getResiduals() {
		int saveFit = fit;
		double[] params = getParams();
		double[] residuals = new double[numPoints];
		for (int i=0; i<numPoints; i++)
			residuals[i] = magData[i] - f(fit, params, trData[i],teData[i]);
		fit = saveFit;
		return residuals;
	}
    
    /* Last "parametre" at each vertex of simplex is sum of residuals
     * for the curve described by that vertex
     */
    public double getSumResidualsSqr() {
        double sumResidualsSqr = (getParams())[getNumParams(fit)];
        return sumResidualsSqr;
    }
    
    /**  Returns the standard deviation of the residuals. */
    public double getSD() {
    	double[] residuals = getResiduals();
		int n = residuals.length;
		double sum=0.0, sum2=0.0;
		for (int i=0; i<n; i++) {
			sum += residuals[i];
			sum2 += residuals[i]*residuals[i];
		}
		double stdDev = (n*sum2-sum*sum)/n;
		return Math.sqrt(stdDev/(n-1.0));
    }
    
    /** Returns R^2, where 1.0 is best.
    <pre>
     r^2 = 1 - SSE/SSD
     
     where:	 SSE = sum of the squares of the errors
                 SSD = sum of the squares of the deviations about the mean.
    </pre>
    */
    public double getRSquared() {
        double sumY = 0.0;
        for (int i=0; i<numPoints; i++) sumY += magData[i];
        double mean = sumY/numPoints;
        double sumMeanDiffSqr = 0.0;
        for (int i=0; i<numPoints; i++)
            sumMeanDiffSqr += sqr(magData[i]-mean);
        double rSquared = 0.0;
        if (sumMeanDiffSqr>0.0)
            rSquared = 1.0 - getSumResidualsSqr()/sumMeanDiffSqr;
        return rSquared;
    }

    /**  Get a measure of "goodness of fit" where 1.0 is best. */
    public double getFitGoodness() {
        double sumY = 0.0;
        for (int i = 0; i < numPoints; i++) sumY += magData[i];
        double mean = sumY / numPoints;
        double sumMeanDiffSqr = 0.0;
        int degreesOfFreedom = numPoints - getNumParams(fit);
        double fitGoodness = 0.0;
        for (int i = 0; i < numPoints; i++) {
            sumMeanDiffSqr += sqr(magData[i] - mean);
        }
        if (sumMeanDiffSqr > 0.0 && degreesOfFreedom != 0)
            fitGoodness = 1.0 - (getSumResidualsSqr() / degreesOfFreedom) * ((numPoints) / sumMeanDiffSqr);
        
        return fitGoodness;
    }
    
    /** Get a string description of the curve fitting results
     * for easy output.
     */
        
    double sqr(double d) { return d * d; }
    
	/** Adds sum of square of residuals to end of array of parameters */
	void sumResiduals (double[] x) {
		x[numParams] = 0.0;
		for (int i=0; i<numPoints; i++) {
			x[numParams] = x[numParams] + sqr(f(fit,x,trData[i],teData[i])-magData[i]);
		}
	}

    /** Keep the "next" vertex */
    void newVertex() {
    	for (int i = 0; i < numVertices; i++)
            simp[worst][i] = next[i];
    }
    
    /** Find the worst, nextWorst and best current set of parameter estimates */
    void order() {
        for (int i = 0; i < numVertices; i++) {
            if (simp[i][numParams] < simp[best][numParams])	best = i;
            if (simp[i][numParams] > simp[worst][numParams]) worst = i;
        }
        nextWorst = best;
        for (int i = 0; i < numVertices; i++) {
            if (i != worst) {
                if (simp[i][numParams] > simp[nextWorst][numParams]) nextWorst = i;
            }
        }
        //        IJ.log("B: " + simp[best][numParams] + " 2ndW: " + simp[nextWorst][numParams] + " W: " + simp[worst][numParams]);
    }

    /** Get number of iterations performed */
    public int getIterations() {
        return numIter;
    }
    
    /** Get maximum number of iterations allowed */
    public int getMaxIterations() {
        return maxIter;
    }
    
    /** Set maximum number of iterations allowed */
    public void setMaxIterations(int x) {
        maxIter = x;
    }
    
    /** Get number of simplex restarts to do */
    public int getRestarts() {
        return defaultRestarts;
    }
    
    /** Set number of simplex restarts to do */
    public void setRestarts(int n) {
        defaultRestarts = n;
    }

	/** Sets the initial parameters, which override the default initial parameters. */
	public void setInitialParameters(double[] params) {
		initialParams = params;
	}

    /**
     * Gets index of highest value in an array.
     * 
     * @param              Double array.
     * @return             Index of highest value.
     */
    public static int getMax(double[] array) {
        double max = array[0];
        int index = 0;
        for(int i = 1; i < array.length; i++) {
            if(max < array[i]) {
            	max = array[i];
            	index = i;
            }
        }
        return index;
    }
    
	public double[] getTePoints() {
		return teData;
	}
	
	public double[] getTrPoints() {
		return trData;
	}

	public double[] getMagPoints() {
		return magData;
	}
	
	public int getFit() {
		return fit;
	}

	
}
