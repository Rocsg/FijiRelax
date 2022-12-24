/*
 * 
 */
package io.github.rocsg.fijirelax.curvefit;

import java.util.ArrayList;

import io.github.rocsg.fijiyama.registration.TransformUtils;
import io.github.rocsg.fijiyama.common.VitimageUtils;
import ij.macro.Interpreter;
import io.github.rocsg.fijirelax.mrialgo.MRUtils;


/**
 * A simplex-based solution to curve fitting of exponential functions over MRI observations points
 * This one is preferred as default solution, as it produce way less outliers than Levenberg implementation, is faster, and converges toward the same values in most cases
 *  	Copyright (C) 2022  io.github.rocsg
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>
*/
public class SimplexDualCurveFitterNoBias{
	
	/** The iterations break. */
	public boolean iterationsBreak=false;
	
	/** The debug. */
	private boolean debug=true;
	
	/** The estimate delta te. */
	public boolean estimateDeltaTe=false;
	
	/** The Constant debugBionano. */
	public static final boolean debugBionano=false;
	
	/** The iter factor. */
	public static  int iterFactor = 1000;
	
	/** The sigma. */
	protected double sigma;   
	
	/** The Constant alpha. */
	private static final double alpha = -1.0;	  // reflection coefficient
	
	/** The Constant beta. */
	private static final double beta = 0.5;	  // contraction coefficient
	
	/** The Constant gamma. */
	private static final double gamma = 2.0;	  // expansion coefficient
	
	/** The Constant root2. */
	private static final double root2 = 1.414214; // square root of 2
	
	/** The fit. */
	protected int fit;                // Number of curve type to fit
	
	/** The mag data. */
	protected double[] trData,teData, magData;  // x,y data to fit
	
	/** The num points. */
	protected int numPoints;          // number of data points
	
	/** The num params. */
	protected int numParams;          // number of parametres
	
	/** The num vertices. */
	protected int numVertices;        // numParams+1 (includes sumLocalResiduaalsSqrd)
	
	/** The worst. */
	private int worst;			// worst current parametre estimates
	
	/** The next worst. */
	private int nextWorst;		// 2nd worst current parametre estimates
	
	/** The best. */
	private int best;			// best current parametre estimates
	
	/** The simp. */
	protected double[][] simp; 		// the simplex (the last element of the array at each vertice is the sum of the square of the residuals)
	
	/** The next. */
	protected double[] next;		// new vertex to be tested
	
	/** The num iter. */
	private int numIter;		// number of iterations so far
	
	/** The max iter. */
	private int maxIter; 	// maximum number of iterations per restart
	
	/** The restarts. */
	private int restarts; 	// number of times to restart simplex after first soln.
	
	/** The default restarts. */
	private static int defaultRestarts = 2;  // default number of restarts
	
	/** The n restarts. */
	private int nRestarts;  // the number of restarts that occurred
	
	/** The max error. */
	private static double maxError = 1e-5;    // maximum error tolerance
	
	/** The initial params. */
	private double[] initialParams;  // user specified initial parameters
	
	/** The time. */
	private long time;  //elapsed time in ms
	
	/** The custom formula. */
	private String customFormula;
	
	/** The custom param count. */
	private static int customParamCount;
	
	/** The initial values. */
	private double[] initialValues;
	
	/** The macro. */
	private Interpreter macro;
	
	/** The bionano params. */
	private double[] bionanoParams;
	
	/** The mag data 2 D. */
	private double[][]magData2D;
	
	/** The tab tr vals. */
	private double[]tabTrVals;
	
	/** The tab te vals. */
	private double[]tabTeVals;
	
	/** The tab tr series length. */
	private int []tabTrSeriesLength;
	
	/** The t 2. */
	private double t2;
	
	/** The m 0 t 2. */
	private double m0t2;
	
	/** The r 2. */
	private double r2;
	
	/** The bionano factor. */
	private double bionanoFactor=1+Math.exp(-0.25)+Math.exp(-0.5)+Math.exp(-0.75);
	
	/** The t 1. */
	private double t1;
	
	/** The m 0 t 1. */
	private double m0t1;
	
	/** The parameters boundaries. */
	private double[][] parametersBoundaries;
	
	/** The min delta chi 2. When variation of chi 2 will be lower than it, induce an early break of fit*/
	private double minDeltaChi2;

    /**
     *  Construct a new SimplexCurveFitter.
     *
     * @param trData the recovery times
     * @param teData the echo times
     * @param magData the magnitude data
     * @param fitType the type of fit, chosen among @link io.github.rocsg.fijirelax.mrialgo.MRDataType.Class
     * @param sigma the sigma of Rice noise
     * @param debu a debug flag for getting more verbose output
     */
    public SimplexDualCurveFitterNoBias (double[] trData, double[]teData,double[] magData, int fitType,double sigma,boolean debu) {
    	this.debug=debu;
    	this.sigma=sigma;
        this.trData = trData;
        this.teData = teData;
        this.magData = magData;
        numPoints = trData.length;
        this.fit=fitType;
        initialize(fit);
    }
    
    /**
     * Sets the max number of iterations before early break
     *
     * @param maxIterations the new iter factor
     */
    public void setIterFactor(int maxIterations) {
    	iterFactor=maxIterations;
    }
    
    /**
     * Start the fit with the defined setup in the current object
     */
    public void doFit() {
    	doFit(fit);
    }
    
    /**
     * Starts the fit with the defined setup in the current object, and the given fit type
     *
     * @param fitType the fit type, chosen among @link io.github.rocsg.fijirelax.mrialgo.MRDataType.Class
     */
    public void doFit(int fitType) {
    	if(fitType==MRUtils.T1T2_BIONANO) {computeBioNano();return;}  	
    	doFit(fitType, false);
    }
    
    /**
     * Starts the fit with the defined setup in the current object, and the given fit type
     *
     * @param fitType the fit type, chosen among @link io.github.rocsg.fijirelax.mrialgo.MRDataType.Class
     * @param showSettings a flag (without effect since v4.0.2, deprecated soon)
     */
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
        	if(debug)System.out.println("Iteration "+numIter);
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
                if(next[i]<parametersBoundaries[i][0])next[i]=parametersBoundaries[i][0];
                if(next[i]>parametersBoundaries[i][1])next[i]=parametersBoundaries[i][1];
            }
            sumResiduals(next);
            // if it's better than the best...
            if (next[numParams] <= simp[best][numParams]) {
                newVertex();
                // try expanding it
                for (int i = 0; i < numParams; i++) {
                    next[i] = center[i] + gamma * (simp[worst][i] - center[i]);
	                if(next[i]<parametersBoundaries[i][0])next[i]=parametersBoundaries[i][0];
	                if(next[i]>parametersBoundaries[i][1])next[i]=parametersBoundaries[i][1];
                }
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
                for (int i = 0; i < numParams; i++) {
                    next[i] = center[i] + beta*(simp[worst][i] - center[i]);
                    if(next[i]<parametersBoundaries[i][0])next[i]=parametersBoundaries[i][0];
                    if(next[i]>parametersBoundaries[i][1])next[i]=parametersBoundaries[i][1];
                }
                sumResiduals(next);
                // if this is better than the second worst, keep it.
                if (next[numParams] <= simp[nextWorst][numParams]) {
                    newVertex();
                }
                // if all else fails, contract simplex in on best
                else {
                    for (int i = 0; i < numVertices; i++) {
                        if (i != best) {
                            for (int j = 0; j < numVertices; j++) {
                                simp[i][j] = beta*(simp[i][j]+simp[best][j]);
                                if(simp[i][j]<parametersBoundaries[j][0])simp[i][j]=parametersBoundaries[j][0];
                                if(simp[i][j]>parametersBoundaries[j][1])simp[i][j]=parametersBoundaries[j][1];
                            }
                            sumResiduals(simp[i]);
                        }
                    }
                }
                if(debug)System.out.println(TransformUtils.stringVectorN(simp[best], "Best vertex : "));                
            }
            order();
    
    		double rtol = sigma<=0 ? (Math.abs(simp[best][numParams] - simp[worst][numParams]))/((numPoints-numParams))  : (Math.abs(simp[best][numParams] - simp[worst][numParams]))/(sigma*sigma*(numPoints-numParams));
            if(debug)System.out.println("Delta khi2="+rtol+"\n");
            if (numIter >= maxIter) {
            	iterationsBreak=true;
            	done = true;
            }
            else if (rtol < MRUtils.minDeltaChi2) {
                restarts--;
                if (restarts < 0)
                    done = true;
                else
                    restart(best);
             }
        }
        fitType = saveFitType;
    }
        
    /**
     * Gets the default parameters boundaries i.e. the acceptable boundaries for the exponential curve to estimate.
     *
     * @param fitType the fit type, chosen among @link io.github.rocsg.fijirelax.mrialgo.MRDataType.Class
     * @return the default parameters boundaries
     */
    public double[][] getDefaultParametersBoundaries(int fitType) {
    	double[][]parametersBoundaries=new double[MRUtils.getNparams(fitType)+1][2];
    	for(int i=0;i<parametersBoundaries.length;i++)parametersBoundaries[i]=new double[] {-MRUtils.infinity,MRUtils.infinity};
    	return parametersBoundaries;
    }
    
    /**
     * Sets the min delta chi 2.
     *
     * @param delta the new min delta chi 2
     */
    public void setMinDeltaChi2(double delta) {minDeltaChi2 = delta;}


    /**
     *  Initialise the simplex.
     *
     * @param fitType the fit type, chosen among @link io.github.rocsg.fijirelax.mrialgo.MRDataType.Class
     */
    protected void initialize(int fitType) {
        // Calculate some things that might be useful for predicting parametres
        numParams = MRUtils.getNparams(fitType);
        numVertices = numParams + 1;      // need 1 more vertice than parametres,
		parametersBoundaries=new double[numVertices][2];
		parametersBoundaries[numVertices-1]=new double[] {-MRUtils.infinity,MRUtils.infinity};
        simp = new double[numVertices][numVertices];
        next = new double[numVertices];        
        double minTr = VitimageUtils.min(trData);
        double maxTr = VitimageUtils.max(trData);
        double medTr=minTr/2+maxTr/2;
        double minTe = VitimageUtils.min(teData);
        double maxTe = VitimageUtils.max(teData);
        double medTe=minTe/2+maxTe/2;
        double maxMag=VitimageUtils.max(magData);
        double medMag=maxMag/2;
        double minMag=VitimageUtils.min(magData);
        if(magData.length<2)minMag= magData[magData.length-1];
        else if(magData.length<3)minMag= (magData[magData.length-1]+magData[magData.length-2])/2;
        else minMag= (magData[magData.length-1]+magData[magData.length-2]+magData[magData.length-3])/3;
        maxIter = iterFactor * numParams;
        restarts = defaultRestarts;
        nRestarts = 0;
        switch (fit) {

        case MRUtils.T1_MONO:
     	   simp[0][0] = maxMag;
     	   simp[0][1] = medTr;
           if(MRUtils.useBoundaries) {
					parametersBoundaries[0]=new double[] {MRUtils.epsilon,maxMag*MRUtils.factorT1M0MaxRatio};
					parametersBoundaries[1]=new double[] {minTr*MRUtils.factorT1MinRatio,maxTr*MRUtils.factorT1MaxRatio};
           }
		   else parametersBoundaries=getDefaultParametersBoundaries(fitType);
           break;
        case MRUtils.T1_MONO_BIAS:
     	   simp[0][0] = maxMag;
     	   simp[0][1] = medTr;
     	   simp[0][2] = 0;
            if(MRUtils.useBoundaries) {
					parametersBoundaries[0]=new double[] {MRUtils.epsilon,maxMag*MRUtils.factorT1M0MaxRatio};
					parametersBoundaries[1]=new double[] {minTr*MRUtils.factorT1MinRatio,maxTr*MRUtils.factorT1MaxRatio};
					parametersBoundaries[2]=new double[] {0,maxMag};
            }
			else parametersBoundaries=getDefaultParametersBoundaries(fitType);
            break;
        case MRUtils.T1_MONO_RICE:
     	   simp[0][0] = maxMag;
     	   simp[0][1] = medTr;
            if(MRUtils.useBoundaries) {
					parametersBoundaries[0]=new double[] {MRUtils.epsilon,maxMag*MRUtils.factorT1M0MaxRatio};
					parametersBoundaries[1]=new double[] {minTr*MRUtils.factorT1MinRatio,maxTr*MRUtils.factorT1MaxRatio};
            }
			else parametersBoundaries=getDefaultParametersBoundaries(fitType);
            break;


            
       case MRUtils.T2_MONO:
      	   simp[0][0] = maxMag;
      	   simp[0][1] = medTe;
             if(MRUtils.useBoundaries) {
 					parametersBoundaries[0]=new double[] {MRUtils.epsilon,maxMag*MRUtils.factorT2M0MaxRatio};
 					parametersBoundaries[1]=new double[] {minTe*MRUtils.factorT2MinRatio,maxTe*MRUtils.factorT2MaxRatio};
             }
 			 else parametersBoundaries=getDefaultParametersBoundaries(fitType);
             break;
       case MRUtils.T2_MONO_BIAS:
     	   simp[0][0] = maxMag;
     	   simp[0][1] = medTe;
     	   simp[0][2] = minMag;
            if(MRUtils.useBoundaries) {
					parametersBoundaries[0]=new double[] {MRUtils.epsilon,maxMag*MRUtils.factorT2M0MaxRatio};
					parametersBoundaries[1]=new double[] {minTe*MRUtils.factorT2MinRatio,maxTe*MRUtils.factorT2MaxRatio};
					parametersBoundaries[2]=new double[] {0,maxMag};
            }
			   else parametersBoundaries=getDefaultParametersBoundaries(fitType);
            break;
       case MRUtils.T2_MONO_RICE:
     	   simp[0][0] = maxMag;
     	   simp[0][1] = medTe;
            if(MRUtils.useBoundaries) {
					parametersBoundaries[0]=new double[] {MRUtils.epsilon,maxMag*MRUtils.factorT2M0MaxRatio};
					parametersBoundaries[1]=new double[] {minTe*MRUtils.factorT2MinRatio,maxTe*MRUtils.factorT2MaxRatio};
            }
			else parametersBoundaries=getDefaultParametersBoundaries(fitType);
            break;

            
            
       case MRUtils.T2_MULTI:
    	   simp[0][0] = medMag*1.6;
    	   simp[0][1] = maxMag*0.6;
    	   simp[0][2] = medTe*0.4;
    	   simp[0][3] = medTe*1.3;
           if(MRUtils.useBoundaries) {
				parametersBoundaries[0]=new double[] {MRUtils.epsilon,maxMag*MRUtils.factorT2M0MaxRatio};
				parametersBoundaries[1]=new double[] {MRUtils.epsilon,maxMag*MRUtils.factorT2M0MaxRatio};
				parametersBoundaries[2]=new double[] {minTe*MRUtils.factorT2MinRatio,maxTe*MRUtils.factorT2MaxRatio};
				parametersBoundaries[3]=new double[] {minTe*MRUtils.factorT2MinRatio,maxTe*MRUtils.factorT2MaxRatio};
           }
		   else parametersBoundaries=getDefaultParametersBoundaries(fitType);
           break;
       case MRUtils.T2_MULTI_BIAS:
    	   simp[0][0] = medMag*1.6;
    	   simp[0][1] = maxMag*0.6;
    	   simp[0][2] = medTe*0.4;
    	   simp[0][3] = medTe*1.3;
    	   simp[0][4] = minMag;
           if(MRUtils.useBoundaries) {
				parametersBoundaries[0]=new double[] {MRUtils.epsilon,maxMag*MRUtils.factorT2M0MaxRatio};
				parametersBoundaries[1]=new double[] {MRUtils.epsilon,maxMag*MRUtils.factorT2M0MaxRatio};
				parametersBoundaries[2]=new double[] {minTe*MRUtils.factorT2MinRatio,maxTe*MRUtils.factorT2MaxRatio};
				parametersBoundaries[3]=new double[] {minTe*MRUtils.factorT2MinRatio,maxTe*MRUtils.factorT2MaxRatio};
				parametersBoundaries[4]=new double[] {0,maxMag};
           }
		   else parametersBoundaries=getDefaultParametersBoundaries(fitType);
           break;
       case MRUtils.T2_MULTI_RICE:
    	   simp[0][0] = medMag*1.6;
    	   simp[0][1] = maxMag*0.6;
    	   simp[0][2] = medTe*0.4;
    	   simp[0][3] = medTe*1.3;
           if(MRUtils.useBoundaries) {
				parametersBoundaries[0]=new double[] {MRUtils.epsilon,maxMag*MRUtils.factorT2M0MaxRatio};
				parametersBoundaries[1]=new double[] {MRUtils.epsilon,maxMag*MRUtils.factorT2M0MaxRatio};
				parametersBoundaries[2]=new double[] {minTe*MRUtils.factorT2MinRatio,maxTe*MRUtils.factorT2MaxRatio};
				parametersBoundaries[3]=new double[] {minTe*MRUtils.factorT2MinRatio,maxTe*MRUtils.factorT2MaxRatio};
           }
		   else parametersBoundaries=getDefaultParametersBoundaries(fitType);
           break;
            
        
        case MRUtils.T1T2_MONO:
        	   simp[0][0] = maxMag;
        	   simp[0][1] = medTr;
        	   simp[0][2] = medTe;
               if(MRUtils.useBoundaries) {
					parametersBoundaries[0]=new double[] {MRUtils.epsilon,maxMag*MRUtils.factorT2M0MaxRatio};
					parametersBoundaries[1]=new double[] {minTr*MRUtils.factorT1MinRatio,maxTr*MRUtils.factorT1MaxRatio};
					parametersBoundaries[2]=new double[] {minTe*MRUtils.factorT2MinRatio,maxTe*MRUtils.factorT2MaxRatio};
               }
			   else parametersBoundaries=getDefaultParametersBoundaries(fitType);
               break;
        case MRUtils.T1T2_MONO_BIAS:
     	   simp[0][0] = maxMag;
     	   simp[0][1] = medTr;
     	   simp[0][2] = medTe;
     	   simp[0][3] = minMag;
            if(MRUtils.useBoundaries) {
					parametersBoundaries[0]=new double[] {MRUtils.epsilon,maxMag*MRUtils.factorT2M0MaxRatio};
					parametersBoundaries[1]=new double[] {minTr*MRUtils.factorT1MinRatio,maxTr*MRUtils.factorT1MaxRatio};
					parametersBoundaries[2]=new double[] {minTe*MRUtils.factorT2MinRatio,maxTe*MRUtils.factorT2MaxRatio};
					parametersBoundaries[2]=new double[] {0,maxMag};
            }
			   else parametersBoundaries=getDefaultParametersBoundaries(fitType);
            break;
        case MRUtils.T1T2_MONO_RICE:
     	   simp[0][0] = maxMag;
     	   simp[0][1] = medTr;
     	   simp[0][2] = medTe;
            if(MRUtils.useBoundaries) {
					parametersBoundaries[0]=new double[] {MRUtils.epsilon,maxMag*MRUtils.factorT2M0MaxRatio};
					parametersBoundaries[1]=new double[] {minTr*MRUtils.factorT1MinRatio,maxTr*MRUtils.factorT1MaxRatio};
					parametersBoundaries[2]=new double[] {minTe*MRUtils.factorT2MinRatio,maxTe*MRUtils.factorT2MaxRatio};
            }
			   else parametersBoundaries=getDefaultParametersBoundaries(fitType);
            break;

        
        
        
        case MRUtils.T1T2_MULTI:
        	   simp[0][0] = medMag*1.6;
        	   simp[0][1] = medMag*0.6;
        	   simp[0][2] = medTr;//T1
        	   simp[0][3] = medTe*0.4;//T2
        	   simp[0][4] = medTe*1.3;//T2
               if(MRUtils.useBoundaries) {
					parametersBoundaries[0]=new double[] {MRUtils.epsilon,maxMag*MRUtils.factorT2M0MaxRatio};
					parametersBoundaries[1]=new double[] {MRUtils.epsilon,maxMag*MRUtils.factorT2M0MaxRatio};
					parametersBoundaries[2]=new double[] {minTr*MRUtils.factorT1MinRatio,maxTr*MRUtils.factorT1MaxRatio};
					parametersBoundaries[3]=new double[] {minTe*MRUtils.factorT2MinRatio,maxTe*MRUtils.factorT2MaxRatio};
					parametersBoundaries[4]=new double[] {minTe*MRUtils.factorT2MinRatio,maxTe*MRUtils.factorT2MaxRatio};
               }
			   else parametersBoundaries=getDefaultParametersBoundaries(fitType);
               break;
        case MRUtils.T1T2_MULTI_BIAS:
     	   simp[0][0] = medMag*1.6;
     	   simp[0][1] = medMag*0.6;
     	   simp[0][2] = medTr;//T1
     	   simp[0][3] = medTe*0.4;//T2
     	   simp[0][4] = medTe*1.3;//T2
     	   simp[0][5] = minMag;//Bias
            if(MRUtils.useBoundaries) {
					parametersBoundaries[0]=new double[] {MRUtils.epsilon,maxMag*MRUtils.factorT2M0MaxRatio};
					parametersBoundaries[1]=new double[] {MRUtils.epsilon,maxMag*MRUtils.factorT2M0MaxRatio};
					parametersBoundaries[2]=new double[] {minTr*MRUtils.factorT1MinRatio,maxTr*MRUtils.factorT1MaxRatio};
					parametersBoundaries[3]=new double[] {minTe*MRUtils.factorT2MinRatio,maxTe*MRUtils.factorT2MaxRatio};
					parametersBoundaries[4]=new double[] {minTe*MRUtils.factorT2MinRatio,maxTe*MRUtils.factorT2MaxRatio};
					parametersBoundaries[5]=new double[] {0,maxMag};
            }
			   else parametersBoundaries=getDefaultParametersBoundaries(fitType);
            break;
        case MRUtils.T1T2_MULTI_RICE:
     	   simp[0][0] = medMag*1.6;
     	   simp[0][1] = medMag*0.6;
     	   simp[0][2] = medTr;//T1
     	   simp[0][3] = medTe*0.4;//T2
     	   simp[0][4] = medTe*1.3;//T2
            if(MRUtils.useBoundaries) {
					parametersBoundaries[0]=new double[] {MRUtils.epsilon,maxMag*MRUtils.factorT2M0MaxRatio};
					parametersBoundaries[1]=new double[] {MRUtils.epsilon,maxMag*MRUtils.factorT2M0MaxRatio};
					parametersBoundaries[2]=new double[] {minTr*MRUtils.factorT1MinRatio,maxTr*MRUtils.factorT1MaxRatio};
					parametersBoundaries[3]=new double[] {minTe*MRUtils.factorT2MinRatio,maxTe*MRUtils.factorT2MaxRatio};
					parametersBoundaries[4]=new double[] {minTe*MRUtils.factorT2MinRatio,maxTe*MRUtils.factorT2MaxRatio};
            }
			   else parametersBoundaries=getDefaultParametersBoundaries(fitType);
            break;

        
        case MRUtils.T1T2_DEFAULT_T2_MONO_RICE:
        	   simp[0][0] = maxMag;
        	   simp[0][1] = medTr;//T1
        	   simp[0][2] = medTe;//T2
        	   simp[0][3] = 2*minTe;//DeltaT2
               if(MRUtils.useBoundaries) {
					parametersBoundaries[0]=new double[] {MRUtils.epsilon,maxMag*MRUtils.factorT2M0MaxRatio};
					parametersBoundaries[1]=new double[] {minTr*MRUtils.factorT1MinRatio,maxTr*MRUtils.factorT1MaxRatio};
					parametersBoundaries[2]=new double[] {minTe*MRUtils.factorT2MinRatio,maxTe*MRUtils.factorT2MaxRatio};
					parametersBoundaries[3]=new double[] {0,maxTe*MRUtils.factorT2MaxRatio};
               }
			   else parametersBoundaries=getDefaultParametersBoundaries(fitType);
               break;
           case MRUtils.T1T2_DEFAULT_T2_MULTI_RICE:
        	   simp[0][0] = medMag;
        	   simp[0][1] = medMag;
        	   simp[0][2] = medTr;
        	   simp[0][3] = minTe+MRUtils.epsilon;
        	   simp[0][4] = maxTe+MRUtils.epsilon;
        	   simp[0][5] = 2*minTe;//DeltaT2
               if(MRUtils.useBoundaries) {
					parametersBoundaries[0]=new double[] {MRUtils.epsilon,maxMag*MRUtils.factorT2M0MaxRatio};
					parametersBoundaries[1]=new double[] {MRUtils.epsilon,maxMag*MRUtils.factorT2M0MaxRatio};
					parametersBoundaries[2]=new double[] {minTr*MRUtils.factorT1MinRatio,maxTr*MRUtils.factorT1MaxRatio};
					parametersBoundaries[3]=new double[] {minTe*MRUtils.factorT2MinRatio,maxTe*MRUtils.factorT2MaxRatio};
					parametersBoundaries[4]=new double[] {minTe*MRUtils.factorT2MinRatio,maxTe*MRUtils.factorT2MaxRatio};
					parametersBoundaries[5]=new double[] {0,maxTe*MRUtils.factorT2MaxRatio};
               }
			   else parametersBoundaries=getDefaultParametersBoundaries(fitType);
               break;
               default:break;
        }
    }
        
    /**
     *  Restart the simplex at the nth vertex.
     *
     * @param n the n
     */
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
        	if(debug)System.out.println("Testing boundaries of vertices : ");
        	if(debug)System.out.println(TransformUtils.stringVectorN(simp[i], " "));  
            for (int j = 0; j < numParams; j++) {
                simp[i][j] = simp[i-1][j] + q[j];
            }
            simp[i][i-1] = simp[i][i-1] + p[i-1];
            for (int j = 0; j < numParams; j++) {
            	//System.out.println(parametersBoundaries[j][0]+"  -  "+simp[i][j]+"   -   "+parametersBoundaries[j][1]);
            	//if(simp[i][j]<parametersBoundaries[j][0]) {System.out.println("YES");simp[i][j]=parametersBoundaries[j][0];}
                //if(simp[i][j]>parametersBoundaries[j][1]) {System.out.println("YES");simp[i][j]=parametersBoundaries[j][1];}
             }
        
            sumResiduals(simp[i]);
        }
        // Initialise current lowest/highest parametre estimates to simplex 1
        best = 0;
        worst = 0;
        nextWorst = 0;
        order();
        nRestarts++;
    }
        
    /**
     * Show simplex. Display it in the log window, as text
     *
     * @param iter the iter
     */
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
        
        
	/**
	 *  Returns formula value for parameters 'p' at 'x'.
	 *
	 * @param p the p
	 * @param tr the tr
	 * @param te the te
	 * @return the double
	 */
	public double f(double[] p, double tr,double te) {
		return f(fit, p, tr,te);
	}

	
   /**
    *  Returns 'fit' formula value for parameters "p" at "x".
    *
    * @param fit the fit
    * @param p the p
    * @param tr the tr
    * @param te the te
    * @return the double
    */
    public double f(int fit, double[] p, double tr,double te) {
    	return MRUtils.getFitFunctionValue(tr, te, p, sigma, fit);
    }
    
    /**
     *  Get the set of parameter values from the best corner of the simplex.
     *
     * @return the params
     */
    public double[] getParams() {
        if(fit==MRUtils.T1T2_BIONANO)return bionanoParams;
    	order();
        return simp[best];
    }
    
	/**
	 *  Returns residuals array ie. differences between data and curve.
	 *
	 * @return the residuals
	 */
	public double[] getResiduals() {
		int saveFit = fit;
		double[] params = getParams();
		double[] residuals = new double[numPoints];
		for (int i=0; i<numPoints; i++)
			residuals[i] = magData[i] - f(fit, params, trData[i],teData[i]);
		fit = saveFit;
		return residuals;
	}
    
    /**
     * Gets the sum residuals sqr.
     *
     * @return the sum residuals sqr
     */
    /* Last "parametre" at each vertex of simplex is sum of residuals
     * for the curve described by that vertex
     */
    public double getSumResidualsSqr() {
        double sumResidualsSqr = (getParams())[MRUtils.getNparams(fit)];
        return sumResidualsSqr;
    }
    
    /**
     *   Returns the standard deviation of the residuals.
     *
     * @return the sd
     */
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
    
    /**
     *  Returns R^2, where 1.0 is best.
     *     <pre>
     *      r^2 = 1 - SSE/SSD
     *      
     *      where:	 SSE = sum of the squares of the errors
     *                  SSD = sum of the squares of the deviations about the mean.
     *     </pre>
     *
     * @return the r squared
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

    /**
     *   Get a measure of "goodness of fit" where 1.0 is best.
     *
     * @return the fit goodness
     */
    public double getFitGoodness() {
        double sumY = 0.0;
        for (int i = 0; i < numPoints; i++) sumY += magData[i];
        double mean = sumY / numPoints;
        double sumMeanDiffSqr = 0.0;
        int degreesOfFreedom = numPoints - MRUtils.getNparams(fit);
        double fitGoodness = 0.0;
        for (int i = 0; i < numPoints; i++) {
            sumMeanDiffSqr += sqr(magData[i] - mean);
        }
        if (sumMeanDiffSqr > 0.0 && degreesOfFreedom != 0)
            fitGoodness = 1.0 - (getSumResidualsSqr() / degreesOfFreedom) * ((numPoints) / sumMeanDiffSqr);
        
        return fitGoodness;
    }
    
    /**
     *  Get a string description of the curve fitting results
     * for easy output.
     *
     * @param d the d
     * @return the double
     */
        
    double sqr(double d) { return d * d; }
    
	/**
	 *  Adds sum of square of residuals to end of array of parameters.
	 *
	 * @param x the x
	 */
	void sumResiduals (double[] x) {
		x[numParams] = 0.0;
		for (int i=0; i<numPoints; i++) {
			x[numParams] = x[numParams] + sqr(f(fit,x,trData[i],teData[i])-magData[i]);
		}
	}

    /**
     *  Keep the "next" vertex.
     */
    void newVertex() {
    	for (int i = 0; i < numVertices; i++)
            simp[worst][i] = next[i];
    }
    
    /**
     *  Find the worst, nextWorst and best current set of parameter estimates.
     */
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
    }

    /**
     *  Get number of iterations performed.
     *
     * @return the iterations
     */
    public int getIterations() {
        return numIter;
    }
    
    /**
     *  Get maximum number of iterations allowed.
     *
     * @return the max iterations
     */
    public int getMaxIterations() {
        return maxIter;
    }
    
    /**
     *  Set maximum number of iterations allowed.
     *
     * @param x the new max iterations
     */
    public void setMaxIterations(int x) {
        maxIter = x;
    }
    
    /**
     *  Get number of simplex restarts to do.
     *
     * @return the restarts
     */
    public int getRestarts() {
        return defaultRestarts;
    }
    
    /**
     *  Set number of simplex restarts to do.
     *
     * @param n the new restarts
     */
    public void setRestarts(int n) {
        defaultRestarts = n;
    }

	/**
	 *  Sets the initial parameters, which override the default initial parameters.
	 *
	 * @param params the new initial parameters
	 */
	public void setInitialParameters(double[] params) {
		initialParams = params;
	}

    /**
     * Gets index of highest value in an array.
     *
     * @param array the array
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
    
	/**
	 * Gets the te points.
	 *
	 * @return the te points
	 */
	public double[] getTePoints() {
		return teData;
	}
	
	/**
	 * Gets the tr points.
	 *
	 * @return the tr points
	 */
	public double[] getTrPoints() {
		return trData;
	}

	/**
	 * Gets the mag points.
	 *
	 * @return the mag points
	 */
	public double[] getMagPoints() {
		return magData;
	}
	
	/**
	 * Gets the fit.
	 *
	 * @return the fit
	 */
	public int getFit() {
		return fit;
	}

	
	
	
	
    
    /**
     * Compute approximative parameters using the three first echoes. Used as benchmark 
     */
    public void computeBioNano(){
    	initialize2DTab();
    	inspect2DTab();
    	computeT2Bionano();
    	inspectT2EstimatedValues(); 
    	computeT1Bionano();
    }
    
	/**
	 * Initialize 2 D tab for storing successive values of Repetition times and Echo times
	 */
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
	
	/**
	 * Inspect 2 D tab storing successive values of Repetition times and Echo times
	 */
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

	/**
	 * Inspect values when computing approximative parameters using the three first echoes. Used as benchmark 
	 */
	public void inspectT2EstimatedValues() {
		if(debugBionano)System.out.println("T2 vas computed. Values : T2="+this.t2+"  , R2="+this.r2+"  , M0="+this.m0t2);
	}
	
	/**
	 * Computing approximative value of T2 using the three first echoes. Used as benchmark 
	 */
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
	
	
	/**
	 * Computing approximative value of T1 using the three first recovery times. Used as benchmark 
	 */
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
    
	
	/**
	 * Used in benchmark for estimative approximation of parrameters. A replicate of the Hampel function of Matlab
	 *
	 * @param tab the tab
	 * @return the double
	 */
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
	

	
	
	
	
	
}