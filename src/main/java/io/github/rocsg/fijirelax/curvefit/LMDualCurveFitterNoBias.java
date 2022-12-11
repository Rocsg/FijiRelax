package io.github.rocsg.fijirelax.curvefit;
import io.github.rocsg.fijiyama.common.VitimageUtils;
import io.github.rocsg.fijirelax.lma.LMAMatrix.InvertException;
import io.github.rocsg.fijirelax.mrialgo.MRUtils;
import io.github.rocsg.fijirelax.lma.LMAMultiDimFunction;
import io.github.rocsg.fijirelax.lma.LMA;


/**
 *  Implements the MRLMA capabilities for fitting exponential curves to match observation points
 * 	Copyright (C) 2022  <io.github.rocsg>
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
public class LMDualCurveFitterNoBias {	
	
	/** The debug LM. */
	public boolean debugLM=false;
	
	/** The Constant timeunits. */
	public final static String[] timeunits={"ms", "s"};
	
	/** The Constant timeitems. */
	public final static int[] timeitems={MRUtils.MSEC, MRUtils.SEC};
	
	/** The Constant fititems2. */
	public final static String[] fititems2={"Simplex","Levenberg-Marquardt"};
	
	/** The Constant constitems2. */
	public final static int[] constitems2={MRUtils.SIMPLEX,MRUtils.LM};
	
	/** The fit. */
	protected int fit;                // Number of curve type to fit

	/** The data. */
	protected double[][] data; // x,y data to fit
	
	/** The lma. */
	private MRLMA lma;
	
	/** The parameters. */
	private double[] parameters;
	
	/** The parameters boundaries. */
	private float[][] parametersBoundaries;
	
	/** The gfit. */
	private float gfit;
	
	/** The iterations break. */
	public boolean iterationsBreak=false;
	
	/** The lambda. */
	public static double lambda = 0.001;//lambda = 0.01 in original version;
	
	/** The min delta chi 2. */
	public static double minDeltaChi2 = 1e-7;
	
	/** The max iterations. */
	public static int maxIterations=200;

		
	
    /**
     * Gets the default parameters boundaries i.e. the acceptable boundaries for the exponential curve to estimate.
     *
     * @param fitType the fit type
     * @return the default parameters boundaries
     */
    public float[][] getDefaultParametersBoundaries(int fitType) {
    	float[][]parametersBoundaries=new float[MRUtils.getNparams(fitType)+1][2];
    	for(int i=0;i<parametersBoundaries.length;i++)parametersBoundaries[i]=new float[] {(float) -MRUtils.infinity,(float) MRUtils.infinity};
    	return parametersBoundaries;
    }


	/**
	 *  Construct a new CurveFitter object to proceed to exponential fit
	 *
	 * @param trData Recovery times (along T1 relax)
	 * @param teData the te times (along T2 relax)
	 * @param yData the magnitude observed
	 * @param fitType the fit type, chosen among @link io.github.rocsg.fijirelax.mrialgo.MRDataType.Class
	 * @param sigma the estimated sigma of the Rice distribution
	 * @param debugLM a debug flag producing additional verbose output data
	 */
	public LMDualCurveFitterNoBias (double[] trData, double[] teData, double[] yData, int fitType,double sigma,boolean debugLM) {
		this.fit=fitType;
		this.debugLM=debugLM;
		int xlength=trData.length;
		int ylength=yData.length;
		if (xlength!= ylength)
			throw new IllegalArgumentException("Arrays' lengths do not match");

		data=new double[xlength][3];
		for(int i=0;i<xlength;i++) data[i]=new double[]{yData[i],trData[i],teData[i]};
        double minTr=VitimageUtils.min(trData);
        double maxTr=VitimageUtils.max(trData);
        double minTe=VitimageUtils.min(teData);
        double maxTe=VitimageUtils.max(teData);
        double medTr=minTr/2+maxTr/2;
        double medTe=minTe/2+maxTe/2;
        double maxMag=VitimageUtils.max(yData);
        double minMag=VitimageUtils.min(yData);
        double medMag=maxMag/2;

        int np=MRUtils.getNparams(fitType);
		double[] inparameters=new double[np];
		parameters=new double[np];
		parametersBoundaries=new float[np][2];
		switch (fit) {
		case MRUtils.T1_MONO:
			inparameters[0] = maxMag*1.1;
			inparameters[1] = medTr;
			if(MRUtils.useBoundaries) {
				parametersBoundaries[0]=new float[] {(float) MRUtils.epsilon,(float) ((maxMag)*MRUtils.factorT1M0MaxRatio)};
				parametersBoundaries[1]=new float[] {(float) (minTr*MRUtils.factorT1MinRatio),(float) (maxTr*MRUtils.factorT1MaxRatio)};
			}
			else parametersBoundaries=getDefaultParametersBoundaries(fitType);
			T1Mono t1mrmono=new T1Mono();
			t1mrmono.setSigma(sigma);
			t1mrmono.setLM(this);
			
			lma = new MRLMA(t1mrmono,inparameters,data,parametersBoundaries,debugLM,sigma);
			break;
		case MRUtils.T1_MONO_BIAS:
			inparameters[0] = maxMag*1.1;
			inparameters[1] = medTr;
			inparameters[2] = 0;
			if(MRUtils.useBoundaries) {
				parametersBoundaries[0]=new float[] {(float) MRUtils.epsilon,(float) ((maxMag)*MRUtils.factorT1M0MaxRatio)};
				parametersBoundaries[1]=new float[] {(float) (minTr*MRUtils.factorT1MinRatio),(float) (maxTr*MRUtils.factorT1MaxRatio)};
				parametersBoundaries[2]=new float[] {0,(float) maxMag};
			}
			else parametersBoundaries=getDefaultParametersBoundaries(fitType);
			T1MonoBias t1mrbias=new T1MonoBias();
			t1mrbias.setSigma(sigma);
			t1mrbias.setLM(this);
			
			lma = new MRLMA(t1mrbias,inparameters,data,parametersBoundaries,debugLM,sigma);
			break;
		case MRUtils.T1_MONO_RICE:
			inparameters[0] = maxMag*1.1;
			inparameters[1] = medTr;
			if(MRUtils.useBoundaries) {
				parametersBoundaries[0]=new float[] {(float) MRUtils.epsilon,(float) ((maxMag)*MRUtils.factorT1M0MaxRatio)};
				parametersBoundaries[1]=new float[] {(float) (minTr*MRUtils.factorT1MinRatio),(float) (maxTr*MRUtils.factorT1MaxRatio)};
			}
			else parametersBoundaries=getDefaultParametersBoundaries(fitType);
			T1MonoRice t1mr=new T1MonoRice();
			t1mr.setSigma(sigma);
			t1mr.setLM(this);
			
			lma = new MRLMA(t1mr,inparameters,data,parametersBoundaries,debugLM,sigma);
			break;
			
			
			
			
			

		case MRUtils.T2_MONO:
			inparameters[0] = maxMag*1.1;
			inparameters[1] = medTe;
			if(MRUtils.useBoundaries) {
				parametersBoundaries[0]=new float[] {(float) MRUtils.epsilon,(float) ((maxMag)*MRUtils.factorT2M0MaxRatio)};
				parametersBoundaries[1]=new float[] {(float) MRUtils.epsilon,(float) ((maxMag)*MRUtils.factorT2M0MaxRatio)};
			}
			else parametersBoundaries=getDefaultParametersBoundaries(fitType);
			T2Mono t2mrmono=new T2Mono();
			t2mrmono.setSigma(sigma);
			t2mrmono.setLM(this);
			lma = new MRLMA(t2mrmono,inparameters,data,parametersBoundaries,debugLM,sigma);
			break;
		case MRUtils.T2_MONO_BIAS:
			inparameters[0] = maxMag*1.1;
			inparameters[1] = medTe;
			inparameters[2] = minMag;
			if(MRUtils.useBoundaries) {
				parametersBoundaries[0]=new float[] {(float) MRUtils.epsilon,(float) ((maxMag)*MRUtils.factorT2M0MaxRatio)};
				parametersBoundaries[1]=new float[] {(float) MRUtils.epsilon,(float) ((maxMag)*MRUtils.factorT2M0MaxRatio)};
				parametersBoundaries[2]=new float[] {0,(float) (maxMag)};
			}
			else parametersBoundaries=getDefaultParametersBoundaries(fitType);
			T2MonoBias t2mrb=new T2MonoBias();
			t2mrb.setSigma(sigma);
			t2mrb.setLM(this);
			lma = new MRLMA(t2mrb,inparameters,data,parametersBoundaries,debugLM,sigma);
			break;
		case MRUtils.T2_MONO_RICE:
			inparameters[0] = maxMag*1.1;
			inparameters[1] = medTe;
			if(MRUtils.useBoundaries) {
				parametersBoundaries[0]=new float[] {(float) MRUtils.epsilon,(float) ((maxMag)*MRUtils.factorT2M0MaxRatio)};
				parametersBoundaries[1]=new float[] {(float) MRUtils.epsilon,(float) ((maxMag)*MRUtils.factorT2M0MaxRatio)};
			}
			else parametersBoundaries=getDefaultParametersBoundaries(fitType);
			T2MonoRice t2mr=new T2MonoRice();
			t2mr.setSigma(sigma);
			t2mr.setLM(this);
			lma = new MRLMA(t2mr,inparameters,data,parametersBoundaries,debugLM,sigma);
			break;

			
			
			
		case MRUtils.T2_MULTI:
			inparameters[0] = medMag*1.6;
			inparameters[1] = medMag*0.6;
			inparameters[2] = medTe*0.4;
			inparameters[3] = medTe*1.3;
			if(MRUtils.useBoundaries) {
				parametersBoundaries[0]=new float[] {(float) MRUtils.epsilon,(float) ((maxMag)*MRUtils.factorT1M0MaxRatio)};
				parametersBoundaries[1]=new float[] {(float) MRUtils.epsilon,(float) ((maxMag)*MRUtils.factorT1M0MaxRatio)};
				parametersBoundaries[2]=new float[] {(float) (minTe*MRUtils.factorT2MinRatio),(float) (maxTe*MRUtils.factorT2MaxRatio)};
				parametersBoundaries[3]=new float[] {(float) (minTe*MRUtils.factorT2MinRatio),(float) (maxTe*MRUtils.factorT2MaxRatio)};
			}
			else parametersBoundaries=getDefaultParametersBoundaries(fitType);
			T2Multi t2Mr=new T2Multi();
			t2Mr.setSigma(sigma);
			t2Mr.setLM(this);
			lma = new MRLMA(t2Mr,inparameters,data,parametersBoundaries,debugLM,sigma);
			break;
		case MRUtils.T2_MULTI_BIAS:
			inparameters[0] = medMag*1.6;
			inparameters[1] = medMag*0.6;
			inparameters[2] = medTe*0.4;
			inparameters[3] = medTe*1.3;
			inparameters[4] = minMag;
			if(MRUtils.useBoundaries) {
				parametersBoundaries[0]=new float[] {(float) MRUtils.epsilon,(float) ((maxMag)*MRUtils.factorT1M0MaxRatio)};
				parametersBoundaries[1]=new float[] {(float) MRUtils.epsilon,(float) ((maxMag)*MRUtils.factorT1M0MaxRatio)};
				parametersBoundaries[2]=new float[] {(float) (minTe*MRUtils.factorT2MinRatio),(float) (maxTe*MRUtils.factorT2MaxRatio)};
				parametersBoundaries[3]=new float[] {(float) (minTe*MRUtils.factorT2MinRatio),(float) (maxTe*MRUtils.factorT2MaxRatio)};
				parametersBoundaries[4]=new float[] {(float) MRUtils.epsilon,(float) ((maxMag)*MRUtils.factorT1M0MaxRatio)};
			}
			else parametersBoundaries=getDefaultParametersBoundaries(fitType);
			T2MultiBias t2MrBias=new T2MultiBias();
			t2MrBias.setSigma(sigma);
			t2MrBias.setLM(this);
			lma = new MRLMA(t2MrBias,inparameters,data,parametersBoundaries,debugLM,sigma);
			break;
		case MRUtils.T2_MULTI_RICE:
			inparameters[0] = medMag*1.6;
			inparameters[1] = medMag*0.6;
			inparameters[2] = medTe*0.4;
			inparameters[3] = medTe*1.3;
			if(MRUtils.useBoundaries) {
				parametersBoundaries[0]=new float[] {(float) MRUtils.epsilon,(float) ((maxMag)*MRUtils.factorT1M0MaxRatio)};
				parametersBoundaries[1]=new float[] {(float) MRUtils.epsilon,(float) ((maxMag)*MRUtils.factorT1M0MaxRatio)};
				parametersBoundaries[3]=new float[] {(float) (minTe*MRUtils.factorT2MinRatio),(float) (maxTe*MRUtils.factorT2MaxRatio)};
				parametersBoundaries[4]=new float[] {(float) (minTe*MRUtils.factorT2MinRatio),(float) (maxTe*MRUtils.factorT2MaxRatio)};
			}
			else parametersBoundaries=getDefaultParametersBoundaries(fitType);
			T2MultiRice t2MrRice=new T2MultiRice();
			t2MrRice.setSigma(sigma);
			t2MrRice.setLM(this);
			lma = new MRLMA(t2MrRice,inparameters,data,parametersBoundaries,debugLM,sigma);
			break;
			
			
			
			
		
		case MRUtils.T1T2_MONO:
			inparameters[0] = maxMag*1.1;
			inparameters[1] = medTr;
			inparameters[2] = medTe;
			if(MRUtils.useBoundaries) {
				parametersBoundaries[0]=new float[] {(float) MRUtils.epsilon,(float) ((maxMag)*MRUtils.factorT1M0MaxRatio)};
				parametersBoundaries[1]=new float[] {(float) (minTr*MRUtils.factorT1MinRatio),(float) (maxTr*MRUtils.factorT1MaxRatio)};
				parametersBoundaries[2]=new float[] {(float) (minTe*MRUtils.factorT2MinRatio),(float) (maxTe*MRUtils.factorT2MaxRatio)};
			}
			else parametersBoundaries=getDefaultParametersBoundaries(fitType);
			T1T2Mono t1t2mr=new T1T2Mono();
			t1t2mr.setSigma(sigma);
			t1t2mr.setLM(this);
			
			lma = new MRLMA(t1t2mr,inparameters,data,parametersBoundaries,debugLM,sigma);
			break;
		case MRUtils.T1T2_MONO_BIAS:
			inparameters[0] = maxMag*1.1;
			inparameters[1] = medTr;
			inparameters[2] = medTe;
			inparameters[3] = minMag;
			if(MRUtils.useBoundaries) {
				parametersBoundaries[0]=new float[] {(float) MRUtils.epsilon,(float) ((maxMag)*MRUtils.factorT1M0MaxRatio)};
				parametersBoundaries[1]=new float[] {(float) (minTr*MRUtils.factorT1MinRatio),(float) (maxTr*MRUtils.factorT1MaxRatio)};
				parametersBoundaries[2]=new float[] {(float) (minTe*MRUtils.factorT2MinRatio),(float) (maxTe*MRUtils.factorT2MaxRatio)};
				parametersBoundaries[3]=new float[] {(float) MRUtils.epsilon,(float) ((maxMag)*MRUtils.factorT1M0MaxRatio)};
			}
			else parametersBoundaries=getDefaultParametersBoundaries(fitType);
			T1T2MonoBias t1t2mrbias=new T1T2MonoBias();
			t1t2mrbias.setSigma(sigma);
			t1t2mrbias.setLM(this);
			
			lma = new MRLMA(t1t2mrbias,inparameters,data,parametersBoundaries,debugLM,sigma);
			break;
		case MRUtils.T1T2_MONO_RICE:
			inparameters[0] = maxMag*1.1;
			inparameters[1] = medTr;
			inparameters[2] = medTe;
			if(MRUtils.useBoundaries) {
				parametersBoundaries[0]=new float[] {(float) MRUtils.epsilon,(float) ((maxMag)*MRUtils.factorT1M0MaxRatio)};
				parametersBoundaries[1]=new float[] {(float) (minTr*MRUtils.factorT1MinRatio),(float) (maxTr*MRUtils.factorT1MaxRatio)};
				parametersBoundaries[2]=new float[] {(float) (minTe*MRUtils.factorT2MinRatio),(float) (maxTe*MRUtils.factorT2MaxRatio)};
			}
			else parametersBoundaries=getDefaultParametersBoundaries(fitType);
			T1T2MonoRice t1t2mrRice=new T1T2MonoRice();
			t1t2mrRice.setSigma(sigma);
			t1t2mrRice.setLM(this);
			
			lma = new MRLMA(t1t2mrRice,inparameters,data,parametersBoundaries,debugLM,sigma);
			break;

			
			
			
		case MRUtils.T1T2_MULTI:
			inparameters[0] = medMag*1.6;
			inparameters[1] = medMag*0.6;
			inparameters[2] = medTr;
			inparameters[3] = medTe*0.4;
			inparameters[4] = medTe*1.3;
			if(MRUtils.useBoundaries) {
				parametersBoundaries[0]=new float[] {(float) MRUtils.epsilon,(float) ((maxMag)*MRUtils.factorT1M0MaxRatio)};
				parametersBoundaries[1]=new float[] {(float) MRUtils.epsilon,(float) ((maxMag)*MRUtils.factorT1M0MaxRatio)};
				parametersBoundaries[2]=new float[] {(float) (minTr*MRUtils.factorT1MinRatio),(float) (maxTr*MRUtils.factorT1MaxRatio)};
				parametersBoundaries[3]=new float[] {(float) (minTe*MRUtils.factorT2MinRatio),(float) (maxTe*MRUtils.factorT2MaxRatio)};
				parametersBoundaries[4]=new float[] {(float) (minTe*MRUtils.factorT2MinRatio),(float) (maxTe*MRUtils.factorT2MaxRatio)};
			}
			else parametersBoundaries=getDefaultParametersBoundaries(fitType);
			T1T2Multi t1t2Mr=new T1T2Multi();
			t1t2Mr.setSigma(sigma);
			t1t2Mr.setLM(this);
			lma = new MRLMA(t1t2Mr,inparameters,data,parametersBoundaries,debugLM,sigma);
			break;
		case MRUtils.T1T2_MULTI_BIAS:
			inparameters[0] = medMag*1.6;
			inparameters[1] = medMag*0.6;
			inparameters[2] = medTr;
			inparameters[3] = medTe*0.4;
			inparameters[4] = medTe*1.3;
			inparameters[5] = maxMag;
			if(MRUtils.useBoundaries) {
				parametersBoundaries[0]=new float[] {(float) MRUtils.epsilon,(float) ((maxMag)*MRUtils.factorT1M0MaxRatio)};
				parametersBoundaries[1]=new float[] {(float) MRUtils.epsilon,(float) ((maxMag)*MRUtils.factorT1M0MaxRatio)};
				parametersBoundaries[2]=new float[] {(float) (minTr*MRUtils.factorT1MinRatio),(float) (maxTr*MRUtils.factorT1MaxRatio)};
				parametersBoundaries[3]=new float[] {(float) (minTe*MRUtils.factorT2MinRatio),(float) (maxTe*MRUtils.factorT2MaxRatio)};
				parametersBoundaries[4]=new float[] {(float) (minTe*MRUtils.factorT2MinRatio),(float) (maxTe*MRUtils.factorT2MaxRatio)};
				parametersBoundaries[5]=new float[] {(float) MRUtils.epsilon,(float) ((maxMag)*MRUtils.factorT1M0MaxRatio)};
			}
			else parametersBoundaries=getDefaultParametersBoundaries(fitType);
			T1T2MultiBias t1t2MrBias=new T1T2MultiBias();
			t1t2MrBias.setSigma(sigma);
			t1t2MrBias.setLM(this);
			lma = new MRLMA(t1t2MrBias,inparameters,data,parametersBoundaries,debugLM,sigma);
			break;
		case MRUtils.T1T2_MULTI_RICE:
			inparameters[0] = medMag*1.6;
			inparameters[1] = medMag*0.6;
			inparameters[2] = medTr;
			inparameters[3] = medTe*0.4;
			inparameters[4] = medTe*1.3;
			if(MRUtils.useBoundaries) {
				parametersBoundaries[0]=new float[] {(float) MRUtils.epsilon,(float) ((maxMag)*MRUtils.factorT1M0MaxRatio)};
				parametersBoundaries[1]=new float[] {(float) MRUtils.epsilon,(float) ((maxMag)*MRUtils.factorT1M0MaxRatio)};
				parametersBoundaries[2]=new float[] {(float) (minTr*MRUtils.factorT1MinRatio),(float) (maxTr*MRUtils.factorT1MaxRatio)};
				parametersBoundaries[3]=new float[] {(float) (minTe*MRUtils.factorT2MinRatio),(float) (maxTe*MRUtils.factorT2MaxRatio)};
				parametersBoundaries[4]=new float[] {(float) (minTe*MRUtils.factorT2MinRatio),(float) (maxTe*MRUtils.factorT2MaxRatio)};
			}
			else parametersBoundaries=getDefaultParametersBoundaries(fitType);
			T1T2MultiRice t1t2MrRice=new T1T2MultiRice();
			t1t2MrRice.setSigma(sigma);
			t1t2MrRice.setLM(this);
			lma = new MRLMA(t1t2MrRice,inparameters,data,parametersBoundaries,debugLM,sigma);
			break;

			
			
			
		case MRUtils.T1T2_DEFAULT_T2_MONO_RICE:
			inparameters[0] = maxMag*1.1;
			inparameters[1] = medTr;
			inparameters[2] = medTe;
			inparameters[3] = minTe*2;
			if(MRUtils.useBoundaries) {
				parametersBoundaries[0]=new float[] {(float) MRUtils.epsilon,(float) ((maxMag)*MRUtils.factorT1M0MaxRatio)};
				parametersBoundaries[1]=new float[] {(float) (minTr*MRUtils.factorT1MinRatio),(float) (maxTr*MRUtils.factorT1MaxRatio)};
				parametersBoundaries[2]=new float[] {(float) (minTe*MRUtils.factorT2MinRatio),(float) (maxTe*MRUtils.factorT2MaxRatio)};
				parametersBoundaries[3]=new float[] {(float) (0),(float) (maxTe)};
			}
			else parametersBoundaries=new float[][] {{(float) -MRUtils.infinity,(float) MRUtils.infinity},{(float) -MRUtils.infinity,(float) MRUtils.infinity},{(float) -MRUtils.infinity,(float) MRUtils.infinity},{(float) -MRUtils.infinity,(float) MRUtils.infinity},{(float) -MRUtils.infinity,(float) MRUtils.infinity}};
			T1T2DefaultMonoRice t1t2Defmr=new T1T2DefaultMonoRice();
			t1t2Defmr.setSigma(sigma);
			t1t2Defmr.setLM(this);
			
			lma = new MRLMA(t1t2Defmr,inparameters,data,parametersBoundaries,debugLM,sigma);
			break;

		case MRUtils.T1T2_DEFAULT_T2_MULTI_RICE:
			inparameters[0] = medMag*1.6;
			inparameters[1] = medMag*0.6;
			inparameters[2] = medTr;
			inparameters[3] = medTe*0.4;
			inparameters[4] = medTe*1.3;
			inparameters[5] = minTe*2;
			if(MRUtils.useBoundaries) {
				parametersBoundaries[0]=new float[] {(float) MRUtils.epsilon,(float) ((maxMag)*MRUtils.factorT1M0MaxRatio)};
				parametersBoundaries[1]=new float[] {(float) MRUtils.epsilon,(float) ((maxMag)*MRUtils.factorT1M0MaxRatio)};
				parametersBoundaries[2]=new float[] {(float) (minTr*MRUtils.factorT1MinRatio),(float) (maxTr*MRUtils.factorT1MaxRatio)};
				parametersBoundaries[3]=new float[] {(float) (minTe*MRUtils.factorT2MinRatio),(float) (maxTe*MRUtils.factorT2MaxRatio)};
				parametersBoundaries[4]=new float[] {(float) (minTe*MRUtils.factorT2MinRatio),(float) (maxTe*MRUtils.factorT2MaxRatio)};
				parametersBoundaries[5]=new float[] {(float) (0),(float) (maxTe)};
			}
			else parametersBoundaries=new float[][] {{(float) -MRUtils.infinity,(float) MRUtils.infinity},{(float) -MRUtils.infinity,(float) MRUtils.infinity},{(float) -MRUtils.infinity,(float) MRUtils.infinity},{(float) -MRUtils.infinity,(float) MRUtils.infinity},{(float) -MRUtils.infinity,(float) MRUtils.infinity},{(float) -MRUtils.infinity,(float) MRUtils.infinity},{(float) -MRUtils.infinity,(float) MRUtils.infinity}};
			T1T2DefaultMultiRice t1t2DefMr=new T1T2DefaultMultiRice();
			t1t2DefMr.setSigma(sigma);
			t1t2DefMr.setLM(this);
			lma = new MRLMA(t1t2DefMr,inparameters,data,parametersBoundaries,debugLM,sigma);
			break;
		
		default:
			break;
		}
	}


	/**
	 * Setters for LM parameters
	 *
	 * @param lambda the lambda
	 * @param minDeltaChi2 the min delta chi 2
	 * @param maxIterations the max iterations
	 */
	public void configLMA(double lambda, double minDeltaChi2, int maxIterations) {
		if (lma==null) return;
		lma.lambda=lambda;
		lma.minDeltaChi2=minDeltaChi2;
		lma.maxIterations=maxIterations;
	}

	/**
	 * Run the fit using the declared setup of the current object
	 */
	public void doFit() {
		try {
			lma.fit();
			this.iterationsBreak=lma.iterationsBreak;
			parameters=lma.parameters;
			gfit=lma.chi2Goodness();

		} catch (InvertException e) {
			for(int i=0;i<parameters.length;i++)parameters[i]=MRUtils.ERROR_VALUE;
		}
	}

	/**
	 * Run the fit using the declared setup of the current object and the given LM parameters
	 *
	 * @param lambda the lambda
	 * @param minDeltaChi2 the min delta chi 2
	 * @param maxIterations the max iterations
	 */
	public void doFit(double lambda, double minDeltaChi2, int maxIterations) {
		try {
			lma.fit( lambda,  minDeltaChi2,  maxIterations);
			parameters=lma.parameters;
			gfit=lma.chi2Goodness();
		} catch (InvertException e) {
			for(int i=0;i<parameters.length;i++)parameters[i]=MRUtils.ERROR_VALUE;}
	}

	/**
	 * Gets the estimated params of the fitted function.
	 *
	 * @return the params
	 */
	public double[] getParams() {
		double[]ret=new double[parameters.length+1];
		for(int i=0;i<parameters.length;i++)ret[i]=parameters[i];
		ret[parameters.length]=lma.chi2Goodness();
		return ret;
	}

	/**
	 * Gets the lma object.
	 *
	 * @return the lma
	 */
	public LMA getLMA() {
		return lma;
	}


	
	
	
	
	

	/**
	 * The Class T1Mono. Describes fitting of a T1 relaxation with a single T1 component, without offset and without Rice noise
	 */
	//////////////////////  1-  T1MONO  /////////////////////////////////////////////////////////////////////////////////////////
	public static class T1Mono extends LMAMultiDimFunction {
		
		/** The Constant Nparams. */
		public static final int Nparams=2;
		
		/** The sigma. */
		public double sigma=0;
		
		/** The lm. */
		public LMDualCurveFitterNoBias lm;
		
		/**
		 * Sets the optimiser.
		 *
		 * @param lm the new lm
		 */
		public void setLM(LMDualCurveFitterNoBias lm) {this.lm=lm;}
		
		/**
		 * Sets the Rice sigma
		 *
		 * @param sig the new sigma
		 */
		public void setSigma(double sig){
			sigma=sig;
		}

		/**
		 * Gets estimates of what would be the observed value in absence of Rice noise corruption
		 *
		 * @param Tr the recovery time
		 * @param Te the echo time
		 * @param a the parameters of the equation y=  Rice_sigma [ a[0]* (1 - Math.exp(-(Tr / a[1])) ]
		 * @return the non riced value
		 */
		public double getNonRicedY(double Tr,double Te, double[] a) {
			return a[0]*(1-Math.exp(-Tr/a[1]));
		}
		
		/**
		 * Get the equation currently estimated
		 *
		 * @return the string
		 */
		public String toString() {
			return "y=  [ a[0]* (1 - Math.exp(-(Tr / a[1])) ]";
		}

		/**
		 * Gets estimates of what would be the observed value in absence of Rice noise corruption
		 *
		 * @param x =[Tr,Te]
		 * @param a the parameters of the equation y=  [ a[0]* (1 - Math.exp(-(Tr / a[1])) ]
		 * @return the value in absence of any Rice corruption
		 */
		@Override
		public double getY(double[] x, double[] a) {
			double Tr=x[0];
			double Te=x[1];					
			return  getNonRicedY(Tr,Te,a) ; // a[1] - echo times
		}
		
		/**
		 * Gets the partial derivate of the signal according to one of the parameters
		 *
		 * @param x =[Tr,Te]
		 * @param a the parameters of the equation y=   [ a[0]* (1 - Math.exp(-(Tr / a[1])) ]
		 * @param parameterIndex the index of the parameter in a
		 * @return the partial derivate value around this point
		 */
		@Override
		public double getPartialDerivate(double[] x, double[] a, int parameterIndex) {
			double Tr=x[0];
			double Te=x[1];					
			switch (parameterIndex) {
			case 0: return  (1-Math.exp(-Tr/a[1]));  
			case 1: return  a[0]*Math.exp(-Tr/a[1])*Tr/(-a[1]*a[1]);
			} 
			throw new RuntimeException("No such parameter index: " + parameterIndex);
		}

	} // end class
	
	/**
	 * The Class T1MonoBias. Describes fitting of a T1 relaxation with a single T1 component, with offset and withoutRice noise
	 */
	//////////////////////  1-  T1MONOBIAS  /////////////////////////////////////////////////////////////////////////////////////////
	public static class T1MonoBias extends LMAMultiDimFunction {
		
		/** The Constant Nparams. */
		public static final int Nparams=3;
		
		/** The sigma. */
		public double sigma=0;
		
		/** The lm. */
		public LMDualCurveFitterNoBias lm;
		
		/**
		 * Sets the optimiser.
		 *
		 * @param lm the new lm
		 */
		public void setLM(LMDualCurveFitterNoBias lm) {this.lm=lm;}
		
		/**
		 * Sets the Rice sigma.
		 *
		 * @param sig the new sigma
		 */
		public void setSigma(double sig){
			sigma=sig;
		}

		/**
		 * Gets estimates of what would be the observed value in absence of Rice noise corruption
		 *
		 * @param Tr the recovery time
		 * @param Te the echo time
		 * @param a the parameters of the equation y=  Rice_sigma [ a[0]* (1 - Math.exp(-(Tr / a[1])) ]
		 * @return the non riced value
		 */
		public double getNonRicedY(double Tr,double Te, double[] a) {
			return a[0]*(1-Math.exp(-Tr/a[1]))+a[2];
		}
		
		/**
		 * Get the equation currently estimated
		 *
		 * @return the string
		 */
		public String toString() {
			return "y=   [ a[0]* (1 - Math.exp(-(Tr / a[1])) ] + a[2]";
		}

		/**
		 * Gets estimates of what would be the observed value in absence of Rice noise corruption
		 *
		 * @param x =[Tr,Te]
		 * @param a the parameters of the equation y=  [ a[0]* (1 - Math.exp(-(Tr / a[1])) ] + a[2]
		 * @return the value in absence of any Rice corruption
		 */
		@Override
		public double getY(double[] x, double[] a) {
			double Tr=x[0];
			double Te=x[1];					
			return  getNonRicedY(Tr,Te,a); // a[1] - echo times
		}
		
		/**
		 * Gets the partial derivate of the signal according to one of the parameters
		 *
		 * @param x =[Tr,Te]
		 * @param a the parameters of the equation y=   [ a[0]* (1 - Math.exp(-(Tr / a[1])) ] + a[2]
		 * @param parameterIndex the index of the parameter in a
		 * @return the partial derivate value around this point
		 */
		@Override
		public double getPartialDerivate(double[] x, double[] a, int parameterIndex) {
			double Tr=x[0];
			double Te=x[1];					
			switch (parameterIndex) {
			case 0: return  (1-Math.exp(-Tr/a[1]));  
			case 1: return  a[0]*Math.exp(-Tr/a[1])*Tr/(-a[1]*a[1]);
			case 2: return  1;
			} 
			throw new RuntimeException("No such parameter index: " + parameterIndex);
		}

	} // end class
	
	/**
	 * The Class T1MonoRice. Describes fitting of a T1 relaxation with a single T1 component, without offset and with Rice noise
	 */
	//////////////////////  1-  T1MONORICE  /////////////////////////////////////////////////////////////////////////////////////////
	public static class T1MonoRice extends LMAMultiDimFunction {
		
		/** The Constant Nparams. */
		public static final int Nparams=2;
		
		/** The sigma. */
		public double sigma=0;
		
		/** The lm. */
		public LMDualCurveFitterNoBias lm;
		
		/**
		 * Sets the optimiser.
		 *
		 * @param lm the new lm
		 */
		public void setLM(LMDualCurveFitterNoBias lm) {this.lm=lm;}
		
		/**
		 * Sets the Rice sigma.
		 *
		 * @param sig the new sigma
		 */
		public void setSigma(double sig){
			sigma=sig;
		}

		/**
		 * Gets estimates of what would be the observed value in absence of Rice noise corruption
		 *
		 * @param Tr the recovery time
		 * @param Te the echo time
		 * @param a the parameters of the equation y=  [ a[0]* (1 - Math.exp(-(Tr / a[1])) ]
		 * @return the non riced value
		 */
		public double getNonRicedY(double Tr,double Te, double[] a) {
			return a[0]*(1-Math.exp(-Tr/a[1]));
		}
		
		/**
		 * Get the equation currently estimated
		 *
		 * @return the string
		 */
		public String toString() {
			return "y=  Rice_sigma [ a[0]* (1 - Math.exp(-(Tr / a[1])) ]";
		}

		/**
		 * Gets estimates of what would be the observed value in absence of Rice noise corruption
		 *
		 * @param x =[Tr,Te]
		 * @param a the parameters of the equation y= Ricecorruption_sigma [ a[0]* (1 - Math.exp(-(Tr / a[1])) ]
		 * @return the value in absence of any Rice corruption
		 */
		@Override
		public double getY(double[] x, double[] a) {
			double Tr=x[0];
			double Te=x[1];					
			return  MRUtils.besFunkCost( getNonRicedY(Tr,Te,a) , sigma); // a[1] - echo times
		}
		
		/**
		 * Gets the partial derivate of the signal according to one of the parameters
		 *
		 * @param x =[Tr,Te]
		 * @param a the parameters of the equation y= Rice [ a[0]* (1 - Math.exp(-(Tr / a[1])) ] + a[2]
		 * @param parameterIndex the index of the parameter in a
		 * @return the partial derivate value around this point
		 */
		@Override
		public double getPartialDerivate(double[] x, double[] a, int parameterIndex) {
			double Tr=x[0];
			double Te=x[1];					
			switch (parameterIndex) {
			case 0: return  MRUtils.besFunkDeriv( getNonRicedY(Tr,Te,a) , sigma)*(1-Math.exp(-Tr/a[1]));  
			case 1: return  MRUtils.besFunkDeriv( getNonRicedY(Tr,Te,a) , sigma)*a[0]*Math.exp(-Tr/a[1])*Tr/(-a[1]*a[1]);
			} 
			throw new RuntimeException("No such parameter index: " + parameterIndex);
		}

	} // end class
	
	

	
	
	
	
	/**
	 * The Class T2Mono. Describes fitting of a T2 relaxation with a single T2 component, without offset and without Rice noise
	 */
	//////////////////////  2-  T2MONO  /////////////////////////////////////////////////////////////////////////////////////////
	public static class T2Mono extends LMAMultiDimFunction {
		
		/** The Constant Nparams. */
		public static final int Nparams=2;
		
		/** The sigma. */
		public double sigma=0;
		
		/** The lm. */
		public LMDualCurveFitterNoBias lm;
		
		/**
		 * Sets the optimiser
		 *
		 * @param lm the new lm
		 */
		public void setLM(LMDualCurveFitterNoBias lm) {this.lm=lm;}
		
		/**
		 * Sets the Rice sigma.
		 *
		 * @param sig the new sigma
		 */
		public void setSigma(double sig){
			sigma=sig;
		}

		/**
		 * Gets estimates of what would be the observed value in absence of Rice noise corruption
		 *
		 * @param Tr the recovery time
		 * @param Te the echo time
		 * @param a the parameters of the equation y=  [ a[0]* Math.exp(-(Te / a[1])) ]
		 * @return the non riced value
		 */
		public double getNonRicedY(double Tr,double Te, double[] a) {
			return a[0]* Math.exp(-(Te / a[1]));
		}
		
		/**
		 * Get the equation currently estimated
		 *
		 * @return the string
		 */
		public String toString() {
			return "y=   [ a[0]* Math.exp(-(Te / a[1])) ]";
		}

		/**
		 * Gets estimates of what would be the observed value in absence of Rice noise corruption
		 *
		 * @param x =[Tr,Te]
		 * @param a the parameters of the equation y= [ a[0]* Math.exp(-(Te / a[1])) ]
		 * @return the value in absence of any Rice corruption
		 */
		@Override
		public double getY(double[] x, double[] a) {
			double Tr=x[0];
			double Te=x[1];					
			return  getNonRicedY(Tr,Te,a); // a[1] - echo times
		}
		
		/**
		 * Gets the partial derivate of the signal according to one of the parameters
		 *
		 * @param x =[Tr,Te]
		 * @param a the parameters of the equation y=   [ a[0]* Math.exp(-(Te / a[1])) ]
		 * @param parameterIndex the index of the parameter in a
		 * @return the partial derivate value around this point
		 */
		@Override
		public double getPartialDerivate(double[] x, double[] a, int parameterIndex) {
			double Tr=x[0];
			double Te=x[1];					
			switch (parameterIndex) {
			case 0: return  (Math.exp(-Te/a[1]));  
			case 1: return  a[0]*Math.exp(-Te/a[1])*Te/(a[1]*a[1]);
			} 
			throw new RuntimeException("No such parameter index: " + parameterIndex);
		}

	} // end class
	
	/**
	 * The Class T2MonoBias. Describes fitting of a T2 relaxation with a single T2 component, with offset and without Rice noise
	 */
	//////////////////////  2-  T2MONOBIAS  /////////////////////////////////////////////////////////////////////////////////////////
	public static class T2MonoBias extends LMAMultiDimFunction {
		
		/** The Constant Nparams. */
		public static final int Nparams=3;
		
		/** The sigma. */
		public double sigma=0;
		
		/** The lm. */
		public LMDualCurveFitterNoBias lm;
		
		/**
		 * Sets the optimiser
		 *
		 * @param lm the new lm
		 */
		public void setLM(LMDualCurveFitterNoBias lm) {this.lm=lm;}
		
		/**
		 * Sets the Rice sigma.
		 *
		 * @param sig the new sigma
		 */
		public void setSigma(double sig){
			sigma=sig;
		}

		/**
		 * Gets estimates of what would be the observed value in absence of Rice noise corruption
		 *
		 * @param Tr the recovery time
		 * @param Te the echo time
		 * @param a the parameters of the equation y=   [ a[0]* Math.exp(-(Te / a[1])) ]
		 * @return the non riced value
		 */
		public double getNonRicedY(double Tr,double Te, double[] a) {
			return a[0]* Math.exp(-(Te / a[1]))+a[2];
		}
		
		/**
		 * Get the equation currently estimated
		 *
		 * @return the string
		 */
		public String toString() {
			return "y=  [ a[0]* Math.exp(-(Te / a[1])) ]";
		}

		/**
		 * Gets estimates of what would be the observed value in absence of Rice noise corruption
		 *
		 * @param x =[Tr,Te]
		 * @param a the parameters of the equation y=  [ a[0]* Math.exp(-(Te / a[1])) ]
		 * @return the value in absence of any Rice corruption
		 */
		@Override
		public double getY(double[] x, double[] a) {
			double Tr=x[0];
			double Te=x[1];					
			return  getNonRicedY(Tr,Te,a) ; // a[1] - echo times
		}
		
		/**
		 * Gets the partial derivate of the signal according to one of the parameters
		 *
		 * @param x =[Tr,Te]
		 * @param a the parameters of the equation y=   [ a[0]* Math.exp(-(Te / a[1])) ]
		 * @param parameterIndex the index of the parameter in a
		 * @return the partial derivate value around this point
		 */
		@Override
		public double getPartialDerivate(double[] x, double[] a, int parameterIndex) {
			double Tr=x[0];
			double Te=x[1];					
			switch (parameterIndex) {
			case 0: return  (Math.exp(-Te/a[1]));  
			case 1: return  a[0]*Math.exp(-Te/a[1])*Te/(a[1]*a[1]);
			case 2: return  1;
			} 
			throw new RuntimeException("No such parameter index: " + parameterIndex);
		}

	} // end class
	
	/**
	 * The Class T2MonoRice. Describes fitting of a T2 relaxation with a single T1 component, without offset and with Rice noise
	 */
	//////////////////////  2-  T2MONORICE  /////////////////////////////////////////////////////////////////////////////////////////
	public static class T2MonoRice extends LMAMultiDimFunction {
		
		/** The Constant Nparams. */
		public static final int Nparams=2;
		
		/** The sigma. */
		public double sigma=0;
		
		/** The lm. */
		public LMDualCurveFitterNoBias lm;
		
		/**
		 * Sets the optimiser
		 *
		 * @param lm the new lm
		 */
		public void setLM(LMDualCurveFitterNoBias lm) {this.lm=lm;}
		
		/**
		 * Sets the Rice sigma.
		 *
		 * @param sig the new sigma
		 */
		public void setSigma(double sig){
			sigma=sig;
		}

		/**
		 * Gets estimates of what would be the observed value in absence of Rice noise corruption
		 *
		 * @param Tr the recovery time
		 * @param Te the echo time
		 * @param a the parameters of the equation y=  Rice_sigma [ a[0]* Math.exp(-(Te / a[1])) ]
		 * @return the non riced value
		 */
		public double getNonRicedY(double Tr,double Te, double[] a) {
			return a[0]* Math.exp(-(Te / a[1]));
		}
		
		/**
		 * Get the equation currently estimated
		 *
		 * @return the string
		 */
		public String toString() {
			return "y=  Rice_sigma [ a[0]* Math.exp(-(Te / a[1])) ]";
		}

		/**
		 * Gets estimates of what would be the observed value in absence of Rice noise corruption
		 *
		 * @param x =[Tr,Te]
		 * @param a the parameters of the equation y= Rice_sigma [ a[0]* Math.exp(-(Te / a[1])) ]
		 * @return the value in absence of any Rice corruption
		 */
		@Override
		public double getY(double[] x, double[] a) {
			double Tr=x[0];
			double Te=x[1];					
			return  MRUtils.besFunkCost( getNonRicedY(Tr,Te,a) , sigma); // a[1] - echo times
		}
		
		/**
		 * Gets the partial derivate of the signal according to one of the parameters
		 *
		 * @param x =[Tr,Te]
		 * @param a the parameters of the equation y= Rice_sigma [ a[0]* Math.exp(-(Te / a[1])) ]
		 * @param parameterIndex the index of the parameter in a
		 * @return the partial derivate value around this point
		 */
		@Override
		public double getPartialDerivate(double[] x, double[] a, int parameterIndex) {
			double Tr=x[0];
			double Te=x[1];					
			switch (parameterIndex) {
			case 0: return  MRUtils.besFunkDeriv( getNonRicedY(Tr,Te,a) , sigma)*(Math.exp(-Te/a[1]));  
			case 1: return  MRUtils.besFunkDeriv( getNonRicedY(Tr,Te,a) , sigma)*a[0]*Math.exp(-Te/a[1])*Te/(a[1]*a[1]);
			} 
			throw new RuntimeException("No such parameter index: " + parameterIndex);
		}

	} // end class


	
	
	
	
	
	
	
	
	
	
	/**
	 * The Class T2Multi. Describes fitting of a T2 relaxation with two T2 components, without offset and without Rice noise
	 */
	//////////////////////  3-  T2MULTI  /////////////////////////////////////////////////////////////////////////////////////////
	public static class T2Multi extends LMAMultiDimFunction {
		
		/** The Constant Nparams. */
		public static final int Nparams=4;
		
		/** The sigma. */
		public double sigma=0;
		
		/** The lm. */
		public LMDualCurveFitterNoBias lm;
		
		/**
		 * Sets the optimiser
		 *
		 * @param lm the new lm
		 */
		public void setLM(LMDualCurveFitterNoBias lm) {this.lm=lm;}
		
		/**
		 * Sets the Rice sigma.
		 *
		 * @param sig the new sigma
		 */
		public void setSigma(double sig){
			sigma=sig;
		}

		
		/**
		 * Gets estimates of what would be the observed value in absence of Rice noise corruption
		 *
		 * @param x =[Tr,Te]
		 * @param a the parameters of the equation y=   [ a[0]* Math.exp(-(Te / a[2]))  +  a[1]* Math.exp(-(Te / a[3])) ]
		 * @return the non riced value
		 */
//		@Override
		public double getY(double[] x, double[] a) {
			double Tr=x[0];
			double Te=x[1];					
			return  getNonRicedY(Tr,Te,a); // a[1] - echo times
		}

		/**
		 * Gets estimates of what would be the observed value in absence of Rice noise corruption
		 *
		 * @param Tr the recovery time
		 * @param Te the echo time
		 * @param a the parameters of the equation y=  [ a[0]* Math.exp(-(Te / a[2]))  +  a[1]* Math.exp(-(Te / a[3])) ]
		 * @return the value in absence of any Rice corruption
		 */
		public double getNonRicedY(double Tr,double Te, double[] a) {
			return  ( a[0]* Math.exp(-(Te / a[2])) + a[1] * Math.exp(-(Te / a[3])));
		}
		

		/**
		 * The equation being optimized.
		 *
		 * @return the string
		 */
		public String toString() {
			return "y=   [ a[0]* Math.exp(-(Te / a[2]))  +  a[1]* Math.exp(-(Te / a[3])) ]";
		}
		
		/**
		 * Gets the partial derivate of the signal according to one of the parameters
		 *
		 * @param x =[Tr,Te]
		 * @param a the parameters of the equation y=   [ a[0]* Math.exp(-(Te / a[2]))  +  a[1]* Math.exp(-(Te / a[3])) ]
		 * @param parameterIndex the index of the parameter in a
		 * @return the partial derivate value around this point
		 */
		@Override
		public double getPartialDerivate(double[] x, double[] a, int parameterIndex) {
			double Tr=x[0];
			double Te=x[1];					
			switch (parameterIndex) {
			case 0: return  (Math.exp(-Te/a[2]));  
			case 1: return  (Math.exp(-Te/a[3]));  
			case 2: return  a[0]*Math.exp(-Te/a[2])*Te/(a[2]*a[2]);
			case 3: return  a[1]*Math.exp(-Te/a[3])*Te/(a[3]*a[3]);
			} 
			throw new RuntimeException("No such parameter index: " + parameterIndex);
		}
	} // end class
	
	/**
	 * The Class T2MultiBias. Describes fitting of a T2 relaxation with two T2 components, with offset and without Rice noise
	 */
	//////////////////////  3-  T2MULTIRICE  /////////////////////////////////////////////////////////////////////////////////////////
	public static class T2MultiBias extends LMAMultiDimFunction {
		
		/** The Constant Nparams. */
		public static final int Nparams=5;
		
		/** The sigma. */
		public double sigma=0;
		
		/** The lm. */
		public LMDualCurveFitterNoBias lm;
		
		/**
		 * Sets the optimiser
		 *
		 * @param lm the new lm
		 */
		public void setLM(LMDualCurveFitterNoBias lm) {this.lm=lm;}
		
		/**
		 * Sets the Rice sigma.
		 *
		 * @param sig the new sigma
		 */
		public void setSigma(double sig){
			sigma=sig;
		}

		
		/**
		 * Gets estimates of what would be the observed value in absence of Rice noise corruption
		 *
		 * @param x =[Tr,Te]
		 * @param a the parameters of the equation y=   [ a[0]* Math.exp(-(Te / a[2]))  +  a[1]* Math.exp(-(Te / a[3])) ]
		 * @return the non riced value
		 */
//		@Override
		public double getY(double[] x, double[] a) {
			double Tr=x[0];
			double Te=x[1];					
			return  getNonRicedY(Tr,Te,a); // a[1] - echo times
		}

		/**
		 * Gets estimates of what would be the observed value in absence of Rice noise corruption
		 *
		 * @param Tr the recovery time
		 * @param Te the echo time
		 * @param a the parameters of the equation y=   [ a[0]* Math.exp(-(Te / a[2]))  +  a[1]* Math.exp(-(Te / a[3])) ]
		 * @return the non riced value
		 */
		public double getNonRicedY(double Tr,double Te, double[] a) {
			return  ( a[0]* Math.exp(-(Te / a[2])) + a[1] * Math.exp(-(Te / a[3]))+a[4]);
		}
		

		/**
		 * Get the equation currently estimated
		 *
		 * @return the string
		 */
		public String toString() {
			return "y=  [ a[0]* Math.exp(-(Te / a[2]))  +  a[1]* Math.exp(-(Te / a[3])) ]";
		}
		
		/**
		 * Gets the partial derivate of the signal according to one of the parameters
		 *
		 * @param x =[Tr,Te]
		 * @param a the parameters of the equation y=  [ a[0]* Math.exp(-(Te / a[2]))  +  a[1]* Math.exp(-(Te / a[3])) ]
		 * @param parameterIndex the index of the parameter in a
		 * @return the partial derivate value around this point
		 */
		@Override
		public double getPartialDerivate(double[] x, double[] a, int parameterIndex) {
			double Tr=x[0];
			double Te=x[1];					
			switch (parameterIndex) {
			case 0: return  (Math.exp(-Te/a[2]));  
			case 1: return  (Math.exp(-Te/a[3]));  
			case 2: return  a[0]*Math.exp(-Te/a[2])*Te/(a[2]*a[2]);
			case 3: return  a[1]*Math.exp(-Te/a[3])*Te/(a[3]*a[3]);
			case 4: return  1;
			} 
			throw new RuntimeException("No such parameter index: " + parameterIndex);
		}
	} // end class
	
	/**
	 * The Class T2MultiRice. Describes fitting of a T2 relaxation with two T2 components, without offset and with Rice noise
	 */
	//////////////////////  3-  T2MULTIRICE  /////////////////////////////////////////////////////////////////////////////////////////
	public static class T2MultiRice extends LMAMultiDimFunction {
		
		/** The Constant Nparams. */
		public static final int Nparams=4;
		
		/** The sigma. */
		public double sigma=0;
		
		/** The lm. */
		public LMDualCurveFitterNoBias lm;
		
		/**
		 * Sets the optimiser
		 *
		 * @param lm the new lm
		 */
		public void setLM(LMDualCurveFitterNoBias lm) {this.lm=lm;}
		
		/**
		 * Sets the Rice sigma.
		 *
		 * @param sig the new sigma
		 */
		public void setSigma(double sig){
			sigma=sig;
		}

		
		/**
		 * Gets estimates of what would be the observed value in absence of Rice noise corruption
		 *
		 * @param x =[Tr,Te]
		 * @param a the parameters of the equation y= Rice_sigma [ a[0]* Math.exp(-(Te / a[2]))  +  a[1]* Math.exp(-(Te / a[3])) ]
		 * @return the value in absence of any Rice corruption
		 */
//		@Override
		public double getY(double[] x, double[] a) {
			double Tr=x[0];
			double Te=x[1];					
			return  MRUtils.besFunkCost( getNonRicedY(Tr,Te,a)  , sigma); // a[1] - echo times
		}

		/**
		 * Gets estimates of what would be the observed value in absence of Rice noise corruption
		 *
		 * @param Tr the recovery time
		 * @param Te the echo time
		 * @param a the parameters of the equation y=  Rice_sigma [ a[0]* Math.exp(-(Te / a[2]))  +  a[1]* Math.exp(-(Te / a[3])) ]
		 * @return the non riced value
		 */
		public double getNonRicedY(double Tr,double Te, double[] a) {
			return  ( a[0]* Math.exp(-(Te / a[2])) + a[1] * Math.exp(-(Te / a[3])));
		}
		

		/**
		 * Get the equation currently estimated
		 *
		 * @return the string
		 */
		public String toString() {
			return "y=  Rice_sigma [ a[0]* Math.exp(-(Te / a[2]))  +  a[1]* Math.exp(-(Te / a[3])) ]";
		}
		
		/**
		 * Gets the partial derivate of the signal according to one of the parameters
		 *
		 * @param x =[Tr,Te]
		 * @param a the parameters of the equation y= Rice_sigma [ a[0]* Math.exp(-(Te / a[2]))  +  a[1]* Math.exp(-(Te / a[3])) ]
		 * @param parameterIndex the index of the parameter in a
		 * @return the partial derivate value around this point
		 */
		@Override
		public double getPartialDerivate(double[] x, double[] a, int parameterIndex) {
			double Tr=x[0];
			double Te=x[1];					
			switch (parameterIndex) {
			case 0: return  MRUtils.besFunkDeriv( getNonRicedY(Tr,Te,a)  , sigma)*(Math.exp(-Te/a[2]));  
			case 1: return  MRUtils.besFunkDeriv( getNonRicedY(Tr,Te,a)  , sigma)*(Math.exp(-Te/a[3]));  
			case 2: return  MRUtils.besFunkDeriv( getNonRicedY(Tr,Te,a)  , sigma)*a[0]*Math.exp(-Te/a[2])*Te/(a[2]*a[2]);
			case 3: return  MRUtils.besFunkDeriv( getNonRicedY(Tr,Te,a)  , sigma)*a[1]*Math.exp(-Te/a[3])*Te/(a[3]*a[3]);
			} 
			throw new RuntimeException("No such parameter index: " + parameterIndex);
		}
	} // end class

	
	
	
	
	
	
	
	/**
	 * The Class T1T2Mono. Describes cross-fitting of a T1 and T2 relaxation with a single T2 component, without offset and without Rice noise
	 */
	//////////////////////  4-  T1T2MONO  /////////////////////////////////////////////////////////////////////////////////////////
	public static class T1T2Mono extends LMAMultiDimFunction {
		
		/** The Constant Nparams. */
		public static final int Nparams=3;
		
		/** The sigma. */
		public double sigma=0;
		
		/** The lm. */
		public LMDualCurveFitterNoBias lm;
		
		/**
		 * Sets the optimiser
		 *
		 * @param lm the new lm
		 */
		public void setLM(LMDualCurveFitterNoBias lm) {this.lm=lm;}
		
		/**
		 * Sets the Rice sigma.
		 *
		 * @param sig the new sigma
		 */
		public void setSigma(double sig){
			sigma=sig;
		}

		/**
		 * Gets estimates of what would be the observed value in absence of Rice noise corruption
		 *
		 * @param Tr the recovery time
		 * @param Te the echo time
		 * @param a the parameters of the equation y=  [ a[0]* Math.exp(-(Te / a[2]))* (1 - Math.exp(-(Tr / a[1])) ]
		 * @return the value in absence of any Rice corruption
		 */
		public double getNonRicedY(double Tr,double Te, double[] a) {
			return a[0]* Math.exp(-(Te / a[2]))*(1-Math.exp(-Tr/a[1]));
		}
		
		/**
		 * The equation being optimized.
		 *
		 * @return the string
		 */
		public String toString() {
			return "y=  [ a[0]* Math.exp(-(Te / a[2]))* (1 - Math.exp(-(Tr / a[1])) ]";
		}

		/**
		 * Gets estimates of what would be the observed value in absence of Rice noise corruption
		 *
		 * @param x =[Tr,Te]
		 * @param a the parameters of the equation y=   [ a[0]* Math.exp(-(Te / a[2]))  +  a[1]* Math.exp(-(Te / a[3])) ]
		 * @return the non riced value
		 */
		@Override
		public double getY(double[] x, double[] a) {
			double Tr=x[0];
			double Te=x[1];					
			return   getNonRicedY(Tr,Te,a); // a[1] - echo times
		}
		
		/**
		 * Gets the partial derivate of the signal according to one of the parameters
		 *
		 * @param x =[Tr,Te]
		 * @param a the parameters of the equation y=   [ a[0]* Math.exp(-(Te / a[2]))* (1 - Math.exp(-(Tr / a[1])) ]
		 * @param parameterIndex the index of the parameter in a
		 * @return the partial derivate value around this point
		 */
		@Override
		public double getPartialDerivate(double[] x, double[] a, int parameterIndex) {
			double Tr=x[0];
			double Te=x[1];					
			switch (parameterIndex) {
			case 0: return (1-Math.exp(-Tr/a[1]))*(Math.exp(-Te/a[2]));  
			case 1: return a[0]*Math.exp(-Te/a[2])*Math.exp(-Tr/a[1])*Tr/(-a[1]*a[1]);
			case 2: return a[0]*(1-Math.exp(-Tr/a[1]))*Math.exp(-Te/a[2])*Te/(a[2]*a[2]);
			} 
			throw new RuntimeException("No such parameter index: " + parameterIndex);
		}

	} // end class
	
	/**
	 * The Class T1T2MonoBias. Describes cross-fitting of a T1 and T2 relaxation with a single T2 component, with offset and without Rice noise
	 */
	//////////////////////  4-  T1T2MONOBIAS  /////////////////////////////////////////////////////////////////////////////////////////
	public static class T1T2MonoBias extends LMAMultiDimFunction {
		
		/** The Constant Nparams. */
		public static final int Nparams=4;
		
		/** The sigma. */
		public double sigma=0;
		
		/** The lm. */
		public LMDualCurveFitterNoBias lm;
		
		/**
		 * Sets the optimiser
		 *
		 * @param lm the new lm
		 */
		public void setLM(LMDualCurveFitterNoBias lm) {this.lm=lm;}
		
		/**
		 * Sets the Rice sigma.
		 *
		 * @param sig the new sigma
		 */
		public void setSigma(double sig){
			sigma=sig;
		}

		/**
		 * Gets estimates of what would be the observed value in absence of Rice noise corruption
		 *
		 * @param Tr the recovery time
		 * @param Te the echo time
		 * @param a the parameters of the equation y=  [ a[0]* Math.exp(-(Te / a[2]))* (1 - Math.exp(-(Tr / a[1])) ]
		 * @return the non riced value
		 */
		public double getNonRicedY(double Tr,double Te, double[] a) {
			return a[0]* Math.exp(-(Te / a[2]))*(1-Math.exp(-Tr/a[1]))+a[3];
		}
		
		/**
		 * Get the equation currently estimated
		 *
		 * @return the string
		 */
		public String toString() {
			return "y=  [ a[0]* Math.exp(-(Te / a[2]))* (1 - Math.exp(-(Tr / a[1])) ]";
		}

		/**
		 * Gets estimates of what would be the observed value in absence of Rice noise corruption
		 *
		 * @param x =[Tr, Te]
		 * @param a the parameters of the equation y=  [ a[0]* Math.exp(-(Te / a[2]))* (1 - Math.exp(-(Tr / a[1])) ]
		 * @return the non riced value
		 */
		@Override
		public double getY(double[] x, double[] a) {
			double Tr=x[0];
			double Te=x[1];					
			return  getNonRicedY(Tr,Te,a); // a[1] - echo times
		}
		
		/**
		 * Gets the partial derivate of the signal according to one of the parameters
		 *
		 * @param x =[Tr,Te]
		 * @param a the parameters of the equation y=  [ a[0]* Math.exp(-(Te / a[2]))* (1 - Math.exp(-(Tr / a[1])) ]
		 * @param parameterIndex the index of the parameter in a
		 * @return the partial derivate value around this point
		 */
		@Override
		public double getPartialDerivate(double[] x, double[] a, int parameterIndex) {
			double Tr=x[0];
			double Te=x[1];					
			switch (parameterIndex) {
			case 0: return  (1-Math.exp(-Tr/a[1]))*(Math.exp(-Te/a[2]));  
			case 1: return  a[0]*Math.exp(-Te/a[2])*Math.exp(-Tr/a[1])*Tr/(-a[1]*a[1]);
			case 2: return  a[0]*(1-Math.exp(-Tr/a[1]))*Math.exp(-Te/a[2])*Te/(a[2]*a[2]);
			case 3: return  1;
			} 
			throw new RuntimeException("No such parameter index: " + parameterIndex);
		}

	} // end class
	
	/**
	 * The Class T1T2MonoRice. Describes cross fitting of a T1-T2 relaxation with a single T2 component, without offset and with Rice noise
	 */
	//////////////////////  4-  T1T2MONORICE  /////////////////////////////////////////////////////////////////////////////////////////
	public static class T1T2MonoRice extends LMAMultiDimFunction {
		
		/** The Constant Nparams. */
		public static final int Nparams=3;
		
		/** The sigma. */
		public double sigma=0;
		
		/** The lm. */
		public LMDualCurveFitterNoBias lm;
		
		/**
		 * Sets the optimiser
		 *
		 * @param lm the new lm
		 */
		public void setLM(LMDualCurveFitterNoBias lm) {this.lm=lm;}
		
		/**
		 * Sets the Rice sigma.
		 *
		 * @param sig the new sigma
		 */
		public void setSigma(double sig){
			sigma=sig;
		}

		/**
		 * Gets estimates of what would be the observed value in absence of Rice noise corruption
		 *
		 * @param Tr the recovery time
		 * @param Te the echo time
		 * @param a the parameters of the equation y=   Rice_sigma [ a[0]* Math.exp(-(Te / a[2]))* (1 - Math.exp(-(Tr / a[1])) ]
		 * @return the non riced value
		 */
		public double getNonRicedY(double Tr,double Te, double[] a) {
			return a[0]* Math.exp(-(Te / a[2]))*(1-Math.exp(-Tr/a[1]));
		}
		
		/**
		 * Get the equation currently estimated
		 *
		 * @return the string
		 */
		public String toString() {
			return "y=  Rice_sigma [ a[0]* Math.exp(-(Te / a[2]))* (1 - Math.exp(-(Tr / a[1])) ]";
		}

		/**
		 * Gets estimates of what would be the observed value in absence of Rice noise corruption
		 *
		 * @param x =[Tr,Te]
		 * @param a the parameters of the equation y=  Rice_sigma [ a[0]* Math.exp(-(Te / a[2]))* (1 - Math.exp(-(Tr / a[1])) ]
		 * @return the value in absence of any Rice corruption
		 */
		@Override
		public double getY(double[] x, double[] a) {
			double Tr=x[0];
			double Te=x[1];					
			return  MRUtils.besFunkCost( getNonRicedY(Tr,Te,a) , sigma); // a[1] - echo times
		}
		
		/**
		 * Gets the partial derivate of the signal according to one of the parameters
		 *
		 * @param x =[Tr,Te]
		 * @param a the parameters of the equation y=  Rice_sigma [ a[0]* Math.exp(-(Te / a[2]))* (1 - Math.exp(-(Tr / a[1])) ]
		 * @param parameterIndex the index of the parameter in a
		 * @return the partial derivate value around this point
		 */
		@Override
		public double getPartialDerivate(double[] x, double[] a, int parameterIndex) {
			double Tr=x[0];
			double Te=x[1];					
			switch (parameterIndex) {
			case 0: return  MRUtils.besFunkDeriv( getNonRicedY(Tr,Te,a) , sigma)*(1-Math.exp(-Tr/a[1]))*(Math.exp(-Te/a[2]));  
			case 1: return  MRUtils.besFunkDeriv( getNonRicedY(Tr,Te,a) , sigma)*a[0]*Math.exp(-Te/a[2])*Math.exp(-Tr/a[1])*Tr/(-a[1]*a[1]);
			case 2: return  MRUtils.besFunkDeriv( getNonRicedY(Tr,Te,a) , sigma)*a[0]*(1-Math.exp(-Tr/a[1]))*Math.exp(-Te/a[2])*Te/(a[2]*a[2]);
			} 
			throw new RuntimeException("No such parameter index: " + parameterIndex);
		}

	} // end class

	

	
	
	
	/**
	 * The Class T1T2Multi. Describes cross-fitting of a T1 and T2 relaxation with a single T2 component, without offset and without Rice noise
	 */
	//////////////////////  5-  T1T2MULTI  /////////////////////////////////////////////////////////////////////////////////////////
	public static class T1T2Multi extends LMAMultiDimFunction {
		
		/** The Constant Nparams. */
		public static final int Nparams=5;
		
		/** The sigma. */
		public double sigma=0;
		
		/** The lm. */
		public LMDualCurveFitterNoBias lm;
		
		/**
		 * Sets the optimiser.
		 *
		 * @param lm the new lm
		 */
		public void setLM(LMDualCurveFitterNoBias lm) {this.lm=lm;}
		
		/**
		 * Sets the Rice sigma.
		 *
		 * @param sig the new sigma
		 */
		public void setSigma(double sig){
			sigma=sig;
		}

		
		/**
		 * Gets estimates of what would be the observed value in absence of Rice noise corruption
		 *
		 * @param x =[Tr,Te]
		 * @param a the parameters of the equation y=  [ (      a[0]* Math.exp(-(Te / a[3])) + a[1]* Math.exp(-(Te / a[4]) ) * (1 - Math.exp(-(Tr / a[2])) ]
		 * @return the non riced value
		 */
//		@Override
		public double getY(double[] x, double[] a) {
			double Tr=x[0];
			double Te=x[1];					
			return  getNonRicedY(Tr,Te,a) ; // a[1] - echo times
		}

		/**
		 * Gets estimates of what would be the observed value in absence of Rice noise corruption
		 *
		 * @param Tr the recovery time
		 * @param Te the echo time
		 * @param a the parameters of the equation y=  [ (      a[0]* Math.exp(-(Te / a[3])) + a[1]* Math.exp(-(Te / a[4]) ) * (1 - Math.exp(-(Tr / a[2])) ]
		 * @return the value in absence of any Rice corruption
		 */
		public double getNonRicedY(double Tr,double Te, double[] a) {
			return (1-Math.exp(-Tr/a[2])) * ( a[0]* Math.exp(-(Te / a[3])) + a[1] * Math.exp(-(Te / a[4])));
		}
		

		/**
		 * The equation being optimized.
		 *
		 * @return the string
		 */
		public String toString() {
			return "y=   [ (      a[0]* Math.exp(-(Te / a[3])) + a[1]* Math.exp(-(Te / a[4]) ) * (1 - Math.exp(-(Tr / a[2])) ]";
		}
		
		/**
		 * Gets the partial derivate of the signal according to one of the parameters
		 *
		 * @param x =[Tr,Te]
		 * @param a the parameters of the equation y=  [ (      a[0]* Math.exp(-(Te / a[3])) + a[1]* Math.exp(-(Te / a[4]) ) * (1 - Math.exp(-(Tr / a[2])) ]
		 * @param parameterIndex the index of the parameter in a
		 * @return the partial derivate value around this point
		 */
		@Override
		public double getPartialDerivate(double[] x, double[] a, int parameterIndex) {
			double Tr=x[0];
			double Te=x[1];					
			switch (parameterIndex) {
			case 0: return  (1-Math.exp(-Tr/a[2]))*(Math.exp(-Te/a[3]));  
			case 1: return  (1-Math.exp(-Tr/a[2]))*(Math.exp(-Te/a[4]));  
			case 2: return  (a[0]*Math.exp(-Te/a[3])+a[1]*Math.exp(-Te/a[4]))*Math.exp(-Tr/a[2])*Tr/(-a[2]*a[2]);
			case 3: return  a[0]*(1-Math.exp(-Tr/a[2]))*Math.exp(-Te/a[3])*Te/(a[3]*a[3]);
			case 4: return  (1-Math.exp(-Tr/a[2]))*Math.exp(-Te/a[4])*Te/(a[4]*a[4]);
			} 
			throw new RuntimeException("No such parameter index: " + parameterIndex);
		}
	} // end class

	/**
	 * The Class T1T2MultiBias. Describes cross-fitting of a T1 and T2 relaxation with a single T2 component, with an offset and without Rice noise
	 */
	//////////////////////  5-  T1T2MULTIBIAS  /////////////////////////////////////////////////////////////////////////////////////////
	public static class T1T2MultiBias extends LMAMultiDimFunction {
		
		/** The Constant Nparams. */
		public static final int Nparams=6;
		
		/** The sigma. */
		public double sigma=0;
		
		/** The lm. */
		public LMDualCurveFitterNoBias lm;
		
		/**
		 * Sets the optimiser
		 *
		 * @param lm the new lm
		 */
		public void setLM(LMDualCurveFitterNoBias lm) {this.lm=lm;}
		
		/**
		 * Sets the Rice sigma.
		 *
		 * @param sig the new sigma
		 */
		public void setSigma(double sig){
			sigma=sig;
		}

		
		/**
		 * Gets estimates of what would be the observed value in absence of Rice noise corruption
		 *
		 * @param x =[Tr,Te]
		 * @param a the parameters of the equation y=  [ (      a[0]* Math.exp(-(Te / a[3])) + a[1]* Math.exp(-(Te / a[4]) ) * (1 - Math.exp(-(Tr / a[2])) ]
		 * @return the non riced value
		 */
//		@Override
		public double getY(double[] x, double[] a) {
			double Tr=x[0];
			double Te=x[1];					
			return  getNonRicedY(Tr,Te,a) ; // a[1] - echo times
		}

		/**
		 * Gets estimates of what would be the observed value in absence of Rice noise corruption
		 *
		 * @param Tr the recovery time
		 * @param Te the echo time
		 * @param a the parameters of the equation y=  [ (      a[0]* Math.exp(-(Te / a[3])) + a[1]* Math.exp(-(Te / a[4]) ) * (1 - Math.exp(-(Tr / a[2])) ]
		 * @return the non riced value
		 */
		public double getNonRicedY(double Tr,double Te, double[] a) {
			return (1-Math.exp(-Tr/a[2])) * ( a[0]* Math.exp(-(Te / a[3])) + a[1] * Math.exp(-(Te / a[4])))+a[5];
		}
		

		/**
		 * Get the equation currently estimated
		 *
		 * @return the string
		 */
		public String toString() {
			return "y= [ (      a[0]* Math.exp(-(Te / a[3])) + a[1]* Math.exp(-(Te / a[4]) ) * (1 - Math.exp(-(Tr / a[2])) ]";
		}
		
		/**
		 * Gets the partial derivate of the signal according to one of the parameters
		 *
		 * @param x =[Tr,Te]
		 * @param a the parameters of the equation y=  [ (      a[0]* Math.exp(-(Te / a[3])) + a[1]* Math.exp(-(Te / a[4]) ) * (1 - Math.exp(-(Tr / a[2])) ]
		 * @param parameterIndex the index of the parameter in a
		 * @return the partial derivate value around this point
		 */
		@Override
		public double getPartialDerivate(double[] x, double[] a, int parameterIndex) {
			double Tr=x[0];
			double Te=x[1];					
			switch (parameterIndex) {
			case 0: return  (1-Math.exp(-Tr/a[2]))*(Math.exp(-Te/a[3]));  
			case 1: return  (1-Math.exp(-Tr/a[2]))*(Math.exp(-Te/a[4]));  
			case 2: return  (a[0]*Math.exp(-Te/a[3])+a[1]*Math.exp(-Te/a[4]))*Math.exp(-Tr/a[2])*Tr/(-a[2]*a[2]);
			case 3: return  a[0]*(1-Math.exp(-Tr/a[2]))*Math.exp(-Te/a[3])*Te/(a[3]*a[3]);
			case 4: return  a[1]*(1-Math.exp(-Tr/a[2]))*Math.exp(-Te/a[4])*Te/(a[4]*a[4]);
			case 5: return  1;
			} 
			throw new RuntimeException("No such parameter index: " + parameterIndex);
		}
	} // end class

	/**
	 * The Class T1T2MultiRice. Describes cross fitting of a T1-T2 relaxation with a single T2 component, without offset and with Rice noise
	 */
	//////////////////////  5-  T1T2MULTIRICE  /////////////////////////////////////////////////////////////////////////////////////////
	public static class T1T2MultiRice extends LMAMultiDimFunction {
		
		/** The Constant Nparams. */
		public static final int Nparams=5;
		
		/** The sigma. */
		public double sigma=0;
		
		/** The optimiser */
		public LMDualCurveFitterNoBias lm;
		
		/**
		 * Sets the lm.
		 *
		 * @param lm the new lm
		 */
		public void setLM(LMDualCurveFitterNoBias lm) {this.lm=lm;}
		
		/**
		 * Sets the Rice sigma.
		 *
		 * @param sig the new sigma
		 */
		public void setSigma(double sig){
			sigma=sig;
		}

		
		/**
		 * Gets estimates of what would be the observed value in absence of Rice noise corruption
		 *
		 * @param x =[Tr,Te]
		 * @param a the parameters of the equation y=  Rice_sigma [ (      a[0]* Math.exp(-(Te / a[3])) + a[1]* Math.exp(-(Te / a[4]) ) * (1 - Math.exp(-(Tr / a[2])) ]
		 * @return the value in absence of any Rice corruption
		 */
//		@Override
		public double getY(double[] x, double[] a) {
			double Tr=x[0];
			double Te=x[1];					
			return  MRUtils.besFunkCost( getNonRicedY(Tr,Te,a)  , sigma); // a[1] - echo times
		}

		/**
		 * Gets estimates of what would be the observed value in absence of Rice noise corruption
		 *
		 * @param Tr the recovery time
		 * @param Te the echo time
		 * @param a the parameters of the equation y=   Rice_sigma [ (      a[0]* Math.exp(-(Te / a[3])) + a[1]* Math.exp(-(Te / a[4]) ) * (1 - Math.exp(-(Tr / a[2])) ]
		 * @return the non riced value
		 */
		public double getNonRicedY(double Tr,double Te, double[] a) {
			return (1-Math.exp(-Tr/a[2])) * ( a[0]* Math.exp(-(Te / a[3])) + a[1] * Math.exp(-(Te / a[4])));
		}
		

		/**
		 * Get the equation currently estimated
		 *
		 * @return the string
		 */
		public String toString() {
			return "y=  Rice_sigma [ (      a[0]* Math.exp(-(Te / a[3])) + a[1]* Math.exp(-(Te / a[4]) ) * (1 - Math.exp(-(Tr / a[2])) ]";
		}
		
		/**
		 * Gets the partial derivate of the signal according to one of the parameters
		 *
		 * @param x =[Tr,Te]
		 * @param a the parameters of the equation y=  Rice_sigma [ (      a[0]* Math.exp(-(Te / a[3])) + a[1]* Math.exp(-(Te / a[4]) ) * (1 - Math.exp(-(Tr / a[2])) ]
		 * @param parameterIndex the index of the parameter in a
		 * @return the partial derivate value around this point
		 */
		@Override
		public double getPartialDerivate(double[] x, double[] a, int parameterIndex) {
			double Tr=x[0];
			double Te=x[1];					
			switch (parameterIndex) {
			case 0: return  MRUtils.besFunkDeriv( getNonRicedY(Tr,Te,a)  , sigma)*(1-Math.exp(-Tr/a[2]))*(Math.exp(-Te/a[3]));  
			case 1: return  MRUtils.besFunkDeriv( getNonRicedY(Tr,Te,a)  , sigma)*(1-Math.exp(-Tr/a[2]))*(Math.exp(-Te/a[4]));  
			case 2: return  MRUtils.besFunkDeriv( getNonRicedY(Tr,Te,a)  , sigma)*(a[0]*Math.exp(-Te/a[3])+a[1]*Math.exp(-Te/a[4]))*Math.exp(-Tr/a[2])*Tr/(-a[2]*a[2]);
			case 3: return  MRUtils.besFunkDeriv( getNonRicedY(Tr,Te,a)  , sigma)*a[0]*(1-Math.exp(-Tr/a[2]))*Math.exp(-Te/a[3])*Te/(a[3]*a[3]);
			case 4: return  MRUtils.besFunkDeriv( getNonRicedY(Tr,Te,a)  , sigma)*a[1]*(1-Math.exp(-Tr/a[2]))*Math.exp(-Te/a[4])*Te/(a[4]*a[4]);
			} 
			throw new RuntimeException("No such parameter index: " + parameterIndex);
		}
	} // end class


	
	/**
	 * The Class T1T2DefaultMonoRice.  Describes cross-fitting of a T1 and T2 relaxation with a single T2 component, without offset and without Rice noise
	 */
	//////////////////////  6-  T1T2DEFAULTMONORICE  /////////////////////////////////////////////////////////////////////////////////////////
	public static class T1T2DefaultMonoRice extends LMAMultiDimFunction {
		
		/** The Constant Nparams. */
		public static final int Nparams=4;
		
		/** The sigma. */
		public double sigma=0;
		
		/** The lm. */
		public LMDualCurveFitterNoBias lm;
		
		/**
		 * Sets the optimiser
		 *
		 * @param lm the new lm
		 */
		public void setLM(LMDualCurveFitterNoBias lm) {this.lm=lm;}
		
		/**
		 * Sets the sigma.
		 *
		 * @param sig the new sigma
		 */
		public void setSigma(double sig){
			sigma=sig;
		}

		/**
		 * Gets the non riced Y.
		 *
		 * @param Tr the tr
		 * @param Te the te
		 * @param a the a
		 * @return the non riced Y
		 */
		public double getNonRicedY(double Tr,double Te, double[] a) {
			return a[0]* Math.exp(-((Te+(Tr>9999 ? a[3] : 0)) / a[2]))*(1-Math.exp(-Tr/a[1]));
		}
		
		/**
		 * To string.
		 *
		 * @return the string
		 */
		public String toString() {
			return "y=  Rice_sigma [ a[0]* Math.exp(-((Te+a[3]) / a[2]))* (1 - Math.exp(-(Tr / a[1])) ]";
		}

		/**
		 * Gets the y.
		 *
		 * @param x the x
		 * @param a the a
		 * @return the y
		 */
		@Override
		public double getY(double[] x, double[] a) {
			double Tr=x[0];
			double Te=x[1];					
			return  MRUtils.besFunkCost( getNonRicedY(Tr,Te,a) , sigma); // a[1] - echo times
		}
		
		/**
		 * Gets the partial derivate.
		 *
		 * @param x the x
		 * @param a the a
		 * @param parameterIndex the parameter index
		 * @return the partial derivate
		 */
		@Override
		public double getPartialDerivate(double[] x, double[] a, int parameterIndex) {
			double Tr=x[0];
			double Te=x[1];					
			switch (parameterIndex) {
			case 0: return  MRUtils.besFunkDeriv( getNonRicedY(Tr,Te,a) , sigma)*(1-Math.exp(-Tr/a[1]))*(Math.exp(-(Te+a[3])/a[2]));  
			case 1: return  MRUtils.besFunkDeriv( getNonRicedY(Tr,Te,a) , sigma)*a[0]*Math.exp(-(Te+(Tr>9999 ? a[3] : 0))/a[2])*Math.exp(-Tr/a[1])*Tr/(-a[1]*a[1]);
			case 2: return  MRUtils.besFunkDeriv( getNonRicedY(Tr,Te,a) , sigma)*a[0]*(1-Math.exp(-Tr/a[1]))*Math.exp(-(Te+(Tr>9999 ? a[3] : 0))/a[2])*(Te+(Tr>9999 ? a[3] : 0))/(a[2]*a[2]);
			case 3: return  (Tr<10000 ? 0 : MRUtils.besFunkDeriv( getNonRicedY(Tr,Te,a) , sigma)*a[0]*(1-Math.exp(-Tr/a[1]))*Math.exp(-(Te+(Tr>9999 ? a[3] : 0))/a[2])*(-1/a[2]));
			} 
			throw new RuntimeException("No such parameter index: " + parameterIndex);
		}

	} // end class

	
	/**
	 * The Class T1T2DefaultMultiRice. Describes cross-fitting of a T1 and T2 relaxation with two T2 components, without offset and without Rice noise
	 */
	//////////////////////  7-  T1T2DEFAULTMULTIRICE  /////////////////////////////////////////////////////////////////////////////////////////
	public static class T1T2DefaultMultiRice extends LMAMultiDimFunction {
		
		/** The Constant Nparams. */
		public static final int Nparams=6;
		
		/** The sigma. */
		public double sigma=0;
		
		/** The lm. */
		public LMDualCurveFitterNoBias lm;
		
		/**
		 * Sets the optimiser
		 *
		 * @param lm the new lm
		 */
		public void setLM(LMDualCurveFitterNoBias lm) {this.lm=lm;}
		
		/**
		 * Sets the Rice sigma.
		 *
		 * @param sig the new sigma
		 */
		public void setSigma(double sig){
			sigma=sig;
		}

		
		/**
		 * Gets estimates of what would be the observed value in presence of Rice noise corruption
		 *
		 * @param x the vector of parameters {Tr, Te}
		 * @param a the parameters of the equation y= Rice_sigma [ (      a[0]* Math.exp(-((Te+a[5]) / a[3])) + a[1]* Math.exp(-((Te+a[5]) / a[4]) ) * (1 - Math.exp(-(Tr / a[2])) ]
		 * @return the non riced value
		 */
//		@Override
		public double getY(double[] x, double[] a) {
			double Tr=x[0];
			double Te=x[1];					
			return  MRUtils.besFunkCost( getNonRicedY(Tr,Te,a)  , sigma); // a[1] - echo times
		}

		/**
		 * Gets estimates of what would be the observed value in absence of Rice noise corruption
		 *
		 * @param Tr the recovery time
		 * @param Te the echo time
		 * @param a the parameters of the equation y=  Rice_sigma [ (      a[0]* Math.exp(-((Te+a[5]) / a[3])) + a[1]* Math.exp(-((Te+a[5]) / a[4]) ) * (1 - Math.exp(-(Tr / a[2])) ]
		 * @return the value in absence of any Rice corruption
		 */
		public double getNonRicedY(double Tr,double Te, double[] a) {
			return (1-Math.exp(-Tr/a[2])) * ( a[0]* Math.exp(-((Te+(Tr>9999 ? a[5] : 0)) / a[3])) + a[1] * Math.exp(-((Te+(Tr>9999 ? a[5] : 0)) / a[4])));
		}
		

		/**
		 * The equation being optimized.
		 *
		 * @return the string
		 */
		public String toString() {
			return "y=  Rice_sigma [ (      a[0]* Math.exp(-((Te+a[5]) / a[3])) + a[1]* Math.exp(-((Te+a[5]) / a[4]) ) * (1 - Math.exp(-(Tr / a[2])) ]";
		}
		
		/**
		 * Gets the partial derivate of the signal according to one of the parameters
		 *
		 * @param x =[Tr,Te]
		 * @param a the parameters of the equation y=  Rice_sigma [ (      a[0]* Math.exp(-((Te+a[5]) / a[3])) + a[1]* Math.exp(-((Te+a[5]) / a[4]) ) * (1 - Math.exp(-(Tr / a[2])) ]
		 * @param parameterIndex the index of the parameter in a
		 * @return the partial derivate value around this point
		 */
		@Override
		public double getPartialDerivate(double[] x, double[] a, int parameterIndex) {
			double Tr=x[0];
			double Te=x[1];					
			switch (parameterIndex) {
			case 0: return  MRUtils.besFunkDeriv( getNonRicedY(Tr,Te,a)  , sigma)*(1-Math.exp(-Tr/a[2]))*(Math.exp(-(Te+(Tr>9999 ? a[5] : 0))/a[3]));  
			case 1: return  MRUtils.besFunkDeriv( getNonRicedY(Tr,Te,a)  , sigma)*(1-Math.exp(-Tr/a[2]))*(Math.exp(-(Te+(Tr>9999 ? a[5] : 0))/a[4]));  
			case 2: return  MRUtils.besFunkDeriv( getNonRicedY(Tr,Te,a)  , sigma)*(a[0]*Math.exp(-(Te+(Tr>9999 ? a[5] : 0))/a[3])+a[1]*Math.exp(-(Te+(Tr>9999 ? a[5] : 0))/a[4]))*Math.exp(-Tr/a[2])*Tr/(-a[2]*a[2]);
			case 3: return  MRUtils.besFunkDeriv( getNonRicedY(Tr,Te,a)  , sigma)*a[0]*(1-Math.exp(-Tr/a[2]))*Math.exp(-(Te+(Tr>9999 ? a[5] : 0))/a[3])*(Te+(Tr>9999 ? a[5] : 0))/(a[3]*a[3]);
			case 4: return  MRUtils.besFunkDeriv( getNonRicedY(Tr,Te,a)  , sigma)*a[1]*(1-Math.exp(-Tr/a[2]))*Math.exp(-(Te+(Tr>9999 ? a[5] : 0))/a[4])*(Te+(Tr>9999 ? a[5] : 0))/(a[4]*a[4]);
			case 5: return  (Tr<100000 ? 0 : MRUtils.besFunkDeriv( getNonRicedY(Tr,Te,a) , sigma)*(1-Math.exp(-Tr/a[1]))*( a[0]* Math.exp(-(Te+(Tr>9999 ? a[5] : 0))/a[3])*(-1/a[3]) + a[1]* Math.exp(-(Te+(Tr>9999 ? a[5] : 0))/a[4])*(-1/a[4]) )  ) ;
			} 
			throw new RuntimeException("No such parameter index: " + parameterIndex);
		}
	} // end class


}