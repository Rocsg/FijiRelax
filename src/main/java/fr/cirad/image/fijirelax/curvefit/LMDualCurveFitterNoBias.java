package fr.cirad.image.fijirelax.curvefit;
import fr.cirad.image.common.VitimageUtils;

import fr.cirad.image.fijirelax.mrialgo.MRUtils;
import ij.IJ;
import lma.LMAFunction;
import lma.LMAMatrix.InvertException;
import lma.LMAMultiDimFunction;
import lma.implementations.LMA;

/** Class for building curve fitter */
public class LMDualCurveFitterNoBias {	
	public boolean debugLM=false;
	public final static String[] timeunits={"ms", "s"};
	public final static int[] timeitems={MRUtils.MSEC, MRUtils.SEC};
	public final static String[] fititems2={"Simplex","Levenberg-Marquardt"};
	public final static int[] constitems2={MRUtils.SIMPLEX,MRUtils.LM};
	protected int fit;                // Number of curve type to fit

	protected double[][] data; // x,y data to fit
	private MRLMA lma;
	private double[] parameters;
	private float[][] parametersBoundaries;
	
	private float gfit;
	public boolean iterationsBreak=false;
	public static double lambda = 0.01;
	public static double minDeltaChi2 = 1e-7;
	public static int maxIterations=200;

		
	
    public float[][] getDefaultParametersBoundaries(int fitType) {
    	float[][]parametersBoundaries=new float[MRUtils.getNparams(fitType)+1][2];
    	for(int i=0;i<parametersBoundaries.length;i++)parametersBoundaries[i]=new float[] {(float) -MRUtils.infinity,(float) MRUtils.infinity};
    	return parametersBoundaries;
    }


	/** Construct a new CurveFitter. */
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


	public void configLMA(double lambda, double minDeltaChi2, int maxIterations) {
		if (lma==null) return;
		lma.lambda=lambda;
		lma.minDeltaChi2=minDeltaChi2;
		lma.maxIterations=maxIterations;
	}

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

	public void doFit(double lambda, double minDeltaChi2, int maxIterations) {
		try {
			lma.fit( lambda,  minDeltaChi2,  maxIterations);
			parameters=lma.parameters;
			gfit=lma.chi2Goodness();
		} catch (InvertException e) {
			for(int i=0;i<parameters.length;i++)parameters[i]=MRUtils.ERROR_VALUE;}
	}

	public double[] getParams() {
		double[]ret=new double[parameters.length+1];
		for(int i=0;i<parameters.length;i++)ret[i]=parameters[i];
		ret[parameters.length]=lma.chi2Goodness();
		return ret;
	}

	public LMA getLMA() {
		return lma;
	}


	
	
	
	
	

	//////////////////////  1-  T1MONO  /////////////////////////////////////////////////////////////////////////////////////////
	public static class T1Mono extends LMAMultiDimFunction {
		public static final int Nparams=2;
		public double sigma=0;
		public LMDualCurveFitterNoBias lm;
		public void setLM(LMDualCurveFitterNoBias lm) {this.lm=lm;}
		public void setSigma(double sig){
			sigma=sig;
		}

		public double getNonRicedY(double Tr,double Te, double[] a) {
			return a[0]*(1-Math.exp(-Tr/a[1]));
		}
		
		public String toString() {
			return "y=  Rice_sigma [ a[0]* (1 - Math.exp(-(Tr / a[1])) ]";
		}

		@Override
		public double getY(double[] x, double[] a) {
			double Tr=x[0];
			double Te=x[1];					
			return  getNonRicedY(Tr,Te,a) ; // a[1] - echo times
		}
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
	//////////////////////  1-  T1MONOBIAS  /////////////////////////////////////////////////////////////////////////////////////////
	public static class T1MonoBias extends LMAMultiDimFunction {
		public static final int Nparams=3;
		public double sigma=0;
		public LMDualCurveFitterNoBias lm;
		public void setLM(LMDualCurveFitterNoBias lm) {this.lm=lm;}
		public void setSigma(double sig){
			sigma=sig;
		}

		public double getNonRicedY(double Tr,double Te, double[] a) {
			return a[0]*(1-Math.exp(-Tr/a[1]))+a[2];
		}
		
		public String toString() {
			return "y=  Rice_sigma [ a[0]* (1 - Math.exp(-(Tr / a[1])) ]";
		}

		@Override
		public double getY(double[] x, double[] a) {
			double Tr=x[0];
			double Te=x[1];					
			return  getNonRicedY(Tr,Te,a); // a[1] - echo times
		}
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
	//////////////////////  1-  T1MONORICE  /////////////////////////////////////////////////////////////////////////////////////////
	public static class T1MonoRice extends LMAMultiDimFunction {
		public static final int Nparams=2;
		public double sigma=0;
		public LMDualCurveFitterNoBias lm;
		public void setLM(LMDualCurveFitterNoBias lm) {this.lm=lm;}
		public void setSigma(double sig){
			sigma=sig;
		}

		public double getNonRicedY(double Tr,double Te, double[] a) {
			return a[0]*(1-Math.exp(-Tr/a[1]));
		}
		
		public String toString() {
			return "y=  Rice_sigma [ a[0]* (1 - Math.exp(-(Tr / a[1])) ]";
		}

		@Override
		public double getY(double[] x, double[] a) {
			double Tr=x[0];
			double Te=x[1];					
			return  MRUtils.besFunkCost( getNonRicedY(Tr,Te,a) , sigma); // a[1] - echo times
		}
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
	
	

	
	
	
	
	//////////////////////  2-  T2MONO  /////////////////////////////////////////////////////////////////////////////////////////
	public static class T2Mono extends LMAMultiDimFunction {
		public static final int Nparams=2;
		public double sigma=0;
		public LMDualCurveFitterNoBias lm;
		public void setLM(LMDualCurveFitterNoBias lm) {this.lm=lm;}
		public void setSigma(double sig){
			sigma=sig;
		}

		public double getNonRicedY(double Tr,double Te, double[] a) {
			return a[0]* Math.exp(-(Te / a[1]));
		}
		
		public String toString() {
			return "y=  Rice_sigma [ a[0]* Math.exp(-(Te / a[1])) ]";
		}

		@Override
		public double getY(double[] x, double[] a) {
			double Tr=x[0];
			double Te=x[1];					
			return  getNonRicedY(Tr,Te,a); // a[1] - echo times
		}
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
	//////////////////////  2-  T2MONOBIAS  /////////////////////////////////////////////////////////////////////////////////////////
	public static class T2MonoBias extends LMAMultiDimFunction {
		public static final int Nparams=3;
		public double sigma=0;
		public LMDualCurveFitterNoBias lm;
		public void setLM(LMDualCurveFitterNoBias lm) {this.lm=lm;}
		public void setSigma(double sig){
			sigma=sig;
		}

		public double getNonRicedY(double Tr,double Te, double[] a) {
			return a[0]* Math.exp(-(Te / a[1]))+a[2];
		}
		
		public String toString() {
			return "y=  Rice_sigma [ a[0]* Math.exp(-(Te / a[1])) ]";
		}

		@Override
		public double getY(double[] x, double[] a) {
			double Tr=x[0];
			double Te=x[1];					
			return  getNonRicedY(Tr,Te,a) ; // a[1] - echo times
		}
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
	//////////////////////  2-  T2MONORICE  /////////////////////////////////////////////////////////////////////////////////////////
	public static class T2MonoRice extends LMAMultiDimFunction {
		public static final int Nparams=2;
		public double sigma=0;
		public LMDualCurveFitterNoBias lm;
		public void setLM(LMDualCurveFitterNoBias lm) {this.lm=lm;}
		public void setSigma(double sig){
			sigma=sig;
		}

		public double getNonRicedY(double Tr,double Te, double[] a) {
			return a[0]* Math.exp(-(Te / a[1]));
		}
		
		public String toString() {
			return "y=  Rice_sigma [ a[0]* Math.exp(-(Te / a[1])) ]";
		}

		@Override
		public double getY(double[] x, double[] a) {
			double Tr=x[0];
			double Te=x[1];					
			return  MRUtils.besFunkCost( getNonRicedY(Tr,Te,a) , sigma); // a[1] - echo times
		}
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


	
	
	
	
	
	
	
	
	
	
	//////////////////////  3-  T2MULTI  /////////////////////////////////////////////////////////////////////////////////////////
	public static class T2Multi extends LMAMultiDimFunction {
		public static final int Nparams=4;
		public double sigma=0;
		public LMDualCurveFitterNoBias lm;
		public void setLM(LMDualCurveFitterNoBias lm) {this.lm=lm;}
		public void setSigma(double sig){
			sigma=sig;
		}

		
//		@Override
		public double getY(double[] x, double[] a) {
			double Tr=x[0];
			double Te=x[1];					
			return  getNonRicedY(Tr,Te,a); // a[1] - echo times
		}

		public double getNonRicedY(double Tr,double Te, double[] a) {
			return  ( a[0]* Math.exp(-(Te / a[2])) + a[1] * Math.exp(-(Te / a[3])));
		}
		

		public String toString() {
			return "y=  Rice_sigma [ a[0]* Math.exp(-(Te / a[2]))  +  a[1]* Math.exp(-(Te / a[3])) ]";
		}
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
	//////////////////////  3-  T2MULTIRICE  /////////////////////////////////////////////////////////////////////////////////////////
	public static class T2MultiBias extends LMAMultiDimFunction {
		public static final int Nparams=5;
		public double sigma=0;
		public LMDualCurveFitterNoBias lm;
		public void setLM(LMDualCurveFitterNoBias lm) {this.lm=lm;}
		public void setSigma(double sig){
			sigma=sig;
		}

		
//		@Override
		public double getY(double[] x, double[] a) {
			double Tr=x[0];
			double Te=x[1];					
			return  getNonRicedY(Tr,Te,a); // a[1] - echo times
		}

		public double getNonRicedY(double Tr,double Te, double[] a) {
			return  ( a[0]* Math.exp(-(Te / a[2])) + a[1] * Math.exp(-(Te / a[3]))+a[4]);
		}
		

		public String toString() {
			return "y=  Rice_sigma [ a[0]* Math.exp(-(Te / a[2]))  +  a[1]* Math.exp(-(Te / a[3])) ]";
		}
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
	//////////////////////  3-  T2MULTIRICE  /////////////////////////////////////////////////////////////////////////////////////////
	public static class T2MultiRice extends LMAMultiDimFunction {
		public static final int Nparams=4;
		public double sigma=0;
		public LMDualCurveFitterNoBias lm;
		public void setLM(LMDualCurveFitterNoBias lm) {this.lm=lm;}
		public void setSigma(double sig){
			sigma=sig;
		}

		
//		@Override
		public double getY(double[] x, double[] a) {
			double Tr=x[0];
			double Te=x[1];					
			return  MRUtils.besFunkCost( getNonRicedY(Tr,Te,a)  , sigma); // a[1] - echo times
		}

		public double getNonRicedY(double Tr,double Te, double[] a) {
			return  ( a[0]* Math.exp(-(Te / a[2])) + a[1] * Math.exp(-(Te / a[3])));
		}
		

		public String toString() {
			return "y=  Rice_sigma [ a[0]* Math.exp(-(Te / a[2]))  +  a[1]* Math.exp(-(Te / a[3])) ]";
		}
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

	
	
	
	
	
	
	
	//////////////////////  4-  T1T2MONO  /////////////////////////////////////////////////////////////////////////////////////////
	public static class T1T2Mono extends LMAMultiDimFunction {
		public static final int Nparams=3;
		public double sigma=0;
		public LMDualCurveFitterNoBias lm;
		public void setLM(LMDualCurveFitterNoBias lm) {this.lm=lm;}
		public void setSigma(double sig){
			sigma=sig;
		}

		public double getNonRicedY(double Tr,double Te, double[] a) {
			return a[0]* Math.exp(-(Te / a[2]))*(1-Math.exp(-Tr/a[1]));
		}
		
		public String toString() {
			return "y=  Rice_sigma [ a[0]* Math.exp(-(Te / a[2]))* (1 - Math.exp(-(Tr / a[1])) ]";
		}

		@Override
		public double getY(double[] x, double[] a) {
			double Tr=x[0];
			double Te=x[1];					
			return   getNonRicedY(Tr,Te,a); // a[1] - echo times
		}
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
	//////////////////////  4-  T1T2MONOBIAS  /////////////////////////////////////////////////////////////////////////////////////////
	public static class T1T2MonoBias extends LMAMultiDimFunction {
		public static final int Nparams=4;
		public double sigma=0;
		public LMDualCurveFitterNoBias lm;
		public void setLM(LMDualCurveFitterNoBias lm) {this.lm=lm;}
		public void setSigma(double sig){
			sigma=sig;
		}

		public double getNonRicedY(double Tr,double Te, double[] a) {
			return a[0]* Math.exp(-(Te / a[2]))*(1-Math.exp(-Tr/a[1]))+a[3];
		}
		
		public String toString() {
			return "y=  Rice_sigma [ a[0]* Math.exp(-(Te / a[2]))* (1 - Math.exp(-(Tr / a[1])) ]";
		}

		@Override
		public double getY(double[] x, double[] a) {
			double Tr=x[0];
			double Te=x[1];					
			return  getNonRicedY(Tr,Te,a); // a[1] - echo times
		}
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
	//////////////////////  4-  T1T2MONORICE  /////////////////////////////////////////////////////////////////////////////////////////
	public static class T1T2MonoRice extends LMAMultiDimFunction {
		public static final int Nparams=3;
		public double sigma=0;
		public LMDualCurveFitterNoBias lm;
		public void setLM(LMDualCurveFitterNoBias lm) {this.lm=lm;}
		public void setSigma(double sig){
			sigma=sig;
		}

		public double getNonRicedY(double Tr,double Te, double[] a) {
			return a[0]* Math.exp(-(Te / a[2]))*(1-Math.exp(-Tr/a[1]));
		}
		
		public String toString() {
			return "y=  Rice_sigma [ a[0]* Math.exp(-(Te / a[2]))* (1 - Math.exp(-(Tr / a[1])) ]";
		}

		@Override
		public double getY(double[] x, double[] a) {
			double Tr=x[0];
			double Te=x[1];					
			return  MRUtils.besFunkCost( getNonRicedY(Tr,Te,a) , sigma); // a[1] - echo times
		}
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

	

	
	
	
	//////////////////////  5-  T1T2MULTI  /////////////////////////////////////////////////////////////////////////////////////////
	public static class T1T2Multi extends LMAMultiDimFunction {
		public static final int Nparams=5;
		public double sigma=0;
		public LMDualCurveFitterNoBias lm;
		public void setLM(LMDualCurveFitterNoBias lm) {this.lm=lm;}
		public void setSigma(double sig){
			sigma=sig;
		}

		
//		@Override
		public double getY(double[] x, double[] a) {
			double Tr=x[0];
			double Te=x[1];					
			return  getNonRicedY(Tr,Te,a) ; // a[1] - echo times
		}

		public double getNonRicedY(double Tr,double Te, double[] a) {
			return (1-Math.exp(-Tr/a[2])) * ( a[0]* Math.exp(-(Te / a[3])) + a[1] * Math.exp(-(Te / a[4])));
		}
		

		public String toString() {
			return "y=  Rice_sigma [ (      a[0]* Math.exp(-(Te / a[3])) + a[1]* Math.exp(-(Te / a[4]) ) * (1 - Math.exp(-(Tr / a[2])) ]";
		}
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

	//////////////////////  5-  T1T2MULTIBIAS  /////////////////////////////////////////////////////////////////////////////////////////
	public static class T1T2MultiBias extends LMAMultiDimFunction {
		public static final int Nparams=6;
		public double sigma=0;
		public LMDualCurveFitterNoBias lm;
		public void setLM(LMDualCurveFitterNoBias lm) {this.lm=lm;}
		public void setSigma(double sig){
			sigma=sig;
		}

		
//		@Override
		public double getY(double[] x, double[] a) {
			double Tr=x[0];
			double Te=x[1];					
			return  getNonRicedY(Tr,Te,a) ; // a[1] - echo times
		}

		public double getNonRicedY(double Tr,double Te, double[] a) {
			return (1-Math.exp(-Tr/a[2])) * ( a[0]* Math.exp(-(Te / a[3])) + a[1] * Math.exp(-(Te / a[4])))+a[5];
		}
		

		public String toString() {
			return "y=  Rice_sigma [ (      a[0]* Math.exp(-(Te / a[3])) + a[1]* Math.exp(-(Te / a[4]) ) * (1 - Math.exp(-(Tr / a[2])) ]";
		}
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

	//////////////////////  5-  T1T2MULTIRICE  /////////////////////////////////////////////////////////////////////////////////////////
	public static class T1T2MultiRice extends LMAMultiDimFunction {
		public static final int Nparams=5;
		public double sigma=0;
		public LMDualCurveFitterNoBias lm;
		public void setLM(LMDualCurveFitterNoBias lm) {this.lm=lm;}
		public void setSigma(double sig){
			sigma=sig;
		}

		
//		@Override
		public double getY(double[] x, double[] a) {
			double Tr=x[0];
			double Te=x[1];					
			return  MRUtils.besFunkCost( getNonRicedY(Tr,Te,a)  , sigma); // a[1] - echo times
		}

		public double getNonRicedY(double Tr,double Te, double[] a) {
			return (1-Math.exp(-Tr/a[2])) * ( a[0]* Math.exp(-(Te / a[3])) + a[1] * Math.exp(-(Te / a[4])));
		}
		

		public String toString() {
			return "y=  Rice_sigma [ (      a[0]* Math.exp(-(Te / a[3])) + a[1]* Math.exp(-(Te / a[4]) ) * (1 - Math.exp(-(Tr / a[2])) ]";
		}
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


	
	//////////////////////  6-  T1T2DEFAULTMONORICE  /////////////////////////////////////////////////////////////////////////////////////////
	public static class T1T2DefaultMonoRice extends LMAMultiDimFunction {
		public static final int Nparams=4;
		public double sigma=0;
		public LMDualCurveFitterNoBias lm;
		public void setLM(LMDualCurveFitterNoBias lm) {this.lm=lm;}
		public void setSigma(double sig){
			sigma=sig;
		}

		public double getNonRicedY(double Tr,double Te, double[] a) {
			return a[0]* Math.exp(-((Te+(Tr>9999 ? a[3] : 0)) / a[2]))*(1-Math.exp(-Tr/a[1]));
		}
		
		public String toString() {
			return "y=  Rice_sigma [ a[0]* Math.exp(-((Te+a[3]) / a[2]))* (1 - Math.exp(-(Tr / a[1])) ]";
		}

		@Override
		public double getY(double[] x, double[] a) {
			double Tr=x[0];
			double Te=x[1];					
			return  MRUtils.besFunkCost( getNonRicedY(Tr,Te,a) , sigma); // a[1] - echo times
		}
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

	
	//////////////////////  7-  T1T2DEFAULTMULTIRICE  /////////////////////////////////////////////////////////////////////////////////////////
	public static class T1T2DefaultMultiRice extends LMAMultiDimFunction {
		public static final int Nparams=6;
		public double sigma=0;
		public LMDualCurveFitterNoBias lm;
		public void setLM(LMDualCurveFitterNoBias lm) {this.lm=lm;}
		public void setSigma(double sig){
			sigma=sig;
		}

		
//		@Override
		public double getY(double[] x, double[] a) {
			double Tr=x[0];
			double Te=x[1];					
			return  MRUtils.besFunkCost( getNonRicedY(Tr,Te,a)  , sigma); // a[1] - echo times
		}

		public double getNonRicedY(double Tr,double Te, double[] a) {
			return (1-Math.exp(-Tr/a[2])) * ( a[0]* Math.exp(-((Te+(Tr>9999 ? a[5] : 0)) / a[3])) + a[1] * Math.exp(-((Te+(Tr>9999 ? a[5] : 0)) / a[4])));
		}
		

		public String toString() {
			return "y=  Rice_sigma [ (      a[0]* Math.exp(-((Te+a[5]) / a[3])) + a[1]* Math.exp(-((Te+a[5]) / a[4]) ) * (1 - Math.exp(-(Tr / a[2])) ]";
		}
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