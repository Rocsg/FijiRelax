/*
 * 
 */
package io.github.rocsg.fijirelax.curvefit;
import io.github.rocsg.fijiyama.common.VitimageUtils;
import io.github.rocsg.fijirelax.mrialgo.MRUtils;
import io.github.rocsg.fijirelax.lma.LMAFunction;
import io.github.rocsg.fijirelax.lma.LMAMatrix.InvertException;
import io.github.rocsg.fijirelax.lma.LMA;

// TODO: Auto-generated Javadoc
/**
 * The Class LMCurveFitterNoBias.
 */
public class LMCurveFitterNoBias {
	
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
	
	/** The lambda. */
	public static double lambda = 0.001;
	
	/** The min delta chi 2. */
	public static double minDeltaChi2 = 1e-30;
	
	/** The max iterations. */
	public static int maxIterations=200;

	
	
	 
	//TODO : after test set, remove from there
	/**
	 *  
	 * Return a first order approximation of initial signal leading to value valFunk, after being corrupted by Rice noise with stddev=sigma .
	 *
	 * @param valFunk The observed value after noise alteration
	 * @param sigma The parameter of the Rice noise
	 * @return First order approximation of the signal, as a double
	 */
	static double sigmaWay(double valFunk,double sigma){
		double ret=valFunk*valFunk-2*sigma*sigma;
		return (ret<0 ? 0 : Math.sqrt(ret));
	}
	    
	/**
	 * Bes funk cost.
	 *
	 * @param d the d
	 * @param sigma2 the sigma 2
	 * @return the double
	 */
	static double besFunkCost(double d,double sigma2) {
		double alpha=d*d/(4*sigma2*sigma2);
		return (double)(Math.sqrt(Math.PI*sigma2*sigma2/2.0)*( (1+2*alpha)*bessi0NoExp(alpha) + 2*alpha*bessi1NoExp(alpha) ));
	}

	/**
	 * Bes funk deriv.
	 *
	 * @param d the d
	 * @param sigma the sigma
	 * @return the double
	 */
	static double besFunkDeriv(double d,double sigma) {
		double alpha=d*d/(4*sigma*sigma);
		return (double)(Math.sqrt(Math.PI*alpha/2.0)*(bessi0NoExp(alpha) + bessi1NoExp(alpha) ));
	}
 
	/**
	 * Bessi 0 no exp.
	 *
	 * @param alpha the alpha
	 * @return the double
	 */
	static double bessi0NoExp( double alpha ){
		double x=(double) alpha;
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
	   return (double)ans;
	}

	/**
	 * Bessi 1 no exp.
	 *
	 * @param alpha the alpha
	 * @return the double
	 */
	static double bessi1NoExp( double alpha){
			double ax,ans;
			double y;
			double x=(double)alpha;
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
		   return ((double)(x < 0.0 ? -ans : ans));
		}
		
		
	
	
	
	

	/**
	 *  Construct a new CurveFitter.
	 *
	 * @param xData the x data
	 * @param yData the y data
	 * @param fitType the fit type
	 * @param sigma the sigma
	 * @param debugLM the debug LM
	 */
	public LMCurveFitterNoBias (double[] xData, double[] yData, int fitType,double sigma,boolean debugLM) {
		this.fit=fitType;
		this.debugLM=debugLM;
		int xlength=xData.length;
		int ylength=yData.length;
		if (xlength!= ylength)
			throw new IllegalArgumentException("Arrays' lengths do not match");

		data=new double[2][xlength];
		data[0]=xData;
		data[1]=yData;
        double minT=VitimageUtils.min(xData);
        double maxT=VitimageUtils.max(xData);
        double medT=minT/2+maxT/2;
        double maxMag=VitimageUtils.max(yData);
        double medMag=maxMag/2;


        int np=getNumParams(fitType);
		double[] inparameters=new double[np];
		parameters=new double[np];
		parametersBoundaries=new float[np][2];

		switch (fit) {
		case MRUtils.T1_MONO_RICE:
			inparameters[0] = maxMag;
			inparameters[1] = medT;
			if(MRUtils.useBoundaries) {
				parametersBoundaries[0]=new float[] {(float) MRUtils.epsilon,(float) (maxMag*MRUtils.factorT1M0MaxRatio)};
				parametersBoundaries[1]=new float[] {(float) (minT*MRUtils.factorT1MinRatio),(float) (maxT*MRUtils.factorT1MaxRatio)};
			}
			else parametersBoundaries=new float[][] {{(float) -MRUtils.infinity,(float) MRUtils.infinity},{(float) -MRUtils.infinity,(float) MRUtils.infinity}};
			T1RecoveryRice t1rr=new T1RecoveryRice();
			t1rr.setSigma(sigma);
			t1rr.setLM(this);
			lma = new MRLMA(t1rr,inparameters,data,parametersBoundaries);
			break;

		case MRUtils.T2_MONO_RICE:
			inparameters[0] = maxMag;
			inparameters[1] = medT;//tHope;
			if(MRUtils.useBoundaries) {
				parametersBoundaries[0]=new float[] {(float) MRUtils.epsilon,(float) (maxMag*MRUtils.factorT2M0MaxRatio)};
				parametersBoundaries[1]=new float[] {(float) (minT*MRUtils.factorT2MinRatio),(float) (maxT*MRUtils.factorT2MaxRatio)};
			}
			else parametersBoundaries=new float[][] {{(float) -MRUtils.infinity,(float) MRUtils.infinity},{(float) -MRUtils.infinity,(float) MRUtils.infinity}};
			T2RelaxRice t2rr=new T2RelaxRice();
			t2rr.setSigma(sigma);
			t2rr.setLM(this);
			lma = new MRLMA(t2rr,inparameters,data,parametersBoundaries);
			break;

		case MRUtils.T2_MULTI_RICE:
			inparameters[0] = medMag+MRUtils.epsilon;
			inparameters[1] = minT+MRUtils.epsilon;//tHope;
			inparameters[2] = medMag-MRUtils.epsilon;
			inparameters[3] = maxT-MRUtils.epsilon;//tHope2;
			if(MRUtils.useBoundaries) {
				parametersBoundaries[0]=new float[] {(float) MRUtils.epsilon,(float) ((maxMag)*MRUtils.factorT2M0MaxRatio)};
				parametersBoundaries[1]=new float[] {(float) (minT*MRUtils.factorT2MinRatio),(float) (maxT*MRUtils.factorT2MaxRatio)};
				parametersBoundaries[2]=new float[] {(float) MRUtils.epsilon,(float) ((maxMag)*MRUtils.factorT2M0MaxRatio)};
				parametersBoundaries[3]=new float[] {(float) (minT*MRUtils.factorT2MinRatio),(float) (maxT*MRUtils.factorT2MaxRatio)};
			}
			else parametersBoundaries=new float[][] {{(float) -MRUtils.infinity,(float) MRUtils.infinity},{(float) -MRUtils.infinity,(float) MRUtils.infinity},{(float) -MRUtils.infinity,(float) MRUtils.infinity},{(float) -MRUtils.infinity,(float) MRUtils.infinity}};
			T2MulticompRice t2mcr=new T2MulticompRice();
			t2mcr.setSigma(sigma);
			t2mcr.setLM(this);
			lma = new MRLMA(t2mcr,inparameters,data,parametersBoundaries);
			break;		
		}
	}


	/**
	 * Config LMA.
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
	 * Do fit.
	 */
	public void doFit() {
		try {
			lma.fit();
			parameters=lma.parameters;
			gfit=lma.chi2Goodness();

		} catch (InvertException e) {
			for(int i=0;i<parameters.length;i++)parameters[i]=MRUtils.ERROR_VALUE;
		}
	}

	/**
	 * Do fit.
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
	 * Gets the params.
	 *
	 * @return the params
	 */
	public double[] getParams() {
		return parameters;
	}

	/**
	 * Gets the lma.
	 *
	 * @return the lma
	 */
	public LMA getLMA() {
		return lma;
	}

	/**
	 *  Get number of parameters for current fit formula.
	 *
	 * @param fitType the fit type
	 * @return the num params
	 */
	public static int getNumParams(int fitType) {
		switch (fitType) {
		case MRUtils.T1_MONO_RICE: return T1RecoveryRice.Nparams ;  
		case MRUtils.T2_MONO_RICE:  return  T2RelaxRice.Nparams  ;
		case MRUtils.T2_MULTI_RICE:  return  T2MulticompRice.Nparams  ;
		}
		return 0;
	}




	/**
	 * The Class T1RecoveryRice.
	 */
	/*    	p[0]*(1 - Math.exp(-(x / p[1]))); // p[1] - repetition times
	 */
	public static class T1RecoveryRice extends LMAFunction {
		
		/** The sigma. */
		public double sigma=0;
		
		/** The lm. */
		public LMCurveFitterNoBias lm;
		
		/**
		 * Sets the lm.
		 *
		 * @param lm the new lm
		 */
		public void setLM(LMCurveFitterNoBias lm) {this.lm=lm;}
		
		/**
		 * Sets the sigma.
		 *
		 * @param sig the new sigma
		 */
		public void setSigma(double sig){
			sigma=sig;
		}
		
		/** The Constant Nparams. */
		public static final int Nparams=2;
		 
		/**
		 * Gets the y.
		 *
		 * @param x the x
		 * @param a the a
		 * @return the y
		 */
		public double getY(double x, double[] a) {
			if(lm.debugLM)System.out.println("T1RecovRice sigma="+sigma+"  x="+x+"  a[0]="+a[0]+"  a[1]="+a[1]+"   res="+besFunkCost(a[0]*(1 - Math.exp(-(x / a[1]))),sigma));
			if(lm.debugLM)VitimageUtils.waitFor(10);
			return  besFunkCost(a[0]*(1 - Math.exp(-(x / a[1]))),sigma); // a[1] - repetition times
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
		public double getPartialDerivate(double x, double[] a, int parameterIndex) {
			switch (parameterIndex) {
			case 0: return  besFunkDeriv( a[0]*(1 - Math.exp(-(x / a[1]))), sigma )*(1.0-Math.exp(-x/a[1]));
			case 1: return  besFunkDeriv( a[0]*(1 - Math.exp(-(x / a[1]))), sigma )*(-(a[0]*x*Math.exp(-x/a[1]))/(a[1]*a[1]));

			} throw new RuntimeException("No such parameter index: " + parameterIndex);
		}

		/**
		 * To string.
		 *
		 * @return the string
		 */
		public String toString() {
			return "y=a[0]*(1 - Math.exp(-(x / a[1])))";
		}
	} // end class



	/**
	 * The Class T1MulticompRice.
	 */
	public static class T1MulticompRice extends LMAFunction {
		
		/** The sigma. */
		public double sigma=0;
		
		/** The lm. */
		public LMCurveFitterNoBias lm;
		
		/**
		 * Sets the lm.
		 *
		 * @param lm the new lm
		 */
		public void setLM(LMCurveFitterNoBias lm) {this.lm=lm;}
		
		/**
		 * Sets the sigma.
		 *
		 * @param sig the new sigma
		 */
		public void setSigma(double sig){
			sigma=sig;
		}
		
		/** The Constant Nparams. */
		public static final int Nparams=2;
		 
		/**
		 * Gets the y.
		 *
		 * @param x the x
		 * @param a the a
		 * @return the y
		 */
		public double getY(double x, double[] a) {
			return  besFunkCost(a[0]*(1 - Math.exp(-(x / a[1])))+a[2]*(1 - Math.exp(-(x / a[3]))),sigma); // a[1] - repetition times
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
		public double getPartialDerivate(double x, double[] a, int parameterIndex) {
			switch (parameterIndex) {
			case 0: return  besFunkDeriv( a[0]*(1 - Math.exp(-(x / a[1]))), sigma )*(1.0-Math.exp(-x/a[1]));
			case 1: return  besFunkDeriv( a[0]*(1 - Math.exp(-(x / a[1]))), sigma )*(-(a[0]*x*Math.exp(-x/a[1]))/(a[1]*a[1]));
			case 2: return  besFunkDeriv( a[2]*(1 - Math.exp(-(x / a[3]))), sigma )*(1.0-Math.exp(-x/a[3]));
			case 3: return  besFunkDeriv( a[2]*(1 - Math.exp(-(x / a[3]))), sigma )*(-(a[2]*x*Math.exp(-x/a[3]))/(a[3]*a[3]));

			} throw new RuntimeException("No such parameter index: " + parameterIndex);
		}

		/**
		 * To string.
		 *
		 * @return the string
		 */
		public String toString() {
			return "y=a[0]*(1 - Math.exp(-(x / a[1])))";
		}
	} // end class


	/*
	 * T2_RELAX_RICE:
	                return p[0]*Math.exp(-(x / p[1]) )+ p[2]; // p[1] - echo times
	 */

	/**
	 * The Class T2RelaxRice.
	 */
	public static class T2RelaxRice extends LMAFunction {
		
		/** The Constant Nparams. */
		public static final int Nparams=2;
		
		/** The sigma. */
		public double sigma=0;
		
		/** The lm. */
		public LMCurveFitterNoBias lm;
		
		/**
		 * Sets the lm.
		 *
		 * @param lm the new lm
		 */
		public void setLM(LMCurveFitterNoBias lm) {this.lm=lm;}
		
		/**
		 * Sets the sigma.
		 *
		 * @param sig the new sigma
		 */
		public void setSigma(double sig){
			sigma=sig;
		}

		/**
		 * Gets the y.
		 *
		 * @param x the x
		 * @param a the a
		 * @return the y
		 */
		@Override
		public double getY(double x, double[] a) {
			if(lm.debugLM)System.out.println("T2Rice sigma="+sigma+"  x="+x+"  a[0]="+a[0]+"  a[1]="+a[1]+"   res="+besFunkCost( a[0]* Math.exp(-(x / a[1])),sigma));
			if(lm.debugLM)VitimageUtils.waitFor(10);
			return  besFunkCost( a[0]* Math.exp(-(x / a[1])),sigma);	
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
		public double getPartialDerivate(double x, double[] a, int parameterIndex) {
			switch (parameterIndex) {
			case 0: return  besFunkDeriv( a[0]* Math.exp(-(x / a[1])), sigma ) * Math.exp(-x/a[1]);  
			case 1: return  besFunkDeriv( a[0]* Math.exp(-(x / a[1])), sigma ) * a[0]*x*Math.exp(-x/a[1])/(a[1]*a[1]);
			} 

			throw new RuntimeException("No such parameter index: " + parameterIndex);
		}

		/**
		 * To string.
		 *
		 * @return the string
		 */
		public String toString() {
			return "y=a[0]* Math.exp(-(x / a[1])) +a[2] + rice noise (sigma)";
		}
	} // end class




	/*
	 * T2MULTICOMP_RICE
	                return p[0]*Math.exp(-(x / p[1]) )+ p[2]; // p[1] - echo times
	 */

	/**
	 * The Class T2MulticompRice.
	 */
	public static class T2MulticompRice extends LMAFunction {
		
		/** The Constant Nparams. */
		public static final int Nparams=4;
		
		/** The sigma. */
		public double sigma=0;
		
		/** The lm. */
		public LMCurveFitterNoBias lm;
		
		/**
		 * Sets the lm.
		 *
		 * @param lm the new lm
		 */
		public void setLM(LMCurveFitterNoBias lm) {this.lm=lm;}
		
		/**
		 * Sets the sigma.
		 *
		 * @param sig the new sigma
		 */
		public void setSigma(double sig){
			sigma=sig;
		}

		
		/**
		 * Gets the y.
		 *
		 * @param x the x
		 * @param a the a
		 * @return the y
		 */
		@Override
		public double getY(double x, double[] a) {
			//return  a[0]* Math.exp(-(x / a[1])) +a[2]; // a[1] - echo times
			return  besFunkCost( a[0]* Math.exp(-(x / a[1])) + a[2]* Math.exp(-(x / a[3])), sigma); // a[1] - echo times
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
		public double getPartialDerivate(double x, double[] a, int parameterIndex) {
			switch (parameterIndex) {
			case 0: return  besFunkDeriv( a[0]* Math.exp(-(x / a[1])) + a[2]* Math.exp(-(x / a[3])), sigma)*Math.exp(-x/a[1]);  
			case 1: return  besFunkDeriv( a[0]* Math.exp(-(x / a[1])) + a[2]* Math.exp(-(x / a[3])), sigma)*a[0]*x*Math.exp(-x/a[1])/(a[1]*a[1]);
			case 2: return  besFunkDeriv( a[0]* Math.exp(-(x / a[1])) + a[2]* Math.exp(-(x / a[3])), sigma)*Math.exp(-x/a[3]);  
			case 3: return  besFunkDeriv( a[0]* Math.exp(-(x / a[1])) + a[2]* Math.exp(-(x / a[3])), sigma)*a[2]*x*Math.exp(-x/a[3])/(a[3]*a[3]);
			} 

			throw new RuntimeException("No such parameter index: " + parameterIndex);
		}

		/**
		 * To string.
		 *
		 * @return the string
		 */
		public String toString() {
			return "y=a[0]* Math.exp(-(x / a[1])) + a[2]* Math.exp(-(x / a[3])) + Rice noise estimation";
		}
	} // end class



}