package io.github.rocsg.fijirelax.mrialgo;

import java.util.Random;

import io.github.rocsg.fijiyama.common.VitimageUtils;

/** This class provides utilities to estimate MRI relaxation parameters in presence of a Rice noise
 * Rice noise affects the signal in a way that make it complicated to invert : its moments depend on the value of the unaltered signal
 * For more information, refer to Fernandez et al. 2022 FijiRelax: Fast and noise-corrected estimation of MRI relaxation maps in 3D + t (in prep.)
 * 
 * 
 */
public class RiceEstimator {


	/** The Constant epsilon. */
	private final static double epsilon=0.00001;
	
	/** The Constant debug. */
	private final static boolean debug=false;
	
	/** The sigma rice min. */
	private double sigmaRiceMin;
	
	/** The sigma rice max. */
	private double sigmaRiceMax;
	
	/** The observation min. */
	private double observationMin;
	
	/** The observation max. */
	private double observationMax;
	
	/** The observation range. */
	private double []observationRange;
	
	/** The sigma range. */
	private double []sigmaRange;
	
	/** The n sig. */
	private int nSig;
	
	/** The n obs. */
	private int nObs;
	
	/** The lut. */
	private double[][][]lut;
	
	

	
	/**
	 * Test of the random Rice simulator
	 */
	public static double[] testRandomRiceMeanDrift() {
		int N=100000;
		int SNR=2;
		double signal=700;
		double sigma=signal/SNR;
		double[]tab=new double[N];
		for(int i=0;i<N;i++) {
			tab[i]=getRandomRiceRealization(signal,sigma);
		}
		double mean=VitimageUtils.statistics1D(tab)[0];
		System.out.println(VitimageUtils.statistics1D(tab)[0]);
		return new double[] {mean,795.4,4};
	}
	
	/**
	 * Test of the random Rice simulator
	 */
	public static double[] testSimulateAndEstimateRiceSigma() {
		double sigma=0.05;
		double[]stats=computeSigmaAndMeanBgFromRiceSigmaStatic(sigma);
		System.out.println(stats[0]);
		return new double[] {stats[0],0.0626657068,VitimageUtils.EPSILON};
		//Sigma rice="+sigma+"  equivalent mean in zero signal="+stats[0]+"  equivalent std in zero signal="+stats[1]);
	}

	
	/**
	 * Test of the random Rice simulator and correction scheme
	 */
	public static double []testCorruptAndRecoverSignal() {
		//Build an estimator.
		RiceEstimator rice=getDefaultRiceEstimatorForNormalizedHyperEchoesT1AndT2Images();
		double sigma=0.05;
		for(double obs=0;obs<1;obs+=0.01) {
			double index=rice.getSigmaCoord(sigma);
			double coord=rice.getObservationCoord(0.5);
			System.out.println("Obs="+VitimageUtils.dou(obs)+"  Index="+index+" Coord="+coord);
		}		

		double[]sigs=rice.estimateSigmas(100,new double[] {120,1000});
		System.out.println("Sigmas="+sigs[0]+" , "+sigs[1]);
		
		
		double val=computeRiceSigmaFromBackgroundValuesStatic(2,1.2); 
		val=computeRiceSigmaFromBackgroundValuesStatic(20,1.2); 
		return new double[] {val, 15.957691,VitimageUtils.EPSILON};
		
	}
	
	
	
	/**
	 * Gets the default rice estimator for normalized hyper echoes T 1 and T 2 images, includin a common scale of possible signals and possible rice values
	 *
	 * @return the default rice estimator for normalized hyper echoes T 1 and T 2 images
	 */
	public static RiceEstimator getDefaultRiceEstimatorForNormalizedHyperEchoesT1AndT2Images() {
		RiceEstimator rice=new RiceEstimator(0.0001,1,0.0001,20,1000,1000);
		rice.start();
		return rice;
	}
	
	/**
	 * Instantiates a new rice estimator.
	 *
	 * @param sigmaRiceMin the minimal possible sigma for Rice 
	 * @param sigmaRiceMax the maximal possible sigma for Rice 
	 * @param observationMin the minimal possible observation
	 * @param observationMax the maximal possible observation
	 * @param nSig the number of sigma
	 * @param nObs the number of observations
	 */
	//TODO : document
	public RiceEstimator(double sigmaRiceMin,double sigmaRiceMax,double observationMin,double observationMax,int nSig,int nObs) {
		this.sigmaRiceMax=sigmaRiceMax;
		this.sigmaRiceMin=sigmaRiceMin;
		this.observationMin=observationMin;
		this.observationMax=observationMax;
		this.nSig=nSig;
		this.nObs=nObs;
		this.observationRange=new double[nObs];
		this.sigmaRange=new double[nSig];
	}

	/**
	 * Starting point after instantiating RiceEstimator
	 */
	public void start() {
		buildSigmaRange();
		buildObservationRange();
		buildLookupTable();
	}
	
	
	
	/**
	 * Builds the sigma range tab
	*/
	public void buildSigmaRange() {
		for(int sig=0;sig<nSig;sig++) {
			sigmaRange[sig]=sigmaRiceMin+(sig*1.0)*(sigmaRiceMax-sigmaRiceMin)/(nSig-1);
		}
	}
	
	/**
	 * Gets the coordinate to access the sigma value
	 *
	 * @param sigma the sigma
	 * @return the sigma coord
	 */
	//index of a sigma in the table
	public double getSigmaCoord(double sigma) {
		double coord=(sigma-sigmaRiceMin)/(sigmaRiceMax-sigmaRiceMin)*(nSig-1);
		if(coord<=0)coord=epsilon;
		if(coord>=nSig-1)coord=nSig-1-epsilon;
		return coord;
	}
	
	/**
	 * Builds the observation range.
	 */
	public void buildObservationRange() {
		if(debug)System.out.print("Construction array observations : ");
		double multFactor=this.observationMax/this.observationMin;
		for(int obs=0;obs<nObs;obs++) {		
			observationRange[obs]=this.observationMin*Math.pow(multFactor,obs*1.0/(nObs-1));
			if(debug)System.out.print("obs"+obs+"="+observationRange[obs]+" , ");
		}		
		if(debug)System.out.println();
	}

	/**
	 * Gets the coordinate corresponding to this observation
	 *
	 * @param observation the observation
	 * @return the observation coord
	 */
	//index of an observation in the table
	public double getObservationCoord(double observation) {
		if(observation>=this.observationMax)return nObs-1-epsilon;
		else if(observation<=this.observationMin)return epsilon;
		
		double valMin=this.observationMin;double valMax=this.observationMax;int indMin=0;int indMax=nObs-1;int indMed;double valMed;
		do {
			indMed=(int)Math.round(0.5*(indMax+indMin));
			valMed=observationRange[indMed];
			if(observation<valMed) {indMax=indMed;valMax=valMed;}
			else {indMin=indMed;valMin=valMed;}			
		}
		while((indMax-indMin)>1);
		return indMin+(observation-valMin)/(valMax-valMin);
	}

	
	/**
	 * Builds the lookup table among possible values of sigma and observations
	 */
	public void buildLookupTable() {
		double interMin, interMax,interMed,valMin,valMax,valMed;
		System.out.print("Building lookup table for estimation of initial signal amplitude using max likelyhood reverse estimation from observations X rice noise ...");
		int nTot=this.nObs*this.nSig;
		int incr=0;
		this.lut=new double[nSig][nObs][3];//SigmaRice, valObservation, EquivalentNormalSigma
		for(int sig=0;sig<nSig;sig++) {
			double curSigRice=sigmaRange[sig];
			for(int obs=0;obs<nObs;obs++) {
				if((incr++)%50000==0)System.out.print(""+(incr*100)/nTot+"%  ");
				double curObs=observationRange[obs];

				//Looking for the most likelyhood initial signal which leads to observe curObs, when noised by a Rice noise of parameter curSigRice 
				interMax=curObs;interMin=0;interMed=0;
				if(besFunkCost(interMin,curSigRice)>curObs) {}
				else {
					//Recherche iterative
					valMin=besFunkCost(interMin, curSigRice);
					valMax=besFunkCost(interMax, curSigRice);
					do {
						interMed=(interMax+interMin)/2;
						valMed=besFunkCost(interMed, curSigRice);
						if(valMed>=curObs) {interMax=interMed;valMax=valMed;}
						else {interMin=interMed;valMin=valMed;}
					}
					while(Math.abs(valMed-curObs)>epsilon);
				}
				lut[sig][obs]=new double[] {curSigRice,curObs,interMed};
				if(debug)System.out.println("sig"+sig+"  obs"+obs+"  : sigma="+curSigRice+"  observation="+ curObs+"  initialSignal="+interMed);
			}
		}
		System.out.println();
	}
			
	
	/**
	 * Estimate the besselFunction sigma for these values
	 *
	 * @param sigmaRice the sigma rice
	 * @param observation the observation
	 * @return the double
	 */
	public double estimateSigma(double sigmaRice,double observation) {
		double initialSignal=estimateInitialSignal(sigmaRice,observation);
		return besFunkSigma(initialSignal,sigmaRice);
	}
	
	/**
	 * Estimate the besselFunction sigma for these values
	 *
	 * @param sigmaRice the sigma rice
	 * @param observations the observations
	 * @return the double[]
	 */
	public double []estimateSigmas(double sigmaRice,double []observations) {
		double[]ret=new double[observations.length];
		for(int r=0;r<ret.length;r++)ret[r]=estimateSigma(sigmaRice,observations[r]);
		return ret;
	}
	
	/**
	 * Estimate the initial uncorrupted signal, given the value of sigmaRice, and the value of observation of the magnitude signal
	 *
	 * @param sigmaRice the sigma rice
	 * @param observation the observation
	 * @return the double
	 */
	public double estimateInitialSignal(double sigmaRice,double observation) {
		double coordObs=getObservationCoord(observation);
		double coordSig=getSigmaCoord(sigmaRice);
		int corObs0=(int)Math.floor(coordObs);
		int corSig0=(int)Math.floor(coordSig);
		double dObs=coordObs-corObs0;
		double dSig=coordSig-corSig0;
		double targetValue=  dSig    *    dObs    *this.lut[corSig0+1][corObs0+1][2]  +  
		              		 (1-dSig)*    dObs    *this.lut[corSig0  ][corObs0+1][2]  +  
		              		 dSig    *  (1-dObs)  *this.lut[corSig0+1][corObs0  ][2]  +  
		              		 (1-dSig)*  (1-dObs)  *this.lut[corSig0  ][corObs0  ][2];
	
		double erreur=Math.abs(besFunkCost(targetValue, sigmaRice)-observation)/observation*100.0;
		if(debug)System.out.println("\nESTIMATING INITIAL SIGMA from sigmaRice="+sigmaRice+" and obs="+observation);
		if(debug)System.out.println("Résultat : "+targetValue+ " jitter="+erreur+" % ");
		if(debug) {
			System.out.println("\nESTIMATING INITIAL SIGMA from sigmaRice="+sigmaRice+" and obs="+observation);
			System.out.println("Coordonnees detectees : coordSig="+coordSig+"  coordObs="+coordObs);		
			System.out.println("Valeurs 1 1 pond "+(dSig*dObs)+" : sig="+this.lut[corSig0 +1 ][corObs0  +1  ][0]+"  obs="+this.lut[corSig0  +1  ][corObs0  +1  ][1]+"  signal="+this.lut[corSig0  +1  ][corObs0   +1 ][2]);
			System.out.println("Valeurs 0 1 pond "+((1-dSig)*dObs)+" : sig="+this.lut[corSig0  ][corObs0 +1 ][0]+"  obs="+this.lut[corSig0  ][corObs0 +1 ][1]+"  signal="+this.lut[corSig0  ][corObs0 +1 ][2]);
			System.out.println("Valeurs 1 0 pond "+(dSig*(1-dObs))+" : sig="+this.lut[corSig0 +1 ][corObs0  ][0]+"  obs="+this.lut[corSig0 +1 ][corObs0  ][1]+"  signal="+this.lut[corSig0 +1 ][corObs0  ][2]);
			System.out.println("Valeurs 0 0 pond "+((1-dSig)*(1-dObs))+" : sig="+this.lut[corSig0  ][corObs0  ][0]+"  obs="+this.lut[corSig0  ][corObs0  ][1]+"  signal="+this.lut[corSig0  ][corObs0  ][2]);
			System.out.println("Résultat : "+targetValue);
		}
		return targetValue;
	}
	

	
	/**
	 * Bes funk cost. Helper function to ease computing the derivative of the signal as a function of MRI parameters
	 *
	 * @param d the d
	 * @param sigma2 the sigma 2
	 * @return the double
	 */
	static double besFunkCost(double d,double sigma2) {
		if(sigma2<=0)return d;
		double alpha=d*d/(4*sigma2*sigma2);
		return (double)(Math.sqrt(Math.PI*sigma2*sigma2/2.0)*( (1+2*alpha)*bessi0NoExp(alpha) + 2*alpha*bessi1NoExp(alpha) ));
	}
	

	/**
	 * Bes funk sigma.Helper function to ease computing the derivative of the signal as a function of MRI parameters
	 *
	 * @param d the d
	 * @param sigma2 the sigma 2
	 * @return the double
	 */
	/*
	 * Helper function to ease computing the derivative of the signal as a function of MRI parameters*/
	static double besFunkSigma(double d,double sigma2) {
		double alpha=d*d/(4*sigma2*sigma2);
		return Math.sqrt(2*sigma2*sigma2+d*d  - (Math.PI*sigma2*sigma2/2.0)*Math.pow( (1+2*alpha)*bessi0NoExp(alpha) + 2*alpha*bessi1NoExp(alpha) ,2));
	}
	


	
	/**
	 * Bessi 0 no exp. Function of bessi of 0th order, as given in numerical recipes
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
	 * Bessi 1 no exp. Function of bessi of 1th order, as given in numerical recipes
	 *
	 * @param alpha the alpha
	 * @return the double
	 */
	/*
	 * Function of bessi of 1th order, as given in numerical recipes
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
	 * Estimation of the rice noise stddev from values measured in the background
	 * assuming that there is no object in such background. The value is computed by two methods, using the mean and stdev of the BG, then the first one is returned
	 * If the values estimated by the two methods diverge by more than 30 %, display a Warning
	 *
	 * @param meanBg the mean background value
	 * @param sigmaBg the sigma of background values
	 * @return The estimated sigma of the Rice noise
	 */
	public static double computeRiceSigmaFromBackgroundValuesStatic(double meanBg,double sigmaBg) {
		boolean debug=true;
		double val1=meanBg * Math.sqrt(2.0/Math.PI);
		double val2=sigmaBg * Math.sqrt(2.0/(4.0-Math.PI));
		if(debug && Math.abs((val1-val2)/val2) >0.3) System.out.println("Warning : Acquisition > computeRiceSigmaStatic. Given :M="+meanBg+" , S="+sigmaBg+" gives sm="+val1+" , ss="+val2+"  .sm/ss="+VitimageUtils.dou(val1/val2)+". Using the first one..");
		
		return val1;
	}

	

	/**
	 * Estimation of the rice noise stddev from values measured in the background
	 * assuming that there is no object in such background. The value is computed by two methods, using the mean and stdev of the BG, then the two values are returned
	 *
	 * @param sigmaRice the sigma rice
	 * @return A double array yielding the estimated sigma in the two methods
	 */
	public static double []computeSigmaAndMeanBgFromRiceSigmaStatic(double sigmaRice) {
		boolean debug=true;
		double meanBg=sigmaRice/(Math.sqrt(2.0/Math.PI));
		double stdBg=sigmaRice/(Math.sqrt(2.0/(4.0-Math.PI)));
		return new double[] {meanBg,stdBg};
	}

	
	/**
	 * Simulation of Rice Noise over an unaltered signal.
	 *
	 * @param originalSignal the original signal
	 * @param sigmaRice the sigma rice
	 * @param nRepets the n repets
	 * @return A double corresponding to the value of the altered signal
	*/ 
	public static double getRandomRiceRealization(double originalSignal,double sigmaRice,int nRepets){
		double acc=0;
		for(int i=0;i<nRepets;i++) {
			acc+=getRandomRiceRealization(originalSignal,sigmaRice);
		}
		return acc/nRepets;
	}

	/**
	 * Simulation of Rice Noise over an unaltered signal.
	 *
	 * @param originalSignal the original signal
	 * @param sigmaRice the sigma rice
	 * @return A double corresponding to the value of the altered signal
	 */
	public static double getRandomRiceRealization(double originalSignal,double sigmaRice){
		Random rand=new Random();
		double gaussX=rand.nextGaussian()*sigmaRice;
		double gaussY=rand.nextGaussian()*sigmaRice;
		double sigmaT=rand.nextDouble()*2*Math.PI;
		double varX=originalSignal*Math.cos(sigmaT)+gaussX;
		double varY=originalSignal*Math.sin(sigmaT)+gaussY;
		return Math.sqrt(varX*varX+varY*varY);
	}
	
	
	
	
	
	
	
	
}