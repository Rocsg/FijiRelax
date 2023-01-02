/*
 * 
 */
package io.github.rocsg.fijirelax.mrialgo;


import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;
import java.util.concurrent.atomic.AtomicInteger;

import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import io.github.rocsg.fijiyama.registration.TransformUtils;
import io.github.rocsg.fijiyama.common.VitimageUtils;
import ij.IJ;
import ij.ImagePlus;
import ij.plugin.Duplicator;
import ij.process.ImageProcessor;
import io.github.rocsg.fijirelax.curvefit.LMDualCurveFitterNoBias;
import io.github.rocsg.fijirelax.curvefit.SimplexDualCurveFitterNoBias;

/**
 * This class is a utility class to hold many helpers used for simulating MRI echoes, computing rice noise estimation, reading information into TIFF metadata.
 * 	Copyright (C) 2022  io.github.rocsg
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
 *   along with this program.  If not, see https://www.gnu.org/licenses/
*/
public class MRUtils  {
	
	
	
	
	


	/** The Constant epsilon. */
	public static final double epsilon=10E-10;
	
	/** The Constant infinity. */
	public static final double infinity=1/epsilon;
	
	/** The Constant ERROR_VALUE. */
	public static final int ERROR_VALUE= 0;
	
	/** The Constant ERROR_KHI2. */
	public static final int ERROR_KHI2= 9999;
	
	/** The Constant SIMPLEX. */
	public static final int SIMPLEX = 1;
    
    /** The Constant LM. */
    public static final int LM=2;
    
    /** The Constant fititems2. */
    public final static String[] fititems2={"Simplex","Levenberg-Marquardt"};
    
    /** The Constant constitems2. */
    public final static int[] constitems2={SIMPLEX,LM};
    
    /** The Constant MSEC. */
    public static final int MSEC=11;
    
    /** The Constant SEC. */
    public static final int SEC=12;
	
	/** The Constant timeunits. */
	public final static String[] timeunits={"ms", "s"};
    
    /** The Constant timeitems. */
    public final static int[] timeitems={MSEC, SEC};

    /** The Constant TWOPOINTS. */
    public static final int TWOPOINTS=300;

    /** The Constant T1_MONO. */
    public static final int T1_MONO =1; //T1 exponential : PD  (1- exp(-tr/T1))
    
    /** The Constant T2_MONO. */
    public static final int T2_MONO =2; // T2 exponential : M0 exp(-te/T2)
    
    /** The Constant T2_MULTI. */
    public static final int T2_MULTI=3; // Two components exponential : PD1 exp(-te/T21) + PD2 exp(-te/T22)
	
	/** The Constant T1T2_MONO. */
	public static final int T1T2_MONO = 4; //T1-T2 cycle with a single T2 component : PD exp(-te/T2) (1- exp(-tr/T1))
	
	/** The Constant T1T2_MULTI. */
	public static final int T1T2_MULTI = 5; //T1-T2 cycle with two T2 component : ( PD1exp(-te/T21) + PD2exp(-te/T22) ) (1- exp(-tr/T1))

	/** The Constant T1_MONO_BIAS. */
	public static final int T1_MONO_BIAS =11;//The same, with an additive offset 
    
    /** The Constant T2_MONO_BIAS. */
    public static final int T2_MONO_BIAS =12; 
    
    /** The Constant T2_MULTI_BIAS. */
    public static final int T2_MULTI_BIAS=13;
	
	/** The Constant T1T2_MONO_BIAS. */
	public static final int T1T2_MONO_BIAS = 14; 
	
	/** The Constant T1T2_MULTI_BIAS. */
	public static final int T1T2_MULTI_BIAS = 15; 

	/** The Constant T1_MONO_RICE. */
	public static final int T1_MONO_RICE =21;//The same, with an estimated rice noise contribution
    
    /** The Constant T2_MONO_RICE. */
    public static final int T2_MONO_RICE =22;
    
    /** The Constant T2_MULTI_RICE. */
    public static final int T2_MULTI_RICE=23; 
	
	/** The Constant T1T2_MONO_RICE. */
	public static final int T1T2_MONO_RICE = 24; 
	
	/** The Constant T1T2_MULTI_RICE. */
	public static final int T1T2_MULTI_RICE = 25;
	
	
	/** The Constant T1T2_BIONANO. */
	public static final int T1T2_BIONANO = 30; //Calcul de T1 par la méthode des polynomes, calcul de T2 par la méthode des deux points
	
	/** The Constant T1T2_DEFAULT_T2_MONO_RICE. */
	public static final int T1T2_DEFAULT_T2_MONO_RICE=31;//Amendement de T1T2 avec un delta T2 de 25 ms
	
	/** The Constant T1T2_DEFAULT_T2_MULTI_RICE. */
	public static final int T1T2_DEFAULT_T2_MULTI_RICE=32;//Amendement de T1T2 avec un delta T2 de 25 ms
	
	/** The min delta chi 2. */
	public static double minDeltaChi2=1E-6;

	/** The Constant PARAM_SIGMA. */
	public static final int PARAM_SIGMA=31;
	
	/** The Constant PARAM_TR. */
	public static final int PARAM_TR=32;
	
	/** The Constant PARAM_TE. */
	public static final int PARAM_TE=33;
    
    /** The Constant DELTA_TE_BOUTURE_TRICK. */
    public static final double DELTA_TE_BOUTURE_TRICK=20;
    
 
    /** The Constant factorT1M0MaxRatio. */
    public static final double factorT1M0MaxRatio=3;
    
    /** The Constant factorT1MinRatio. */
    public static final double factorT1MinRatio=1;
    
    /** The Constant factorT1MaxRatio. */
    public static final double factorT1MaxRatio=3;

    /** The Constant factorT2M0MaxRatio. */
    public static final double factorT2M0MaxRatio=3;
    
    /** The Constant factorT2MinRatio. */
    public static final double factorT2MinRatio=1;
    
    /** The Constant factorT2MaxRatio. */
    public static final double factorT2MaxRatio=3;
	
	/** The Constant THRESHOLD_RATIO_BETWEEN_MAX_ECHO_AND_SIGMA_RICE. */
	public static final double THRESHOLD_RATIO_BETWEEN_MAX_ECHO_AND_SIGMA_RICE=5.25;//1.25 stands for mean Rayleigh and 4 for the three std around the mean;
	
	/** The Constant THRESHOLD_RATIO_BETWEEN_MEAN_ECHO_AND_SIGMA_RICE. */
	public static final double THRESHOLD_RATIO_BETWEEN_MEAN_ECHO_AND_SIGMA_RICE=2;
    
	
	/** The max khi 2. */
	protected static double maxKhi2=5;
	
	/** The max khi 2 MONO. */
	protected static double maxKhi2MONO=10E8;
	
	/** The max khi 2 after fit. */
	protected static float maxKhi2AfterFit=10;

    /** The Constant minAcceptableBionanoT1. */
    public static final double minAcceptableBionanoT1=500;
	
	/** The Constant minAcceptableBionanoT1Bouture. */
	public static final double minAcceptableBionanoT1Bouture=100;
	
	/** The Constant minAcceptableBionanoT2. */
	public static final double minAcceptableBionanoT2=10;
	
	/** The Constant N_ITER_T1T2. */
	public static final int N_ITER_T1T2=400;
	
	/** The Constant maxAcceptableBionanoM0. */
	public static final double maxAcceptableBionanoM0=2*MRUtils.maxM0BionanoForNormalization;
	
	/** The Constant maxAcceptableBionanoT2. */
	public static final double maxAcceptableBionanoT2=500;
	
	/** The Constant maxDisplayedBionanoT2. */
	public static final int maxDisplayedBionanoT2=150;
	
	/** The Constant maxDisplayedBionanoT1. */
	public static final int maxDisplayedBionanoT1=3200;
	
	/** The Constant maxDisplayedBionanoM0. */
	public static final int maxDisplayedBionanoM0=8666;
	
	/** The Constant maxDisplayedBionanoGE3D. */
	public static final int maxDisplayedBionanoGE3D=20000;
	
	/** The Constant maxM0BionanoForNormalization. */
	public static final int maxM0BionanoForNormalization=6500;
	
	/** The Constant maxDisplayedPDratio. */
	public static final double maxDisplayedPDratio=1.05;
	
	/** The Constant maxDisplayedT1ratio. */
	public static final double maxDisplayedT1ratio=1.5;
	
	/** The Constant maxDisplayedT2ratio. */
	public static final double maxDisplayedT2ratio=0.7;
	
	/** The Constant THRESHOLD_BOUTURE_FACTOR. */
	public static final double THRESHOLD_BOUTURE_FACTOR=0.5;
	
	/** The Constant maxAcceptableBionanoT1. */
	public static final int maxAcceptableBionanoT1=3500;
    
    /** The use boundaries. */
    public static boolean useBoundaries=false;
	
	/** The Constant DEFAULT_N_ITER_T1T2. */
	protected static final int DEFAULT_N_ITER_T1T2 = 500;
	
	/** The Constant minPossibleMaxT1. */
	public static final int minPossibleMaxT1=3500;
    
    
    
    
	
	
	
	
	/**
	 * Instantiates a new MR utils.
	 */
	public MRUtils() {	}

	
	/**
	 * Eliminate the first echo of the sequence
	 *
	 * @param tab the initial measures
	 * @return the double[][]
	 */
	public static double[][]forgetLastT2(double[][]tab){
		double[][]tab2=new double[tab.length][];
		for(int i=0;i<tab.length;i++) {
			tab2[i]=new double[tab[i].length-1];
			for(int j=1;j<tab[i].length;j++) tab2[i][j-1]=tab[i][j];			
		}
		return tab2;
	}
	
	/**
	 * Eliminate the first echo of the sequence
	 *
	 * @param tab the tab
	 * @return the double[][][]
	 */
	public static double[][][]forgetLastT2(double[][][]tab){
		//System.out.println("Forget last T2 building : "+tab.length +" x  "+tab[0].length+" x "+(tab[0][0].length-1));
		double[][][]tab2=new double[tab.length][][];
		for(int i=0;i<tab.length;i++) {
			tab2[i]=new double[tab[i].length][];
			for(int j=0;j<tab[i].length;j++) {	
				tab2[i][j]=new double[tab[i][j].length-1];
				for(int k=1;k<tab[i][j].length;k++) tab2[i][j][k-1]=tab[i][j][k];			
			}
		}
		return tab2;
	}
	
	
	/**
	 * Simulate echoes of a T2 relaxation sequence.
	 *
	 * @param maps an ImagePlus giving the dimensions of the image to simulate
	 * @param Tr the value of the recovery time
	 * @param Te the double[] giving the succesive values of echoe times
	 * @param nRepetitions the averaging (number of repeated measurements to average for each pixel)
	 * @param typeAcq the acquisition type (MRDataType.T2SEQ)
	 * @param sigma the sigma of Rice noise to simulate
	 * @param nameSimul the name of the simulation, to give a title to the ImagePlus
	 * @return the resulting ImagePlus, with echoes organized along the channel slicer
	 */
	public static ImagePlus simulateEchoesT2Relax(ImagePlus maps,double Tr,double[]Te,int nRepetitions,MRDataType typeAcq,double sigma,String nameSimul) {
		double[]tabTr=new double[Te.length];
		MRDataType[]tab=new MRDataType[Te.length];
		for(int i=0;i<Te.length;i++) {tabTr[i]=Tr;tab[i]=MRDataType.T2SEQ;}
		return simulateEchoesT1T2Relax(maps,tabTr,Te,nRepetitions,tab,sigma,nameSimul) ;
	}

	/**
	 * Simulate echoes of a T2 relaxation sequence.
	 *
	 * @param maps an ImagePlus giving the dimensions of the image to simulate
	 * @param Tr the double[] values of the recovery time
	 * @param Te the value of echo times
	 * @param nRepetitions the averaging (number of repeated measurements to average for each pixel)
	 * @param typeAcq the acquisition type (MRDataType.T2SEQ)
	 * @param sigma the sigma of Rice noise to simulate
	 * @param nameSimul the name of the simulation, to give a title to the ImagePlus
	 * @return the resulting ImagePlus, with echoes organized along the channel slicer
	 */
	public static ImagePlus simulateEchoesT1Relax(ImagePlus maps,double[]Tr,double Te,int nRepetitions,MRDataType typeAcq,double sigma,String nameSimul) {
		double[]tabTe=new double[Tr.length];
		MRDataType[]tab=new MRDataType[Tr.length];
		for(int i=0;i<Tr.length;i++) {tabTe[i]=Te;tab[i]=MRDataType.T1SEQ;}
		return simulateEchoesT1T2Relax(maps,Tr,tabTe,nRepetitions,tab,sigma,nameSimul) ;
	}
	
	/**
	 * Simulate echoes of a T2 relaxation sequence.
	 *
	 * @param hyperMap an ImagePlus giving the dimensions of the image to simulate
	 * @param Tr the double[] value of the recovery times
	 * @param Te the double[] value of the echoe times
	 * @param nRepetitions the averaging (number of repeated measurements to average for each pixel)
	 * @param typeAcq the acquisition type (MRDataType.T2SEQ)
	 * @param sigma the sigma of Rice noise to simulate
	 * @param nameSimul the name of the simulation, to give a title to the ImagePlus
	 * @return the resulting ImagePlus, with echoes organized along the channel slicer
	 */
	public static ImagePlus simulateEchoesT1T2Relax(ImagePlus hyperMap,double[]Tr,double[]Te,int nRepetitions,MRDataType []typeAcq,double sigma,String nameSimul) {
		double maxMag=0;
		double maxTr=0;
		double maxTe=0;
		int nMapsInput=3;
		int nMapsOutput=0;
		boolean hasT1=false;
		boolean hasT2=false;
		boolean hasT1T2=false;
		for(MRDataType m : typeAcq) {
			if(m==MRDataType.T1SEQ)hasT1=true;
			if(m==MRDataType.T2SEQ)hasT2=true;
		}
		hasT1T2=hasT1 && hasT2;
		if(!hasT1T2) nMapsOutput=3;
		else {
			nMapsOutput=4;
			nMapsInput=4;
		}
		int nEchoes=Tr.length;
		ImagePlus[]maps=VitimageUtils.stacksFromHyperstackFastBis(hyperMap);
		ImagePlus[]tab=new ImagePlus[nEchoes+nMapsOutput];

		for(int i=0;i<nMapsOutput-1;i++)tab[i]=maps[i];
		for(int i=nMapsOutput-1;i<tab.length;i++)tab[i]=VitimageUtils.nullImage(maps[0]);
		float[][]echoVals=new float[nEchoes][];
		float[][]mapVals=new float[nMapsInput][];
		int dimZ=maps[0].getNSlices();
		int dimX=maps[0].getWidth();
		int dimY=maps[0].getHeight();
		for(int z=0;z<dimZ;z++) {
			for(int i=0;i<maps.length;i++) {
				mapVals[i]=(float[]) maps[i].getStack().getPixels(z+1);
			}
			for(int i=0;i<tab.length-nMapsOutput;i++) {
				echoVals[i]=(float[]) tab[i+nMapsOutput].getStack().getPixels(z+1);
				for(int x=0;x<dimX;x++) {
					for(int y=0;y<dimY;y++) {
						int pix=y*dimX+x;
						double acc=0;
						for(int r=0;r<nRepetitions;r++) {
							acc+=RiceEstimator.getRandomRiceRealization(
									MRUtils.getFitFunctionValue(
											Tr[i], Te[i],new double[] {mapVals[0][pix],mapVals[1][pix],mapVals[2][pix]}, 0, MRUtils.T1T2_MONO_RICE
									),sigma
								);	
						}
						echoVals[i][pix]=(float)(acc/nRepetitions); 
						if(echoVals[i][pix]>maxMag)maxMag=echoVals[i][pix];
					}
				}
			}
		}

		for(int i=0;i<nEchoes;i++) {
			if(Tr[i]>maxTr && Tr[i]<10000)maxTr=Tr[i];
			if(Te[i]>maxTe)maxTe=Te[i];
		}

		tab[0]=maps[0];//PD
		tab[nMapsOutput-1]=VitimageUtils.nullImage(tab[0]);//Mask
		if(hasT1T2) {tab[1]=maps[1];tab[2]=maps[2];}
		else if(hasT1)tab[1]=maps[1];
		else if(hasT2)tab[1]=maps[2];

		ImagePlus hyper=VitimageUtils.hyperStackingChannels(tab);
		for(int c=nMapsOutput;c<nMapsOutput+nEchoes;c++) {
			int cc=c-nMapsOutput;
			for(int z=0;z<dimZ;z++) {
				hyper.getStack().setSliceLabel(""+typeAcq[0]+"_"+nameSimul+"_TR="+VitimageUtils.dou(Tr[cc])+"_TE="+VitimageUtils.dou(Te[cc])+"_SIGMARICE="+VitimageUtils.dou(sigma),
						VitimageUtils.getCorrespondingSliceInHyperImage(hyper, c, z, 0));
			}
			hyper.setC(c+1);
			hyper.setDisplayRange(0, maxMag);
			//IJ.run(hyper,"Fire","");
		}
		for(int z=0;z<dimZ;z++) {
			hyper.getStack().setSliceLabel("PDMAP_"+nameSimul,VitimageUtils.getCorrespondingSliceInHyperImage(hyper, 0, z, 0));
			if(hasT1T2) {
				hyper.getStack().setSliceLabel("T1MAP_"+nameSimul,VitimageUtils.getCorrespondingSliceInHyperImage(hyper, 1, z, 0));
				hyper.getStack().setSliceLabel("T2MAP_"+nameSimul,VitimageUtils.getCorrespondingSliceInHyperImage(hyper, 2, z, 0));
			}
			else if(hasT1)hyper.getStack().setSliceLabel("T1MAP_"+nameSimul,VitimageUtils.getCorrespondingSliceInHyperImage(hyper, 1, z, 0));
			else if(hasT2)hyper.getStack().setSliceLabel("T2MAP_"+nameSimul,VitimageUtils.getCorrespondingSliceInHyperImage(hyper, 1, z, 0));
			hyper.getStack().setSliceLabel("MASKMAP_"+nameSimul,VitimageUtils.getCorrespondingSliceInHyperImage(hyper, nMapsOutput-1, z, 0));
		}

		hyper.setC(1);
		hyper.setDisplayRange(0, maxMag*MRUtils.maxDisplayedPDratio);
		//IJ.run(hyper,"Fire","");
		if(hasT1T2) {
			hyper.setC(2);
			hyper.setDisplayRange(0, maxTe*MRUtils.maxDisplayedT2ratio);
			//IJ.run(hyper,"Fire","");		
			hyper.setC(3);
			hyper.setDisplayRange(0, maxTr*MRUtils.maxDisplayedT1ratio);
			//IJ.run(hyper,"Fire","");		
		}
		else if(hasT1) {
			hyper.setC(2);
			hyper.setDisplayRange(0, maxTr*MRUtils.maxDisplayedT1ratio);
			//IJ.run(hyper,"Fire","");		
		}
		else if(hasT2) {
			hyper.setC(2);
			hyper.setDisplayRange(0, maxTe*MRUtils.maxDisplayedT2ratio);
			//IJ.run(hyper,"Fire","");		
		}
		
		hyper.setTitle("Hypermap_simulated_"+nameSimul);
		return hyper;
	}	
	
	/**
	 * Compute T1 T2 maps slice per slice, in a multi-threaded way.
	 *
	 * @param imgsTemp the array of succesive images obtained during the acquisition sequence
	 * @param mask the mask identifying which voxels should be computed or not
	 * @param sigmaSmoothing the sigma smoothing if there is a need to "smooth" the data before maps computation
	 * @param fitType the fit type (exemple :MRUtils.T1_MONO_BIAS)
	 * @param algType the alg type (exemple: MRUtils.LM)
	 * @param debug a flag for displaying info about how the thread are going
	 * @param forgetFirstT2 boolean flag to define whether or not the first echo should be taken in consideration
	 * @param silentMode boolean flag to activate the performance mode and make the measurement of optimal performance
	 * @return the tab of computed maps
	 */
	public static ImagePlus[] computeT1T2MapMultiThreadSlices(final ImagePlus[]imgsTemp,ImagePlus mask,double sigmaSmoothing,int fitType,int algType, boolean debug,boolean forgetFirstT2,boolean silentMode) {
		int[]dims=VitimageUtils.getDimensions(imgsTemp[0]);		final int X=dims[0];		final int Y=dims[1];		final int Z=dims[2];  final int YZ=Y*Z;
		int nEc=imgsTemp.length;

		int nCores=VitimageUtils.getNbCores();
		int nVox=X*Y*Z;
		int[][]listForThreads=VitimageUtils.listForThreads(nVox,nCores);
		ImagePlus []imgs=new ImagePlus[nEc];
		for(int i=0;i<nEc;i++) {imgs[i]=imgsTemp[i].duplicate();imgs[i]=VitimageUtils.convertToFloat(imgs[i]);}

		ImagePlus[]maps=new ImagePlus[4];
		maps[0]=VitimageUtils.nullImage(imgsTemp[0]);
		maps[1]=VitimageUtils.nullImage(imgsTemp[0]);
		maps[2]=VitimageUtils.nullImage(imgsTemp[0]);
		maps[3]=VitimageUtils.nullImage(imgsTemp[0]);
		final double [][]listTrForThreads=new double[nCores][];  
		final double [][]listTeForThreads=new double[nCores][];
		final double [][]listSigmaForThreads=new double[nCores][];
		final int [][]listIndex=new int[nCores][];
		final double [][][]listMask=new double[nCores][][];
		final double [][][]listData=new double[nCores][][];
		final double [][][]results=new double[nCores][][];
		final int []fitTypeTab;
		final int[]algTypeTab;
		int []tab1=new int[nCores];
		int[]tab2=new int[nCores];
		boolean[]tab3=new boolean[nCores];
		for(int nC=0;nC<nCores;nC++) {
			listTrForThreads[nC]=getTrFrom3DRelaxationImageTab(imgsTemp);
			listTeForThreads[nC]=getTeFrom3DRelaxationImageTab(imgsTemp)[0];
			listSigmaForThreads[nC]=getSigmaFrom3DRelaxationImageTab(imgsTemp);
			listIndex[nC]=copyTab(listForThreads[nC]);
			listData[nC]=new double[listForThreads[nC].length][];
			listMask[nC]=new double[listForThreads[nC].length][];
			tab1[nC]=fitType;
			tab2[nC]=algType;
			tab3[nC]=forgetFirstT2;
		}
		fitTypeTab=copyTab(tab1);
		algTypeTab=copyTab(tab2);
		final boolean []forgetFirstT2Tab=copyTab(tab3);
		final double [][]listTrForThreadsForget=forgetFirstT2 ? forgetLastT2( listTrForThreads ):null; 
		final double [][]listTeForThreadsForget=forgetFirstT2 ? forgetLastT2( listTeForThreads ):null; 
		
		//Construire le fetch list
		int [][]fetch=new int[X*Y*Z][5];
		for(int i=0;i<nCores;i++) {
			results[i]=new double[listIndex[i].length][];
			for(int j=0;j<listIndex[i].length;j++) {
				int index=listIndex[i][j];
				int x=index/(YZ);
				int y=(index-x*YZ)/Z;
				int z=index-x*YZ-y*Z;
				fetch[index]=new int[] {i,j,x,y,z};
			}
		}

		float[][][]imgData=new float[nEc][Z][];
		float[][][]mapData=new float[4][Z][];
		float[][]maskData=new float[Z][];
		//Ouvrir les voxels des images mask et all
		for(int z=0;z<Z;z++) {
			for(int n=0;n<nEc;n++) {
				imgData[n][z]=(float[])imgs[n].getStack().getPixels(z+1);
			}
			for(int n=0;n<4;n++) {
				mapData[n][z]=(float[])maps[n].getStack().getPixels(z+1);
			}
			maskData[z]=(float[])mask.getStack().getPixels(z+1);
		}

		//Collect data for each thread
		double[]data=new double[nEc];
		for(int ind=0;ind<fetch.length;ind++) {
			int proc=fetch[ind][0];
			int indexProc=fetch[ind][1];
			int x=fetch[ind][2];
			int y=fetch[ind][3];
			int z=fetch[ind][4];
			for(int n=0;n<nEc;n++) {
				data[n]=imgData[n][z][y*X+x];
			}			
			listData[proc][indexProc]=copyTab(data);
			listMask[proc][indexProc]=copyTab(new double[] { maskData[z][y*X+x] } );
		}
		final double [][][]listDataForget=forgetLastT2( listData ); 

		
		final AtomicInteger incrTotalFit = new AtomicInteger(0);
		final AtomicInteger incrProc = new AtomicInteger(0);
		final AtomicInteger nEarlyBreaks = new AtomicInteger(0);
		final int onePercent=1+(X*YZ)/10;
		final int XYZ=X*YZ;
		IJ.log((  ("Multi-threaded with ("+Z+" threads)"))+" T1 and/or T2 map computation.\n Start fit on "+(X*Y*Z)+" voxels with sigma="+listSigmaForThreads[0][0]);
		final Thread[] threads = VitimageUtils.newThreadArray(nCores);    
		for (int ithread = 0; ithread < nCores; ithread++) {  
			threads[ithread] = new Thread() {  { setPriority(Thread.NORM_PRIORITY); }  
			public void run() {  
				int numProc = incrProc.getAndIncrement();
				int nWork=listIndex[numProc].length;

				for(int w=0;w<nWork;w++) {
					
					int incrTot=incrTotalFit.getAndIncrement();
					if(incrTot%onePercent==0 && (!silentMode)) {
						IJ.log("Maps computation : "+(VitimageUtils.dou(incrTot*100.0/XYZ)+" %"));
						IJ.showProgress(incrTot*1.0/XYZ);
					}
					
					if(listMask[numProc][w][0]<=0) {
						results[numProc][w]=new double[] {0,0,0,0,0,0,0};
						continue;
					}

					Object []obj=null;
					if(forgetFirstT2Tab[numProc]) {
						obj=MRUtils.makeFit(listTrForThreadsForget[numProc], listTeForThreadsForget[numProc],  listDataForget[numProc][w],   fitTypeTab[numProc],  algTypeTab[numProc],  400 , listSigmaForThreads[numProc][0],false);						
					}
					else {
						obj=MRUtils.makeFit(listTrForThreads[numProc], listTeForThreads[numProc],  listData[numProc][w],   fitTypeTab[numProc],  algTypeTab[numProc],  400 , listSigmaForThreads[numProc][0],false);
					}
					double []estimatedParams=((double[])(obj[0]));
					if(( (boolean) obj[1]) )nEarlyBreaks.getAndIncrement();
					results[numProc][w]=estimatedParams;
				}				
			} } ;  
		}  		
		VitimageUtils.startAndJoin(threads);  

		System.out.println("End of fit. Early breaks = "+nEarlyBreaks.get()+" / "+(Z*X*Y));
		for(int ind=0;ind<fetch.length;ind++) {
			int proc=fetch[ind][0];
			int indexProc=fetch[ind][1];
			int x=fetch[ind][2];
			int y=fetch[ind][3];
			int z=fetch[ind][4];
			mapData[0][z][y*X+x]=(float) results[proc][indexProc][0];
			if(fitType==MRUtils.T1_MONO)mapData[1][z][y*X+x]=(float) results[proc][indexProc][1];
			else if(fitType==MRUtils.T1_MONO_BIAS)mapData[1][z][y*X+x]=(float) results[proc][indexProc][1];
			else if(fitType==MRUtils.T1_MONO_RICE)mapData[1][z][y*X+x]=(float) results[proc][indexProc][1];
			else if(fitType==MRUtils.T2_MONO)mapData[1][z][y*X+x]=(float) results[proc][indexProc][1];
			else if(fitType==MRUtils.T2_MONO_BIAS)mapData[1][z][y*X+x]=(float) results[proc][indexProc][1];
			else if(fitType==MRUtils.T2_MONO_RICE)mapData[1][z][y*X+x]=(float) results[proc][indexProc][1];
			else if(fitType==MRUtils.T1T2_MONO) {
				mapData[1][z][y*X+x]=(float) results[proc][indexProc][1];
				mapData[2][z][y*X+x]=(float) results[proc][indexProc][2];
			}
			else if(fitType==MRUtils.T1T2_MONO_BIAS) {
				mapData[1][z][y*X+x]=(float) results[proc][indexProc][1];
				mapData[2][z][y*X+x]=(float) results[proc][indexProc][2];
			}
			else if(fitType==MRUtils.T1T2_MONO_RICE) {
				mapData[1][z][y*X+x]=(float) results[proc][indexProc][1];
				mapData[2][z][y*X+x]=(float) results[proc][indexProc][2];
			}
			else {
				IJ.showMessage("Critical fail in MRUtils : fit not expected here : "+fitType);
			}
//			if(indexProc==246)System.out.println("Going to make the stuff "+proc+"-"+indexProc+" . And results is "+TransformUtils.stringVectorN(results[proc][indexProc], "")+" and setting in mapData"+MRUtils.getNparams(fitType%10)+" have leng"+mapData.length);
			mapData[MRUtils.getNparams(fitType%10)][z][y*X+x]=(float) results[proc][indexProc][MRUtils.getNparams(fitType)];
		}
		
		maps[0].setDisplayRange(0, maxDisplayedBionanoM0);
		maps[1].setDisplayRange(0, maxDisplayedBionanoT1);
		maps[2].setDisplayRange(0, maxDisplayedBionanoT2);
		maps[3].setDisplayRange(0, 2);

		ImagePlus []tabRet=null;
		switch(fitType) {
			case T1_MONO : tabRet=new ImagePlus[] {maps[0],maps[1],maps[2]};break;
			case T1_MONO_BIAS : tabRet=new ImagePlus[] {maps[0],maps[1],maps[2]};break;
			case T1_MONO_RICE : tabRet=new ImagePlus[] {maps[0],maps[1],maps[2]};break;

			case T2_MONO : tabRet=new ImagePlus[] {maps[0],maps[1],maps[2]};break;
			case T2_MONO_BIAS : tabRet=new ImagePlus[] {maps[0],maps[1],maps[2]};break;
			case T2_MONO_RICE : tabRet=new ImagePlus[] {maps[0],maps[1],maps[2]};break;

			case T1T2_MONO : tabRet=new ImagePlus[] {maps[0],maps[1],maps[2],maps[3]};break;
			case T1T2_MONO_BIAS : tabRet=new ImagePlus[] {maps[0],maps[1],maps[2],maps[3]};break;
			case T1T2_MONO_RICE : tabRet=new ImagePlus[] {maps[0],maps[1],maps[2],maps[3]};break;
			default : break;
		}
		return tabRet;
	}
	
	
	



	/**
	 *  Sort an array of double by increasing order
	 *
	 * @param tab the tab
	 * @return the double[]
	 */
	public static double[] doubleArraySort(double[]tab) {
		double[]tabRet=new double[tab.length];
		ArrayList<Double> list=new ArrayList<Double>();
		for(double dou : tab)list.add(dou);
		Collections.sort(list);
		for(int i=0;i<list.size();i++)tabRet[i]=list.get(i);
		return tabRet;
	}

	/**
	 *  Check data type in slice labels, by testing if it contains patterns (ex : "T1SEQ").
	 *
	 * @param img the image to test
	 * @param c the channel
	 * @param z the zslice
	 * @param f the timeframe
	 * @return the data type of this magnetic resonance slice
	 */
	public static MRDataType getDataTypeOfThisMagneticResonanceSlice(ImagePlus img,int c,int z,int f) {
		String label=img.getStack().getSliceLabel(VitimageUtils.getCorrespondingSliceInHyperImage(img,c,z,f) );		
		int nbCat=0;
		MRDataType dataT=null;
		if(label==null)return MRDataType.OTHER;
		if(label.contains("T1SEQ")) {dataT=MRDataType.T1SEQ;nbCat++;}
		if(label.contains("T1T2SEQ")) {dataT=MRDataType.T1T2SEQ;nbCat++;}
		if(label.contains("T2SEQ") && (!label.contains("T1T2SEQ")) ) {dataT=MRDataType.T2SEQ;nbCat++;}
		if(label.contains("T1MAP")) {dataT=MRDataType.T1MAP;nbCat++;}
		if(label.contains("T2MAP")) {dataT=MRDataType.T2MAP;nbCat++;}
		if(label.contains("PDMAP")) {dataT=MRDataType.PDMAP;nbCat++;}
		if(label.contains("M0MAP")) {dataT=MRDataType.PDMAP;nbCat++;}
		if(nbCat>1) {
			IJ.showMessage("Critical fail in MRUtils : get "+nbCat+" categories instead of 1 \nwhen calling getDataTypeOfThisMagneticResonanceSlice(ImagePlus img,"+c+","+z+","+f+"\nLabel was : "+label);
		}
		if(nbCat<1)return MRDataType.OTHER;
		return dataT;
	}
	
	/**
	 * Read values of the capillary in slices label
	 *
	 * @param img the image to test
	 * @param c the channel
	 * @param z the zslice
	 * @param t the timeframe
	 * @return the double[]
	 */
	public static double[] readCapValuesInSliceLabel(ImagePlus img,int c, int z, int t) {
		double []ret=new double[2];
		String label=img.getStack().getSliceLabel(VitimageUtils.getCorrespondingSliceInHyperImage(img, c, z, t));
		String[]strTab=label.split("_");
		for(String str:strTab) {
			if(str.contains("CAP")) {
				ret[0]=Double.parseDouble((str.split("=")[1]).split("\\+-")[0]);
				ret[1]=Double.parseDouble((str.split("=")[1]).split("\\+-")[1]);
			}
		}
		return ret;
	}

	/**
	 * Compute the mean of the double tab
	 *
	 * @param tab the tab
	 * @return the double
	 */
	protected static double mean(double[] tab) {
		double m=0;
		for(int i=0;i<tab.length;i++) m+=tab[i];
		return (m/tab.length);
	}

	

	
	
	
	/**
	 *  Estimate rice noise in each slice from the background values.
	 *
	 * @param img the img
	 * @param keepOldString flag indicating whether or not the slice label should be kept or updated
	 * @param optionalAdditionalPrefix the optional additional prefix to add before the measured value in slice label
	 * @param makeBoutureTrick flag for specific kind of data
	 */
	public static void computeTeTrAndRiceSigmaOfEachSliceAndWriteItInTheLabels(ImagePlus img,boolean keepOldString,String optionalAdditionalPrefix,boolean makeBoutureTrick) {
		int nbC=img.getNChannels();
		int nbZ=img.getNSlices();
		int nbF=img.getNFrames();
		ImageProcessor imgP;
		double tr=VitimageUtils.getRepetitionTime(img);
		double te=VitimageUtils.getEchoTime(img);
//		System.out.println("Make bouture trick ?"+makeBoutureTrick);
		for(int nc=0;nc<nbC;nc++) {
			for(int nz=0;nz<nbZ;nz++) {
				for(int nf=0;nf<nbF;nf++) {
					int sli=VitimageUtils.getCorrespondingSliceInHyperImage(img,nc, nz, nf);
					imgP=img.getStack().getProcessor(sli);
					double[]stats=getBackgroundStatsFromProcessor(imgP);
					double sigmaRice=RiceEstimator.computeRiceSigmaFromBackgroundValuesStatic(stats[0],stats[1]);
					String chain="_TR="+VitimageUtils.dou(tr)+"_TE="+VitimageUtils.dou(te+ (((makeBoutureTrick && tr>9000) ? (nz==0 ? 0 : 0) : 0 )));
					img.getStack().setSliceLabel((keepOldString ? img.getStack().getSliceLabel(sli) : "") + "_"+optionalAdditionalPrefix+chain+"_SIGMARICE="+VitimageUtils.dou(sigmaRice), sli);
				}
			}
		}
	}
	
	/**
	 * 	This routine compute the background stats in eight little square near the border
	 * Before estimation, the 3 ones with the highest values are discarded
	 * because there is often the object somewhere on the border, or the capillary.
	 *
	 * @param imgP the image processor where to estimate
	 * @return the background stats from processor
	 */
	public static double[]getBackgroundStatsFromProcessor(ImageProcessor imgP) {
		int dimX=imgP.getWidth();
		int dimY=imgP.getHeight();
		int samplSize=Math.min(10+20,dimX/10);
		if(dimX<100)samplSize=12;
		if(dimX>500)samplSize=40;

		int x0=(3*samplSize)/2;
		int y0=(3*samplSize)/2;
		int x1=dimX/2;
		int y1=dimY/2;
		int x2=dimX-(3*samplSize)/2;
		int y2=dimY-(3*samplSize)/2;
		double[][] vals=new double[8][];
		vals[0]=VitimageUtils.valuesOfImageProcessor(imgP,x0,y0,samplSize/2);
		vals[1]=VitimageUtils.valuesOfImageProcessor(imgP,x0,y2,samplSize/2);
		vals[2]=VitimageUtils.valuesOfImageProcessor(imgP,x2,y0,samplSize/2);
		vals[3]=VitimageUtils.valuesOfImageProcessor(imgP,x2,y2,samplSize/2);		
		vals[4]=VitimageUtils.valuesOfImageProcessor(imgP,x0,y1,samplSize/4);
		vals[5]=VitimageUtils.valuesOfImageProcessor(imgP,x1,y0,samplSize/4);
		vals[6]=VitimageUtils.valuesOfImageProcessor(imgP,x2,y1,samplSize/4);
		vals[7]=VitimageUtils.valuesOfImageProcessor(imgP,x2,y2,samplSize/4);

		//Compute the global mean over all these squares
		double[]tempStatsMeanVar=null;
		double[]tempStats=new double[vals.length];
		
		//Measure local stats, and guess that three of them can host the object
		for(int i=0;i<vals.length;i++) {
			tempStatsMeanVar=VitimageUtils.statistics1D(vals[i]);
			tempStats[i]=tempStatsMeanVar[0];
		}

		double[]tempStatsCorrected=null;
		int incr=0;
		double[][]valsBis=null;
		tempStatsCorrected=new double[5];//Suppress the 3 maximum, that should be the border or corner where the object lies

		double[]tempStats2=doubleArraySort(tempStats);
		for(int i=0;i<5;i++)tempStatsCorrected[i]=tempStats2[i];
		valsBis=new double[5][];			
		for(int i=0;i<8 && incr<5;i++) {if(tempStats[i]<=tempStatsCorrected[4]) {valsBis[incr++]=vals[i];}}
			
		double []valsRetBis=VitimageUtils.statistics2D(valsBis);		
		return valsRetBis;
	}	
	
	/**
	 * 	This routine compute the background stats in eight very little square near the border
	 * Before estimation, the 3 ones with the highest values are discarded
	 * because there is often the object somewhere on the border, or the capillary.
	 *
	 * @param imgP the image processor where to estimate
	 * @return the background stats from processor
	 */
	public static double[]getBackgroundStatsFromProcessorTight(ImageProcessor imgP) {
		int dimX=imgP.getWidth();
		int dimY=imgP.getHeight();
		int samplSize=Math.min(18,dimX/10);
		if(dimX<100)samplSize=12;
		if(dimX>500)samplSize=40;

		int x0=(samplSize)/2;
		int y0=(samplSize)/2;
		int x1=dimX/2;
		int y1=dimY/2;
		int x2=dimX-(samplSize)/2;
		int y2=dimY-(samplSize)/2;
		double[][] vals=new double[8][];
		vals[0]=VitimageUtils.valuesOfImageProcessor(imgP,x0,y0,samplSize/2);
		vals[1]=VitimageUtils.valuesOfImageProcessor(imgP,x0,y2,samplSize/2);
		vals[2]=VitimageUtils.valuesOfImageProcessor(imgP,x2,y0,samplSize/2);
		vals[3]=VitimageUtils.valuesOfImageProcessor(imgP,x2,y2,samplSize/2);		
		vals[4]=VitimageUtils.valuesOfImageProcessor(imgP,x0,y1,samplSize/4);
		vals[5]=VitimageUtils.valuesOfImageProcessor(imgP,x1,y0,samplSize/4);
		vals[6]=VitimageUtils.valuesOfImageProcessor(imgP,x2,y1,samplSize/4);
		vals[7]=VitimageUtils.valuesOfImageProcessor(imgP,x2,y2,samplSize/4);

		//Compute the global mean over all these squares
		double[]tempStatsMeanVar=null;
		double[]tempStats=new double[vals.length];
		
		//Measure local stats, and guess that three of them can host the object
		for(int i=0;i<vals.length;i++) {
			tempStatsMeanVar=VitimageUtils.statistics1D(vals[i]);
			tempStats[i]=tempStatsMeanVar[0];
		}

		double[]tempStatsCorrected=null;
		int incr=0;
		double[][]valsBis=null;
		tempStatsCorrected=new double[5];//Suppress the 3 maximum, that should be the border or corner where the object lies

		double[]tempStats2=doubleArraySort(tempStats);
		for(int i=0;i<5;i++)tempStatsCorrected[i]=tempStats2[i];
		valsBis=new double[5][];			
		for(int i=0;i<8 && incr<5;i++) {if(tempStats[i]<=tempStatsCorrected[4]) {valsBis[incr++]=vals[i];}}
			
		double []valsRetBis=VitimageUtils.statistics2D(valsBis);		
		return valsRetBis;
	}	
	
	
	

	
	
	/**
	 *  Extract Tr values from the relaxation images by reading this in slice labels
	 *
	 * @param imgTr the images where to look for the values
	 * @return the recovery time.
	 */
	public static double[]getTrFrom3DRelaxationImageTab(ImagePlus[]imgTr){
		double[]tab=new double[imgTr.length];
		for(int i=0;i<imgTr.length;i++) {
			tab[i]=readValueOfSigmaTrTeInSliceLabel(imgTr[i],PARAM_TR,0,0,0);
		}
		return tab;
	}

	/**
	 *  Extract Te values from the relaxation images 
	 *
	 * @param imgTr the img tr
	 * @return the echoes times
	 */
	public static double[][]getTeFrom3DRelaxationImageTab(ImagePlus[]imgTr){
		double[][]tab=new double[imgTr[0].getNSlices()][imgTr.length];
		for(int i=0;i<imgTr.length;i++) {
			for(int z=0;z<imgTr[0].getNSlices();z++) {
				tab[z][i]=readValueOfSigmaTrTeInSliceLabel(imgTr[i],PARAM_TE,0,z,0);
			}
		}
		return tab;
	}
	
	/**
	 * Read recovery time in slice label
	 *
	 * @param img the img
	 * @param c the channel
	 * @param z the zslice
	 * @param f the timeframe
	 * @return the Tr, as a double
	 */
	public static double readTrInSliceLabel(ImagePlus img,int c, int z, int f) {
		return readValueOfSigmaTrTeInSliceLabel(img, PARAM_TR, c, z, f);
	}

	/**
	 * Read Rice noise sigma in slice label
	 *
	 * @param img the img
	 * @param c the channel
	 * @param z the zslice
	 * @param f the timeframe
	 * @return the Sigma, as a double
	 */
	public static double readSigmaInSliceLabel(ImagePlus img,int c, int z, int f) {
		return readValueOfSigmaTrTeInSliceLabel(img, PARAM_SIGMA, c, z, f);
	}

	/**
	 * Read the echoe time in slice label
	 *
	 * @param img the img
	 * @param c the channel
	 * @param z the zslice
	 * @param f the timeframe
	 * @return the Te as a double
	 */
	public static double readTeInSliceLabel(ImagePlus img,int c, int z, int f) {
		return readValueOfSigmaTrTeInSliceLabel(img, PARAM_TE, c, z, f);
	}
	
	/**
	 * Read value of sigma tr te in slice label.
	 *
	 * @param img the img
	 * @param prm0Sigma_1Tr_2Te the prm 0 sigma 1 tr 2 te
	 * @param c the channel
	 * @param z the zslice
	 * @param f the timeframe
	 * @return the value as a double
	 */
	public static double readValueOfSigmaTrTeInSliceLabel(ImagePlus img,int prm0Sigma_1Tr_2Te,int c,int z,int f) {
		String tit=img.getStack().getSliceLabel(VitimageUtils.getCorrespondingSliceInHyperImage(img,c,z,f) );
		String[]datas=tit.replace("__","_").split("_");
		for(int i=0;i<datas.length;i++) {
			String []prm=datas[i].split("=");
			if(prm.length>1) {
				if(prm0Sigma_1Tr_2Te==PARAM_SIGMA && prm[0].equals("SIGMARICE"))return Double.parseDouble(prm[1]);
				if(prm0Sigma_1Tr_2Te==PARAM_TR && prm[0].equals("TR"))return Double.parseDouble(prm[1]);
				if(prm0Sigma_1Tr_2Te==PARAM_TE && prm[0].equals("TE"))return Double.parseDouble(prm[1]);
			}
		}
		return 0;
	}

	/**
	 * Modify Rice noise sigma, and write it in the slice label
	 *
	 * @param hyperImg the hyper img
	 * @param c the channel
	 * @param z the zslice
	 * @param t the timeframe
	 * @param newSigma the new sigma to write
	 */
	public static void modifySigma(ImagePlus hyperImg,int c,int z,int t,double newSigma) {
		String label=hyperImg.getStack().getSliceLabel(VitimageUtils.getCorrespondingSliceInHyperImage(hyperImg, c, z, t));
		String[]strTab=label.split("_");
		String newLab="";
		for(int i=0;i<strTab.length;i++) {
			String str=strTab[i];
			if(str.contains("SIGMARICE"))newLab+=("SIGMARICE="+VitimageUtils.dou(newSigma));
			else newLab+=str;
			if(i<strTab.length-1)newLab+="_";
		}
		hyperImg.getStack().setSliceLabel(newLab, VitimageUtils.getCorrespondingSliceInHyperImage(hyperImg, c, z, t));
	}


	
	
	
	
	
	/**
	 *  Helper to compute sigma Rice mean over the successive echoes image
	 * If one of the sigma is +-30% from the mean, a warning is sent to the log console.
	 *
	 * @param imgTr the img tab for estimating the sigmas
	 * @return the sigma Rice tab estimation in each slice
	 */
	public static double[]getSigmaFrom3DRelaxationImageTab(ImagePlus[]imgTr){
		double[]tab=new double[imgTr[0].getNSlices()];
		double[]tabTemp=new double[imgTr.length];
		for(int i=0;i<tab.length;i++) {
			for(int j=0;j<imgTr.length;j++) {
				tabTemp[j]=readValueOfSigmaTrTeInSliceLabel(imgTr[j],PARAM_SIGMA,0,i,0);
			}
			tab[i]=VitimageUtils.statistics1D(tabTemp)[0];
			for(int j=0;j<imgTr.length;j++) {
				if(tabTemp[j]>tab[i]*1.3 || tabTemp[j]<tab[i]*0.7)IJ.log("Warning in MRUtils getSigmaFrom3DRelaxation. Values="+TransformUtils.stringVectorN(tabTemp,"")); 
			}
		}			
		return tab;
	}


	
	/**
	 *  Helper to make a linspace along relaxation times in order to prepare display of the relaxation curves.
	 *
	 * @param valMin the val min
	 * @param valMax the val max
	 * @param step the step
	 * @return the proportional times
	 */
	public static double []getProportionalTimes(double valMin, double valMax,double step){
		double[]tab=new double[(int)(Math.ceil((double)0.00001+((valMax-valMin)/step)))];
		for (int i=0;i<tab.length;i++){
			tab[i]=valMin+(double)((i*step));
		}
		return tab;
	}
	
	
	
	
	
	///////////////// Helpers for accessing signal in magnitude image with successive echoes, in a given roi ////////////////////
	/**
	 *  Give the average decay signal of the voxels in the neighbouring of v(xCor,yCor,curZ,curTime
	 * The neighbouring is a parallepiped of size (2*crossWidth+1) along X-Y and (2*crossThick+1) along Z.
	 *
	 * @param imgIn the img in
	 * @param xCor the x0
	 * @param yCor the y0
	 * @param curT the current time
	 * @param totalT the total T
	 * @param curZ the cur Z
	 * @param totalZ the total Z
	 * @param nEchoesExpected the n echoes expected
	 * @param crossWidth the cross width
	 * @param crossThick the cross thick
	 * @param gaussianWeighting the gaussian weighting
	 * @return the data for voxel
	 */
	public static double[]getDataForVoxel(ImagePlus imgIn,int xCor,int yCor,int curT,int totalT,int curZ,int totalZ,int nEchoesExpected,int crossWidth,int crossThick,boolean gaussianWeighting){
		double[][]tabData=getFullDataForVoxel(imgIn, xCor, yCor, curT, totalT, curZ, totalZ, nEchoesExpected, crossWidth, crossThick);
		double[]tabRet=new double[nEchoesExpected];
		for(int ec=0;ec<nEchoesExpected;ec++) 			tabRet[ec]=VitimageUtils.mean(tabData[ec]);
		return tabRet;
	}

	/**
	 *  Give the magnitude decay signal of the voxels in the neighbouring of v(xCor,yCor,curZ,curTime
	 * The neighbouring is a parallepiped of size (2*crossWidth+1) along X-Y and (2*crossThick+1) along Z.
	 *
	 * @param imgIn the img in
	 * @param xCor the x cor
	 * @param yCor the y cor
	 * @param curT the cur T
	 * @param totalT the total T
	 * @param curZ the cur Z
	 * @param totalZ the total Z
	 * @param nEchoesExpected the n echoes expected
	 * @param crossWidth the cross width
	 * @param crossThick the cross thick
	 * @return the full data for voxel
	 */
	public static double[][]getFullDataForVoxel(ImagePlus imgIn,int xCor,int yCor,int curT,int totalT,int curZ,int totalZ,int nEchoesExpected,int crossWidth,int crossThick){
		int xMax=imgIn.getWidth();
		int yMax=imgIn.getHeight();
		int zMax=imgIn.getNSlices();
		int xm,ym,xM,yM,zm,zM;
		xm=xCor-crossWidth;
		xM=xCor+crossWidth;
		ym=yCor-crossWidth;
		yM=yCor+crossWidth;
		zm=curZ-crossThick;
		zM=curZ+crossThick;
		xm=Math.max(xm, 0);
		xM=Math.min(xMax-1, xM);
		ym=Math.max(ym, 0);
		yM=Math.min(yMax-1, yM);
		zm=Math.max(zm, 0);
		zM=Math.min(zMax-1, zM);
		int nHits=(xM-xm+1)*(yM-ym+1)*(zM-zm+1);
		int incr=0;
		double[][]data= new double[nEchoesExpected][nHits];
		if( (xCor>xMax-1) || (yCor>yMax-1)) {IJ.log("Bad coordinates. Data set to 0"); return null;}
		for(int z=zm;z<=zM;z++) {
			for(int x=xm;x<=xM;x++) {
				for(int y=ym;y<=yM;y++) {
					for(int ec=0;ec<nEchoesExpected;ec++) {			
						int indexZ=nEchoesExpected*totalZ*curT+nEchoesExpected*z+1+ec;
						data[ec][incr]=(double)((float[])(imgIn.getStack().getProcessor(indexZ).getPixels()))[x + xMax * y];
					}
					incr++;
				}
			}
		}
		return data;
	}
	
	/**
	 *  Give the coords of the voxels in the neighbouring of v(xCor,yCor,curZ,curTime
	 * The neighbouring is a parallepiped of size (2*crossWidth+1) along X-Y and (2*crossThick+1) along Z.
	 *
	 * @param imgIn the img in
	 * @param xCor the x cor
	 * @param yCor the y cor
	 * @param curT the cur T
	 * @param totalT the total T
	 * @param curZ the cur Z
	 * @param totalZ the total Z
	 * @param nEchoesExpected the n echoes expected
	 * @param crossWidth the cross width
	 * @param crossThick the cross thick
	 * @param gaussianWeighting the gaussian weighting
	 * @return the coords of corresponding voxels used in estimation around this point
	 */
	public static int[][]getCoordsOfCorrespondingVoxelsUsedInEstimationAroundThisPoint(ImagePlus imgIn,int xCor,int yCor,int curT,int totalT,int curZ,int totalZ,int nEchoesExpected,int crossWidth,int crossThick,boolean gaussianWeighting){
		int xMax=imgIn.getWidth();
		int yMax=imgIn.getHeight();
		int zMax=imgIn.getNSlices();
		int xm,ym,xM,yM,zm,zM;
		xm=xCor-crossWidth;
		xM=xCor+crossWidth;
		ym=yCor-crossWidth;
		yM=yCor+crossWidth;
		zm=curZ-crossThick;
		zM=curZ+crossThick;
		xm=Math.max(xm, 0);
		xM=Math.min(xMax-1, xM);
		ym=Math.max(ym, 0);
		yM=Math.min(yMax-1, yM);
		zm=Math.max(zm, 0);
		zM=Math.min(zMax-1, zM);
		int nHits=(xM-xm+1)*(yM-ym+1)*(zM-zm+1);
		int incr=0;
		int[][]data= new int[nHits][8];
		if( (xCor>xMax-1) || (yCor>yMax-1)) {IJ.log("Bad coordinates. Data set to 0"); return null;}
		for(int z=zm;z<=zM;z++) {
			for(int x=xm;x<=xM;x++) {
				for(int y=ym;y<=yM;y++) {
					for(int ec=0;ec<nEchoesExpected;ec++) {			
						data[incr]=new int[] {x,y,0,0,0,0,0,0};
					}
					incr++;
				}
			}
		}
		return data;
	}

	/**
	 *  Give the average decay signal of the voxels in the neighbouring of v(xCor,yCor,curZ,curTime
	 * The neighbouring is a parallepiped of size (2*crossWidth+1) along X-Y and (2*crossThick+1) along Z
	 * 
	 * This method is surnumerary : with assumption of non using of sigma, it can be done with the other one.
	 *
	 * @param imgIn the img in
	 * @param xCor the x cor
	 * @param yCor the y cor
	 * @param curT the cur T
	 * @param totalT the total T
	 * @param curZ the cur Z
	 * @param totalZ the total Z
	 * @param nEchoesExpected the n echoes expected
	 * @param crossWidth the cross width
	 * @param crossThick the cross thick
	 * @param gaussianWeighting the gaussian weighting
	 * @return the data for voxel gaussian
	 */	
	public static double[]getDataForVoxelGaussian(ImagePlus imgIn,int xCor,int yCor,int curT,int totalT,int curZ,int totalZ,int nEchoesExpected,int crossWidth,int crossThick,boolean gaussianWeighting){
		if(gaussianWeighting) {IJ.showMessage("Not applicable : gaussian weighting. Abort");gaussianWeighting=false;}
		int xMax=imgIn.getWidth();
		int yMax=imgIn.getHeight();
		int zMax=imgIn.getNSlices();
		int zCor=curZ;
		int xm,ym,xM,yM,zm,zM;
		xm=xCor-crossWidth;
		xM=xCor+crossWidth;
		ym=yCor-crossWidth;
		yM=yCor+crossWidth;
		zm=curZ-crossThick;
		zM=curZ+crossThick;
		xm=Math.max(xm, 0);
		xM=Math.min(xMax-1, xM);
		ym=Math.max(ym, 0);
		yM=Math.min(yMax-1, yM);
		zm=Math.max(zm, 0);
		zM=Math.min(zMax-1, zM);
		double[]data= new double[nEchoesExpected];
		int nHits=(xM-xm+1)*(yM-ym+1)*(zM-zm+1);
		if( (xCor>xMax-1) || (yCor>yMax-1)) {IJ.log("Bad coordinates. Data set to 0"); return data;}
		double[][][]weights=new double[xM-xm+1][yM-ym+1][zM-zm+1];
		double sigmaX=crossWidth;
		double sigmaY=crossWidth;
		double sigmaZ=crossThick;
		
		
		double sum=0;
		for(int x=xm;x<=xM;x++) {
			for(int y=ym;y<=yM;y++) {
				for(int z=zm;z<=zM;z++) {	
					if(!gaussianWeighting)weights[x-xm][y-ym][z-zm]=1.0/nHits;
					else {
						if(sigmaX<=0 && sigmaZ<=0) {//One point
							weights[x-xm][y-ym][z-zm]=1;
						}
						else if(sigmaX<=0 && sigmaZ>0) {//Mono dimensional along z
							weights[x-xm][y-ym][z-zm]=1/(Math.pow(2*Math.PI,0.5)*sigmaZ) * Math.exp(- (z-zCor)*(z-zCor)/(sigmaZ*sigmaZ));
						}
						else if(sigmaX>0 && sigmaZ<=0) {//Two dimensional along x and y
							weights[x-xm][y-ym][z-zm]=1/(Math.pow(2*Math.PI,1)*sigmaX*sigmaY) * Math.exp(-  (x-xCor)*(x-xCor)/(sigmaX*sigmaX)  -  (y-yCor)*(y-yCor)/(sigmaY*sigmaY));
						}
						else {//Three dimensional gaussian
							weights[x-xm][y-ym][z-zm]=1/(Math.pow(2*Math.PI,1.5)*sigmaX*sigmaY*sigmaZ) * Math.exp(-  (x-xCor)*(x-xCor)/(sigmaX*sigmaX)  -  (y-yCor)*(y-yCor)/(sigmaY*sigmaY)  -  (z-zCor)*(z-zCor)/(sigmaZ*sigmaZ));
						}
					}
					sum+=weights[x-xm][y-ym][z-zm];
				}
			}
		}
		if(gaussianWeighting) {
			for(int x=xm;x<=xM;x++) {
				for(int y=ym;y<=yM;y++) {
					for(int z=zm;z<=zM;z++) {	
						weights[x-xm][y-ym][z-zm]/=sum;
					}
				}
			}
		}
		for(int ec=0;ec<nEchoesExpected;ec++) {			
			for(int z=zm;z<=zM;z++) {
				int indexZ=nEchoesExpected*totalZ*curT+nEchoesExpected*z+1+ec;
				for(int x=xm;x<=xM;x++) {
					for(int y=ym;y<=yM;y++) {
						data[ec]+=weights[x-xm][y-ym][z-zm]*(double)((float[])(imgIn.getStack().getProcessor(indexZ).getPixels()))[x + xMax * y];
					}
				}
			}
		}
		return data;
		
	}

    
	

	
	
	/**
	 *  Compute a simple T1 or T2 fit, or a T1T2 crossfit, using Simplex or Levenberg-Marquardt, 
	 *  and estimate a MONO or BI exponential, corrupted by rice Noise.
	 *
	 * @param tabTrTimesTemp the tab of recovery times
	 * @param tabTeTimesTemp the tab of echo times
	 * @param tabData the tab of magnitude data observed
	 * @param fitType the fit type to use, Ex MRUtils.T1_MONO_BIAS
	 * @param algType the alg type; Ex: MRUtils.SIMPLEX
	 * @param nbIter the nb max of iterations before early break if simplex is in use
	 * @param sigma the estimated value of Rice sigma
	 * @param debug the debug flag to display verbose output
	 * @return the object[]{estimatedParams(double[]),earlyBreak(bool)}
	 */ 
	public static Object[]makeFit(final double[]tabTrTimesTemp, final double[]tabTeTimesTemp, double[]tabData,int fitType,int algType,int nbIter,double sigma,boolean debug){
		double[]estimatedParams;
		double[]tabTrTimes=new double[tabTrTimesTemp.length];
		for(int i=0;i<tabTrTimes.length;i++)tabTrTimes[i]=tabTrTimesTemp[i];
		double[]tabTeTimes=new double[tabTeTimesTemp.length];
		for(int i=0;i<tabTeTimes.length;i++)tabTeTimes[i]=tabTeTimesTemp[i];
		boolean earlyBreak=false;
		if(algType==SIMPLEX) {
			SimplexDualCurveFitterNoBias simpfitter=new SimplexDualCurveFitterNoBias(tabTrTimes,tabTeTimes,tabData,fitType,sigma,debug);
		 	simpfitter.doFit();
			estimatedParams=simpfitter.getParams();
			estimatedParams[estimatedParams.length-1]/=( (sigma<=0 ? 1 : (sigma*sigma))*(tabTrTimes.length-getNparams(fitType)));
			earlyBreak=simpfitter.iterationsBreak;
		}	                                                
		else{
			LMDualCurveFitterNoBias lm=new LMDualCurveFitterNoBias(tabTrTimes, tabTeTimes, tabData,fitType, sigma, debug);
		 	lm.configLMA(LMDualCurveFitterNoBias.lambda,minDeltaChi2,nbIter);
			lm.doFit();
			estimatedParams=lm.getParams();
			estimatedParams[estimatedParams.length-1]/=(sigma<=0 ? 1 : (sigma*sigma));
			earlyBreak=lm.iterationsBreak;
		}

		return new Object[] {estimatedParams,earlyBreak};
	}	
	
	/**
	 *  Compute a simple T1 or T2 fit, or a T1T2 crossfit,
	 *  but compute also the distribution of estimated parameters when applying randomized modifications of the relaxation data, to check the sensibility of the operation.
	 *
	 * @param tabTrTimes the tab of recovery times
	 * @param tabTeTimes the tab of echo times
	 * @param tabData the tab of magnitude data observed
	 * @param fitType the fit type to use, Ex MRUtils.T1_MONO_BIAS
	 * @param algType the alg type; Ex: MRUtils.SIMPLEX
	 * @param nbIter the nb max of iterations before early break if simplex is in use
	 * @param sigma the estimated value of Rice sigma
	 * @param debug the debug flag to display verbose output
	 * @param nRepetMonteCarlo the n repet monte carlo
	 * @param nAveragePts the n average pts
	 * @param riceEstimator the rice estimator
	 * @return the double[][]
	 */	
	public static double[][]makeFitMonteCarlo(double[]tabTrTimes, double[]tabTeTimes, double[]tabData,int fitType,int algType,int nbIter,double sigma,int nRepetMonteCarlo,int nAveragePts,RiceEstimator riceEstimator,boolean debug){
		Random rand=new Random();
		int nParams=getNparams(fitType);
		double[]tmpSimulatedData=new double[tabData.length];
		double[][]resultsMonteCarlo=new double[nParams][nRepetMonteCarlo];
		double[]estimatedParams;
		double[]simEstimatedParams;
		double[]estimatedSigmas=new double[nParams];
		double[]estimatedMeans=new double[nParams];
		double sigmaZero;
		if(algType==LM){
			LMDualCurveFitterNoBias lmfitter=new LMDualCurveFitterNoBias(tabTrTimes,tabTeTimes,tabData,fitType,sigma,debug);
		 	lmfitter.configLMA(LMDualCurveFitterNoBias.lambda,LMDualCurveFitterNoBias.minDeltaChi2,nbIter);
		 	lmfitter.doFit();
			estimatedParams=lmfitter.getParams();
			for(int n=0;n<nRepetMonteCarlo;n++) {
				for(int dat=0;dat<tabData.length;dat++) {
					sigmaZero=riceEstimator.estimateSigma(sigma,tabData[dat])/Math.sqrt(nAveragePts); 
					tmpSimulatedData[dat]=tabData[dat]+rand.nextGaussian()*sigmaZero;
				}
				lmfitter=new LMDualCurveFitterNoBias(tabTrTimes,tabTeTimes,tmpSimulatedData,fitType,sigma,debug);
				lmfitter.doFit();
				simEstimatedParams=lmfitter.getParams();
				for(int p=0;p<nParams;p++)resultsMonteCarlo[p][n]=simEstimatedParams[p];
			}
		}
		else {
			SimplexDualCurveFitterNoBias simpfitter=new SimplexDualCurveFitterNoBias(tabTrTimes,tabTeTimes,tabData,fitType,sigma,debug);
		 	simpfitter.setIterFactor(nbIter);
		 	simpfitter.doFit();
			estimatedParams=simpfitter.getParams();
			for(int n=0;n<nRepetMonteCarlo;n++) {
				for(int dat=0;dat<tabData.length;dat++) {
					sigmaZero=riceEstimator.estimateSigma(sigma,tabData[dat])/Math.sqrt(nAveragePts); 
					tmpSimulatedData[dat]=tabData[dat]+rand.nextGaussian()*sigmaZero;
				}
				simpfitter=new SimplexDualCurveFitterNoBias(tabTrTimes,tabTeTimes,tmpSimulatedData,fitType,sigma,debug);
			 	simpfitter.setIterFactor(nbIter);
				simpfitter.doFit();
				simEstimatedParams=simpfitter.getParams();
				for(int p=0;p<nParams;p++)resultsMonteCarlo[p][n]=simEstimatedParams[p];
			}
		}
		if(nRepetMonteCarlo>0)for(int p=0;p<nParams;p++) {estimatedSigmas[p]=VitimageUtils.statistics1D(resultsMonteCarlo[p])[1];estimatedMeans[p]=VitimageUtils.statistics1D(resultsMonteCarlo[p])[0];}
		double [][]ret=new double[2][nParams];
		ret[0]=estimatedParams;
		ret[1]=estimatedSigmas;
		return ret;
	}
	
	/**
	 *  Estimate the curve corresponding to the parameters estimated, and compute the accuracy of fit, khi2, and corresponding pvalue.
	 *
	 * @param Tr the tr
	 * @param Te the te
	 * @param param the param
	 * @param sigma the sigma
	 * @param fitType the fit type
	 * @return the double[]
	 */	
	public static double[] fittenRelaxationCurve(double[]Tr,double[]Te,double []param,double sigma,int fitType){
		double[]tab=new double[Te.length];
		for(int indT=0;indT<Te.length;indT++) tab[indT]=getFitFunctionValue(Tr[indT],Te[indT],param,sigma,fitType);
		return tab;
	}

	/**
	 *  Get the number of parameters used by a model described by a fit Type,(ex: MRUTils.T1_MONO)
	 *
	 * @param fitType the fit type
	 * @return the nparams
	 */
	public static int getNparams(int fitType){
		int ret=0;
		switch(fitType) {
		case T1_MONO:ret=2;break;
		case T1_MONO_BIAS:ret=3;break;
		case T1_MONO_RICE:ret=2;break;

		case T2_MONO:ret=2;break;
		case T2_MONO_BIAS:ret=3;break;
		case T2_MONO_RICE:ret=2;break;

		case T2_MULTI:ret=4;break;
		case T2_MULTI_BIAS:ret=5;break;
		case T2_MULTI_RICE:ret=4;break;

		case T1T2_MONO:ret=3;break;
		case T1T2_MONO_BIAS:ret=4;break;
		case T1T2_MONO_RICE:ret=3;break;

		case T1T2_MULTI:ret=5;break;
		case T1T2_MULTI_BIAS:ret=6;break;
		case T1T2_MULTI_RICE:ret=5;break;

		case T1T2_DEFAULT_T2_MONO_RICE:ret=4;break;
		case T1T2_DEFAULT_T2_MULTI_RICE:ret=6;break;
		default:IJ.showMessage("In MRUtils.getNParams, unexpected fit type : "+fitType+" . \nBrutal stop now.");ret=1/0;
		}
		return ret;
	}
	
	/**
	 *  Get the value of the fitting equation of the model described by a fit Type,(ex: MRUTils.T1_MONO), given the recovery time and the echo time
	 *
	 * @param trecov the trecov
	 * @param techo the techo
	 * @param param the double[] describing the parameters of the model (see io.github.rocsg.fijirelax.curvit.LMDualCurveFitterNoBias)
	 * @param sigma the sigma of Rice noise
	 * @param fitType the fit type in use ex: MRUTils.T1_MONO)
	 * @return the fit function value
	 */
	public static double getFitFunctionValue(double trecov,double techo,double []param,double sigma,int fitType){
		double ret=0;
			switch(fitType) {
			case T1_MONO:ret=param[0]* (1-Math.exp(-(trecov / param[1])));break;
			case T1_MONO_BIAS:ret=param[0]* (1-Math.exp(-(trecov / param[1])))+param[2];break;
			case T1_MONO_RICE:ret=RiceEstimator.besFunkCost(param[0]* (1-Math.exp(-(trecov / param[1]))),sigma);break;

			case T2_MONO:ret=param[0]* Math.exp(-(techo / param[1]));break;
			case T2_MONO_BIAS:ret=param[2]+param[0]* Math.exp(-(techo / param[1]));break;
			case T2_MONO_RICE:ret=RiceEstimator.besFunkCost(param[0]* Math.exp(-(techo / param[1])),sigma);break;

			case T2_MULTI:ret= (param[0]* Math.exp(-(techo / param[2])) + param[1]* Math.exp(-(techo / param[3])  ));break;
			case T2_MULTI_BIAS:ret=param[4] + (param[0]* Math.exp(-(techo / param[2])) + param[1]* Math.exp(-(techo / param[3])  )) ;break;
			case T2_MULTI_RICE:ret=RiceEstimator.besFunkCost(  (param[0]* Math.exp(-(techo / param[2])) + param[1]* Math.exp(-(techo / param[3])  )) ,sigma);break;

			case T1T2_MONO:ret=param[0]* (1-Math.exp(-(trecov / param[1]))) * Math.exp(-(techo / param[2]));break;
			case T1T2_MONO_BIAS:ret=param[3]+ param[0]* (1-Math.exp(-(trecov / param[1]))) * Math.exp(-(techo / param[2]));break;
			case T1T2_MONO_RICE:ret=RiceEstimator.besFunkCost(param[0]* (1-Math.exp(-(trecov / param[1]))) * Math.exp(-(techo / param[2])),sigma);break;

			case T1T2_MULTI:ret= (param[0]*Math.exp(-(techo / param[3]) )+param[1]*Math.exp(-(techo / param[4])) ) * (1-Math.exp(-(trecov / param[2])) ) ;break;
			case T1T2_MULTI_BIAS:ret=param[5]+(param[0]*Math.exp(-(techo / param[3]) )+param[1]*Math.exp(-(techo / param[4])) ) * (1-Math.exp(-(trecov / param[2])) ) ;break;
			case T1T2_MULTI_RICE:ret=RiceEstimator.besFunkCost(  (param[0]*Math.exp(-(techo / param[3]) )+param[1]*Math.exp(-(techo / param[4])) ) * (1-Math.exp(-(trecov / param[2])) )  , sigma);break;

			case T1T2_DEFAULT_T2_MONO_RICE:ret=RiceEstimator.besFunkCost(param[0]* (1-Math.exp(-(trecov / param[1]))) * Math.exp(-((techo+(trecov>9999 ? param[3] : 0)) / param[2])  ),sigma);break;
			case T1T2_DEFAULT_T2_MULTI_RICE:ret=RiceEstimator.besFunkCost(  (param[0]*Math.exp(-(  (techo+(trecov>9999 ? param[5] : 0))   / param[3]) )+param[1]*Math.exp(-((techo+(trecov>9999 ? param[5] : 0)) / param[4])) * (1-Math.exp(-((trecov) / param[2])) ) ) , sigma);break;
			default:IJ.showMessage("In MRUtils.getFitFunctionValue, unexpected fit type : "+fitType+" . \nBrutal stop now.");ret=1/0;
			}
		return ret;
	}

	/**
	 * Compute normalised khi2 and p values according to the sigma parameter of the Rice noise
	 * the number of points averaged for the estimation is used to estimate sigma of the sum
	 * Another interesting point is that sigma is reestimated locally for each magnitude value
	 * as the rice sigma is not equal to the observed sigma for different values of the magnitude.
	 *
	 * @param tabData the tab data
	 * @param tabTrTimes the tab of recovery times
	 * @param tabTeTimes the tab of echo times
	 * @param sigma the sigma of Rice noise
	 * @param estimatedParams the estimated parameters of the fitting function
	 * @param fitType the fit type (ex: T2_MONO_BIAS)
	 * @param debug the debug flag in order to display some debugging infformations
	 * @param riceEstimator the RiceEstimator object in use to handle the analysis
	 * @param nbPoints the nb of points averaged for getting this data
	 * @return the double[]{khi2, p}
	 */
	public static double[] fittingAccuracies(double[] tabData,double[] tabTrTimes,double[] tabTeTimes,double sigma,double[] estimatedParams,int fitType,boolean debug,RiceEstimator riceEstimator,int nbPoints){		
		int N_PARAMS=getNparams(fitType);
		double []realP=fittenRelaxationCurve(tabTrTimes,tabTeTimes,estimatedParams,sigma,fitType);
		double khi2=0;		double diff=0;		double sigmaEst=0;
		for(int indT=0;indT<tabData.length;indT++) {
			sigmaEst=riceEstimator.estimateSigma(sigma, realP[indT])/Math.sqrt(nbPoints);
			diff=(realP[indT]-tabData[indT]);
			khi2+=diff*diff/(sigmaEst*sigmaEst);
		}
		double pVal=getPvalue(khi2,tabTrTimes.length-N_PARAMS);
		double khi2Ret=khi2/(tabTrTimes.length-N_PARAMS);
 		return new double[] {khi2Ret,pVal*100};
	}
	
	/**
	 * Gets the pvalue associated with the khi2 and the number of freedom degrees
	 *
	 * @param khi2 the khi 2
	 * @param freedomDegrees the freedom degrees
	 * @return the pvalue
	 */
	public static double getPvalue(double khi2,int freedomDegrees){
		if (freedomDegrees<1)return MRUtils.infinity;
		ChiSquaredDistribution x2 = new ChiSquaredDistribution( freedomDegrees );
		try {return(x2.cumulativeProbability(khi2));} catch (Exception e) {e.printStackTrace();}
		return -1;
	}
	
	/**
	 * Convert data from custom MRImaging device in percentage of maximal magnitude
	 *
	 * @param m0 the m 0
	 * @param doIt the do it
	 * @return the double
	 */
	public static double convertBionanoM0InPercentage(double m0,boolean doIt) {
		if(doIt)return (m0/maxM0BionanoForNormalization);
		else return m0;
	}
	
	/**
	 * Copy a boolean tab
	 *
	 * @param tab the tab
	 * @return the boolean[]
	 */
	public static boolean[]copyTab(boolean[]tab){
		boolean[]tab2=new boolean[tab.length];
		for(int i=0;i<tab.length;i++)tab2[i]=tab[i];
		return tab2;
	}

	/**
	 * Copy double tab
	 *
	 * @param tab the tab
	 * @return the double[]
	 */
	public static double[]copyTab(double[]tab){
		double[]tab2=new double[tab.length];
		for(int i=0;i<tab.length;i++)tab2[i]=tab[i];
		return tab2;
	}
	
	/**
	 * Copy int tab
	 *
	 * @param tab the tab
	 * @return the int[]
	 */
	public static int[]copyTab(int[]tab){
		int[]tab2=new int[tab.length];
		for(int i=0;i<tab.length;i++)tab2[i]=tab[i];
		return tab2;
	}
	

	
	
	/**
	 * Bes funk cost. Application of the rice noise corruption to the initial signal (see FijiRelax paper)
	 *
	 * @param value the value
	 * @param sigma the sigma
	 * @return the double
	 */

	public static double besFunkCost(double value,double sigma) {
		double alpha=value*value/(4*sigma*sigma);
		return Math.sqrt(Math.PI*sigma*sigma/2.0)*( (1+2*alpha)*bessi0NoExp(alpha) + 2*alpha*bessi1NoExp(alpha) );
	}	

	/**
	 * Bes funk deriv. derivative of the bessel function on the corresponding value
	 *
	 * @param value the value
	 * @param sigma the sigma
	 * @return the double
	 */
	public static double besFunkDeriv(double value,double sigma) {
		double alpha=value*value/(4*sigma*sigma);
		return Math.sqrt(Math.PI*alpha/2.0)*(bessi0NoExp(alpha) + bessi1NoExp(alpha) );
	}

	/**
	 * Bessi 0 no exp. Fonction de bessel modifiee de premiere espece d'ordre 0
	 *
	 * @param x the x
	 * @return the double
	 */
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

	/**
	 * Bessi 1 no exp. Fonction de bessel modifiee de premiere espece d'ordre 1
	 *
	 * @param x the x
	 * @return the double
	 */
	//
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

	
	
	
	
	
	
	
	
	/**
	 *  A routine to normalize a z-stack along the z-axis, to compensate for the non-uniformity of the magnetic field 
	 * This method detects the capillary and do 
	 * the normalization of the capillary area and the rest of the image separately
	 *
	 * @param imgsEchoes the succesive echo images
	 * @param imgMask the image of mask where to estimate the fit or not
	 * @return the image plus[]
	 */
	public static ImagePlus[] normalizeBeforeComputation(ImagePlus []imgsEchoes,ImagePlus imgMask) {
		ImagePlus []tabRet=VitimageUtils.imageTabCopy(imgsEchoes);
		ImagePlus imgSum=VitimageUtils.sumOfImages(imgsEchoes);
		imgSum=VitimageUtils.makeOperationBetweenTwoImages(imgSum, imgMask, 2, true);
		int Z=imgSum.getNSlices();
		int X=VitimageUtils.getDimensions(imgSum)[0];
		int Y=VitimageUtils.getDimensions(imgSum)[1];
		int[][]coords=new int[Z][3];
		double[]meanOverAll=new double[Z];
		double[]meanOverCap=new double[Z];
		double[]meanOverRest=new double[Z];
		double globalMeanOverCap;
		double globalMeanOverRest;
		int radiusAroundCapillary=(int)Math.round(VitimageUtils.bionanoCapillaryRadius*1.5/VitimageUtils.getVoxelSizes(imgSum)[0]);
		//Pour chaque Z
		for(int z=0;z<Z;z++) {
			ImagePlus sliceTemp=new Duplicator().run(imgSum,1,1,z+1,z+1,1,1);
			//Pour cette image, localiser le capillaire, sur chaque Z
			coords[z]=VitimageUtils.findCapillaryCenterInSlice(sliceTemp, VitimageUtils.bionanoCapillaryRadius);
			if(coords[z][0]<radiusAroundCapillary+5 || coords[z][1]<radiusAroundCapillary+5 || coords[z][0]>(X-radiusAroundCapillary-5) || coords[z][0]>(Y-radiusAroundCapillary-5) ) {
				IJ.showMessage("Warning in normalizeBeforeComputation : unsteady mean computation at z="+z+" because capillary is detected near the border, at coordinates "+coords[z][0]+" , "+coords[z][1]+" with an indicated security radius of "+radiusAroundCapillary);
			}
			
			//Mesurer la moyenne autour du capillaire
			meanOverCap[z]=VitimageUtils.meanValueofImageAround(sliceTemp, coords[z][0], coords[z][1], 0,radiusAroundCapillary);
			
			
			//Mesurer la moyenne de l'image
			int radiusImg=Math.min(X/2,Y/2)-2;
			meanOverAll[z]=VitimageUtils.meanValueofImageAround(sliceTemp, X/2, Y/2, 0,radiusImg);
			
			//Oter l'un à l'autre
			meanOverRest[z]=meanOverAll[z]*(2*radiusImg+1)*(2*radiusImg+1) - meanOverCap[z]*(2*radiusAroundCapillary+1)*(2*radiusAroundCapillary+1);
			meanOverRest[z]=meanOverRest[z]/( (2*radiusImg+1)*(2*radiusImg+1) - (2*radiusAroundCapillary+1)*(2*radiusAroundCapillary+1) );

		}
		//Faire la moyenne des valeurs mesurées along Z
		globalMeanOverRest=VitimageUtils.statistics1D(meanOverRest)[0];
		globalMeanOverCap=VitimageUtils.statistics1D(meanOverCap)[0];

		int radiusSq=radiusAroundCapillary*radiusAroundCapillary;
		double distSq;
		//Pour chaque Z
		for(int im=0;im<tabRet.length;im++) {
			tabRet[im]=VitimageUtils.convertToFloat(tabRet[im]);
			for(int z=0;z<Z;z++) {
				float[]tabData=(float[])tabRet[im].getStack().getProcessor(z+1).getPixels();
				double factorMultCap=globalMeanOverCap/meanOverCap[z];
				double factorMultRest=globalMeanOverRest/meanOverRest[z];
				for(int x=0;x<X;x++) {
					for(int y=0;y<Y;y++) {
						distSq=(coords[z][0]-x)*(coords[z][0]-x) + (coords[z][1]-y)*(coords[z][1]-y);
						//Multiplier par moy / facteur 1(Z) sauf dans cap ou on multiplie par moy / facteur 2(Z) 
						if(distSq<radiusSq)tabData[X*y+x]*=(factorMultCap);
						else tabData[X*y+x]*=(factorMultRest);
					}
				}
			}
		}					
		return tabRet;
	}
	
	
    
    
}