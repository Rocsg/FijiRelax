package com.vitimage.mrutils;


import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;
import java.util.concurrent.atomic.AtomicInteger;

import org.apache.commons.math3.distribution.ChiSquaredDistribution;

import com.vitimage.common.TransformUtils;
import com.vitimage.common.VitiDialogs;
import com.vitimage.common.VitimageUtils;

import ij.IJ;
import ij.ImagePlus;
import ij.plugin.Duplicator;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
public class MRUtils  {
//	static float stdValMaxIRM=50000;
	//	public static float stdSigmaIRM=159;
	
	public MRUtils() {
		
	}


	/** Top level entry points : compute maps from the successive images of the echo or repetition times*/
	public static ImagePlus[] computeT2Map(final ImagePlus[]imgTe,double sigmaSmoothing) {
		final int nEch=imgTe.length;
		ImagePlus[]imgTeTemp=VitimageUtils.imageTabCopy(imgTe);
		double[]voxs=VitimageUtils.getVoxelSizes(imgTeTemp[0]);
		int[]dims=VitimageUtils.getDimensions(imgTeTemp[0]);
		
		for(int i=0;i<nEch ; i++) {
			imgTeTemp[i]=VitimageUtils.gaussianFiltering(imgTeTemp[i],voxs[0]*sigmaSmoothing,voxs[1]*sigmaSmoothing,voxs[0]*sigmaSmoothing);//It's no error : no "big smoothing" over Z, due to misalignment
		}
		ImagePlus imgM0=VitimageUtils.imageCopy(imgTeTemp[0]);
		IJ.run(imgM0,"32-bit","");
		ImagePlus imgT2=VitimageUtils.imageCopy(imgM0);
		double factorSeconds=8E-5;
		final int fitType=MRUtils.T2_RELAX_RICE;
		final int algType=MRUtils.SIMPLEX;
		final int X=dims[0];		final int Y=dims[1];		final int Z=dims[2];
		final double []listTeForThreads=getTeFrom3DRelaxationImageTab(imgTe)[0];
		final double []listSigmaForThreads=getSigmaFrom3DRelaxationImageTab(imgTe);
		final float[][][]tabData=new float[nEch][Z][];
		for(int i=0;i<nEch ;i++)for(int z=0;z<Z ; z++) tabData[i][z]=(float[])imgTeTemp[i].getStack().getProcessor(z+1).getPixels();
		final FloatProcessor[] tempComputedT2=new FloatProcessor[Z];
		final FloatProcessor[] tempComputedM0=new FloatProcessor[Z];
		final AtomicInteger incrZ = new AtomicInteger(0);
		final AtomicInteger incrProcess = new AtomicInteger(0);
		final int totLines=Y*Z;
		final int twopercent=totLines/50;
		IJ.log("Multi-threaded T2 map computation. start fit on "+(X*Y*Z)+" voxels.\n--> Estimated time  @2.5Ghz @12 cores : "+VitimageUtils.dou(factorSeconds*(X*Y*Z) )+" s");
		final Thread[] threads = VitimageUtils.newThreadArray(Z);    
		for (int ithread = 0; ithread < Z; ithread++) {  
			threads[ithread] = new Thread() {  { setPriority(Thread.NORM_PRIORITY); }  
			public void run() {  
				FloatProcessor threadCompT2=new FloatProcessor(X,Y);
				FloatProcessor threadCompM0=new FloatProcessor(X,Y);
				float[]tabThreadT2=(float[])threadCompT2.getPixels();
				float[]tabThreadM0=(float[])threadCompM0.getPixels();
				
				double[]echoesForThisVoxel=new double[nEch];		
				double[]estimatedParams;
				int z = incrZ.getAndIncrement();

				for(int y=0;y<Y;y++) {
					int incrTot=incrProcess.getAndIncrement();
					if(incrTot%twopercent==0)IJ.log("Processing map :, "+VitimageUtils.dou(incrTot*100.0/totLines)+" %");
					for(int x=0;x<X;x++) {
						int index=y*X+x;
						for(int ech=0;ech<nEch;ech++)echoesForThisVoxel[ech]= tabData[ech][z][index];
						estimatedParams=MRUtils.makeFitSimple(listTeForThreads, echoesForThisVoxel,fitType,algType,400,listSigmaForThreads[z]);
						tabThreadT2[index]=(float)estimatedParams[1];
						tabThreadM0[index]=(float)estimatedParams[0];
						if(tabThreadM0[index]<0 || tabThreadM0[index]>(2*echoesForThisVoxel[0]) 
								|| ((echoesForThisVoxel[0])<3*listSigmaForThreads[z])
								|| tabThreadT2[index]<5 || tabThreadT2[index]>1000  || ((tabThreadT2[index]>200) && (tabThreadM0[index]<6*listSigmaForThreads[z] ))) {
							tabThreadM0[index]=(float)(echoesForThisVoxel[0]);
							tabThreadT2[index]=0;
						}
					}
				}
				tempComputedT2[z]=threadCompT2;
				tempComputedM0[z]=threadCompM0;
			}};  
		}  		
		VitimageUtils.startAndJoin(threads);  
		for (int z=0; z< Z; z++) {  
			imgM0.getStack().setProcessor(tempComputedM0[z], z+1);
			imgT2.getStack().setProcessor(tempComputedT2[z], z+1);
		}  
		imgM0.setDisplayRange(0, 500);
		imgT2.setDisplayRange(0, 200);
		return new ImagePlus[] {imgM0,imgT2};
	}  
	
	public static ImagePlus[] computeT1Map(final ImagePlus[]imgTr,double sigmaSmoothing) {
		final int nEch=imgTr.length;
		ImagePlus[]imgTrTemp=VitimageUtils.imageTabCopy(imgTr);
		double[]voxs=VitimageUtils.getVoxelSizes(imgTrTemp[0]);
		int[]dims=VitimageUtils.getDimensions(imgTrTemp[0]);
		
		for(int i=0;i<nEch ; i++) {
			imgTrTemp[i]=VitimageUtils.gaussianFiltering(imgTrTemp[i],voxs[0]*sigmaSmoothing,voxs[1]*sigmaSmoothing,voxs[0]*sigmaSmoothing);//It's no error : no "big smoothing" over Z, due to misalignment
		}
		ImagePlus imgM0=VitimageUtils.imageCopy(imgTrTemp[0]);
		IJ.run(imgM0,"32-bit","");
		ImagePlus imgT1=VitimageUtils.imageCopy(imgM0);
		final int fitType=MRUtils.T1_RECOVERY_RICE;
		final int algType=MRUtils.SIMPLEX;
		final int X=dims[0];		final int Y=dims[1];		final int Z=dims[2];
		final double []listTrForThreads=getTrFrom3DRelaxationImageTab(imgTr);
		final double []listSigmaForThreads=getSigmaFrom3DRelaxationImageTab(imgTr);
				
		final float[][][]tabData=new float[nEch][Z][];
		for(int i=0;i<nEch ;i++)for(int z=0;z<Z ; z++) tabData[i][z]=(float[])imgTrTemp[i].getStack().getProcessor(z+1).getPixels();
		final FloatProcessor[] tempComputedT1=new FloatProcessor[Z];
		final FloatProcessor[] tempComputedM0=new FloatProcessor[Z];
		final AtomicInteger incrZ = new AtomicInteger(0);
		final AtomicInteger incrProcess = new AtomicInteger(0);
		final int totLines=Y*Z;
		final int twoPercent=totLines/50;
		IJ.log("Multi-threaded T1 map computation. start fit on "+(X*Y*Z)+" voxels with sigma="+listSigmaForThreads[0]+" and times="+TransformUtils.stringVectorN(listTrForThreads,""));
		final Thread[] threads = VitimageUtils.newThreadArray(Z);    
		for (int ithread = 0; ithread < Z; ithread++) {  
			threads[ithread] = new Thread() {  { setPriority(Thread.NORM_PRIORITY); }  
			public void run() {  
				FloatProcessor threadCompT1=new FloatProcessor(X,Y);
				FloatProcessor threadCompM0=new FloatProcessor(X,Y);
				float[]tabThreadT1=(float[])threadCompT1.getPixels();
				float[]tabThreadM0=(float[])threadCompM0.getPixels();
				
				double[]echoesForThisVoxel=new double[nEch];		
				double[]estimatedParams;
				int z = incrZ.getAndIncrement();

				for(int y=0;y<Y;y++) {
					int incrTot=incrProcess.getAndIncrement();
					if(incrTot%twoPercent==0)IJ.log("Processing map :, "+VitimageUtils.dou(incrTot*100.0/totLines)+" %");
					for(int x=0;x<X;x++) {
						boolean debug=false;
//						if(x==262 && y==160)debug=true;
						int index=y*X+x;
						for(int ech=0;ech<nEch;ech++)echoesForThisVoxel[ech]= tabData[ech][z][index];
						estimatedParams=MRUtils.makeFitSimple(listTrForThreads, echoesForThisVoxel,fitType,algType,400,listSigmaForThreads[z]);

						
						tabThreadT1[index]=(float)estimatedParams[1];
						tabThreadM0[index]=(float)estimatedParams[0];
						if(tabThreadM0[index]<0 || tabThreadM0[index]>10*echoesForThisVoxel[nEch-1] || ((echoesForThisVoxel[nEch-3]+echoesForThisVoxel[nEch-1]+echoesForThisVoxel[nEch-2])<9*listSigmaForThreads[z])
								|| tabThreadT1[index]<0 || tabThreadT1[index]>3300) {
							tabThreadM0[index]=(float)(echoesForThisVoxel[nEch-1]);
							tabThreadT1[index]=0;
						}
					}
				}
				tempComputedT1[z]=threadCompT1;
				tempComputedM0[z]=threadCompM0;
			}};  
		}  		
		VitimageUtils.startAndJoin(threads);  
		for (int z=0; z< Z; z++) {  
			imgM0.getStack().setProcessor(tempComputedM0[z], z+1);
			imgT1.getStack().setProcessor(tempComputedT1[z], z+1);
		}  
		imgM0.setDisplayRange(0, 5000);
		imgT1.setDisplayRange(0, 4000);
		return new ImagePlus[] {imgM0,imgT1};
	}  



	
	
	
	
	
	
	
	public static ImagePlus[] computeT1T2Map(final ImagePlus[]imgs,double sigmaSmoothing,int fitType) {
		double[]voxs=VitimageUtils.getVoxelSizes(imgs[0]);
		int[]dims=VitimageUtils.getDimensions(imgs[0]);
		int n=imgs.length;
		ImagePlus []imgsTemp=new ImagePlus[n];
		for(int i=0;i<n;i++)imgsTemp[i]=imgs[i].duplicate();
		
		for(int i=0;i<n ; i++) {
			imgsTemp[i]=VitimageUtils.gaussianFiltering(imgsTemp[i],voxs[0]*sigmaSmoothing,voxs[1]*sigmaSmoothing,voxs[0]*sigmaSmoothing);//It's no error : no "big smoothing" over Z, due to misalignment
		}
		ImagePlus imgM0=imgsTemp[0].duplicate();
		IJ.run(imgM0,"32-bit","");
		ImagePlus imgT1=imgM0.duplicate();
		ImagePlus imgT2=imgM0.duplicate();
		ImagePlus imgT12=imgM0.duplicate();
		ImagePlus imgT22=imgM0.duplicate();
		ImagePlus imgM02=imgM0.duplicate();
		final int algType=MRUtils.SIMPLEX;
		final int X=dims[0];		final int Y=dims[1];		final int Z=dims[2];
		final double []listTrForThreads=getTrFrom3DRelaxationImageTab(imgsTemp);
		final double [][]listTeForThreads=getTeFrom3DRelaxationImageTab(imgsTemp);
		final double []listSigmaForThreads=getSigmaFrom3DRelaxationImageTab(imgsTemp);
				
		final float[][][]tabData=new float[n][Z][];
		final float[][][]tabDataNonSmooth=new float[n][Z][];
		for(int i=0;i<n ;i++)for(int z=0;z<Z ; z++) {
			tabData[i][z]=(float[])imgsTemp[i].getStack().getProcessor(z+1).getPixels();
			tabDataNonSmooth[i][z]=(float[])imgs[i].getStack().getProcessor(z+1).getPixels();
		}
		final FloatProcessor[] tempComputedM0=new FloatProcessor[Z];
		final FloatProcessor[] tempComputedT1=new FloatProcessor[Z];
		final FloatProcessor[] tempComputedT2=new FloatProcessor[Z];
		final FloatProcessor[] tempComputedM02=new FloatProcessor[Z];
		final FloatProcessor[] tempComputedT12=new FloatProcessor[Z];
		final FloatProcessor[] tempComputedT22=new FloatProcessor[Z];
		final AtomicInteger incrZ = new AtomicInteger(0);
		final AtomicInteger incrProcess = new AtomicInteger(0);
		final int totLines=Y*Z;
		final int nThread=n;
		final int onePercent=totLines/100;
		IJ.log("Multi-threaded T1 and T2 map computation. start fit on "+(X*Y*Z)+" voxels with sigma="+listSigmaForThreads[0]+" and times="+TransformUtils.stringVectorN(listTrForThreads,""));
		final Thread[] threads = VitimageUtils.newThreadArray(Z);    
		for (int ithread = 0; ithread < Z; ithread++) {  
			threads[ithread] = new Thread() {  { setPriority(Thread.NORM_PRIORITY); }  
			public void run() {  
				FloatProcessor threadCompT2=new FloatProcessor(X,Y);
				FloatProcessor threadCompT1=new FloatProcessor(X,Y);
				FloatProcessor threadCompM0=new FloatProcessor(X,Y);
				float[]tabThreadT1=(float[])threadCompT1.getPixels();
				float[]tabThreadM0=(float[])threadCompM0.getPixels();
				float[]tabThreadT2=(float[])threadCompT2.getPixels();
				FloatProcessor threadCompT22=new FloatProcessor(X,Y);
				FloatProcessor threadCompT12=new FloatProcessor(X,Y);
				FloatProcessor threadCompM02=new FloatProcessor(X,Y);
				float[]tabThreadT12=(float[])threadCompT12.getPixels();
				float[]tabThreadM02=(float[])threadCompM02.getPixels();
				float[]tabThreadT22=(float[])threadCompT22.getPixels();
				
				double[]echoesForThisVoxel=new double[nThread];		
				double[]echoesNonGaussForThisVoxel=new double[nThread];		
				double[]estimatedParams;
				int z = incrZ.getAndIncrement();

				for(int y=0;y<Y;y++) {
					int incrTot=incrProcess.getAndIncrement();
					if(incrTot%onePercent==0)IJ.log("Processing map :, "+VitimageUtils.dou(incrTot*100.0/totLines)+" %");
					for(int x=0;x<X;x++) {
						int index=y*X+x;
						for(int ech=0;ech<nThread;ech++) {
							echoesForThisVoxel[ech]= tabData[ech][z][index];
							echoesNonGaussForThisVoxel[ech]= tabDataNonSmooth[ech][z][index];
						}
						if(VitimageUtils.moyTwoMax(echoesNonGaussForThisVoxel)<THRESHOLD_RATIO_BETWEEN_MAX_ECHO_AND_SIGMA_RICE*listSigmaForThreads[z]) {
							tabThreadT2[index]=0;
							tabThreadT1[index]=0;
							tabThreadM0[index]=0;
							if(fitType==T1T2_MULTI_RICE) {
								tabThreadT22[index]=0;
								tabThreadM02[index]=0;
							}
							if(fitType==T1T2_MULTIMULTI_RICE) {
								tabThreadT22[index]=0;
								tabThreadM02[index]=0;
								tabThreadT12[index]=0;
							}
							continue;
						}
						estimatedParams=MRUtils.makeFitDouble(listTrForThreads, listTeForThreads[z],echoesForThisVoxel,fitType,algType,N_ITER_T1T2,listSigmaForThreads[z]);

						double max=VitimageUtils.max(echoesForThisVoxel);
						

						
						if(fitType==T1T2_MULTIMULTI_RICE) {
							tabThreadM0[index]=(float)estimatedParams[0];
							tabThreadM02[index]=(float)estimatedParams[1];
							tabThreadT1[index]=(float)estimatedParams[2];
							tabThreadT12[index]=(float)estimatedParams[3];
							tabThreadT2[index]=(float)estimatedParams[4];
							tabThreadT22[index]=(float)estimatedParams[5];
							if(tabThreadM0[index]<0 || tabThreadM0[index]>3*max || tabThreadM02[index]<0 || tabThreadM02[index]>3*max  || (VitimageUtils.max(echoesForThisVoxel)<THRESHOLD_RATIO_BETWEEN_MAX_ECHO_AND_SIGMA_RICE*listSigmaForThreads[z])
									|| tabThreadT1[index]<0 || tabThreadT1[index]>4000 || tabThreadT2[index]<5 || tabThreadT2[index]>500 
									|| tabThreadT12[index]<0 || tabThreadT12[index]>4000 || tabThreadT22[index]<5 || tabThreadT22[index]>500) {
								tabThreadM0[index]=(float) max;
								tabThreadM02[index]=(float) max;
								tabThreadT1[index]=0;
								tabThreadT2[index]=0;
								tabThreadT12[index]=0;
								tabThreadT22[index]=0;
							}
						}
						else if(fitType==T1T2_MULTI_RICE) {
							tabThreadM0[index]=(float)estimatedParams[0];
							tabThreadM02[index]=(float)estimatedParams[1];
							tabThreadT1[index]=(float)estimatedParams[2];
							tabThreadT2[index]=(float)estimatedParams[3];
							tabThreadT22[index]=(float)estimatedParams[4];
							if(tabThreadM0[index]<0 || tabThreadM0[index]>3*max || tabThreadM02[index]<0 || tabThreadM02[index]>3*max  || (VitimageUtils.max(echoesForThisVoxel)<THRESHOLD_RATIO_BETWEEN_MAX_ECHO_AND_SIGMA_RICE*listSigmaForThreads[z])
									|| tabThreadT1[index]<0 || tabThreadT1[index]>4000 || tabThreadT2[index]<5 || tabThreadT2[index]>500 
									|| tabThreadT22[index]<5 || tabThreadT22[index]>500) {
								tabThreadM0[index]=(float) max;
								tabThreadM02[index]=(float) max;
								tabThreadT1[index]=0;
								tabThreadT2[index]=0;
								tabThreadT12[index]=0;
								tabThreadT22[index]=0;
							}
						}
						else {
							tabThreadM0[index]=(float)estimatedParams[0];
							tabThreadT1[index]=(float)estimatedParams[1];
							tabThreadT2[index]=(float)estimatedParams[2];
							if(tabThreadM0[index]<0 || tabThreadM0[index]>3*max || (VitimageUtils.max(echoesForThisVoxel)<THRESHOLD_RATIO_BETWEEN_MAX_ECHO_AND_SIGMA_RICE*listSigmaForThreads[z])|| tabThreadM0[index]<THRESHOLD_RATIO_BETWEEN_MAX_ECHO_AND_SIGMA_RICE*listSigmaForThreads[z]
								|| tabThreadT1[index]<0 || tabThreadT1[index]>4000 || tabThreadT2[index]<5 || tabThreadT2[index]>500) {
							tabThreadM0[index]=(float) max;
							tabThreadT1[index]=0;
							tabThreadT2[index]=0;
							}
						}
					}
				}
				tempComputedT2[z]=threadCompT2;
				tempComputedT1[z]=threadCompT1;
				tempComputedM0[z]=threadCompM0;
				if(fitType==T1T2_MULTIMULTI_RICE) {
					tempComputedT22[z]=threadCompT22;
					tempComputedT12[z]=threadCompT12;
					tempComputedM02[z]=threadCompM02;					
				}
				if(fitType==T1T2_MULTI_RICE) {
					tempComputedT22[z]=threadCompT22;
					tempComputedM02[z]=threadCompM02;
					
				}
			}};  
		}  		
		VitimageUtils.startAndJoin(threads);  
		for (int z=0; z< Z; z++) {  
			imgM0.getStack().setProcessor(tempComputedM0[z], z+1);
			imgT1.getStack().setProcessor(tempComputedT1[z], z+1);
			imgT2.getStack().setProcessor(tempComputedT2[z], z+1);
			if(fitType==T1T2_MULTIMULTI_RICE) {
				imgM02.getStack().setProcessor(tempComputedM02[z], z+1);
				imgT12.getStack().setProcessor(tempComputedT12[z], z+1);
				imgT22.getStack().setProcessor(tempComputedT22[z], z+1);
			}
			if(fitType==T1T2_MULTI_RICE) {
				imgM02.getStack().setProcessor(tempComputedM02[z], z+1);
				imgT22.getStack().setProcessor(tempComputedT22[z], z+1);
			}
		}  
		imgM0.setDisplayRange(0, maxDisplayedM0);
		imgT1.setDisplayRange(0, maxDisplayedT1);
		imgT2.setDisplayRange(0, maxDisplayedT2);
		if(fitType==T1T2_MULTIMULTI_RICE) {
			imgM02.setDisplayRange(0, maxDisplayedM0);
			imgT12.setDisplayRange(0, maxDisplayedT1);
			imgT22.setDisplayRange(0, maxDisplayedT2);
			return new ImagePlus[] {imgM0,imgM02,imgT1,imgT12,imgT2,imgT22};
		}
		else if(fitType==T1T2_MULTI_RICE) {
			imgM02.setDisplayRange(0, maxDisplayedM0);
			imgT22.setDisplayRange(0, maxDisplayedT2);
			return new ImagePlus[] {imgM0,imgM02,imgT1,imgT2,imgT22};
		}
		else return new ImagePlus[] {imgM0,imgT1,imgT2};
	}  

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	public static double[] doubleArraySort(double[]tab) {
		double[]tabRet=new double[tab.length];
		ArrayList<Double> list=new ArrayList<Double>();
		for(double dou : tab)list.add(dou);
		Collections.sort(list);
		for(int i=0;i<list.size();i++)tabRet[i]=list.get(i);
		return tabRet;
	}


	
	
	/** Estimate rice noise in each slice from the background values*/
	public static void computeTeTrAndRiceSigmaOfEachSliceAndWriteItInTheLabels(ImagePlus img,boolean keepOldString,String optionalAdditionalString,boolean makeBoutureTrick) {
		int nbC=img.getNChannels();
		int nbZ=img.getNSlices();
		int nbF=img.getNFrames();
		ImageProcessor imgP;
		double tr=VitimageUtils.getRepetitionTime(img);
		double te=VitimageUtils.getEchoTime(img);
		boolean isBoutureT2=VitimageUtils.isBouture(img) && tr>9000;
		System.out.println("Make bouture trick ?"+makeBoutureTrick);
		for(int nc=0;nc<nbC;nc++) {
			for(int nz=0;nz<nbZ;nz++) {
				for(int nf=0;nf<nbF;nf++) {
					int sli=VitimageUtils.getCorrespondingSliceInHyperImage(img,nc, nz, nf);
					imgP=img.getStack().getProcessor(sli);
					double[]stats=getBackgroundStatsFromProcessor(imgP,3);
					double sigmaRice=RiceEstimator.computeRiceSigmaFromBackgroundValuesStatic(stats[0],stats[1]);
					String chain="_TR="+VitimageUtils.dou(tr)+"_TE="+VitimageUtils.dou(te+ (((makeBoutureTrick && tr>9000) ? (nz==0 ? 0 : MRUtils.DELTA_TE_BOUTURE_TRICK) : 0 )));
					System.out.println("Setting chain "+chain);
					img.getStack().setSliceLabel((keepOldString ? img.getStack().getSliceLabel(sli) : "") + "_"+optionalAdditionalString+chain+"_SIGMARICE="+VitimageUtils.dou(sigmaRice), sli);
				}
			}
		}
	}
	
	public static double[]getBackgroundStatsFromProcessor(ImageProcessor imgP,int config0firstCorner_config1Corners_config2Cross_config3Both) {
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
		double[][] vals=null;
		if(config0firstCorner_config1Corners_config2Cross_config3Both==0) {
			vals=new double[1][];
		}
		else if(config0firstCorner_config1Corners_config2Cross_config3Both==1) {
			vals=new double[4][];
			vals[0]=VitimageUtils.valuesOfImageProcessor(imgP,x0,y0,samplSize/2);
			vals[1]=VitimageUtils.valuesOfImageProcessor(imgP,x0,y2,samplSize/2);
			vals[2]=VitimageUtils.valuesOfImageProcessor(imgP,x2,y0,samplSize/2);
			vals[3]=VitimageUtils.valuesOfImageProcessor(imgP,x2,y2,samplSize/2);		
		}
		else if(config0firstCorner_config1Corners_config2Cross_config3Both==2) {
			vals=new double[4][];
			vals[0]=VitimageUtils.valuesOfImageProcessor(imgP,x0,y1,samplSize/4);
			vals[1]=VitimageUtils.valuesOfImageProcessor(imgP,x1,y0,samplSize/4);
			vals[2]=VitimageUtils.valuesOfImageProcessor(imgP,x2,y1,samplSize/4);
			vals[3]=VitimageUtils.valuesOfImageProcessor(imgP,x1,y2,samplSize/4);		
		}
		else if(config0firstCorner_config1Corners_config2Cross_config3Both==3) {
			vals=new double[8][];
			vals[0]=VitimageUtils.valuesOfImageProcessor(imgP,x0,y0,samplSize/2);
			vals[1]=VitimageUtils.valuesOfImageProcessor(imgP,x0,y2,samplSize/2);
			vals[2]=VitimageUtils.valuesOfImageProcessor(imgP,x2,y0,samplSize/2);
			vals[3]=VitimageUtils.valuesOfImageProcessor(imgP,x2,y2,samplSize/2);		
			vals[4]=VitimageUtils.valuesOfImageProcessor(imgP,x0,y1,samplSize/4);
			vals[5]=VitimageUtils.valuesOfImageProcessor(imgP,x1,y0,samplSize/4);
			vals[6]=VitimageUtils.valuesOfImageProcessor(imgP,x2,y1,samplSize/4);
			vals[7]=VitimageUtils.valuesOfImageProcessor(imgP,x2,y2,samplSize/4);
		}
		double []valsRet=VitimageUtils.statistics2D(vals);
		double[]tempStatsMeanVar=null;
		boolean objectIsOnTheBorder=false;
		double[]tempStats=new double[vals.length];
		for(int i=0;i<vals.length;i++) {
			tempStatsMeanVar=VitimageUtils.statistics1D(vals[i]);
			tempStats[i]=tempStatsMeanVar[0];
			if(true || tempStats[0]>valsRet[0]*1.3 || tempStats[0]<valsRet[0]*0.7 || tempStats[1]>valsRet[1]*3 || tempStats[0]<valsRet[0]*0.3) objectIsOnTheBorder=true;
		}
		if(objectIsOnTheBorder==false)return valsRet;

		
//		IJ.log("\nWarning in MRUtils at rice estimation. Divergence exists between stats of the background computed on the borders of image\nObject should be on a border, what can affect noise estimation");
//		IJ.log("\nGlobal mean value="+valsRet[0]+" and details mean values in borders="+TransformUtils.stringVectorN(tempStats, ""));
		double[]tempStatsCorrected=null;
		int incr=0;
		double[][]valsBis=null;
		if(config0firstCorner_config1Corners_config2Cross_config3Both==3) {
			tempStatsCorrected=new double[5];//Suppress the 3 maximum, that should be the border or corner where the object lies
			double[]tempStats2=doubleArraySort(tempStats);
			for(int i=0;i<5;i++)tempStatsCorrected[i]=tempStats2[i];
			double[]valsCorrected=new double[5];
			valsBis=new double[5][];
			for(int i=0;i<8;i++)if(tempStats[i]<=tempStatsCorrected[4])valsBis[incr++]=vals[i];
		}
		else{
			tempStatsCorrected=new double[3];//Suppress the 3 maximum, that should be the border or corner where the object lies
			double[]tempStats2=doubleArraySort(tempStats);
			for(int i=0;i<3;i++)tempStatsCorrected[i]=tempStats2[i];				
			double[]valsCorrected=new double[3];
			valsBis=new double[3][];
			for(int i=0;i<4;i++)if(tempStats[i]<=tempStatsCorrected[2])valsBis[incr++]=vals[i];
		}
		double []valsRetBis=VitimageUtils.statistics2D(valsBis);		
		IJ.log("MR slice import. Noise statistics before correction="+TransformUtils.stringVectorN(valsRet, "")+" / after correction="+TransformUtils.stringVectorN(valsRetBis, ""));
		return valsRetBis;
	}	
	
	
	

	/** Check data type in slice labels. */
	public static MRDataType getDataTypeOfThisMagneticResonanceSlice(ImagePlus img,int c,int z,int f) {
		String label=img.getStack().getSliceLabel(VitimageUtils.getCorrespondingSliceInHyperImage(img,c,z,f) );		
		int nbCat=0;
		MRDataType dataT=null;
		if(label.contains("T1SEQ")) {dataT=MRDataType.T1SEQ;nbCat++;}
		if(label.contains("T1T2SEQ")) {dataT=MRDataType.T1T2SEQ;nbCat++;}
		if(label.contains("T2SEQ") && !label.contains("T2SEQ") ) {dataT=MRDataType.T2SEQ;nbCat++;}
		if(label.contains("T1MAP")) {dataT=MRDataType.T1MAP;nbCat++;}
		if(label.contains("T2MAP")) {dataT=MRDataType.T2MAP;nbCat++;}
		if(label.contains("M0MAP")) {dataT=MRDataType.M0MAP;nbCat++;}
		if(nbCat>1) {
			IJ.showMessage("Critical fail in MRUtils : get "+nbCat+" categories instead of 1 \nwhen calling getDataTypeOfThisMagneticResonanceSlice(ImagePlus img,"+c+","+z+","+f);
		}
		if(nbCat<1)return null;
		return dataT;
	}
	
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

	
	

	/** Extract Tr, Te or Sigma values from the relaxation images */
	public static double[]getTrFrom3DRelaxationImageTab(ImagePlus[]imgTr){
		double[]tab=new double[imgTr.length];
		for(int i=0;i<imgTr.length;i++) {
			tab[i]=readValueOfSigmaTrTeInSliceLabel(imgTr[i],PARAM_TR,0,0,0);
		}
		return tab;
	}

	public static double[][]getTeFrom3DRelaxationImageTab(ImagePlus[]imgTr){
		double[][]tab=new double[imgTr[0].getNSlices()][imgTr.length];
		for(int i=0;i<imgTr.length;i++) {
			for(int z=0;z<imgTr[0].getNSlices();z++) {
				tab[z][i]=readValueOfSigmaTrTeInSliceLabel(imgTr[i],PARAM_TE,0,z,0);
			}
		}
		return tab;
	}
	
	public static double readTeInSliceLabel(ImagePlus img,int c, int z, int f) {
		return readValueOfSigmaTrTeInSliceLabel(img, PARAM_TE, c, z, f);
	}

	public static double readTrInSliceLabel(ImagePlus img,int c, int z, int f) {
		return readValueOfSigmaTrTeInSliceLabel(img, PARAM_TR, c, z, f);
	}

	public static double readSigmaInSliceLabel(ImagePlus img,int c, int z, int f) {
		return readValueOfSigmaTrTeInSliceLabel(img, PARAM_SIGMA, c, z, f);
	}

	
	
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


	
	
	
	
	/** Extract relaxation values for a voxel, or a voxel and its neighbours, get the corresponding coordinates in the images */

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	public static double []getProportionalTimes(double valMin, double valMax,double step){
		double[]tab=new double[(int)(Math.ceil((double)0.00001+((valMax-valMin)/step)))];
		for (int i=0;i<tab.length;i++){
			tab[i]=valMin+(double)((i*step));
		}
		return tab;
	}
	
	
	
	
	public static double[]getDataForVoxel(ImagePlus imgIn,int xCor,int yCor,int curT,int totalT,int curZ,int totalZ,int nEchoesExpected,int crossWidth,int crossThick,boolean gaussianWeighting){
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
	
	public static double[][]getFullDataForVoxel(ImagePlus imgIn,int xCor,int yCor,int curT,int totalT,int curZ,int totalZ,int nEchoesExpected,int crossWidth,int crossThick,boolean gaussianWeighting){
		if(gaussianWeighting) {IJ.showMessage("Not applicable : gaussian weighting. Abort");gaussianWeighting=false;}
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

	

	
	/** Compute a simple T1 or T2 fit, using Simplex or Levenberg-Marquardt, and estimate a MONO or BI exponential, corrupted by rice Noise.*/ 
	public static double[]makeFitDouble(final double[]tabTrTimesTemp, final double[]tabTeTimesTemp, double[]tabData,int fitType,int algType,int nbIter,double sigma){
		double[]estimatedParams;
		double[]tabTrTimes=new double[tabTrTimesTemp.length];
		for(int i=0;i<tabTrTimes.length;i++)tabTrTimes[i]=tabTrTimesTemp[i];
		double[]tabTeTimes=new double[tabTeTimesTemp.length];
		for(int i=0;i<tabTeTimes.length;i++)tabTeTimes[i]=tabTeTimesTemp[i];
		SimplexDualCurveFitter simpfitter=new SimplexDualCurveFitter(tabTrTimes,tabTeTimes,tabData,fitType,sigma);
	 	simpfitter.config(nbIter);
	 	simpfitter.doFit();
		estimatedParams=simpfitter.getParams();
		return estimatedParams;
	}	
	
	
	
	/** Compute a simple T1 or T2 fit, using Simplex or Levenberg-Marquardt, and estimate a MONO or BI exponential, corrupted by rice Noise.*/ 
	public static double[]makeFitSimple(double[]tabTimesTemp, double[]tabData,int fitType,int algType,int nbIter,double sigma){
		double[]estimatedParams;
		double[]tabTimes=new double[tabTimesTemp.length];
		for(int i=0;i<tabTimes.length;i++)tabTimes[i]=tabTimesTemp[i];
		if(algType==LM){
			LMCurveFitterNoBias lmfitter=new LMCurveFitterNoBias(tabTimes,tabData,fitType,sigma,false);
		 	lmfitter.configLMA(LMCurveFitterNoBias.lambda,LMCurveFitterNoBias.minDeltaChi2,nbIter);
		 	lmfitter.doFit();
			estimatedParams=lmfitter.getParams();
		}
		else {
			SimplexCurveFitterNoBias simpfitter=new SimplexCurveFitterNoBias(tabTimes,tabData,fitType,sigma);
		 	simpfitter.config(nbIter);
		 	simpfitter.doFit();
			estimatedParams=simpfitter.getParams();
		}
		return estimatedParams;
	}	
	

	/** Same as above, but compute also the distribution of estimated parameters when applying randomized modifications of the relaxation data */	
	public static double[][]makeFit(double[]tabTimes, double[]tabData,int fitType,int algType,int nbIter,double sigma,int nRepetMonteCarlo,int nAveragePts,RiceEstimator riceEstimator){
		boolean debugLM=false;
		Random rand=new Random();
		int nParams=(fitType==MRUtils.T2_RELAX_RICE || fitType==MRUtils.T1_RECOVERY_RICE ) ? 2 : 4;
		double[]tmpSimulatedData=new double[tabData.length];
		double[][]resultsMonteCarlo=new double[nParams][nRepetMonteCarlo];
		double[]estimatedParams;
		double[]simEstimatedParams;
		double[]estimatedSigmas=new double[nParams];
		double[]estimatedMeans=new double[nParams];
		double sigmaZero;
		if(algType==LM){
			LMCurveFitterNoBias lmfitter=new LMCurveFitterNoBias(tabTimes,tabData,fitType,sigma,debugLM);
		 	lmfitter.configLMA(LMCurveFitterNoBias.lambda,LMCurveFitterNoBias.minDeltaChi2,nbIter);
		 	lmfitter.doFit();
			estimatedParams=lmfitter.getParams();
			for(int n=0;n<nRepetMonteCarlo;n++) {
				for(int dat=0;dat<tabData.length;dat++) {
					sigmaZero=riceEstimator.estimateSigma(sigma,tabData[dat])/Math.sqrt(nAveragePts); 
					tmpSimulatedData[dat]=tabData[dat]+rand.nextGaussian()*sigmaZero;
				}
				lmfitter=new LMCurveFitterNoBias(tabTimes,tmpSimulatedData,fitType,sigma,debugLM);
				lmfitter.doFit();
				simEstimatedParams=lmfitter.getParams();
				for(int p=0;p<nParams;p++)resultsMonteCarlo[p][n]=simEstimatedParams[p];
			}
		}
		else {
			SimplexCurveFitterNoBias simpfitter=new SimplexCurveFitterNoBias(tabTimes,tabData,fitType,sigma);
		 	simpfitter.config(nbIter);
		 	simpfitter.doFit();
			estimatedParams=simpfitter.getParams();
			for(int n=0;n<nRepetMonteCarlo;n++) {
				for(int dat=0;dat<tabData.length;dat++) {
					sigmaZero=riceEstimator.estimateSigma(sigma,tabData[dat])/Math.sqrt(nAveragePts); 
					tmpSimulatedData[dat]=tabData[dat]+rand.nextGaussian()*sigmaZero;
				}
				simpfitter=new SimplexCurveFitterNoBias(tabTimes,tmpSimulatedData,fitType,sigma);
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


	
	
	public static double[] fittenRelaxationCurveT1T2(double[]tr,double []te, double []param,double sigma,int fitType){
		double[]tabRes=new double[tr.length];
		for(int indT=0;indT<tr.length;indT++){
			switch(fitType){
				case MRUtils.T1T2_MONO_RICE: tabRes[indT]=(double)(RiceEstimator.besFunkCost(param[0]* (double)Math.exp(-(te[indT] / param[2])) * (1-Math.exp(-(tr[indT] / param[1]))),sigma));break;
				case MRUtils.T1T2_MULTI_RICE: tabRes[indT]=(double)(RiceEstimator.besFunkCost(
						(1-Math.exp(-(tr[indT] / param[2]))) *
						(    param[0]* (double)Math.exp(-(te[indT] / param[3]))  +
						     param[1]* (double)Math.exp(-(te[indT] / param[4]))    )
					,sigma) );break;
				case MRUtils.T1T2_MULTIMULTI_RICE: tabRes[indT]=(double)(RiceEstimator.besFunkCost(
						param[0]* (double)Math.exp(-(te[indT] / param[4])) * (1-Math.exp(-(tr[indT] / param[2])))  +
						param[1]* (double)Math.exp(-(te[indT] / param[5])) * (1-Math.exp(-(tr[indT] / param[3])))
					,sigma) ) ;break;
				case MRUtils.T1T2_BIONANO: tabRes[indT]=  ( (tr[indT]>=10000) && (te[indT]>0) ) ? 
						(double)(RiceEstimator.besFunkCost(param[2]* (double)Math.exp(-(te[indT] / param[3])) * (1-Math.exp(-(tr[indT] / param[1]))),sigma)) : 
						(double)(RiceEstimator.besFunkCost(param[0]* (double)Math.exp(-(te[indT] / param[3])) * (1-Math.exp(-(tr[indT] / param[1]))),sigma));
				break;
			}
		}
		return tabRes;
	}

	
	/** Compute fitting results : estimated curve corresponding to the parameters estimated, accuracy of fit, khi2, and corresponding pvalue */	
	public static double[] fittenRelaxationCurve(double[]tEchos,double []param,double sigma,int fitType){
		double[]tab=new double[tEchos.length];
		double techo;
		for(int indT=0;indT<tEchos.length;indT++){
			techo=tEchos[indT];
			switch(fitType){
				case MRUtils.T2_RELAX_RICE: tab[indT]=(double)(RiceEstimator.besFunkCost(param[0]* (double)Math.exp(-(techo / param[1])),sigma));break;
				case MRUtils.T2_MULTICOMP_RICE: tab[indT]=(double)(RiceEstimator.besFunkCost(param[0]* (double)Math.exp(-(techo / param[1]))+param[2]* (double)Math.exp(-(techo / param[3])),sigma));break;
				case MRUtils.T1_RECOVERY_RICE: tab[indT]=(double)(RiceEstimator.besFunkCost(param[0]* (double)(1-Math.exp(-(techo / param[1]))) , sigma));break;
				case MRUtils.T1_MULTICOMP_RICE: tab[indT]=(double)(RiceEstimator.besFunkCost(param[0]* (double)(1-Math.exp(-(techo / param[1])) ) +param[2]* (double)(1-Math.exp(-(techo / param[3])) ) , sigma));break;
			}
		}
		return tab;
	}

	public static double[] fittingAccuracies(double[] tabData,double[] tabTimes,double sigma,double[] estimatedParams,int fitType,boolean debug,RiceEstimator riceEstimator,int nbPoints){		
		double []realP=fittenRelaxationCurve(tabTimes,estimatedParams,sigma,fitType);
		double khi2=0;
		double diff=0;
		double sigmaEst=0;
		for(int indT=0;indT<tabData.length;indT++) {
			sigmaEst=riceEstimator.estimateSigma(sigma, realP[indT])/Math.sqrt(nbPoints);
			diff=(realP[indT]-tabData[indT]);
			khi2+=diff*diff/(sigmaEst*sigmaEst);
		}
		int N_PARAMS=( (fitType==MRUtils.T1_RECOVERY_RICE || fitType==MRUtils.T2_RELAX_RICE) ? 2 : 4);
		double pVal=getPvalue(khi2,tabTimes.length-N_PARAMS);
		double khi2Ret=khi2/(tabTimes.length-N_PARAMS);
 		return new double[] {khi2Ret,pVal*100};
	}	
	

	public static double[] fittingAccuraciesT1T2(double[] tabData,double[] tabTrTimes,double[] tabTeTimes,double sigma,double[] estimatedParams,int fitType,boolean debug,RiceEstimator riceEstimator,int nbPoints){		
		int N_PARAMS=( (fitType==MRUtils.T1T2_MONO_RICE) ? 3 : (fitType==MRUtils.T1T2_MULTI_RICE ?  5 : 6 ));
		if(fitType==MRUtils.T1T2_BIONANO)N_PARAMS=4;
		double []realP=fittenRelaxationCurveT1T2(tabTrTimes,tabTeTimes,estimatedParams,sigma,fitType);
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
	

	
	
	
	
	public static double[] fittingAccuraciesFullData(double[][] tabData,double[] tabTimes,double sigma,double[] estimatedParams,int fitType,boolean debug,RiceEstimator riceEstimator){		
		double []realP=fittenRelaxationCurve(tabTimes,estimatedParams,sigma,fitType);
		double khi2=0;
		double diff=0;
		double sigmaEst=0;
		for(int indT=0;indT<tabData.length;indT++) {
			sigmaEst=riceEstimator.estimateSigma(sigma, realP[indT]);
			for(int indP=0;indP<tabData[indT].length;indP++) {
				diff=(realP[indT]-tabData[indT][indP]);
				khi2+=diff*diff/(sigmaEst*sigmaEst);
			}
		}
		int N_PARAMS=( (fitType==MRUtils.T1_RECOVERY_RICE || fitType==MRUtils.T2_RELAX_RICE) ? 2  : 4);
		double pVal=getPvalue(khi2,tabTimes.length*tabData[0].length-N_PARAMS);
		double khi2Ret=khi2/(tabTimes.length*tabData[0].length-N_PARAMS);
 		return new double[] {khi2Ret,pVal*100};
	}	

	public static double getPvalue(double khi2,int freedomDegrees){
		ChiSquaredDistribution x2 = new ChiSquaredDistribution( freedomDegrees );
		try {
			return(x2.cumulativeProbability(khi2));
		} catch (Exception e) {
			e.printStackTrace();
		}
		return -1;
	}

	
	public static double convertM0InPercentage(double m0) {
		return (m0/maxM0ForNormalization);
	}
	
	public static ImagePlus[] normalizeBeforeComputation(ImagePlus []imgsEchoes,ImagePlus imgMask) {
//		System.out.println("Normalisation d une entree\n");
		//Entrees : imgT1T2[][];
		//Pour chaque i, sommer tous les echos disponibles en une seule image float
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
		double globalMeanOverAll;
		double globalMeanOverCap;
		double globalMeanOverRest;
		int radiusAroundCapillary=(int)Math.round(HyperMRIT1T2.bionanoCapillaryRadius*1.5/VitimageUtils.getVoxelSizes(imgSum)[0]);
		//Pour chaque Z
		for(int z=0;z<Z;z++) {
			ImagePlus sliceTemp=new Duplicator().run(imgSum,1,1,z+1,z+1,1,1);
			ImagePlus sliceMask=new Duplicator().run(imgMask,1,1,z+1,z+1,1,1);
			//Pour cette image, localiser le capillaire, sur chaque Z
			coords[z]=HyperMRIT1T2.findCapillaryCenterInSlice(sliceTemp, HyperMRIT1T2.bionanoCapillaryRadius);
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
		globalMeanOverAll=VitimageUtils.statistics1D(meanOverAll)[0];
		globalMeanOverRest=VitimageUtils.statistics1D(meanOverRest)[0];
		globalMeanOverCap=VitimageUtils.statistics1D(meanOverCap)[0];

		int radiusSq=radiusAroundCapillary*radiusAroundCapillary;
		double distSq;
		//Pour chaque Z
		for(int im=0;im<tabRet.length;im++) {
			IJ.run(tabRet[im],"32-bit","");
			for(int z=0;z<Z;z++) {
				float[]tabData=(float[])tabRet[im].getStack().getProcessor(z+1).getPixels();
				double factorMultCap=globalMeanOverCap/meanOverCap[z];
				double factorMultRest=globalMeanOverRest/meanOverRest[z];
//				System.out.println("Moyennes detectees à z="+z+" : All="+meanOverAll[z]+"/"+VitimageUtils.dou(globalMeanOverAll)+" , Cap="+meanOverCap[z]+"/"+VitimageUtils.dou(globalMeanOverCap)+" , Rest="+meanOverRest[z]+"/"+VitimageUtils.dou(globalMeanOverRest)+
//						"   -> factor Cap="+factorMultCap+"   -> factor Rest="+factorMultRest);
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
	
	

	public static final int ERROR_VALUE= 0;
	public static final int SIMPLEX = 1;
    public static final int LM=2;
    public static final int TWOPOINTS=3; 	   

    public static final int MSEC=11;
    public static final int SEC=12;

    public static final int T2_RELAX_RICE =21; //offset 3
    public static final int T2_MULTICOMP_RICE=22;
	public static final int T1_RECOVERY_RICE = 23; //offset 3
	public static final int T1_MULTICOMP_RICE = 24; //offset 3
	public static final int T1T2_MONO_RICE = 25; //offset 3
	public static final int T1T2_MULTI_RICE = 26; //offset 3
	public static final int T1T2_MULTIMULTI_RICE = 27; //offset 3
	public static final int T1T2_BIONANO = 29; //offset 3

	public final static String[] timeunits={"ms", "s"};
    public final static int[] timeitems={MSEC, SEC};
    public final static String[] fititems2={"Simplex","Levenberg-Marquardt"};
    public final static int[] constitems2={SIMPLEX,LM};

	public static final int PARAM_SIGMA=31;
	public static final int PARAM_TR=32;
	public static final int PARAM_TE=33;
    public static final double DELTA_TE_BOUTURE_TRICK=33.453;
    
	public static final int maxDisplayedT2=150;
	public static final int maxDisplayedT1=3200;
	public static final int maxDisplayedM0=8000;
	public static final int maxM0ForNormalization=7000;
	public static final int THRESHOLD_RATIO_BETWEEN_MAX_ECHO_AND_SIGMA_RICE=7;
  
    
    public static int N_ITER_T1T2=30;
}
