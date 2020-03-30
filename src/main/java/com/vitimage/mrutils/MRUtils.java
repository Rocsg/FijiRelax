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
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
public class MRUtils  {
//	static float stdValMaxIRM=50000;
	//	public static float stdSigmaIRM=159;
	
	public MRUtils() {
		
	}


	/** Top level entry points : compute maps from the successive images of the echo or repetition times*/
	public static ImagePlus[] computeT2Map(final ImagePlus[]imgTe,double sigmaSmoothing) {
		if(imgTe[0].getType() != ImagePlus.GRAY16) {VitiDialogs.notYet("computeT2map != 16 in computeMaps in MRI_T1_Seq : "+imgTe[0].getType());return null;}
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
		final double []listTeForThreads=getTeFrom3DRelaxationImageTab(imgTe);
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


	
	public static double[] doubleArraySort(double[]tab) {
		double[]tabRet=new double[tab.length];
		ArrayList<Double> list=new ArrayList<Double>();
		for(double dou : tab)list.add(dou);
		Collections.sort(list);
		for(int i=0;i<list.size();i++)tabRet[i]=list.get(i);
		return tabRet;
	}


	
	
	
	/** Estimate rice noise in each slice from the background values*/
	public static void computeRiceSigmaOfEachSliceAndWriteItInTheLabels(ImagePlus img,boolean keepOldString,String optionalAdditionalString) {
		int nbC=img.getNChannels();
		int nbZ=img.getNSlices();
		int nbF=img.getNFrames();
		ImageProcessor imgP;
		for(int nc=0;nc<nbC;nc++) {
			for(int nz=0;nz<nbZ;nz++) {
				for(int nf=0;nf<nbF;nf++) {
					int sli=nbC*nbZ*nf+nbC*nz+nc+1;
					imgP=img.getStack().getProcessor(sli);
					double[]stats=getBackgroundStatsFromProcessor(imgP,3);
					double sigmaRice=RiceEstimator.computeRiceSigmaFromBackgroundValuesStatic(stats[0],stats[1]);
					img.getStack().setSliceLabel((keepOldString ? img.getStack().getSliceLabel(sli) : "") + "_"+optionalAdditionalString+"_SIGMARICE="+VitimageUtils.dou(sigmaRice), sli);
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
		if(label.contains("T2SEQ")) {dataT=MRDataType.T2SEQ;nbCat++;}
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
				if(prm0Sigma_1Tr_2Te==PARAM_TR && prm[0].equals("TR"))return Integer.parseInt(prm[1]);
				if(prm0Sigma_1Tr_2Te==PARAM_TE && prm[0].equals("TE"))return Integer.parseInt(prm[1]);
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

	public static double[]getTeFrom3DRelaxationImageTab(ImagePlus[]imgTr){
		double[]tab=new double[imgTr.length];
		for(int i=0;i<imgTr.length;i++) {
			tab[i]=readValueOfSigmaTrTeInSliceLabel(imgTr[i],PARAM_TE,0,0,0);
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
	public static double[]makeFitSimple(double[]tabTimesTemp, double[]tabData,int fitType,int algType,int nbIter,double sigma){
		double[]estimatedParams;
		double[]tabTimes=new double[tabTimesTemp.length];
		for(int i=0;i<tabTimes.length;i++)tabTimes[i]=tabTimesTemp[i];
		if(algType==MRUtils.TWOPOINTS){
			TwoPointsCurveFitterNoBias twopointsfitter=new TwoPointsCurveFitterNoBias(tabTimes, tabData,fitType, sigma);
			twopointsfitter.doFit();
			estimatedParams=twopointsfitter.getParams();
		}
		else if(algType==LM){
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
		int nParams=(fitType==MRUtils.T2_RELAX_RICE || fitType==MRUtils.T1_RECOVERY_RICE || fitType==MRUtils.TWOPOINTS) ? 2 : (fitType==MRUtils.TRICOMP_RICE ) ? 6 : 4;
		double[]tmpSimulatedData=new double[tabData.length];
		double[][]resultsMonteCarlo=new double[nParams][nRepetMonteCarlo];
		double[]estimatedParams;
		double[]simEstimatedParams;
		double[]estimatedSigmas=new double[nParams];
		double[]estimatedMeans=new double[nParams];
		double sigmaZero;
		if(algType==MRUtils.TWOPOINTS){
			TwoPointsCurveFitterNoBias twopointsfitter=new TwoPointsCurveFitterNoBias(tabTimes, tabData,fitType, sigma);
			twopointsfitter.doFit();
			estimatedParams=twopointsfitter.getParams();
			for(int n=0;n<nRepetMonteCarlo;n++) {
				for(int dat=0;dat<tabData.length;dat++) {
					sigmaZero=riceEstimator.estimateSigma(sigma,tabData[dat])/Math.sqrt(nAveragePts); 
					tmpSimulatedData[dat]=tabData[dat]+rand.nextGaussian()*sigmaZero;
				}
				twopointsfitter=new TwoPointsCurveFitterNoBias(tabTimes, tmpSimulatedData,fitType, sigma);
				twopointsfitter.doFit();
				simEstimatedParams=twopointsfitter.getParams();
				for(int p=0;p<nParams;p++)resultsMonteCarlo[p][n]=simEstimatedParams[p];
			}
		}
		else if(algType==LM){
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
		for(int p=0;p<nParams;p++) {estimatedSigmas[p]=VitimageUtils.statistics1D(resultsMonteCarlo[p])[1];estimatedMeans[p]=VitimageUtils.statistics1D(resultsMonteCarlo[p])[0];}
		double [][]ret=new double[2][nParams];
		ret[0]=estimatedParams;
		ret[1]=estimatedSigmas;
		return ret;
	}


	
	
	
	
	/** Compute fitting results : estimated curve corresponding to the parameters estimated, accuracy of fit, khi2, and corresponding pvalue */	
	public static double[] fittenRelaxationCurve(double[]tEchos,double []param,double sigma,int fitType){
		double[]tab=new double[tEchos.length];
		double techo;
		for(int indT=0;indT<tEchos.length;indT++){
			techo=tEchos[indT];
			switch(fitType){
				case MRUtils.T2_RELAX: tab[indT]=(double)(param[0]* (double)Math.exp(-(techo / param[1])));break;
				case MRUtils.T2_RELAX_SIGMA: tab[indT]=(double)(param[0]* (double)Math.exp(-(techo / param[1])));break;
			 	case MRUtils.T2_RELAX_BIAS: tab[indT]=(double)(param[0]* (double)Math.exp(-(techo / param[1]))+param[2]);break;
				case MRUtils.T2_RELAX_RICE: tab[indT]=(double)(RiceEstimator.besFunkCost(param[0]* (double)Math.exp(-(techo / param[1])),sigma));break;

				case MRUtils.T1_RECOVERY: tab[indT]=(double)(param[0]* (double)(1-Math.exp(-(techo / param[1]))));break;
				case MRUtils.T1_RECOVERY_RICE: tab[indT]=(double)(RiceEstimator.besFunkCost(param[0]* (double)(1-Math.exp(-(techo / param[1]))) , sigma));break;
			 	case MRUtils.T1_RECOVERY_RICE_NORMALIZED: tab[indT]=(double)(RiceEstimator.besFunkCost( (double)(1-Math.exp(-(techo / param[0]))) , sigma));break;

				case MRUtils.MULTICOMP: tab[indT]=(double)(param[0]* (double)Math.exp(-(techo / param[1]))+param[2]* (double)Math.exp(-(techo / param[3])));break;
				case MRUtils.MULTICOMP_BIAS: tab[indT]=(double)(param[0]* (double)Math.exp(-(techo / param[1]))+param[2]* (double)Math.exp(-(techo / param[3]))+param[4]);break;
				case MRUtils.MULTICOMP_RICE: tab[indT]=(double)(RiceEstimator.besFunkCost(param[0]* (double)Math.exp(-(techo / param[1]))+param[2]* (double)Math.exp(-(techo / param[3])),sigma));break;
				case MRUtils.TRICOMP_RICE: tab[indT]=(double)(RiceEstimator.besFunkCost(param[0]* (double)Math.exp(-(techo / param[1]))+param[2]* (double)Math.exp(-(techo / param[3])+param[4]* (double)Math.exp(-(techo / param[5]))),sigma) );break;
			}
		}
		return tab;
	}
/*
	public static double fittingAccuracy(double[] tabData,double[] tabTimes,double sigma,double[] estimatedParams,int fitType,boolean debug){		
		double []realP=fittenRelaxationCurve(tabTimes,estimatedParams,sigma,fitType);
		double cumulator=0;
		double meanBg=RiceEstimator.computeSigmaAndMeanBgFromRiceSigmaStatic(sigma)[0];
		for(int indT=0;indT<tabTimes.length;indT++)cumulator+=(double)((realP[indT]-tabData[indT]) * (realP[indT]-tabData[indT]));

		//Case T2 relaxation
		double valRef=0;
		for(int t=0;t<tabData.length;t++)valRef+=(1.0/tabData.length)*tabData[t];
		cumulator=Math.sqrt(cumulator/tabTimes.length);
		cumulator=cumulator*100.0/((Math.abs(valRef-meanBg)));
		if(cumulator>99)cumulator=99;
		return cumulator;
	}	
*/
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
		int N_PARAMS=( (fitType==MRUtils.T1_RECOVERY_RICE || fitType==MRUtils.T2_RELAX_RICE) ? 2 : (fitType==MRUtils.TRICOMP_RICE) ? 6 : 4);
		double pVal=getPvalue(khi2,tabTimes.length-N_PARAMS);
		double khi2Ret=khi2/(tabTimes.length-N_PARAMS);
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
		int N_PARAMS=( (fitType==MRUtils.T1_RECOVERY_RICE || fitType==MRUtils.T2_RELAX_RICE) ? 2 : (fitType==MRUtils.TRICOMP_RICE) ? 6 : 4);
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

	

	
	
	
	public static final int RICE=100;
	public static final int BIAS=1000;
	public static final int SIGMA=10000;
	public static final int MULTI=100000;
    public static final int ALL_AVAILABLE_FIT=1000000;
	public static final int TRI=100000000;

	public static final int STRAIGHT_LINE=0,EXPONENTIAL=STRAIGHT_LINE+4,EXP_RECOVERY=STRAIGHT_LINE+13;
	public static final int T2_RELAX = EXP_RECOVERY +3; //offset 3
    public static final int T2_RELAX_BIAS = T2_RELAX+BIAS; //offset 3
    public static final int T2_RELAX_SIGMA = T2_RELAX+SIGMA; //offset 3
    public static final int T2_RELAX_RICE = T2_RELAX+RICE; //offset 3
	public static final int MULTICOMP=T2_RELAX+MULTI;
    public static final int MULTICOMP_BIAS=T2_RELAX+MULTI+BIAS;
    public static final int MULTICOMP_SIGMA=T2_RELAX+MULTI+SIGMA;
    public static final int MULTICOMP_RICE=T2_RELAX+MULTI+RICE;
    public static final int TRICOMP_RICE=T2_RELAX+TRI+RICE;
 	public static final int T1_RECOVERY = 500; //offset 3
	public static final int T1_RECOVERY_RICE = 501; //offset 3
	public static final int T1_RECOVERY_RICE_NORMALIZED = 504; //offset 3
	public static final int GAUSSIAN=17;
	public static final int ERROR_VALUE= 0;
	public static final int SIMPLEX = 1;
    public static final int LM=2;
    public static final int TWOPOINTS=3; 	   
    public static final int MSEC=3;
    public static final int SEC=0;

    public final static String[] timeunits={"ms", "s"};
    public final static int[] timeitems={MSEC, SEC};
    public final static String[] fititems2={"Simplex","Levenberg-Marquardt"};
    public final static int[] constitems2={SIMPLEX,LM};

	public static final int PARAM_SIGMA=0;
	public static final int PARAM_TR=1;
	public static final int PARAM_TE=2;
    
    
    
  
    
    
}
