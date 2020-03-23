package com.vitimage.mrutils;

import java.io.File;
import com.vitimage.common.VitiDialogs;
import com.vitimage.common.VitimageUtils;
import com.vitimage.fijiyama.RegistrationAction;
import com.vitimage.fijiyama.RegistrationManager;
import com.vitimage.registration.BlockMatchingRegistration;
import com.vitimage.registration.ItkTransform;
import com.vitimage.registration.MetricType;
import com.vitimage.registration.Transform3DType;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.plugin.Concatenator;
import ij.plugin.FolderOpener;
import ij.plugin.HyperStackConverter;
import ij.plugin.frame.PlugInFrame;
import math3d.Point3d;

public class Bionano_T1T2M0_Importer extends PlugInFrame{
	boolean debug=false;
	
	boolean debugCompute=false;
	private final int maxDisplayedT2=125;
	private final int maxDisplayedT1=3000;
	private final int maxDisplayedM0=7000;
	private static final long serialVersionUID = 1L;
	private String inputDir="";
	private String outputDir="";
	private int test;
	private double[]voxs;
	private int[]dims;
	private int []valT1Trs;
	private int []valT2Tes;
	private boolean haveT2Seq;
	ImagePlus []imgT2Seq;
	ImagePlus []imgT1Seq;
	ImagePlus []imgT2SeqTest;
	ImagePlus []imgT1SeqTest;
	ImagePlus M0map;
	ImagePlus T1map;
	ImagePlus T2map;
	private boolean haveT1Seq;
	private static final int TrValueUsedForT2Seq=10000;
	
	public Bionano_T1T2M0_Importer() {
		super("");
	}

	public static void main(String[] args) {
		@SuppressWarnings("unused")
		ImageJ ij=new ImageJ();
		ItkTransform tr=new ItkTransform();
		IJ.log("Running T1-T2-M0 import");
		Bionano_T1T2M0_Importer bn=new Bionano_T1T2M0_Importer();
		bn.debugCompute=false;
		bn.run(bn.debug ? new Object[]{"/home/fernandr/Bureau/Traitements/Sorgho/Test SSM1/Source/SSM1_F1_0923_SL_RAW","/home/fernandr/Bureau/Traitements/Sorgho/Test SSM1/Out_F1"} : null);
	}

	public void run(String arg) {
		IJ.log("Starting MRI Importer");
		run(new Object[] {});
	}
	
	public void run(Object[]obj) {
		//Ask for the directory input, and the directory output
		IJ.log("Opening input and output dirs");
		if(obj != null && obj.length==2) {
			inputDir=(String)obj[0];
			outputDir=(String)obj[1];
		}
		else{
			inputDir=VitiDialogs.chooseDirectoryUI("Localize inputDir (SL_RAW)", "Select this directory");
			outputDir=VitiDialogs.chooseDirectoryUI("Select outputDir", "Select this directory");
		}

		//Read the voxel sizes and the image sizes		
		IJ.log("Gathering image informations");
		Object []objs=getFirstImageOfT1T2Serie(inputDir);
		ImagePlus testDcm=(ImagePlus)objs[0];
		int nSlices=(int)objs[1];
		this.voxs=VitimageUtils.getVoxelSizes(testDcm);
		dims=VitimageUtils.getDimensions(testDcm);
		dims[2]=nSlices;


		//Look for a T2 sequence, take the successives Te
		IJ.log("Gathering T2 sequence");
		if(new File(inputDir,"TR010000").exists()) {
			haveT2Seq=true;
			readT2Echoes();
			displayT2Sequence();
		}

		//Look for a T1 sequence, take the successives Tr
		IJ.log("Gathering T1 sequence");
		if(new File(inputDir).list().length>2) {
			haveT1Seq=true;
			readT1Reps();
			displayT1Sequence();
		}

		//Register T1 seq if needed
		IJ.log("Registering T1 sequence");
		if(haveT1Seq && (debug || VitiDialogs.getYesNoUI("Register T1 repetition images ?","Register T1 repetition images ??"))){
			registerT1Reps();
		}
		displayT1Sequence();
			

		//Show T1 estimation and T1 M0 for multiple sigma, then compute T1 and M0
		IJ.log("Computing T1 map");
		computeT1(new Object[] {1.0});
		T1map.show();
		T1map.setDisplayRange(0, maxDisplayedT1);
		T1map.setTitle("T1 map");


		
		//Show T2 estimation and T2 M0 for multiple sigma, then compute T2 and M0
		IJ.log("Computing T2 map");
		computeT2(new Object[] {1.0});
		T2map.show();
		T2map.setDisplayRange(0, maxDisplayedT2);
		M0map.show();
		M0map.setDisplayRange(0, maxDisplayedM0);
		T2map.setTitle("T2 map");
		M0map.setTitle("M0 map");
		
		
		//Stack M0, T1, T2, T1 reps, T2 echoes and save the result
		IJ.log("Stacking results");
		stackSaveAndShowResults();
		IJ.log("Finished !");
	}
	
	
	public Object[] getFirstImageOfT1T2Serie(String inputDir) {
		String imgPath=inputDir;
		String[]subFile=null;
		while(!imgPath.substring(imgPath.length()-4, imgPath.length()).equals(".dcm")){
			IJ.log("Not yet an img :"+imgPath);
			subFile=VitimageUtils.stringArraySort(new File(imgPath).list());
			imgPath=new File(imgPath,subFile[0]).getAbsolutePath();
		}
		return new Object[] {IJ.openImage(imgPath),subFile.length};
	}
	
	
	public void displayT2Sequence() {
		ImagePlus[]imgs=new ImagePlus[imgT2Seq.length];int nF=imgs.length;int nZ=imgT2Seq[0].getNSlices();
		for(int i=0;i<imgT2Seq.length;i++) imgs[i]=VitimageUtils.imageCopy(imgT2Seq[i]);	

		ImagePlus test=HyperStackConverter.toHyperStack(Concatenator.run(imgs), 1,nZ,nF,"xyztc","Fire");
		test.setSlice(dims[2]/2);
		test.setTitle("Stacked T2 echoes in original position");
		test.resetDisplayRange();
		test.show();
	}
	
	public void displayT1Sequence() {
		ImagePlus[]imgs=new ImagePlus[imgT1Seq.length];int nF=imgs.length;int nZ=imgT1Seq[0].getNSlices();
		for(int i=0;i<imgT1Seq.length;i++) imgs[i]=VitimageUtils.imageCopy(imgT1Seq[i]);

		ImagePlus test=HyperStackConverter.toHyperStack(Concatenator.run(imgs), 1,nZ,nF,"xyztc","Fire");
		test.setTitle("Stacked T1 rep times in original position");
		test.setSlice(test.getStackSize()-1-dims[2]/2);
		test.resetDisplayRange();
		test.show();
	}

	
	
	
	
	public void readT2Echoes() {
		String[]listFiles=VitimageUtils.stringArraySort(new File(inputDir,"TR010000").list());
		this.valT2Tes=new int[listFiles.length];
		this.imgT2Seq=new ImagePlus[listFiles.length];
		this.imgT2SeqTest=new ImagePlus[listFiles.length];
		for(int i=0;i<this.valT2Tes.length;i++) {
			String str=new File(new File(inputDir,"TR0"+TrValueUsedForT2Seq),listFiles[i]).getAbsolutePath();//TODO : genericity among TRvalue : TrValueUsedForT2Seq is fixed
			imgT2Seq[i] = FolderOpener.open(str, "");
			this.valT2Tes[i]=Integer.parseInt(listFiles[i].replace("TE",""));
			MRUtils.computeRiceSigmaOfEachSliceAndWriteItInTheLabels(imgT2Seq[i],false,"T2SEQ_TR="+TrValueUsedForT2Seq+"_TE="+this.valT2Tes[i]);
			//imgT2SeqTest[i]=VitimageUtils.cropImageShort(imgT2Seq[i], (dims[0]*1)/3,(dims[1]*1)/3, 0,(dims[0])/3,(dims[1])/3,4);
			//if(debug)imgT2Seq[i]=VitimageUtils.imageCopy(imgT2SeqTest[i]);
		}
	}		
			
	public void readT1Reps() {
		String[]listFiles=VitimageUtils.stringArraySort(new File(inputDir).list());
		String[]strTes=VitimageUtils.stringArraySort(new File(inputDir,listFiles[0]).list());
		String strChosenTe=strTes[0];
		int teValue=Integer.parseInt(strChosenTe.replace("TE",""));
		this.valT1Trs=new int[listFiles.length-(haveT2Seq ? 1 : 0)];
		this.imgT1Seq=new ImagePlus[this.valT1Trs.length];
		this.imgT1SeqTest=new ImagePlus[this.valT1Trs.length];
		for(int i=0;i<this.valT1Trs.length;i++) {
			String str=new File(new File(inputDir,listFiles[i]),strChosenTe).getAbsolutePath();
			imgT1Seq[i] = FolderOpener.open(str, "");
			this.valT1Trs[i]=Integer.parseInt(listFiles[i].replace("TR",""));
			MRUtils.computeRiceSigmaOfEachSliceAndWriteItInTheLabels(imgT1Seq[i],false,"T1SEQ_TR="+this.valT1Trs[i]+"_TE="+teValue);
			//imgT1SeqTest[i]=VitimageUtils.cropImageShort(imgT1Seq[i],  (dims[0]*1)/3,(dims[1]*1)/3, 0,(dims[0])/3,(dims[1])/3, 4);
			//if(debug)imgT1Seq[i]=VitimageUtils.imageCopy(imgT1SeqTest[i]);
		}
	}		

	public void registerT1Reps(){		
		boolean crossFit=false;		
		int nT1=this.imgT1Seq.length;
		ItkTransform  []trs=new ItkTransform[nT1];
		trs[nT1-1]=new ItkTransform();
		ImagePlus imgMov,imgRef;
		ItkTransform tr=new ItkTransform();
		RegistrationAction regAct;
		BlockMatchingRegistration bmRegistration;

		//Registration of last repetition time of T1 seq (mov) onto the first echo time of T2 seq (ref)
		if(haveT1Seq && haveT2Seq) {
			ImagePlus comp=VitimageUtils.compositeOf(imgT2Seq[0],imgT1Seq[nT1-1],"red=T2 first echo , green = T1 last rep time");
			comp.show();
			if( VitiDialogs.getYesNoUI("Register T1 and T2 rigid ?", "Register T1 and T2 rigid ?")){
				comp.close();
				imgMov=VitimageUtils.imageCopy(imgT1Seq[nT1-1]);
				imgRef=VitimageUtils.imageCopy(imgT2Seq[0]);
				if(VitiDialogs.getYesNoUI("Manual registration first ?", "Manual registration first ?\nNeeded if T1 stack and T2 stack are largely disaligned")) {
					Point3d[][]pts=VitiDialogs.registrationPointsUI(5,imgRef,imgMov,true);
					tr=ItkTransform.estimateBestRigid3D(pts[1],pts[0]); 
				}
				else tr=new ItkTransform();
				imgRef.show();
				imgMov.show();
				regAct=new RegistrationAction();
				regAct.defineSettingsSimplyFromTwoImages(imgRef,imgMov);
				regAct.higherAcc=1;
				bmRegistration=BlockMatchingRegistration.setupBlockMatchingRegistration(imgRef,imgMov,regAct);
				bmRegistration.consoleOutputActivated=false;
				bmRegistration.timingMeasurement=true;
				bmRegistration.refRange=new double[] {imgRef.getDisplayRangeMin(),imgRef.getDisplayRangeMax()};
				bmRegistration.movRange=new double[] {imgMov.getDisplayRangeMin(),imgMov.getDisplayRangeMax()};
				bmRegistration.flagRange=true;
				bmRegistration.minBlockVariance=10;
				bmRegistration.displayRegistration=2;
				bmRegistration.displayR2=false;
				bmRegistration.levelMin=-1;
				bmRegistration.returnComposedTransformationIncludingTheInitialTransformationGiven=true;
				tr=bmRegistration.runBlockMatching(tr);
				bmRegistration.closeLastImages();

				//Ask for using the TR10000 to make the T1 fit. In this case, the registration between both should be really accurate, implying potentially
				ImagePlus temp=tr.transformImage( imgT2Seq[0], imgT1Seq[nT1-1], false);
				comp=VitimageUtils.compositeOf(imgT2Seq[0],temp,"red=T2 first echo , green = T1 last rep time");
				comp.show();
				if(VitiDialogs.getYesNoUI("Register T1 and T2 dense ? ", "Register T1 and T2 dense ? \n(to obtain a very high accuracy in positioning)")) {
					comp.close();
					regAct=new RegistrationAction();
					regAct.defineSettingsSimplyFromTwoImages(imgRef,imgMov);
					regAct.higherAcc=1;
					regAct.typeTrans=Transform3DType.DENSE;
					bmRegistration=BlockMatchingRegistration.setupBlockMatchingRegistration(imgRef,imgMov,regAct);
					bmRegistration.consoleOutputActivated=false;
					bmRegistration.timingMeasurement=true;
					bmRegistration.refRange=new double[] {imgRef.getDisplayRangeMin(),imgRef.getDisplayRangeMax()};
					bmRegistration.movRange=new double[] {imgMov.getDisplayRangeMin(),imgMov.getDisplayRangeMax()};
					bmRegistration.flagRange=true;
					bmRegistration.minBlockVariance=10;
					bmRegistration.displayRegistration=2;
					bmRegistration.displayR2=false;
					bmRegistration.levelMin=-1;
					bmRegistration.returnComposedTransformationIncludingTheInitialTransformationGiven=true;
					tr=bmRegistration.runBlockMatching(tr);
					bmRegistration.closeLastImages();
				}
			}
			if(comp.isVisible())comp.close();
			ImagePlus temp=tr.transformImage( imgT2Seq[0], imgT1Seq[nT1-1], false);
			comp=VitimageUtils.compositeOf(imgT2Seq[0],temp,"red=T2 first echo , green = T1 last rep time");
			comp.show();
			if(VitiDialogs.getYesNoUI("Use the first T2 echo for T1 seq computation ?", "Use the TR10000 of T2 seq as a T1 seq point")) crossFit=true;	
			comp.close();
		}
		

		
		
		//Registration of T1 successive images, if needed
		ImagePlus comp=VitimageUtils.compositeOf(imgT1Seq[0],imgT1Seq[imgT1Seq.length-1],"T1 seq (red=first rep time , green = last rep time)");
		comp.show();
		if(VitiDialogs.getYesNoUI("Register T1 images before T1 map estimation ? ", "Register T1 images before T1 map estimation ? (rigid reg.)")) {
			comp.close();
			imgRef=imgT1Seq[nT1-1];
			for(int iMov=0;iMov<nT1;iMov++) {
				if(iMov<nT1-1) {
					imgMov=imgT1Seq[iMov];
					regAct=new RegistrationAction();
					regAct.defineSettingsSimplyFromTwoImages(imgRef,imgMov);
					regAct.higherAcc=1;
					regAct.typeTrans=Transform3DType.EULER2D;
					bmRegistration=BlockMatchingRegistration.setupBlockMatchingRegistration(imgRef,imgMov,regAct);
					bmRegistration.consoleOutputActivated=false;
					bmRegistration.timingMeasurement=true;
					bmRegistration.refRange=new double[] {imgRef.getDisplayRangeMin(),imgRef.getDisplayRangeMax()};
					bmRegistration.movRange=new double[] {imgMov.getDisplayRangeMin(),imgMov.getDisplayRangeMax()};
					bmRegistration.flagRange=true;
					bmRegistration.minBlockVariance=10;
					bmRegistration.displayRegistration=2;
					bmRegistration.displayR2=false;
					bmRegistration.returnComposedTransformationIncludingTheInitialTransformationGiven=true;
					bmRegistration.levelMin=-1;
					trs[iMov]=bmRegistration.runBlockMatching(null);
					bmRegistration.closeLastImages();
				}
				else trs[iMov]=new ItkTransform();
				trs[iMov].addTransform(tr);
				IJ.log("Registration of T1 times finished !");
				String []labels=new String[imgT1Seq[iMov].getStackSize()];
				for(int i=0;i<imgT1Seq[iMov].getStackSize();i++)labels[i]=imgT1Seq[iMov].getStack().getSliceLabel(i+1);
				imgT1Seq[iMov]=trs[iMov].transformImage(imgT2Seq[0],imgT1Seq[iMov],false);
				imgT1Seq[iMov]=VitimageUtils.convertFloatToShortWithoutDynamicChanges(imgT1Seq[iMov]);
				for(int i=0;i<imgT1Seq[iMov].getStackSize();i++)imgT1Seq[iMov].getStack().setSliceLabel(labels[i],i+1);
			}
		}
		if(comp.isVisible())comp.close();

		
		
		//Adding the first image of T2 seq into T1 seq, as the last repetition time
		if(crossFit) {
			int nOld=nT1;
			ImagePlus[]imSeqNew=new ImagePlus[nT1+1];//TODO : verifier que les T1 succ ont les memes dims, pareil pour les T2
			for(int i=0;i<nOld;i++)imSeqNew[i]=VitimageUtils.imageCopy(imgT1Seq[i]);
			imSeqNew[nOld]=VitimageUtils.imageCopy(imgT2Seq[0]);
			for(int z=0;z<imSeqNew[nOld].getStackSize();z++)imSeqNew[nOld].getStack().setSliceLabel(imSeqNew[nOld].getStack().getSliceLabel(z+1).replace("T2SEQ", "T1SEQ"), z+1);
			imgT1Seq=imSeqNew;			
			int[]cpTab=new int[nOld+1];
			for(int i=0;i<nOld;i++)cpTab[i]=this.valT1Trs[i];
			cpTab[nOld]=TrValueUsedForT2Seq;
			this.valT1Trs=cpTab;
		}
				
		IJ.log("Registration finished !");
	}
	
	
	
	public void computeT1(Object[]params) {
		double sigmaSmoothing=0;
		if(params!=null) sigmaSmoothing=(double)params[0];
		else {
			double[]successiveSigmasSmoothing=new double[] {1};//{0.3,0.5,0.67,1.0,1.5,2.0,3.0,5.0};
			int nSucSig=successiveSigmasSmoothing.length;
			ImagePlus []resultsTab=new ImagePlus[nSucSig*2];
			for(int i=0;i<nSucSig;i++) {
				double valSig=successiveSigmasSmoothing[i];
				ImagePlus[]mapsComputed=MRUtils.computeT1Map(imgT1SeqTest, valSig);
				resultsTab[i*2]=mapsComputed[0];
				resultsTab[i*2].getStack().setSliceLabel("M0 map from T1 seq\nWith sigmaSmooth="+valSig,1);
				resultsTab[i*2+1]=mapsComputed[1];
				resultsTab[i*2+1].getStack().setSliceLabel("T1 map from T1 seq\nWith sigmaSmooth="+valSig,1);
			}
			ImagePlus test=HyperStackConverter.toHyperStack(new Concatenator().concatenate(resultsTab,false), nSucSig, 1,2,"xyztc","Fire");
			test.setTitle("T1 map test");
			test.show();
			int count=30;
			while((count--) >0 && test.isVisible()){
				test.setTitle("T1 seq results, "+count+" seconds before closing");
				VitimageUtils.waitFor(1000);
			}
			test.close();
			sigmaSmoothing=VitiDialogs.getDoubleUI("Choose sigma for smoothing", "",1);
		}
		if(debugCompute)T1map=imgT1Seq[imgT1Seq.length/2];
		else T1map=MRUtils.computeT1Map(imgT1Seq, sigmaSmoothing)[1];
	}
	
	
	public void computeT2(Object[]params) {
		double sigmaSmoothing=0;
		if(params!=null) sigmaSmoothing=(double)params[0];
		else {
			double[]successiveSigmasSmoothing=new double[] {1};//{0.3,0.5,0.67,1.0,1.5,2.0,3.0,5.0};
			int nSucSig=successiveSigmasSmoothing.length;
			ImagePlus []resultsTab=new ImagePlus[nSucSig*2];
			for(int i=0;i<nSucSig;i++) {
				double valSig=successiveSigmasSmoothing[i];
				ImagePlus[]mapsComputed=MRUtils.computeT2Map(imgT2SeqTest, valSig);
				resultsTab[i*2]=VitimageUtils.writeWhiteOnImage("M0 map from T2 seq\nWith sigmaSmooth="+valSig, mapsComputed[0], 7,0);
				resultsTab[i*2+1]=VitimageUtils.writeWhiteOnImage("T2 map from T2 seq\nWith sigmaSmooth="+valSig, mapsComputed[1], 7,0);
			}
			ImagePlus test=HyperStackConverter.toHyperStack(new Concatenator().concatenate(resultsTab,false), nSucSig, 1,2,"xyztc","Fire");
			test.show();
			int count=30;
			while((count--) >0 && test.isVisible()){
				test.setTitle("T1 seq results, "+count+" seconds before closing");
				VitimageUtils.waitFor(1000);
			}
			test.close();
			sigmaSmoothing=VitiDialogs.getDoubleUI("Choose sigma for smoothing", "",1);
		}
		ImagePlus []results=null;
		if(debugCompute) {
			T2map=imgT2Seq[imgT2Seq.length-1];
			M0map=imgT2Seq[0];
		}
		else {
			results= MRUtils.computeT2Map(imgT2Seq, sigmaSmoothing);
			T2map=results[1];
			M0map=results[0];
		}
	}
	
	//TODO : si T1 et T2, verifier que le tout matche en terme de XYZ
	public void stackSaveAndShowResults() {
		ImagePlus []tabIm=new ImagePlus[ (haveT1Seq ? (1+this.valT1Trs.length) : 0) +  (haveT2Seq ? (2+ this.valT2Tes.length) : 0) ];  
		int nZ=0;
		if(haveT2Seq)nZ=imgT2Seq[0].getNSlices();
		else if(haveT1Seq)nZ=imgT1Seq[0].getNSlices();
		int curInd=0;
		int nMaps=3;
		int nT1=this.valT1Trs.length;
		int nT2=this.valT2Tes.length;
		if(haveT2Seq) {//Write M0 map
			tabIm[curInd]=VitimageUtils.imageCopy(M0map);
			tabIm[curInd].setDisplayRange(0,maxDisplayedM0 );
			IJ.run(tabIm[curInd],"32-bit","");
			for(int z=1;z<=nZ;z++) {
				tabIm[curInd].getStack().setSliceLabel("M0MAP", VitimageUtils.getCorrespondingSliceInHyperImage(tabIm[curInd], 0, z-1, 0));
			}
			curInd++;
		}
		if(haveT1Seq) {//Write T1 map
			tabIm[curInd]=VitimageUtils.imageCopy(T1map);
			tabIm[curInd].setDisplayRange(0,maxDisplayedT1 );
			IJ.run(tabIm[curInd],"32-bit","");
			for(int z=1;z<=nZ;z++) {
				tabIm[curInd].getStack().setSliceLabel("T1MAP", VitimageUtils.getCorrespondingSliceInHyperImage(tabIm[curInd], 0, z-1, 0));
			}
			curInd++;
		}
		if(haveT2Seq) {//Write T2 map
			tabIm[curInd]=VitimageUtils.imageCopy(T2map);
			tabIm[curInd].setDisplayRange(0,maxDisplayedT2 );
			IJ.run(tabIm[curInd],"32-bit","");
			for(int z=1;z<=nZ;z++) {
				tabIm[curInd].getStack().setSliceLabel("T2MAP", VitimageUtils.getCorrespondingSliceInHyperImage(tabIm[curInd], 0, z-1, 0));
			}
			curInd++;
		}

		if(haveT1Seq) {//Write T1seq
			for(int i=0;i<this.valT1Trs.length;i++) {
				tabIm[curInd]=VitimageUtils.imageCopy(imgT1Seq[i]);
				tabIm[curInd].setDisplayRange(0,maxDisplayedM0 );
				IJ.run(tabIm[curInd++],"32-bit","");
			}
		}
		if(haveT2Seq) {//Write T2 seq
			for(int i=0;i<this.valT2Tes.length;i++) {
				tabIm[curInd]=VitimageUtils.imageCopy(imgT2Seq[i]);
				tabIm[curInd].setDisplayRange(0,maxDisplayedM0 );
				IJ.run(tabIm[curInd++],"32-bit","");
			}
		}

		int nC=tabIm.length;
		ImagePlus hyperImg=Concatenator.run(tabIm);
		hyperImg=HyperStackConverter.toHyperStack(hyperImg, nC,nZ,1,"xyztc","Grayscale");

		curInd=1;
		if(haveT2Seq) {
			hyperImg.setC(curInd++);
			hyperImg.setDisplayRange(0,maxDisplayedM0 );
			IJ.run(hyperImg,"Fire","");
		}

		if(haveT1Seq) {
			hyperImg.setC(curInd++);
			hyperImg.setDisplayRange(0,maxDisplayedT1 );
			IJ.run(hyperImg,"Fire","");
		}

		if(haveT2Seq) {
			hyperImg.setC(curInd++);
			hyperImg.setDisplayRange(0,maxDisplayedT2 );
			IJ.run(hyperImg,"Fire","");
		}

		for(int c=curInd;c<=nC;c++) {
			hyperImg.setC(c);
			hyperImg.setDisplayRange(0,maxDisplayedM0 );				
			IJ.run(hyperImg,"Fire","");
		}
		hyperImg.setTitle("Hyperimage");
		IJ.saveAsTiff(hyperImg,new File(outputDir,"hyperImage.tif").getAbsolutePath());
		hyperImg.close();
		ImagePlus result=IJ.openImage(new File(outputDir,"hyperImage.tif").getAbsolutePath());
		result.show();
	}
}
