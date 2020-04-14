package com.vitimage.mrutils;

import com.vitimage.aplimtools.Bionano_T1T2M0_Importer;
import com.vitimage.common.TransformUtils;
import com.vitimage.common.VitiDialogs;
import com.vitimage.common.VitimageUtils;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.plugin.Duplicator;
import ij.plugin.filter.Convolver;
import ij.process.ImageProcessor;


public class HyperMRIT1T2 {
	final static double bionanoCapillaryRadius=1.15/2;//size in mm in the bionanoNMRI lab
	ImagePlus hyperImg;
	public int nTimepoints;
	public int nChannels;
	final int m0MapChannel=0;
	final int t1MapChannel=1;
	final int t2mapChannel=2;
	public boolean hasT2sequence;
	public boolean hasT1sequence;
	public boolean hasT1T2sequence;
	public int[]dims;
	public double[]voxs;
	public double[][][]Te;
	public double[][][]Tr;
	public double [][][]sigmaRice;
	public MRDataType [][] mrDataType;
	public String[]actualDay;
	public double[][]tabSigmasT1Seq;
	public double[][]tabSigmasT2Seq;
	public double[][]tabSigmasT1T2Seq;
	private double TR_ALPHA=0;
	private double TE_ALPHA=0;

	public static void main(String[]args) {
		ImageJ ij=new ImageJ();
		HyperMRIT1T2 hyp=new HyperMRIT1T2(IJ.openImage("/home/fernandr/Bureau/Traitements/Sorgho/Test SSM1/Out_input_Fijiyama/Output/Exported_data/Data_combined.tif"));
		hyp.showCopy("Initial Hyper MRI");
		hyp.normalizeMRSignal();
		hyp.showCopy("Normalized Hyper MRI");
	}

	
	public void showCopy(String text) {
		ImagePlus temp=hyperImg.duplicate();
		temp.show();
		temp.setTitle(text);
	}
	
	
	public HyperMRIT1T2(ImagePlus hyperImg) {
		this.hyperImg=hyperImg;
		setupStructures();
		fillInInformationFields();
		computeSigmaTabs();
	}

	public void setupStructures() {
		this.nTimepoints=hyperImg.getNFrames();
		this.actualDay=new String[this.nTimepoints];
		for(int t=0;t<this.nTimepoints;t++)this.actualDay[t]=""+t;
		this.nChannels=hyperImg.getNChannels();
		this.dims=VitimageUtils.getDimensions(hyperImg);
		this.voxs=VitimageUtils.getVoxelSizes(hyperImg);
		
		this.Te=new double[this.nTimepoints][dims[2]][this.nChannels];
		this.Tr=new double[this.nTimepoints][dims[2]][this.nChannels];
		this.sigmaRice=new double[this.nTimepoints][dims[2]][this.nChannels];
		this.mrDataType=new MRDataType[this.nTimepoints][this.nChannels];
	}
	
	public void fillInInformationFields() {
		for(int t=0;t<this.nTimepoints;t++) {
			hasT1sequence=false;
			hasT2sequence=false;
			hasT1T2sequence=false;
			for(int c=0;c<this.nChannels;c++) {
				for(int z=0;z<dims[2];z++) {
					mrDataType[t][c]=MRUtils.getDataTypeOfThisMagneticResonanceSlice(hyperImg, c, z, t);
					if(mrDataType[t][c]!=MRDataType.T1SEQ  && mrDataType[t][c]!=MRDataType.T2SEQ  && mrDataType[t][c]!=MRDataType.T1T2SEQ) {
						sigmaRice[t][z][c]=Float.NaN;
						Te[t][z][c]=Float.NaN;
						Tr[t][z][c]=Float.NaN;
					}
					else {
						if(mrDataType[t][c]==MRDataType.T1SEQ)hasT1sequence=true;
						if(mrDataType[t][c]==MRDataType.T1T2SEQ)hasT1T2sequence=true;
						if(mrDataType[t][c]==MRDataType.T2SEQ)hasT2sequence=true;
						sigmaRice[t][z][c]=MRUtils.readSigmaInSliceLabel(hyperImg, c, z,t);
						Te[t][z][c]=MRUtils.readTeInSliceLabel(hyperImg, c, z, t);
						Tr[t][z][c]=MRUtils.readTrInSliceLabel(hyperImg, c, z, t);		
					}
				}
			}
			if(((!hasT1sequence) || (!hasT2sequence)) && (!hasT1T2sequence)) {
				VitiDialogs.notYet("Handling only T1 or only T2 data : not yet coded in mrutils/HyperMRIT1T2.java");
				System.exit(0);
			}
		}
	}

	
	public void computeSigmaTabs() {
		if(hasT1T2sequence) {
			this.tabSigmasT1T2Seq=new double[this.nTimepoints][dims[2]];			
			for(int t=0;t<this.nTimepoints;t++) {
				for(int z=0;z<dims[2];z++) {
					double average=0;
					for(int ind=3;ind<nChannels;ind++)average+=MRUtils.readSigmaInSliceLabel(hyperImg, ind, z, t);
					this.tabSigmasT1T2Seq[t][z]=average/(nChannels-3);
				}
			}
			IJ.log(TransformUtils.stringMatrixMN("Tab sigmas T1T2 seq", tabSigmasT1T2Seq));
		}
		else{
			this.tabSigmasT1Seq=new double[this.nTimepoints][dims[2]];
			this.tabSigmasT2Seq=new double[this.nTimepoints][dims[2]];
			for(int t=0;t<this.nTimepoints;t++) {
				for(int z=0;z<dims[2];z++) {
					double average=0;
					int[]indices=getT1Indices(t);
					for(int ind:indices)average+=MRUtils.readSigmaInSliceLabel(hyperImg, ind, z, t);
					this.tabSigmasT1Seq[t][z]=average/indices.length;
					
					average=0;
					indices=getT2Indices(t);
					for(int ind:indices)average+=MRUtils.readSigmaInSliceLabel(hyperImg, ind, z, t);
					this.tabSigmasT2Seq[t][z]=average/indices.length;
				}
			}
			IJ.log(TransformUtils.stringMatrixMN("Tab sigmas T1 seq", tabSigmasT1Seq));
			IJ.log(TransformUtils.stringMatrixMN("Tab sigmas T2 seq", tabSigmasT2Seq));
		}
	
	}
	
	
	/** This function normalize MR signal along Z and T axis, using the capillary values*/
	public void normalizeMRSignal() {
		double[][]capillaryValuesInM0Map=new double[nTimepoints][dims[2]];
		double factor=1;
		//Measure capillary values in T1 and T2 sequences, at each time
		//Adjust image dynamics of M0 map, T1seq, and T2 seq in order that value being 10000;
		for(int t=0;t<this.nTimepoints;t++) {
			capillaryValuesInM0Map[t]=measureCapillaryValuesInM0Map(t);
			for(int c=0;c<this.nChannels;c++) {
				if((c)!=t1MapChannel && (c)!=t2mapChannel) {
					for(int z=0;z<dims[2];z++) {
						int slice=VitimageUtils.getCorrespondingSliceInHyperImage(hyperImg, c, z, t);
						factor=MRUtils.maxM0ForNormalization/capillaryValuesInM0Map[t][z];
						ImageProcessor ip=hyperImg.getStack().getProcessor(slice);
						ip.multiply(factor);
						if((c)!=m0MapChannel)MRUtils.modifySigma(hyperImg,c,z,t,factor*MRUtils.readSigmaInSliceLabel(hyperImg, c, z, t));
					}
					hyperImg.setC(c+1);
					hyperImg.setDisplayRange(0, MRUtils.maxDisplayedM0);
				}
			}
		}		
	}
	
	
	public static double measureMeanCapillaryValueAlongZ(ImagePlus imgM0) {
		//Find capillary in the central slice
		ImagePlus img=new Duplicator().run(imgM0,1,1,imgM0.getNSlices()/2+1,imgM0.getNSlices()/2+1,1,1);
		int []coordsCentral=findCapillaryCenterInSlice(img,bionanoCapillaryRadius);
		
		//Find best match for the capillary in other slices, in a neighbourhood around, and measure its value
		double[]vals=capillaryValuesAlongZStatic(imgM0,coordsCentral,bionanoCapillaryRadius);
		System.out.println("Mean capillary values="+TransformUtils.stringVectorN(vals, ""));
		return VitimageUtils.statistics1D(vals)[0];
	}

	public double[]measureCapillaryValuesInM0Map(int time){

		//Find capillary in the central slice
		ImagePlus img=new Duplicator().run(hyperImg,m0MapChannel+1,m0MapChannel+1,dims[2]/2+1,dims[2]/2+1,time+1,time+1);
		int []coordsCentral=findCapillaryCenterInSlice(img,bionanoCapillaryRadius);
		
		//Find best match for the capillary in other slices, in a neighbourhood around, and measure its value
		return capillaryValuesAlongZ(time,coordsCentral,bionanoCapillaryRadius);
	}

	
	public double[]capillaryValuesAlongZ(int time,int[]coords,double capillaryRadius){
		IJ.log("\nCapillary detection at time="+time+".\nStarting from coordinates :"+coords[0]+", "+coords[1]);
		Duplicator dup=new Duplicator();
		int rayPix=(int)Math.round(capillaryRadius/(voxs[0]));
		int semiRayPix=rayPix/2;
		int lookupRadius=(int)Math.round(1.5*rayPix);

		int xMed=coords[0];int yMed=coords[1];int zMed=dims[2]/2;
		double[]capVals=new double[dims[2]];
		ImagePlus imgTemp,img3D;
		int xLast=xMed,yLast=yMed;
		long t0= System.currentTimeMillis();

		img3D=dup.run(hyperImg,m0MapChannel+1,m0MapChannel+1,1,dims[2],time+1,time+1);
		for(int z=zMed;z>=0;z--) {
			//Extract a patch around the finding of upside
			imgTemp=VitimageUtils.cropImageFloat(img3D,xLast-lookupRadius,yLast-lookupRadius,z, lookupRadius*2,lookupRadius*2,1);
			//Find cap center in it
			int[]coordsNew=findCapillaryCenterInSlice(imgTemp,bionanoCapillaryRadius);
			
			//Update coordinates of last
			xLast=xLast-lookupRadius+coordsNew[0];
			yLast=yLast-lookupRadius+coordsNew[1];

			//Gather information of M0 value
			double[]stats=VitimageUtils.statistics1D(VitimageUtils.valuesOfBlock(img3D,xLast-semiRayPix, yLast-semiRayPix, z-1,xLast+semiRayPix, yLast+semiRayPix, z+1));
			IJ.log("Capillary detected at z="+z+" at coordinates "+xLast+", "+yLast+" with M0="+stats[0]+" std="+stats[1]);
			capVals[z]=stats[0];
		}
		xLast=xMed;
		yLast=yMed;
		for(int z=zMed;z<dims[2];z++) {
			//Extract a patch around the finding of upside
			imgTemp=VitimageUtils.cropImageFloat(img3D,xLast-lookupRadius,yLast-lookupRadius, z, lookupRadius*2,lookupRadius*2,1);
			
			//Find cap center in it
			int[]coordsNew=findCapillaryCenterInSlice(imgTemp,bionanoCapillaryRadius);
			
			//Update coordinates of last
			xLast=xLast-lookupRadius+coordsNew[0];
			yLast=yLast-lookupRadius+coordsNew[1];
			
			//Gather information of M0 value
			double[]stats=VitimageUtils.statistics1D(VitimageUtils.valuesOfBlock(img3D,xLast-semiRayPix, yLast-semiRayPix, z-1,xLast+semiRayPix, yLast+semiRayPix, z+1));
			System.out.println("Capillary detected at z="+z+" at coordinates "+xLast+", "+yLast+" with M0="+stats[0]+" std="+stats[1]);
			capVals[z]=stats[0];
		}
		return capVals;
	}

	
	public static double[]capillaryValuesAlongZStatic(ImagePlus img3D,int[]coords,double capillaryRadius){
		double[]voxs=VitimageUtils.getVoxelSizes(img3D);
		int[]dims=VitimageUtils.getDimensions(img3D);
		int rayPix=(int)Math.round(capillaryRadius/(voxs[0]));
		int semiRayPix=rayPix/2;
		int lookupRadius=(int)Math.round(1.5*rayPix);

		int xMed=coords[0];int yMed=coords[1];int zMed=dims[2]/2;
		double[]capVals=new double[dims[2]];
		ImagePlus imgTemp;
		int xLast=xMed,yLast=yMed;

		for(int z=zMed;z>=0;z--) {
			//Extract a patch around the finding of upside
			imgTemp=VitimageUtils.cropImageFloat(img3D,xLast-lookupRadius,yLast-lookupRadius,z, lookupRadius*2,lookupRadius*2,1);
			//Find cap center in it
			int[]coordsNew=findCapillaryCenterInSlice(imgTemp,bionanoCapillaryRadius);
			
			//Update coordinates of last
			xLast=xLast-lookupRadius+coordsNew[0];
			yLast=yLast-lookupRadius+coordsNew[1];

			//Gather information of M0 value
			double[]stats=VitimageUtils.statistics1D(VitimageUtils.valuesOfBlock(img3D,xLast-semiRayPix, yLast-semiRayPix, z-1,xLast+semiRayPix, yLast+semiRayPix, z+1));
			IJ.log("Capillary detected at z="+z+" at coordinates "+xLast+", "+yLast+" with M0="+stats[0]+" std="+stats[1]);
			capVals[z]=stats[0];
		}
		xLast=xMed;
		yLast=yMed;
		for(int z=zMed;z<dims[2];z++) {
			//Extract a patch around the finding of upside
			imgTemp=VitimageUtils.cropImageFloat(img3D,xLast-lookupRadius,yLast-lookupRadius, z, lookupRadius*2,lookupRadius*2,1);
			
			//Find cap center in it
			int[]coordsNew=findCapillaryCenterInSlice(imgTemp,bionanoCapillaryRadius);
			
			//Update coordinates of last
			xLast=xLast-lookupRadius+coordsNew[0];
			yLast=yLast-lookupRadius+coordsNew[1];
			
			//Gather information of M0 value
			double[]stats=VitimageUtils.statistics1D(VitimageUtils.valuesOfBlock(img3D,xLast-semiRayPix, yLast-semiRayPix, z-1,xLast+semiRayPix, yLast+semiRayPix, z+1));
			System.out.println("Capillary detected at z="+z+" at coordinates "+xLast+", "+yLast+" with M0="+stats[0]+" std="+stats[1]);
			capVals[z]=stats[0];
		}
		return capVals;
	}

	
	public static int[] findCapillaryCenterInSlice(ImagePlus imgTemp,double capillaryRadius) {
		ImagePlus img=imgTemp.duplicate();
		double[]voxs=VitimageUtils.getVoxelSizes(img);		
		int rayPix=(int)Math.round(capillaryRadius/(voxs[0]));
		int raySquare=rayPix*rayPix;
		int kernelSize=3*rayPix;
		if (kernelSize%2==0)kernelSize++;
		int distX,distY;
		float[]kernel=new float[kernelSize*kernelSize];
		for(int i=0;i<kernelSize;i++)for(int j=0;j<kernelSize;j++) {
			distX=(i-(kernelSize/2));
			distY=(j-(kernelSize/2));
			kernel[i*kernelSize+j]=((distX*distX+distY*distY)<raySquare) ? 1 : -1;
		}
		new Convolver().convolve(img.getProcessor(), kernel, kernelSize, kernelSize);
		img=VitimageUtils.makeOperationOnOneImage(img,2,-1, false);
		int[]coordsOfMax=VitimageUtils.indMaxOfImage(img);
		return coordsOfMax;
	}

	
	public Object[] getMapsImage(boolean writeParametersForExplorer) {
		int fontSize=15;
		ImagePlus maps= new Duplicator().run(hyperImg,1,3,1,dims[2],1,nTimepoints);
		String[][][]texts=new String[3][dims[2]][nTimepoints];
		for(int z=0;z<dims[2];z++)for(int t=0;t<nTimepoints;t++)for(int c=0;c<3;c++) {
			texts[c][z][t]="   "+MRUtils.getDataTypeOfThisMagneticResonanceSlice(maps, c, z, t)+"  Day="+actualDay[t]+"    Z="+z;
			if(writeParametersForExplorer) VitimageUtils.writeTextOnGivenImage(texts[c][z][t],maps ,fontSize ,30 ,dims[1]-30,c, z, t);
		}  
		maps.setPosition(VitimageUtils.getCorrespondingSliceInHyperImage(maps, 0, dims[2]/2, 0));
		return new Object[] {maps,texts};
	}


	public Object[] getEchoesImage(boolean writeParametersForExplorer) {
		int fontSize=15;
		ImagePlus echoes= new Duplicator().run(hyperImg,4,nChannels,1,dims[2],1,nTimepoints);
		String[][][]texts=new String[nChannels-3][dims[2]][nTimepoints];
		for(int z=0;z<dims[2];z++)for(int t=0;t<nTimepoints;t++) for(int c=0;c<nChannels-3;c++) {
			texts[c][z][t]="   "+MRUtils.getDataTypeOfThisMagneticResonanceSlice(echoes, c, z, t)+"  Day="+actualDay[t]+"   Z="+z+
						"  TR="+(MRUtils.readTrInSliceLabel(echoes, c, z, t))+"  TE="+(MRUtils.readTeInSliceLabel(echoes, c, z, t));
			if(writeParametersForExplorer) VitimageUtils.writeTextOnGivenImage(texts[c][z][t],echoes ,fontSize ,30 ,10,c , z, t);
			//System.out.println("Just read tr="+MRUtils.readTrInSliceLabel(echoes, c, z, t)+"te="+MRUtils.readTeInSliceLabel(echoes, c, z, t));
		}
		echoes.setPosition(VitimageUtils.getCorrespondingSliceInHyperImage(echoes, 1, dims[2]/2, 0));
		return new Object[] {echoes,texts};
	}

	
	public int getT1T2SeqNumberReps(int t){
		int nb=0;
		for(int c=0;c<this.nChannels;c++)if(mrDataType[t][c]==MRDataType.T1T2SEQ)nb++;
		return nb;
	}

	
	public int getT1SeqNumberReps(int t){
		int nb=0;
		for(int c=0;c<this.nChannels;c++)if(mrDataType[t][c]==MRDataType.T1SEQ)nb++;
		return nb;
	}
		
	public int getT2SeqNumberReps(int t){
		int nb=0;
		for(int c=0;c<this.nChannels;c++)if(mrDataType[t][c]==MRDataType.T2SEQ)nb++;
		return nb;
	}


	
	public int[]getT1Indices(int t){
		int[]ret=new int[getT1SeqNumberReps(t)];
		int incr=0;
		for(int c=0;c<this.nChannels;c++)if(mrDataType[t][c]==MRDataType.T1SEQ)ret[incr++]=c;		
		return ret;
	}

	public int[]getT2Indices(int t){
		int[]ret=new int[getT2SeqNumberReps(t)];
		int incr=0;
		for(int c=0;c<this.nChannels;c++)if(mrDataType[t][c]==MRDataType.T2SEQ)ret[incr++]=c;		
		return ret;
	}

	public int[]getT1T2Indices(int t){
		int[]ret=new int[getT1T2SeqNumberReps(t)];
		int incr=0;
		for(int c=0;c<this.nChannels;c++)if(mrDataType[t][c]==MRDataType.T1T2SEQ)ret[incr++]=c;		
		return ret;
	}

	public int[][]getT1Indices(){
		int[][]ret=new int[this.nTimepoints][];
		for(int t=0;t<this.nTimepoints;t++)ret[t]=this.getT1Indices(t);
		return ret;
	}

	public int[][]getT2Indices(){
		int[][]ret=new int[this.nTimepoints][];
		for(int t=0;t<this.nTimepoints;t++)ret[t]=this.getT2Indices(t);
		return ret;
	}

	public int[][]getT1T2Indices(){
		int[][]ret=new int[this.nTimepoints][];
		for(int t=0;t<this.nTimepoints;t++)ret[t]=this.getT1T2Indices(t);
		return ret;
	}

/*
	public double[]getT1TrTimes(int t){
		int incr=0;
		double[]ret=new double[getT1SeqNumberReps(t)];
		for(int c=0;c<this.nChannels;c++)if(mrDataType[t][c]==MRDataType.T1SEQ)ret[incr++]=this.Tr[t][c]-TR_ALPHA;		
		return ret;
	}

	public double[][]getT1TrTimes(){
		double[][]ret=new double[this.nTimepoints][];
		for(int t=0;t<this.nTimepoints;t++)ret[t]=getT1TrTimes(t);
		return ret;
	}

	public double[]getT1T2TrTimes(int t){
		int incr=0;
		double[]ret=new double[getT1T2SeqNumberReps(t)];
		for(int c=0;c<this.nChannels;c++)if(mrDataType[t][c]==MRDataType.T1T2SEQ)ret[incr++]=this.Tr[t][c]-TR_ALPHA;		
		return ret;
	}

	public double[][]getT1T2TrTimes(){
		double[][]ret=new double[this.nTimepoints][];
		for(int t=0;t<this.nTimepoints;t++)ret[t]=getT1T2TrTimes(t);
		return ret;
	}
*/
	public double[][][]getT1T2TeTimes(){
		double[][][][]init=getT1T2TrTeTimes();
		double[][][]ret=new double[this.nTimepoints][dims[2]][];
		for(int t=0;t<this.nTimepoints;t++)for(int z=0;z<dims[2];z++) {
			ret[t][z]=new double[init[t][z].length];
			for(int c=0;c<init[t][z].length;c++) {
				ret[t][z][c]=init[t][z][c][1];
			}
		}
		return ret;
	}

	public double[][][]getT1T2TrTimes(){
		double[][][][]init=getT1T2TrTeTimes();
		double[][][]ret=new double[this.nTimepoints][dims[2]][];
		for(int t=0;t<this.nTimepoints;t++)for(int z=0;z<dims[2];z++) {
			ret[t][z]=new double[init[t][z].length];
			for(int c=0;c<init[t][z].length;c++) {
				ret[t][z][c]=init[t][z][c][0];
			}
		}
		return ret;
	}
	
	
	public double[][][]getT1T2TrTeTimes(int t){
		int incr=0;
		double[][][]ret=new double[dims[2]][getT1T2SeqNumberReps(t)][2];
		for(int c=0;c<this.nChannels;c++)if(mrDataType[t][c]==MRDataType.T1T2SEQ) {
			for(int z=0;z<dims[2];z++)
				ret[z][incr]=new double[] {this.Tr[t][z][c],this.Te[t][z][c]};
			incr++;
		}
		return ret;
	}
	
	public double[][][][]getT1T2TrTeTimes(){
		double[][][][]ret=new double[nTimepoints][][][];
		for(int t=0;t<this.nTimepoints;t++)ret[t]=getT1T2TrTeTimes(t);
		return ret;
	}
	/*
	public double[]getT2TeTimes(int t){
		double[]ret=new double[getT2SeqNumberReps(t)];
		int incr=0;
		for(int c=0;c<this.nChannels;c++)if(mrDataType[t][c]==MRDataType.T2SEQ)ret[incr++]=this.Te[t][c]+TR_ALPHA;		
		return ret;
	}

	public double[][]getT2TeTimes(){
		double[][]ret=new double[this.nTimepoints][];
		for(int t=0;t<this.nTimepoints;t++)ret[t]=getT2TeTimes(t);
		return ret;
	}

	public double[]getT1T2TeTimes(int t){
		double[]ret=new double[getT1T2SeqNumberReps(t)];
		int incr=0;
		for(int c=0;c<this.nChannels;c++)if(mrDataType[t][c]==MRDataType.T1T2SEQ)ret[incr++]=this.Te[t][c]+(this.Tr[t][c]>9000 ? TE_ALPHA : 0);		
		return ret;
	}

*/
	/*
	public double[][][]getMRISignalForThisVoxel(int x,int y, int z,int t){
		return getMeanMRISignalAroundThisVoxel(x, y, z, t, 0,0);
	}*/
/*	
	public double[][][]getMeanMRISignalAroundThisVoxel(int x,int y, int z,int t,int crossWidth,int crossThick){
		int xm,ym,xM,yM,zm,zM;
		xm=Math.max(0, x-crossWidth);
		xM=Math.min(dims[0]-1, x+crossWidth);
		ym=Math.max(0, y-crossWidth);
		yM=Math.min(dims[1]-1, y+crossWidth);
		zm=Math.max(0, z-crossThick);
		zM=Math.min(dims[2]-1, z+crossThick);

		double[]t1Times=this.getT1T2TrTeTimes(t);
		double[]t2Times=this.getT2TeTimes(t);
		int[]t1Indices=this.getT1Indices(t);
		int[]t2Indices=this.getT2Indices(t);
		int nHits=(xM-xm+1)*(yM-ym+1)*(zM-zm+1);
		int currentChan;
		double[]mrDataT1Seq= new double[t1Times.length];
		double[]mrDataT2Seq= new double[t2Times.length];
		double weightForEachVoxel=1.0/nHits;
		float[]pixels;
		for(int t1trInd=0;t1trInd<t1Indices.length;t1trInd++) {			
			currentChan=t1Indices[t1trInd];
			for(int zz=zm;zz<=zM;zz++) {
				int indexSlice=this.nChannels*dims[2]*t + this.nChannels*zz + currentChan + 1;
				pixels=((float[])(this.hyperImg.getStack().getProcessor(indexSlice).getPixels()));
				for(int xx=xm;xx<=xM;xx++) {
					for(int yy=ym;yy<=yM;yy++) {
						mrDataT1Seq[t1trInd]+=weightForEachVoxel*(double)pixels[xx + this.dims[0] * yy];
					}
				}
			}
		}
		for(int t2teInd=0;t2teInd<t1Indices.length;t2teInd++) {			
			currentChan=t2Indices[t2teInd];
			for(int zz=zm;zz<=zM;zz++) {
				int indexSlice=this.nChannels*dims[2]*t + this.nChannels*zz + currentChan + 1;
				pixels=((float[])(this.hyperImg.getStack().getProcessor(indexSlice).getPixels()));
				for(int xx=xm;xx<=xM;xx++) {
					for(int yy=ym;yy<=yM;yy++) {
						mrDataT2Seq[t2teInd]+=weightForEachVoxel*(double)pixels[xx + this.dims[0] * yy];
					}
				}
			}
		}
		
		return new double[][][] {{t1Times,mrDataT1Seq},{t2Times,mrDataT2Seq}};		
	}
*/
/*
	public double[][][]getMeanMRISignalAroundThisVoxelT1T2(int x,int y, int z,int t,int crossWidth,int crossThick){
		int xm,ym,xM,yM,zm,zM;
		xm=Math.max(0, x-crossWidth);
		xM=Math.min(dims[0]-1, x+crossWidth);
		ym=Math.max(0, y-crossWidth);
		yM=Math.min(dims[1]-1, y+crossWidth);
		zm=Math.max(0, z-crossThick);
		zM=Math.min(dims[2]-1, z+crossThick);

		double[]t1t2TrTimes=this.getT1T2TrTimes(t);
		double[]t1t2TeTimes=this.getT1T2TeTimes(t);
		int[]t1t2Indices=this.getT1T2Indices(t);
		int nHits=(xM-xm+1)*(yM-ym+1)*(zM-zm+1);
		int currentChan;
		double[]mrDataT1T2Seq= new double[t1t2TrTimes.length];
		double weightForEachVoxel=1.0/nHits;
		float[]pixels;
		for(int echoInd=0;echoInd<t1t2Indices.length;echoInd++) {			
			currentChan=t1t2Indices[echoInd];
			for(int zz=zm;zz<=zM;zz++) {
				int indexSlice=this.nChannels*dims[2]*t + this.nChannels*zz + currentChan + 1;
				pixels=((float[])(this.hyperImg.getStack().getProcessor(indexSlice).getPixels()));
				for(int xx=xm;xx<=xM;xx++) {
					for(int yy=ym;yy<=yM;yy++) {
						mrDataT1T2Seq[echoInd]+=weightForEachVoxel*(double)pixels[xx + this.dims[0] * yy];
					}
				}
			}
		}
		
		return new double[][][] {{t1t2TrTimes},{t1t2TeTimes},{mrDataT1T2Seq}};		
	}
	
	public double[][][][]getFullMRISignalAroundThisVoxel(int x,int y, int z,int crossWidth,int crossThick){
		int xm,ym,xM,yM,zm,zM;
		xm=Math.max(0, x-crossWidth);
		xM=Math.min(dims[0]-1, x+crossWidth);
		ym=Math.max(0, y-crossWidth);
		yM=Math.min(dims[1]-1, y+crossWidth);
		zm=Math.max(0, z-crossThick);
		zM=Math.min(dims[2]-1, z+crossThick);

		int[][]t1Indices=this.getT1Indices();
		int[][]t2Indices=this.getT2Indices();
		int nHits=(xM-xm+1)*(yM-ym+1)*(zM-zm+1);
		double[][][][]data= new double[nTimepoints][nHits][2][];//[times][vox][Seq][echos]
		
		int currentChan;
		int index;
		
		for(int t=0;t<this.nTimepoints;t++) {			
			for(int zz=zm;zz<=zM;zz++) {
				for(int xx=xm;xx<=xM;xx++) {
					for(int yy=ym;yy<=yM;yy++) {
						index=(xx-xm) + (xM-xm+1)*(yy-ym) + (xM-xm+1)*(yM-ym+1)*(zz-zm);
						data[t][index][0]=new double[t1Indices[t].length];
						for(int t1trInd=0;t1trInd<t1Indices[t].length;t1trInd++) {			
							currentChan=t1Indices[t][t1trInd];
							int indexSlice=this.nChannels*dims[2]*t + this.nChannels*zz + currentChan + 1;
							data[t][index][0][t1trInd]=hyperImg.getStack().getProcessor(indexSlice).getPixelValue(xx, yy);
						}
						data[t][index][1]=new double[t2Indices[t].length];
						for(int t2teInd=0;t2teInd<t2Indices[t].length;t2teInd++) {			
							currentChan=t2Indices[t][t2teInd];
							int indexSlice=this.nChannels*dims[2]*t + this.nChannels*zz + currentChan + 1;
							data[t][index][1][t2teInd]=hyperImg.getStack().getProcessor(indexSlice).getPixelValue(xx, yy);
						}
					}
				}	
			}
		}
		return data;
	}
*/
	
	
	public double[][][][]getFullMRISignalAroundThisVoxelT1T2(int x,int y, int z,int crossWidth,int crossThick,double sigmaXYInVoxels){
		int xm,ym,xM,yM,zm,zM;
		xm=Math.max(0, x-crossWidth);
		xM=Math.min(dims[0]-1, x+crossWidth);
		ym=Math.max(0, y-crossWidth);
		yM=Math.min(dims[1]-1, y+crossWidth);
		zm=Math.max(0, z-crossThick);
		zM=Math.min(dims[2]-1, z+crossThick);

		int xm2,ym2,xM2,yM2,zm2,zM2;
		xm2=Math.max(0, x-crossWidth-5);
		xM2=Math.min(dims[0]-1, x+crossWidth+5);
		ym2=Math.max(0, y-crossWidth-5);
		yM2=Math.min(dims[1]-1, y+crossWidth+5);

		int deltaX=xm-xm2;
		int deltaY=ym-ym2;
		
		int[][]t1t2Indices=this.getT1T2Indices();
		int nHits=(xM-xm+1)*(yM-ym+1)*(zM-zm+1);
		double[][][][]data= new double[nTimepoints][nHits][1][];//[times][vox][Seq][echos]
		
		float[]pixels;
		int currentChan;
		int index;
		
		ImagePlus temp=VitimageUtils.cropMultiChannelFloatImage(hyperImg,xm2,xM2,ym2,yM2,zm,zM);
		temp=VitimageUtils.gaussianFilteringMultiChannel(temp,sigmaXYInVoxels,sigmaXYInVoxels,0);//It's no error : no "big smoothing" over Z, due to misalignment
		
		for(int t=0;t<this.nTimepoints;t++) {			
			for(int zz=zm;zz<=zM;zz++) {
				for(int xx=xm;xx<=xM;xx++) {
					for(int yy=ym;yy<=yM;yy++) {
						index=(xx-xm) + (xM-xm+1)*(yy-ym) + (xM-xm+1)*(yM-ym+1)*(zz-zm);
						data[t][index][0]=new double[t1t2Indices[t].length];
						for(int t1trInd=0;t1trInd<t1t2Indices[t].length;t1trInd++) {			
							currentChan=t1t2Indices[t][t1trInd];
							//int indexSlice=this.nChannels*dims[2]*t + this.nChannels*zz + currentChan + 1;
							//data[t][index][0][t1trInd]=hyperImg.getStack().getProcessor(indexSlice).getPixelValue(xx, yy);
							int indexSlice=this.nChannels*(zM-zm+1)*t + this.nChannels*(zz-zm) + currentChan + 1;
							data[t][index][0][t1trInd]=temp.getStack().getProcessor(indexSlice).getPixelValue(xx-xm2, yy-ym2);
							System.out.println("data["+t+"]["+index+"][0]["+t1trInd+"]="+data[t][index][0][t1trInd]);
						}
					}
				}	
			}
		}
		return data;
	}


	
	
	/*
	
	public double[][][][]getFullMRISignalInTheseCoordinates(int xCor,int yCor, int zCor,int[][]coordinates){
		double[][]t1Times=this.getT1T2TrTeTimes();
		double[][]t2Times=this.getT2TeTimes();
		int[][]t1Indices=this.getT1Indices();
		int[][]t2Indices=this.getT2Indices();
		int nHits=coordinates.length;
		double[][][][]data= new double[t1Times.length][nHits][2][];//[times][vox][Seq][echos]		
		int currentChan;
		int index,xx,yy,zz;		
		for(int t=0;t<this.nTimepoints;t++) {			
			for(int vo=0;vo<nHits;vo++) {
				xx=coordinates[vo][0];
				yy=coordinates[vo][1];
				zz=zCor;
				data[t][vo][0]=new double[t1Indices[t].length];
				data[t][vo][1]=new double[t2Indices[t].length];
				//System.out.println("Processing ROI coordinates : zCor="+zCor+" t="+t+", pts="+xx+","+yy);
				for(int t1trInd=0;t1trInd<t1Indices[t].length;t1trInd++) {			
					currentChan=t1Indices[t][t1trInd];
					int indexSlice=this.nChannels*dims[2]*t + this.nChannels*zz + currentChan + 1;
					data[t][vo][0][t1trInd]=hyperImg.getStack().getProcessor(indexSlice).getPixelValue(xx, yy);
				}
				//System.out.println(TransformUtils.stringVectorN(data[t][vo][0],"T1seq"));
				for(int t2teInd=0;t2teInd<t2Indices[t].length;t2teInd++) {			
					currentChan=t2Indices[t][t2teInd];
					int indexSlice=this.nChannels*dims[2]*t + this.nChannels*zz + currentChan + 1;
					data[t][vo][1][t2teInd]=hyperImg.getStack().getProcessor(indexSlice).getPixelValue(xx, yy);
				}
				//System.out.println(TransformUtils.stringVectorN(data[t][vo][1],"T2seq")+"\n");
			}
		}
		return data;
	}*/

	public double[][][][]getFullMRISignalInTheseCoordinatesT1T2(int xCor,int yCor, int zCor,int[][]coordinates,double sigmaXYInVoxels){
		int[][]t1t2Indices=this.getT2Indices();
		int nHits=coordinates.length;
		int bboxX0=100000000;
		int bboxXf=0;
		int bboxY0=100000000;
		int bboxYf=0;
		for(int vo=0;vo<nHits;vo++) {
			if(coordinates[vo][0]<bboxX0)bboxX0=coordinates[vo][0];
			if(coordinates[vo][0]>bboxXf)bboxXf=coordinates[vo][0];
			if(coordinates[vo][1]<bboxY0)bboxY0=coordinates[vo][1];
			if(coordinates[vo][1]>bboxYf)bboxYf=coordinates[vo][1];
		}
		bboxX0=Math.max(0, bboxX0-5);
		bboxXf=Math.min(dims[0]-1, bboxXf+5);
		bboxY0=Math.max(0, bboxY0-5);
		bboxYf=Math.min(dims[1]-1, bboxYf+5);
	
		ImagePlus temp=VitimageUtils.cropMultiChannelFloatImage(hyperImg,bboxX0,bboxXf,bboxY0,bboxYf,zCor,zCor);
		temp=VitimageUtils.gaussianFiltering(temp,sigmaXYInVoxels,sigmaXYInVoxels,0);//It's no error : no "big smoothing" over Z, due to misalignment
		
		double[][][][]data= new double[t1t2Indices.length][nHits][1][];//[times][vox][Seq][echos]		
		int currentChan;
		int index,xx,yy,zz;		
		for(int t=0;t<this.nTimepoints;t++) {			
			for(int vo=0;vo<nHits;vo++) {
				xx=coordinates[vo][0]-bboxX0;
				yy=coordinates[vo][1]-bboxY0;
				zz=0;
				data[t][vo][1]=new double[t1t2Indices[t].length];
				for(int t1trInd=0;t1trInd<t1t2Indices[t].length;t1trInd++) {			
					currentChan=t1t2Indices[t][t1trInd];
					int indexSlice=this.nChannels*1*t + this.nChannels*0 + currentChan + 1;
					data[t][vo][0][t1trInd]=temp.getStack().getProcessor(indexSlice).getPixelValue(xx, yy);
				}
			}
		}
		return data;
	}

	public int[][]getCoordinatesAroundThisVoxel(int x,int y, int z,int crossWidth,int crossThick){
		int xm,ym,xM,yM,zm,zM;
		xm=Math.max(0, x-crossWidth);
		xM=Math.min(dims[0]-1, x+crossWidth);
		ym=Math.max(0, y-crossWidth);
		yM=Math.min(dims[1]-1, y+crossWidth);
		zm=Math.max(0, z-crossThick);
		zM=Math.min(dims[2]-1, z+crossThick);

		int nHits=(xM-xm+1)*(yM-ym+1)*(zM-zm+1);
		int[][]coords= new int[nHits][3];//[vox][dim]
		for(int zz=zm;zz<=zM;zz++) {
			for(int xx=xm;xx<=xM;xx++) {
				for(int yy=ym;yy<=yM;yy++) {
					int index=(xx-xm) + (xM-xm+1)*(yy-ym) + (xM-xm+1)*(yM-ym+1)*(zz-zm);
					coords[index]=new int[] {xx,yy,zz};
				}
			}
		}
		return coords;
	}	
	
	
	
	
}
