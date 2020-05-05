package com.vitimage.mrutils;

import java.util.ArrayList;
import com.vitimage.aplimtools.T1T2Seq_Importer;
import com.vitimage.common.TransformUtils;
import com.vitimage.common.VitiDialogs;
import com.vitimage.common.VitimageUtils;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.plugin.Concatenator;
import ij.plugin.Duplicator;
import ij.plugin.HyperStackConverter;
import ij.plugin.filter.Convolver;
import ij.process.ImageProcessor;


public class HyperMRIT1T2 {
	private double sigmaInUse=0.75;
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
					for(int ind=getNumberOfMaps();ind<nChannels;ind++)average+=MRUtils.readSigmaInSliceLabel(hyperImg, ind, z, t);
					this.tabSigmasT1T2Seq[t][z]=average/(nChannels-getNumberOfMaps());
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
	
	
	/** Capillary utilities. The first normalize MR signal along Z and T axis, using the capillary values*/
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

	public static ImagePlus getAntiCapillaryMask(ImagePlus img) {
		//Find capillary in the central slice
		ImagePlus imgSlice=new Duplicator().run(img,1,1,img.getNSlices()/2+1,img.getNSlices()/2+1,1,1);
		int []coordsCentral=findCapillaryCenterInSlice(imgSlice,bionanoCapillaryRadius);
		int[][]capillaryCenters=capillaryCoordinatesAlongZ(img,coordsCentral,bionanoCapillaryRadius);

		ImagePlus ret=new Duplicator().run(img,1,1,1,img.getNSlices(),1,1);
		IJ.run(ret,"32-bit","");
		ret=VitimageUtils.makeOperationOnOneImage(ret, 2, 0, true);
		int radiusCapLarge=(int)Math.round(bionanoCapillaryRadius*2.0/VitimageUtils.getVoxelSizes(img)[0]);
		int radiusSquare=radiusCapLarge*radiusCapLarge;
		int radiusCapShort=(int)Math.round(bionanoCapillaryRadius*0.7/VitimageUtils.getVoxelSizes(img)[0]);
		int radiusSquareShort=radiusCapShort*radiusCapShort;
		
		int dimX=ret.getWidth();
		int dimY=ret.getHeight();
		for(int z=0;z<img.getNSlices();z++) {
			System.out.println("Pour z="+z+" found "+capillaryCenters[z][0]+" , "+capillaryCenters[z][1]);
			float[]retVal=(float[])ret.getStack().getProcessor(z+1).getPixels();
			int xm=Math.max(0,capillaryCenters[z][0]-radiusCapLarge);
			int ym=Math.max(0,capillaryCenters[z][1]-radiusCapLarge);
			int xM=Math.min(dimX-1,capillaryCenters[z][0]+radiusCapLarge);
			int yM=Math.min(dimY-1,capillaryCenters[z][1]+radiusCapLarge);
			for(int x=xm;x<=xM;x++)			for(int y=ym;y<=yM;y++) {
				if( (x-capillaryCenters[z][0])*(x-capillaryCenters[z][0]) + (y-capillaryCenters[z][1]) * (y-capillaryCenters[z][1]) < radiusSquareShort )retVal[dimX*y+x]=4;
				else if( (x-capillaryCenters[z][0])*(x-capillaryCenters[z][0]) + (y-capillaryCenters[z][1]) * (y-capillaryCenters[z][1]) < radiusSquare )retVal[dimX*y+x]=2;
				else retVal[dimX*y+x]=0;
			}
		}
		return ret;
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

	public static int[][]capillaryCoordinatesAlongZ(ImagePlus img3D,int[]coords,double capillaryRadius){
		int[][]ret=new int[img3D.getNSlices()][2];
		int[]dims=VitimageUtils.getDimensions(img3D);
		double[]voxs=VitimageUtils.getVoxelSizes(img3D);
		
		Duplicator dup=new Duplicator();
		int rayPix=(int)Math.round(capillaryRadius/(voxs[0]));
		int semiRayPix=rayPix/2;
		int lookupRadius=(int)Math.round(1.5*rayPix);

		int xMed=coords[0];int yMed=coords[1];int zMed=dims[2]/2;
		double[]capVals=new double[dims[2]];
		ImagePlus imgTemp;
		int xLast=xMed,yLast=yMed;
		long t0= System.currentTimeMillis();

		for(int z=zMed;z>=0;z--) {
			//Extract a patch around the finding of upside
			imgTemp=VitimageUtils.cropImageFloat(img3D,xLast-lookupRadius,yLast-lookupRadius,z, lookupRadius*2,lookupRadius*2,1);
			//Find cap center in it
			int[]coordsNew=findCapillaryCenterInSlice(imgTemp,bionanoCapillaryRadius);
			
			//Update coordinates of last
			xLast=xLast-lookupRadius+coordsNew[0];
			yLast=yLast-lookupRadius+coordsNew[1];

			ret[z][0]=xLast;
			ret[z][1]=yLast;
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
			
			ret[z][0]=xLast;
			ret[z][1]=yLast;
		}
		return ret;
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


	
	
	
	
	/** Access to maps, echoes and Tr/Te values*/
	public Object[] getMapsImage(boolean writeParametersForExplorer) {
		int fontSize=15;
		int nMaps=getNumberOfMaps();
		ImagePlus maps= new Duplicator().run(hyperImg,1,nMaps,1,dims[2],1,nTimepoints);
		String[][][]texts=new String[nMaps][dims[2]][nTimepoints];
		for(int z=0;z<dims[2];z++)for(int t=0;t<nTimepoints;t++)for(int c=0;c<nMaps;c++) {
			texts[c][z][t]="   "+MRUtils.getDataTypeOfThisMagneticResonanceSlice(maps, c, z, t)+"  Day="+actualDay[t]+"    Z="+z;
			if(writeParametersForExplorer) VitimageUtils.writeTextOnGivenImage(texts[c][z][t],maps ,fontSize ,30 ,dims[1]-30,c, z, t);
		}  
		maps.setPosition(VitimageUtils.getCorrespondingSliceInHyperImage(maps, 0, dims[2]/2, 0));
		return new Object[] {maps,texts};
	}

	public int getNumberOfMaps() {
		int nMaps=0;
		for(int i=0;i<hyperImg.getNChannels();i++) {
			if(hyperImg.getStack().getSliceLabel(VitimageUtils.getCorrespondingSliceInHyperImage(hyperImg, i, 0, 0)).contains("MAP")) {
				nMaps=i+1;		
			}
		}
		return nMaps;
	}

	public Object[] getEchoesImage(boolean writeParametersForExplorer) {
		int fontSize=15;
		int nMaps=getNumberOfMaps();
		System.out.println("Nmaps="+nMaps);
		System.out.println("new Duplicator().run(hyperImg,"+(nMaps+1)+","+nChannels+",1,"+dims[2]+",1,"+nTimepoints);
		VitimageUtils.printImageResume(hyperImg);
		ImagePlus echoes= new Duplicator().run(hyperImg,nMaps+1,nChannels,1,dims[2],1,nTimepoints);
		String[][][]texts=new String[nChannels-nMaps][dims[2]][nTimepoints];
		for(int z=0;z<dims[2];z++)for(int t=0;t<nTimepoints;t++) for(int c=0;c<nChannels-nMaps;c++) {
			texts[c][z][t]="   "+MRUtils.getDataTypeOfThisMagneticResonanceSlice(echoes, c, z, t)+"  Day="+actualDay[t]+"   Z="+z+
						"  TR="+(MRUtils.readTrInSliceLabel(echoes, c, z, t))+"  TE="+(MRUtils.readTeInSliceLabel(echoes, c, z, t));
			if(writeParametersForExplorer) VitimageUtils.writeTextOnGivenImage(texts[c][z][t],echoes ,fontSize ,30 ,10,c , z, t);
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
		for(int c=0;c<this.nChannels;c++) {
			System.out.print("Inspecting channel t"+t+" c"+c+" with "+this.Tr[t][0][c]+" , "+this.Te[t][0][c]+"...");
			if(mrDataType[t][c]==MRDataType.T1T2SEQ) {
				System.out.println("Ok.");
				for(int z=0;z<dims[2];z++)
					ret[z][incr]=new double[] {this.Tr[t][z][c],this.Te[t][z][c]};
				incr++;
			}
			else System.out.println("No.");
		}
		return ret;
	}
		
	public double[][][][]getT1T2TrTeTimes(){
		double[][][][]ret=new double[nTimepoints][][][];
		for(int t=0;t<this.nTimepoints;t++)ret[t]=getT1T2TrTeTimes(t);
		return ret;
	}
	
	
	/**Access to magnitude data*/
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
							//System.out.println("data["+t+"]["+index+"][0]["+t1trInd+"]="+data[t][index][0][t1trInd]);
						}
					}
				}	
			}
		}
		return data;
	}

	public double[][][][]getFullMRISignalInTheseCoordinatesT1T2(int xCor,int yCor, int zCor,int[][]coordinates,double sigmaXYInVoxels){
		int[][]t1t2Indices=this.getT1T2Indices();
		//System.out.println("GOT T1T2"+t1t2TrTeIndices[0][0].length+" , "+t1t2TrTeIndices[0][0][0].length);
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
		temp=VitimageUtils.gaussianFilteringMultiChannel(temp,sigmaXYInVoxels,sigmaXYInVoxels,0);//It's no error : no "big smootphrahing" over Z, due to misalignment
		double[][][][]data= new double[t1t2Indices.length][nHits][1][];//[times][vox][Seq][echos]		
		int currentChan;
		int index,xx,yy;	
		for(int t=0;t<this.nTimepoints;t++) {			
			for(int vo=0;vo<nHits;vo++) {
				xx=coordinates[vo][0]-bboxX0;
				yy=coordinates[vo][1]-bboxY0;
				System.out.println("\nProcessing voxel "+vo+" at coordinates "+xx+" , "+yy);
				data[t][vo][0]=new double[t1t2Indices[t].length];
				for(int t1trInd=0;t1trInd<t1t2Indices[t].length;t1trInd++) {			
					currentChan=t1t2Indices[t][t1trInd];
					int indexSlice=this.nChannels*t + currentChan + 1;
					data[t][vo][0][t1trInd]=temp.getStack().getProcessor(indexSlice).getPixelValue(xx, yy);
				}
				System.out.println(" ->"+TransformUtils.stringVectorN(data[t][vo][0], ""));
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
	
	
	
	
	
	
	public static HyperMRIT1T2 importHyperMapFromRawDataDir(String inputDir,String nameObservation) {
		T1T2Seq_Importer t1t2=new T1T2Seq_Importer();
		ImagePlus hyperMap=t1t2.run(inputDir,nameObservation);
		return new HyperMRIT1T2(hyperMap);
	}
	
	
	public void computeMapsAgain() {
		hyperImg.show();
		VitimageUtils.waitFor(10000);
		if(hyperImg.getNFrames()>1) {VitiDialogs.notYet("ComputeMapsAgain of 5D hypermaps in HyperMRIT1T2");return;}
		int Z=hyperImg.getNSlices();
		int C=hyperImg.getNChannels();
		int T=hyperImg.getNFrames();
		int X=hyperImg.getWidth();
		int Y=hyperImg.getHeight();
		ImagePlus imgMask=new Duplicator().run(hyperImg,4,4,1,Z,1,T);
		ImagePlus []imgT1T2Line=VitimageUtils.stacksFromHyperstackFastBis((ImagePlus)(getEchoesImage(false)[0]));
		boolean makeBoutureTrick=hyperImg.getStack().getSliceLabel(1).contains("BOUTURE");
		IJ.run(imgMask,"32-bit","");
		for(int i=0;i<imgT1T2Line.length;i++)IJ.run(imgT1T2Line[i],"32-bit","");
		
		//Compute T1T2mono
		ImagePlus[]	maps=MRUtils.computeT1T2Map(imgT1T2Line, sigmaInUse,makeBoutureTrick ? MRUtils.T1T2_MONO_RICE : MRUtils.T1T2_MONO_RICE);		
		ImagePlus M0Mono=maps[0];
		ImagePlus T1Mono=maps[0];
		ImagePlus T2Mono=maps[0];
		ImagePlus khi2Mono=maps[3];

		
		//Compute T1T2bi
		maps=MRUtils.computeT1T2Map(imgT1T2Line, sigmaInUse,makeBoutureTrick ? MRUtils.T1T2_MULTI_RICE : MRUtils.T1T2_MULTI_RICE);
		ImagePlus M0Multi1=maps[0];
		ImagePlus M0Multi2=maps[1];
		ImagePlus T1Multi=maps[2];
		ImagePlus T2Multi1=maps[3];
		ImagePlus T2Multi2=maps[4];
		ImagePlus khi2Multi=maps[5];

		
		//Combine
		ImagePlus retM0=M0Mono.duplicate();
		ImagePlus retT2=T2Mono.duplicate();
		ImagePlus retT1=T1Mono.duplicate();
		ImagePlus selectedModel=imgMask.duplicate();
		ImagePlus portionMap=imgMask.duplicate();
		int index;
		for(int z=0;z<dims[2];z++) {
			float[]valsRetM0=(float[]) retM0.getStack().getProcessor(z+1).getPixels();
			float[]valsRetT2=(float[]) retT2.getStack().getProcessor(z+1).getPixels();
			float[]valsRetT1=(float[]) retT1.getStack().getProcessor(z+1).getPixels();
			float[]valsM0=(float[]) M0Mono.getStack().getProcessor(z+1).getPixels();
			float[]valsM0Multi1=(float[]) M0Multi1.getStack().getProcessor(z+1).getPixels();
			float[]valsM0Multi2=(float[]) M0Multi2.getStack().getProcessor(z+1).getPixels();
			float[]valsT1=(float[]) T1Mono.getStack().getProcessor(z+1).getPixels();
			float[]valsT1Multi=(float[]) T1Multi.getStack().getProcessor(z+1).getPixels();
			float[]valsT2=(float[]) T2Mono.getStack().getProcessor(z+1).getPixels();
			float[]valsT2Multi1=(float[]) T2Multi1.getStack().getProcessor(z+1).getPixels();
			float[]valsT2Multi2=(float[]) T2Multi2.getStack().getProcessor(z+1).getPixels();
			float[]valMask=(float[]) imgMask.getStack().getProcessor(z+1).getPixels();
			float[]valKhi2=(float[]) khi2Mono.getStack().getProcessor(z+1).getPixels();
			float[]valKhi2Multi=(float[]) khi2Multi.getStack().getProcessor(z+1).getPixels();
			float[]valSelectedModel=(float[]) selectedModel.getStack().getProcessor(z+1).getPixels();
			float[]valPortionMap=(float[]) portionMap.getStack().getProcessor(z+1).getPixels();
			for(int y=0;y<dims[1];y++) {
				for(int x=0;x<dims[0];x++) {
					index=y*dims[0]+x;
					if((int)Math.round(valMask[index])%2==0) {
						valsRetM0[index]=0;
						valsRetT1[index]=0;
						valsRetT2[index]=0;
						continue;//If out of mask, set to 0  (work already done)
					}
					if(valKhi2Multi[index]>=MRUtils.ERROR_KHI2) {
						if(valKhi2[index]<MRUtils.ERROR_KHI2) {
							//if multicomp is out, and monocomp khi < errorkhi, set to monocomp (already done)
						}
						else {
							//else set to 0
							valsRetM0[index]=0;
							valsRetT1[index]=0;
							valsRetT2[index]=0;
						}
					}
					else {
						if(valKhi2Multi[index]<valKhi2[index]) {
							//if T2 bi less Khi, // eventually : and T2 vals > 10  , < 500 and ratio between components > 0.02, < 50
							valsRetM0[index]=valsM0Multi1[index]+valsM0Multi2[index];
							valsRetT1[index]=valsT1Multi[index];
							valsRetT2[index]=(valsT2Multi1[index]*valsM0Multi1[index]+valsT2Multi2[index]*valsM0Multi2[index])/(valsM0Multi1[index]+valsM0Multi2[index]);
							valSelectedModel[index]=2;
							valPortionMap[index]=valsM0Multi1[index]/(valsM0Multi1[index]+valsM0Multi2[index]);
						}
						else { 
							if(valKhi2[index]<MRUtils.ERROR_KHI2) {
								//if T1 khi < errorkhi, set to T1
								valsRetM0[index]=valsM0[index];
								valsRetT1[index]=valsT1[index];
								valsRetT2[index]=valsT2[index];
							}
							else {
								//else set to 0
								valsRetM0[index]=0;
								valsRetT1[index]=0;
								valsRetT2[index]=0;
							}
						}
					}
				}
			}
		}
		
		ImagePlus []tempRes=new ImagePlus[C];
		tempRes[0]=retM0;
		tempRes[1]=retT1;
		tempRes[2]=retT2;
		tempRes[3]=imgMask;
		for(int c=4;c<C;c++) {tempRes[c]=imgT1T2Line[c-4];}
		for(int c=0;c<C;c++)tempRes[c]=VitimageUtils.convertFloatToShortWithoutDynamicChanges(tempRes[c]);			

		ImagePlus hyperImg=Concatenator.run(tempRes);
		hyperImg=HyperStackConverter.toHyperStack(hyperImg, C,Z,T,"xyztc","Grayscale");

		
		//Update capillaryvalues
		short[]capVals;
		int nHits;
		int incr;
		int slice=0;
		double sum;
		for(int t=0;t<T;t++) {
			for(int z=0;z<Z;z++) {
				ArrayList<int[]> coordsPointsCap=new ArrayList<int[]>();
				capVals=(short[])imgMask.getStack().getProcessor(z+1).getPixels();
				nHits=0;
				for(int x=0;x<X;x++)for(int y=0;y<Y;y++) if(( capVals[VitimageUtils.getPixelIndex(x, X, y)] & 0xffff ) >4.5) {nHits++; coordsPointsCap.add(new int[] {x,y});};	
				for(int c=0;c<C;c++) {
					double[]vals=new double[nHits];
					incr=0;
					slice=VitimageUtils.getCorrespondingSliceInHyperImage(hyperImg, c, z, t);
					capVals=(short[])hyperImg.getStack().getProcessor(slice).getPixels();				
					for(int i=0;i<nHits;i++) {
						vals[incr++]=(int)((capVals[VitimageUtils.getPixelIndex(coordsPointsCap.get(i)[0], X, coordsPointsCap.get(i)[1])]) & 0xffff);
					}
					double[]stats=VitimageUtils.statistics1D(vals);
					String label=hyperImg.getStack().getSliceLabel(slice);
					String []infos=label.split("_");
					label="";
					for(int i=0;i<infos.length;i++) {
						if(infos[i].contains("CAP="))infos[i]=("CAP="+VitimageUtils.dou(stats[0])+"+-"+VitimageUtils.dou(stats[1]));
						label+=(infos[i])+( (i==infos.length-1) ?  "" : "_");
					}
					hyperImg.getStack().setSliceLabel(label,slice);
				}
			}
		}
		
		
		//Update ranges and luts
		hyperImg.setC(1);
		IJ.run(hyperImg,"Fire","");
		hyperImg.setDisplayRange(0,MRUtils.maxDisplayedM0 );

		hyperImg.setC(2);
		IJ.run(hyperImg,"Fire","");
		hyperImg.setDisplayRange(0,MRUtils.maxDisplayedT1 );

		hyperImg.setC(3);
		IJ.run(hyperImg,"Fire","");
		hyperImg.setDisplayRange(0,MRUtils.maxDisplayedT2 );

		hyperImg.setC(4);
		IJ.run(hyperImg,"Fire","");
		hyperImg.setDisplayRange(0,1 );

		hyperImg.setC(5);
		IJ.run(hyperImg,"Fire","");
		hyperImg.setDisplayRange(-1,1 );

		for(int c=6;c<=C;c++) {
			hyperImg.setC(c);
			IJ.run(hyperImg,"Fire","");
			hyperImg.setDisplayRange(0,MRUtils.maxDisplayedM0 );
		}

		hyperImg.show();

		
	}
	

	
	
	
	
	
	/**Random helpers*/
	public void showCopy(String text) {
		ImagePlus temp=hyperImg.duplicate();
		temp.show();
		temp.setTitle(text);
	}
	
	public ImagePlus getCopy() {
		return hyperImg.duplicate();
	}
	
	
	
	
	
	
	

}
