package com.vitimage.mrutils;

import java.io.File;

import com.vitimage.aplimtools.MRI_HyperCurvesExplorer;
import com.vitimage.common.VitimageUtils;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.plugin.Concatenator;
import ij.plugin.Duplicator;
import ij.plugin.HyperStackConverter;

public class TestingMRUtilsPackage {


	public static void testDims() {
		String pathSource="/home/fernandr/Bureau/Traitements/Sorgho/Donnees_brutes_export_Romain";
		String tempStr=pathSource;
		for (String str1 : VitimageUtils.stringArraySort( new File(tempStr).list())) {
			System.out.println("\n\n\nSpecimen="+str1);
			tempStr=pathSource+"/"+str1;
			for (String str2 : VitimageUtils.stringArraySort( new File(tempStr).list())) {
				System.out.println("|\n|--Experiment="+str2);
				tempStr=pathSource+"/"+str1+"/"+str2;
				for (String str3 : VitimageUtils.stringArraySort( new File(tempStr).list())) {
					tempStr=pathSource+"/"+str1+"/"+str2+"/"+str3;
					String str4 = VitimageUtils.stringArraySort( new File(tempStr).list())[0];
					tempStr=pathSource+"/"+str1+"/"+str2+"/"+str3+"/"+str4;
					String str5 = VitimageUtils.stringArraySort( new File(tempStr).list())[0];
					tempStr=pathSource+"/"+str1+"/"+str2+"/"+str3+"/"+str4+"/"+str5;
//					VitimageUtils.printImageResume(IJ.openImage(tempStr));
					double []dims=VitimageUtils.getVoxelSizes(IJ.openImage(tempStr));
					double vol=VitimageUtils.getVoxelVolume(IJ.openImage(tempStr));
					System.out.println("|----"+str2+" "+str3+" "+str4+" VoxelVolume="+VitimageUtils.dou(1000*vol)+" .10^-3 mm3     "+"Voxelsizes = "+dims[0]+" X "+dims[1]+" X "+dims[1]);
				}	
			}			
		}
	}
	
	
	public static void testNorm() {
		
		ImagePlus[]tab=new ImagePlus[4];
		tab[0]=IJ.openImage("/home/fernandr/Bureau/Temp/m0.tif");
		HyperMRIT1T2.measureMeanCapillaryValueAlongZ(tab[0]);
		VitimageUtils.waitFor(1000000);
		tab[1]=IJ.openImage("/home/fernandr/Bureau/Temp/img22.tif");
		tab[2]=IJ.openImage("/home/fernandr/Bureau/Temp/img33.tif");
		tab[3]=IJ.openImage("/home/fernandr/Bureau/Temp/img44.tif");
		VitimageUtils.waitFor(10000000);
	}
	
	public static void chiasse() {
		String dirIn="/home/fernandr/Bureau/Traitements/Sorgho/Cartes_calculees_methode_Romain";
		boolean makeShort=false;
		String dirOut="";
		dirOut="/home/fernandr/Bureau/Traitements/Sorgho/Cartes_rangees_par_specimen";
		if(makeShort)dirOut="/home/fernandr/Bureau/Traitements/Sorgho/Cartes_short";
		String[]specs=new File(dirIn).list();
		for(String sp : specs) {
			System.out.println("Processing "+sp);
			String[]infoTab=sp.split("_");
			String strNew=infoTab[0]+"_"+infoTab[2]+".tif";
			ImagePlus img=IJ.openImage(new File(dirIn,sp).getAbsolutePath());
			img.setC(1);img.setDisplayRange(0, MRUtils.maxDisplayedM0);
			img.setC(2);img.setDisplayRange(0, MRUtils.maxDisplayedT1);
			img.setC(3);img.setDisplayRange(0, MRUtils.maxDisplayedT2);
			img.setC(4);img.setDisplayRange(0, 1);
			img.setC(5);img.setDisplayRange(-1, 1);
			for(int c=6;c<=img.getNChannels();c++) {img.setC(c);img.setDisplayRange(0, MRUtils.maxDisplayedM0);}
			
			String fileOut=new File(dirOut,infoTab[0]).getAbsolutePath();
			fileOut=new File(fileOut,strNew).getAbsolutePath();
			System.out.println("Sauvegarde in "+fileOut);
			if(makeShort)img=new Duplicator().run(img,1,3,1,4,1,1);
			IJ.saveAsTiff(img,fileOut);			
		}
		System.exit(0);
	}
	
	
	public static void testKhi2() {
		int[][]vals=new int[][] {
			{5,10 },
			{6,12 },
			{10,20 },
			{10,10 },
			{20,10 },
			{50,10 },
			{100,10 },
			{500,10 }
		};
		for(int i=0;i<vals.length;i++) {		System.out.println(" Khi="+vals[i][0] + " Nprms="+vals[i][1]+"  pval="+MRUtils.getPvalue(vals[i][0],vals[i][1]) );}
	}

	
	public static void main(String[]args) {
		ImageJ ij=new ImageJ();
		testDims();
		System.exit(0);
		//testKhi2();
		//testNorm();
		chiasse();
		String path="/home/fernandr/Bureau/Test/Ghetto/fdfddfd/";
		ImagePlus imgT1=IJ.openImage(path+"T1short.tif");
		ImagePlus imgT2=IJ.openImage(path+"T2short.tif");
		ImagePlus imgMaps=IJ.openImage(path+"Maps.tif");
		int[]tr=new int[] {600,1200,2400};
		for(int t=1;t<=imgT1.getNFrames();t++) {
			for(int z=1;z<=imgT1.getNSlices();z++) {
				System.out.println("t="+t+" , z="+z);
				for(int c=1;c<=imgMaps.getNChannels();c++)imgMaps.getStack().setSliceLabel("_SIGMARICE=130_MAPS_TR="+tr[c-1]+"_TE=11", VitimageUtils.getCorrespondingSliceInHyperImage(imgMaps, c-1, z-1, t-1));
				for(int c=1;c<=imgT1.getNChannels();c++)imgT1.getStack().setSliceLabel("_SIGMARICE=130_T1SEQ_TR="+tr[c-1]+"_TE=11", VitimageUtils.getCorrespondingSliceInHyperImage(imgT1, c-1, z-1, t-1));
				for(int c=1;c<=imgT2.getNChannels();c++)imgT2.getStack().setSliceLabel("_SIGMARICE=130_T2SEQ_TR=10000_TE="+(11*c), VitimageUtils.getCorrespondingSliceInHyperImage(imgT2, c-1, z-1, t-1));
			}
		}
		ImagePlus []tabT1=stacksFromHyperstackFastBis(imgT1);
		ImagePlus []tabT2=stacksFromHyperstackFastBis(imgT2);
		ImagePlus []tabMaps=stacksFromHyperstackFastBis(imgMaps);
		ImagePlus []imgTab=new ImagePlus[3*5+3*5+16*5];
		for(int i=0;i<3*5;i++)imgTab[i]=tabMaps[i];
		for(int i=0;i<3*5;i++)imgTab[3*5+i]=tabT1[i];
		for(int i=0;i<16*5;i++)imgTab[3*5+3*5+i]=tabT2[i];
		ImagePlus res=Concatenator.run(imgTab);
		VitimageUtils.printImageResume(res);
		System.out.println("22 , " +imgT1.getNSlices()+" , "+imgT1.getNFrames());
		ImagePlus res2=HyperStackConverter.toHyperStack(res, 22,imgT1.getNSlices(),imgT1.getNFrames(),"xyztc","Fire");
		IJ.run(res2,"32-bit","");
		res2.show();
	}
	public static ImagePlus[]stacksFromHyperstackFastBis(ImagePlus hyper){
		int nbZ=hyper.getNSlices();
		int nbT=hyper.getNFrames();
		int nbC=hyper.getNChannels();
		int nb=nbT*nbC;
		ImagePlus []ret=new ImagePlus[nb];
		for(int ic=0;ic<nbC;ic++) {
			for(int it=0;it<nbT;it++) {
				int i=ic*nbT+it;
				System.out.println(ic+"/"+nbC+" ,  "+it+"/"+nbT);
				ret[i] = new Duplicator().run(hyper, 1+ic, 1+ic, 1, nbZ, 1+it, 1+it);
				VitimageUtils.adjustImageCalibration(ret[i],hyper);
				IJ.run(ret[i],"Grays","");
			}
		}
		return ret;
	}

	public TestingMRUtilsPackage() {
		MRI_HyperCurvesExplorer explorer=new MRI_HyperCurvesExplorer();
	}

}
