package com.vitimage.mrutils;

import java.io.File;

import com.vitimage.common.VitimageUtils;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.plugin.Concatenator;
import ij.plugin.Duplicator;
import ij.plugin.HyperStackConverter;

public class TestingMRUtilsPackage {

	//ajout
	//ajout2
	//Test hotfix
	public static void main(String[]args) {
		ImageJ ij=new ImageJ();
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
