package io.github.rocsg.fijirelax.test;

import static org.junit.jupiter.api.Assertions.assertEquals;

import org.junit.jupiter.api.Test;
import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.plugin.Duplicator;
import io.github.rocsg.fijiyama.common.VitimageUtils;
import io.github.rocsg.fijiyama.registration.ItkTransform;
import io.github.rocsg.fijirelax.mrialgo.HyperMap;
import io.github.rocsg.fijirelax.mrialgo.MRDataType;
import io.github.rocsg.fijirelax.mrialgo.MRUtils;
import io.github.rocsg.fijirelax.mrialgo.NoiseManagement;
import io.github.rocsg.fijirelax.mrialgo.RiceEstimator;

import java.io.File;

			
public class Test_FijiRelaxPackage {
	boolean completeTest=true;

	

	public Test_FijiRelaxPackage() {		
	}
	
	
	
	public static void main(String[]args) throws Exception {
		ImageJ ij=new ImageJ();
		Test_FijiRelaxPackage te=new Test_FijiRelaxPackage();		
		te.test_03_RegisterAndComputeMapsCircle();
	}
	
	

	// Test main computation features of the plugin
	//@Test
	public void test_03_RegisterAndComputeMapsCircle()throws Exception{
		System.out.println("\n\n\nRUNNING TEST 03\n\n\n");
		if(!completeTest)return;
		try {ItkTransform tr=new ItkTransform();}
		catch (UnsatisfiedLinkError error) {
			System.out.println(error.getMessage());
			System.out.println("\n\n\nNo test for registration cause dependance is missing\n\n\n");
			System.out.println("LD PATH="+System.getProperty("java.library.path"));
			System.out.println("Content of lib path=");
			File f=new File(System.getProperty("java.library.path"));
			if(f!=null) {
				String[]ss=f.list();
				if(ss!=null)for(String s: f.list() )System.out.println(s);
			}
			System.out.println("PROP="+System.getProperties().toString());
			VitimageUtils.waitFor(5000);
			return;
		}
		
		String mainDir=new File(getClass().getClassLoader().getResource("data/test_1/").getFile()).getAbsolutePath();
		
		//Open the source data
		File f1=new File(mainDir,"SimulatedData_echoes.tif");
		String sourcePathSimulatedEchoes=f1.getAbsolutePath();
		HyperMap hyperFantom=new HyperMap(IJ.openImage(sourcePathSimulatedEchoes)); 

		


 
		//Run of the registration routine
		try{hyperFantom.registerEchoes();}catch(java.awt.HeadlessException he) {he.printStackTrace();return;}

		//Run copy constructor
		HyperMap hyperFantom2=HyperMap.hyperMapFactory(hyperFantom);

		//Run of maps computation, with default arguments
		hyperFantom2.computeMaps();
		
		
		//Run the ranging routine
		hyperFantom2.setDisplayRange();
		//Test of the obtained map range
		String sComp=hyperFantom2.getRangeInfos();
		double maxComp=VitimageUtils.maxOfImage(new Duplicator().run(hyperFantom2.getMapsImage(),1,1,1,1,1,1));
		double maxExp=3000.0021972656;
		assertEquals(maxComp,maxExp,VitimageUtils.EPSILON);
		
		//Test if the computed map is equal to the expected one
		ImagePlus computedMaps=hyperFantom2.getMapsImage();		
		File f2=new File(mainDir,"SimulatedData_expected_maps.tif");
		String sourcePathExpectedMaps=f2.getAbsolutePath();
		ImagePlus expectedMaps=IJ.openImage(sourcePathExpectedMaps);
		ImagePlus diff=VitimageUtils.substract(computedMaps,expectedMaps, true);
		double lowerBound=VitimageUtils.minOfImage(diff);
		double upperBound=VitimageUtils.maxOfImage(diff);		
		assertEquals(lowerBound, 0,VitimageUtils.EPSILON);
		assertEquals(upperBound, 0,VitimageUtils.EPSILON);

		int[]tabAlg=new int[] {MRUtils.SIMPLEX/*,MRUtils.LM*/};
		NoiseManagement[]tabNoise=new NoiseManagement[] {NoiseManagement.RICE,NoiseManagement.NOTHING,NoiseManagement.OFFSET};
		boolean []tabForget=new boolean[] {true,false};
		int i=0;
		for(int alg:tabAlg) {
			for(NoiseManagement noi:tabNoise) {
				for(boolean forg:tabForget) {
					System.out.println("Test "+(i++)+"/12");
					HyperMap hyperFantom3=HyperMap.hyperMapFactory(hyperFantom);
					hyperFantom3.computeMapsAgainAndMask(alg,false,noi,forg,null,4,false,0);
				}
			}
		}
	}
	
	
	// Test main computation features of the plugin
	//@Test
	public void test_04_all_fits() {
		System.out.println("\n\n\nRUNNING TEST 04\n\n\n");
		if(!completeTest)return;
		ImagePlus t1map=IJ.createImage("T1MAP", 100, 100, 1, 32);
		VitimageUtils.set32bitToValue(t1map, 1600, 0);
		ImagePlus t2map=IJ.createImage("T2MAP", 100, 100, 1, 32);
		VitimageUtils.set32bitToValue(t1map, 3000, 0);
		ImagePlus pdmap=IJ.createImage("PDMAP", 100, 100, 1, 32);
		VitimageUtils.set32bitToValue(t1map, 5000, 0);
		ImagePlus maskmap=IJ.createImage("MASKMAP", 100, 100, 1, 32);
		VitimageUtils.set32bitToValue(t1map, 255, 0);
		ImagePlus testMaps;
		ImagePlus testEchoes;
		
		testMaps=VitimageUtils.hyperStackingChannels(new ImagePlus[] {pdmap,t1map,maskmap});
		try {
		testEchoes=MRUtils.simulateEchoesT1Relax(
				testMaps, new double[] {600,1200,2400,4000}, 11, 12,MRDataType.T1SEQ,100,"TestT1" 
		);
		
		testMaps=VitimageUtils.hyperStackingChannels(new ImagePlus[] {pdmap,t2map,maskmap});
		testEchoes=MRUtils.simulateEchoesT2Relax(
				testMaps, 2400,  new double[] {11,11,11,11}, 12,MRDataType.T1SEQ,100,"TestT2" 
		);

		testMaps=VitimageUtils.hyperStackingChannels(new ImagePlus[] {pdmap,t1map,t2map,maskmap});
		testEchoes=MRUtils.simulateEchoesT1T2Relax(
				testMaps, new double[] {600,1200,2400,4000,4000,4000,4000,4000,4000,4000},  new double[] {11,11,11,11,22,33,44,55,66,77}, 12,new MRDataType[] {MRDataType.T1SEQ,MRDataType.T2SEQ,MRDataType.T1T2SEQ},100,"TestT1T2" 
		);
		}catch(java.awt.HeadlessException he) {System.out.println("A Headless exception occured. Test early break");return;}
	}

	
	// Test the rice estimator
	//@Test
	public void test_05_rice() {
		if(!completeTest)return;
		RiceEstimator.testRandomRiceMeanDrift();
		RiceEstimator.testCorruptAndRecoverSignal();	
		RiceEstimator.testSimulateAndEstimateRiceSigma();
	}


	//@Test
	public void test_06_Importers() throws Exception{
		System.out.println("\n\n\nRUNNING TEST 06\n\n\n");
		if(!completeTest)return;
		String path=new File(getClass().getClassLoader().getResource("data/test_2/").getFile()).getAbsolutePath();
		ImagePlus imgMap=IJ.openImage(
				new File(getClass().getClassLoader().getResource("data/test_2/TestDataHyperMap/smallSorghoSerieRaw.tif").getFile()).getAbsolutePath());

		HyperMap hyp=HyperMap.importHyperMapFromCustomData(new File(path,"TestDataImportRaw").getAbsolutePath(),"img_TR{TR}_TE{TE}.tif");
		ImagePlus resultCustom=hyp.getEchoesImage();
		//IJ.saveAsTiff(resultCustom, new File(pathTemp,"testCust.tif").getAbsolutePath());
		
		hyp=HyperMap.importHyperMapFromRawDicomData(new File(path,"TestDataImportDicom").getAbsolutePath(),"testDicom");
		ImagePlus resultDicom=hyp.getEchoesImage();
		resultDicom=new Duplicator().run(resultDicom,1,resultDicom.getNChannels(),1,1,1,1);
		//IJ.saveAsTiff(resultDicom, new File(pathTemp,"testDicom.tif").getAbsolutePath());
		
		ImagePlus resultExpected=IJ.openImage(new File(path,"TestDataHyperMap/smallSorghoSerieRaw.tif").getAbsolutePath());
		resultExpected=new Duplicator().run(resultExpected,1,resultExpected.getNChannels(),1,1,1,1);
		//IJ.saveAsTiff(resultExpected, new File(pathTemp,"testExp.tif").getAbsolutePath());

		
		ImagePlus diff=VitimageUtils.substract(resultExpected,resultCustom, true);
		double lowerBound=VitimageUtils.minOfImage(diff);
		double upperBound=VitimageUtils.maxOfImage(diff);		
		assertEquals(lowerBound, 0,VitimageUtils.EPSILON);
		assertEquals(upperBound, 0,VitimageUtils.EPSILON);
		//IJ.saveAsTiff(diff, "/home/rfernandez/Bureau/testDiff1.tif");
		
		
		diff=VitimageUtils.substract(resultExpected,resultDicom, true);
		//IJ.saveAsTiff(diff, "/home/rfernandez/Bureau/testDiff2.tif");
		lowerBound=VitimageUtils.minOfImage(diff);
		upperBound=VitimageUtils.maxOfImage(diff);		
		assertEquals(lowerBound, 0,VitimageUtils.EPSILON);
		assertEquals(upperBound, 0,VitimageUtils.EPSILON);
	}

	//@Test
	public void test_07_Estimators() throws Exception{
		System.out.println("\n\n\nRUNNING TEST 07\n\n\n");
		if(!completeTest)return;
		ImagePlus imgMap=IJ.openImage(
				new File(getClass().getClassLoader().getResource("data/test_2/TestDataHyperMap/smallSorghoSerieRaw.tif").getFile()).getAbsolutePath());

		
		HyperMap map=new HyperMap(imgMap);
		
		int[]algTypes=new int[] {MRUtils.SIMPLEX,MRUtils.LM};
		boolean[]separateds=new boolean[] {false,true};
		NoiseManagement[]noises=new NoiseManagement[] {NoiseManagement.OFFSET,NoiseManagement.NOTHING,NoiseManagement.RICE};
		boolean[]forgets=new boolean[] {false,true};
		double[]nbStds=new double[] {3,5};
		int incr=0;
		for(int alg:algTypes)for(boolean separated:separateds)for(NoiseManagement noi: noises)for(boolean forget:forgets)for(double std:nbStds) {
			IJ.log(""+(incr++)+" / 48");
			map.computeMapsAgainAndMask(alg, separated, noi, forget, null, std,false,0);
		}
	}

	//@Test
	public void test_08_Outliers() throws Exception{
		System.out.println("\n\n\nRUNNING TEST 08\n\n\n");
		if(!completeTest)return;
		ImagePlus imgMap=IJ.openImage(
				new File(getClass().getClassLoader().getResource("data/test_2/TestDataHyperMap/smallSorghoSerieRaw.tif").getFile()).getAbsolutePath());
		HyperMap map=new HyperMap(imgMap);
		map.computeMaps();
		map.simulateOutlierRemoval(10, 10, 20, 20, 0, 0);
	
	}


}
