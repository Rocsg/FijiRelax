package io.github.rocsg.fijirelax.test;

import static org.junit.jupiter.api.Assertions.assertEquals;

import org.junit.Before;
import org.junit.jupiter.api.Test;
import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.plugin.Duplicator;
import io.github.rocsg.fijiyama.common.VitimageUtils;
import io.github.rocsg.fijirelax.mrialgo.HyperMap;
import io.github.rocsg.fijiyama.fijiyamaplugin.RegistrationAction;
import ij.IJ;
import ij.ImagePlus;
import java.io.File;


public class Test_FijiRelaxPackage {
	@Before
	public void initialize() {
		ImageJ ij=new ImageJ();
	}	

	public Test_FijiRelaxPackage() {		
	}
	
	
	
	public static void main(String[]args) throws Exception {
		ImageJ ij=new ImageJ();
		Test_FijiRelaxPackage te=new Test_FijiRelaxPackage();		
		te.test_03_RegisterAndComputeMapsCircle();
	}
	
	@Test
	public void test_01_Trial()throws Exception{
		assertEquals(0, 0);
	}

	@Test
	public void test_02_OpenImage() throws Exception{
		String path="src/test/resources/data/test_1/SimulatedData_echoes.tif";
		ImagePlus img=IJ.openImage(path);
		img.show();
		img.close();
	}


	// Test main computation features of the plugin
	@Test
	public void test_03_RegisterAndComputeMapsCircle()throws Exception{
		//switch variable to use or not use visual debug information
//		boolean useVisualDebugThatIsVeryInterestingButImplySlowerComputation=true;
		//Get FijiRelax archive directory path
		String mainDir="src/test/resources/data/test_1/";
		File f1=new File(mainDir,"SimulatedData_echoes.tif");
		String sourcePathSimulatedEchoes=f1.getAbsolutePath();
		File f2=new File(mainDir,"SimulatedData_expected_maps.tif");
		String sourcePathExpectedMaps=f2.getAbsolutePath();
		ImagePlus imgEchoesInit=IJ.openImage(sourcePathSimulatedEchoes);
		HyperMap hyperFantom=new HyperMap(imgEchoesInit); 

		//Compute maps
		//hyperFantom.computeMaps();
		HyperMap hyperFantom2=HyperMap.hyperMapFactory(hyperFantom);
		hyperFantom2.registerEchoes();
		hyperFantom2.computeMaps();
		ImagePlus computedMaps=hyperFantom2.getMapsImage();
		computedMaps=new Duplicator().run(computedMaps,1,2,1,1,1,1);
		ImagePlus expectedMaps=IJ.openImage(sourcePathExpectedMaps);
		VitimageUtils.printImageResume(computedMaps);
		VitimageUtils.printImageResume(expectedMaps);
		ImagePlus diff=VitimageUtils.substract(computedMaps,expectedMaps, true);
		computedMaps.show();
		expectedMaps.show();
		diff.show();
//		VitimageUtils.waitFor(1000000000);
		int val=VitimageUtils.isNullFloatImage(diff) ? 0 : 1;
		System.out.println("Val="+val);
		//assertEquals(val, 0);
	}




}
