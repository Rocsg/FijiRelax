package com.vitimage.fijirelax;

import ij.IJ;
import ij.plugin.frame.PlugInFrame;

public class FijiRelax_Gui   extends PlugInFrame{

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	public static void main(String[]args) {		
		FijiRelax_Gui fj=new FijiRelax_Gui();
		fj.run("");
	}
	
	public FijiRelax_Gui() {
		super("");
		// TODO Auto-generated constructor stub
	}

	public void run(String arg) {
		IJ.showMessage("Plugin is started !");

	}
	
	
}
