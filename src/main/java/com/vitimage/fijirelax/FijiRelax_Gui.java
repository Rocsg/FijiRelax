package com.vitimage.fijirelax;

import java.awt.Color;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.GridLayout;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.File;
import java.lang.management.ManagementFactory;
import java.net.URISyntaxException;
import java.util.ArrayList;
import java.util.Random;

import javax.swing.BorderFactory;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSeparator;
import javax.swing.JTextArea;

import com.vitimage.common.VitiDialogs;
import com.vitimage.common.VitimageUtils;
import com.vitimage.fijiyama.RegistrationManager;
import com.vitimage.registration.OptimizerType;
import com.vitimage.registration.Transform3DType;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.plugin.Memory;
import ij.plugin.frame.PlugInFrame;
import ij.plugin.frame.RoiManager;

public class FijiRelax_Gui extends PlugInFrame  implements ActionListener {

	//Buttons, Frames, Panels
	//Registration Frame buttons
	public JButton importButton=new JButton("Import Nifti / Dicom data");
	private JButton openButton=new JButton("Open existing hypermap");
	private JButton registerButton=new JButton("Register T1 / T2 sequences");
	private JButton computeButton=new JButton("Compute PD / T1 / T2 maps");
	private JButton sosButton = new JButton("Sos");
	private JButton explorerButton=new JButton("Explorer");
	private JButton undoButton=new JButton("Undo");
	private JButton abortButton=new JButton("Abort");
	private JButton stressButton=new JButton("I'm feeling stressed...");

	//Some more Gui objects and constants
	public JFrame registrationFrame;
	private JTextArea logArea=new JTextArea("", 11,10);
	private Color colorIdle;
	private int screenHeight=0;
	private int screenWidth=0;
	public int[] lastViewSizes=new int[] {700,700};//Only useful for serie running
	public final String displayedNameImage1="Image 1";
	public final String displayedNameImage2="Image 2";
	public final String displayedNameCombinedImage="Data_combined";
	public final String displayedNameHyperImage="Data_combined";
	private final int waitingTimeHyperImage=30;

	private Font mySurvivalFontForLittleDisplays=null;
	private boolean undoButtonHasBeenPressed=false;

	public String versionName="Handsome honeysuckle";
	public String timeVersionFlag="  Release time : 2020-12-16 - 21:39 PM";
	public String versionFlag=versionName+timeVersionFlag;
	public ImagePlus imgView;
	private Color colorStdActivatedButton;
	private Color colorGreenRunButton;
	public boolean interfaceIsRunning=false;

	private static final String saut="<div style=\"height:1px;display:block;\"> </div>";
	private static final String startPar="<p  width=\"650\" >";
	private static final String nextPar="</p><p  width=\"650\" >";

	//Identifiers for buttons
	private static final int IMPORT=91;
	private static final int OPEN=92;
	private static final int REGISTER=93;
	private static final int COMPUTE=94;
	private static final int EXPLORER=95;
	private static final int SOS=96;
	private static final int ABORT=97;
	private static final int STRESS=98;
	private static final int UNDO=99;

	private volatile boolean actionAborted=false;
	boolean developerMode=false;
	String spaces="                                                                                                                                   ";
	public int viewSlice;

	private boolean debugMode=true;
	private boolean isSurvivorVncTunnelLittleDisplay=false;
	private int nbCpu;
	private int jvmMemory;
	private long memoryFullSize;
	private HyperMap currentHypermap;
	private ArrayList<HyperMap>hypermapList=new ArrayList<HyperMap>();
	private ArrayList<ImagePlus>imgShowedList=new ArrayList<ImagePlus>();

	/**
	 * 
	 * 
	 */
	private static final long serialVersionUID = 1L;

	public static void main(String[]args) {		
		ImageJ ij=new ImageJ();
		FijiRelax_Gui fj=new FijiRelax_Gui();
		fj.run("");
	}
	
	public FijiRelax_Gui() {
		super("");
		// TODO Auto-generated constructor stub
	}

	public void run(String arg) {
		if(new File("/users/bionanonmri/").exists()) {
			this.isSurvivorVncTunnelLittleDisplay=true;
			IJ.showMessage("Detected Bionano server. \nSurvival display, but numerous cores");
		}
		startFijiRelaxInterface();
		welcomeAndInformAboutComputerCapabilities();
		
	}
	
	public void welcomeAndInformAboutComputerCapabilities() {		
		String[]str=checkComputerCapacity(true);
		addLog(str[0],0);
		addLog(str[1],0);		
	}

	
	public void addLog(String t,int level) {
		logArea.append((level==0 ? "\n > ": (level==1) ? "\n " : " ")+t);
		logArea.setCaretPosition(logArea.getDocument().getLength());
	}

	public void displaySosMessage(int context){
		IJ.showMessage(
		"Start menu\n"+
		"* First trial ? Choose register two images and test Fijiyama with demo images or with your own data.\n"+
		"* Ready to process a N-time M-modalities series ? Choose Register 3d series\n"+
		"* Go on a previous experiment (two images or series) ? Choose Open a previous study\n"+
		"More information ? Visit the webpage : www.imagej.net/Fijiyama\n"
		);
	}

	
	@SuppressWarnings("restriction")
	public String []checkComputerCapacity(boolean verbose) {
		this.nbCpu=Runtime.getRuntime().availableProcessors();
		this.jvmMemory=(int)((new Memory().maxMemory() /(1024*1024)));//Java virtual machine available memory (in Megabytes)
		this.memoryFullSize=0;
		String []str=new String[] {"",""};

		str[0]="Welcome to FijiRelax "+versionFlag+" ! \nFirst trial ? Click on \"Contextual help\" to get started. \n User requests, issues ? Contact : romain.fernandezATcirad.fr\n";
		str[1]="System check. Available memory in JVM="+this.jvmMemory+" MB. #Available processor cores="+this.nbCpu+".";
		try {
			this.memoryFullSize =  (( ((com.sun.management.OperatingSystemMXBean) ManagementFactory
		        .getOperatingSystemMXBean()).getTotalPhysicalMemorySize() )/(1024*1024));
		}	catch(Exception e) {return str;}		
		if((this.memoryFullSize>this.jvmMemory*2) && (this.memoryFullSize-this.jvmMemory>4000))  {
			str[1]+="\nIt seems that your computer have more memory : total memory="+(this.memoryFullSize)+" MB.\n"+
					"Registration and map computation are time and memory consuming. To give more memory and computation power to Fiji,"+
					" close the plugin then use the Fiji menu \"Edit / Options / Memory & threads\". "+
					"Let at least "+VitimageUtils.getSystemNeededMemory()+" unused to keep your "+VitimageUtils.getSystemName()+" stable and responsive.";
		}
		if(verbose)		return str;
		else return new String[] {"",""};
	}

	
	
	/* Registration Manager gui  and launching interface gui ************************************************************************************************/
	public void startFijiRelaxInterface() {
		IJ.log("Starting FijiRelax registration interface");

		//Panel with console-style log informations and requests
		JPanel consolePanel=new JPanel();
		consolePanel.setBorder(BorderFactory.createEmptyBorder(10,10,10,10));
		consolePanel.setLayout(new GridLayout(1,1,0,0));
		logArea.setSize(isSurvivorVncTunnelLittleDisplay ? 400 : 600,isSurvivorVncTunnelLittleDisplay ?  57 : 80);
		logArea.setBackground(new Color(10,10,10));
		logArea.setForeground(new Color(245,255,245));
		logArea.setFont(new Font(Font.DIALOG,Font.PLAIN,isSurvivorVncTunnelLittleDisplay ? 10 : 14));
		JScrollPane jscroll=new JScrollPane(logArea);
        logArea.setLineWrap(true);
        logArea.setWrapStyleWord(true);
        logArea.setEditable(false);	
        consolePanel.add(logArea);

       //Panel with step settings, used for registration of two images, and when programming registration pipelines for series
		JPanel stepSettingsPanel=new JPanel();
		mySurvivalFontForLittleDisplays=null;
		if(isSurvivorVncTunnelLittleDisplay ) {
			String name=openButton.getFont().getFamily();
			mySurvivalFontForLittleDisplays=new Font(name,Font.BOLD,11);
			importButton.setFont(mySurvivalFontForLittleDisplays);
			openButton.setFont(mySurvivalFontForLittleDisplays);
			registerButton.setFont(mySurvivalFontForLittleDisplays);
			computeButton.setFont(mySurvivalFontForLittleDisplays);
			sosButton.setFont(mySurvivalFontForLittleDisplays);
			undoButton.setFont(mySurvivalFontForLittleDisplays);
			stressButton.setFont(mySurvivalFontForLittleDisplays);
		}


		//settingsButton.setToolTipText("<html><p width=\"500\">" +"Advanced settings let you manage more parameters of the automatic registration algorithms"+"</p></html>");
//		setRunToolTip(MAIN);

		
		//Panel with buttons for the context of two image registration
		JPanel []buttonsPanel=new JPanel[4];
		for(int pan=0;pan<4;pan++) {
			buttonsPanel[pan]=new JPanel();
			if(isSurvivorVncTunnelLittleDisplay ) {
				buttonsPanel[pan].setBorder(BorderFactory.createEmptyBorder(0,0,5,0));
				buttonsPanel[pan].setLayout(new GridLayout(pan==1 ? 3 : 1,2,10,10));
			}
			else {
				buttonsPanel[pan].setBorder(BorderFactory.createEmptyBorder(25,25,25,25));
				buttonsPanel[pan].setLayout(new GridLayout(pan==1 ? 3 : 1,2,70,5));
			}
		}
		buttonsPanel[0].add(importButton);
		buttonsPanel[0].add(openButton);
		buttonsPanel[1].add(registerButton);
		buttonsPanel[1].add(computeButton);
		buttonsPanel[1].add(new JLabel(""));
		buttonsPanel[1].add(new JLabel(""));
		buttonsPanel[1].add(undoButton);
		buttonsPanel[1].add(abortButton);
		buttonsPanel[2].add(explorerButton);
		buttonsPanel[3].add(sosButton);
		buttonsPanel[3].add(stressButton);
			
		importButton.addActionListener(this);
		openButton.addActionListener(this);
		registerButton.addActionListener(this);
		computeButton.addActionListener(this);
		sosButton.addActionListener(this);
		undoButton.addActionListener(this);
		stressButton.addActionListener(this);
		abortButton.addActionListener(this);

		colorStdActivatedButton=openButton.getBackground();
		colorGreenRunButton=new Color(100,255,100);
			
		int width=isSurvivorVncTunnelLittleDisplay ? 350 : 500;
			
		this.colorIdle=abortButton.getBackground();
//		disable(ABORT);
//		enable(new int[] {RUN,UNDO,SAVE,FINISH,SETTINGS});



		//Main frame and main panel
		registrationFrame=new JFrame();
		JPanel registrationPanelGlobal=new JPanel();
		if(isSurvivorVncTunnelLittleDisplay ) {
			registrationFrame.setSize(600,680);
			registrationPanelGlobal.setBorder(BorderFactory.createEmptyBorder(5,5,5,5));
		}
		else {
			registrationFrame.setSize(750,850);
			registrationPanelGlobal.setBorder(BorderFactory.createEmptyBorder(5,5,5,5));			
		}
		registrationPanelGlobal.setLayout(new BoxLayout(registrationPanelGlobal, BoxLayout.Y_AXIS));
		String[]panelText=new String[] {"Open data               ","Process data","Graphical explorer","Need help ?"};
		for(int pan=0;pan<4;pan++) {
			registrationPanelGlobal.add(new JSeparator());
			JPanel tmpPanel=new JPanel();
			tmpPanel.setLayout(new GridLayout(2,1,10,10));
			tmpPanel.add(new JLabel(panelText[pan]));
			tmpPanel.add(new JLabel(""));
			registrationPanelGlobal.add(tmpPanel);			
			registrationPanelGlobal.add(buttonsPanel[pan]);
			registrationPanelGlobal.add(new JLabel(""));
			registrationPanelGlobal.add(new JLabel(""));
		}
		registrationFrame.setLayout(new BoxLayout(registrationFrame.getContentPane(), BoxLayout.Y_AXIS));
		//registrationFrame.setLayout(new BoxLayout(registrationFrame, BoxLayout.Y_AXIS));
		registrationFrame.add(consolePanel);
		registrationFrame.add(registrationPanelGlobal);
		registrationFrame.setTitle("FijiRelax ");
		registrationFrame.pack();
		if(isSurvivorVncTunnelLittleDisplay ) {
			registrationFrame.setSize(600,680);
		}
		else {
			registrationFrame.setSize(750,850);
		}
		
		registrationFrame.setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
		registrationFrame.addWindowListener(new WindowAdapter(){
             public void windowClosing(WindowEvent e){
                   IJ.showMessage("See you next time !");
                   registrationFrame.setVisible(false);
                   closeAllViews();
               }
		});
		registrationFrame.setVisible(true);
		registrationFrame.repaint();
		VitimageUtils.adjustFrameOnScreen(registrationFrame,2,0);		
		logArea.setVisible(true);
		logArea.repaint();
		
		
		//TODO : setup the tooltips
		importButton.setToolTipText("<html><p width=\"500\">" +	"Import dicom dir, nifti 4D or series of nifti 3D."+"</p></html>");
		openButton.setToolTipText("<html><p width=\"500\">" +	"Open a previously imported Hypermap from 4D nifti or 4D tif file."+"</p></html>");
		registerButton.setToolTipText("<html><p width=\"500\">" +	"Import dicom dir, nifti 4D or series of nifti 3D."+"</p></html>");
		computeButton.setToolTipText("<html><p width=\"500\">" +	"Import dicom dir, nifti 4D or series of nifti 3D."+"</p></html>");
		abortButton.setToolTipText("<html><p width=\"500\">" +	"Import dicom dir, nifti 4D or series of nifti 3D."+"</p></html>");
		undoButton.setToolTipText("<html><p width=\"500\">" +	"Import dicom dir, nifti 4D or series of nifti 3D."+"</p></html>");
		explorerButton.setToolTipText("<html><p width=\"500\">" +	"Import dicom dir, nifti 4D or series of nifti 3D."+"</p></html>");
		sosButton.setToolTipText("<html><p width=\"500\">" +	"Import dicom dir, nifti 4D or series of nifti 3D."+"</p></html>");
		stressButton.setToolTipText("<html><p width=\"500\">" +	"Import dicom dir, nifti 4D or series of nifti 3D."+"</p></html>");
		disable(new int[] {REGISTER,COMPUTE,ABORT,UNDO,EXPLORER});
		enable(OPEN);
		enable(IMPORT);
		enable(SOS);
		enable(STRESS);
	}

	@Override
	public void actionPerformed(ActionEvent e) {
		// TODO Auto-generated method stub
		if(e.getSource()==this.importButton)runActionImport();
		if(e.getSource()==this.openButton)runActionOpen();
		if(e.getSource()==this.registerButton)runActionRegistration();
		if(e.getSource()==this.computeButton)runActionCompute();
		if(e.getSource()==this.abortButton)runActionAbort();
		if(e.getSource()==this.undoButton)runActionUndo();
		if(e.getSource()==this.explorerButton)runActionExplorer();
		if(e.getSource()==this.sosButton)runActionSos();
		if(e.getSource()==this.stressButton)runActionStress();
	}


	
	

	//TODO : Open from nifti data
	public void runActionOpen() {
		ImagePlus imgNew=VitiDialogs.chooseOneImageUI("Hypermap selection", "Choose nifti or tif hypermap");
		if(imgNew!=null) {
			this.currentHypermap=new HyperMap(imgNew);
			this.hypermapList.add(this.currentHypermap);
		}
		updateView();
		enable(new int[] {REGISTER,COMPUTE,ABORT,UNDO,EXPLORER,SOS,STRESS} );
	}
	

	//TODO : Open from dicom data and nifti 3D
	public void runActionImport() {		
		VitiDialogs.notYet("runActionImport");		
	}

	
	
	
	
	
	//TODO : Test registration
		public void runActionRegistration() {
		HyperMap tmp=HyperMap.hyperMapFactory(this.currentHypermap);
		tmp.registerEchoes();
		hypermapList.add(tmp);
		updateView();
	}

	//TODO : Adapt thresholds to data
	public void runActionCompute() {
		HyperMap tmp=HyperMap.hyperMapFactory(this.currentHypermap);
		tmp.computeMaps();
		hypermapList.add(tmp);
		updateView();
	}
	
	//TODO : threading
	public void runActionAbort() {
		VitiDialogs.notYet("runActionAbort");
	}
	
	public void runActionUndo() {
		if(hypermapList.size()<2)return;
		this.currentHypermap=hypermapList.get(hypermapList.size()-2);
		this.hypermapList.remove(hypermapList.size()-1);
		this.imgShowedList.remove(imgShowedList.size()-1);
		updateView();
	}

	//TODO : set up back the explorer
	public void runActionExplorer() {		
		VitiDialogs.notYet("runActionExplorer");
	}
	


	
	
	
	
	public int getRelativeOptimalPositionFor2DView() {
		java.awt.Dimension currentScreen = Toolkit.getDefaultToolkit().getScreenSize();
        int screenX=(int)Math.round(currentScreen.width);
        if(screenX>1920)screenX/=2;        
		if(registrationFrame.getLocationOnScreen().x+registrationFrame.getSize().getWidth()/2 > screenX/2) return 0;
		else return 2;
	}

	
	
	public void updateView() {
		if(this.hypermapList.size()==0)return;
		ImagePlus temp=null;
		boolean firstView=false;
		if(this.imgView==null || (!this.imgView.isVisible()))firstView=true;
		
		if(!firstView) {
			this.viewSlice=this.imgView.getZ();
			this.imgView.hide();
		}
		else this.viewSlice=(this.hypermapList.get(this.hypermapList.size()-1).getAsImagePlus().getNSlices())/2;
		
		
		this.imgView=this.hypermapList.get(this.hypermapList.size()-1).getAsImagePlus();
		if(firstView)this.viewSlice = this.imgView.getNSlices()/2;

		
		this.imgView.setZ(this.viewSlice);
		imgView.show();
		this.imgShowedList.add(imgView);
		VitimageUtils.adjustFrameOnScreenRelative(imgView.getWindow(),registrationFrame,getRelativeOptimalPositionFor2DView(),0,10);

		java.awt.Rectangle w = imgView.getWindow().getBounds();
		int max=0;

		
		//If little image, enlarge it until its size is between half screen and full screen
		while(imgView.getWindow().getWidth()<(screenWidth/2) && imgView.getWindow().getHeight()<(screenHeight/2) && (max++)<4) {
			int sx=imgView.getCanvas().screenX((int) (w.x+w.width));
			int sy=imgView.getCanvas().screenY((int) (w.y+w.height));
			imgView.getCanvas().zoomIn(sx, sy);
			VitimageUtils.waitFor(50);
			this.lastViewSizes=new int[] {imgView.getWindow().getWidth(),imgView.getWindow().getHeight()};
		}

		
		//If big image, reduce it until its size is between half screen and full screen
		while(imgView.getWindow().getWidth()>(screenWidth) || imgView.getWindow().getHeight()>(screenHeight)) {
			int sx=imgView.getCanvas().screenX((int) (w.x+w.width));
			int sy=imgView.getCanvas().screenY((int) (w.y+w.height));
			imgView.getCanvas().zoomOut(sx, sy);
			this.lastViewSizes=new int[] {imgView.getWindow().getWidth(),imgView.getWindow().getHeight()};
		}
		VitimageUtils.adjustFrameOnScreenRelative(imgView.getWindow(),registrationFrame,0,0,10);

		imgView.setSlice(viewSlice);
		imgView.updateAndRepaintWindow();
	}
	
	public static String getImgViewText(int st){
		return ( (st==0) ? "Superimposition before registration" :  ("Registration results after "+(st)+" step"+((st>1)? "s" : "")) );
	}
	
	public void closeAllViews() {
		for(ImagePlus img : imgShowedList) {if(img!=null)img.close();};
	}	

	
	
	
	
	
	
	
	public void enable(int but) {
		setState(new int[] {but},true);
	}
	public void disable(int but) {
		setState(new int[] {but},false);
	}

	public void enable(int[]tabBut) {
		setState(tabBut,true);
	}
	public void disable(int[]tabBut) {
		setState(tabBut,false);
	}
			
	public void setState(int[]tabBut,boolean state) {
		for(int but:tabBut) {
			switch(but) {
			case OPEN:this.openButton.setEnabled(state);break;
			case IMPORT:this.importButton.setEnabled(state);break;
			case REGISTER:this.registerButton.setEnabled(state);break;
			case COMPUTE:this.computeButton.setEnabled(state);break;
			case ABORT:this.abortButton.setEnabled(state);break;
			case UNDO:this.undoButton.setEnabled(state);break;
			case EXPLORER:this.explorerButton.setEnabled(state);break;
			case SOS:this.sosButton.setEnabled(state);break;
			case STRESS:this.stressButton.setEnabled(state);break;
			}	
		}	
	}


	
	public void runActionStress(){
		getRelaxingPopup("Feeling stressed with FijiRelax ?",true);
	}

	
	
	public void runActionSos(){
		String textToDisplay="";
		String basicFijiText="<b> </b>FijiRelax is a tool for.....,"+
				" (texte 2 ........) ";
		String mainWindowText=basicFijiText+
				"Texte3 ........ :"+startPar+
				" <b>1)   Texte4........ "+nextPar+
				" <b>2-a) Texte5 ........ ."+nextPar+
				"automatic algorithms settings can be modified using the settings dialog ( <b>\"Advanced settings\"</b> button )."+saut+saut; 

		
	
		disable(SOS);
		textToDisplay+="<b>Citing this work :</b> R. Fernandez and C. Moisy, <i>FijiRelax : fast and noise-corrected estimation of NMR relaxation maps in 3D+t</i> (under review)"+saut+
				"<b>Credits :</b> this work was supported by the \"Plan deperissement du vignoble\"   </p>";
		IJ.showMessage("FijiRelax help", textToDisplay);
		enable(SOS);
	}
	
	
	
	
	public static void getRelaxingPopup(String introductionSentence,boolean external) {
		File testf=null;
		try {		testf = new File( FijiRelax_Gui.class.getResource( File.separator+"Zen_Quotes.txt" ).toURI() );	} catch (URISyntaxException e) {			e.printStackTrace();		}
		String []str=VitimageUtils.readStringFromFile(testf.getAbsolutePath()).split("\n");
		Random rand=new Random();
		int line=rand.nextInt(str.length-1);
		String message=""+introductionSentence+"\n";
		message+="\nRead these words carefully :\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n#"+(""+line)+" "+str[line]+"\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"+
				"\n\nKeep this sentence in mind while taking three deep breaths.\n When you feel confident to go back to work, click OK"; 
		IJ.showMessage(message);
	}

	
	
}
