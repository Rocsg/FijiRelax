package com.vitimage.aplimtools;

import java.io.File;
import java.text.DateFormat;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.Date;

import com.vitimage.common.VitimageUtils;
import com.vitimage.mrutils.MRDataType;

import ij.IJ;

public class Fid {
	public static final String DAY_NO_TIME_LAPSE="NODAY";
	public String path;
	public String fileName;
	public String fileParent;
	public String name;
	public LocalDateTime dateAcq;
	public LocalDateTime dateSubmitted;
	public int versionNumber;
	public boolean isSelected=false;
	public String specimen;
	public String fungus;
	public String day;
	public double daysAfterStart=0;
	public MRDataType typeAcq;
	public String nameAcq;
	public String dayCountFromStart;
	public int dayIndexInExp;

	public Fid(String path) {
		this.path=path;
		this.fileName=new File(path).getName();
		this.fileParent=new File(path).getParent();
		String[]tab=VitimageUtils.readStringFromFile(this.path+"/procpar").split("\n");
		boolean foundName=false,foundDate=false,foundSubmit=false;String s3=null,s4=null,s5=null,s6=null;
		for(String s2 : tab) {
			if(foundName) {
				this.name=s2.split("\\\"")[1];
				foundName=false;
			}
			if(foundSubmit) {
				s5=s2.split("\\\"")[1].replace("T","").substring(0,14);
				foundSubmit=false;
			}
			if(foundDate) {
				s3=s2.split("\\\"")[1].replace("T","").substring(0,14);
				foundDate=false;
			}
			if(s2.startsWith("name"))foundName=true;
			if(s2.startsWith("time_run"))foundDate=true;
			if(s2.startsWith("time_submitted"))foundSubmit=true;
		}
		DateTimeFormatter formatter = DateTimeFormatter.ofPattern("yyyyMMddHHmmss");
		this.dateAcq=LocalDateTime.parse(s3,formatter);
		this.dateSubmitted=LocalDateTime.parse(s5,formatter);
		this.versionNumber =Integer.parseInt(this.fileName.substring(this.fileName.length()-6).replace(".fid",""));
		this.name=this.name.replace("__","_");
		this.name=this.name.replace("3D_HR","3D");
		if(this.name.contains("_bis")){this.versionNumber++;this.name=this.name.replace("_bis", "");}
		if(this.name.equals("B051_CT_J133_B051_CT_J133_TR10000"))this.name="B051_CT_J133_TR10000";
		this.specimen=this.name.split("_")[0];
		this.fungus=this.name.split("_")[1];
		this.day=this.name.split("_")[2];
		if(this.day.equals("J1"))this.day="J0";
		
		if(this.name.split("_").length==3) {
			this.day="J0";
			this.typeAcq=parseMRType(this.name.split("_")[2]);
		}
		else if(this.name.split("_").length>4 && this.name.split("_")[4].contains("HR")){
			this.typeAcq=parseMRType(this.name.split("_")[3]+"_HR");
			this.nameAcq=this.name.split("_")[3]+"_HR";
		}
		else {			
			this.typeAcq=parseMRType(this.name.split("_")[3]);
			this.nameAcq=this.name.split("_")[3];
		}
		if(this.typeAcq==MRDataType.OTHER)IJ.showMessage("Oula : "+this.toFullString());
		this.correctKnownMistakeIfNeeded();
	}
	
	
	
	public void correctKnownMistakeIfNeeded() {
		DateTimeFormatter formatter = DateTimeFormatter.ofPattern("yyyyMMddHHmmss");
		if(this.dateAcq.equals(LocalDateTime.parse("20180410190641",formatter))){this.specimen="B013";}
		if(this.dateAcq.equals(LocalDateTime.parse("20180410191701",formatter))){this.specimen="B013";}
		if(this.dateAcq.equals(LocalDateTime.parse("20180410190127",formatter))){this.specimen="B013";}
		if(this.dateAcq.equals(LocalDateTime.parse("20180410193738",formatter))){this.specimen="B013";}
		if(this.dateAcq.equals(LocalDateTime.parse("20180410210342",formatter))){this.specimen="B013";}
		if(this.dateAcq.equals(LocalDateTime.parse("20190722072023",formatter))){this.day="J133";}
		if(this.dateAcq.equals(LocalDateTime.parse("20200219125049",formatter))){this.day="J100";}
	}
	

	public static MRDataType parseMRType(String s) {
		if(s.contains("TR10000_HR"))return MRDataType.T2SEQHR;
		if(s.contains("3D"))return MRDataType.GE3D;
		if(!s.contains("TR")){IJ.showMessage("Warning in Fid : unexpected type : "+s);return MRDataType.OTHER;}
		if(s.contains("TR600"))return MRDataType.TR600;
		if(s.contains("TR1200"))return MRDataType.TR1200;
		if(s.contains("TR2400"))return MRDataType.TR2400;
		if(s.contains("TR4800"))return MRDataType.TR4800;
		if(s.contains("TR10000"))return MRDataType.TR10000;
		return MRDataType.OTHER;
	}
		
	
	public String toString() {
		String s="";
		s+="Fid "+(isSelected ? "X" : " ")+"[ Date = "+(dateAcqString())+"(d0+"+VitimageUtils.dou(daysAfterStart)+") | Specimen = "+specimen+" | Fungus = "+fungus+" | DayExp = "+
				day+(day.length()==2 ? "  " : day.length()==3 ? " " : "")+" | TypeAcq = "+typeAcq+" | Version = "+versionNumber+" ]";
		return s;
	}
	
	public String toFullString() {
		String s="";
		s+="Fid "+(isSelected ? "X" : " ")+" DateSubmit= "+(dateSubString())+" | DateRun= "+(dateAcqString())+"(d0+"+VitimageUtils.dou(daysAfterStart)+") | Specimen = "+specimen+" | Fungus = "+fungus+" | DayExp = "+
				day+(day.length()==2 ? "  " : day.length()==3 ? " " : "")+" \n     | TypeAcq = "+typeAcq+" | Version = "+versionNumber+" |  name="+this.name+
				"\n     "+" | Path = "+this.path;
		return s;
	}
	
	
	public String dateSubString() {
		DateTimeFormatter formatter = DateTimeFormatter.ofPattern("yyyy-MM-dd HH:mm:ss");
	    return this.dateSubmitted.format(formatter);
	}
	
	public String dateAcqString() {
		DateTimeFormatter formatter = DateTimeFormatter.ofPattern("yyyy-MM-dd HH:mm:ss");
	    return this.dateAcq.format(formatter);
	}

	
	public boolean equals(Fid f) {
		if(!  specimen.equals(f.specimen ))return false;		
		if(!   fungus.equals(f.fungus ))return false;
		if(!   day.equals(f.day ))return false;
		if(!   typeAcq.equals(f.typeAcq ))return false;
		if(!   (versionNumber==f.versionNumber))return false;
		return true;
	}
	
}
