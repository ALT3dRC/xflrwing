#include <istream>
#include <iterator>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <stdio.h>
//#include <cmath>
#include <sstream>
#include <math.h>      
#include <ctime>
#include <cstdlib>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include "spline.h"
#define XXX 512
#define PI 3.14159265
#define SSTR( x ) static_cast< std::ostringstream & >(( std::ostringstream() << std::dec << x ) ).str()

using namespace std;

double sx[XXX][501],sxe[XXX][501],sy[XXX],sz[XXX][501],sze[XXX][501],thl[XXX],tht[XXX],tll[XXX],tlt[XXX];
double lx0[XXX],ly0[XXX],tx0[XXX],ty0[XXX],lx1[XXX],ly1[XXX],tx1[XXX],ty1[XXX],lx2[XXX],ly2[XXX],tx2[XXX],ty2[XXX],lx3[XXX],ly3[XXX],tx3[XXX],ty3[XXX];
double cx[XXX][30],cy[XXX][30],cz[XXX][30],mx[251],mz[251];
double ww,rr,hh,ltd,lts,led,ttd,tts,ted,lew,tew,teh,leh,zbase,hingeprc,xspace=25.,yspace=25.,dh,sh,pmloz,pmhiz;
double rcjoiner1h,rcjoiner2h,rcjoiner3h,rcjoiner1w,rcjoiner2w,rcjoiner3w,rcjoiner1p,rcjoiner2p,rcjoiner3p;
double cmjoiner1h,cmjoiner2h,cmjoiner3h,cmjoiner1w,cmjoiner2w,cmjoiner3w,cmjoiner1p,cmjoiner2p,cmjoiner3p;
double mtjoiner1h,mtjoiner2h,mtjoiner3h,mtjoiner1w,mtjoiner2w,mtjoiner3w,mtjoiner1p,mtjoiner2p,mtjoiner3p;
int maxidx,cmbrk,mtbrk,mybrk,bs;

string line,txt_line,proj_name,file_name,file_name1,file_name2,airfoilt,airfoil[25],cfg_par,seg_text,lba_text,alt_text,tbf_text,htf_text,hla_text,alh_text,ath_text;
int imax,nr_seg,nrl,trigger,i,ii,j,jj,k,l,ll,ls,iy,iu,jy,ju,fo,wo,co,mo,po,to,last_k,nr_slice,cmin,cmax,cpos,f,wbrk,piy,cdih,mdih,tdih,ymax,zbreak,wtip,ztip,ypos,spos,semn,kla,ktf;
int auw,khlf, khlb,khtf,khtb,kh,khfa,mr,Q,wing_cfg_idx,twistsel,cpmoutput,dxfoutput,idxle[25],idxte[25],idxu[25],idxl[25],nodigits,flatmodel,ctrl1typ,ctrl2typ,ctrl3typ,ctrlprc;
int elliptic,sparprc,wingspan,rootspan,spos1,spos2,spos3,fpos,apos,dtype,iprc,spoiler,flaps,aileron,fspartype,bspartype,rib,smin,smax,amin,amax;
double airfoilx,sf,wp[51][501],ax[50][500],ay[50][500],aposy[50],achord[50],aoffset[50],aoffsetd[50],adihedral[50],adjust[50],xscale,yscale,xpower,cmpos,mtpos;
double tubemin,tubemax,posmax,atwist[50],aunk1[50],aunk2[50],dboxprc,twistl,twistr,posl,posr,maximh,camber,ctrldim,hingedim;
double xpow,chordl,chordr,offsetl,offsetr,offsetdl,offsetdr,srv1d1,srv1d2,srv1p1,srv1p2,srv1z,srv2d1,srv2d2,srv2p1,srv2p2,srv2z,srv3d1,srv3d2,srv3p1,srv3p2,srv3z;
double aunk3[50],aunk4[50],adat[50],px[2][500],pz[2][500],dx[501],dy,dz[501],ps[501],sv[501],x,y,z,r,nx,ny,x1,x2,y2,z1,z2,fi,posy,chord,offset,offsetd,la,lb,lc,ld,le,lp;
double px1,px2,px3,px4,py1,py2,py3,py4,dihedral,twist,unk1,unk2,unk3,unk4,deltax,deltaz,surface,csurface,volume,weight,xb,xf,xxb,xxf,minim,W,cfg_val,ytx;
double zoffb,zofft,zoffset1,zoffset2,zoffsetl,zoffsetr,zoffset,zoffsetf,zoffsetm,zoffsetc,xchscale,ywsscale,zhtscale,ymaxnew,flapsprc,flapsmin,flapsmax,flapsdim,ailerprc,yc,yt,yr;
double ailermin,ailermax,ailerdim,resinfiber,moldpos,zedim,moldmargin,moldheight,xc1,xc2,yc1,yc2,zc1,zc2,xt1,xt2,yt1,yt2,zt1,hinge,posydx,estfin,estelev,estwing,estfuse;
double estauw,maxim,ctrl1d1,ctrl1d2,ctrl1dim,ctrl1prc,ctrl2d1,ctrl2d2,ctrl2dim,ctrl2prc,ctrl3d1,ctrl3d2,ctrl3dim,ctrl3prc,cfweight,xpsweight,deadweight,rcy,zt2,zch,zth;
double xch,xth,xh,yh,zh,uhx,uhy,uhz,zlo,zhi,xtrans,ytrans,tipextra,xcnc[500],zcnc[500],zc1c,zc2c,zt1c,zt2c,ribspace,balsaxmax,balsaymax,sbalsasize,sparprc1,sparprc2;
double tesw,tesh,lesw,lesh,ort,orl,wingpos,cx0,cy0,cx1,cx2,cx3,cx4,cy1,cy2,cy3,cy4,dcx,dcy,dr,wingmax,sparcaph,scw,ribdist,sparwidth1,sparwidth2,xe,ye,xo,yo,xes,ddx,dbox;
double cut,ribwidth,tubeh,tubez,tubew,tuberef,dboxe[401],dboxi[401],lboxi,lboxe,oxn,oyn,wxn,wyn,exn,eyn,oxs,oys,wxs,wys,exs,eys,oxl,oyl,wxl,wyl,exl,eye,oxr,oyr,wxr,wyr,exr,eyr;
double xne[251],xni[251],yne[251],yni[251],xle[251],xli[25],yle[251],yli[251],xre[251],xri[251],yre[251],yri[251];
double xse[251],xsi[251],yse[251],ysi[251],cdx,compfoil,pos,cntsp,tipsp,hi,lo,hl,compfoill,compfoils,compfoilm,compfoilt,cf,xlc,zlc,xtc,ztc,xlb,zlb,xtb,ztb,xlt,zlt,xtt,ztt;	
tk::spline sch,sxo,szo,sda,swa,ssp,spr[501]; 
ofstream df,stl;

int between(double x, double x1, double x2)
{
	return (((x>=x1)&&(x<=x2))||((x>=x2)&&(x<=x1))?1:0);; 
};

int asy_sort()
{
	double temp;
	int min;
	for (int i=0; i<maxidx; i++) 
		{
			min=i;
			for (j=i+1; j<=maxidx; j++) if (sy[j]<sy[min]) min=j;
			temp=sy[i]; 	sy[i]=sy[min];		sy[min]=temp;
        };
	return maxidx;
};

int asy_compact()
{ 
	int i,j=0; double y[999];
	asy_sort();
	y[0]=sy[0]; j=1;
	for (i=1;i<=maxidx;i++) if (!between(sy[i],sy[i-1]-0.001,sy[i-1]+0.001)) y[j++]=sy[i];
	maxidx=j-1;
	for (i=0;i<=maxidx;i++) sy[i]=y[i];
	return maxidx;
};

int asy_find (double yval)
{
	int pos=-1;
	for (int i=0;i<=maxidx;i++) pos=(between(yval,sy[i]-0.001,sy[i]+0.001)?i:pos);
	return (pos);
};

int asy_add (double yval)
{
	if (yval>=0.0) {sy[++maxidx]=yval;};
	return (maxidx);
};

double limited(double v, double w)
{
	return (v<=0.0?max(v,-w):min(v,w));
};
	
int stlv1(int cw, ofstream &stl_file,double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3)
 {
	stl_file << "  facet normal 0.0 0.0 0.0" << endl; 
	stl_file << "    outer loop" << endl;
	stl_file << "      vertex " << fixed << setprecision(3) << x1           << " " << y1           << " " << z1			  << endl;
	stl_file << "      vertex " << fixed << setprecision(3) << (cw>0?x2:x3) << " " << (cw>0?y2:y3) << " " << (cw>0?z2:z3) << endl;
	stl_file << "      vertex " << fixed << setprecision(3) << (cw>0?x3:x2) << " " << (cw>0?y3:y2) << " " << (cw>0?z3:z2) << endl;
	stl_file << "    endloop" << endl; stl_file << "  endfacet" << endl; 
	return 0;
};

int stlv2(int cw, ofstream &stl_file,double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3, double x4, double y4, double z4)
{ 
	stlv1(cw,stl_file,x1,y1,z1,x2,y2,z2,x4,y4,z4);	
	stlv1(cw,stl_file,x3,y3,z3,x1,y1,z1,x4,y4,z4);
	return 0;
};

int stlvc(int cw, ofstream &stl_file,double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3, double x4, double y4, double z4)
{ 
	stlv1(cw,stl_file,x1,y1,z1,x2,y2,z2,x3,y3,z3);	
	stlv1(cw,stl_file,x1,y1,z1,x3,y3,z3,x4,y4,z4);
	return 0;
};

int stl_surface(int cw, int my, int dw, ofstream &stl_file, int biy, int eiy)
{
	double xx1,xx2,xx3,xx4,yy1,yy2,yy3,yy4,zz1,zz2,zz3,zz4;
	for (int i=biy; i<eiy; i++)
		{
			yy1=yy2=(my>0?sy[i]:-sy[i]);
			yy3=yy4=(my>0?sy[i+1]:-sy[i+1]);
			for (int j=((dw==-1)?250:0); j<((dw==1)?250:500); j++)
				{
					xx1=sx[i][j];
					xx2=sx[i][j+1];
					xx3=sx[i+1][j];
					xx4=sx[i+1][j+1];	
					zz1=sz[i][j];
					zz2=sz[i][j+1];
					zz3=sz[i+1][j];
					zz4=sz[i+1][j+1];
					stlv2(cw,stl_file,xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3,xx4,yy4,zz4);
				};
		};
	return 0;	
};

int is_hole(int i)
{
	double py,pyh;
	py=sy[i];
	pyh=round(py/dh)*dh;
	return between(py,pyh-sh/2.+0.001,pyh+sh/2.-0.001);
};	

int stl_contour(int cw, int my, int db, int dw, ofstream &stl_file, int biy, int eiy, int exb, int exe)
// cw clockwise, my y_mirrored, db hi/lo base, dw hi/lo wing
{
	double xle,yle,zle,xlw,zlw,xlt,zlt,xlh,zlh,xte,zte,xtw,ztw,xtt,ztt,xth,zth;
	double xx1,xx2,xx3,xx4,yy1,yy2,yy3,yy4,zz1,zz2,zz3,zz4;

    for (int i=biy; i<=eiy; i++)
		{
			xte=sx[i][0];
			xle=sx[i][250];
			zte=sz[i][0];
			zle=sz[i][250];
			ztt=zte+dw*tts/2.;
			zlt=zle+dw*lts/2.;
			ztw=zte-dw*teh;
			zlw=zle-dw*leh;
			zlh=zle+db*sh;
			zth=zte+db*sh;
			// point 0 is LE / aft leveled leading wall
			cx[i][0]=xle;
			cy[i][0]=sy[i];
			cz[i][0]=zle;
			// point 1 is aft displaced leading wall
			cx[i][1]=xle;
			cy[i][1]=sy[i];
			cz[i][1]=zlw;
			// point 2 is fore displaced leading wall
			cx[i][2]=lx0[i];
 			cy[i][2]=ly0[i];
			cz[i][2]=zlw;
			// point 3 is fore leveled leading wall			
			cx[i][3]=lx0[i]; 
			cy[i][3]=ly0[i];
			cz[i][3]=zle;
			// point 4 is aft leveled leading trench
			cx[i][4]=lx1[i]; 
			cy[i][4]=ly1[i]; 
			cz[i][4]=zle;
			// point 5 is aft displaced leading trench			
			cx[i][5]=lx1[i];
 			cy[i][5]=ly1[i]; 
			cz[i][5]=zlt;
			// point 6 is fore displaced leading trench			
			cx[i][6]=lx2[i];
 			cy[i][6]=ly2[i]; 
			cz[i][6]=zlt;
			// point 7 is fore leveled leading trench				
			cx[i][7]=lx2[i];
 			cy[i][7]=ly2[i];
			cz[i][7]=zle;
			// point 8 is aft leveled centering
			cx[i][8]=(lx2[i]+lx3[i]+sh)/2.; 
			cy[i][8]=(ly2[i]+ly3[i])/2.; 
			cz[i][8]=zle;
			// point 9 is aft displaced centering			
			cx[i][9]=(lx2[i]+lx3[i]+sh)/2.-1.; 
			cy[i][9]=(ly2[i]+ly3[i])/2.; 
			cz[i][9]=(is_hole(i)?zlh:zle);
			// point 10 is fore displaced centering			
			cx[i][10]=(lx2[i]+lx3[i]-sh)/2.+1.; 
			cy[i][10]=(ly2[i]+ly3[i])/2.; 
			cz[i][10]=(is_hole(i)?zlh:zle);
			// point 11 is fore leveled centering			
			cx[i][11]=(lx2[i]+lx3[i]-sh)/2.; 
			cy[i][11]=(ly2[i]+ly3[i])/2.; 
			cz[i][11]=zle;
			// point 12 is leveled leading point				
			cx[i][12]=lx3[i]; 
			cy[i][12]=ly3[i];
			cz[i][12]=zle;
			// point 13 is base / displaced leading point
			cx[i][13]=lx3[i];
 			cy[i][13]=ly3[i]; 
			cz[i][13]=(db==1?pmhiz:-pmloz);
			// point 14 is base / displaced leading edge			
			cx[i][14]=xle;
			cy[i][14]=sy[i];
			cz[i][14]=(db==1?pmhiz:-pmloz);
			// point 15 is base / displaced trailing edge			
			cx[i][15]=xte;
			cy[i][15]=sy[i]; 
			cz[i][15]=(db==1?pmhiz:-pmloz);
			// point 16 is base / displaced trailing point				
			cx[i][16]=tx3[i]; 
			cy[i][16]=ty3[i];
			cz[i][16]=(db==1?pmhiz:-pmloz);
			// point 17 is leveled trailing point			
			cx[i][17]=tx3[i]; 
			cy[i][17]=ty3[i]; 
			cz[i][17]=zte ;
			// point 18 is aft leveled trailing centering			
			cx[i][18]=(tx2[i]+tx3[i]+sh)/2.; 
			cy[i][18]=(ty2[i]+ty3[i])/2.;
			cz[i][18]=zte;
			// point 19 is aft displaced trailing centering			
			cx[i][19]=(tx2[i]+tx3[i]+sh)/2.-1.; 
			cy[i][19]=(ty2[i]+ty3[i])/2.;
			cz[i][19]=(is_hole(i)?zth:zte);
			// point 20 is fore displaced trailing centering			
			cx[i][20]=(tx2[i]+tx3[i]-sh)/2.+1.;	
			cy[i][20]=(ty2[i]+ty3[i])/2.;
			cz[i][20]=(is_hole(i)?zth:zte);
			// point 21 is fore leveled trailing centering  			
			cx[i][21]=(tx2[i]+tx3[i]-sh)/2.; 
			cy[i][21]=(ty2[i]+ty3[i])/2.; 
			cz[i][21]=zte;
			// point 22 is aft leveled trailing trench			
			cx[i][22]=tx2[i]; 
			cy[i][22]=ty2[i]; 
			cz[i][22]=zte;
			// point 23 is aft displaced trailing trench			
			cx[i][23]=tx2[i]; 
			cy[i][23]=ty2[i]; 
			cz[i][23]=ztt;
			// point 24 is fore displaced trailing trench			
			cx[i][24]=tx1[i];
			cy[i][24]=ty1[i];
			cz[i][24]=ztt;
			// point 25 is fore leveled trailing trench			
			cx[i][25]=tx1[i]; 
			cy[i][25]=ty1[i]; 
			cz[i][25]=zte;
			// point 26 is aft leveled trailing wall			
			cx[i][26]=tx0[i]; 
			cy[i][26]=ty0[i]; 
			cz[i][26]=zte;
			// point 27 is aft displaced trailing wall			
			cx[i][27]=tx0[i]; 
			cy[i][27]=ty0[i];
			cz[i][27]=ztw;
			// point 28 is fore displaced trailing wall			
			cx[i][28]=xte; 
			cy[i][28]=sy[i];
			cz[i][28]=ztw;
			// point 29 is TE / fore leveled trailing wall 				
			cx[i][29]=xte; 
			cy[i][29]=sy[i]; 
			cz[i][29]=zte;

		};	
	for (int i=biy; i<eiy; i++)
		{
			for (int j=0; j<=28; j++)
				{
					xx1=cx[i][j];
					yy1=(my>0?cy[i][j]:-cy[i][j]);
					zz1=cz[i][j];	
					xx2=cx[i][j+1];
					yy2=(my>0?cy[i][j+1]:-cy[i][j+1]);
					zz2=cz[i][j+1];				
					xx3=cx[i+1][j];
					yy3=(my>0?cy[i+1][j]:-cy[i+1][j]);
					zz3=cz[i+1][j];				
					xx4=cx[i+1][j+1];
					yy4=(my>0?cy[i+1][j+1]:-cy[i+1][j+1]);
					zz4=cz[i+1][j+1];
					stlv2(cw,stl_file,xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3,xx4,yy4,zz4);
				};
		};
	// center 
	if (exb!=0)
	{
		xx1=xx3=cx[biy][0];	
		xx2=xx4=cx[biy][12];
		yy1=yy3=yy2=yy4=my*(cy[biy][0]-exb*yspace);
		zz1=zz2=sz[biy][250];
		zz3=zz4=(db==1?pmhiz:-pmloz);  	
 	    stlv2(-cw,stl_file,xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3,xx4,yy4,zz4);
		
		xx1=xx3=cx[biy][17];
		xx2=xx4=cx[biy][29];
		yy1=yy3=yy2=yy4=my*(cy[biy][29]-exb*yspace); 
		zz1=zz2=sz[biy][0];
		zz3=zz4=(db==1?pmhiz:-pmloz);  	
		stlv2(-cw,stl_file,xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3,xx4,yy4,zz4);
		
		xx1=xx3=cx[biy][0];
		xx2=xx4=cx[biy][2];
		yy1=yy3=my*cy[biy][0];
		yy2=yy4=my*cy[biy][2];
		zz1=zz2=cz[biy][1];	zz3=zz4=cz[biy][0];  	
		stlv2(cw,stl_file,xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3,xx4,yy4,zz4);
		
		xx1=xx3=cx[biy][4];
		xx2=xx4=cx[biy][6];
		yy1=yy3=my*cy[biy][4];
		yy2=yy4=my*cy[biy][6];	
		zz1=zz2=cz[biy][5];
		zz3=zz4=cz[biy][4];  	
		stlv2(cw,stl_file,xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3,xx4,yy4,zz4);
		
		xx1=xx3=cx[biy][8];
		xx2=xx4=cx[biy][10];	
		yy1=yy3=my*cy[biy][8];
		yy2=yy4=my*cy[biy][10];
		zz1=zz2=cz[biy][9];
		zz3=zz4=cz[biy][8];  	
		stlv2(cw,stl_file,xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3,xx4,yy4,zz4);
		
		xx1=xx3=cx[biy][29];
		xx2=xx4=cx[biy][27];
		yy1=yy3=my*cy[biy][29];
		yy2=yy4=my*cy[biy][27];
		zz1=zz2=cz[biy][28];
		zz3=zz4=cz[biy][29];  	
		stlv2(cw,stl_file,xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3,xx4,yy4,zz4);

		xx1=xx3=cx[biy][25];
		xx2=xx4=cx[biy][23];
		yy1=yy3=my*cy[biy][25];
		yy2=yy4=my*cy[biy][23];
		zz1=zz2=cz[biy][24];
		zz3=zz4=cz[biy][25];  	
		stlv2(cw,stl_file,xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3,xx4,yy4,zz4);
		
		xx1=xx3=cx[biy][21];
		xx2=xx4=cx[biy][19];
		yy1=yy3=my*cy[biy][21]; 
		yy2=yy4=my*cy[biy][19];
		zz1=zz2=cz[biy][20];
		zz3=zz4=cz[biy][21];  	
		stlv2(cw,stl_file,xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3,xx4,yy4,zz4);
	};		
	if (exe!=0)
	{
		xx1=xx3=cx[eiy][0];	
		xx2=xx4=cx[eiy][12];
		yy1=yy3=yy2=yy4=my*(cy[eiy][0]+exe*yspace);
		zz1=zz2=sz[eiy][250];
		zz3=zz4=(db==1?pmhiz:-pmloz);  	
		stlv2(cw,stl_file,xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3,xx4,yy4,zz4);
	
		xx1=xx3=cx[eiy][17];
		xx2=xx4=cx[eiy][29];	
		yy1=yy3=yy2=yy4=my*(cy[eiy][29]+exe*yspace);
		zz1=zz2=sz[eiy][0];	
		zz3=zz4=(db==1?pmhiz:-pmloz);  	
		stlv2(cw,stl_file,xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3,xx4,yy4,zz4);
	
		xx1=xx3=cx[eiy][0];	
		xx2=xx4=cx[eiy][2];	
		yy1=yy3=my*cy[eiy][0]; 
		yy2=yy4=my*cy[eiy][2];	
		zz1=zz2=cz[eiy][1];	
		zz3=zz4=cz[eiy][0];  	
		stlv2(cw,stl_file,xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3,xx4,yy4,zz4);
	
		xx1=xx3=cx[eiy][4];
		xx2=xx4=cx[eiy][6];
		yy1=yy3=my*cy[eiy][4]; 
		yy2=yy4=my*cy[eiy][6];
		zz1=zz2=cz[eiy][5];	
		zz3=zz4=cz[eiy][4];  	
		stlv2(cw,stl_file,xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3,xx4,yy4,zz4);
	
		xx1=xx3=cx[eiy][8];
		xx2=xx4=cx[eiy][10];
		yy1=yy3=my*cy[eiy][8]; 
		yy2=yy4=my*cy[eiy][10];	
		zz1=zz2=cz[biy][9];	
		zz3=zz4=cz[eiy][8];  	
		stlv2(cw,stl_file,xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3,xx4,yy4,zz4);
	
		xx1=xx3=cx[eiy][29];
		xx2=xx4=cx[eiy][27];
		yy1=yy3=my*cy[eiy][29]; 
		yy2=yy4=my*cy[eiy][27];
		zz1=zz2=cz[eiy][28];
		zz3=zz4=cz[eiy][29];  	
		stlv2(cw,stl_file,xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3,xx4,yy4,zz4);
	
		xx1=xx3=cx[eiy][25];
		xx2=xx4=cx[eiy][23];
		yy1=yy3=my*cy[eiy][25]; 
		yy2=yy4=my*cy[eiy][23];
		zz1=zz2=cz[eiy][24];
		zz3=zz4=cz[eiy][25];  	
		stlv2(cw,stl_file,xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3,xx4,yy4,zz4);
	
		xx1=xx3=cx[eiy][21];
		xx2=xx4=cx[eiy][19];
		yy1=yy3=my*cy[eiy][21]; 
		yy2=yy4=my*cy[eiy][19];
		zz1=zz2=cz[eiy][20];
		zz3=zz4=cz[eiy][21];  	
		stlv2(cw,stl_file,xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3,xx4,yy4,zz4);
	};
	if (exb!=0)
		for (int j=0; j<=28; j++)
			{
				xx1=cx[biy][j];	
				xx2=cx[biy][j+1];
				xx3=cx[biy][j];
				xx4=cx[biy][j+1];	
				yy1=my*(sy[biy]-yspace);
				yy2=my*(sy[biy]-yspace);
				yy3=my*cy[biy][j];
				yy4=my*cy[biy][j+1];
				zz1=zz3=cz[biy][j];	
				if (j<12) zz1=zz3=sz[biy][250];
				if (j>17) zz1=zz3=sz[biy][0];
				zz2=zz4=cz[biy][j+1];
				if (j<12) zz2=zz4=sz[biy][250];	
				if (j>17) zz2=zz4=sz[biy][0];
				stlv2(cw,stl_file,xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3,xx4,yy4,zz4);
			};
	if (exe!=0)
		for (int j=0; j<=28; j++)
			{
				xx1=cx[eiy][j];
				xx2=cx[eiy][j+1];
				xx3=cx[eiy][j];	
				xx4=cx[eiy][j+1];	
				yy1=my*(sy[eiy]+yspace);
				yy2=my*(sy[eiy]+yspace);
				yy3=my*cy[eiy][j];
				yy4=my*cy[eiy][j+1];
				zz1=zz3=cz[eiy][j];	
				if (j<12) zz1=zz3=sz[eiy][250]; 
				if (j>17) zz1=zz3=sz[eiy][0];
				zz2=zz4=cz[eiy][j+1]; 
				if (j<12) zz2=zz4=sz[eiy][250];	
				if (j>17) zz2=zz4=sz[eiy][0];				
				stlv2(-cw,stl_file,xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3,xx4,yy4,zz4);
			};
	return 0;	
}	

int stl_slicecut(int cw, int my, ofstream &stl_file, int piy)
{
	double x1,x2,x3,x4,z1,z2,z3,z4,y;
	y=(my>0?sy[piy]:-sy[piy]);
	for (int i=0; i<=249; i++)
		{
			x1=sx[piy][250+i];	
			x2=sx[piy][250-i];	
			x3=sx[piy][251+i];
			x4=sx[piy][249-i];
			z1=sz[piy][250+i];	
			z2=sz[piy][250-i];	
			z3=sz[piy][251+i];	
			z4=sz[piy][249-i];
			stlv2(cw,stl_file,x1,y,z1,x2,y,z2,x3,y,z3,x4,y,z4);
		};
	return 0;		
};

int stl_tipclose(int cw, int my, int db, ofstream &stl_file)
{	
	double x1,x2,x3,x4,z1,z2,z3,z4,y1,y2,y3,y4;
	int lmybrk=mybrk;
	cz[lmybrk][13]=(db==1?pmhiz:-pmloz);	
	cz[lmybrk][14]=(db==1?pmhiz:-pmloz); 
	cz[lmybrk][15]=(db==1?pmhiz:-pmloz);
	cz[lmybrk][16]=(db==1?pmhiz:-pmloz);
	for (int tc=0;tc<=13;tc++)
		{
			x1=cx[lmybrk][tc];	
			x2=cx[lmybrk][tc+1];
			x3=cx[lmybrk][28-tc];
			x4=cx[lmybrk][29-tc];
			y1=my*cy[lmybrk][tc];
			y2=my*cy[lmybrk][tc+1];	
			y3=my*cy[lmybrk][28-tc];
			y4=my*cy[lmybrk][29-tc];
			z1=cz[lmybrk][tc];	
			z2=cz[lmybrk][tc+1];
			z3=cz[lmybrk][28-tc];	
			z4=cz[lmybrk][29-tc];
			stlv1(cw,stl_file,x1,y1,z1,x2,y2,z2,x3,y3,z3);
			stlv1(cw,stl_file,x1,y1,z1,x3,y3,z3,x4,y4,z4);
		};	
	return 0;	
}	

int stl_vertical(int cw, int my, ofstream &stl_file, int piy, int src, int dst,int dirz, int yadj)
{
	double linz[251],linx[251],linxu[251],x1,z1,x2,z2,x3,z3,x4,z4,y,posx,posz,dx,dxs,dzrc1,dzrc2,dzrc3,dzcm1,dzcm2,dzcm3,dzmt1,dzmt2,dzmt3;
	int rcidx1,rcidx2,rcidx3,cmidx1,cmidx2,cmidx3,mtidx1,mtidx2,mtidx3,nx;
	
	y=(my>0?sy[piy]+yadj*yspace:-sy[piy]-yadj*yspace); 
	dxs=(sx[piy][0]-sx[piy][250])/200.; 
	rcidx1=(int)(rcjoiner1p*2.);
	rcidx2=(int)(rcjoiner2p*2.);
	rcidx3=(int)(rcjoiner3p*2.);
	cmidx1=(int)(cmjoiner1p*2.);
	cmidx2=(int)(cmjoiner2p*2.);
	cmidx3=(int)(cmjoiner3p*2.);	
	mtidx1=(int)(mtjoiner1p*2.);
	mtidx2=(int)(mtjoiner2p*2.);
	mtidx3=(int)(mtjoiner3p*2.);
	dzrc1=(dirz>0?1.:-1.)*rcjoiner1h/2.;
	dzrc2=(dirz>0?1.:-1.)*rcjoiner2h/2.;
	dzrc3=(dirz>0?1.:-1.)*rcjoiner3h/2.;	
	dzcm1=(dirz>0?1.:-1.)*cmjoiner1h/2.;
	dzcm2=(dirz>0?1.:-1.)*cmjoiner2h/2.;
	dzcm3=(dirz>0?1.:-1.)*cmjoiner3h/2.;	
	dzmt1=(dirz>0?1.:-1.)*mtjoiner1h/2.;
	dzmt2=(dirz>0?1.:-1.)*mtjoiner2h/2.;		dzmt3=(dirz>0?1.:-1.)*mtjoiner3h/2.;
	int rcjoiner=((piy==0)?1:0);
	int cmjoiner=((piy==cmbrk)?1:0);
	int mtjoiner=((piy==mtbrk)?1:0);
	for (int i=0;i<=250;i++) 
	{	
		linx[i]=sx[piy][i];
		linxu[i]=sx[piy][i];
		linz[i]=(sz[piy][i]+sz[piy][500-i])/2;	
	};	
	linz[(int)(hingeprc*2.)]=(linz[(int)(hingeprc*2.)-1]+linz[(int)(hingeprc*2.)+1])/2.;
	
	if (rcjoiner)
		{
			if (rcjoiner1h>0)
				{
					posx=linx[rcidx1];	
					posz=linz[rcidx1];	
					nx=(int)(ceil(rcjoiner1w/dxs/2.)); 
					for (int k=(rcidx1-nx); k<=(rcidx1+nx); k++)
						linz[k]=posz+dzrc1;
					linx[rcidx1-nx-1]=linx[rcidx1-nx]=posx+rcjoiner1w/2.;
					linx[rcidx1+nx+1]=linx[rcidx1+nx]=posx-rcjoiner1w/2.; 
				};
			if (rcjoiner2h>0)
				{	
					posx=linx[rcidx2];	
					posz=linz[rcidx2];
					nx=(int)(ceil(rcjoiner2w/dxs/2.));
					for (int k=(rcidx2-nx); k<=(rcidx2+nx); k++)
						linz[k]=posz+dzrc2;
					linx[rcidx2-nx-1]=linx[rcidx2-nx]=posx+rcjoiner2w/2.;
					linx[rcidx2+nx+1]=linx[rcidx2+nx]=posx-rcjoiner2w/2.;
				};
			if (rcjoiner3h>0)
				{	
					posx=linx[rcidx3];
					posz=linz[rcidx3];	
					nx=(int)(ceil(rcjoiner3w/dxs/2.)); 
					for (int k=(rcidx3-nx); k<=(rcidx3+nx); k++)
						linz[k]=posz+dzrc3;
					linx[rcidx3-nx-1]=linx[rcidx3-nx]=posx+rcjoiner3w/2.;
					linx[rcidx3+nx+1]=linx[rcidx3+nx]=posx-rcjoiner3w/2.;
				};	
		};
	
	if (cmjoiner)
		{
			if (cmjoiner1h>0)
				{
					posx=linx[cmidx1];
					posz=linz[cmidx1];
					nx=(int)(ceil(cmjoiner1w/dxs/2.)); 
					for (int k=(cmidx1-nx); k<=(cmidx1+nx); k++)
						linz[k]=posz+dzcm1;
					linx[cmidx1-nx-1]=linx[cmidx1-nx]=posx+cmjoiner1w/2.;
					linx[cmidx1+nx+1]=linx[cmidx1+nx]=posx-cmjoiner1w/2.; 
				};
			if (cmjoiner2h>0)
				{	
					posx=linx[cmidx2];
					posz=linz[cmidx2];	
					nx=(int)(ceil(cmjoiner2w/dxs/2.));
					for (int k=(cmidx2-nx); k<=(cmidx2+nx); k++)
						linz[k]=posz+dzcm2;
					linx[cmidx2-nx-1]=linx[cmidx2-nx]=posx+cmjoiner2w/2.;
					linx[cmidx2+nx+1]=linx[cmidx2+nx]=posx-cmjoiner2w/2.;
				};
			if (cmjoiner3h>0)
				{	
					posx=linx[cmidx3];
					posz=linz[cmidx3];	
					nx=(int)(ceil(cmjoiner3w/dxs/2.)); 
					for (int k=(cmidx3-nx); k<=(cmidx3+nx); k++)
						linz[k]=posz+dzcm3;
					linx[cmidx3-nx-1]=linx[cmidx3-nx]=posx+cmjoiner3w/2.;
					linx[cmidx3+nx+1]=linx[cmidx3+nx]=posx-cmjoiner3w/2.;
				};	
		};
	
	if (mtjoiner)
		{
			if (mtjoiner1h>0)
				{
					posx=linx[mtidx1];
					posz=linz[mtidx1];	
					nx=(int)(ceil(mtjoiner1w/dxs/2.)); 
					for (int k=(mtidx1-nx); k<=(mtidx1+nx); k++)
						linz[k]=posz+dzmt1;
					linx[mtidx1-nx-1]=linx[mtidx1-nx]=posx+mtjoiner1w/2.;
					linx[mtidx1+nx+1]=linx[mtidx1+nx]=posx-mtjoiner1w/2.;
				};
			if (mtjoiner2h>0)
				{
					posx=linx[mtidx2];
					posz=linz[mtidx2];	
					nx=(int)(ceil(mtjoiner2w/dxs/2.));
					for (int k=(mtidx2-nx); k<=(mtidx2+nx); k++)
						linz[k]=posz+dzmt2;
					linx[mtidx2-nx-1]=linx[mtidx2-nx]=posx+mtjoiner2w/2.;
					linx[mtidx2+nx+1]=linx[mtidx2+nx]=posx-mtjoiner2w/2.;
				};
			if (mtjoiner3h>0)
				{
					posx=linx[mtidx3];
					posz=linz[mtidx3];
					nx=(int)(ceil(mtjoiner3w/dxs/2.));
					for (int k=(mtidx3-nx); k<=(mtidx3+nx); k++)
						linz[k]=posz+dzmt3;
					linx[mtidx3-nx-1]=linx[mtidx3-nx]=posx+mtjoiner3w/2.;
					linx[mtidx3+nx+1]=linx[mtidx3+nx]=posx-mtjoiner3w/2.;
				};	
		};	
	for (int i=0; i<=249; i++)
		{
			switch (src)
				{
					case 0:		
								x1=linx[i];	
								x3=linx[i+1];
								z1=linz[i]; 
								z3=linz[i+1]; 
								break;
					case 1:		
								x1=linxu[i];
								x3=linxu[i+1];
								z1=sz[piy][i];
								z3=sz[piy][i+1]; 
								break;
					case -1:	
								x1=linxu[i];
								x3=linxu[i+1];
								z1=sz[piy][500-i];
								z3=sz[piy][499-i];
								break;
					case -2:	
								x1=linxu[i];
								x3=linxu[i+1];
								z1=-pmloz;
								z3=-pmloz;
								break;
					case 2:		
								x1=linxu[i];
								x3=linxu[i+1];
								z1=+pmhiz;
								z3=+pmhiz;	
								break;
				};
			switch (dst)
				{
					case 0:	
								x2=linx[i];	
								x4=linx[i+1];
								z2=linz[i]; 
								z4=linz[i+1]; 
								break;
					case 1:		
								x2=linxu[i];
								x4=linxu[i+1];
								z2=sz[piy][i];
								z4=sz[piy][i+1];
								break;
					case -1:	
								x2=linxu[i];
								x4=linxu[i+1];
								z2=sz[piy][500-i];
								z4=sz[piy][499-i];
								break;
					case -2:	
								x2=linxu[i];
								x4=linxu[i+1];
								z2=-pmloz;
								z4=-pmloz;	
								break;
					case 2:		
								x2=linxu[i];
								x4=linxu[i+1];
								z2=+pmhiz;
								z4=+pmhiz;
								break;
				};
			stlv2(cw,stl_file,x1,y,z1,x2,y,z2,x3,y,z3,x4,y,z4);
		};
	return 0;	
};

int stl_terminal(int cw, int my, ofstream &stl_file, int piy, int dir, int dirz)
{
	double linx1[251],linx2[251],linz1[251],linz2[251],x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,posx,posz,dxs1,dxs2,dzrc1,dzrc2,dzrc3,dzcm1,dzcm2,dzcm3,dzmt1,dzmt2,dzmt3;
	int rcidx1,rcidx2,rcidx3,cmidx1,cmidx2,cmidx3,mtidx1,mtidx2,mtidx3,nx,piy2;
	
	rcidx1=(int)(rcjoiner1p*2.);
	rcidx2=(int)(rcjoiner2p*2.);
	rcidx3=(int)(rcjoiner3p*2.);	
	cmidx1=(int)(cmjoiner1p*2.);
	cmidx2=(int)(cmjoiner2p*2.);
	cmidx3=(int)(cmjoiner3p*2.);	
	mtidx1=(int)(mtjoiner1p*2.);
	mtidx2=(int)(mtjoiner2p*2.);
	mtidx3=(int)(mtjoiner3p*2.);
	int mys=(my>0?1:-1);
	dxs1=(sx[piy][0]-sx[piy][250])/200.; 
	dxs2=(sx[piy][0]-sx[piy][250])/200.;  
	int rcjoiner=((piy==0)?1:0);
	int cmjoiner=((piy==cmbrk)?1:0);
	int mtjoiner=((piy==mtbrk)?1:0);
	dzrc1=(dirz>0?1.:-1.)*rcjoiner1h/2.;
	dzrc2=(dirz>0?1.:-1.)*rcjoiner2h/2.;
	dzrc3=(dirz>0?1.:-1.)*rcjoiner3h/2.;	
	dzcm1=(dirz>0?1.:-1.)*cmjoiner1h/2.;
	dzcm2=(dirz>0?1.:-1.)*cmjoiner2h/2.;
	dzcm3=(dirz>0?1.:-1.)*cmjoiner3h/2.;	
	dzmt1=(dirz>0?1.:-1.)*mtjoiner1h/2.;
	dzmt2=(dirz>0?1.:-1.)*mtjoiner2h/2.;
	dzmt3=(dirz>0?1.:-1.)*mtjoiner3h/2.;
	for (int i=0;i<=250;i++) 
	{ 	
		linx1[i]=sx[piy][i];
		linx2[i]=sx[piy][i]; 
		linz1[i]=linz2[i]=(sz[piy][i]+sz[piy][500-i])/2; 
	};	
	linz1[(int)(hingeprc*2.)]=(linz1[(int)(hingeprc*2.)-1]+linz1[(int)(hingeprc*2.)+1])/2.; 
	linz2[(int)(hingeprc*2.)]=(linz2[(int)(hingeprc*2.)-1]+linz2[(int)(hingeprc*2.)+1])/2.;

	if (rcjoiner)
		{
			if (rcjoiner1h>0)
				{
					posx=linx1[rcidx1];	
					posz=linz1[rcidx1];	
					nx=(int)(ceil(rcjoiner1w/dxs1/2.));
					for (int k=(rcidx1-nx); k<=(rcidx1+nx); k++) 
						linz1[k]=posz+dzrc1;
					linx1[rcidx1-nx-1]=linx1[rcidx1-nx]=posx+rcjoiner1w/2.;	
					linx1[rcidx1+nx+1]=linx1[rcidx1+nx]=posx-rcjoiner1w/2.;
				};
			if (rcjoiner2h>0)
				{
					posx=linx1[rcidx2];	
					posz=linz1[rcidx2];	
					nx=(int)(ceil(rcjoiner2w/dxs1/2.));
					for (int k=(rcidx2-nx); k<=(rcidx2+nx); k++)
						linz1[k]=posz+dzrc2;
					linx1[rcidx2-nx-1]=linx1[rcidx2-nx]=posx+rcjoiner2w/2.;	
					linx1[rcidx2+nx+1]=linx1[rcidx2+nx]=posx-rcjoiner2w/2.;
				};
			if (rcjoiner3h>0)
				{
					posx=linx1[rcidx3];	
					posz=linz1[rcidx3];	
					nx=(int)(ceil(rcjoiner3w/dxs1/2.)); 
					for (int k=(rcidx3-nx); k<=(rcidx3+nx); k++)
						linz1[k]=posz+dzrc3;
					linx1[rcidx3-nx-1]=linx1[rcidx3-nx]=posx+rcjoiner3w/2.;
					linx1[rcidx3+nx+1]=linx1[rcidx3+nx]=posx-rcjoiner3w/2.;
				};
		};

	if (cmjoiner)
		{
			if (cmjoiner1h>0)
				{
					posx=linx1[cmidx1];	
					posz=linz1[cmidx1];	
					nx=(int)(ceil(cmjoiner1w/dxs1/2.));
					for (int k=(cmidx1-nx); k<=(cmidx1+nx); k++) 
						linz1[k]=posz+dzcm1;
					linx1[cmidx1-nx-1]=linx1[cmidx1-nx]=posx+cmjoiner1w/2.;	
					linx1[cmidx1+nx+1]=linx1[cmidx1+nx]=posx-cmjoiner1w/2.;
				};
			if (cmjoiner2h>0)
				{
					posx=linx1[cmidx2];	
					posz=linz1[cmidx2];	
					nx=(int)(ceil(cmjoiner2w/dxs1/2.));
					for (int k=(cmidx2-nx); k<=(cmidx2+nx); k++) 
						linz1[k]=posz+dzcm2;
					linx1[cmidx2-nx-1]=linx1[cmidx2-nx]=posx+cmjoiner2w/2.;	
					linx1[cmidx2+nx+1]=linx1[cmidx2+nx]=posx-cmjoiner2w/2.;
				};
			if (cmjoiner3h>0)
				{
					posx=linx1[cmidx3];
					posz=linz1[cmidx3];	
					nx=(int)(ceil(cmjoiner3w/dxs1/2.)); 
					for (int k=(cmidx3-nx); k<=(cmidx3+nx); k++)
						linz1[k]=posz+dzcm3;
					linx1[cmidx3-nx-1]=linx1[cmidx3-nx]=posx+cmjoiner3w/2.;	
					linx1[cmidx3+nx+1]=linx1[cmidx3+nx]=posx-cmjoiner3w/2.;
				};
		};
	
	if (mtjoiner)
		{
			if (mtjoiner1h>0)
				{
					posx=linx1[mtidx1];	
					posz=linz1[mtidx1];	
					nx=(int)(ceil(mtjoiner1w/dxs1/2.));
					for (int k=(mtidx1-nx); k<=(mtidx1+nx); k++) 
						linz1[k]=posz+dzmt1;
					linx1[mtidx1-nx-1]=linx1[mtidx1-nx]=posx+mtjoiner1w/2.;	
					linx1[mtidx1+nx+1]=linx1[mtidx1+nx]=posx-mtjoiner1w/2.;
				};
			if (mtjoiner2h>0)
				{
					posx=linx1[mtidx2];	
					posz=linz1[mtidx2];	
					nx=(int)(ceil(mtjoiner2w/dxs1/2.));
					for (int k=(mtidx2-nx); k<=(mtidx2+nx); k++) 
						linz1[k]=posz+dzmt2;
					linx1[mtidx2-nx-1]=linx1[mtidx2-nx]=posx+mtjoiner2w/2.;	
					linx1[mtidx2+nx+1]=linx1[mtidx2+nx]=posx-mtjoiner2w/2.;
				};
			if (mtjoiner3h>0)
				{
					posx=linx1[mtidx3];	
					posz=linz1[mtidx3];	
					nx=(int)(ceil(mtjoiner3w/dxs1/2.));
					for (int k=(mtidx3-nx); k<=(mtidx3+nx); k++) 
						linz1[k]=posz+dzmt3;
					linx1[mtidx3-nx-1]=linx1[mtidx3-nx]=posx+mtjoiner3w/2.;	
					linx1[mtidx3+nx+1]=linx1[mtidx3+nx]=posx-mtjoiner3w/2.;
				};
		};
		
	if (rcjoiner)
		{
			if (rcjoiner1h>0)
				{
					posx=linx2[rcidx1];	
					posz=linz2[rcidx1];	
					nx=(int)(ceil(rcjoiner1w/dxs2/2.));
					for (int k=(rcidx1-nx); k<=(rcidx1+nx); k++) 
						linz2[k]=posz+dzrc1;
					linx2[rcidx1-nx-1]=linx2[rcidx1-nx]=posx+rcjoiner1w/2.;	
					linx2[rcidx1+nx+1]=linx2[rcidx1+nx]=posx-rcjoiner1w/2.;
				};
			if (rcjoiner2h>0)
				{
					posx=linx2[rcidx2];	
					posz=linz2[rcidx2];	
					nx=(int)(ceil(rcjoiner2w/dxs2/2.));
					for (int k=(rcidx2-nx); k<=(rcidx2+nx); k++) 
						linz2[k]=posz+dzrc2;
					linx2[rcidx2-nx-1]=linx2[rcidx2-nx]=posx+rcjoiner2w/2.;	
					linx2[rcidx2+nx+1]=linx2[rcidx2+nx]=posx-rcjoiner2w/2.;
				};
			if (rcjoiner3h>0)
				{
					posx=linx2[rcidx3];	
					posz=linz2[rcidx3];	
					nx=(int)(ceil(rcjoiner3w/dxs2/2.));
					for (int k=(rcidx3-nx); k<=(rcidx3+nx); k++) 
						linz2[k]=posz+dzrc3;
					linx2[rcidx3-nx-1]=linx2[rcidx3-nx]=posx+rcjoiner3w/2.;	
					linx2[rcidx3+nx+1]=linx2[rcidx3+nx]=posx-rcjoiner3w/2.;
				};	
		};
	if (cmjoiner)
		{
			if (cmjoiner1h>0)
				{
					posx=linx2[cmidx1];	
					posz=linz2[cmidx1];	
					nx=(int)(ceil(cmjoiner1w/dxs2/2.));
					for (int k=(cmidx1-nx); k<=(cmidx1+nx); k++) 
						linz2[k]=posz+dzcm1;
					linx2[cmidx1-nx-1]=linx2[cmidx1-nx]=posx+cmjoiner1w/2.;	
					linx2[cmidx1+nx+1]=linx2[cmidx1+nx]=posx-cmjoiner1w/2.;
				};
			if (cmjoiner2h>0)
				{
					posx=linx2[cmidx2];
					posz=linz2[cmidx2];	
					nx=(int)(ceil(cmjoiner2w/dxs2/2.));
					for (int k=(cmidx2-nx); k<=(cmidx2+nx); k++)
						linz2[k]=posz+dzcm2;
					linx2[cmidx2-nx-1]=linx2[cmidx2-nx]=posx+cmjoiner2w/2.;	
					linx2[cmidx2+nx+1]=linx2[cmidx2+nx]=posx-cmjoiner2w/2.;
				};
			if (cmjoiner3h>0)
				{
					posx=linx2[cmidx3];	
					posz=linz2[cmidx3];	
					nx=(int)(ceil(cmjoiner3w/dxs2/2.));
					for (int k=(cmidx3-nx); k<=(cmidx3+nx); k++) 
						linz2[k]=posz+dzcm3;
					linx2[cmidx3-nx-1]=linx2[cmidx3-nx]=posx+cmjoiner3w/2.;	
					linx2[cmidx3+nx+1]=linx2[cmidx3+nx]=posx-cmjoiner3w/2.;
				};	
		};
	if (mtjoiner)
		{	
			if (mtjoiner1h>0)
				{
					posx=linx2[mtidx1];	
					posz=linz2[mtidx1];	
					nx=(int)(ceil(mtjoiner1w/dxs2/2.));
					for (int k=(mtidx1-nx); k<=(mtidx1+nx); k++) 
						linz2[k]=posz+dzmt1;
					linx2[mtidx1-nx-1]=linx2[mtidx1-nx]=posx+mtjoiner1w/2.;	
					linx2[mtidx1+nx+1]=linx2[mtidx1+nx]=posx-mtjoiner1w/2.;
				};
			if (mtjoiner2h>0)
				{
					posx=linx2[mtidx2];	
					posz=linz2[mtidx2];	
					nx=(int)(ceil(mtjoiner2w/dxs2/2.));
					for (int k=(mtidx2-nx); k<=(mtidx2+nx); k++) 
						linz2[k]=posz+dzmt2;
					linx2[mtidx2-nx-1]=linx2[mtidx2-nx]=posx+mtjoiner2w/2.;	
					linx2[mtidx2+nx+1]=linx2[mtidx2+nx]=posx-mtjoiner2w/2.;
				};
			if (mtjoiner3h>0)
				{
					posx=linx2[mtidx3];	
					posz=linz2[mtidx3];	
					nx=(int)(ceil(mtjoiner3w/dxs2/2.));
					for (int k=(mtidx3-nx); k<=(mtidx3+nx); k++) 
						linz2[k]=posz+dzmt3;
					linx2[mtidx3-nx-1]=linx2[mtidx3-nx]=posx+mtjoiner3w/2.;	
					linx2[mtidx3+nx+1]=linx2[mtidx3+nx]=posx-mtjoiner3w/2.;
				};	
		};
	for (int i=0; i<=249; i++)
		{
			x1=linx1[i];
			x2=linx2[i];
			x3=linx1[i+1];	
			x4=linx2[i+1];
			y1=sy[piy]*mys;	
			y2=y1+mys*yspace*dir;
			y3=sy[piy]*mys;	
			y4=y3+mys*yspace*dir;
			z1=linz1[i];
			z2=linz2[i];
			z3=linz1[i+1];
			z4=linz2[i+1];
			stlv2(cw,stl_file,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4);
		};
	return 0;	
};


string dxf_layer (int c)
{
	string layer;
	switch(c)
		{	case 7:
			case 8:
			case 0: layer="lasercut"; 	break;
			case 1: layer="engrave"; 	break;
			case 2: layer="joiner";		break;
			case 3: layer="spoiler";    break;
			case 4: layer="flaps"; 		break;
			case 6: layer="aileron";    break;
			case 9: layer="auxiliar";   break;
		};
	return layer;
}

int dxf_line(ofstream &dxf_file,double x1, double y1, double x2, double y2, int line_color)
{
	dxf_file << "LINE" << endl << "8" << endl << dxf_layer(line_color) << endl << "6" << endl << "CONTINUOUS" << endl;
	dxf_file << "10" << endl << "   " << x1 << endl << "20" << endl << "   " << y1 << endl;
	dxf_file << "11" << endl << "   " << x2 << endl << "21" << endl << "   " << y2 << endl;
	dxf_file << "39" << endl << "0" << endl << "62" << endl << line_color << endl << "0" << endl;
	return 0;
}

int main () 
{
	ifstream def_file("xflrwing.def");
    ifstream cfg_file("xflrwing.cfg");
 	nodigits=2; cout.setf(ios::fixed,ios::floatfield); cout.precision(nodigits);

	cout << "XFLRWING - STL GENERATOR FOR XFLR WINGS " << endl;

	if (def_file.is_open())
		{
			for (i=0;i<500;i++)		
				wp[0][i]=1.000000000;
			for (i=0;i<500;i++) 	
				for (j=1;j<21;j++) 	
					wp[j][i]=0.000000000;
			for (i=1; i<=190; i++) 	
				wp[0][i]=wp[0][i-1]-0.005;			
			for (i=191;i<=230;i++) 	
				wp[0][i]=wp[0][i-1]-0.001; 
			for (i=231;i<=250;i++) 
				wp[0][i]=wp[0][i-1]-0.0005;			
			wp[0][250]=0.0;			
			for (i=251;i<=500;i++)	
				wp[0][i]=wp[0][500-i];
			nr_seg=0;  			
			posy=0.0;

			cout << endl << "STAGE 1:reading wing configuration / definition and airfoils data" << endl << endl<< "Following arameters have non-default values :" << endl;

			//some default values
			wingspan=3600.0;	
			rootspan=254;		
			sf=0.;				
			cmpos=400.;		
			mtpos=1200.;	
			cdih=0.0;		
			mdih=2.0;		
			tdih=7.0;		
			Q=1;
			W=25.;
			pmloz=25.;
			pmhiz=30.;	
			fo=0;
			auw=1;
			wo=1;
			co=1;
			mo=1;	
			po=1;	
			cfweight=80;
			xpsweight=35;
			deadweight=400;	
			resinfiber=2.0;
			xchscale=1.0;
			ywsscale=1.0;
			zhtscale=1.0;
			elliptic=1;	
			xpow=2.4128;
			mr=1;
			flatmodel=0;
			dh=100.;
			sh=6.;				
			hingeprc=30;
			dboxprc=75;	
			sparprc=70;	
			zbase=25.;	
			tts=10.;
			lts=10.;	
			ttd=5.;	
			ltd=5.;	
			led=25.;
			ted=25.;
			yspace=10.;	
			lew=2.;	
			tew=2.0;
			leh=2.0;
			teh=2.0;		
			compfoil=1.;
			compfoill=1.; 
			compfoils=2.;	
			compfoilm=1.;	
			compfoilt=0.;	
			twistsel=-1;
			bs=0;	
			to=0;
			dh=100.;
			sh=6.0;	
			cntsp=0.;	
			tipsp=0.0;
			ctrl1typ=1;	
			ctrl1d1=0.;	
			ctrl1d2=400.;
			ctrl1prc=30.;
			ctrl1dim=0.3;		
			ctrl2typ=1;	
			ctrl2d1=400.;
			ctrl2d2=1200.;
			ctrl2prc=30.;
			ctrl2dim=0.2;		
			ctrl3typ=0.;
			ctrl3d1=1200.;	
			ctrl3d2=2000.;	
			ctrl3prc=30.;	
			ctrl3dim=0.1;		
			srv1d1=145.;
			srv1d2=195.;	
			srv1p1=45.;		
			srv1p2=65.;	
			srv1z=0.0;			
			srv2d1=775.;
			srv2d2=825.;
			srv2p1=42.5;	
			srv2p2=65.;		
			srv2z=0.0;			
			srv3d1=1275.;	
			srv3d2=1325.;	
			srv3p1=40.;		
			srv3p2=65.;		
			srv3z=0.0;			
			rcjoiner1p=72.;	
			rcjoiner1h=0.;	
			rcjoiner1w=4.;	
			rcjoiner2p=34.;	
			rcjoiner2h=0.;
			rcjoiner2w=2.;	
			rcjoiner3p=53.;	
			rcjoiner3h=0.;	
			rcjoiner3w=3.;		
			cmjoiner1p=72.;	
			cmjoiner1h=0.;	
			cmjoiner1w=4.;	
			cmjoiner2p=34.;	
			cmjoiner2h=0.;	
			cmjoiner2w=2.;	
			cmjoiner3p=53.;	
			cmjoiner3h=0.;	
			cmjoiner3w=3.;		
			mtjoiner1p=72.;	
			mtjoiner1h=0.; 	
			mtjoiner1w=4.;	
			mtjoiner2p=34.;	
			mtjoiner2h=0.; 	
			mtjoiner2w=2.;	
			mtjoiner3p=53.;	
			mtjoiner3h=0.;	
			mtjoiner3w=3.;		
			nodigits=3;			
			//now let's read user defined values
			if (cfg_file.is_open())
				{
					while (getline (cfg_file, txt_line)) 
						{
							cfg_file >> cfg_par >> cfg_val;
							if (cfg_par=="WINGSPAN") 	
								{ cout << "WINGSPAN " 	<< cfg_val << endl; 	wingspan=cfg_val; }
							if (cfg_par=="ROOTSPAN") 	
								{ cout << "ROOTSPAN " 	<< cfg_val << endl; 	rootspan=cfg_val; }
							if (cfg_par=="STABFIN") 	
								{ cout << "STABFIN " 	<< cfg_val << endl; 	sf=cfg_val; }	
							if (cfg_par=="MOLDRECT") 	
								{ cout << "MOLDRECT " 	<< cfg_val << endl; 	mr=(int)cfg_val; }	
							if (cfg_par=="FLATMODEL") 	
								{ cout << "FLATMODEL " 	<< cfg_val << endl; 	flatmodel=(int)cfg_val; }	
							if (cfg_par=="CNTMIDPOS") 	
								{ cout << "CNTMIDPOS " 	<< cfg_val << endl; 	cmpos=cfg_val; }
							if (cfg_par=="MIDTIPPOS")
								{ cout << "MIDTIPPOS " 	<< cfg_val << endl; 	mtpos=cfg_val; }
							if (cfg_par=="CNTANGLE") 
								{ cout << "CNTANGLE "	<< cfg_val << endl; 	cdih=cfg_val; }
							if (cfg_par=="MIDANGLE") 
								{ cout << "MIDANGLE "	<< cfg_val << endl; 	mdih=cfg_val; }
							if (cfg_par=="TIPANGLE") 
								{ cout << "TIPANGLE "	<< cfg_val << endl; 	tdih=cfg_val; }
							if (cfg_par=="CHORDRES") 
								{ cout << "CHORDRES " 	<< cfg_val << endl; 	Q=(int) cfg_val; }
							if (cfg_par=="WSPANRES") 
								{ cout << "WSPANRES " 	<< cfg_val << endl; 	W=cfg_val; }
							if (cfg_par=="FULLOUTPUT")
								{ cout << "FULOUTPUT " 	<< cfg_val << endl; 	fo=(int)cfg_val; }
							if (cfg_par=="AUWOUTPUT") 
								{ cout << "AUWOUTPUT " 	<< cfg_val << endl; 	auw=(int)cfg_val; }
							if (cfg_par=="WINGOUTPUT") 
								{ cout << "WINGOUTPUT " << cfg_val << endl; 	wo=(int)cfg_val; }
							if (cfg_par=="COREOUTPUT") 
								{ cout << "COREOUTPUT " << cfg_val << endl; 	co=(int)cfg_val; }
							if (cfg_par=="MOLDOUTPUT")
								{ cout << "MOLDOUTPUT " << cfg_val << endl; 	mo=(int)cfg_val; }
							if (cfg_par=="PLUGOUTPUT")
								{ cout << "PLUGOUTPUT " << cfg_val << endl; 	po=(int)cfg_val; }
							if (cfg_par=="TESTOUTPUT")
								{ cout << "TESTOUTPUT " << cfg_val << endl; 	to=(int)cfg_val; }
							if (cfg_par=="CFWEIGHT") 
								{ cout << "CFWEIGHT " 	<< cfg_val << endl; 	cfweight=(int)cfg_val; }
							if (cfg_par=="XPSWEIGHT") 
								{ cout << "XPSWEIGHT " 	<< cfg_val << endl; 	xpsweight=(int)cfg_val; }
							if (cfg_par=="DEADWEIGHT") 
								{ cout << "DEADWEIGHT " << cfg_val << endl; 	deadweight=(int)cfg_val; }
							if (cfg_par=="RESINFIBER")
								{ cout << "RESINFIBER " << cfg_val << endl; 	resinfiber=(int)cfg_val; }
							if (cfg_par=="XCHSCALE")
								{ cout << "XCHSCALE " 	<< cfg_val << endl; 	xchscale=cfg_val; }	
							if (cfg_par=="YWSSCALE") 
								{ cout << "YWSSCALE " 	<< cfg_val << endl; 	ywsscale=cfg_val; }	
							if (cfg_par=="ZHTSCALE") 
								{ cout << "ZHTSCALE " 	<< cfg_val << endl; 	zhtscale=cfg_val; }	
							if (cfg_par=="ELLIPTIC") 
								{ cout << "ELLIPTIC " 	<< cfg_val << endl; 	elliptic=(int)cfg_val; }
							if (cfg_par=="XPOWER") 
								{ cout << "XPOWER " 	<< cfg_val << endl; 	xpow=cfg_val; }
							if (cfg_par=="DISTHOLES") 
								{ cout << "DISTHOLES " 	<< cfg_val << endl; 	dh=cfg_val; }	
							if (cfg_par=="SIZEHOLES") 
								{ cout << "SIZEHOLES " 	<< cfg_val << endl; 	sh=cfg_val; }	
							if (cfg_par=="HINGEPRC") 
								{ cout << "HINGEPRC " 	<< cfg_val << endl; 	hingeprc=(int)cfg_val; ctrlprc=ctrl1prc=ctrl2prc=ctrl3prc=hingeprc;}
							if (cfg_par=="CTRLPRC") 
								{ cout << "CTRLPRC " 	<< cfg_val << endl; 	hingeprc=(int)cfg_val; ctrlprc=ctrl1prc=ctrl2prc=ctrl3prc=hingeprc;}
							if (cfg_par=="CTRLDIM") 
								{ cout << "CTRLDIM " 	<< cfg_val << endl; 	ctrldim=ctrl1dim=ctrl2dim=ctrl3dim=cfg_val;}
							if (cfg_par=="DBOXPRC") 
								{ cout << "DBOXPRC " 	<< cfg_val << endl; 	dboxprc=(int)cfg_val; /* dboxprc=(int)(200-2*cfg_val);*/} 
							if (cfg_par=="SPARPRC") 
								{ cout << "SPARPRC " 	<< cfg_val << endl; 	sparprc=(int)cfg_val;  /*(int)(200-2*cfg_val);*/} 
							if (cfg_par=="MOLDLOSIZE")	
								{ cout << "MOLDLOSIZE " << cfg_val << endl; 	pmloz = cfg_val; }
							if (cfg_par=="MOLDHISIZE")
								{ cout << "MOLDHISIZE " << cfg_val << endl; 	pmhiz = cfg_val; }
							if (cfg_par=="TETRENCH") 
								{ cout << "TETRENCH " 	<< cfg_val << endl; 	tts=cfg_val; }
							if (cfg_par=="LETRENCH") 
								{ cout << "LETRENCH " 	<< cfg_val << endl; 	lts=cfg_val; }
							if (cfg_par=="TRENCH") 	
								{ cout << "TRENCH " 	<< cfg_val << endl; 	tts=lts=cfg_val; }
							if (cfg_par=="TETRDIST") 
								{ cout << "TETRDIST " 	<< cfg_val << endl; 	ttd=cfg_val;}
							if (cfg_par=="LETRDIST") 
								{ cout << "LETRDIST " 	<< cfg_val << endl; 	ltd=cfg_val; }
							if (cfg_par=="TRENCHDIST") 
								{ cout << "TRENCHDIST " << cfg_val << endl; 	ttd=ltd=cfg_val; }
							if (cfg_par=="LESPACE") 
								{ cout << "LESPACE "	<< cfg_val << endl; 	led=cfg_val; }
							if (cfg_par=="TESPACE") 
								{ cout << "TESPACE "	<< cfg_val << endl; 	ted=cfg_val; }	
							if (cfg_par=="XSPACE") 	
								{ cout << "XSPACE "		<< cfg_val << endl; 	led=ted=cfg_val; }
							if (cfg_par=="YSPACE") 	
								{ cout << "YSPACE "		<< cfg_val << endl; 	yspace=cfg_val; }
							if (cfg_par=="LEWALLW")	
								{ cout << "LEWALLW " 	<< cfg_val << endl; 	lew=cfg_val; }	
							if (cfg_par=="TEWALLW")	
								{ cout << "TEWALLW "	<< cfg_val << endl; 	tew=cfg_val; }	
							if (cfg_par=="WALLW")	
								{ cout << "WALLW " 		<< cfg_val << endl; 	lew=tew=cfg_val; }	
							if (cfg_par=="LEWALLH")	
								{ cout << "LEWALLH " 	<< cfg_val << endl; 	leh=cfg_val; }	
							if (cfg_par=="TEWALLH")	
								{ cout << "TEWALLH "	<< cfg_val << endl; 	teh=cfg_val; }	
							if (cfg_par=="WALLH")	
								{ cout << "WALLH " 		<< cfg_val << endl; 	leh=teh=cfg_val; }	
							if (cfg_par=="OFFSET") 	
								{ cout << "OFFSET " 	<< cfg_val << endl; 	compfoil=cfg_val; }
							if (cfg_par=="OFFSETL") 
								{ cout << "OFFSETL " 	<< cfg_val << endl; 	compfoill=cfg_val; }
							if (cfg_par=="OFFSETS") 
								{ cout << "OFFSETS " 	<< cfg_val << endl; 	compfoils=cfg_val; }
							if (cfg_par=="OFFSETM") 
								{ cout << "OFFSETM " 	<< cfg_val << endl; 	compfoilm=cfg_val; }
							if (cfg_par=="OFFSETT") 
								{ cout << "OFFSETT " 	<< cfg_val << endl; 	compfoilt=cfg_val; }
							if (cfg_par=="TWISTSEL") 
								{ cout << "TWISTSEL " 	<< cfg_val << endl; 	twistsel=cfg_val; }	
							if (cfg_par=="BIGSPACE") 
								{ cout << "BIGSPACE " 	<< cfg_val << endl; 	bs=(int)cfg_val; }	
							if (cfg_par=="CTRL1TYPE") 
								{ cout << "CTRL1T " 	<< cfg_val << endl; 	ctrl1typ=(int)cfg_val; }
							if (cfg_par=="CTRL1D1") 
								{ cout << "CTRL1D1 " 	<< cfg_val << endl; 	ctrl1d1=cfg_val; }
							if (cfg_par=="CTRL1D2") 
								{ cout << "CTRL1D2 " 	<< cfg_val << endl; 	ctrl1d2=cfg_val;}
							if (cfg_par=="CTRL1PRC") 
								{ cout << "CTRL1PRC " 	<< cfg_val << endl; 	ctrl1prc=cfg_val; }
							if (cfg_par=="CTRL1DIM") 
								{ cout << "CTRL1DIM " 	<< cfg_val << endl; 	ctrl1dim=cfg_val; }
							if (cfg_par=="CTRL2TYPE") 
								{ cout << "CTRL2TYPE " 	<< cfg_val << endl; 	ctrl2typ=(int)cfg_val; }
							if (cfg_par=="CTRL2D1") 
								{ cout << "CTRL2D1 " 	<< cfg_val << endl; 	ctrl2d1=cfg_val; }
							if (cfg_par=="CTRL2D2") 
								{ cout << "CTRL2D2 " 	<< cfg_val << endl; 	ctrl2d2=cfg_val;}
							if (cfg_par=="CTRL2PRC") 
								{ cout << "CTRL2PRC " 	<< cfg_val << endl; 	ctrl2prc=cfg_val; }
							if (cfg_par=="CTRL2DIM") 
								{ cout << "CTRL2DIM " 	<< cfg_val << endl; 	ctrl2dim=cfg_val; }
							if (cfg_par=="CTRL3TYPE") 
								{ cout << "CTRL3TYPE " 	<< cfg_val << endl; 	ctrl3typ=(int)cfg_val; }
							if (cfg_par=="CTRL3D1") 
								{ cout << "CTRL3D1 " 	<< cfg_val << endl; 	ctrl3d1=cfg_val; }
							if (cfg_par=="CTRL3D2") 	
								{ cout << "CTRL3D2 " 	<< cfg_val << endl; 	ctrl3d2=cfg_val;}
							if (cfg_par=="CTRL3PRC") 
								{ cout << "CTRL3PRC " 	<< cfg_val << endl; 	ctrl3prc=cfg_val; }
							if (cfg_par=="CTRL3DIM") 
								{ cout << "CTRL3DIM " 	<< cfg_val << endl; 	ctrl3dim=cfg_val; }
							if (cfg_par=="SERVO1D1") 
								{ cout << "SERVO1D1 "	<< cfg_val << endl; 	srv1d1=cfg_val; }
							if (cfg_par=="SERVO1D2") 
								{ cout << "SERVO1D2 "	<< cfg_val << endl; 	srv1d2=cfg_val; }
							if (cfg_par=="SERVO1P1") 
								{ cout << "SERVO1P1 "	<< cfg_val << endl; 	srv1p1=cfg_val; }
							if (cfg_par=="SERVO1P2") 
								{ cout << "SERVO1P2 "	<< cfg_val << endl; 	srv1p2=cfg_val; }
							if (cfg_par=="SERVO1Z") 
								{ cout << "SERVO1Z "	<< cfg_val << endl; 	srv1z=cfg_val; }
							if (cfg_par=="SERVO2D1") 
								{ cout << "SERVO2D1 "	<< cfg_val << endl; 	srv2d1=cfg_val; }
							if (cfg_par=="SERVO2D2") 
								{ cout << "SERVO2D2 "	<< cfg_val << endl; 	srv2d2=cfg_val; }
							if (cfg_par=="SERVO2P1") 
								{ cout << "SERVO2P1 "	<< cfg_val << endl; 	srv2p1=cfg_val; }
							if (cfg_par=="SERVO2P2") 
								{ cout << "SERVO2P2 "	<< cfg_val << endl; 	srv2p2=cfg_val; }
							if (cfg_par=="SERVO2Z") 
								{ cout << "SERVO2Z "	<< cfg_val << endl; 	srv2z=cfg_val; }
							if (cfg_par=="SERVO3D1") 
								{ cout << "SERVO3D1 "	<< cfg_val << endl; 	srv3d1=cfg_val; }
							if (cfg_par=="SERVO3D2") 
								{ cout << "SERVO3D2 "	<< cfg_val << endl; 	srv3d2=cfg_val; }
							if (cfg_par=="SERVO3P1") 
								{ cout << "SERVO3P1 "	<< cfg_val << endl; 	srv3p1=cfg_val; }
							if (cfg_par=="SERVO3P2") 
								{ cout << "SERVO3P2 "	<< cfg_val << endl; 	srv3p2=cfg_val; }
							if (cfg_par=="SERVO3Z") 
								{ cout << "SERVO3Z "	<< cfg_val << endl; 	srv3z=cfg_val; }
							if (cfg_par=="RCJOINER1P") 
								{ cout << "RCJOINER1P " << cfg_val << endl; 	rcjoiner1p=cfg_val; }
							if (cfg_par=="RCJOINER1H") 
								{ cout << "RCJOINER1H " << cfg_val << endl; 	rcjoiner1h=cfg_val; }
							if (cfg_par=="RCJOINER1W") 
								{ cout << "RCJOINER1W " << cfg_val << endl; 	rcjoiner1w=cfg_val; }
							if (cfg_par=="RCJOINER2P") 
								{ cout << "RCJOINER2P " << cfg_val << endl; 	rcjoiner2p=cfg_val; }
							if (cfg_par=="RCJOINER2H") 
								{ cout << "RCJOINER2H " << cfg_val << endl; 	rcjoiner2h=cfg_val; }
							if (cfg_par=="RCJOINER2W") 
								{ cout << "RCJOINER2W " << cfg_val << endl; 	rcjoiner2w=cfg_val; }
							if (cfg_par=="RCJOINER3P") 
								{ cout << "RCJOINER3P " << cfg_val << endl; 	rcjoiner3p=cfg_val; }
							if (cfg_par=="RCJOINER3H") 
								{ cout << "RCJOINER3H " << cfg_val << endl; 	rcjoiner3h=cfg_val; }
							if (cfg_par=="RCJOINER3W") 
								{ cout << "RCJOINER3W " << cfg_val << endl; 	rcjoiner3w=cfg_val; }
							if (cfg_par=="CMJOINER1P") 
								{ cout << "CMJOINER1P " << cfg_val << endl; 	cmjoiner1p=cfg_val; }
							if (cfg_par=="CMJOINER1H") 
								{ cout << "CMJOINER1H " << cfg_val << endl; 	cmjoiner1h=cfg_val; }
							if (cfg_par=="CMJOINER1W") 
								{ cout << "CMJOINER1W " << cfg_val << endl; 	cmjoiner1w=cfg_val; }
							if (cfg_par=="CMJOINER2P") 
								{ cout << "CMJOINER2P " << cfg_val << endl; 	cmjoiner2p=cfg_val; }
							if (cfg_par=="CMJOINER2H") 
								{ cout << "CMJOINER2H " << cfg_val << endl; 	cmjoiner2h=cfg_val; }
							if (cfg_par=="CMJOINER2W") 
								{ cout << "CMJOINER2W " << cfg_val << endl; 	cmjoiner2w=cfg_val; }
							if (cfg_par=="CMJOINER3P") 
								{ cout << "CMJOINER3P " << cfg_val << endl; 	cmjoiner3p=cfg_val; }
							if (cfg_par=="CMJOINER3H") 
								{ cout << "CMJOINER3H " << cfg_val << endl; 	cmjoiner3h=cfg_val; }
							if (cfg_par=="CMJOINER3W") 
								{ cout << "CMJOINER3W " << cfg_val << endl; 	cmjoiner3w=cfg_val; }
							if (cfg_par=="MTJOINER1P") 
								{ cout << "MTJOINER1P " << cfg_val << endl; 	mtjoiner1p=cfg_val; }
							if (cfg_par=="MTJOINER1H") 
								{ cout << "MTJOINER1H " << cfg_val << endl; 	mtjoiner1h=cfg_val; }
							if (cfg_par=="MTJOINER1W") 
								{ cout << "MTJOINER1W " << cfg_val << endl; 	mtjoiner1w=cfg_val; }
							if (cfg_par=="MTJOINER2H") 
								{ cout << "MTJOINER2H " << cfg_val << endl; 	mtjoiner2h=cfg_val; }
							if (cfg_par=="MTJOINER2W") 
								{ cout << "MTJOINER2W " << cfg_val << endl; 	mtjoiner2w=cfg_val; }
							if (cfg_par=="MTJOINER3P") 
								{ cout << "MTJOINER3P " << cfg_val << endl; 	mtjoiner3p=cfg_val; }
							if (cfg_par=="MTJOINER3H") 
								{ cout << "MTJOINER3H " << cfg_val << endl; 	mtjoiner3h=cfg_val; }
							if (cfg_par=="MTJOINER3W") 
								{ cout << "MTJOINER3W " << cfg_val << endl; 	mtjoiner3w=cfg_val; }
							if (cfg_par=="NODIGITS") 
								{ cout << "NODIGITS " 	<< cfg_val << endl; 	nodigits=(int)cfg_val; cout.setf(ios::fixed,ios::floatfield); cout.precision(nodigits);}
						};
				};	
			if (fo==1) 
				wo=co=mo=po=1;	
			W=(W<2.?2.:W);
			Q=(Q>2?2:Q);	
			cfg_file.close(); 

			while (getline (def_file, txt_line)) 
				{
					def_file >> posy >> chord >> offset >> dihedral >> twist >> unk1 >> unk2 >> unk3 >> unk4 >> airfoilt >> airfoilx;
					aposy[nr_seg]=posy*ywsscale;
					achord[nr_seg]=chord*xchscale; 
					aoffset[nr_seg]=offset*xchscale; 
					adihedral[nr_seg]=dihedral; 
					atwist[nr_seg]=twist;	
					airfoil[nr_seg]=airfoilt;
					adjust[nr_seg]=airfoilx;
					nr_seg++;
				};
			cout << "Wing definition has " << --nr_seg << " sections" << endl;
			std::vector<double> VPY(nr_seg),VCH(nr_seg),VXO(nr_seg),VDA(nr_seg),VWA(nr_seg),VZO(nr_seg),VPR(nr_seg),VSP(nr_seg);
			def_file.close();
			ymax=(wingspan==0?posy*1000.0*ywsscale:wingspan/2.0);
			aoffsetd[0]=0.0;
			yscale=wingspan/2000.0/aposy[nr_seg-1];
			xscale=rootspan/1000.0/achord[0];
			wingmax=ymax;
			ww=ymax*1.0;
			rr=(rootspan==0?achord[0]*1000.0*xchscale:rootspan);
			cout << "xscale " << xscale << " yscale " << yscale << endl;
			for (i=0;i<nr_seg;i++)
				{
					aposy[i]=aposy[i]*1000.0*yscale;	
					VPY[i]=aposy[i];	
					achord[i]=achord[i]*1000.0*xscale;	
					VCH[i]=achord[i];	
					aoffset[i]=aoffset[i]*1000.0*xscale;	
					VXO[i]=aoffset[i];
					VSP[i]=VXO[i]+0.25*VCH[i];	
					VDA[i]=adihedral[i];	
					VWA[i]=atwist[i];		
					VZO[i]=(i==0?0.0:VZO[i-1]+(VPY[i]-VPY[i-1])*sin(PI*VDA[i-1]/180));
					cout << "Y " << VPY[i] << " CH " << VCH[i] << " XO " << VXO[i] << " DA " << VDA[i] << " WA " << VWA[i] << " ZO " << VZO[i] << " airfoil " << airfoil[i] << endl;
					file_name=airfoil[i]+".dat";	
					ifstream dat_file(file_name.c_str());
					nrl=0;	
					cout << " reading " << file_name << " adjusted with " << adjust[nrl] ;
					idxle[i]=idxte[i]=idxu[i]=idxl[i]=0;
					while (getline(dat_file,txt_line))
						{ 
							dat_file>>x>>y; 
							ax[i][nrl]=x; 
							ay[i][nrl]=y*zhtscale*adjust[i]; 
							if ((nrl>0)&&(ax[i][nrl]<ax[i][nrl-1])) idxu[i]++; nrl++; 
						};
					dat_file.close();	
					cout << "  ... " << nrl << " lines " << endl;
					adat[i]=nrl-1;
				};

			cout << endl << "STAGE 2:computing high-resolution internal airfoil model" << endl << endl;

			for (int i=0;i<nr_seg;i++) 
				{
					for (j=0;j<250;j++)	
					{
						k=0;
						while(!between(wp[0][j],ax[i][k],ax[i][k+1])) 
							k++; 
						wp[i+1][j]=ay[i][k]+(ay[i][k+1]-ay[i][k])/(ax[i][k+1]-ax[i][k])*(wp[0][j]-ax[i][k]);
					};
					last_k=k+1;
					k=0;
					minim=ax[i][k];	
					while (ax[i][k+1]<minim)	
						minim=ax[i][++k];
					wp[i+1][250]=ay[i][k];
					for (j=250;j<=500;j++) 
					{
						k=last_k;
						while(!between(wp[0][j],ax[i][k],ax[i][k+1])) 
							k++;
						wp[i+1][j]=ay[i][k]+(ay[i][k+1]-ay[i][k])/(ax[i][k+1]-ax[i][k])*(wp[0][j]-ax[i][k]);
					};
                };	
			if (sf!=0) 
				for (int i=0;i<nr_seg;i++)
					for (j=0;j<=250;j++) 
						{
							wp[i+1][j]=(wp[i+1][j]-wp[i+1][500-j])*sf/2.; 
							wp[i+1][500-j]=-wp[i+1][j]; 
						};
			for (int i=1;i<=nr_seg;i++)
				{
					maximh=0.0; 
					for (int j=0;j<=250; j++) 
						maximh=(((wp[i][j]-wp[i][500-j])>maximh)?(wp[i][j]-wp[i][500-j]):maximh); 
					camber=0.0;	
					for (int j=0;j<=250; j++) 
						camber=((((wp[i][j]+wp[i][500-j])/2.)>camber)?(wp[i][j]+wp[i][500-j])/2.:camber);
					cout << "Section " << i << " airfoil " << airfoil[i-1] << "  " << maximh*100.0 << "%  camber " << camber*100.0 << "%" << endl;
				};	
			posmax=1;
			maxim=0;

			cout << endl << "STAGE 3: computing spline interpolation for chord, offsets, dihedral, washout and profiles";

			ymaxnew=(elliptic==1?wingspan/2.0:ymax);
			ymax=ymaxnew;
			ssp.set_points(VPY,VSP); 
			sch.set_points(VPY,VCH); 
			sxo.set_points(VPY,VXO); 
			szo.set_points(VPY,VZO); 
			sda.set_points(VPY,VDA); 
			swa.set_points(VPY,VWA);	
			for (int i=0;i<501;i++) 
				for (int j=0;j<nr_seg;j++)
					{ 	
						VPR[j]=wp[j+1][i];
						spr[i].set_points(VPY,VPR);	
					};

			cout << endl << endl << "STAGE 4:computing space coordinates for wing surface" << endl << endl;

			sy[0]=cmpos; 
			sy[1]=mtpos; 
			sy[2]=ymaxnew; 
			maxidx=2;
			for (posy=0.0;					posy<(ymaxnew-5.*W);		posy=posy+W) 
				asy_add(posy);
			for (posy=(ymaxnew-5.*W);		posy<(ymaxnew-5.*W/5.);		posy=posy+W/5.) 
				asy_add(posy);
			for (posy=(ymaxnew-5.*W/5.);	posy<(ymaxnew-5.*W/25.);	posy=posy+W/25.)
				asy_add(posy);
			for (posy=(ymaxnew-5.*W/25.);	posy<(ymaxnew-5.*W/100.);	posy=posy+W/100.)
				asy_add(posy);
			for (posy=(ymaxnew-5.*W/100.);	posy<ymaxnew;				posy=posy+W/500.)
				asy_add(posy);

			cout << "Initial " << maxidx << endl;	for (int i=0;i<=maxidx;i++) cout << i << "@"<<sy[i] << "|"; cout << endl;

			asy_add(srv1d1-0.2);	
			asy_add(srv1d1-0.1); 
			asy_add(srv1d1);	
			asy_add(srv1d1+0.1);	
			asy_add(srv1d1+0.2);
			
			asy_add(srv1d2-0.2);
			asy_add(srv1d2-0.1); 
			asy_add(srv1d2);
			asy_add(srv1d2+0.1);
			asy_add(srv1d2+0.2);
			
			asy_add(srv2d1-0.2);
			asy_add(srv2d1-0.1); 
			asy_add(srv2d1);
			asy_add(srv2d1+0.1);
			asy_add(srv2d1+0.2);	
			
			asy_add(srv2d2-0.2);
			asy_add(srv2d2-0.1); 
			asy_add(srv2d2);
			asy_add(srv2d2+0.1);
			asy_add(srv2d2+0.2);
			
			asy_add(srv3d1-0.2);
			asy_add(srv3d1-0.1); 
			asy_add(srv3d1);	
			asy_add(srv3d1+0.1);
			asy_add(srv3d1+0.2);		
			asy_add(srv3d2-0.2);	
			asy_add(srv3d2-0.1); 
			asy_add(srv3d2);	
			asy_add(srv3d2+0.1);	
			asy_add(srv3d2+0.2);	

			asy_add(ctrl1d1);	
			asy_add(ctrl1d1+0.1);
			asy_add(ctrl1d2-0.1); 
			asy_add(ctrl1d2);	
			asy_add(ctrl1d2+0.1);	
			
			asy_add(ctrl2d1);
			asy_add(ctrl2d1+0.1);
			asy_add(ctrl2d2-0.1); 
			asy_add(ctrl2d2);
			asy_add(ctrl2d2+0.1);
			
			asy_add(ctrl3d1);
			asy_add(ctrl3d1+0.1);
			asy_add(ctrl3d2-0.1); 
			asy_add(ctrl3d2);	
			asy_add(ctrl3d2+0.1);
			
			for (float f=0.;f<ymaxnew;f=f+dh) 
				{
					if (f>0) 
						{
							asy_add(f-sh/2.);
							asy_add(f-sh/2.+0.5);
						};
					asy_add(f+sh/2.-0.5); 
					asy_add(f+sh/2.);
				};
			maxidx--;
			asy_sort();		
			cout << "Sorted " << asy_sort() << endl;
			for (int i=0;i<=maxidx;i++) 
				cout << i << "@"<<sy[i] << "|"; cout << endl;
			asy_compact();
			cout << "Points of Interest " << asy_compact() << endl;	
			for (int i=0;i<=maxidx;i++) 
				cout << i << "@"<<sy[i] << "|"; cout << endl;
			cmbrk=asy_find(cmpos);
			cout << " Central / mid cut " << cmpos << " found @ " << cmbrk << endl;//cmbrk=(int)(cmpos/W);
			mtbrk=asy_find(mtpos);	
			cout << " Mid / tip cut "<< mtpos <<" found @ " << mtbrk << endl;//mtbrk=(int)(mtpos/W);
			mybrk=asy_find(ymax);   
			cout << " Tip " << ymax << " found @ " << mybrk << endl; //
			for (int i=0; i<=mybrk; i++)
				{	
					xpower=2.0+(2.0*sy[i]/wingspan)*(xpow-2.0);
					posy=sy[i]; 
					chord=(elliptic==1?2.0*rootspan*pow(pow(0.5,xpower)-pow(posy/wingspan,xpower),1/xpower)+0.01:sch(posy))+0.001;
					offset=(elliptic==1?(-(1-hingeprc/100.0)*2*rootspan*pow(pow(0.5,xpower)-pow(posy/wingspan,xpower),1/xpower)+(1-hingeprc/100.0)*rootspan):sxo(posy));
					offsetd=((flatmodel==1)?0.:((elliptic==1)?((i<=cmbrk)?(posy*sin(cdih*PI/180.)):((i<=mtbrk)?(cmpos*sin(cdih*PI/180.)+(posy-cmpos)*sin(mdih*PI/180.)):(cmpos*sin(cdih*PI/180.)+(mtpos-cmpos)*sin(mdih*PI/180.)+(posy-mtpos)*sin(tdih*PI/180.)))):szo(posy)));
					twist=(sf==1?0:swa(posy));
					switch(twistsel)
						{
							case -1:	
										xtrans=wp[0][(int)(hingeprc*2.0)]; 
										ytrans=(spr[(int)(hingeprc*2.0)](posy)+spr[(int)(500-hingeprc*2.0)](posy))/2.0;	
										break;
							case 0: 	
										xtrans=0.0; 
										ytrans=0.0;	
										break;
							case 1: 	
										xtrans=wp[0][250]; 	
										ytrans=spr[250](posy);	
										break;
							case 2: 	
										xtrans=wp[0][150]; 
										ytrans=spr[150](posy); 
										break;
						};
					for (k=0;k<=500;k++)
						{ 	
							x=wp[0][k]-xtrans;
							y=spr[k](posy)-ytrans;
							r=sqrt(x*x+y*y); 
							fi=atan2(y,x)*180/PI-twist;
							nx=r*cos(fi*PI/180)+xtrans;	
							ny=r*sin(fi*PI/180)+ytrans;	
							sxe[i][k]=nx*chord+offset;
							sze[i][k]=ny*chord+offsetd;
						}; 
				 };
			cout << "... ready!" << endl;
			posydx=0.0; 
			csurface=0.0;
			for (int i=0;i<mybrk;i++)
				{
					ps[i]=0.0; 
					sv[i]=0.0;
					csurface=csurface+(sxe[i][0]-sxe[i][250]+sxe[i+1][0]-sxe[i+1][250])*(sy[i+1]-sy[i])/10000.;
					for (k=0;k<=249;k++)
						{
							px1=sxe[i][k]; 
							px2=sxe[i][k+1];
							px3=sxe[i][499-k]; 
							px4=sxe[i][500-k]; 
							py1=sze[i][k]; 
							py2=sze[i][k+1]; 
							py3=sze[i][499-k];	
							py4=sze[i][500-k];
							ps[i]=ps[i]+sqrt((px1-px2)*(px1-px2)+(py1-py2)*(py1-py2))+sqrt((px4-px3)*(px4-px3)+(py4-py3)*(py4-py3));
							la=sqrt((px1-px2)*(px1-px2)+(py1-py2)*(py1-py2));
							lb=sqrt((px3-px2)*(px3-px2)+(py3-py2)*(py3-py2));	
							lc=sqrt((px4-px3)*(px4-px3)+(py4-py3)*(py4-py3));
							ld=sqrt((px1-px4)*(px1-px4)+(py1-py4)*(py1-py4));
							le=sqrt((px4-px2)*(px4-px2)+(py4-py2)*(py4-py2));	
							lp=(lb+lc+le)/2.0;	
							sv[i]=sv[i]+sqrt(lp*(lp-lb)*(lp-lc)*(lp-le));	
							lp=(la+ld+le)/2.0;	
							sv[i]=sv[i]+sqrt(lp*(lp-la)*(lp-ld)*(lp-le));
						}
					cout << "slice " << i << " Position " << sy[i] << " Perimeter " << ps[i] << " mm   Surface " << sv[i] << " mm2 " << csurface << endl;
				};
			surface=0.0; 
			volume=0.0; 
			for (int i=0; i<maxidx; i++) 
				{ 
					surface=surface+(ps[i]+ps[i+1])*(sy[i+1]-sy[i])/1000000.; 
					volume=volume+(sv[i]+sv[i+1])*(sy[i+1]-sy[i])/1000.; 
				};
			cout << endl << "Wingspan " << wingspan << "mm Root chord " << rootspan << " mm  Aspect Ratio " << wingspan*wingspan/csurface/10000.0 << endl ;
			cout << "Wing surface " << csurface << " dm2 Total surface " << surface << " mp Volume " << volume << " cm3" << endl;	
			if (auw==1)
				{
					estwing=volume*xpsweight/1000.0 + surface*cfweight*resinfiber+wingspan/25;  
					cout << "Estimated wing weight " << estwing << " grams " << endl;
					estfin=surface*cfweight*resinfiber*0.06 +10; 
					cout << "Estimated fin weight " << estfin << " grams" << endl;
					estelev=surface*cfweight*resinfiber*0.11+10;	
					cout << "Estimated elevator weight " << estelev << " grams" << endl;
					estfuse=wingspan/20.0-20.0;	
					cout << "Estimated fuselage weight " << estfuse << " grams" << endl;
					estauw=estwing+estfin+estelev+estfuse+deadweight;
					cout << "Estimated AUW " << estauw << " grams" << endl;  
					cout << "Estimated wing load " << estauw/csurface/1.11 << "-" << estauw/csurface << "gr/dmp"  << endl;
				};	
			for (int i=0;i<=maxidx;i++)
				{
					thl[i]=tht[i]=tll[i]=tlt[i]=0.0;
					for (int j=(int)(2*hingeprc); j<250;j++)
						thl[i]+=sqrt((sxe[i][j]-sxe[i][j+1])*(sxe[i][j]-sxe[i][j+1])+(sze[i][j]-sze[i][j+1])*(sze[i][j]-sze[i][j+1]));
					for (int j=0; j<(int)(2*hingeprc); j++)	
						tht[i]+=sqrt((sxe[i][j]-sxe[i][j+1])*(sxe[i][j]-sxe[i][j+1])+(sze[i][j]-sze[i][j+1])*(sze[i][j]-sze[i][j+1]));
					for (int j=250;j<(500-(int)(2*hingeprc));j++)
						tll[i]+=sqrt((sxe[i][j]-sxe[i][j+1])*(sxe[i][j]-sxe[i][j+1])+(sze[i][j]-sze[i][j+1])*(sze[i][j]-sze[i][j+1]));
					for (int j=500-(int)(2*hingeprc);j<500;j++)	
						tlt[i]+=sqrt((sxe[i][j]-sxe[i][j+1])*(sxe[i][j]-sxe[i][j+1])+(sze[i][j]-sze[i][j+1])*(sze[i][j]-sze[i][j+1]));
				};
			remove("template.dxf");		df.open("template.dxf");
			df << "0" << endl << "SECTION" << endl << "2" << endl << "HEADER" << endl ;
			df << "9" << endl << "$EXTMAX" << endl << "10" << endl << "	-900." << endl << "20" << endl << "	90." << endl;
			df << "9" << endl << "$EXTMIN" << endl << "10" << endl << "	5000." << endl << "20" << endl << "	-3000." << endl;
			df << "0" << endl << "ENDSEC" << endl << "0" << endl << "SECTION" << endl << "2" << endl << "ENTITIES" << endl << "0" << endl;
			dxf_line(df,-300,+cmpos,300,+cmpos,1);	dxf_line(df,-300,-cmpos,300,-cmpos,1);	dxf_line(df,-300,+mtpos,300,+mtpos,1);	dxf_line(df,-300,-mtpos,300,-mtpos,1);
			for (int i=0;i<maxidx;i++)
				{
					dxf_line(df,-thl[i],+sy[i],-thl[i+1],+sy[i+1],8);	
					dxf_line(df,+tht[i],+sy[i],+tht[i+1],+sy[i+1],8);
					dxf_line(df,-tll[i],-sy[i],-tll[i+1],-sy[i+1],8);
					dxf_line(df,+tlt[i],-sy[i],+tlt[i+1],-sy[i+1],8);
				};
			df << "0" << endl << "ENDSEC" << endl << "0" << endl << "EOF" << endl ;		df.close(); 
			//compute inner normal
			for (int i=0;i<=maxidx;i++)
				{
					sx[i][0]=sxe[i][0]-compfoilt;
					sx[i][250]=sxe[i][250]+compfoill;
					sx[i][500]=sx[i][0];
					sz[i][0]=sze[i][0];	
					sz[i][250]=sze[i][250];	
					sz[i][500]=sze[i][500];
					for (int j=1;j<250;j++)
						{	
							cf=0.0;
							if (j<=2*ctrlprc)		
								cf=compfoilt;
							else 
								if (j<=2*sparprc)		
									cf=compfoilm;
								else 
									if (j<=2*dboxprc)	
										cf=compfoils;
									else			
										cf=compfoill;
							cx1=sxe[i][j-1];
							cx2=sxe[i][j+1];	
							dcx=cx2-cx1;	
							cy1=sze[i][j-1];
							cy2=sze[i][j+1];
							dcy=cy2-cy1;
							dr=sqrt(dcx*dcx+dcy*dcy);
							sx[i][j]=sxe[i][j]-dcy/dr*cf;
							sz[i][j]=sze[i][j]+dcx/dr*cf;
						};
					for (int j=251;j<500;j++)
						{	
							cf=0.0;
							if ((500-j)>=2*dboxprc)			
								cf=compfoill;
							else 
								if ((500-j)>=2*sparprc)  
									cf=compfoils;
								else 
									if ((500-j)>=2*ctrlprc)
										cf=compfoilm;
									else 	
										cf=compfoilt;
							cx1=sxe[i][j-1];	
							cx2=sxe[i][j+1];	
							dcx=cx2-cx1;
							cy1=sze[i][j-1];
							cy2=sze[i][j+1];	
							dcy=cy2-cy1;
							dr=sqrt(dcx*dcx+dcy*dcy);
							sx[i][j]=(cx1+cx2)/2-dcy/dr*cf;	
							sz[i][j]=(cy1+cy2)/2+dcx/dr*cf;
						};
					//adjust inner normal
					for (int j=0;j<=250;j++) 
						{	
							hi=sz[i][j];
							lo=sz[i][500-j];
							hl=(hi+lo)/2;	
							sz[i][j]=(hi>lo?sz[i][j]:hl); 
							sz[i][500-j]=(hi>lo?sz[i][500-j]:hl);	
						};
				};

			cout << endl << "STAGE 5:computing servo pockets and control surfaces hinges" << endl ;

			for (int i1=asy_find(srv1d1);i1<=asy_find(srv1d2); i1++)
				for (int j1=(int)(srv1p1*2.);j1<=(int)(srv1p2*2.);j1++) 
					sz[i1][500-j1]=sz[i1][500-j1]+srv1z;
			for (int i2=asy_find(srv2d1);i2<=asy_find(srv2d2); i2++) 
				for (int j2=(int)(srv2p1*2.);j2<=(int)(srv2p2*2.);j2++) 
					sz[i2][500-j2]=sz[i2][500-j2]+srv2z;
			for (int i3=asy_find(srv3d1);i3<=asy_find(srv3d2); i3++)
				for (int j3=(int)(srv3p1*2.);j3<=(int)(srv3p2*2.);j3++) 
					sz[i3][500-j3]=sz[i3][500-j3]+srv3z;

			cmin=asy_find(ctrl1d1)+1; 
			cmax=asy_find(ctrl1d2)-1; 
			cpos=(elliptic==0?(int)(2.*ctrl1prc):(int)(2.*hingeprc));
			switch(ctrl1typ)
				{
					case 1: 
								for (int f=cmin;f<=cmax;f++) 
									{ 
										sz[f][cpos]= sz[f][500-cpos]+ctrl1dim; 
										if (Q==2)	
											sz[f][cpos+1]= sz[f][500-cpos-1]+ctrl1dim;
									};
								break;
					case 2: 
								for (int f=cmin;f<=cmax;f++)
									{ 
										sz[f][500-cpos]= sz[f][cpos]-ctrl1dim;
										if (Q==2)
											sz[f][500-cpos-1]= sz[f][cpos+1]-ctrl1dim;
									};
								break;
					case 3: 
								for (int f=cmin;f<=cmax;f++)
									{
										hinge=(sz[f][500-cpos]+sz[f][cpos])/2; 
										sz[f][cpos]=hinge+ctrl1dim/2;
										if (Q==2)
											sz[f][cpos+1]=hinge+ctrl1dim/2;
										sz[f][500-cpos]=hinge-ctrl1dim/2;
										if (Q==2)
											sz[f][500-cpos-1]=hinge-ctrl1dim/2;
									};
								break;
				};
			cmin=asy_find(ctrl2d1)+1;
			cmax=asy_find(ctrl2d2)-1;
			cpos=(elliptic==0?(int)(2.*ctrl2prc):(int)(2.*hingeprc));
			switch(ctrl2typ)
				{
					case 1: 
								for (int f=cmin;f<=cmax;f++)
									{ 
										sz[f][cpos]= sz[f][500-cpos]+ctrl2dim; 
										if (Q==2)
											sz[f][cpos+1]= sz[f][500-cpos-1]+ctrl2dim; 
									}; 
								break;
					case 2: 
								for (int f=cmin;f<=cmax;f++)
									{ 
										sz[f][500-cpos]= sz[f][cpos]-ctrl2dim;
										if (Q==2)
											sz[f][500-cpos-1]= sz[f][cpos+1]-ctrl2dim;
									};
								break;
					case 3: 
								for (int f=cmin;f<=cmax;f++) 
									{
										hinge=(sz[f][500-cpos]+ sz[f][cpos])/2; 
										sz[f][cpos]=hinge+ctrl2dim/2;
										if (Q==2)
											sz[f][cpos+1]=hinge+ctrl2dim/2;
										sz[f][500-cpos]=hinge-ctrl2dim/2;
										if (Q==2)
											sz[f][500-cpos-1]=hinge-ctrl2dim/2;
									};
								break;
				};
			cmin=asy_find(ctrl3d1)+1;
			cmax=asy_find(ctrl3d2)-1; 
			cpos=(elliptic==0?(int)(2.*ctrl3prc):(int)(2.*hingeprc));
			switch(ctrl3typ)
				{
					case 1: 
								for (int f=cmin;f<=cmax;f++) 
									{
										sz[f][cpos]= sz[f][500-cpos]+ctrl3dim; 
										if (Q==2) 
											sz[f][cpos+1]= sz[f][500-cpos-1]+ctrl3dim; 
									}; 
								break;
					case 2:	
								for (int f=cmin;f<=cmax;f++)	
									{
										sz[f][500-cpos]= sz[f][cpos]-ctrl3dim;
										if (Q==2)
											sz[f][500-cpos-1]= sz[f][cpos+1]-ctrl3dim; 
									};
								break;
					case 3: 
								for (int f=cmin;f<=cmax;f++) 
									{
										hinge=(sz[f][500-cpos]+ sz[f][cpos])/2;
										sz[f][cpos]=hinge+ctrl3dim/2; 	
										if (Q==2) 
											sz[f][cpos+1]=hinge+ctrl3dim/2;
										sz[f][500-cpos]=hinge-ctrl3dim/2;	
										if (Q==2) 
											sz[f][500-cpos-1]=hinge-ctrl3dim/2;
									};	
								break;
				};

			tx0[0]=sx[0][0]+tew;
			tx1[0]=sx[0][0]+ttd;
			tx2[0]=sx[0][0]+ttd+tts;
			tx3[0]=sx[0][0]+ted;
			
			ty0[0]=sy[0];
			ty1[0]=sy[0];	
			ty2[0]=sy[0];
			ty3[0]=sy[0];
			
			lx0[0]=sx[0][250]-lew;	
			lx1[0]=sx[0][250]-ltd;	
			lx2[0]=sx[0][250]-ltd-lts;	
			lx3[0]=sx[0][250]-led;
			
			ly0[0]=sy[0];
			ly1[0]=sy[0];
			ly2[0]=sy[0];
			ly3[0]=sy[0];
			
			for (int i=1;i<=mybrk;i++)
				{
					cx1=sx[i-1][0];
					cx2=sx[i][0];	
					cy1=sy[i-1];
					cy2=sy[i];
					dcx=sqrt((cx2-cx1)*(cx2-cx1));	
					dcy=sqrt((cy2-cy1)*(cy2-cy1));
					dr=sqrt(dcx*dcx+dcy*dcy);
					tx0[i]=cx2+dcy/dr*tew;	
					ty0[i]=cy2+dcx/dr*tew;
					tx1[i]=cx2+dcy/dr*ttd;	
					ty1[i]=cy2+dcx/dr*ttd;
					tx2[i]=cx2+dcy/dr*(ttd+tts);
					ty2[i]=cy2+dcx/dr*(ttd+tts);
					tx3[i]=((mr==0)?(cx2+dcy/dr*ted):tx3[0]);	
					ty3[i]=cy2+dcx/dr*ted;	
					cx1=sx[i-1][250];	
					cx2=sx[i][250];	
					cy1=sy[i-1];
					cy2=sy[i];
					dcx=sqrt((cx2-cx1)*(cx2-cx1));	
					dcy=sqrt((cy2-cy1)*(cy2-cy1));
					dr=sqrt(dcx*dcx+dcy*dcy);		
					lx0[i]=cx2-dcy/dr*lew;	
					ly0[i]=cy2+dcx/dr*lew;
					lx1[i]=cx2-dcy/dr*ltd;	
					ly1[i]=cy2+dcx/dr*ltd;
					lx2[i]=cx2-dcy/dr*(ltd+lts);
					ly2[i]=cy2+dcx/dr*(ltd+lts);
					lx3[i]=((mr==0)?(cx2-dcy/dr*led):lx3[0]);
					ly3[i]=cy2+dcx/dr*led;			
				};	

			cout << endl << "STAGE 6:computing mesh model for wing segments" << endl ;

			for (l=0;l<120;l++) 
				{ 
					sx[l][500]=sx[l][0];	
					sz[l][500]=sz[l][0]; 
				};
/*
int stl_surface(int cw, int my, int dw, ofstream &stl_file, int biy, int eiy)
int stl_contour(int cw, int my, int db, int dw, ofstream &stl_file, int biy, int eiy, int exb, int exe)
int stl_slicecut(int cw, int my, ofstream &stl_file, int piy)
int stl_tipclose(int cw, int my, int db, ofstream &stl_file)
int stl_vertical(int cw, int my, ofstream &stl_file, int piy, int src, int dst,int dirz, int yadj)
int stl_terminal(int cw, int my, ofstream &stl_file, int piy, int dir, int dirz)
*/			
			if (to!=0)
				{
					remove("test.stl");	
					stl.open("test.stl");	
					stl << "solid test " << endl;
					cout << "test ";
					//mold_lo_tip_right
					stl_surface(+1,+1,-1,stl,mtbrk,mybrk);
					stl_contour(-1,+1,-1,-1,stl,mtbrk,mybrk,1,0);
					stl_vertical(+1,+1,stl,mtbrk,0,-1,-1,0);
					stl_vertical(-1,+1,stl,mtbrk,0,-2,-1,-1);
					stl_terminal(-1,+1,stl,mtbrk,-1,-1);		
					stl_tipclose(-1,+1,-1,stl);
					stl << "endsolid test" << endl;
					stl.close(); 
					cout << " ... done" << endl;
				}
			if (wo!=0)
				{
					//ok
					cout << "wing_left ";
					remove("wing_left.stl");
					stl.open("wing_left.stl");
					stl << "solid wing_left" << endl;	
					stl_surface(+1,-1,0,stl,0,mybrk);
					stl_slicecut(-1,-1,stl,0);
					stl	<< "endsolid wing_left" << endl;
					stl.close();
					cout << " ... done" << endl ;
					
					//ok
					cout   << "wing_right ";
					remove("wing_right.stl");
					stl.open("wing_right.stl");	
					stl << "solid wing_right" << endl;
					stl_surface(+1,+1,0,stl,0,mybrk);
					stl_slicecut(-1,-1,stl,0);
					stl << "endsolid wing_right" << endl;
					stl.close();
					cout << " ... done" << endl ;
					
					//ok
					cout << "wing_full ";
					remove("wing_full.stl");
					stl.open("wing_full.stl");	
					stl << "solid wing_full" << endl;		
					stl_surface(-1,+1,0,stl,0,mybrk);		
					stl_surface(+1,-1,0,stl,0,mybrk);
					stl << "endsolid wing_full" << endl;
					stl.close();
					cout << " ... done" << endl ;
					
					if (sf==0)
						{
							//ok
							cout << "wing_central " ;
							remove("wing_central.stl");	
							stl.open("wing_central.stl");
							stl << "solid wing_central" << endl;	
							stl_slicecut(-1,-1,stl,cmbrk);
							stl_surface(-1,+1,0,stl,0,cmbrk);		
							stl_surface(+1,-1,0,stl,0,cmbrk);
							stl_slicecut(+1,+1,stl,cmbrk);
							stl  << "endsolid wing_central" << endl;
							stl.close();
							cout << " ... done" << endl ;
							
							//ok
							cout << "wing_cnt_left ";	
							remove("wing_cnt_left.stl");
							stl.open("wing_cnt_left.stl");	
							stl << "solid wing_cnt_left" << endl;		 
							stl_slicecut(+1,-1,stl,0);
							stl_surface(+1,-1,0,stl,0,cmbrk);
							stl_slicecut(-1,-1,stl,cmbrk);
							stl << "endsolid wing_cnt_left" << endl;
							stl.close();
							cout << " ... done" << endl ;
							
							//ok
							cout << "wing_cnt_right ";
							remove("wing_cnt_right.stl");
							stl.open("wing_cnt_right.stl");
							stl  << "solid wing_cnt_right" << endl;	 
							stl_slicecut(-1,+1,stl,0);
							stl_surface(-1,+1,0,stl,0,cmbrk);
							stl_slicecut(+1,+1,stl,cmbrk);
							stl << "endsolid wing_cnt_right" << endl;
							stl.close();
							cout << " ... done" << endl ;
							
							//ok
							cout << "wing_mid_left ";
							remove("wing_mid_left.stl");
							stl.open("wing_mid_left.stl");
							stl << "solid wing_mid_left" << endl;		 
							stl_slicecut(+1,-1,stl,cmbrk);
							stl_surface(+1,-1,0,stl,cmbrk,mtbrk);
							stl_slicecut(-1,-1,stl,mtbrk);
							stl << "endsolid wing_mid_left" << endl;
							stl.close();
							cout << " ... done" << endl ;
							
							//ok
							cout << "wing_mid_right ";
							remove("wing_mid_right.stl");
							stl.open("wing_mid_right.stl");
							stl  << "solid wing_mid_right" << endl;	 
							stl_slicecut(-1,+1,stl,cmbrk);
							stl_surface(-1,+1,0,stl,cmbrk,mtbrk);
							stl_slicecut(+1,+1,stl,mtbrk);
							stl << "endsolid wing_mid_right" << endl;
							stl.close();
							cout << " ... done" << endl ;
							
							//ok
							cout << "wing_tip_left ";
							remove("wing_tip_left.stl");
							stl.open("wing_tip_left.stl");
							stl << "solid wing_tip_left" << endl;		
							stl_slicecut(+1,-1,stl,mtbrk);
							stl_surface(+1,-1,0,stl,mtbrk,mybrk);
							stl_slicecut(-1,-1,stl,mybrk);
							stl << "endsolid wing_tip_left" << endl;
							stl.close();
							cout << " ... done" << endl ;
							
							//ok
							cout << "wing_tip_right ";	
							remove("wing_tip_right.stl");
							stl.open("wing_tip_right.stl");	
							stl << "solid wing_tip_right" << endl;		
							stl_slicecut(-1,+1,stl,mtbrk);
							stl_surface(-1,+1,0,stl,mtbrk,mybrk);
							stl_slicecut(+1,+1,stl,mybrk);
							stl << "endsolid wing_tip_right" << endl;
							stl.close();
							cout << " ... done" << endl ;
						};	
				};
			if (co!=0)
				{
					if (sf!=0)
						{
							//ok
							cout << "core_hi_left ";
							remove("core_hi_left.stl");	
							stl.open("core_hi_left.stl");
							stl << "solid core_hi_left" << endl;		
							stl_surface(+1,-1,+1,stl,0,mybrk);
							stl << "endsolid core_hi_left" << endl; 
							stl.close();
							cout << " ... done" << endl ;
							
							//ok
							cout << "core_lo_left ";
							remove("core_lo_left.stl");	
							stl.open("core_lo_left.stl");
							stl << "solid core_lo_left" << endl;			
							stl_surface(+1,-1,-1,stl,0,mybrk);
							stl << "endsolid core_lo_left" << endl;
							stl.close();
							cout << " ... done" << endl ;
							
							//ok
							cout << "core_hi_right ";	
							remove("core_hi_right.stl");
							stl.open("core_hi_right.stl");
							stl << "solid core_hi_right" << endl;		
							stl_surface(+1,+1,+1,stl,0,mybrk);
							stl << "endsolid core_hi_right" << endl; 
							stl.close();
							cout << " ... done" << endl ;
							
							//ok
							cout << "core_lo_right ";
							remove("core_lo_right.stl");
							stl.open("core_lo_right.stl");
							stl << "solid core_lo_right" << endl;		
							stl_surface(+1,+1,-1,stl,0,mybrk);
							stl   << "endsolid core_lo_right" << endl;
							stl.close();
							cout << " ... done" << endl ;
							
							//ok
							cout << "core_hi_full ";
							remove("core_hi_full.stl");
							stl.open("core_hi_full.stl");
							stl << "solid core_hi_full" << endl;			
							stl_surface(-1,+1,+1,stl,0,mybrk);		
							stl << "endsolid core_hi_full" << endl; 
							stl.close();
							cout << " ... done" << endl ;
							
							//ok
							cout << "core_lo_full ";
							remove("core_lo_full.stl");	
							stl.open("core_lo_full.stl");
							stl << "solid core_lo_full" << endl;			
							stl_surface(-1,+1,-1,stl,0,mybrk);		
							stl_surface(+1,-1,-1,stl,0,mybrk);
							stl << "endsolid core_lo_full" << endl; 
							stl.close();
							cout << " ... done" << endl ;
						}	
					else
						{
							//ok
							cout << "core_hi_central ";
							remove("core_hi_central.stl");
							stl.open("core_hi_central.stl");
							stl << "solid core_hi_central" << endl;						
							stl_surface(-1,+1,+1,stl,0,cmbrk);	
							stl_surface(+1,-1,+1,stl,0,cmbrk);
							stl  << "endsolid core_hi_central" << endl; 
							stl.close();
							cout << " ... done" << endl ;			
							
							//ok
							cout << "core_hi_cnt_left "; 
							remove("core_hi_cnt_left.stl");
							stl.open("core_hi_cnt_left.stl");
							stl << "solid core_hi_cnt_left" << endl;		
							stl_surface(+1,-1,+1,stl,0,cmbrk);	
							stl << "endsolid core_hi_cnt_left" << endl;	
							stl.close();
							cout << " ... done" << endl ;
							
							//ok
							cout << "core_hi_cnt_right ";
							remove("core_hi_cnt_right.stl");
							stl.open("core_hi_cnt_right.stl");	
							stl << "solid core_hi_cnt_right" << endl;		
							stl_surface(-1,+1,+1,stl,0,cmbrk);	
							stl << "endsolid core_hi_cnt_right" << endl;
							stl.close();
							cout << " ... done" << endl ;			
							
							//ok
							cout << "core_hi_mid_left ";
							remove("core_hi_mid_left.stl");
							stl.open("core_hi_mid_left.stl");
							stl << "solid core_hi_mid_left" << endl;		
							stl_surface(+1,-1,+1,stl,cmbrk,mtbrk);	
							stl << "endsolid core_hi_mid_left" << endl;
							stl.close();
							cout << " ... done" << endl ;
							
							//ok
							cout << "core_hi_mid_right ";
							remove("core_hi_mid_right.stl");
							stl.open("core_hi_mid_right.stl");
							stl << "solid core_hi_mid_right" << endl;		
							stl_surface(-1,+1,+1,stl,cmbrk,mtbrk);	
							stl << "endsolid core_hi_mid_right" << endl;
							stl.close();
							cout << " ... done" << endl ;			
							
							//ok
							cout << "core_hi_tip_left ";
							remove("core_hi_tip_left.stl");	
							stl.open("core_hi_tip_left.stl");
							stl << "solid core_hi_tip_left" << endl;		
							stl_surface(+1,-1,+1,stl,mtbrk,mybrk);					
							stl << "endsolid core_hi_tip_left" << endl;
							stl.close();
							cout << " ... done" << endl ;
							
							//ok
							cout << "core_hi_tip_right ";
							remove("core_hi_tip_right.stl");
							stl.open("core_hi_tip_right.stl");
							stl << "solid core_hi_tip_right" << endl;			
							stl_surface(-1,+1,+1,stl,mtbrk,mybrk);
							stl << "endsolid core_hi_tip_right" << endl;
							stl.close();
							cout << " ... done" << endl ;
							
							//ok
							cout << "core_lo_central ";	
							remove("core_lo_central.stl");
							stl.open("core_lo_central.stl"); 
							stl  << "solid core_lo_central" << endl;		 
							stl_surface(+1,+1,-1,stl,0,cmbrk);
							stl_surface(-1,-1,-1,stl,0,cmbrk);
							stl  << "endsolid core_lo_central" << endl;
							stl.close();
							cout << " ... done" << endl ;
							
							//ok
							cout << "core_lo_cnt_left ";
							remove("core_lo_cnt_left.stl");	
							stl.open("core_lo_cnt_left.stl");
							stl << "solid core_lo_cnt_left" << endl;		
							stl_surface(-1,-1,-1,stl,0,cmbrk);	
							stl << "endsolid core_lo_cnt_left" << endl;
							stl.close();
							cout << " ... done" << endl ;	
							
							//ok
							cout << "core_lo_cnt_right ";
							remove("core_lo_cnt_right.stl");
							stl.open("core_lo_cnt_right.stl");
							stl << "solid core_lo_cnt_right" << endl;		
							stl_surface(+1,+1,-1,stl,0,cmbrk);	
							stl << "endsolid core_lo_cnt_right" << endl;
							stl.close();
							cout << " ... done" << endl ;					
							
							//ok
							cout << "core_lo_mid_left ";
							remove("core_lo_mid_left.stl");
							stl.open("core_lo_mid_left.stl");
							stl << "solid core_lo_mid_left" << endl;		
							stl_surface(-1,-1,-1,stl,cmbrk,mtbrk);	
							stl << "endsolid core_lo_mid_left" << endl;
							stl.close();
							cout << " ... done" << endl ;	
							
							//ok
							cout << "core_lo_mid_right ";
							remove("core_lo_mid_right.stl");
							stl.open("core_lo_mid_right.stl");
							stl << "solid core_lo_mid_right" << endl;		
							stl_surface(+1,+1,-1,stl,cmbrk,mtbrk);	
							stl << "endsolid core_lo_mid_right" << endl;
							stl.close();
							cout << " ... done" << endl ;					
							
							//ok
							cout << "core_lo_tip_left ";
							remove("core_lo_tip_left.stl"); 
							stl.open("core_lo_tip_left.stl");
							stl << "solid core_lo_tip_left" << endl;		
							stl_surface(-1,-1,-1,stl,mtbrk,mybrk);	
							stl << "endsolid core_lo_tip_left" << endl;
							stl.close();
							cout << " ... done" << endl ;
							
							//ok
							cout << "core_lo_tip_right " ;
							remove("core_lo_tip_right.stl");
							stl.open("core_lo_tip_right.stl");
							stl << "solid core_lo_tip_right" << endl;		
							stl_surface(+1,+1,-1,stl,mtbrk,mybrk);	
							stl << "endsolid core_lo_tip_right" << endl;
							stl.close();
							cout << " ... done" << endl;	
						};
				};	
			if (mo!=0)
				{
							//ok
							cout << "mold_hi_full ";
							remove("mold_hi_full.stl");
							stl.open("mold_hi_full.stl");
							stl << "solid mold_hi_full" << endl;		
							stl_surface(+1,+1,+1,stl,0,mybrk);
							stl_surface(-1,-1,+1,stl,0,mybrk);
							stl_contour(+1,+1,+1,+1,stl,0,mybrk,0,0);
							stl_contour(-1,-1,+1,+1,stl,0,mybrk,0,0);
							stl_tipclose(+1,+1,+1,stl);
							stl_tipclose(-1,-1,+1,stl);
							stl  << "endsolid mold_hi_full" << endl;
							stl.close();
							cout << " ... done" << endl;	
							 
							//ok
							cout << "mold_hi_left ";
							remove("mold_hi_left.stl");	
							stl.open("mold_hi_left.stl");
							stl << "solid mold_hi_left" << endl;		
							stl_surface(-1,-1,+1,stl,0,mybrk);
							stl_contour(-1,-1,+1,+1,stl,0,mybrk,1,0);
							stl_vertical(+1,-1,stl,0,0,+1,+1,0);
							stl_vertical(-1,-1,stl,0,0,+2,+1,-1);
							stl_terminal(-1,-1,stl,0,-1,+1);		
							stl_tipclose(-1,-1,+1,stl);
							stl  << "endsolid mold_hi_left" << endl;
							stl.close();
							cout << " ... done" << endl;	
							 
							//ok
							cout << "mold_hi_right ";
							remove("mold_hi_right.stl");
							stl.open("mold_hi_right.stl");
							stl << "solid mold_hi_right" << endl;	
							stl_surface(+1,+1,+1,stl,0,mybrk);
							stl_contour(+1,+1,+1,+1,stl,0,mybrk,1,0);
							stl_vertical(+1,+1,stl,0,0,+1,+1,0);
							stl_vertical(+1,+1,stl,0,0,+2,+1,-1);
							stl_terminal(+1,+1,stl,0,-1,+1);		
							stl_tipclose(+1,+1,+1,stl);
							stl  << "endsolid mold_hi_right" << endl;
							stl.close();
							cout << " ... done" << endl;	
							 
							//ok
							cout << "mold_lo_full ";
							remove("mold_lo_full.stl");	
							stl.open("mold_lo_full.stl");
							stl << "solid mold_lo_full" << endl;		 
							stl_surface(+1,+1,-1,stl,0,mybrk);
							stl_surface(-1,-1,-1,stl,0,mybrk);
							stl_contour(-1,+1,-1,-1,stl,0,mybrk,0,0);
							stl_contour(+1,-1,-1,-1,stl,0,mybrk,0,0);
							stl_tipclose(-1,+1,-1,stl);
							stl_tipclose(+1,-1,-1,stl);
							stl  << "endsolid mold_lo_full" << endl;
							stl.close();
							cout << " ... done" << endl;	
							 
							//ok
							cout << "mold_lo_left ";
							remove("mold_lo_left.stl");	
							stl.open("mold_lo_left.stl");
							stl << "solid mold_lo_left" << endl;		 
							stl_surface(-1,-1,-1,stl,0,mybrk);
							stl_contour(+1,-1,-1,-1,stl,0,mybrk,1,0);
							stl_vertical(-1,-1,stl,0,0,-1,-1,0);
							stl_vertical(+1,-1,stl,0,0,-2,-1,-1);
							stl_terminal(+1,-1,stl,0,-1,-1);		
							stl_tipclose(+1,-1,-1,stl);
							stl  << "endsolid mold_lo_left" << endl;
							stl.close();
							cout << " ... done" << endl;	
							 
							//ok
							cout << "mold_lo_right ";
							remove("mold_lo_right.stl");
							stl.open("mold_lo_right.stl");
							stl << "solid mold_lo_right" << endl;	 
							stl_surface(+1,+1,-1,stl,0,mybrk);
							stl_contour(-1,+1,-1,-1,stl,0,mybrk,1,0);
							stl_vertical(+1,+1,stl,0,0,-1,-1,0);
							stl_vertical(-1,+1,stl,0,0,-2,-1,-1);
							stl_terminal(-1,+1,stl,0,-1,-1);		
							stl_tipclose(-1,+1,-1,stl);
							stl  << "endsolid mold_lo_right" << endl;
							stl.close();
							cout << " ... done" << endl;	
							 
							//ok
							cout << "mold_hi_central ";
							remove("mold_hi_central.stl");
							stl.open("mold_hi_central.stl");
							stl  << "solid mold_hi_central" << endl;
							stl_surface(+1,+1,+1,stl,0,cmbrk);
							stl_surface(-1,-1,+1,stl,0,cmbrk);
							stl_contour(+1,+1,+1,+1,stl,0,cmbrk,0,1);
							stl_contour(-1,-1,+1,+1,stl,0,cmbrk,1,1);
							stl_vertical(-1,-1,stl,cmbrk,0,+1,+1,0);
							stl_vertical(+1,+1,stl,cmbrk,0,+1,+1,0);
							stl_vertical(+1,-1,stl,cmbrk,0,+2,+1,+1);
							stl_vertical(-1,+1,stl,cmbrk,0,+2,+1,+1);
							stl_terminal(+1,-1,stl,cmbrk,+1,+1);
							stl_terminal(-1,+1,stl,cmbrk,+1,+1);
							stl << "endsolid mold_hi_central" << endl;
							stl.close();
							cout << " ... done" << endl;	
							 
							//ok
							cout << "mold_hi_cnt_left ";
							remove("mold_hi_cnt_left.stl");
							stl.open("mold_hi_cnt_left.stl");
							stl << "solid mold_hi_cnt_left" << endl;
							stl_surface(-1,-1,+1,stl,0,cmbrk);
							stl_contour(-1,-1,+1,+1,stl,0,cmbrk,1,1);
							stl_vertical(+1,-1,stl,0,0,+1,+1,0);
							stl_vertical(-1,-1,stl,0,0,+2,+1,-1);
							stl_vertical(-1,-1,stl,cmbrk,0,+1,+1,0);
							stl_vertical(+1,-1,stl,cmbrk,0,+2,+1,+1);
							stl_terminal(-1,-1,stl,0,-1,+1);
							stl_terminal(+1,-1,stl,cmbrk,+1,+1);
							stl << "endsolid mold_hi_cnt_left" << endl;	
							stl.close();
							cout << " ... done" << endl;				
							 
							//ok
							cout << "mold_hi_cnt_right ";
							remove("mold_hi_cnt_right.stl");
							stl.open("mold_hi_cnt_right.stl");
							stl << "solid mold_hi_cnt_right" << endl; 
							stl_surface(+1,+1,+1,stl,0,cmbrk);
							stl_contour(+1,+1,+1,+1,stl,0,cmbrk,1,1);
							stl_vertical(-1,+1,stl,0,0,+1,+1,0);
							stl_vertical(+1,+1,stl,0,0,+2,+1,-1);
							stl_vertical(+1,+1,stl,cmbrk,0,+1,+1,0);
							stl_vertical(-1,+1,stl,cmbrk,0,+2,+1,+1);
							stl_terminal(+1,+1,stl,0,-1,+1);
							stl_terminal(-1,+1,stl,cmbrk,+1,+1);
							stl << "endsolid mold_hi_cnt_right" << endl;
							stl.close();
							cout << " ... done" << endl;	
							 
							//ok
							cout << "mold_hi_mid_left ";
							remove("mold_hi_mid_left.stl");	
							stl.open("mold_hi_mid_left.stl");
							stl << "solid mold_hi_mid_left" << endl;
							stl_surface(-1,-1,+1,stl,cmbrk,mtbrk);
							stl_contour(-1,-1,+1,+1,stl,cmbrk,mtbrk,1,1);
							stl_vertical(+1,-1,stl,cmbrk,0,+1,+1,0);
							stl_vertical(-1,-1,stl,cmbrk,0,+2,+1,-1);
							stl_vertical(+1,-1,stl,mtbrk,0,+1,+1,0);
							stl_vertical(+1,-1,stl,mtbrk,0,+2,+1,+1);
							stl_terminal(-1,-1,stl,cmbrk,-1,+1);
							stl_terminal(+1,-1,stl,mtbrk,+1,+1);
							stl << "endsolid mold_hi_mid_left" << endl;
							stl.close();
							cout << " ... done" << endl;				
							 
							//ok
							cout << "mold_hi_mid_right ";
							remove("mold_hi_mid_right.stl");
							stl.open("mold_hi_mid_right.stl");
							stl << "solid mold_hi_mid_right" << endl; 
							stl_surface(+1,+1,+1,stl,cmbrk,mtbrk);
							stl_contour(+1,+1,+1,+1,stl,cmbrk,mtbrk,1,1);
							stl_vertical(-1,+1,stl,cmbrk,0,+1,+1,0);
							stl_vertical(+1,+1,stl,cmbrk,0,+2,+1,-1);
							stl_vertical(+1,+1,stl,mtbrk,0,+1,+1,0);
							stl_vertical(-1,+1,stl,mtbrk,0,+2,+1,+1);
							stl_terminal(+1,+1,stl,cmbrk,-1,+1);
							stl_terminal(-1,+1,stl,mtbrk,+1,+1);
							stl << "endsolid mold_hi_mid_right" << endl;
							stl.close();
							cout << " ... done" << endl;	
							 
							//ok
							cout << "mold_hi_tip_left ";
							remove("mold_hi_tip_left.stl");	
							stl.open("mold_hi_tip_left.stl");
							stl << "solid mold_hi_tip_left" << endl;
							stl_surface(-1,-1,+1,stl,mtbrk,mybrk);
							stl_contour(-1,-1,+1,+1,stl,mtbrk,mybrk,1,0);
							stl_tipclose(-1,-1,+1,stl);
							stl_vertical(-1,-1,stl,mtbrk,0,+1,+1,0);
							stl_vertical(-1,-1,stl,mtbrk,0,+2,+1,-1);
							stl_terminal(-1,-1,stl,mtbrk,-1,+1);
							stl << "endsolid mold_hi_tip_left" << endl;
							stl.close();
							cout << " ... done" << endl;		
							 
							//ok
							cout  << "mold_hi_tip_right ";
							remove("mold_hi_tip_right.stl");
							stl.open("mold_hi_tip_right.stl");
							stl << "solid mold_hi_tip_right" << endl;		
							stl_surface(+1,+1,+1,stl,mtbrk,mybrk);
							stl_contour(+1,+1,+1,+1,stl,mtbrk,mybrk,1,0);
							stl_tipclose(+1,+1,+1,stl);
							stl_vertical(-1,+1,stl,mtbrk,0,+1,+1,0);
							stl_vertical(+1,+1,stl,mtbrk,0,+2,+1,-1);
							stl_terminal(+1,+1,stl,mtbrk,-1,+1);
							stl << "endsolid mold_hi_tip_right" << endl;
							stl.close();
							cout << " ... done" << endl;	
							 
							//ok
							cout << "mold_lo_central ";	
							remove("mold_lo_central.stl");
							stl.open("mold_lo_central.stl");
							stl << "solid mold_lo_central" << endl;		 
							stl_surface(+1,+1,-1,stl,0,cmbrk);
							stl_surface(-1,-1,-1,stl,0,cmbrk);
							stl_contour(-1,+1,-1,-1,stl,0,cmbrk,0,1);
							stl_contour(+1,-1,-1,-1,stl,0,cmbrk,1,1);
							stl_vertical(+1,-1,stl,cmbrk,0,-1,-1,0);
							stl_vertical(-1,+1,stl,cmbrk,0,-1,-1,0);
							stl_vertical(-1,-1,stl,cmbrk,0,-2,-1,+1);
							stl_vertical(+1,+1,stl,cmbrk,0,-2,-1,+1);
							stl_terminal(-1,-1,stl,cmbrk,+1,-1);
							stl_terminal(+1,+1,stl,cmbrk,+1,-1);		
							stl  << "endsolid mold_lo_central" << endl;
							stl.close();
							cout << " ... done" << endl;	
							 				
							//ok
							cout << "mold_lo_cnt_left ";
							remove("mold_lo_cnt_left.stl");
							stl.open("mold_lo_cnt_left.stl");
							stl << "solid mold_lo_cnt_left" << endl;		
							stl_surface(-1,-1,-1,stl,0,cmbrk);
							stl_contour(+1,-1,-1,-1,stl,0,cmbrk,1,1);
							stl_vertical(-1,-1,stl,0,0,-1,-1,0);
							stl_vertical(+1,-1,stl,0,0,-2,-1,-1);
							stl_vertical(+1,-1,stl,cmbrk,0,-1,-1,0);
							stl_vertical(-1,-1,stl,cmbrk,0,-2,-1,+1);
							stl_terminal(+1,-1,stl,0,-1,-1);
							stl_terminal(-1,-1,stl,cmbrk,+1,-1);		
							stl << "endsolid mold_lo_cnt_left" << endl;
							stl.close();
							cout << " ... done" << endl;	
							 					
							//ok
							cout << "mold_lo_cnt_right ";
							remove("mold_lo_cnt_right.stl"); 
							stl.open("mold_lo_cnt_right.stl");
							stl << "solid mold_lo_cnt_right" << endl;		
							stl_surface(+1,+1,-1,stl,0,cmbrk);
							stl_contour(-1,+1,-1,-1,stl,0,cmbrk,1,1);
							stl_vertical(+1,+1,stl,0,0,-1,-1,0);
							stl_vertical(-1,+1,stl,0,0,-2,-1,-1);
							stl_vertical(-1,+1,stl,cmbrk,0,-1,-1,0);
							stl_vertical(+1,+1,stl,cmbrk,0,-2,-1,+1);
							stl_terminal(-1,+1,stl,0,-1,-1);
							stl_terminal(+1,+1,stl,cmbrk,+1,-1);		
							stl << "endsolid mold_lo_cnt_right" << endl;
							stl.close();
							cout << " ... done" << endl;	
							 					
							//ok
							cout << "mold_lo_mid_left ";
							remove("mold_lo_mid_left.stl");
							stl.open("mold_lo_mid_left.stl");
							stl << "solid mold_lo_mid_left" << endl;		
							stl_surface(-1,-1,-1,stl,cmbrk,mtbrk);
							stl_contour(+1,-1,-1,-1,stl,cmbrk,mtbrk,1,1);
							stl_vertical(-1,-1,stl,cmbrk,0,-1,-1,0);
							stl_vertical(+1,-1,stl,cmbrk,0,-2,-1,-1);
							stl_vertical(+1,-1,stl,mtbrk,0,-1,-1,0);
							stl_vertical(-1,-1,stl,mtbrk,0,-2,-1,+1);
							stl_terminal(+1,-1,stl,cmbrk,-1,-1);
							stl_terminal(-1,-1,stl,mtbrk,+1,-1);		
							stl << "endsolid mold_lo_mid_left" << endl;
							stl.close();
							cout << " ... done" << endl;	
							 
							//ok
							cout << "mold_lo_mid_right ";
							remove("mold_lo_mid_right.stl");
							stl.open("mold_lo_mid_right.stl");
							stl << "solid mold_lo_mid_right" << endl;		
							stl_surface(+1,+1,-1,stl,cmbrk,mtbrk);
							stl_contour(-1,+1,-1,-1,stl,cmbrk,mtbrk,1,1);
							stl_vertical(+1,+1,stl,cmbrk,0,-1,-1,0);
							stl_vertical(-1,+1,stl,cmbrk,0,-2,-1,-1);
							stl_vertical(-1,+1,stl,mtbrk,0,-1,-1,0);
							stl_vertical(+1,+1,stl,mtbrk,0,-2,-1,+1);
							stl_terminal(-1,+1,stl,cmbrk,-1,-1);
							stl_terminal(+1,+1,stl,mtbrk,+1,-1);	
							stl << "endsolid mold_lo_mid_right" << endl;
							stl.close();
							cout << " ... done" << endl;	
							
							//ok
							cout << "mold_lo_tip_left ";
							remove("mold_lo_tip_left.stl"); 
							stl.open("mold_lo_tip_left.stl");
							stl << "solid mold_lo_tip_left" << endl;		
							stl_surface(-1,-1,-1,stl,mtbrk,mybrk);
							stl_contour(+1,-1,-1,-1,stl,mtbrk,mybrk,1,0);
							stl_vertical(-1,-1,stl,mtbrk,0,-1,-1,0);
							stl_vertical(+1,-1,stl,mtbrk,0,-2,-1,-1);
							stl_terminal(+1,-1,stl,mtbrk,-1,-1);		
							stl_tipclose(+1,-1,-1,stl);
							stl << "endsolid mold_lo_tip_left" << endl;
							stl.close();
							cout << " ... done" << endl;	
							
							//ok
							cout << "mold_lo_tip_right ";
							remove("mold_lo_tip_right.stl");
							stl.open("mold_lo_tip_right.stl");
							stl << "solid mold_lo_tip_right" << endl;					
							stl_surface(+1,+1,-1,stl,mtbrk,mybrk);
							stl_contour(-1,+1,-1,-1,stl,mtbrk,mybrk,1,0);
							stl_vertical(+1,+1,stl,mtbrk,0,-1,-1,0);
							stl_vertical(-1,+1,stl,mtbrk,0,-2,-1,-1);
							stl_terminal(-1,+1,stl,mtbrk,-1,-1);		
							stl_tipclose(-1,+1,-1,stl);
							stl << "endsolid mold_lo_tip_right" << endl;
							stl.close(); 
							cout << " ... done" << endl;	
				};	
			if (po!=0)
				{
					        //ok
							cout << "plug_hi_full ";
							remove("plug_hi_full.stl");
							stl.open("plug_hi_full.stl");
							stl << "solid plug_hi_full" << endl;		
							stl_surface(-1,+1,-1,stl,0,mybrk);
							stl_surface(+1,-1,-1,stl,0,mybrk);
							stl_contour(+1,+1,+1,-1,stl,0,mybrk,0,0);
							stl_contour(-1,-1,+1,-1,stl,0,mybrk,0,0);
							stl_tipclose(+1,+1,+1,stl);
							stl_tipclose(-1,-1,+1,stl);
							stl  << "endsolid plug_hi_full" << endl;
							stl.close();
							cout << " ... done" << endl;	
							
							//ok
							cout << "plug_hi_left ";
							remove("plug_hi_left.stl");
							stl.open("plug_hi_left.stl");
							stl << "solid plug_hi_left" << endl;
							stl_surface(+1,-1,-1,stl,0,mybrk);
							stl_contour(-1,-1,+1,-1,stl,0,mybrk,1,0);
							stl_vertical(+1,-1,stl,0,0,-1,-1,0);
							stl_vertical(-1,-1,stl,0,0,+2,-1,-1);
							stl_terminal(-1,-1,stl,0,-1,-1);		
							stl_tipclose(-1,-1,+1,stl);
							stl  << "endsolid plug_hi_left" << endl;
							stl.close();
							cout << " ... done" << endl;	
							
							//ok
							cout << "plug_hi_right ";
							remove("plug_hi_right.stl");
							stl.open("plug_hi_right.stl");	
							stl << "solid plug_hi_right" << endl;
							stl_surface(-1,+1,-1,stl,0,mybrk);
							stl_contour(+1,+1,+1,-1,stl,0,mybrk,1,0);
							stl_vertical(-1,+1,stl,0,0,-1,-1,0);
							stl_vertical(+1,+1,stl,0,0,+2,-1,-1);
							stl_terminal(+1,+1,stl,0,-1,-1);		
							stl_tipclose(+1,+1,+1,stl);
							stl  << "endsolid plug_hi_right" << endl;
							stl.close();
							cout << " ... done" << endl;	
							
							//ok
							cout << "plug_lo_full ";
							remove("plug_lo_full.stl");
							stl.open("plug_lo_full.stl");
							stl << "solid plug_lo_full" << endl;		 
							stl_surface(-1,+1,+1,stl,0,mybrk);
							stl_surface(+1,-1,+1,stl,0,mybrk);
							stl_contour(-1,+1,-1,-1,stl,0,mybrk,0,0);
							stl_contour(+1,-1,-1,-1,stl,0,mybrk,0,0);
							stl_tipclose(-1,+1,-1,stl);
							stl_tipclose(+1,-1,-1,stl);
							stl  << "endsolid plug_lo_full" << endl;
							stl.close();
							cout << " ... done" << endl;	
							
							//ok
							cout << "plug_lo_left ";
							remove("plug_lo_left.stl");	
							stl.open("plug_lo_left.stl");
							stl << "solid plug_lo_left" << endl;		 
							stl_surface(+1,-1,+1,stl,0,mybrk);
							stl_contour(+1,-1,-1,-1,stl,0,mybrk,1,0);
							stl_vertical(-1,-1,stl,0,0,+1,+1,0);
							stl_vertical(+1,-1,stl,0,0,-2,+1,-1);
							stl_terminal(+1,-1,stl,0,-1,+1);		
							stl_tipclose(+1,-1,-1,stl);
							stl  << "endsolid plug_lo_left" << endl;
							stl.close();
							cout << " ... done" << endl;	
							
							//ok
							cout << "plug_lo_right ";
							remove("plug_lo_right.stl");
							stl.open("plug_lo_right.stl");
							stl << "solid plug_lo_right" << endl;	 
							stl_surface(-1,+1,+1,stl,0,mybrk);
							stl_contour(-1,+1,-1,-1,stl,0,mybrk,1,0);
							stl_vertical(+1,+1,stl,0,0,+1,+1,0);
							stl_vertical(-1,+1,stl,0,0,-2,+1,-1);
							stl_terminal(-1,+1,stl,0,-1,+1);		
							stl_tipclose(-1,+1,-1,stl);
							stl  << "endsolid plug_lo_right" << endl;
							stl.close();
							cout << " ... done" << endl;	
							
							//ok
							cout << "plug_hi_central ";	
							remove("plug_hi_central.stl");	
							stl.open("plug_hi_central.stl");
							stl  << "solid plug_hi_central" << endl;
							stl_surface(-1,+1,-1,stl,0,cmbrk);
							stl_surface(+1,-1,-1,stl,0,cmbrk);
							stl_contour(+1,+1,+1,-1,stl,0,cmbrk,0,1);
							stl_contour(-1,-1,+1,-1,stl,0,cmbrk,1,1);
							stl_vertical(+1,+1,stl,cmbrk,0,-1,-1,0);
							stl_vertical(-1,-1,stl,cmbrk,0,-1,-1,0);
							stl_vertical(-1,+1,stl,cmbrk,0,+2,-1,+1);
							stl_vertical(+1,-1,stl,cmbrk,0,+2,-1,+1);
							stl_terminal(-1,+1,stl,cmbrk,+1,-1);
							stl_terminal(+1,-1,stl,cmbrk,+1,-1);
							stl  << "endsolid plug_hi_central" << endl;
							stl.close();		
							cout << " ... done" << endl;	
							
							//ok
							cout << "plug_hi_cnt_left ";
							remove("plug_hi_cnt_left.stl");	
							stl.open("plug_hi_cnt_left.stl");
							stl << "solid plug_hi_cnt_left" << endl;
							stl_surface(+1,-1,-1,stl,0,cmbrk);
							stl_contour(-1,-1,+1,-1,stl,0,cmbrk,1,1);
							stl_vertical(+1,-1,stl,0,0,-1,-1,0);
							stl_vertical(-1,-1,stl,0,0,+2,-1,-1);
							stl_vertical(-1,-1,stl,cmbrk,0,-1,-1,0);
							stl_vertical(+1,-1,stl,cmbrk,0,+2,-1,+1);
							stl_terminal(-1,-1,stl,0,-1,-1);
							stl_terminal(+1,-1,stl,cmbrk,+1,-1);
							stl << "endsolid plug_hi_cnt_left" << endl;
							stl.close();
							cout << " ... done" << endl;				
							
							//ok
							cout << "plug_hi_cnt_right ";
							remove("plug_hi_cnt_right.stl");
							stl.open("plug_hi_cnt_right.stl");
							stl << "solid plug_hi_cnt_right" << endl; 
							stl_surface(-1,+1,-1,stl,0,cmbrk);
							stl_contour(+1,+1,+1,-1,stl,0,cmbrk,1,1);
							stl_vertical(-1,+1,stl,0,0,-1,-1,0);
							stl_vertical(+1,+1,stl,0,0,+2,-1,-1);
							stl_vertical(+1,+1,stl,cmbrk,0,-1,-1,0);
							stl_vertical(-1,+1,stl,cmbrk,0,+2,-1,+1);
							stl_terminal(+1,+1,stl,0,-1,-1);
							stl_terminal(-1,+1,stl,cmbrk,+1,-1);
							stl << "endsolid plug_hi_cnt_right" << endl;
							stl.close();
							cout << " ... done" << endl;	
							
							//ok
							cout << "plug_hi_mid_left ";
							remove("plug_hi_mid_left.stl");
							stl.open("plug_hi_mid_left.stl");
							stl << "solid plug_hi_mid_left" << endl;
							stl_surface(+1,-1,-1,stl,cmbrk,mtbrk);
							stl_contour(-1,-1,+1,-1,stl,cmbrk,mtbrk,1,1);
							stl_vertical(+1,-1,stl,cmbrk,0,-1,-1,0);
							stl_vertical(-1,-1,stl,cmbrk,0,+2,-1,-1);
							stl_vertical(-1,-1,stl,mtbrk,0,-1,-1,0);
							stl_vertical(+1,-1,stl,mtbrk,0,+2,-1,+1);
							stl_terminal(-1,-1,stl,cmbrk,-1,-1);
							stl_terminal(+1,-1,stl,mtbrk,+1,-1);
							stl << "endsolid plug_hi_mid_left" << endl;	
							stl.close();
							cout << " ... done" << endl;				
							
							//ok
							cout << "plug_hi_mid_right ";
							remove("plug_hi_mid_right.stl");
							stl.open("plug_hi_mid_right.stl");
							stl << "solid plug_hi_mid_right" << endl; 
							stl_surface(-1,+1,-1,stl,cmbrk,mtbrk);
							stl_contour(+1,+1,+1,-1,stl,cmbrk,mtbrk,1,1);
							stl_vertical(-1,+1,stl,cmbrk,0,-1,-1,0);
							stl_vertical(+1,+1,stl,cmbrk,0,+2,-1,-1);
							stl_vertical(+1,+1,stl,mtbrk,0,-1,-1,0);
							stl_vertical(-1,+1,stl,mtbrk,0,+2,-1,+1);
							stl_terminal(+1,+1,stl,cmbrk,-1,-1);
							stl_terminal(-1,+1,stl,mtbrk,+1,-1);
							stl << "endsolid plug_hi_mid_right" << endl;
							stl.close();
							cout << " ... done" << endl;	
							
							//ok
							cout << "plug_hi_tip_left ";
							remove("plug_hi_tip_left.stl");	
							stl.open("plug_hi_tip_left.stl");
							stl << "solid plug_hi_tip_left" << endl;
							stl_surface(+1,-1,-1,stl,mtbrk,mybrk);
							stl_contour(-1,-1,+1,-1,stl,mtbrk,mybrk,1,0);
							stl_tipclose(-1,-1,+1,stl);
							stl_vertical(+1,-1,stl,mtbrk,0,-1,-1,0);
							stl_vertical(-1,-1,stl,mtbrk,0,+2,-1,-1);
							stl_terminal(-1,-1,stl,mtbrk,-1,-1);
							stl << "endsolid plug_hi_tip_left" << endl;
							stl.close();
							cout << " ... done" << endl;		
							
							//ok
							cout  << "plug_hi_tip_right ";
							remove("plug_hi_tip_right.stl");
							stl.open("plug_hi_tip_right.stl");
							stl << "solid plug_hi_tip_right" << endl;		
							stl_surface(-1,+1,-1,stl,mtbrk,mybrk);
							stl_contour(+1,+1,+1,-1,stl,mtbrk,mybrk,1,0);
							stl_tipclose(+1,+1,+1,stl);
							stl_vertical(-1,+1,stl,mtbrk,0,-1,-1,0);
							stl_vertical(+1,+1,stl,mtbrk,0,+2,-1,-1);
							stl_terminal(+1,+1,stl,mtbrk,-1,-1);
							stl << "endsolid plug_hi_tip_right" << endl;
							stl.close();
							cout << " ... done" << endl;	
							
							//ok
							cout << "plug_lo_central ";
							remove("plug_lo_central.stl");
							stl.open("plug_lo_central.stl");
							stl << "solid plug_lo_central" << endl;		 
							stl_surface(-1,+1,+1,stl,0,cmbrk);
							stl_surface(+1,-1,+1,stl,0,cmbrk);
							stl_contour(-1,+1,-1,-1,stl,0,cmbrk,0,1);
							stl_contour(+1,-1,-1,-1,stl,0,cmbrk,1,1);
							stl_vertical(-1,+1,stl,cmbrk,0,+1,+1,0);
							stl_vertical(+1,-1,stl,cmbrk,0,+1,+1,0);
							stl_vertical(+1,+1,stl,cmbrk,0,-2,+1,+1);
							stl_vertical(-1,-1,stl,cmbrk,0,-2,+1,+1);
							stl_terminal(+1,+1,stl,cmbrk,+1,+1);
							stl_terminal(-1,-1,stl,cmbrk,+1,+1);		
							stl  << "endsolid plug_lo_central" << endl;
							stl.close();
							cout << " ... done" << endl;	
												
							//ok
							cout << "plug_lo_cnt_left ";
							remove("plug_lo_cnt_left.stl"); 
							stl.open("plug_lo_cnt_left.stl");
							stl << "solid plug_lo_cnt_left" << endl;		
							stl_surface(+1,-1,+1,stl,0,cmbrk);
							stl_contour(+1,-1,-1,-1,stl,0,cmbrk,1,1);
							stl_vertical(-1,-1,stl,0,0,+1,+1,0);
							stl_vertical(+1,-1,stl,0,0,-2,+1,-1);
							stl_vertical(+1,-1,stl,cmbrk,0,+1,+1,0);
							stl_vertical(-1,-1,stl,cmbrk,0,-2,+1,+1);
							stl_terminal(+1,-1,stl,0,-1,+1);
							stl_terminal(-1,-1,stl,cmbrk,+1,+1);		
							stl << "endsolid plug_lo_cnt_left" << endl;
							stl.close();
							cout << " ... done" << endl;	
												
							//ok
							cout << "plug_lo_cnt_right ";
							remove("plug_lo_cnt_right.stl"); 
							stl.open("plug_lo_cnt_right.stl");	
							stl << "solid plug_lo_cnt_right" << endl;		
							stl_surface(-1,+1,+1,stl,0,cmbrk);
							stl_contour(-1,+1,-1,-1,stl,0,cmbrk,1,1);
							stl_vertical(+1,+1,stl,0,0,+1,+1,0);
							stl_vertical(-1,+1,stl,0,0,-2,+1,-1);
							stl_vertical(-1,+1,stl,cmbrk,0,+1,+1,0);
							stl_vertical(+1,+1,stl,cmbrk,0,-2,+1,+1);
							stl_terminal(-1,+1,stl,0,-1,+1);
							stl_terminal(+1,+1,stl,cmbrk,+1,+1);		
							stl << "endsolid plug_lo_cnt_right" << endl;
							stl.close();
							cout << " ... done" << endl;	
												
							//ok
							cout << "plug_lo_mid_left ";
							remove("plug_lo_mid_left.stl");
							stl.open("plug_lo_mid_left.stl");
							stl << "solid plug_lo_mid_left" << endl;		
							stl_surface(+1,-1,+1,stl,cmbrk,mtbrk);
							stl_contour(+1,-1,-1,-1,stl,cmbrk,mtbrk,1,1);
							stl_vertical(-1,-1,stl,cmbrk,0,+1,+1,0);
							stl_vertical(+1,-1,stl,cmbrk,0,-2,+1,-1);
							stl_vertical(+1,-1,stl,mtbrk,0,+1,+1,0);
							stl_vertical(-1,-1,stl,mtbrk,0,-2,+1,+1);
							stl_terminal(+1,-1,stl,cmbrk,-1,+1);
							stl_terminal(-1,-1,stl,mtbrk,+1,+1);		
							stl << "endsolid plug_lo_mid_left" << endl;
							stl.close();
							cout << " ... done" << endl;	
							
							//ok
							cout << "plug_lo_mid_right ";
							remove("plug_lo_mid_right.stl");
							stl.open("plug_lo_mid_right.stl");
							stl << "solid plug_lo_mid_right" << endl;		
							stl_surface(-1,+1,+1,stl,cmbrk,mtbrk);
							stl_contour(-1,+1,-1,-1,stl,cmbrk,mtbrk,1,1);
							stl_vertical(+1,+1,stl,cmbrk,0,+1,+1,0);
							stl_vertical(-1,+1,stl,cmbrk,0,-2,+1,-1);
							stl_vertical(-1,+1,stl,mtbrk,0,+1,+1,0);
							stl_vertical(+1,+1,stl,mtbrk,0,-2,+1,+1);
							stl_terminal(-1,+1,stl,cmbrk,-1,+1);
							stl_terminal(+1,+1,stl,mtbrk,+1,+1);	
							stl << "endsolid plug_lo_mid_right" << endl;
							stl.close();
							cout << " ... done" << endl;	
							
							//ok
							cout << "plug_lo_tip_left ";
							remove("plug_lo_tip_left.stl");
							stl.open("plug_lo_tip_left.stl");
							stl << "solid plug_lo_tip_left" << endl;		
							stl_surface(+1,-1,+1,stl,mtbrk,mybrk);
							stl_contour(+1,-1,-1,-1,stl,mtbrk,mybrk,1,0);
							stl_vertical(-1,-1,stl,mtbrk,0,+1,+1,0);
							stl_vertical(+1,-1,stl,mtbrk,0,-2,+1,-1);
							stl_terminal(+1,-1,stl,mtbrk,-1,+1);		
							stl_tipclose(+1,-1,-1,stl);
							stl << "endsolid plug_lo_tip_left" << endl;
							stl.close();
							cout << " ... done" << endl;	
							
							//ok
							cout << "plug_lo_tip_right ";
							remove("plug_lo_tip_right.stl");
							stl.open("plug_lo_tip_right.stl");
							stl << "solid plug_lo_tip_right" << endl;					
							stl_surface(-1,+1,+1,stl,mtbrk,mybrk);
							stl_contour(-1,+1,-1,-1,stl,mtbrk,mybrk,1,0);
							stl_vertical(+1,+1,stl,mtbrk,0,+1,+1,0);
							stl_vertical(-1,+1,stl,mtbrk,0,-2,+1,-1);
							stl_terminal(-1,+1,stl,mtbrk,-1,+1);		
							stl_tipclose(-1,+1,-1,stl);
							stl << "endsolid plug_lo_tip_right" << endl;
							stl.close(); 
							cout << " ... done" << endl;	
		
				};
			cout << endl << endl << "All done ... have fun with the STL files ... " << endl;
		}
	else 
		cout << "Error on opening xflrwing.def and / or xflrwing.cfg"; 
	return 0;
};
