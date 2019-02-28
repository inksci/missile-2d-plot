//
//       通用弹道仿真计算程序(经典版)V1.0
//
//
/*

通用弹道仿真计算程序V1.0
本程序包为vxLand提供的经典版的通用导弹、炮弹、火箭弹等的弹道仿真计算程序，适用于无控弹药， 
用户可以需要根据加入三点法、比例导引等导引方法，修改、扩充成相应的有控、制导模拟仿真弹道程 
序。

程序包中MissileTrajectorySimulationRk.h的rk()是经典的4阶龙格库塔方法求常微分方程的数值算法 
的函数子程序，大家可参考使用。

由于时间关系，相应的matlab语言版本暂时没有整理出来。

关键词：导弹 弹道 炮弹 火箭弹 三点法 比例导引 无控 有控 制导 仿真 模拟 计算 通用

MissileTrajectorySimulation_Input.dat是输入数据文件

*/

/******************************************
*   (c) Copyright vxLand 2006,2007,2008
*	y(0)	-	t	(s)
*	y(1)	-	v	(m/s)
*	y(2)	-	Theta (degree)
*	y(3)	-	x	(m)
*	y(4)	-	y	(m)
*	y(5)	-	m   (kg)
*******************************************/

#include "stdafx.h"
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <conio.h>
#include "MissileTrajectorySimulation.h"

double acx[22][2],acy[10][2],dy[NumbVars];

double y[NumbVars],y00[NumbVars];

double s,alpha;
void read_array(double *array,int n1,int n2,FILE *fp);
double insert(double w,int n,double ax[][2]);
void dery(double xy[NumbVars]);

#include "MissileTrajectorySimulationRk.h"

void dery(double xy[NumbVars]){
	double xf,yf,mach,q,rho;
	const double T0=289.05,p0=99991.5,R=287.053;
	double T,ps,sonic;

	T = T0 - 0.0063 * xy[4];
	ps = p0 * pow( T/T0 , 5.4246);
	sonic = sqrt( 1.4*R*T );
	rho = ps / ( R * T );
	q = 0.5*rho*xy[1]*xy[1];
	mach = xy[1] / sonic ;

	alpha = 0.0;
	xf=q*s*insert(mach,22,acx);
	yf=q*s*alpha*insert(mach,10,acy);
	dy[0] = 1.0 ;
	dy[1] = (0.0-xf-G*xy[5]*sin(xy[2]))/xy[5];
	dy[2] = yf/(xy[5]*xy[1])-G*cos(xy[2])/xy[1];
	dy[3] = xy[1]*cos(xy[2]);
	dy[4] = xy[1]*sin(xy[2]);
	dy[5] = 0.0;
}

void main_missile() {
  int i;
  double v_time=0.0;
  FILE *fp;

  fp=fopen("MissileTrajectorySimulation_Output.dat","w"); // 输出结果文件
  read_data();
  y00[2] = y00[2] * ToRadian;
  for(i=0;i<NumbVars;i++) y[i]=y00[i];
  fprintf(fp,"  T(s)   V(m/s)   Theta     X(m)     Y(m)\n");
  while( y[4] >= 0. )  { // 当高度大于0
	fprintf(fp,"%6.2f%9.2f%8.2f%10.3f%10.3f\n",
	  v_time,y[1],y[2]*ToDegree,y[3],y[4]); // 写入到输出数据文件
	rk(RUN_PACE);
	v_time+=RUN_PACE; // 累积时间   			}

	fprintf(fp,"%6.2f%9.2f%8.2f%10.3f%10.3f\n",
	  v_time,y[1],y[2]*ToDegree,y[3],y[4]);
	fclose(fp);
	printf("通用弹道仿真计算程序(经典版)V1.0    BY vxLand   2008\n\n计算完毕!\n");getch();
	}
}
double insert(double w,int n,double ax[][2]) {
  double x;int i;
  if( w < ax[0][0] || w > ax[n-1][0])
{ printf("\tERROR: Mach is out! \n\tMach = %10.3f , V = %9.2f, t = %6.2f\n",w,y[1],y[0]);
		exit(-1); }
  for(i=0;w>ax[i][0];) i++;
  if(i==0) i=1;
  if(i==n-1) i=n-2;
  x=ax[i-1][1]*(w-ax[i][0])*(w-ax[i+1][0])
		 /((ax[i-1][0]-ax[i][0])*(ax[i-1][0]-ax[i+1][0])) +
	 ax[i][1]*(w-ax[i-1][0])*(w-ax[i+1][0])
		 /((ax[i][0]-ax[i-1][0])*(ax[i][0]-ax[i+1][0])) +
	 ax[i+1][1]*(w-ax[i-1][0])*(w-ax[i][0])
		 /((ax[i+1][0]-ax[i-1][0])*(ax[i+1][0]-ax[i][0]));
  return(x); }


void read_array(double *array,int n1,int n2,FILE *fp){
  int i,j;float q;
  for(i=0;i<n1;i++)
  for(j=0;j<n2;j++) { fscanf(fp,"%f",&q); *array++ = (double)q;} }

void read_data(void){
  double d;
	float v1;
  FILE *fp;
  fp=fopen("MissileTrajectorySimulation_Input.dat","r"); // 输入数据文件
  read_array((double *)acx,22,2,fp);
  read_array((double *)acy,10,2,fp);
  fscanf(fp,"%f",&v1); d= (float)v1;
  s=d*d*3.1415926/4.;
  read_array(y00,1,NumbVars,fp);
  fclose(fp); }

