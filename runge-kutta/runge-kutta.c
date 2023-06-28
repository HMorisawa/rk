#include<stdio.h>
#include<math.h>

#define M_PI 3.14159265358979323846	//pi
#define gc 6.67408e-11				//gravitational constant
#define ec 1.6021766208e-19			//elementary charge
#define mp 1.672621898e-27			//weight of proton
#define pe 8.85418782e-12			//permittivity of vacuum
#define moe 5.972e24				//mass of Earth
#define gae 9.80665e0				//gravitational acceleration of Earth
#define GMe 3.986004418e14			//earth's center constant of gravitation (=gc*moe)
#define NN 2000

int mod,nop;		//mode, number of particle
double gdt,mass;	//time step, weight of particle
double tt[NN];
double rx0[NN],ry0[NN],rz0[NN],px0[NN],py0[NN],pz0[NN];
double rx1[NN],ry1[NN],rz1[NN],px1[NN],py1[NN],pz1[NN];

void init_setting(void);
int Runge_Kutta_RK4(int bn,double bdt,double bh0[],double bh1[]);
void func_RK4(double ba0[],double ba1[]);

int main(void)
{
	int i,j,bn;
	double bh0[15],bh1[15];
	FILE *fp;

	gdt=1.0e-9;
	mod=3;
	nop=1;
	bn=15;

	init_setting();
	for(i=0;i<bn;i++){bh0[i]=0.0e0; bh1[i]=0.0e0;}

	bh0[0]=tt[0];
	bh0[1]=rx0[0]; bh0[2]=ry0[0]; bh0[3]=rz0[0];
	bh0[4]=px0[0]; bh0[5]=py0[0]; bh0[6]=pz0[0];
	bh0[7]=rx1[0]; bh0[8]=ry1[0]; bh0[9]=rz1[0];
	bh0[10]=px1[0]; bh0[11]=py1[0]; bh0[12]=pz1[0];

	if(nop==1){//single body problem
		for(i=1;i<NN;i++){
			Runge_Kutta_RK4(7,gdt,bh0,bh1);
			tt[i]=bh1[0];
			rx0[i]=bh1[1]; ry0[i]=bh1[2]; rz0[i]=bh1[3];
			px0[i]=bh1[4]; py0[i]=bh1[5]; pz0[i]=bh1[6];
			for(j=0;j<7;j++) bh0[j]=bh1[j];
		}
		fp=fopen("data_r0.csv","w");
		for(i=0;i<NN;i++) fprintf(fp,"%12.4e,%12.4e,%12.4e,%12.4e\n",tt[i],rx0[i],ry0[i],rz0[i]);
		fclose(fp);
		fp=fopen("data_p0.csv","w");
		for(i=0;i<NN;i++) fprintf(fp,"%12.4e,%12.4e,%12.4e,%12.4e\n",tt[i],px0[i],py0[i],pz0[i]);
		fclose(fp);
	}

	if(nop==2){//two body problem
		for(i=1;i<NN;i++){
			Runge_Kutta_RK4(13,gdt,bh0,bh1);
			tt[i]=bh1[0];
			rx0[i]=bh1[1]; ry0[i]=bh1[2]; rz0[i]=bh1[3];
			px0[i]=bh1[4]; py0[i]=bh1[5]; pz0[i]=bh1[6];
			rx1[i]=bh1[7]; ry1[i]=bh1[8]; rz1[i]=bh1[9];
			px1[i]=bh1[10]; py1[i]=bh1[11]; pz1[i]=bh1[12];
			for(j=0;j<13;j++) bh0[j]=bh1[j];
		}
		fp=fopen("data_r0.csv","w");
		for(i=0;i<NN;i++) fprintf(fp,"%12.4e,%12.4e,%12.4e,%12.4e\n",tt[i],rx0[i],ry0[i],rz0[i]);
		fclose(fp);
		fp=fopen("data_p0.csv","w");
		for(i=0;i<NN;i++) fprintf(fp,"%12.4e,%12.4e,%12.4e,%12.4e\n",tt[i],px0[i],py0[i],pz0[i]);
		fclose(fp);
		fp=fopen("data_r1.csv","w");
		for(i=0;i<NN;i++) fprintf(fp,"%12.4e,%12.4e,%12.4e,%12.4e\n",tt[i],rx1[i],ry1[i],rz1[i]);
		fclose(fp);
		fp=fopen("data_p1.csv","w");
		for(i=0;i<NN;i++) fprintf(fp,"%12.4e,%12.4e,%12.4e,%12.4e\n",tt[i],px1[i],py1[i],pz1[i]);
		fclose(fp);
		//add
		fp=fopen("position.csv","w");
		for(i=0;i<NN;i++) fprintf(fp,"%12.4e,%12.4e,%12.4e,%12.4e\n",rx0[i],ry0[i],rx1[i],ry1[i]);
		fclose(fp);
	}

	return 0;
}

void init_setting(void)
{
	int i;
	double buf0;

	mass=1.0e0;
	for(i=0;i<NN;i++){
		tt[i]=0.0e0;
		rx0[i]=0.0e0; ry0[i]=0.0e0; rz0[i]=0.0e0;
		px0[i]=0.0e0; py0[i]=0.0e0; pz0[i]=0.0e0;
		rx1[i]=0.0e0; ry1[i]=0.0e0; rz1[i]=0.0e0;
		px1[i]=0.0e0; py1[i]=0.0e0; pz1[i]=0.0e0;
	}

	switch (mod)
	{
	case 1://parabolic motion
		mass=3.0e0;
		rx0[0]=0.0e0; rz0[0]=2.0e1;
		buf0=M_PI/3.0e0;
		px0[0]=mass*6.0e1*cos(buf0); pz0[0]=mass*6.0e1*sin(buf0);
		break;

	case 2://air resistance
		mass=5.0e0;
		rx0[0]=0.0e0; rz0[0]=8.0e2;
		buf0=M_PI/4.0e0;
		px0[0]=mass*2.0e1*cos(buf0); pz0[0]=mass*2.0e1*sin(buf0);
		break;

	case 3://charged particle
		mass=mp;
		rx0[0]=0.0e0; ry0[0]=0.0e0;
		px0[0]=mass*1.0e6; py0[0]=0.0e0;
		break;

	case 4://rocket launch
		mass=5.0e2;
		rz0[0]=0.0e0;
		pz0[0]=0.0e0;
		break;

	case 5://gravity assisted acceleration
		mass=1.5e3;
		rx0[0]=5.0e8; ry0[0]=4.0e7; rz0[0]=0.0e0;
		px0[0]=-mass*5.0e2; py0[0]=-mass*1.0e3; pz0[0]=0.0e0;

		rx1[0]=0.0e0; ry1[0]=0.0e0; rz1[0]=0.0e0;
		px1[0]=moe*3.0e4; py1[0]=0.0e0; pz1[0]=0.0e0;
		break;
	}

	return;
}

int Runge_Kutta_RK4(int bn,double bdt,double bh0[],double bh1[])
{
	int i,j,bmax;
	double bdt2,bdt6;			//dt/2, dt/6
	double bufs[50],bufr[50];	//send/recieve data array
	double bufg[4][50];			//[No.1~4] data array

	bmax=50;
	if(bn>bmax) return -1;

	bdt2=bdt/2.0e0;
	bdt6=bdt/6.0e0;
	for(i=0;i<bmax;i++){bufs[i]=0.0e0; bufr[i]=0.0e0;}
	for(i=0;i<4;i++){for(j=0;j<bmax;j++) bufg[i][j]=0.0e0;}

	bufs[0]=bh0[0];									//1st step
	for(i=1;i<bn;i++) bufs[i]=bh0[i];
	func_RK4(bufs,bufr);
	for(i=0;i<bn;i++) bufg[0][i]=bufr[i];

	bufs[0]=bh0[0]+bdt2;							//2nd step
	for(i=1;i<bn;i++) bufs[i]=bh0[i]+bdt2*bufr[i];
	func_RK4(bufs,bufr);
	for(i=0;i<bn;i++) bufg[1][i]=bufr[i];

	bufs[0]=bh0[0]+bdt2;							//3rd step
	for(i=1;i<bn;i++) bufs[i]=bh0[i]+bdt2*bufr[i];
	func_RK4(bufs,bufr);
	for(i=0;i<bn;i++) bufg[2][i]=bufr[i];

	bufs[0]=bh0[0]+bdt;								//4th step
	for(i=1;i<bn;i++) bufs[i]=bh0[i]+bdt*bufr[i];
	func_RK4(bufs,bufr);
	for(i=0;i<bn;i++) bufg[3][i]=bufr[i];

	bh1[0]=bh0[0]+bdt;
	for(i=1;i<bn;i++){
		bh1[i]=bh0[i]+bdt6*(bufg[0][i]+2.0e0*bufg[1][i]+2.0e0*bufg[2][i]+bufg[3][i]);
	}

	return 0;
}

void func_RK4(double ba0[],double ba1[])
{
	int i;
	double bms0,bms1;
	double bex,bey,bbz,bm0,btm,bfz,bx,by,bz,bd;
	double buf0;

	switch (mod)
	{
	case 1://parabolic motion
		bms0=1.0e0/mass;
		ba1[0]=ba0[0];
		ba1[1]=bms0*ba0[4];
		ba1[2]=bms0*ba0[5];
		ba1[3]=bms0*ba0[6];
		ba1[4]=0.0e0;
		ba1[5]=0.0e0;
		ba1[6]=-gae*mass;
		break;

	case 2://air resistance
		buf0=2.0e0;
		bms0=1.0e0/mass;
		ba1[0]=ba0[0];
		ba1[1]=bms0*ba0[4];
		ba1[2]=bms0*ba0[5];
		ba1[3]=bms0*ba0[6];
		ba1[4]=-buf0*ba1[1];
		ba1[5]=-buf0*ba1[2];
		ba1[6]=-buf0*ba1[3]-gae*mass;
		break;

	case 3://charged particle
		bex=2.0e4;
		bey=5.0e4;
		bbz=1.0e0;
		bms0=1.0e0/mass;
		ba1[0]=ba0[0];
		ba1[1]=bms0*ba0[4];
		ba1[2]=bms0*ba0[5];
		ba1[3]=bms0*ba0[6];
		ba1[4]=ec*(bex+bbz*ba1[2]);
		ba1[5]=ec*(bey-bbz*ba1[1]);
		ba1[6]=0.0e0;
		break;

	case 4://rocket launch
		btm=5.0e0;
		bm0=2.0e2;
		bfz=6.0e3;
		if(ba0[0]<=btm){
			buf0=ba0[0]*(bm0-mass)/btm+mass;
			bms0=1.0e0/buf0;
		}else{
			buf0=bm0;
			bms0=1.0e0/bm0;
			bfz=0.0e0;
		}
		ba1[0]=ba0[0];
		ba1[1]=bms0*ba0[4];
		ba1[2]=bms0*ba0[5];
		ba1[3]=bms0*ba0[6];
		ba1[4]=0.0e0;
		ba1[5]=0.0e0;
		ba1[6]=bfz-gae*buf0;
		break;

	case 5://gravity assisted acceleration
		bms0=1.0e0/mass;
		bms1=1.0e0/moe;
		bx=ba0[7]-ba0[1];
		by=ba0[8]-ba0[2];
		bz=ba0[9]-ba0[3];
		bd=sqrt(bx*bx+by*by+bz*bz);
		if(bd<=0.0e0){
			for(i=0;i<13;i++) ba1[i]=0.0e0;
			break;
		}
		bd=GMe*mass/(bd*bd*bd);
		ba1[0]=ba0[0];
		ba1[1]=bms0*ba0[4];
		ba1[2]=bms0*ba0[5];
		ba1[3]=bms0*ba0[6];
		ba1[4]=bd*bx;
		ba1[5]=bd*by;
		ba1[6]=bd*bz;
		ba1[7]=bms1*ba0[10];
		ba1[8]=bms1*ba0[11];
		ba1[9]=bms1*ba0[12];
		ba1[10]=-bd*bx;
		ba1[11]=-bd*by;
		ba1[12]=-bd*bz;
		break;
	}

	return;
}
