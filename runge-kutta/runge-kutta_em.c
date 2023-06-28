#include<stdio.h>
#include<math.h>

#define M_PI 3.14159265358979323846	//pi
#define ec 1.6021766208e-19			//elementary charge
#define mass 1.672621898e-27			//weight of proton
#define NN 1000

int Runge_Kutta_RK4(double bdt,double bh0[],double bh1[]);
void func_RK4(double ba0[],double ba1[]);

int main(void)
{
	int i,j;
	double bh0[7],bh1[7];
	double gdt;	//time step, weight of particle
	double t[NN];
	double rx0[NN],ry0[NN],rz0[NN],px0[NN],py0[NN],pz0[NN];
	FILE *fp;

	gdt=1.0e-9; //time step

	for(i=0;i<NN;i++){
		t[i]=0.0e0;
		rx0[i]=0.0e0; ry0[i]=0.0e0; rz0[i]=0.0e0;
		px0[i]=0.0e0; py0[i]=0.0e0; pz0[i]=0.0e0;
	}

	rx0[0]=0.0e0; ry0[0]=0.0e0; rz0[0]=0.0e0;
	px0[0]=mass*1.0e6; py0[0]=0.0e0; pz0[0]=0.0e0;

	for(i=0;i<7;i++){bh0[i]=0.0e0; bh1[i]=0.0e0;}

	bh0[0]=t[0];
	bh0[1]=rx0[0]; bh0[2]=ry0[0]; bh0[3]=rz0[0];
	bh0[4]=px0[0]; bh0[5]=py0[0]; bh0[6]=pz0[0];

	for(i=1;i<NN;i++){
		Runge_Kutta_RK4(gdt,bh0,bh1);
		t[i]=bh1[0];
		rx0[i]=bh1[1]; ry0[i]=bh1[2]; rz0[i]=bh1[3];
		px0[i]=bh1[4]; py0[i]=bh1[5]; pz0[i]=bh1[6];
		for(j=0;j<7;j++) bh0[j]=bh1[j];
	}
	fp=fopen("data_r0.csv","w");
	for(i=0;i<NN;i++) fprintf(fp,"%12.4e,%12.4e,%12.4e,%12.4e\n",t[i],rx0[i],ry0[i],rz0[i]);
	fclose(fp);
	fp=fopen("data_p0.csv","w");
	for(i=0;i<NN;i++) fprintf(fp,"%12.4e,%12.4e,%12.4e,%12.4e\n",t[i],px0[i],py0[i],pz0[i]);
	fclose(fp);

	return 0;
}

int Runge_Kutta_RK4(double dt,double bh0[],double bh1[]) //Main function of Runge-Kutta
{
	int i,j;
	double dt2,dt6;			//dt/2, dt/6
	double bufs[7],bufr[7];	//send/recieve data array
	double bufk[4][7];		//[No.1~4] data array for kx

	dt2=dt/2.0e0;
	dt6=dt/6.0e0;
	for(i=0;i<7;i++){bufs[i]=0.0e0; bufr[i]=0.0e0;}
	for(i=0;i<4;i++){for(j=0;j<7;j++) bufk[i][j]=0.0e0;}

	bufs[0]=bh0[0];									//1st step for k1
	for(i=1;i<7;i++) bufs[i]=bh0[i];
	func_RK4(bufs,bufr);
	for(i=0;i<7;i++) bufk[0][i]=bufr[i];

	bufs[0]=bh0[0]+dt2;							//2nd step for k2
	for(i=1;i<7;i++) bufs[i]=bh0[i]+dt2*bufr[i];
	func_RK4(bufs,bufr);
	for(i=0;i<7;i++) bufk[1][i]=bufr[i];

	bufs[0]=bh0[0]+dt2;							//3rd step for k3
	for(i=1;i<7;i++) bufs[i]=bh0[i]+dt2*bufr[i];
	func_RK4(bufs,bufr);
	for(i=0;i<7;i++) bufk[2][i]=bufr[i];

	bufs[0]=bh0[0]+dt;								//4th step for k4
	for(i=1;i<7;i++) bufs[i]=bh0[i]+dt*bufr[i];
	func_RK4(bufs,bufr);
	for(i=0;i<7;i++) bufk[3][i]=bufr[i];

	bh1[0]=bh0[0]+dt;
	for(i=1;i<7;i++){
		bh1[i]=bh0[i]+dt6*(bufk[0][i]+2.0e0*bufk[1][i]+2.0e0*bufk[2][i]+bufk[3][i]);
	}

	return 0;
}

void func_RK4(double ba0[],double ba1[]) //function for update
{
	double ex,ey,ez,bx,by,bz,m_inv;

	ex=2.0e4;
	ey=5.0e4;
	ez=2.0e4;
	bx=0.0e0;
	by=0.0e0;
	bz=1.0e0;
	m_inv=1.0e0/mass;
	ba1[0]=ba0[0];
	ba1[1]=m_inv*ba0[4];
	ba1[2]=m_inv*ba0[5];
	ba1[3]=m_inv*ba0[6];
	ba1[4]=ec*(ex+bz*ba1[2]-by*ba1[3]);
	ba1[5]=ec*(ey+bx*ba1[3]-bz*ba1[1]);
	ba1[6]=ec*(ez+by*ba1[1]-bx*ba1[2]);

	return;
}
