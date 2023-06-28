#include<stdio.h>
#include<math.h>

//#define M_PI 3.14159265358979323846	//pi
#define ec 1.6021766208e-19			//elementary charge
#define mass 1.672621898e-27			//weight of proton
#define NN 1000

int Runge_Kutta_RK4(double dt,double data_arr[]);
void func_RK4(double buf[]);

int main(void)
{
	int i,j;
	double data_arr[7];
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

	for(i=0;i<7;i++){data_arr[i]=0.0e0;}

	data_arr[0]=t[0];
	data_arr[1]=rx0[0]; data_arr[2]=ry0[0]; data_arr[3]=rz0[0];
	data_arr[4]=px0[0]; data_arr[5]=py0[0]; data_arr[6]=pz0[0];

	for(i=1;i<NN;i++){
		Runge_Kutta_RK4(gdt,data_arr);
		t[i]=data_arr[0];
		rx0[i]=data_arr[1]; ry0[i]=data_arr[2]; rz0[i]=data_arr[3];
		px0[i]=data_arr[4]; py0[i]=data_arr[5]; pz0[i]=data_arr[6];
	}
	fp=fopen("data_r0.csv","w");
	for(i=0;i<NN;i++) fprintf(fp,"%12.4e,%12.4e,%12.4e,%12.4e\n",t[i],rx0[i],ry0[i],rz0[i]);
	fclose(fp);
	fp=fopen("data_p0.csv","w");
	for(i=0;i<NN;i++) fprintf(fp,"%12.4e,%12.4e,%12.4e,%12.4e\n",t[i],px0[i],py0[i],pz0[i]);
	fclose(fp);

	return 0;
}

int Runge_Kutta_RK4(double dt,double data_arr[]) //Main function of Runge-Kutta
{
	int i,j;
	double dt2,dt6;			//dt/2, dt/6
	double buf[7];	//send/recieve data array
	double bufk[4][7];		//[No.1~4] data array for kx

	dt2=dt/2.0e0;
	dt6=dt/6.0e0;
	for(i=0;i<7;i++){buf[i]=0.0e0;}
	for(i=0;i<4;i++){for(j=0;j<7;j++) bufk[i][j]=0.0e0;}

	buf[0]=data_arr[0];									//1st step for k1
	for(i=1;i<7;i++) buf[i]=data_arr[i];
	func_RK4(buf);
	for(i=0;i<7;i++) bufk[0][i]=buf[i];

	buf[0]=data_arr[0]+dt2;							//2nd step for k2
	for(i=1;i<7;i++) buf[i]=data_arr[i]+dt2*buf[i];
	func_RK4(buf);
	for(i=0;i<7;i++) bufk[1][i]=buf[i];

	buf[0]=data_arr[0]+dt2;							//3rd step for k3
	for(i=1;i<7;i++) buf[i]=data_arr[i]+dt2*buf[i];
	func_RK4(buf);
	for(i=0;i<7;i++) bufk[2][i]=buf[i];

	buf[0]=data_arr[0]+dt;								//4th step for k4
	for(i=1;i<7;i++) buf[i]=data_arr[i]+dt*buf[i];
	func_RK4(buf);
	for(i=0;i<7;i++) bufk[3][i]=buf[i];

	data_arr[0]=data_arr[0]+dt;
	for(i=1;i<7;i++){
		data_arr[i]=data_arr[i]+dt6*(bufk[0][i]+2.0e0*bufk[1][i]+2.0e0*bufk[2][i]+bufk[3][i]);
	}

	return 0;
}

void func_RK4(double buf[]) //function for update
{
	double ex,ey,ez,bx,by,bz,m_inv;

	ex=2.0e4;
	ey=5.0e4;
	ez=2.0e4;
	bx=0.0e0;
	by=0.0e0;
	bz=1.0e0;
	m_inv=1.0e0/mass;
	buf[0]=buf[0];
	buf[1]=m_inv*buf[4];
	buf[2]=m_inv*buf[5];
	buf[3]=m_inv*buf[6];
	buf[4]=ec*(ex+bz*buf[2]-by*buf[3]);
	buf[5]=ec*(ey+bx*buf[3]-bz*buf[1]);
	buf[6]=ec*(ez+by*buf[1]-bx*buf[2]);

	return;
}
