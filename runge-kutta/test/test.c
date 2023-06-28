/*
 * rk2ex2.c (Runge-Kutta method for free fall)
 */

#include <stdio.h>

#define M_PI 3.14159265358979323846	//pi
#define ec 1.6021766208e-19			//elementary charge
#define mass 1.672621898e-27			//weight of proton
#define NN 2000

double g = 9.8;
double runge_kutta(double *r[][NN], double *v[][NN], double *t, double dt, double N);
//int main (void);
double f(double x, double y, double t);


void func_RK4(double t, double v[], double r[])
{
	double ex,ey,ez,bx,by,bz,bm0,btm,bfz,bx,by,bz,bd;

	ex=2.0e4;
	ey=0.0e4;
	ez=1.0e4;
  bx=0.0e0;
  by=1.0e0;
	bz=0.0e0;

	ba1[4]=ec*(ex+bz*v[2]-by*v[3]);
	ba1[5]=ec*(ey+bx*v[3]-bz*v[1]);
	ba1[6]=ec*(ez+by*v[1]-bz*v[2]);

	return;
}

double runge_ketta(double *r, double *v, double *t, double dt, double N)
{
  double k1[3], l1[3], k2[3], l2[3], k3[3], l3[3], k4[3], l4[3];
  // Runge-Kutta 法

  for (int i = 1; i < NN-1; i++) {
    for (int j=0; j<3; j++) {
      k1[j] = dt * f(*r[j][i-1], *v[j][i-1], *t[j][i-1]);
      l1[j] = dt * (*v[j][i-1]);
      k2[j] = dt * f(*r[j][i-1] + k1[]j / 2, *v[j][i-1] + l1[j] / 2, *t[j][i-1] + dt / 2);
      l2[j] = dt * (*v[j][i-1] + l1[j] / 2);
      k3[j] = dt * f(*r[j][i-1] + k2[j] / 2, *v[j][i-1] + l2[j] / 2, *t[j][i-1] + dt / 2);
      l3[j] = dt * (*v[j][i-1] + l2[j] / 2);
      k4[j] = dt * f(*r[j][i-1] + k3[j], *v[j][i-1] + l3[j], *t[j][i-1] + dt);
      l4[j] = dt * (*v + l3);
    }
    for (int j=0; j<3; j++) {
    *v[j][i] = *v[j][i-1]+(k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]) / 6;
    *r[j][i] = *r[j][i-1](l1[j] + 2 * l2[j] + 2 * l3[j] + l4[j]) / 6;
    *t[j][i] = *t[j][i-1]+ dt;
    }
  }
  return 0;
}

void init_setting(void)
{
	int i;

	for(i=0;i<NN;i++){
		tt[i]=0.0e0;
		rx0[i]=0.0e0; ry0[i]=0.0e0; rz0[i]=0.0e0;
		px0[i]=0.0e0; py0[i]=0.0e0; pz0[i]=0.0e0;
	}

	rx0[0]=0.0e0; ry0[0]=0.0e0; rz0[0]=0.0e0
	px0[0]=mass*1.0e6; py0[0]=0.0e0; pz0[0]=0.0e0

	return;
}

int main(void)
{
  int i, N;
  double t, r, v, dt, k1, l1, k2, l2, k3, l3, k4, l4;
  double f(double, double, double), x0, y0;
  double Tmax;
  double gdt;	//time step, weight of particle
  double t[NN];
  double r[3][NN], v[3][NN];

  // 初期値 (100m の高さから初速 0 で降りる)
  x0 = 100.0; y0 = 0.0;
  // 分割数, 最終時刻 -> 時間刻み
  printf("# N, Tmax: "); scanf("%d%lf", &N, &Tmax);
  dt = Tmax / N;
  // 初期値
  t = 0.0;
  
  runge_ketta(&r[][NN], &v[][NN], &t, dt, N);

  printf("# t v r\n");
  printf("%f %f %f\n", t, v, r);
  
  return 0;
}
