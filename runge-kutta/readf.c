#include <stdio.h>
#include <stdlib.h>
 
int main (void)
{
int i,j,row,col;

row=501;
col=1002;
double** val = NULL;
double** ex_amp = NULL;
double** ex_phase = NULL;
val = (double**)malloc(sizeof(double*) * row);
ex_amp = (double**)malloc(sizeof(double*) * row);
ex_phase = (double**)malloc(sizeof(double*) * row);
if(val == NULL){
    return -1;
}
for(i = 0; i < row; i++){
    val[i] = (double*)malloc(sizeof(double) * col);
    if(val[i] == NULL){
        return -1;
    }
    }
for(i = 0; i < row; i++){
    ex_amp[i] = (double*)malloc(sizeof(double) * row);
    ex_phase[i] = (double*)malloc(sizeof(double) * row);
    if(val[i] == NULL){
        return -1;
    }
    }   
FILE *fp;

fp=fopen("ex2.dat","r");

if(fp==NULL){
printf("ファイルがありません\n");
return -1;
}
for(i=0;i<501;i++)for(j=0;j<1002;j++)fscanf(fp,"%lf",&val[i][j]);

//for(j=0;j<1002;j++)printf("%lf ",val[0][j]);
//printf("\n");

for(i=0; i<row; i++){
    for(j=0; i<row; i++){
        ex_amp[i][j]= val[i][2*j];
        ex_phase[i][j]= val[i][2*j+1];
    }
}

fclose(fp);
free(val);
free(ex_amp);
free(ex_phase);


printf(sizeof(*ex_amp)/sizeof(*ex_amp[0]));
printf(sizeof(*ex_amp[0])/sizeof(*ex_amp[0][0]));

return 0;


}

 

