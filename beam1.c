#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <Python.h>
double Q2(double x, double *q, int **index_matrix, double hh);
double NQ(double z, int i, double hh);
double M2(double x, double *q, int **index_matrix, double hh);
double NM(double z, int i, double hh);
double Ntheta(double z, int i, double hh);
double theta2(double x, double *q, int **index_matrix, double hh);
double Nw(double z, int i, double hh);
double w2(double x, double *q, int **index_matrix, double hh);
int fill_matrix(double **A, double hh, double v, double EJ);
int fill_t(double **t, double k, double EJ);
int LUdecompose(double ***A, double ***L, double ***U, int M);
int LUsolve(double ***L, double ***U, double *B, double *X, int M);
double h(double x);
double _M_o(double x);
double _Q_o(double x);
double _M_1(double x, double hh);
double _P(double x, double hh);
double _q_start(double x, double hh);
double _q_end(double x, double hh);
double _R_a(double x, double hh);
double _R_b(double x, double hh);
double _R_c(double x, double hh);
double w(double x, double hh, double EJ, double M_o, double Q_o, double R_a, double R_b, double R_c, double q, double P, double M_1);
double theta_q_end(double x, double hh);
double theta_q_start(double x, double hh);
double theta_p(double x, double hh);
double theta_m(double x, double hh);
double theta_a(double x, double hh);
double theta_b(double x, double hh);
double theta_c(double x, double hh);
double theta(double x, double hh, double EJ, double M_o, double Q_o, double R_a, double R_b, double R_c, double q, double P, double M_1);
double M_q_end(double x, double hh);
double M_q_start(double x, double hh);
double M(double x, double hh, double M_o, double Q_o, double R_a, double R_b, double R_c, double q, double P, double M_1);
double Q(double x, double hh, double M_o, double Q_o, double R_a, double R_b, double R_c, double q, double P, double M_1);
int main(void)
{
  int i, j, N = 5, **index_matrix = NULL;
  double **A = NULL, **L = NULL, **U = NULL, *B = NULL, *X = NULL;
  double *data = NULL;
  double E , J, M_1, P, q, k, l, EJ, hh, qh;
  double R_a, R_b, R_c, M_o, Q_o;
  double t;
  FILE  *input, *output, *test;
  input = fopen("input.txt", "r");
  data = (double *)calloc(7, sizeof(double));
  for(i = 0; i < 7; ++i) fscanf(input, "%lf", &data[i]);
  //for(i = 0; i < 7; ++i) printf("%e ", data[i]); printf("\n");
  E = data[0]; J = data[1]; M_1 = data[2]; P = data[3]; q = data[4]; k = data[5]; l = data[6];
  fclose(input);
  A = (double **) calloc(N, sizeof(double *)); 
  L = (double **) calloc(N, sizeof(double *)); 
  U = (double **) calloc(N, sizeof(double *)); 
  X = (double *) calloc(N, sizeof(double));
  B = (double *) calloc(N, sizeof(double));
  for(i = 0; i < N; i++)
  {
    A[i] = (double *) calloc(N, sizeof(double)); 
    L[i] = (double *) calloc(N, sizeof(double));
    U[i] = (double *) calloc(N, sizeof(double));
  }
  hh = l/20.0;
  
  B[0] = M_1 - P*(14.0*hh) - q *0.5 *(9.0*hh)*(15.0*hh); 
  B[1] = (-1.0) * P - q * (9.0*hh) - (-1.0)*k/(E*J) * ((-1.0)*M_1*0.5*pow(9.0*l/10.0,2) + P*(1.0/6.0)*pow(14.0*hh, 3) + q*(1.0/24.0)*(pow(12.0*hh,4) - pow(3.0*hh,4)));
  B[2] = M_1*0.5*pow(2.0*hh,2); 
  B[3] = M_1*0.5*pow(9.0*hh,2) - P * (1.0/6.0)* pow(5.0*hh,3) - q * (1.0/24.0)*pow(3.0*hh,4);
  B[4] = M_1*0.5*pow(14.0*hh,2) - P * (1.0/6.0)* pow(10.0*hh,3)- q * (1.0/24.0)*pow(8.0*hh,4);
  
  A[0][0] = 1.0; A[0][1] = l; A[0][2] = 16.0*hh; A[0][3] = 9.0*hh; A[0][4] = 4.0*hh;
  A[1][0] = (-1.0)*k/(E*J)*(l*l/2.0); A[1][1] = 1 - k/(E*J)*(pow(l,3)/6.0); A[1][2] = 1 - k/(E*J)*(pow(16.0*hh,3)/6.0); A[1][3] = 1 - k/(E*J)*(pow(9.0*hh,3)/6.0); A[1][4] = 1 - k/(E*J)*(pow(4.0*hh,3)/6.0);
  A[2][0] = 0.5*pow(4.0*hh, 2); A[2][1] = (1.0/6.0)*pow(4.0*hh, 3);
  A[3][0] = 0.5*pow(11.0*hh, 2); A[3][1] = (1.0/6.0)*pow(11.0*hh, 3); A[3][2] = (1.0/6.0)*pow(7.0*hh, 3);
  A[4][0] = 0.5*pow(16.0*hh, 2); A[4][1] = (1.0/6.0)*pow(16.0*hh, 3); A[4][2] = (1.0/6.0)*pow(12.0*hh, 3); A[4][3] = (1.0/6.0)*pow(5.0*hh, 3);
  LUdecompose(&A, &L, &U, N);
  LUsolve(&L, &U, B, X, N);
  EJ = E*J;
  M_o = X[0]; Q_o = X[1]; R_a = X[2]; R_b = X[3]; R_c = X[4]; free(X);
  output = fopen("w_data.txt", "w");
  for(t = 0.0; t <= l; t = t + 0.01) fprintf(output, "%e %.15e\n", t, w(t,hh,EJ,M_o,Q_o,R_a,R_b,R_c,q,P,M_1));
  fclose(output);
  output = fopen("theta_data.txt", "w");
  for(t = 0.0; t <= l; t = t + 0.01) fprintf(output, "%e %.15e\n", t,theta(t,hh,EJ,M_o,Q_o,R_a,R_b,R_c,q,P,M_1));
  fclose(output);
  output = fopen("M_data.txt", "w");
  for(t = 0.0; t <= l; t = t + 0.01) fprintf(output, "%e %.15e\n", t, M(t,hh,M_o,Q_o,R_a,R_b,R_c,q,P,M_1));
  fclose(output);
  output = fopen("Q_data.txt", "w");
  for(t = 0.0; t <= l; t = t + 0.01) fprintf(output, "%e %.15e\n", t, Q(t,hh,M_o,Q_o,R_a,R_b,R_c,q,P,M_1));
  fclose(output);
  free(B);
  for(i = 0; i < N; i++)
  {
    free((A)[i]);
    free((L)[i]);
    free((U)[i]);
  }
  free(A); free(U); free(L);
  N = 42;
  A = (double **) calloc(N, sizeof(double *)); 
  L = (double **) calloc(N, sizeof(double *)); 
  U = (double **) calloc(N, sizeof(double *)); 
  X = (double *) calloc(N, sizeof(double));
  B = (double *) calloc(N, sizeof(double));
  for(i = 0; i < N; i++)
  {
    A[i] = (double *) calloc(N, sizeof(double)); 
    L[i] = (double *) calloc(N, sizeof(double));
    U[i] = (double *) calloc(N, sizeof(double));
  }
  fill_matrix(A, hh, k, EJ);
  qh = q * hh;
  B[5] =  M_1;
  B[12] = P;
  B[16] = 0.5 * qh;
  B[17] = 1.0/12.0 * qh * hh;
  for(i = 0; i < 8; ++i) B[18 + 2 * i] = qh;
  B[22] = 0.0;
  B[32] = 0.0;
  B[34] = 0.5 * qh;
  B[35] = 1.0/12.0 * qh * hh * (-1.0);
  LUdecompose(&A, &L, &U, N);
  LUsolve(&L, &U, B, X, N);
  index_matrix = (int **) calloc(20, sizeof(int *));
  for(i = 0; i < 20; i++)
  {
    index_matrix[i] = (int *) calloc(4, sizeof(int));
    for(j = 0; j < 4; j++)
    {
        index_matrix[i][j] = 2*i + j;
    }
  }
  output = fopen("w2_data.txt", "w");
  for(t = 0.0; t <= l; t = t + 0.01) fprintf(output, "%e %.15e\n", t, w2(t, X, index_matrix, hh));
  fclose(output);
  output = fopen("theta2_data.txt", "w");
  for(t = 0.0; t <= l; t = t + 0.01) fprintf(output, "%e %.15e\n", t, theta2(t, X, index_matrix, hh));
  fclose(output);
  output = fopen("M2_data.txt", "w");
  for(t = 0.0; t <= l; t = t + 0.01) fprintf(output, "%e %.15e\n", t, EJ * M2(t, X, index_matrix, hh));
  fclose(output);
  output = fopen("Q2_data.txt", "w");
  for(t = 0.0; t <= l; t = t + 0.01) fprintf(output, "%e %.15e\n", t, EJ * Q2(t, X, index_matrix, hh));
  fclose(output);
  Py_Initialize();
  test = _Py_fopen("test.py", "r");
  PyRun_SimpleFile(test, "test.py");
  Py_Finalize();
  return 0;
}
double Nw(double z, int i, double hh)
{
    if(i == 0) return 1 - 3*pow(z/hh,2) +2 * pow(z/hh,3);
    if(i == 1) return z - 2/hh * z*z + z * pow(z/hh,2);
    if(i == 2) return 3*pow(z/hh,2) - 2*pow(z/hh,3);
    if(i == 3) return - z*z/hh + z*pow(z/hh,2);
}
double Ntheta(double z, int i, double hh)
{
    if(i == 0) return -6.0/hh * (z/hh) + 6.0/hh * pow(z/hh,2);
    if(i == 1) return 1.0 - 4.0/hh * z + 3.0 * pow(z/hh,2);
    if(i == 2) return 6.0/(hh*hh) * z - 6.0/(hh*hh*hh) * (z * z);
    if(i == 3) return -2.0/hh * z + 3.0/(hh*hh) * (z * z);
}
double NM(double z, int i, double hh)
{
    if(i == 0) return -6.0/(hh*hh) + 12.0/(hh*hh*hh) * z;
    if(i == 1) return  -4.0/hh + 6.0/(hh*hh) * z;
    if(i == 2) return 6.0/(hh*hh) - 12.0/(hh*hh*hh) * z;
    if(i == 3) return -2.0/hh + 6.0/(hh*hh)  * z;
}
double NQ(double z, int i, double hh)
{
    if(i == 0) return 12.0/(hh*hh*hh);
    if(i == 1) return  6.0/(hh*hh);
    if(i == 2) return -12.0/(hh*hh*hh);
    if(i == 3) return 6.0/(hh*hh);;
}
double Q2(double x, double *q, int **index_matrix, double hh)
{
    double sum = 0.0, a = 0.0;
    int i = 0, j;
    
    if(x - 20*hh == 0.0)
    {
        for(j = 0; j < 4; j++) sum = sum + q[index_matrix[19][j]] * NQ(hh, j, hh);
        return sum;
    }
    while (!((a <= x) && (x < a + hh)))
    {
        a = a + hh;
        i = i + 1;
    }
    for(j = 0; j < 4; j++) sum = sum + q[index_matrix[i][j]] * NQ(x - a, j, hh);
    return sum;
}
double M2(double x, double *q, int **index_matrix, double hh)
{
    double sum = 0.0, a = 0.0;
    int i = 0, j;
    
    if(x - 20*hh == 0.0)
    {
        for(j = 0; j < 4; j++) sum = sum + q[index_matrix[19][j]] * NM(hh, j, hh);
        return sum;
    }
    while (!((a <= x) && (x < a + hh)))
    {
        a = a + hh;
        i = i + 1;
    }
    for(j = 0; j < 4; j++) sum = sum + q[index_matrix[i][j]] * NM(x - a, j, hh);
    return sum;
}
double theta2(double x, double *q, int **index_matrix, double hh)
{
    double sum = 0.0, a = 0.0;
    int i = 0, j;
    if(x == 0.0) return q[1];
    if(x - 20*hh == 0.0) return q[41];
    while (!((a <= x) && (x < a + hh)))
    {
        a = a + hh;
        i = i + 1;
    }
    for(j = 0; j < 4; j++) sum = sum + q[index_matrix[i][j]] * Ntheta(x - a, j, hh);
    return sum;
}
double w2(double x, double *q, int **index_matrix, double hh)
{
    double sum = 0.0, a = 0.0;
    int i = 0, j;
    if(x == 0.0) return q[0];
    if(x - 20*hh == 0.0) return q[40];
    while (!((a <= x) && (x < a + hh)))
    {
        a = a + hh;
        i = i + 1;
    }
    for(j = 0; j < 4; j++) sum = sum + q[index_matrix[i][j]] * Nw(x - a, j, hh);
    return sum;
}
int fill_t(double **t, double k, double EJ)
{
    t[0][0] = EJ*   12.0/pow(k, 3);
    t[0][1] =  EJ* 6.0/pow(k, 2);
    t[0][2] =  EJ*(-12.0)/pow(k, 3);
    t[0][3] =  EJ*6.0/pow(k, 2);
    
    t[1][0] =  EJ*6.0/pow(k, 2);
    t[1][1] =  EJ*4.0/pow(k, 1);
    t[1][2] =  EJ*(-6.0)/pow(k, 2);
    t[1][3] =  EJ*2.0/pow(k, 1);
    
    t[2][0] =  EJ*(-12.0)/pow(k, 3);
    t[2][1] =  EJ*(-6.0)/pow(k, 2);
    t[2][2] =  EJ*12.0/pow(k, 3); 
    t[2][3] =  EJ*(-6.0)/pow(k, 2);
    
    t[3][0] =  EJ*6.0/pow(k, 2);
    t[3][1] =  EJ*2.0/pow(k, 1);
    t[3][2] =  EJ*(-6.0)/pow(k, 2);
    t[3][3] =  EJ*4.0/pow(k, 1);
    return 0;
}
int fill_matrix(double **A, double hh, double v, double EJ)
{
    double **t = NULL;
    int i,j, k, le;
    t = (double **) calloc(4, sizeof(double *));
    for(i = 0; i < 4; ++i)
  {
    t[i] = (double *) calloc(4, sizeof(double));
  }
    //printf("1here\n");
    fill_t(t,hh, EJ);
    for(i = 0; i < 20; ++i)
    {
        //printf("i : %d\n", i);
        for(k = 0; k < 4; ++k)
            for(le = 0; le < 4; ++le) {A[2*i + k][2*i + le] = A[2*i + k][2*i + le] + t[k][le];}
        //for(k = 0; k < 4; ++k)
        //{for(le = 0; le < 4; ++le) printf("%lf ", A[2*i + k][2*i + le]); printf("\n");}
        //for(k = 0; k < 4 + 2*i; ++k)
        //{for(le = 0; le < 4 + 2*i; ++le) printf("%.0lf ", A[k][le]); printf("\n");}
        //getchar();
    }
    //for(i = 0; i < 42; ++i)
    //{for(j = 0; j < 42; ++j) printf("%lf  ", A[i][j]); printf("\n");}
    //getchar();
    A[40][40] = A[40][40] + v;
    for(j = 0; j < 42; j++) {A[0][j] = 0.0; A[j][0] = 0.0;}
    A[0][0] = 1.0;
    for(j = 0; j < 42; j++) {A[1][j] = 0.0; A[j][1] = 0.0;}
    A[1][1] = 1.0;
    for(j = 0; j < 42; j++) {A[8][j] = 0.0; A[j][8] = 0.0;}
    A[8][8] = 1.0;
    for(j = 0; j < 42; j++) {A[22][j] = 0.0; A[j][22] = 0.0;}
    A[22][22] = 1.0;
    for(j = 0; j < 42; j++) {A[32][j] = 0.0; A[j][32] = 0.0;}
    A[32][32] = 1.0;
    //for(k = 0; k < 42; ++k)
      //  {for(le = 0; le < 42; ++le) printf("%.0lf ", A[k][le]); printf("\n");}
    //getchar();
    return 0;
}
int LUsolve(double ***L, double ***U, double *B, double *X, int M)
{
  double *Y = NULL, sum;
  int i, j, k, n;
  Y = (double *) calloc(M, sizeof(double));
  Y[0] = B[0]/(*L)[0][0];
  for(k = 1; k < M; k++)
  {
    sum = 0.0;
    for(n = 0; n < k; n++) sum = sum + Y[n] * (*L)[k][n];
    
    Y[k] = (B[k] - sum)/(*L)[k][k];
  }
  X[M - 1] = Y[M - 1]/(*U)[M - 1][M - 1];
  for(i = M - 2; i >= 0; i--)
  {
    sum = 0.0;
    for(j = i + 1; j < M; j++) sum = sum + (*U)[i][j] * X[j];
    
    X[i] = (Y[i] - sum)/(*U)[i][i];
  }
  free(Y);
  return 0;
}
int LUdecompose(double ***A, double ***L, double ***U, int M)
{
  int i, j, k;
  double sum;
  for(i = 0; i < M; i++)
  {
    for(k = i; k < M; k++)
    {
      sum = 0.0;
      for(j = 0; j < i; j++) sum = sum + (*L)[i][j] * (*U)[j][k];
      (*U)[i][k] = (*A)[i][k] - sum;
    }
    for(k = i; k < M; k++)
    {
      if(i == k) (*L)[i][i] = 1.0;
      else
      {
	sum = 0.0;
	for(j = 0; j < i; j++) sum = sum + (*L)[k][j] * (*U)[j][i];
	(*L)[k][i] = ((*A)[k][i] - sum)/(*U)[i][i];
      }
    }
  }
 return 0;
}
double h(double x)
{
    if(x <= 0.0) return 0.0;
    else return 1.0;
}
double _M_o(double x)
{
    return (x*x)/2.0;
}
double _Q_o(double x)
{
    return pow(x,3)/6.0;
}
double _M_1(double x, double hh)
{
    return pow(x - 2.0*hh,2)/2.0 * h(x - 2.0*hh);
}
double _P(double x, double hh)
{
    return pow(x - 6.0*hh,3)/6.0 * h(x - 6*hh);
}
double _q_start(double x, double hh)
{
    return pow(x - 8*hh,4)/24.0 * h(x - 8*hh);
}
double _q_end(double x, double hh)
{
    return pow(x - 17*hh,4)/24.0 * h(x - 17*hh);
}
double _R_a(double x, double hh)
{
    return pow(x - 4*hh,3)/6.0 * h(x - 4*hh);
}
double _R_b(double x, double hh)
{
    return pow(x - 11*hh,3)/6.0 * h(x - 11*hh);
}
double _R_c(double x, double hh)
{
    return pow(x - 16*hh,3)/6.0 * h(x - 16*hh);
}
double w(double x, double hh, double EJ, double M_o, double Q_o, double R_a, double R_b, double R_c, double q, double P, double M_1)
{
    return (1.0/EJ)*(M_o * _M_o(x) + Q_o * _Q_o(x) - M_1 * _M_1(x,hh) + P * _P(x,hh) + q * _q_start(x,hh) - q * _q_end(x,hh) + R_a * _R_a(x,hh) + R_b * _R_b(x,hh) + R_c * _R_c(x,hh));
}
double theta_q_end(double x, double hh)
{
    return pow(x - 17*hh,3)/6.0 * h(x - 17*hh);
}
double theta_q_start(double x, double hh)
{
    return pow(x - 8*hh,3)/6.0 * h(x - 8*hh);
}
double theta_p(double x, double hh)
{
    return pow(x - 6*hh,2)/2.0 * h(x - 6*hh);
}
double theta_m(double x, double hh)
{
    return (x - 2*hh) * h(x - 2*hh);
}
double theta_a(double x, double hh)
{
    return pow(x - 4*hh,2)/2.0 * h(x - 4*hh);
}
double theta_b(double x, double hh)
{
    return pow(x - 11*hh,2)/2.0 * h(x - 11*hh);
}
double theta_c(double x, double hh)
{
    return pow(x - 16*hh,2)/2.0 * h(x - 16*hh);
}
double theta(double x, double hh, double EJ, double M_o, double Q_o, double R_a, double R_b, double R_c, double q, double P, double M_1)
{
    return (1.0/EJ)* (M_o * x  + Q_o *(x*x)/2.0 - M_1 * theta_m(x, hh) + P * theta_p(x,hh) + q * theta_q_start(x,hh) - q * theta_q_end(x,hh) + R_a * theta_a(x,hh) + R_b * theta_b(x,hh) + R_c * theta_c(x,hh));
}
double M_q_end(double x, double hh)
{
    return pow(x - 17*hh,2)/2.0 * h(x - 17*hh);
}
double M_q_start(double x, double hh)
{
    return pow(x - 8*hh,2)/2.0 * h(x - 8*hh);
}
double M(double x, double hh, double M_o, double Q_o, double R_a, double R_b, double R_c, double q, double P, double M_1)
{
    return (M_o + Q_o * x - M_1 * h(x - 2*hh) + P * (x - 6*hh) * h(x - 6*hh) + q * M_q_start(x,hh) - q * M_q_end(x,hh) +R_a * (x - 4*hh) * h(x - 4*hh) + R_b * (x - 11*hh) * h(x - 11*hh) + R_c * (x - 16*hh) * h(x - 16*hh));
}
double Q(double x, double hh, double M_o, double Q_o, double R_a, double R_b, double R_c, double q, double P, double M_1)
{
    return Q_o + P * h(x - 6*hh) + q * (x - 8*hh) * h(x - 8*hh) - q * (x - 17*hh) * h(x - 17*hh) +R_a *h(x - 4*hh) + R_b * h(x - 11*hh) + R_c * h(x - 16*hh);
}
