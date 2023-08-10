/*********************************************************************
 UDF for applying the Carreau Non-newtonian fluid model
 By Rahul Ghosh, Cornell University
 05-07-2023
 **********************************************************************/
 #include "udf.h"
 #define lambda 8.668176047960703 
 #define n 1.05 
 #define mu_0 0.4062852988445957 
 #define mu_inf 0.001 
 DEFINE_PROPERTY(carreau_viscosity,c,t)
 {
 real gamma = C_STRAIN_RATE_MAG(c,t);
 real mu_carreau = mu_inf+(mu_0-mu_inf)*pow((1+pow(gamma*lambda,2)),(n-1)/2);
 if(mu_carreau<=1){
 mu_carreau = mu_carreau;}
 else{
 mu_carreau = 1;
 }
 return mu_carreau;
 }