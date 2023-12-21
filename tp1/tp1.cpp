#include <math.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_cdf.h>
#include <time.h>

#include <gsl/gsl_randist.h>
#define PI 3.14159265358979323846

/*************************************/
double price_Black(double vol, double K, double T, double rho, double St) {
    double d1, d2;

    d1 = (log(St / K) + (rho + vol * vol) * T) / (vol * sqrt(T));
    d2 = d1 - vol * sqrt(T);

    return St * gsl_cdf_ugaussian_P(d1) - K * exp(-rho * T) * gsl_cdf_ugaussian_P(d2);
}
double fplus(double X){
	
	if(X>=0) return X;
	if(X<0) return 0;
}

/*************************************/
//Variable de controle               *
/*************************************/
double calacul_Var_Controle(double K, double r, double T, int n, double S0, double Vol){
	double A, ic=0,ic1,ic2; int i;
	A =S0*exp((r-Vol*Vol/2)*T);
	
	const gsl_rng_type *L;
    gsl_rng *rng;
    gsl_rng_env_setup();
    L = gsl_rng_default;
    rng = gsl_rng_alloc(L);
    
    double s=0;
   
    
	for(i=0;i<=n;i++){
		double G = gsl_ran_gaussian(rng, 1.0);
		s= exp(-r*T)*fplus( K-A *exp(sqrt(T)*Vol*G) )/n+s;
		
	}

	return s+S0-K*exp(-r*T);
	
	
}
/*************************************/
//Variable anthetique                *
/*************************************/

double calcul_var_anth(double K, double r, double T, int n, double S0, double Vol) {
    // Declaration des variables
    double ic, ic1, ic2;
    ic=0;
    double s = 0;
    double A = S0 * exp((r - Vol * Vol / 2) * T);
    double tab[n]; // Initialize the entire array
    tab[0]=0;
    const gsl_rng_type *L;
    gsl_rng *rng;
    gsl_rng_env_setup();
    L = gsl_rng_default;
    rng = gsl_rng_alloc(L);

    // Effectuons la boucle for
    for (int i = 0; i < n; i++) {
        double G = gsl_ran_gaussian(rng, 1.0);
        s = exp(-r * T) * fplus(A * exp(sqrt(T) * Vol * G) - K) / (2 * n) + exp(-r * T) * fplus(A * exp(sqrt(T) * Vol * (-G)) - K) / (2 * n) + s;
		tab[i]=G;
        
    }
    
    // Lib?rer le generateur aleatoire GSL
    gsl_rng_free(rng);

    return s;
}
//*****************************************************************************************************
void intervalle_de_confiance(double K, double r, double T, int n, double S0, double Vol) {
    // Declaration des variables
    double ic, ic1, ic2;
    ic=0;
    double s = 0;
    double A = S0 * exp((r - Vol * Vol / 2) * T);
    double tab[n]; // Initialize the entire array
    tab[0]=0;
    const gsl_rng_type *L;
    gsl_rng *rng;
    gsl_rng_env_setup();
    L = gsl_rng_default;
    rng = gsl_rng_alloc(L);

    // Effectuons la boucle for
    for (int i = 0; i < n; i++) {
        double G = gsl_ran_gaussian(rng, 1.0);
        s = exp(-r * T) * fplus(A * exp(sqrt(T) * Vol * G) - K) / (2 * n) + exp(-r * T) * fplus(A * exp(sqrt(T) * Vol * (-G)) - K) / (2 * n) + s;
		tab[i]=G;
        
    }
    for (int j = 0; j < n; j++) {
    	
        ic = ic + pow(fplus(exp(-r * T) * (A * exp(sqrt(T) * Vol * tab[j]) - K))/2 +fplus(exp(-r * T) * (A * exp(-sqrt(T) * Vol *tab[j]) - K))/2  - s, 2) / (n - 1);
    }
    // Affichage des r?sultats
    ic1 = s - 1.96 * sqrt(ic) / sqrt(n);
    ic2 = s + 1.96 * sqrt(ic) / sqrt(n);
    
    
    printf("l'intervalle de confiance est : [%f, %f]\n",ic1, ic2);

    // Liberer le generateur aleatoire GSL
    gsl_rng_free(rng);


}
///*************************************/
int main(){

	printf("Le prix du call en utilisant la formule de black & scholes : %f\n", price_Black(0.1,100,1,0.1,111));
	printf("**********************************************************\n");
	double res=calacul_Var_Controle(100, 0.1, 1, 10000, 111, 0.1);
	printf("Le prix du call en utilisant les variables de controles %f\n", res);
	printf("**********************************************************\n");
	double x=calcul_var_anth(100, 0.1, 1, 10000, 111, 0.1);
	printf("Le prix du call en utilisant les variables anthetiques %f\n", x);
	printf("****************************L'intervalle de confiance*************************************\n");
	intervalle_de_confiance(100, 0.1, 1, 10000, 111, 0.1);
	return 0;
}




	

