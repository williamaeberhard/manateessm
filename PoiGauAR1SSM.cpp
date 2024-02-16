// Poisson-Gaussian SSM, stationary AR(1) states | v0.1
// allow for NAs in obs vec, obsind id them as 0, but states defined for all t
//   Y_t|X_t indep~ Poi with E[Y_t|X_t] = exp(beta0 + z*beta + X_t), t=1,...,n
//   X_t = phi*X_{t-1} + eps_t, where eps_t iid~ N(0,sigma^2), t=2,...,n
//   X_1 ~ N(0,sigma^2/(1-phi^2)), stationary dist since assuming |phi|<1

#include <TMB.hpp>

template <class Type> 
Type sqr(Type x){
	return x*x;
}

// parameter transformation, R -> (-1,+1)
template <class Type>
Type bound11(Type x){
	return Type(2.0)/(Type(1.0)+exp(-x)) - Type(1.0);
}

template<class Type>
Type objective_function<Type>::operator() () {

	//--------------------------------------------------------------------------
	// Inputs
	//--------------------------------------------------------------------------

	// Data
	DATA_VECTOR(obs); // counts vector, dim n
	// ^ incl NA entries from padding in R, obsind sorts them
	DATA_IVECTOR(obsind); // 1 if count available, 0 if missing, dim n
	DATA_VECTOR(logoffset); // measure of effort, coef 1 in lin pred, dim n
	DATA_MATRIX(Zmat); // matrix of det covariates, dim (n x p)
	DATA_MATRIX(Zmattemp); // piecewise temp lin fx, dim (n x 4)

	// Parameters
	PARAMETER_VECTOR(beta); // det fx, incl intercept, dim p
	PARAMETER_VECTOR(betatemp); // temp piecewise lin fx, dim 2
	PARAMETER(tphi); // transformed AR(1) coef
	PARAMETER(logsigma); // sd of Gaussian proc error for AR(1) on X

	// Random effects
	PARAMETER_VECTOR(X); // unobserved states AR(1), dim n

	// Misc
	DATA_SCALAR(tempthresh); // temperature threshold piecewise lin fx junction


	//--------------------------------------------------------------------------
	// Setup, procedures and init
	//--------------------------------------------------------------------------

	int n = obs.size(); // t = 1, ..., n
	
	Type phi = bound11(tphi);
	Type sigma = exp(logsigma);

	vector<Type> betatempvec(4); // vector of all pw lin fx coeff, incl constr
	betatempvec(0) = Type(0.0);
	// ^ intercept temp <= thresh, constrained to overall intercept
	betatempvec(1) = betatemp(0); // slope temp <= thresh, free param
	betatempvec(2) = tempthresh*(betatemp(0)-betatemp(1));
	// ^ intercept temp > thresh, constrained for continuity at thresh
	betatempvec(3) = betatemp(1); // slope temp > thresh, free param

	vector<Type> zb = Zmat*beta + Zmattemp*betatempvec;
	// ^ all deterministic effects as two linear combinations, dim n

	Type nll = 0.0; // init neg loglik

	//--------------------------------------------------------------------------
	// Random effects
	//--------------------------------------------------------------------------

	Type inisd = sigma/sqrt(Type(1.0) - sqr(phi)); // AR(1) stationary sd
	nll -= dnorm(X(0), Type(0.0), inisd, true);
	// ^ stationary dist for ini states

	for (int i=1; i<n; i++){
		Type Xcondexpec = phi*X(i-1); // AR(1) cond expec
		nll -= dnorm(X(i), Xcondexpec, sigma, true); // Gaussian loglik
	}
	

	//--------------------------------------------------------------------------
	// Observation equations
	//--------------------------------------------------------------------------

	vector<Type> fitted(n); // Poisson cond expec, dim n
	for (int i=0; i<n; i++){
		fitted(i) = exp(logoffset(i) + zb(i) + X(i)); // Poi cond expec
		// ^ fitted = exp(linear predictor)*offset
		if (obsind(i)==1){ // lkhd contrib only when obs available
			nll -= dpois(obs(i), fitted(i), true);
		}
	}


	//--------------------------------------------------------------------------
	// Outputs
	//--------------------------------------------------------------------------

	REPORT(beta);
	REPORT(betatempvec);
	REPORT(phi);
	REPORT(sigma);
	REPORT(inisd);

	REPORT(zb); // detfx linear combination of covariate effects
	REPORT(X); // predicted states
	REPORT(fitted); // Poisson cond expec, fitted values on response scale

	ADREPORT(beta); // for detfx se, needed for z-test stepwise var sel
	ADREPORT(betatemp); // for detfx se, needed for z-test stepwise var sel

	return nll;
}
