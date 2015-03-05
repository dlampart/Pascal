/*******************************************************************************
 * Copyright (c) 2015 IBM Corporation and others.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Eclipse Public License v1.0
 * which accompanies this distribution, and is available at
 * http://www.eclipse.org/legal/epl-v10.html
 *
 * Contributors:
 *     IBM Corporation - initial API and implementation
 *******************************************************************************/
package ch.unil.genescore.vegas;

import ch.unil.genescore.main.Settings;



/**
 * Refactored from R package CompQuadForm by Pierre Lafaye de Micheaux and Pierre Duchesne. See bottom of file for original code.
 * 
 * Title: Imhof (1961) algorithm
 * Ref. (book or article): {J. P. Imhof, Computing the Distribution of Quadratic Forms in Normal Variables, Biometrika, Volume 48, Issue 3/4 (Dec., 1961), 419-426
 * 
 * Description:
 * Distribution function (survival function in fact) of quadratic forms in normal variables using Imhof's method.
 *
 * How absolute and relative errors are defined (GSL documentation):
 * 
 * Each algorithm computes an approximation to a definite integral of the form,
 *		I = \int_a^b f(x) w(x) dx
 *	where w(x) is a weight function (for general integrands w(x)=1). The user provides absolute and relative error bounds 
 * 	(epsabs, epsrel) which specify the following accuracy requirement,
 *		|RESULT - I|  <= max(epsabs, epsrel |I|)
 *	where RESULT is the numerical approximation obtained by the algorithm. The algorithms attempt to estimate the absolute error ABSERR = |RESULT - I| in such a way that the following inequality holds,
 *		|RESULT - I| <= ABSERR <= max(epsabs, epsrel |I|)
 */
public class Imhof implements WeightedChisquareAlgorithm {

	// ARGUMENTS
	/** Distinct non-zero characteristic roots of A\Sigma */
	private double[] lambda_ = null;
	/** Absolute accuracy requested */
	private double epsabs_ = -1;
	/** Relative accuracy requested */
	private double epsrel_ = -1;
	/** Maximum number of subintervals in the partition of the given integration interval */
	private int limit_ = -1;	

	/** The result of the last call to imhof() (set by the native method in the c code) */
	private double result_ = -1;
	/** Estimate of the modulus of the absolute error, which should equal or exceed abs(i-result) (set by the native method in the c code) */
	private double estimatedErrror_ = 1; // Note, if you change the name, you have to adapt the c code, so let's stick with the errrror
	/** The status (gsl error) of the last call to imhof() (is set by the native method in the c code) */
	private int status_ = -41;
	
	
	// ============================================================================
	// NATIVE METHODS

	/** C implementation of imhof method called by JNI */
    private native void imhof(double x, double[] lambda, double epsabs, double epsrel, int limit);

    
	// ============================================================================
	// PUBLIC METHODS
	
	/** Constructor */
	public Imhof(double[] lambda) {

		lambda_ = lambda;
		
		epsabs_ = Settings.requestedAbsolutePrecision_;
		epsrel_ = Settings.requestedRelativePrecision_;	
		limit_ = Settings.gslIntegrationLimit_;
	}

	
	// ----------------------------------------------------------------------------

	/** Compute P(Q > q) */
	public double probQsupx(double q) {
		
		// Commons math did not work well
		//RombergIntegrator integrator = new RombergIntegrator(epsrel_, epsabs_, RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT, RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT);
		//IterativeLegendreGaussIntegrator integrator = new IterativeLegendreGaussIntegrator(21, epsrel_, epsabs_);
		//double result = integrator.integrate(maxEval_, new IntegrandFunction(), 0, 1);
		//numEval_ += integrator.getEvaluations();

		// Sets result_, status_, estimatedError_
		imhof(q, lambda_, epsabs_, epsrel_, limit_);
		return result_;
	}

	
	// ----------------------------------------------------------------------------
 
	/** Get the name of the given GSL error/status */
	public static String getGslErrorName(int err) {
		
		switch (err) {
		case 0: return "GSL_SUCCESS";   
		case -1: return "GSL_FAILURE";  
		case -2: return "GSL_CONTINUE";  /* iteration has not converged */
		case 1: return "GSL_EDOM";       /* input domain error, e.g sqrt(-1) */
		case 2: return "GSL_ERANGE";     /* output range error, e.g. exp(1e100) */
		case 3: return "GSL_EFAULT";     /* invalid pointer */
		case 4: return "GSL_EINVAL";     /* invalid argument supplied by user */
		case 5: return "GSL_EFAILED";    /* generic failure */
		case 6: return "GSL_EFACTOR";    /* factorization failed */
		case 7: return "GSL_ESANITY";    /* sanity check failed - shouldn't happen */
		case 8: return "GSL_ENOMEM";     /* malloc failed */
		case 9: return "GSL_EBADFUNC";   /* problem with user-supplied function */
		case 10: return "GSL_ERUNAWAY";  /* iterative process is out of control */
		case 11: return "GSL_EMAXITER";  /* exceeded max number of iterations */
		case 12: return "GSL_EZERODIV";  /* tried to divide by zero */
		case 13: return "GSL_EBADTOL";   /* user specified an invalid tolerance */
		case 14: return "GSL_ETOL";      /* failed to reach the specified tolerance */
		case 15: return "GSL_EUNDRFLW";  /* underflow */
		case 16: return "GSL_EOVRFLW";   /* overflow  */
		case 17: return "GSL_ELOSS";     /* loss of accuracy */
		case 18: return "GSL_EROUND";    /* failed because of roundoff error */
		case 19: return "GSL_EBADLEN";   /* matrix, vector lengths are not conformant */
		case 20: return "GSL_ENOTSQR";   /* matrix not square */
		case 21: return "GSL_ESING";     /* apparent singularity detected */
		case 22: return "GSL_EDIVERGE";  /* integral or series is divergent */
		case 23: return "GSL_EUNSUP";    /* requested feature is not supported by the hardware */
		case 24: return "GSL_EUNIMPL";   /* requested feature not (yet) implemented */
		case 25: return "GSL_ECACHE";    /* cache limit exceeded */
		case 26: return "GSL_ETABLE";    /* table limit exceeded */
		case 27: return "GSL_ENOPROG";   /* iteration is not making progress towards solution */
		case 28: return "GSL_ENOPROGJ";  /* jacobian evaluations are not improving the solution */
		case 29: return "GSL_ETOLF";     /* cannot reach the specified tolerance in F */
		case 30: return "GSL_ETOLX";     /* cannot reach the specified tolerance in X */
		case 31: return "GSL_ETOLG";     /* cannot reach the specified tolerance in gradient */
		case 32: return "GSL_EOF";       /* end of file */
		default: throw new RuntimeException("Unknown GSL error: " + err);
		}
	}
	

	
	// ============================================================================
	// PRIVATE FUNCTIONS
 
	// ============================================================================
	// SETTERS AND GETTERS

	public double getError() {return estimatedErrror_;}
	public double getResult() { return result_; }
	public int getIfault() { return status_; }
	 
}



			


