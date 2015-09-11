/*******************************************************************************
 * Copyright (c) 2015 David Lamparter, Daniel Marbach
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *******************************************************************************/
package ch.unil.genescore.vegas;

import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;

public class DistributionMethods {
	

private static ChiSquaredDistribution chiSquared1df_ = new ChiSquaredDistribution(1);
private static NormalDistribution normalDist_ = new NormalDistribution();

	
	public static double chiSquared1dfCumulativeProbabilityUpperTail(double q){
		double p;
		if (q > 50){
			double q2 = Math.sqrt(q);
			p=normalCumulativeProbabilityUpperTailApprox(q2)*2;
		}	
		else{
			p=1-chiSquared1df_.cumulativeProbability(q);
		}	
		return(p);
	}

	public static double chiSquared1dfInverseCumulativeProbabilityUpperTail(double p){
		
		double p2=p/2;
		double q;				
		if (p2 < 1E-14){	
			double upper=normalInversionUpperTailApprox(p2);
			q=Math.pow(upper,2);
			
		}
		else{			
			q=chiSquared1df_.inverseCumulativeProbability(1-p);			
		}
		return(q);
		

	}

	
	public static double normalCumulativeProbability(double q){
		double p;
		if (q < -8){	
			double upper=normalCumulativeProbabilityUpperTailApprox(q);
			p=upper;
			
		}
		else{
			p=normalDist_.cumulativeProbability(q);			
		}
		return(p);
	}

	public static double normalCumulativeProbabilityUpperTailApprox(double q){
	
		
		q=Math.abs(q);
		double aa=-(q*q)/2-Math.log(q)-0.5*Math.log(2*Math.PI);
		return(Math.exp(aa));
	}
		
public static double normalInverseCumulativeProbability(double p){
		double q;
		if (p < 10E-14){	
			double upper=normalInversionUpperTailApprox(p);
			q=(-1)*upper;
			
		}
		else{
			q=normalDist_.inverseCumulativeProbability(p);			
		}
		return(q);
	}
	
	
	public static double normalInversionUpperTailApprox(double p){
		// approximates tail integral of normal distribution function:: only use for very low values; below 10^-14
		double lp = Math.log(p);
		double diff=1;
		double a1=1;
		double a=1;
			while(diff>0.001){
			
				
				a=Math.sqrt((-lp-Math.log(Math.sqrt(2*Math.PI))-Math.log(a1))*2);
						diff=Math.abs(a-a1);	
					a1=a;
				
			}
		return(a);
	}

}
