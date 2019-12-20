#pragma once
#define _USE_MATH_DEFINES
#include <math.h>
#include "../../WDL/denormal.h"

class Biquad
	{
	public:
		enum BQType {
			bq_type_lowpass = 0,
			bq_type_highpass,
			bq_type_bandpass,
			bq_type_notch,
			bq_type_peak,
			bq_type_lowshelf,
			bq_type_highshelf,
			typecount
			};
		BQType mType = bq_type_peak;
	private:
		double l1, l2;
		double r1, r2;

		double a0, a1, a2;
		double b0, b1, b2;
		double csl, psl, ppsl;
		double csr, psr, ppsr;



	public:
		double mGain;
		double mFreq;
		double mQ;
		int mParam;
	public:

		Biquad()
			{
			mParam = -1;
			a0 = a1 = a2 = b0 = b1 = b2 = 0;
			csl = csr =  0;
			l1 = l2 = r1 = r2 = 0;

			CalcCoeffs(2000, 0, 1, 48000,-1);
			}
		Biquad(float f, float g, float q, float sr, BQType t = bq_type_lowpass)
			{
			a0 = a1 = a2 = b0 = b1 = b2 = 0;
			csl = csr = 0;
			l1 = l2 = r1 = r2 = 0;

			mType = t;
			CalcCoeffs(f, g, q, sr,-1);
		
			}
		// f=0-22khz
		// g=-12 to 12 db
		// q=0.33 to 12
		// sr= samplerate
		char* Types[typecount] = { "Lowpass","Hipass","Bandpass","Notch","Peak","LowShelf","HiShelf" };
		char* GetType()
			{
			return Types[mType];
			}
		void IncType()
			{
			mType = (BQType)((mType + 1) % typecount);
			Recalc();
			}
		void Recalc()
			{
			CalcCoeffs(mFreq, mGain, mQ, 48000,-1);
			//	x1 = x2 = y1 = y2 = 0;
			}
		void CalcCoeffs(double f, double g, double q, double srate,int p)
			{
			if (p != -1)
				mParam = p;
			if (f < 40) f = 40;
			mQ = q;
			if (q <= 0) q = 0.25f;
			q *= q;
			mFreq = f;
			mGain = g;
	
			double norm;
			double V = pow(10.0, fabs(mGain) / 20.0);
			double K = tan(M_PI * mFreq / srate);
			switch (mType)
				{
				case bq_type_lowpass:
					norm = 1.0 / (1 + K / q + K * K);
					b0 = K * K * norm;
					b1 = 2.0 * b0;
					b2 = b0;
					a1 = 2 * (K * K - 1) * norm;
					a2 = (1 - K / q + K * K) * norm;
					break;

				case bq_type_highpass:
					norm = 1 / (1 + K / q + K * K);
					b0 = 1 * norm;
					b1 = -2 * b0;
					b2 = b0;
					a1 = 2 * (K * K - 1) * norm;
					a2 = (1 - K / q + K * K) * norm;
					break;

				case bq_type_bandpass:
					norm = 1.0 / (1.0 + K / q + K * K);
					b0 = K / q * norm;
					b1 = 0.0;
					b2 = -b0;
					a1 = 2.0 * (K * K - 1.0) * norm;
					a2 = (1.0 - K / q + K * K) * norm;
					break;

				case bq_type_notch:
					norm = 1 / (1 + K / q + K * K);
					b0 = (1 + K * K) * norm;
					b1 = 2 * (K * K - 1) * norm;
					b2 = b0;
					a1 = b1;
					a2 = (1 - K / q + K * K) * norm;
					break;

				case bq_type_peak:
					if (mGain >= 0) {    // boost
						norm = 1.0 / (1.0 + 1.0 / q * K + K * K);
						b0 = (1.0 + V / q * K + K * K) * norm;
						b1 = 2.0 * (K * K - 1.0) * norm;
						b2 = (1.0 - V / q * K + K * K) * norm;
						a1 = b1;
						a2 = (1.0 - 1.0 / q * K + K * K) * norm;
						}
					else {    // cut
						norm = 1.0 / (1.0 + V / q * K + K * K);
						b0 = (1.0 + 1.0 / q * K + K * K) * norm;
						b1 = 2.0 * (K * K - 1) * norm;
						b2 = (1.0 - 1.0 / q * K + K * K) * norm;
						a1 = b1;
						a2 = (1.0 - V / q * K + K * K) * norm;
						}
					break;
				case bq_type_lowshelf:
					if (mGain >= 0) {    // boost
						norm = 1 / (1 + sqrt(2) * K + K * K);
						b0 = (1 + sqrt(2 * V) * K + V * K * K) * norm;
						b1 = 2 * (V * K * K - 1) * norm;
						b2 = (1 - sqrt(2 * V) * K + V * K * K) * norm;
						a1 = 2 * (K * K - 1) * norm;
						a2 = (1 - sqrt(2) * K + K * K) * norm;
						}
					else {    // cut
						norm = 1 / (1 + sqrt(2 * V) * K + V * K * K);
						b0 = (1 + sqrt(2) * K + K * K) * norm;
						b1 = 2 * (K * K - 1) * norm;
						b2 = (1 - sqrt(2) * K + K * K) * norm;
						a1 = 2 * (V * K * K - 1) * norm;
						a2 = (1 - sqrt(2 * V) * K + V * K * K) * norm;
						}
					break;
				case bq_type_highshelf:
					if (mGain >= 0) {    // boost
						norm = 1 / (1 + sqrt(2) * K + K * K);
						b0 = (V + sqrt(2 * V) * K + K * K) * norm;
						b1 = 2 * (K * K - V) * norm;
						b2 = (V - sqrt(2 * V) * K + K * K) * norm;
						a1 = 2 * (K * K - 1) * norm;
						a2 = (1 - sqrt(2) * K + K * K) * norm;
						}
					else {    // cut
						norm = 1 / (V + sqrt(2 * V) * K + K * K);
						b0 = (1 + sqrt(2) * K + K * K) * norm;
						b1 = 2 * (K * K - 1) * norm;
						b2 = (1 - sqrt(2) * K + K * K) * norm;
						a1 = 2 * (K * K - V) * norm;
						a2 = (V - sqrt(2 * V) * K + K * K) * norm;
						}
					break;
				}

			return;
			}
		double ZTransMag(double freq, double mag = 1)
			{
			double w0 = 2 * M_PI * freq;
			double cosW0 = cos(w0)*mag;
			double cos2W0 = cos(2.0 * w0)*mag;
			double	num = (b0*b0) + (b1*b1) + (b2*b2) + (2.0 * cosW0 * ((b0 * b1) + (b1 * b2))) + (2.0 * cos2W0 * b0 * b2);
			double	den = 1.0 + (a1*a1) + (a2*a2) + (2.0 * cosW0 * (a1 + (a1 * a2))) + (2.0 * cos2W0 * a2);
			if (den == 0) return 0;
			if (num == 0) return 0;

			return  (abs(num / den));
			}
	
		inline double process(double in)
			{
			csl = in;
			double out = in * b0 + l1;
			l1 = in * b1 + l2 - a1 * out;
			l2 = in * b2 - a2 * out;
			return out;
			}
		inline void process(double &inl, double &inr)
			{
			csl = inl;
			csr = inr;
			double out = inl * b0 + l1;
			l1 = inl * b1 + l2 - a1 * out;
			l2 = inl * b2 - a2 * out;
			inl = out;
		
			out = inr * b0 + r1;
			r1 = inr * b1 + r2 - a1 * out;
			r2 = inr * b2 - a2 * out;
			inr = out;

			if (WDL_DENORMAL_OR_ZERO_DOUBLE_AGGRESSIVE(&l1))
				l1 = 0.;
			if (WDL_DENORMAL_OR_ZERO_DOUBLE_AGGRESSIVE(&l2))
				l2 = 0.;
			if (WDL_DENORMAL_OR_ZERO_DOUBLE_AGGRESSIVE(&r1))
				r1 = 0.;
			if (WDL_DENORMAL_OR_ZERO_DOUBLE_AGGRESSIVE(&r2))
				r2 = 0.;

			}

	};