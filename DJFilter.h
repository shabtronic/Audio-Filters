#pragma once

#ifndef FILTERTYPE
#define FILTERTYPE float
#endif

/* S.D.Smith (C) 2019 - all rights reserved
 Super simple IIR DJ Filter
 24db only

 No freq tuning
 Just simple 0-1 mapping
 0-0.5 = lowpass
 0.5-1 = highpass
 no DN checking
*/
class DJFilter
	{
	private:
		FILTERTYPE a1L, a1R, a2L, a2R, a3L, a3R, a4L, a4R;
		FILTERTYPE hg, lg;
		FILTERTYPE fc;
	public:
		DJFilter() { a1L = a1R = a2L = a2R = a3L = a3R = a4R = a4L = 0; fc = 0.5f; hg = 0; lg = 1; }
		// 0-0.5< = lowpass
		// >0.5-1 = highpass
		void SetParam(float p)
			{
			if (p <= 0.5)
				{
				lg = 1;
				hg = 0;
				fc = p * 2;
				}
			else
				{
				lg = -1;
				hg = 1;
				fc = ((p - 0.5) * 2);
				}
			fc = pow(fc, 1.0f / 2);
			}
		// 24db Lo/Hi pass - no res
		inline void Process(FILTERTYPE&l, FILTERTYPE&r)
			{
			a1L += (l - a1L)*fc;
			a1R += (r - a1R)*fc;
			a2L += (a1L - a2L)*fc;
			a2R += (a1R - a2R)*fc;
			a3L += (a2L - a3L)*fc;
			a3R += (a2R - a3R)*fc;
			a4L += (a3L - a4L)*fc;
			a4R += (a3R - a4R)*fc;
			l = l * hg + lg * a4L;
			r = r * hg + lg * a4R;
			}
	};