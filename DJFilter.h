class DJFilter
	{
	private:
		float a1L, a1R;
		float a2L, a2R;
		float a3L, a3R;
		float a4L, a4R;
		float hg, lg;
		float fc = 1;
	public:
		DJFilter() {}
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
		inline void Process(float&l, float&r)
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