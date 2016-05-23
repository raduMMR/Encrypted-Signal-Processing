#include "HE_Utils.h"
#include "HE_Integer.h"
#include <assert.h>

extern Evaluator* eval;
extern EncryptionParameters* parms;
extern BigPoly* public_key;
extern BigPoly* secret_key;
extern BigPoly* ctxt_of_1;                  // encryption of 1 (constant)
extern int		t_bits;

/*************************************************************************************/
BigPoly compute_z(int i, int j, vector<BigPoly>& ct_x, vector<BigPoly>& ct_y)
{
	assert(ct_x.size() > 0 && ct_y.size() > 0);
	if (j == 1)
	{
		// ret = ctxt_of_1 + ct_x[i] + ct_y[i];
		BigPoly ret = eval->add( (*ctxt_of_1), ct_x[i]);
		ret = eval->add(ret, ct_y[i]);
		return ret;
	}

	int l;
	l = (j % 2 == 0) ? j / 2 : j / 2 + 1;

	BigPoly ret = compute_z(i + l, j - l, ct_x, ct_y);
	BigPoly ct = compute_z(i, l, ct_x, ct_y);

	ret = eval->multiply(ret, ct);                // ret *= ct;	
	return ret;
}

/*************************************************************************************/
BigPoly compute_t(int i, int j, vector<BigPoly>& ct_x, vector<BigPoly>& ct_y)
{
	assert(ct_x.size() > 0 && ct_y.size() > 0);
	if (j == 1)
	{
		BigPoly ret(ct_x[i]);
		ret = eval->multiply(ret, ct_y[i]);  // ret *= ct_y[i];
		ret = eval->add(ret, ct_x[i]);       // ret += ct_x[i];
		return ret;
	}

	int l;
	l = (j % 2 == 0) ? j / 2 : j / 2 + 1;

	BigPoly ret = compute_t(i + l, j - l, ct_x, ct_y);
	BigPoly ct_z = compute_z(i + l, j - l, ct_x, ct_y);
	BigPoly ct_t = compute_t(i, l, ct_x, ct_y);

	ct_z = eval->multiply(ct_z, ct_t);  // ct_z *= ct_t;
	ret = eval->add(ret, ct_z);		// ret += ct_z;

	return ret;
}

/*************************************************************************************/
BigPoly compute_s(int i, int j, vector<BigPoly>& ct_x, vector<BigPoly>& ct_y)
{
	assert(ct_x.size() > 0 && ct_y.size() > 0);
	if (j == 1)
	{
		BigPoly ret(ct_x[i]);
		ret = eval->multiply(ret, ct_y[i]);		// ret *= ct_y[i];
		ret = eval->add(ret, ct_y[i]);			// ret += ct_y[i];
		ret = eval->add(ret, (*ctxt_of_1));		// ret += *BigPoly_of_1;
		return ret;
	}

	int l;
	l = (j % 2 == 0) ? j / 2 : j / 2 + 1;

	BigPoly ret = compute_t(i + l, j - l, ct_x, ct_y);
	BigPoly ct_z = compute_z(i + l, j - l, ct_x, ct_y);
	BigPoly ct_s = compute_s(i, l, ct_x, ct_y);

	ct_z = eval->multiply(ct_z, ct_s);	// ct_z *= ct_s;
	ret = eval->add(ret, ct_z);		    // ret += ct_z;

	return ret;
}

/*************************************************************************************/
void evaluate_X_gt_Y(vector<BigPoly>& ct_x, vector<BigPoly>& ct_y, int t_bits)
{
	HE_Integer heInt(*parms, t_bits);

	cout << endl << "Evaluating (X>Y)...." << std::flush;
	cout << "Evaluating " << heInt.decryptIntValue(ct_x, *secret_key) << " > " << 
		heInt.decryptIntValue(ct_y, *secret_key) << " compute_t: " << std::flush;
	// t1 = clock();
	BigPoly ct_t = compute_t(0, t_bits, ct_x, ct_y);
	// t2 = clock();

	int dec_t = heInt.decryptBit(ct_t, *secret_key);
	// cout << endl << dec_t << " (" << clock_diff(t1, t2) << "ms)" << endl << std::flush;

}

/*************************************************************************************/
/*void evaluate_X_ge_Y(vector<BigPoly>& ct_x, vector<BigPoly>& ct_y, int t_bits)
{
	HE_Integer heInt(*context, t_bits);

	cout << endl << "Evaluating (X>=Y)..." << std::flush;
	t1 = clock();
	BigPoly ct_s = compute_s(0, t_bits, ct_x, ct_y);
	t2 = clock();
	cout << endl << "after compute_s: findBaseLevel(ct_s): " << ct_s.findBaseLevel() << std::flush;
	cout << endl << "after compute_s: findBaseLevel(ct_x): " << ct_x[0].findBaseLevel() << std::flush;
	cout << endl << "after compute_s: findBaseLevel(ct_y): " << ct_y[0].findBaseLevel() << std::flush;


	int dec_s = heInt.decryptBit(ct_s, *secretKey);
	cout << endl << dec_s << " (" << clock_diff(t1, t2) << "ms)" << endl << std::flush;
}*/

/*************************************************************************************/
/*void evaluate_X_eq_Y(vector<BigPoly>& ct_x, vector<BigPoly>& ct_y, int t_bits)
{
	HE_Integer heInt(*context, t_bits);

	cout << endl << "Evaluating (X=Y)...." << std::flush;
	
	BigPoly ct_z = compute_z(0, t_bits, ct_x, ct_y);
	

	cout << endl << "after compute_z: findBaseLevel(ct_z): " << ct_z.findBaseLevel() << std::flush;
	cout << endl << "after compute_z: findBaseLevel(ct_x): " << ct_x[0].findBaseLevel() << std::flush;
	cout << endl << "after compute_z: findBaseLevel(ct_y): " << ct_y[0].findBaseLevel() << std::flush;

	int dec_z = heInt.decryptBit(ct_z, *secretKey);
	cout << endl << dec_z << " (" << clock_diff(t1, t2) << "ms)" << endl << std::flush;
}*/

/*************************************************************************************/
vector<BigPoly> select(BigPoly& c, vector<BigPoly>& a, vector<BigPoly>& b)
{
	vector<BigPoly> ret;

	vector<BigPoly> vt1, vt2;
	for (int i = 0; i < a.size(); i++)
	{
		vt1.push_back(BigPoly(c));
		vt2.push_back(BigPoly(c));
	}

	for (int i = 0; i < a.size(); i++)
	{
		vt1[i] = eval->multiply(vt1[i], a[i]);     // vt1[i] *= a[i];
		vt2[i] = eval->add(vt2[i], *ctxt_of_1);    // vt2[i] += ENC(1);

		vt2[i] = eval->multiply(vt2[i], b[i]);     // vt2[i] *= b[i];
		vt1[i] = eval->add(vt1[i], vt2[i]);        // vt1[i] += vt2[i];

		ret.push_back(vt1[i]);
	}

	// 	if (ret[0].findBaseLevel() < 2) 
	// 	{
	// 		cout << endl << "select: recryption need =>findBaseLevel: " << ret[0].findBaseLevel() << std::flush;
	// 		batchRecrypt(ret);
	// 	}

	return ret;
}

/*************************************************************************************/
/*vector<BigPoly> getmax(vector<vector<BigPoly> >& vvct)
{
	HE_Integer heInt(*context, t_bits);
	vector<BigPoly>  ct_max = vvct[0];

	for (int i = 1; i < vvct.size(); i++)
	{
		cout << endl << endl << "Step# " << i << "...";
		cout << endl << "Evaluation of " << heInt.decryptIntValue(ct_max, *secretKey) << " > " << heInt.decryptIntValue(vvct[i], *secretKey) << std::flush;
		cout << endl << "getmax (before enter to round) =>findBaseLevel: " << ct_max[0].findBaseLevel() << std::flush;

		t1 = clock();
		if (ct_max[0].findBaseLevel() < 5)
		{
			cout << endl << "before compute_t: recryption need =>findBaseLevel: " << ct_max[0].findBaseLevel() << std::flush;
			batchRecrypt(ct_max);
			cout << endl << "before compute_t: recryption done =>findBaseLevel: " << ct_max[0].findBaseLevel() << std::flush;
		}

		BigPoly ct_t = compute_t(0, t_bits, ct_max, vvct[i]);
		t2 = clock();
		cout << endl << "compute_t: " << heInt.decryptBit(ct_t, *secretKey) << " (" << clock_diff(t1, t2) << " ms)"
			<< "; (baseLevel of ct_t:" << ct_t.findBaseLevel() << ")" << std::flush;

		if (ct_t.findBaseLevel() < 2)
		{
			cout << endl << "after compute_t: recryption need =>findBaseLevel: " << ct_t.findBaseLevel() << std::flush;
			BigPolyRecrypt(ct_t);
			cout << endl << "after compute_t: recryption done =>findBaseLevel: " << ct_t.findBaseLevel() << std::flush;
		}

		t1 = clock();
		ct_max = select(ct_t, ct_max, vvct[i]);
		t2 = clock();
		cout << endl << "select_ct_max: " << heInt.decryptIntValue(ct_max, *secretKey) << " (" << clock_diff(t1, t2) << " ms)"
			<< "; (baseLevel of ct_max:" << ct_max[0].findBaseLevel() << ")" << std::flush;
	}

	return ct_max;
}*/


/*************************************************************************************/

vector<BigPoly> getmax(vector<vector<BigPoly> >& vvct, int start, int n)
{
	assert(n >= 1);

	HE_Integer heInt(*parms, t_bits);
	if (n == 1) return vvct[start];

	vector<BigPoly> ct_max_1 = getmax(vvct, start, n / 2);
	vector<BigPoly> ct_max_2 = getmax(vvct, start + n / 2, n % 2 == 0 ? n / 2 : n / 2 + 1);

	cout << "\n\nEvaluation of " << heInt.decryptIntValue(ct_max_1, *secret_key) << " > ";
	cout << heInt.decryptIntValue(ct_max_2, *secret_key) << endl;

	// FHE_NTIMER_START(getmax);
	// 	if (ct_max_1[0].findBaseLevel() < 6)
	// 	{
	// 		cout << endl << "recryption need =>findBaseLevel(ct_max_1): " << ct_max_1[0].findBaseLevel() << std::flush;
	// 		batchRecrypt(ct_max_1);
	// 		cout << endl << "recryption done =>findBaseLevel (ct_max_1): " << ct_max_1[0].findBaseLevel() << std::flush;			
	// 	}
	// 	if (ct_max_2[0].findBaseLevel() < 6)
	// 	{
	// 		cout << endl << "recryption need =>findBaseLevel(ct_max_2): " << ct_max_2[0].findBaseLevel() << std::flush;
	// 		batchRecrypt(ct_max_2);
	// 		cout << endl << "recryption done =>findBaseLevel (ct_max_2): " << ct_max_2[0].findBaseLevel() << std::flush;			
	// 	}

	// FHE_NTIMER_START(compute_t);
	// t1 = clock();
	BigPoly ct_t = compute_t(0, t_bits, ct_max_1, ct_max_2);
	// t2 = clock();
	// FHE_NTIMER_STOP(compute_t);
	// cout << endl << "ct_t: " << heInt.decryptBit(ct_t, *secretKey) << " (" << clock_diff(t1, t2) << " sec); (baseLevel ct_t:" << ct_t.findBaseLevel() << ")" << std::flush;


	// FHE_NTIMER_START(select);
	// t1 = clock();
	vector<BigPoly> ct_max = select(ct_t, ct_max_1, ct_max_2);
	// t2 = clock();
	// FHE_NTIMER_STOP(select);
	// FHE_NTIMER_STOP(getmax);
	// cout << endl << "ct_max: " << heInt.decryptIntValue(ct_max, *secret_key);; // << " (" <<
		//clock_diff(t1, t2) << " sec); (baseLevel ct_max: " << ct_max[0].findBaseLevel() << ")" << std::flush;

	return ct_max;

	//	cout << endl << "ct_t: " << heInt.decryptBit(ct_t, *secretKey) << "; (baseLevel ct_t:" << ct_t.findBaseLevel() << ")" << std::flush;		
	//	cout << "\nct_max: " << heInt.decryptIntValue(ct_max, *secretKey) << "; (baseLevel ct_max: " << ct_max[0].findBaseLevel() <<  ")" << std::flush;
	//	cout << endl << "ct_t: " << heInt.decryptBit(ct_t, *secretKey) << " (" << tm << " sec); (baseLevel ct_t:" << ct_t.findBaseLevel() << ")" << std::flush;	
	//	cout << endl << "ct_max: " << heInt.decryptIntValue(ct_max, *secretKey) << " (" << tm << " sec); (baseLevel ct_max: " << ct_max[0].findBaseLevel() <<  ")" << std::flush;
}

int gmax(vector<int>& v, int start, int n)
{
	assert(n >= 1);
	if (n == 1) return v[start];

	int max_1 = gmax(v, start, n / 2);
	int max_2 = gmax(v, start + n / 2, n % 2 == 0 ? n / 2 : n / 2 + 1);

	return max_1 > max_2 ? max_1 : max_2;
}
