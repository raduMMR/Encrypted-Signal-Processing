#pragma once
#include "OMP-DGHV\Code\Flat_DGHV.h"

class HE_Signal
{
	Flat_DGHV he_context;

	int Q;	// default 128

	/*int M;
	int R;
	int B;
	int Q;
	int w_Q;*/

public:
	HE_Signal(int lambda = 140, int radix = (int)pow(2, 15));

	void cripteaza_semnal(vector<int> &semnal, vector<Mat_ZZ> &semnal_criptat)const;

	void decripteaza_semnal(vector<int> &semnal, vector<Mat_ZZ> &semnal_criptat)const;

	void convolutie_semnale(vector<Mat_ZZ> &s1, vector<Mat_ZZ> &s2, vector<Mat_ZZ> &s_out)const;

	~HE_Signal();
};

int getmax(vector<int> &semnal, int start, int end);

void convolutie_semnale(vector<int> &s1, vector<int> &s2, vector<int> &s_out);

void citeste_semnal(const char *filename, vector<int> &semnal);

void scrie_semnal(const char *filename, vector<int> &semnal);

void test_speed();

int get_max(vector<int> &semnal);

void testare_convolutie(const char *fisier_semnal, const char *fisier_filtru);