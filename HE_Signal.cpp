#include "HE_Signal.h"
#include <iostream>
#include <fstream>
#include <assert.h>
#include <Windows.h>
#include <time.h>
#include "MyTimer.h"
using namespace std;

extern ZZ *baza;

HE_Signal::HE_Signal(int lambda, int radix):
	he_context(lambda, ZZ(radix))
{
	Q = 128;
}

HE_Signal::~HE_Signal()
{
	delete baza;
	baza = nullptr;
}

void HE_Signal::cripteaza_semnal(vector<int> &semnal, vector<Mat_ZZ> &semnal_criptat)const
{
	semnal_criptat.clear();
	semnal_criptat.reserve(semnal.size());
	for (int i = 0; i < semnal.size(); i++)
	{
		Mat_ZZ C = he_context.encrypt(semnal[i] + Q);
		semnal_criptat.push_back(C);
	}
}

void HE_Signal::decripteaza_semnal(vector<int> &semnal, vector<Mat_ZZ> &semnal_criptat)const
{
	semnal.clear();
	semnal.reserve(semnal_criptat.size());
	for (int i = 0; i < semnal_criptat.size(); i++)
	{
		int s_n;
		ZZ dec = he_context.decrypt(semnal_criptat[i]);
		conv(s_n, dec);
		semnal.push_back(s_n - Q);
	}
}

void HE_Signal::convolutie_semnale(vector<Mat_ZZ> &s1, vector<Mat_ZZ> &s2, vector<Mat_ZZ> &s_out)const
{
	assert(s1.size() != 0);
	assert(s2.size() != 0);

	int len = s1.size() + s2.size() - 1;
	s_out.clear();
	s_out.reserve(len);

	Mat_ZZ zero = he_context.encrypt(0);
	for (int i = 0; i < len; i++)
	{
		s_out.push_back(zero);
	}

	for (int i = 0; i < s1.size(); i++)
	{
		for (int j = 0; j < s2.size(); j++)
		{
			Mat_ZZ C = he_context.hom_mult_opt(s1[i], s2[j]);
			s_out[i + j] = he_context.omp_hom_add(C, s_out[i + j]);
		}
	}
}

void test_speed()
{
	vector<int> v(10000000);
	for (int i = 0; i < 10000000; i++)
	{
		v[i] = rand() % 3200000;
	}
	cout << endl;

	LARGE_INTEGER frequency;
	LARGE_INTEGER start;
	LARGE_INTEGER end;
	double interval;

	QueryPerformanceFrequency(&frequency);
	QueryPerformanceCounter(&start);

	cout << "max2 = " << get_max(v) << endl;

	QueryPerformanceCounter(&end);
	interval = (double)(end.QuadPart - start.QuadPart) / frequency.QuadPart;

	// cout << "Timp abordare liniara : " << interval << endl << endl;

	QueryPerformanceFrequency(&frequency);
	QueryPerformanceCounter(&start);

	cout << "max1 = " << getmax(v, 0, v.size()) << endl;

	QueryPerformanceCounter(&end);
	double interval2 = (double)(end.QuadPart - start.QuadPart) / frequency.QuadPart;

	cout << "arb/liniar = " << interval2 / interval << endl << endl;

	// cout << "Timp abordare arborescenta : " << interval << endl << endl;
}

int getmax(vector<int> &semnal, int start, int end)
{
	assert(end >= 1);
	if (end == 1)
	{
		return semnal[start];
	}

	int max1 = getmax(semnal, start, end / 2);
	int max2 = getmax(semnal, start + end / 2, end % 2 == 0 ? end / 2 : end / 2 + 1);

	return max1 > max2 ? max1 : max2;
}

int get_max(vector<int> &semnal)
{
	assert(semnal.size() != 0);
	int max = semnal[0];

	int len = semnal.size();
	for (int i = 1; i < len; i++)
	{
		if (max < semnal[i])
		{
			max = semnal[i];
		}
	}
	return max;
}

void convolutie_semnale(vector<int> &s1, vector<int> &s2, vector<int> &s_out)
{
	s_out.clear();
	int length = s1.size() + s2.size() - 1;
	vector<int> mult(length);
	for (int i = 0; i < length; i++)
	{
		mult[i] = 0;
	}

	for (int i = 0; i < s1.size(); i++)
	{
		for (int j = 0; j < s2.size(); j++)
		{
			mult[i + j] += s1[i] * s2[j];
		}
	}

	s_out = mult;
}

void citeste_semnal(const char *filename, vector<int> &semnal)
{
	ifstream in;
	in.open(filename);
	semnal.clear();
	int pondere;
	while (in >> pondere)
	{
		semnal.push_back(pondere);
	}
	in.close();
}

void scrie_semnal(const char *filename, vector<int> &semnal)
{
	ofstream out;
	out.open(filename);
	for (int i = 0; i < semnal.size(); i++)
	{
		out << (int)semnal[i] << endl;
	}
	out.close();
}

void testare_convolutie(const char *fisier_semnal, const char *fisier_filtru)
{
	MyTimer timer;
	int w_baza = (int)pow(2, 17);
	baza = new ZZ(w_baza);
	int lambda = 112;

	cout << "Securitate = " << lambda << " biti" << endl;
	cout << "Baza = " << w_baza << endl;

	vector<int> semnal;
	vector<int> filtru;
	vector<int> convolutie;

	citeste_semnal(fisier_semnal, semnal);
	citeste_semnal(fisier_filtru, filtru);

	/*srand(time(NULL));
	for (int i = 0; i < 10; i++)
	{
		semnal[i] = rand() % 256 - 128;
		filtru[i] = rand() % 256 - 128;
	}

	cout << "Semnal = [ ";
	for (int i = 0; i < semnal.size(); i++)
	{
		cout << semnal[i] << " ";
	}
	cout << " ]" << endl << "Filtru = [ ";
	for (int i = 0; i < filtru.size(); i++)
	{
		cout << filtru[i] << " ";
	}
	cout << " ] " << endl;*/

	vector<Mat_ZZ> semnal_criptat;
	vector<Mat_ZZ> filtru_criptat;
	vector<Mat_ZZ> conv_criptata;

	HE_Signal he_signal(lambda, w_baza);
	cout << "creare he_signal - terminata.\n";

	cout << "criptare semnale ...\n";
	timer.start_timer();
	he_signal.cripteaza_semnal(semnal, semnal_criptat);
	cout << "Timp criptare semnal : " << timer.stop_timer << endl;
	timer.start_timer();
	he_signal.cripteaza_semnal(filtru, filtru_criptat);
	cout << "Timp criptare filtru : " << timer.stop_timer() << endl;
	cout << "criptare semnale - terminat.\n";

	/*vector<int> dec_s1;
	vector<int> dec_s2;
	he_signal.decripteaza_semnal(dec_s1, semnal_criptat);
	he_signal.decripteaza_semnal(dec_s2, filtru_criptat);
	scrie_semnal("dec1.dat", dec_s1);
	scrie_semnal("dec2.dat", dec_s2);
	cout << "scriere semnale decriptate -terminat .\n";*/

	timer.start_timer();
	cout << "incepe convolutia semnalelor criptate\n";
	he_signal.convolutie_semnale(semnal_criptat, filtru_criptat, conv_criptata);
	cout << "timp convolutie criptata = " << timer.stop_timer() << endl;

	timer.start_timer();
	convolutie_semnale(semnal, filtru, convolutie);
	cout << "timp convolutie in clar = " << timer.stop_timer() << endl;

	vector<int> conv_dec;

	he_signal.decripteaza_semnal(conv_dec, conv_criptata);

	assert(convolutie.size() == conv_dec.size());
	for (int i = 0; i < convolutie.size(); i++)
	{
		if (convolutie[i] != conv_dec[i])
		{
			cout << "Semnale diferite\n";
			break;
		}
	}

}