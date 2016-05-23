#pragma once
#include<iostream>
#include <windows.h>
#include <ctime>
#define _USE_MATH_DEFINES // for C++
#include <cmath>
#define PI 3.14159265358979323846
using namespace std;

class FFT
{
public:
	int length = 512*16; // (int)pow(2, 10);
	int powers = length / 2;

	// A class for complex numbers
private:
	class Complex
	{
		// the real and imaginary parts
	private:
		double real;
		double imaginary;

		//Constructs a complex number given the real and imaginary parts.
	public:
		Complex() :real(0), imaginary(0) {}

		Complex(double real, double imaginary)
		{
			this->real = real;
			this->imaginary = imaginary;
		}

		//  sum of two complex numbers
		Complex add(Complex c)
		{
			return Complex(real + c.real, imaginary + c.imaginary);
		}

		// difference of two complex numbers
		Complex subtract(Complex c)
		{
			return Complex(real - c.real, imaginary - c.imaginary);
		}

		//product of two complex numbers
		Complex multiply(Complex c)
		{
			return Complex(real*c.real - imaginary*c.imaginary,
				real*c.imaginary + imaginary*c.real);
		}

		Complex multiply_3(Complex c) {
			double ac = real*c.real;
			double bd = imaginary*c.imaginary;
			return Complex(ac - bd, ((real + imaginary)*(c.real + c.imaginary)) - ac - bd);
		}

		// get the real part of complex number
		double getReal()
		{
			return real;
		}

		// To display the complex number as string
		char* toString() {
			return "neimplementata\n";
		}
	};

public:

	void test_fft() {

		LARGE_INTEGER frequency;
		LARGE_INTEGER start;
		LARGE_INTEGER end;

		for (int l = 0; l < 1; l++) { // Number of trails
			Complex *x = new Complex[2 * length];
			// double *p = new double[length];  // P polynomial co- efficients
			// double *q = new double[length];    // Q polynomial co- efficients
			// double p[] = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0 }; // P polynomial co- efficients
			// double q[] = { 1.0, 2.0, 3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0 };// Q polynomial co- efficients

			double *p = generate_random_polynomial(length);
			double *q = generate_random_polynomial(length);

			int n = length * 2;
			int logn = (int)(log(n) /log(2) + 1e-10);

			int* shuffle = precomputed_rbs(n, logn);
			Complex *omega = new Complex[n];
			Complex *omega_inv = new Complex[n];
			omega = omega_computation(omega, n);
			omega_inv = omega_inverse_computation(omega_inv, n);

			QueryPerformanceFrequency(&frequency);
			QueryPerformanceCounter(&start);

			double* result = multiply((double*)p, (double*)q, n, omega, omega_inv);

			QueryPerformanceCounter(&end);
			double interval = (double)(end.QuadPart - start.QuadPart) / frequency.QuadPart;
			cout << "Timp fft: " << interval << endl;

			QueryPerformanceFrequency(&frequency);
			QueryPerformanceCounter(&start);

			double* dp_result = dp_multiply((double*)p, (double*)q, n, omega, shuffle, omega_inv);

			QueryPerformanceCounter(&end);
			double interval2 = (double)(end.QuadPart - start.QuadPart) / frequency.QuadPart;
			cout << "Timp dp_fft : " << interval2 << endl;


			QueryPerformanceFrequency(&frequency);
			QueryPerformanceCounter(&start);
			double *r = naive_poly_mult(p, q);
			QueryPerformanceCounter(&end);
			double interval3 = (double)(end.QuadPart - start.QuadPart) / frequency.QuadPart;
			cout << "simple multiplication : " << interval3 << endl;

			// cout << "dp_result = [ ";
			for (int i = 0; i < length + length - 1; i++)
			{
				if ( abs(r[i] - dp_result[i]) > 0.1 )
				{
					cout << r[i] << " ";
					cout << dp_result[i] << " " << endl;
					// break;
				}
			}
			// cout << "] " << endl;*/
			
			delete[] x;
			delete[] shuffle;
			delete[] omega;
			delete[] omega_inv;
			delete[] result;
			delete[] dp_result;
			delete[] p;
			delete[] q;
		}

	}

	double* naive_poly_mult(double* p, double*q)
	{
		double *r = new double[length * 2-1];
		for (int i = 0; i < length * 2 - 1; i++)
		{
			r[i] = 0;
		}

		for (int i = 0; i < length; i++)
		{
			for (int j = 0; j < length; j++)
			{
				r[i + j] += p[i] * q[j];
			}
		}

		return r;
	}

	double* generate_random_polynomial(int l)
	{
		double *p = new double[l];

		srand(time(NULL));
		for (int i = 0; i < l; i++)
		{
			// p[i] = rand() % 128;

			if (i % 5 == 0)
			{
				p[i] = (-1) * rand() % 128;
			}
			else
			{
				p[i] = rand() % 128;
			}
		}

		return p;
	}

	// Function to multiply two numbers.
	double* multiply(double* p, double* q, int n, Complex* omega, Complex* omega_inv)
	{
		//Generate complex objects with 2 times the size
		Complex* p_double = new Complex[n];
		Complex* q_double = new Complex[n];

		// Copy the coefficents into the real part of complex number and pad zeros to remaining terms.
		for (int i = 0; i<n / 2; i++)
			p_double[i] = Complex(p[i], 0);
		for (int i = n / 2; i<n; i++)
			p_double[i] = Complex(0, 0);
		for (int i = 0; i<n / 2; i++)
			q_double[i] = Complex(q[i], 0);
		for (int i = n / 2; i<n; i++)
			q_double[i] = Complex(0, 0);

		int pow = 1;

		// Apply the FFT to the two factors
		Complex *solp = fft(p_double, omega, n, pow);
		Complex *solq = fft(q_double, omega, n, pow);

		// Multiply the results pointwise recursive
		Complex *finalsol = new Complex[n];
		for (int i = 0; i<n; i++)
			finalsol[i] = solp[i].multiply(solq[i]);

		// Apply the FFT to the pointwise product
		Complex *poly = fft(finalsol, omega_inv, n, pow);

		//
		//   get the results by normalising the results
		//
		double *result = new double[n - 1];
		for (int i = 0; i<n - 1; i++)
			result[i] = poly[i].getReal() / n;

		/*result[0] = poly[0].getReal()/n;
		for (int i=1; i<n-1; i++)
		result[i] = poly[n-i].getReal()/n;*/

		//
		//   get the results by normalising the results.
		//
		delete[] p_double;
		delete[] q_double;
		delete[] finalsol;
		delete[] poly;

		return result;
	}

	double *dp_multiply(double p[], double q[], int n, Complex omega[], int shuffle[], Complex omega_inv[])
	{
		//Generate complex objects with 2 times the size
		Complex* p_double = new Complex[n];
		Complex* q_double = new Complex[n];

		// Copy the coefficents into the real part of complex number and pad zeros to remainaing terms.
		for (int i = 0; i<n / 2; i++)
			p_double[i] = Complex(p[i], 0);
		for (int i = n / 2; i<n; i++)
			p_double[i] = Complex(0, 0);
		for (int i = 0; i<n / 2; i++)
			q_double[i] = Complex(q[i], 0);
		for (int i = n / 2; i<n; i++)
			q_double[i] = Complex(0, 0);

		// Apply the Dynamic Programming FFT
		Complex *dp_solp = dpfft(p_double, omega, n, shuffle);
		Complex *dp_solq = dpfft(q_double, omega, n, shuffle);

		// Multiply the results pointwise for dynamic Programming fft
		Complex *dp_finalsol = new Complex[n];
		for (int i = 0; i<n; i++)
			dp_finalsol[i] = dp_solp[i].multiply_3(dp_solq[i]);

		Complex *dp_poly = dpfft(dp_finalsol, omega_inv, n, shuffle);

		//
		//   get the results by normalising the results for dynamic programming
		//


		double *dp_result = new double[n - 1];
		for (int i = 0; i<n - 1; i++)
			dp_result[i] = dp_poly[i].getReal() / n;
		/* dp_result[0] = dp_poly[0].getReal()/n;
		for (int i=1; i<n-1; i++)
		dp_result[i] = dp_poly[n-i].getReal()/n;
		*/

		delete[] p_double;
		delete[] q_double;
		delete[] dp_solp;
		delete[] dp_solq;
		delete[] dp_finalsol;
		delete[] dp_poly;

		return dp_result;
	}

	// Function to obtain fft
	Complex* fft(Complex pol[], Complex omega[], int n, int pow) {

		// the base case --
		if (n == 1)
		{
			Complex *copy = new Complex[1];
			copy[0] = pol[0];
			return copy;
		}

		// split  coefficients into two pieces - even and idd
		Complex* poly_even = new Complex[n];
		Complex* poly_odd = new Complex[n];

		// send them into two different array of objects
		for (int k = 0; k<n; k++) {
			if (k % 2 == 0)
				poly_even[k / 2] = pol[k];
			else
				poly_odd[k / 2] = pol[k];
		}

		// get the squares of omega values
		/*Complex[] xsquare = new Complex[n/2];
		for(int i=0;i<n/2;i++)
		xsquare[i]=omega[i].multiply(omega[i]);*/

		// call the even and odd objects arrays recursively
		Complex* solution_even = fft(poly_even, omega, n / 2, pow * 2);
		Complex* solution_odd = fft(poly_odd, omega, n / 2, pow * 2);

		// create space for the results
		Complex* poly_sol = new Complex[n];
		for (int i = 0; i < n; i++)
			poly_sol[i] = Complex(0, 0);

		// combine the pieces
		for (int i = 0; i<n / 2; i++)
		{
			poly_sol[i] = solution_even[i].add(omega[i*pow].multiply(solution_odd[i]));
			poly_sol[i + n / 2] = solution_even[i].subtract(omega[i*pow].multiply(solution_odd[i]));
		}

		int a = 2;

		delete[] poly_even;
		delete[] poly_odd;
		delete[] solution_even;
		delete[] solution_odd;

		// return the final solution of fft
		return poly_sol;
	}


	// Pre compute the omega values
	Complex* omega_computation(Complex x[], int length) {
		for (int i = 0; i<length; i++) {
			x[i] = Complex(cos((2 * i* PI) / length), sin((2 * i* PI) / length));
		}
		return x;
	}

	// Pre compute the omega values
	Complex* omega_inverse_computation(Complex x[], int length) {
		for (int i = 0; i<length; i++) {
			x[i] = Complex(cos((2 * i* PI) / length), -1 * sin((2 * i*PI) / length));
		}
		return x;
	}



	// Random values generation
	void random(double p[], double min, double max, int length) {
		double diff = max - min;
		for (int i = 0; i<length; i++) {
			p[i] = min + rand() * diff;
		}
	}

	// Dynamic FFT
	Complex* dpfft(Complex poly[], Complex x[], int n, int shuffle[]) {
		int logn = (int)(log(n) / log(2) + 1e-10);

		//Complex[][] sol = new Complex[logn+1][n];
		Complex **sol = new Complex*[2];
		sol[0] = new Complex[n];
		sol[1] = new Complex[n];

		for (int i = 0; i<n; i++) {
			sol[0][shuffle[i]] = poly[i];
		}
		int power = n / 2;
		int size = 2;
		int k=1;

		for (k = 1; k <= logn; k++) {
			for (int i = 0; i<n; i = i + size) {
				for (int j = 0; j<size / 2; j++) {
					Complex odd = x[j*power].multiply_3(sol[(k - 1) % 2][(i + j + size / 2)]);
					sol[k % 2][i + j] = sol[(k - 1) % 2][i + j].add(odd);
					sol[k % 2][i + j + size / 2] = sol[(k - 1) % 2][i + j].subtract(odd);
				}
			}
			power = power / 2;
			size = size * 2;
		}

		Complex *temp=new Complex[n];
		for (int i = 0; i <n; i++)
			temp[i] = sol[(k - 1) % 2][i];

		delete[] sol[0];
		delete[] sol[1];
		delete[] sol;

		return temp;
	}

	// Calculate  the precomputed RBS value.
	int* precomputed_rbs(int length, int logn) {

		int* shuffled = new int[length];
		for (int i = 0; i< length; i++)
			shuffled[i] = RBS(i, logn);
		return shuffled;
	}

	// Calculate the RBS values
	int RBS(int i, int k) {
		if (k == 0)
			return i;
		if (i % 2 == 1)
			return (int)pow(2, k - 1) + RBS(i / 2, k - 1);
		else
			return RBS(i / 2, k - 1);
	}

	double* copyOfRange(double *p, int i, int j)
	{
		double *copy = new double[j - i+1];

		for (int k = i; k < j; k++)
			copy[k] = p[k];

		return copy;

	}

	//LAST ASSIGMNMWNENT
	double* three_recursive(double p[], double q[], int length) {

		//base case
		if (length == 1) {
			double* pq = new double[1];
			pq[0] = p[0] * q[0];
			return pq;
		}
		else {
			//array declarations for pl,ph,ql,qh
			double* pq = new double[(2 * length)];
			double* pl = copyOfRange(p, 0, (length / 2));
			double* ph = copyOfRange(p, (length / 2), length);
			//System.out.println(Arrays.toString(ph));
			double* ql = copyOfRange(q, 0, (length / 2));
			double* qh = copyOfRange(q, (length / 2), length);

			//recursive calls
			double* plql = three_recursive(pl, ql, length / 2);
			double* phqh = three_recursive(ph, qh, length / 2);

			// pl+ph , ql+qh
			double* plandph = new double[length];
			double* qlandqh = new double[length];

			for (int i = 0; i<length / 2; i++) {
				plandph[i] = pl[i] + ph[i];
				qlandqh[i] = ql[i] + qh[i];
			}

			//recursive call - (pl+ph)(ql+qh)
			double* plmulqh = three_recursive(plandph, qlandqh, length / 2);

			// (pl+qh)(ql+qh) - plql - phqh
			double* plmulqh_f = new double[(length)];
			for (int k = 0; k<length - 1; k++) {
				plmulqh_f[k] = plmulqh[k] - plql[k] - phqh[k];
			}

			// Adding the solutions as plql + ((pl+ph)(ql+qh)-plql-qlqh)x^n/2 +phqhx^n

			int mid = ((length));
			int x = 1;
			while (x <= ((length * 2) - 1)) {
				if (x <= length / 2) {
					pq[x] = plql[x - 1];
				}
				else if (x >= ((length / 2) + 1) && x < length) {
					pq[x] = plql[x - 1] + plmulqh_f[x - ((length / 2) + 1)];
				}
				else if (x == mid) {
					pq[x] = plmulqh_f[x - ((length / 2) + 1)];
				}

				else if (x > mid && x <= ((3 * (length) / 2) - 1)) {
					pq[x] = phqh[x - (mid + 1)] + plmulqh_f[x - ((length / 2) + 1)];
				}
				else {
					pq[x] = phqh[x - (mid + 1)];
				}
				x++;
			}
			// System.out.println(Arrays.toString(pq));
			return copyOfRange(pq, 1, 2 * length);

		}
	}
	
};