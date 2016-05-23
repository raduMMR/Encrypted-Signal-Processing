#include <iostream>
#include <time.h>
#include <fstream>
#include "HE_Signal.h"

#include "HE_Utils.h"
#include "HE_Integer.h"

#include "FFT.h"
#include "MyTimer.h"

Evaluator* eval;
EncryptionParameters* parms;
BigPoly* public_key;
BigPoly* secret_key;
BigPoly* ctxt_of_1;                  // encryption of 1 (constant)
int		t_bits;

void print_example_banner(string title);

void example_basics();
void example_weighted_average();
void example_parameter_selection();
void example_batching();

void greaterThan();
void example_getmax();

void test_encoding_binary_data()
{
	print_example_banner("test_encoding_binary_data ");
	EncryptionParameters parms;
	parms.poly_modulus() = "1x^1024 + 1";
	parms.coeff_modulus() = ChooserEvaluator::default_parameter_options().at(1024);
	parms.plain_modulus() = 1 << 1;
	parms.decomposition_bit_count() = 8;
	parms.noise_standard_deviation() = ChooserEvaluator::default_noise_standard_deviation();
	parms.noise_max_deviation() = ChooserEvaluator::default_noise_max_deviation();

	cout << "Encryption parameters specify " << parms.poly_modulus().significant_coeff_count() << " coefficients with "
		<< parms.coeff_modulus().significant_bit_count() << " bits per coefficient" << endl;


	// Generate keys.
	cout << "Generating keys..." << endl;
	KeyGenerator generator(parms);
	generator.generate();
	cout << "... key generation complete" << endl;
	BigPoly public_key = generator.public_key();
	BigPoly secret_key = generator.secret_key();
	EvaluationKeys evaluation_keys = generator.evaluation_keys();

	// const int values[] = { 0,0, 0,1, 1,0, 1,1 };
	std::uint64_t value_0[] = { 0 };
	std::uint64_t value_1[] = { 1 };
	std::uint64_t *values[8];
	values[0] = value_0;
	values[1] = value_0;
	values[2] = value_0;
	values[3] = value_1;
	values[4] = value_1;
	values[5] = value_0;
	values[6] = value_1;
	values[7] = value_1;

	int i;
	for (i = 0; i < 8; i += 2)
	{
		// BalancedEncoder encoder(parms.plain_modulus());
		// BigPoly encoded1(1, 0, values);
		// BigPoly encoded2 = encoder.encode(values[i + 1]);

		BigPoly encoded1(1, 1, values[i]);
		BigPoly encoded2(1, 1, values[i+1]);

		Encryptor encryptor(parms, public_key);
		BigPoly encrypted1 = encryptor.encrypt(encoded1);
		BigPoly encrypted2 = encryptor.encrypt(encoded2);
		Evaluator evaluator(parms, evaluation_keys);

		BigPoly encryptedsum = evaluator.add(encrypted1, encrypted2);
		BigPoly encryptedproduct = evaluator.multiply(encrypted1, encrypted2);

		Decryptor decryptor(parms, secret_key);
		BigPoly decrypted1 = decryptor.decrypt(encrypted1);
		BigPoly decrypted2 = decryptor.decrypt(encrypted2);

		BigPoly decryptedsum = decryptor.decrypt(encryptedsum);
		BigPoly decryptedproduct = decryptor.decrypt(encryptedproduct);


		cout << decryptedsum.to_string() << " " << decryptedproduct.to_string() << endl;
		/*int decoded1 = encoder.decode_int32(decrypted1);
		int decoded2 = encoder.decode_int32(decrypted2);
		int decodedsum = encoder.decode_int32(decryptedsum);
		int decodedproduct = encoder.decode_int32(decryptedproduct);

		if (values[i][0] != decoded1 || values[i + 1][0] != decoded2)
		{
			cout << "Error at encryption/decryption : ";
			cout << values[i] << " " << decoded1 << endl;
			cout << values[i + 1] << " " << decoded2 << endl;
			break;
		

		if (((values[i][0] + values[i + 1][0]) % 2) != decodedsum)
		{
			cout << "Error at addition : " << values[i][0] << ", " << values[i + 1][0];
			cout << ", " << decodedsum << "\n";
			// break;
		}
		if ((values[i][0] * values[i + 1][0]) != decodedproduct)
		{
			cout << "Error at multiplication : " << values[i][0] << ", " << values[i + 1][0];
			cout << ", " << decodedproduct << "\n";
			// break;
		}*/
	}

	cout << "\nIteration " << i << "\n";
}


/*void getmax_binary()
{
	print_example_banner("Example: getmax Binary");
	EncryptionParameters parms;
	parms.poly_modulus() = "1x^1024 + 1";
	parms.coeff_modulus() = ChooserEvaluator::default_parameter_options().at(1024);
	parms.plain_modulus() = 1 << 1;
	parms.decomposition_bit_count() = 32;
	parms.noise_standard_deviation() = ChooserEvaluator::default_noise_standard_deviation();
	parms.noise_max_deviation() = ChooserEvaluator::default_noise_max_deviation();

	cout << "Encryption parameters specify " << parms.poly_modulus().significant_coeff_count() << " coefficients with "
		<< parms.coeff_modulus().significant_bit_count() << " bits per coefficient" << endl;


	// Generate keys.
	cout << "Generating keys..." << endl;
	KeyGenerator generator(parms);
	generator.generate();
	cout << "... key generation complete" << endl;
	BigPoly public_key = generator.public_key();
	BigPoly secret_key = generator.secret_key();
	EvaluationKeys evaluation_keys = generator.evaluation_keys();

	const int values[] = {0,0, 0,1, 1,0, 1,1};
	for (int i = 0; i < 8; i+=2)
	{
		BalancedEncoder encoder(parms.plain_modulus());
		BigPoly encoded1 = encoder.encode(values[i]);
		BigPoly encoded2 = encoder.encode(values[i+1]);

		Encryptor encryptor(parms, public_key);
		BigPoly encrypted1 = encryptor.encrypt(encoded1);
		BigPoly encrypted2 = encryptor.encrypt(encoded2);
		Evaluator evaluator(parms, evaluation_keys);

		BigPoly encryptedsum = evaluator.add(encrypted1, encrypted2);
		BigPoly encryptedproduct = evaluator.multiply(encrypted1, encrypted2);

		Decryptor decryptor(parms, secret_key);
		BigPoly decrypted1 = decryptor.decrypt(encrypted1);
		BigPoly decrypted2 = decryptor.decrypt(encrypted2);

		BigPoly decryptedsum = decryptor.decrypt(encryptedsum);
		BigPoly decryptedproduct = decryptor.decrypt(encryptedproduct);
		
		int decoded1 = encoder.decode_int32(decrypted1);
		int decoded2 = encoder.decode_int32(decrypted2);
		int decodedsum = encoder.decode_int32(decryptedsum);
		int decodedproduct = encoder.decode_int32(decryptedproduct);

		if ( values[i] != decoded1 || values[i+1] != decoded2 )
		{
			cout << "Eroare la criptare/decriptarea\n";
			cout << values[i] << " " << decoded1 << endl;
			cout << values[i+1] << " " << decoded2 << endl;
			break;
		}

		if (((values[i] + values[i + 1]) % 2) != decodedsum)
		{
			cout << "Eroare la adunare : " << values[i] << ", " << values[i + 1];
			cout<<", "<<decodedsum<<"\n";
			// break;
		}
		if ( (values[i] * values[i + 1]) != decodedproduct )
		{
			cout << "Eroare la inmultire : " << values[i] << ", " << values[i + 1];
			cout << ", " << decodedproduct << "\n";
			// break;
		}
	}
}*/

// void AES_encrypt();
// void AES_decrypt();

void test(const EncryptionParameters& m_parms, BigPoly& ct, int bit,
	const BigPoly& public_key, const BigPoly secret_key)
{
	Encryptor encryptor(m_parms, public_key);
	std::uint64_t value[] = { bit };
	BigPoly ptxt(1, 1, value);
	BigPoly c = encryptor.encrypt(ptxt);

	Decryptor decryptor(m_parms, secret_key);
	BigPoly decrypted = decryptor.decrypt(c);
	bit = atoi(decrypted.to_string().c_str());

	cout << "bit = " << bit << endl;
}

int main2()
{
	example_batching();

	return 0;

	parms = new EncryptionParameters();
	parms->poly_modulus() = "1x^256 + 1";
	parms->coeff_modulus() = ChooserEvaluator::default_parameter_options().at(1024);
	parms->plain_modulus() = 1 << 1;
	parms->decomposition_bit_count() = 8;
	parms->noise_standard_deviation() = ChooserEvaluator::default_noise_standard_deviation();
	parms->noise_max_deviation() = ChooserEvaluator::default_noise_max_deviation();
	cout << "Encryption parameters specify " << parms->poly_modulus().significant_coeff_count() << " coefficients with "
		<< parms->coeff_modulus().significant_bit_count() << " bits per coefficient" << endl;

	cout << "Generating keys..." << endl;
	KeyGenerator generator(*parms);
	generator.generate();
	cout << "... key generation complete" << endl;
	public_key = new BigPoly(generator.public_key());
	secret_key = new BigPoly(generator.secret_key());
	EvaluationKeys evaluation_keys = generator.evaluation_keys();

	eval = new Evaluator(*parms, evaluation_keys);

	Encryptor encryptor(*parms, *public_key);
	std::uint64_t value_1[] = { 1 };
	BigPoly encoded1(1, 1, value_1);
	ctxt_of_1 = new BigPoly(encryptor.encrypt(encoded1));

	t_bits = 3;

	HE_Integer hint(*parms, t_bits);

	// vector<BigPoly> v;
	// hint.encryptIntValue(v, 2, *public_key);
	// cout << hint.decryptIntValue(v, *secret_key) << endl;
	// hint.encryptIntValue(v, 7, *public_key);
	// cout << hint.decryptIntValue(v, *secret_key) << endl;

	size_t vec_size = 2;
	vector<int> values(vec_size);
	for (int i = 0; i < vec_size; i++)
	{
		values[i] = rand() % t_bits;
	}

	vector<vector<BigPoly> > vvct;

	hint.encryptIntVector(vvct, values, (*public_key));
	
	vector<BigPoly> enc_max = getmax(vvct, 0, vvct.size());
	int max = gmax(values, 0, values.size());

	if (max == hint.decryptIntValue(enc_max, (*secret_key)) )
	{
		cout << "Evaluare corecta.\n";
	}
	else
	{
		cout << "Evaluare gresita.\n";
	}

	for (int i = 0; i < values.size(); i++)
	{
		cout << values[i] << " ";
	}
	cout << endl << "Enc_max = ";
	cout << hint.decryptIntValue(enc_max, (*secret_key)) << endl;


	// cleanup
	delete parms; parms = nullptr;
	delete eval; eval = nullptr;
	delete ctxt_of_1; ctxt_of_1 = nullptr;
	delete public_key; public_key = nullptr;
	delete secret_key; secret_key = nullptr;
    return 0;
}

void greaterThan()
{
	print_example_banner("Example: greaterThan");
	EncryptionParameters parms;

	parms.poly_modulus() = "1x^2048 + 1";

	parms.coeff_modulus() = ChooserEvaluator::default_parameter_options().at(2048);

	parms.plain_modulus() = 1 << 8; 

	parms.decomposition_bit_count() = 32;

	parms.noise_standard_deviation() = ChooserEvaluator::default_noise_standard_deviation();

	parms.noise_max_deviation() = ChooserEvaluator::default_noise_max_deviation();

	cout << "Encryption parameters specify " << parms.poly_modulus().significant_coeff_count() << " coefficients with "
		<< parms.coeff_modulus().significant_bit_count() << " bits per coefficient" << endl;

	// Encode two integers as polynomials.
	const int value1 = -2;
	const int value2 = 2;
	const int max_value = -5;
	const int max_value2 = 5;

	BalancedEncoder encoder(parms.plain_modulus());
	BigPoly encoded1 = encoder.encode(value1);
	BigPoly encoded2 = encoder.encode(value2);
	BigPoly encoded3 = encoder.encode(max_value);
	BigPoly encoded4 = encoder.encode(max_value2);
	cout << "Encoded " << value1 << " as polynomial " << encoded1.to_string() << endl;
	cout << "Encoded " << value2 << " as polynomial " << encoded2.to_string() << endl;
	cout << "Encoded " << max_value << " as polynomial " << encoded3.to_string() << endl;
	cout << "Encoded " << max_value2 << " as polynomial " << encoded4.to_string() << endl;

	return;

	// Generate keys.
	cout << "Generating keys..." << endl;
	KeyGenerator generator(parms);
	generator.generate();
	cout << "... key generation complete" << endl;
	BigPoly public_key = generator.public_key();
	BigPoly secret_key = generator.secret_key();
	EvaluationKeys evaluation_keys = generator.evaluation_keys();

	cout << "Encrypting values..." << endl;
	Encryptor encryptor(parms, public_key);
	BigPoly encrypted1 = encryptor.encrypt(encoded1);
	BigPoly encrypted2 = encryptor.encrypt(encoded2);
	BigPoly encrypted3 = encryptor.encrypt(encoded3);

	Evaluator evaluator(parms, evaluation_keys);

	BigPoly gt = evaluator.sub(encrypted1, encrypted2);
	gt = evaluator.add(gt, encrypted3);


	// BigPoly relinproduct = evaluator.relinearize(encryptedproduct);

	cout << "Decrypting results..." << endl;
	Decryptor decryptor(parms, secret_key);
	BigPoly decryptedproduct = decryptor.decrypt(gt);
	int decodedproduct = encoder.decode_int32(decryptedproduct);	

	//BigPoly decryptedrelin = decryptor.decrypt(relinproduct);
	//int decodedrelin = encoder.decode_int32(decryptedrelin);

	cout << "encrypted result of " << decodedproduct << endl;
		//<< value1 << " * " << value2 << " + " << value1 << " = " << decodedproduct << endl;
	//cout << "encrypted relinearization of " << value1 << " and " << value2 << " = " << decodedrelin << endl;
	//int max_noise_bit_count = inherent_noise_max(parms).significant_bit_count();
	printf("\a\a\a");
}

void example_getmax()
{
	print_example_banner("Example: getmax Test");

	// In this example we demonstrate using some of the basic arithmetic operations on integers.

	// Create encryption parameters.
	EncryptionParameters parms;

	/*
	First choose the polynomial modulus. This must be a power-of-2 cyclotomic polynomial,
	i.e. a polynomial of the form "1x^(power-of-2) + 1". We recommend using polynomials of
	degree at least 1024.
	*/
	parms.poly_modulus() = "1x^2048 + 1";

	/*
	Next choose the coefficient modulus. The values we recommend to be used are:

	[ degree(poly_modulus), coeff_modulus ]
	[ 1024, "FFFFFFF00001" ],
	[ 2048, "3FFFFFFFFFFFFFFFFFF00001"],
	[ 4096, "3FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFC0000001"],
	[ 8192, "7FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFE00000001"],
	[ 16384, "7FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00000000000001"].

	These can be conveniently accessed using ChooserEvaluator::default_parameter_options(),
	which returns the above list of options as an std::map, keyed by the degree of the polynomial modulus.

	The user can also relatively easily choose their custom coefficient modulus. It should be a prime number
	of the form 2^A - 2^B + 1, where A > B > degree(poly_modulus). Moreover, B should be as small as possible
	for improved efficiency in modular reduction. For security, we recommend strictly adhering to the following
	size bounds: (see Lepoint-Naehrig (2014) [https://eprint.iacr.org/2014/062])
	/------------------------------------\
	| poly_modulus | coeff_modulus bound |
	| -------------|---------------------|
	| 1x^1024 + 1  | 48 bits             |
	| 1x^2048 + 1  | 96 bits             |
	| 1x^4096 + 1  | 192 bits            |
	| 1x^8192 + 1  | 384 bits            |
	| 1x^16384 + 1 | 768 bits            |
	\------------------------------------/
	*/
	parms.coeff_modulus() = ChooserEvaluator::default_parameter_options().at(2048);

	/*
	Now we set the plaintext modulus. This can be any integer, even though here we take it to be a power of two.
	A larger plaintext modulus causes the noise to grow faster in homomorphic multiplication, and
	also lowers the maximum amount of noise in ciphertexts that the system can tolerate.
	On the other hand, a larger plaintext modulus typically allows for better homomorphic integer arithmetic,
	although this depends strongly on which encoder is used to encode integers into plaintext polynomials.
	*/
	parms.plain_modulus() = 1 << 8;

	/*
	The decomposition bit count affects the behavior of the relinearization (key switch) operation,
	which is typically performed after each homomorphic multiplication. A smaller decomposition
	bit count makes relinearization slower, but improves the noise growth behavior on multiplication.
	Conversely, a larger decomposition bit count makes homomorphic multiplication faster at the cost
	of increased noise growth.
	*/
	parms.decomposition_bit_count() = 32;

	/*
	We use a constant standard deviation for the error distribution. Using a larger standard
	deviation will result in larger noise growth, but in theory should make the system more secure.
	*/
	parms.noise_standard_deviation() = ChooserEvaluator::default_noise_standard_deviation();

	/*
	For the bound on the error distribution we can also use a constant default value
	which is in fact 5 * ChooserEvaluator::default_noise_standard_deviation()
	*/
	parms.noise_max_deviation() = ChooserEvaluator::default_noise_max_deviation();

	cout << "Encryption parameters specify " << parms.poly_modulus().significant_coeff_count() << " coefficients with "
		<< parms.coeff_modulus().significant_bit_count() << " bits per coefficient" << endl;

	// Encode two integers as polynomials.
	const int value1 = 5;
	const int value2 = -7;
	BalancedEncoder encoder(parms.plain_modulus());
	BigPoly encoded1 = encoder.encode(value1);
	BigPoly encoded2 = encoder.encode(value2);
	cout << "Encoded " << value1 << " as polynomial " << encoded1.to_string() << endl;
	cout << "Encoded " << value2 << " as polynomial " << encoded2.to_string() << endl;

	// Generate keys.
	cout << "Generating keys..." << endl;
	KeyGenerator generator(parms);
	generator.generate();
	cout << "... key generation complete" << endl;
	BigPoly public_key = generator.public_key();
	BigPoly secret_key = generator.secret_key();
	EvaluationKeys evaluation_keys = generator.evaluation_keys();
	//cout << "Public Key = " << public_key.to_string() << endl;
	//cout << "Secret Key = " << secret_key.to_string() << endl;

	// Encrypt values.
	cout << "Encrypting values..." << endl;
	Encryptor encryptor(parms, public_key);
	BigPoly encrypted1 = encryptor.encrypt(encoded1);
	BigPoly encrypted2 = encryptor.encrypt(encoded2);

	// Perform arithmetic on encrypted values.
	// cout << "Performing encrypted arithmetic..." << endl;
	Evaluator evaluator(parms, evaluation_keys);
	/*cout << "... Performing negation..." << endl;
	BigPoly encryptednegated1 = evaluator.negate(encrypted1);
	cout << "... Performing addition..." << endl;
	BigPoly encryptedsum = evaluator.add(encrypted1, encrypted2);
	cout << "... Performing subtraction..." << endl;
	BigPoly encrypteddiff = evaluator.sub(encrypted1, encrypted2);
	cout << "... Performing multiplication..." << endl;
	BigPoly encryptedproduct = evaluator.multiply(encrypted1, encrypted2);*/


	BigPoly encryptedproduct = evaluator.multiply(encrypted1, encrypted2);

	// Decrypt results.
	cout << "Decrypting results..." << endl;
	Decryptor decryptor(parms, secret_key);
	/*BigPoly decrypted1 = decryptor.decrypt(encrypted1);
	BigPoly decrypted2 = decryptor.decrypt(encrypted2);
	BigPoly decryptednegated1 = decryptor.decrypt(encryptednegated1);
	BigPoly decryptedsum = decryptor.decrypt(encryptedsum);
	BigPoly decrypteddiff = decryptor.decrypt(encrypteddiff);*/
	BigPoly decryptedproduct = decryptor.decrypt(encryptedproduct);

	// Decode results.
	/*int decoded1 = encoder.decode_int32(decrypted1);
	int decoded2 = encoder.decode_int32(decrypted2);
	int decodednegated1 = encoder.decode_int32(decryptednegated1);
	int decodedsum = encoder.decode_int32(decryptedsum);
	int decodeddiff = encoder.decode_int32(decrypteddiff);*/
	int decodedproduct = encoder.decode_int32(decryptedproduct);

	// Display results.
	/*cout << value1 << " after encryption/decryption = " << decoded1 << endl;
	cout << value2 << " after encryption/decryption = " << decoded2 << endl;
	cout << "encrypted negate of " << value1 << " = " << decodednegated1 << endl;
	cout << "encrypted addition of " << value1 << " and " << value2 << " = " << decodedsum << endl;
	cout << "encrypted subtraction of " << value1 << " and " << value2 << " = " << decodeddiff << endl;
	cout << "encrypted multiplication of " << value1 << " and " << value2 << " = " << decodedproduct << endl;*/


	// How did the noise grow in these operations?
	/*int max_noise_bit_count = inherent_noise_max(parms).significant_bit_count();
	cout << "Noise in encryption of " << value1 << ": " << inherent_noise(encrypted1, parms, secret_key).significant_bit_count()
		<< "/" << max_noise_bit_count << " bits" << endl;
	cout << "Noise in encryption of " << value2 << ": " << inherent_noise(encrypted2, parms, secret_key).significant_bit_count()
		<< "/" << max_noise_bit_count << " bits" << endl;
	cout << "Noise in the sum: " << inherent_noise(encryptedsum, parms, secret_key).significant_bit_count()
		<< "/" << max_noise_bit_count << " bits" << endl;
	cout << "Noise in the product: " << inherent_noise(encryptedproduct, parms, secret_key).significant_bit_count()
		<< "/" << max_noise_bit_count << " bits" << endl;*/
}

/******************************************************************/

void print_example_banner(string title)
{
    if (!title.empty())
    {
        size_t title_length = title.length();
        size_t banner_length = title_length + 2 + 2 * 10;
        string banner_top(banner_length, '*');
        string banner_middle = string(10, '*') + " " + title + " " + string(10, '*');

        cout << endl
            << banner_top << endl
            << banner_middle << endl
            << banner_top << endl
            << endl;
    }
}


void example_signal_processing()
{
	HE_Signal signal2(1, 1, 1);
	// cout << "\n========================================\n";
	// HE_Signal signal2;

	int N = 10;
	vector<int> s1(N);
	vector<int> s2(N);

	srand(time(NULL));
	for (int i = 0; i < N; i++)
	{
		s1[i] = rand() % 256 - 128;
		s2[i] = rand() % 256 - 128;
	}

	vector<BigPoly> enc_s1;
	vector<BigPoly> enc_s2;

	signal2.encrypt_signal(s1, enc_s1);
	signal2.encrypt_signal(s2, enc_s2);

	vector<int> dec_s1;
	vector<int> dec_s2;
	signal2.decrypt_signal(dec_s1, enc_s1);
	signal2.decrypt_signal(dec_s2, enc_s2);

	for (int i = 0; i < dec_s1.size(); i++)
	{
		if (dec_s1[i] != s1[i])
		{
			cout << "Eroare la s1\n";
			break;
		}

		if (dec_s2[i] != s2[i])
		{
			cout << "Eroare la s2\n";
			break;
		}
	}

	vector<int> s_mult;
	vector<BigPoly> enc_s_mult;

	mult_plain_signals(s1, s2, s_mult);

	cout << "s1 = [ ";
	for (int i = 0; i < s1.size(); i++)
	{
		cout << s1[i] << " ";
	}
	cout << "]"<< endl;
	cout << "s2 = [ ";
	for (int i = 0; i < s2.size(); i++)
	{
		cout << s2[i] << " ";
	}
	cout << "]" << endl << endl;
	cout << "s_mult = [ ";
	for (int i = 0; i < s_mult.size(); i++)
	{
		cout << s_mult[i] << " ";
	}
	cout << "]" << endl;

	MyTimer timer;
	timer.start_timer();
	signal2.mult_enc_signals(enc_s1, enc_s2, enc_s_mult);
	cout << "\nTimp convolutie homomorfica : " << timer.stop_timer() << endl << endl;

	vector<int> dec_s_mult;
	signal2.decrypt_no_relin(dec_s_mult, enc_s_mult);
	cout << "DECs_m = [ ";
	for (int i = 0; i < dec_s_mult.size(); i++)
	{
		if (s_mult[i] != dec_s_mult[i] )
		{
			cout << "\n!!! Valorile convolutiei in clar difera de rezultatele convolutiei homomorfice.\n";
			break;
		}
		cout << (int)dec_s_mult[i] << " ";
	}
	cout << " ]\n";

	cout << "\nFinal prelucrare semnale\n\n";

}

void example();

int main()
{
	// FFT fft;
	// fft.test_fft();

	example_signal_processing();

	// example();

	// example_basics();

	// example_batching(); // batching-ul nu suporta nr. negative

	return 0;
}

void example()
{
	EncryptionParameters parms;
	BigPoly public_key, secret_key;
	EvaluationKeys evaluation_keys;

	// parms.poly_modulus() = "1x^16384 + 1";
	// parms.coeff_modulus() = "7FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00000000000001";
	// parms.coeff_modulus() = ChooserEvaluator::default_parameter_options().at(4096);
	// parms.plain_modulus() = 1 << (int)pow(2,14);
	// parms.decomposition_bit_count() = 32;
	parms.poly_modulus() = "1x^2048 + 1";
	parms.coeff_modulus() = ChooserEvaluator::default_parameter_options().at(2048);
	parms.plain_modulus() = 1 << 8;
	parms.decomposition_bit_count() = 32;
	parms.noise_standard_deviation() = ChooserEvaluator::default_noise_standard_deviation();
	parms.noise_max_deviation() = ChooserEvaluator::default_noise_max_deviation();
	cout << "Generating keys..." << endl;
	KeyGenerator generator(parms);
	generator.generate();
	cout << "... key generation complete" << endl;
	public_key = generator.public_key();
	secret_key = generator.secret_key();
	evaluation_keys = generator.evaluation_keys();

	/*ifstream in;
	in.open("parms2048.out");
	parms.load(in);
	in.close();

	in.open("sk2048.out");
	secret_key.load(in);
	in.close();

	in.open("pk2048.out");
	public_key.load(in);
	in.close();

	in.open("ek2048.out");
	evaluation_keys.load(in);
	in.close();*/

	const int value1 = 2;
	const int value2 = 2;
	BalancedEncoder encoder(parms.plain_modulus());
	BigPoly encoded1 = encoder.encode(value1);
	BigPoly encoded2 = encoder.encode(value2);

	cout << "Codificare cu succes.\n";

	Encryptor encryptor(parms, public_key);
	BigPoly encrypted1 = encryptor.encrypt(encoded1);
	BigPoly encrypted2 = encryptor.encrypt(encoded2);

	cout << "Criptare cu succes\n";
	Evaluator evaluator(parms, evaluation_keys);
	cout << "... Performing multiplication..." << endl;

	vector<BigPoly> copies(100);

	// "Mini convolutie fara reliniarizare dupa inmultire"
	MyTimer timer;

	timer.start_timer();
	BigPoly encryptedproduct = evaluator.multiply_norelin(encrypted1, encrypted2);
	encryptedproduct = evaluator.multiply_norelin(encrypted1, encrypted2);
	encryptedproduct = evaluator.multiply_norelin(encrypted1, encrypted2);
	encryptedproduct = evaluator.multiply_norelin(encrypted1, encrypted2);
	encryptedproduct = evaluator.multiply_norelin(encrypted1, encrypted2);
	encryptedproduct = evaluator.multiply_norelin(encrypted1, encrypted2);
	cout << "*** TIMP 6 (SASE) INMULTIRI : " << timer.stop_timer() << endl;

	for (int i = 0; i < 100; i++)
	{
		copies[i] = encryptedproduct;
	}
	BigPoly conv = evaluator.add_many(copies);
	// cout << "timp conv NO RELIN = " << timer.stop_timer() << endl;

	timer.start_timer();
	encryptedproduct = evaluator.multiply(encrypted1, encrypted2);
	for (int i = 0; i < 100; i++)
	{
		copies[i] = encryptedproduct;
	}

	BigPoly conv_relin = evaluator.add_many(copies);
	cout << "timp conv CU RELIN = " << timer.stop_timer() << endl;

	try
	{
		Decryptor decryptor_relin(parms, secret_key);
		BigPoly decoded_relin = decryptor_relin.decrypt(conv_relin);
		int decodedproduct_relin = encoder.decode_int32(decoded_relin);
		cout << "RELIN = " << decodedproduct_relin << endl;
	}
	catch (std::invalid_argument)
	{
		cout << "Invalid argument exception thrown\n";
	}
	

	Decryptor decryptor(parms, secret_key, 2);
	BigPoly decoded = decryptor.decrypt(conv);

	

	int decodedproduct = encoder.decode_int32(decoded);
	
	// int decodedsum = encoder.decode_int32(decryptedsum);
	// int s = encoder.decode_int32(sum);

	cout << "FARA RELIN " << decodedproduct << endl;
	
	// cout << "encrypted add many of " << value1 << " and " << value2 << " = " << s << endl;

	/*ofstream out;
	out.open("parms2048.out");
	parms.save(out);
	out.close();

	out.open("pk2048.out");
	public_key.save(out);
	out.close();

	out.open("sk2048.out");
	secret_key.save(out);
	out.close();

	out.open("ek2048.out");
	evaluation_keys.save(out);
	out.close();*/

	cout << "Final test SEAL.\n";
}

void example_basics()
{
	print_example_banner("Example: Basics");

	// In this example we demonstrate using some of the basic arithmetic operations on integers.

	// Create encryption parameters.
	EncryptionParameters parms;

	/*
	First choose the polynomial modulus. This must be a power-of-2 cyclotomic polynomial,
	i.e. a polynomial of the form "1x^(power-of-2) + 1". We recommend using polynomials of
	degree at least 1024.
	*/
	parms.poly_modulus() = "1x^1024 + 1";

	/*
	Next choose the coefficient modulus. The values we recommend to be used are:

	[ degree(poly_modulus), coeff_modulus ]
	[ 1024, "FFFFFFF00001" ],
	[ 2048, "3FFFFFFFFFFFFFFFFFF00001"],
	[ 4096, "3FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFC0000001"],
	[ 8192, "7FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFE00000001"],
	[ 16384, "7FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00000000000001"].

	These can be conveniently accessed using ChooserEvaluator::default_parameter_options(),
	which returns the above list of options as an std::map, keyed by the degree of the polynomial modulus.

	The user can also relatively easily choose their custom coefficient modulus. It should be a prime number
	of the form 2^A - 2^B + 1, where A > B > degree(poly_modulus). Moreover, B should be as small as possible
	for improved efficiency in modular reduction. For security, we recommend strictly adhering to the following
	size bounds: (see Lepoint-Naehrig (2014) [https://eprint.iacr.org/2014/062])
	/------------------------------------\
	| poly_modulus | coeff_modulus bound |
	| -------------|---------------------|
	| 1x^1024 + 1  | 48 bits             |
	| 1x^2048 + 1  | 96 bits             |
	| 1x^4096 + 1  | 192 bits            |
	| 1x^8192 + 1  | 384 bits            |
	| 1x^16384 + 1 | 768 bits            |
	\------------------------------------/
	*/
	parms.coeff_modulus() = ChooserEvaluator::default_parameter_options().at(1024);

	/*
	Now we set the plaintext modulus. This can be any integer, even though here we take it to be a power of two.
	A larger plaintext modulus causes the noise to grow faster in homomorphic multiplication, and
	also lowers the maximum amount of noise in ciphertexts that the system can tolerate.
	On the other hand, a larger plaintext modulus typically allows for better homomorphic integer arithmetic,
	although this depends strongly on which encoder is used to encode integers into plaintext polynomials.
	*/
	parms.plain_modulus() = 1 << 8;

	/*
	The decomposition bit count affects the behavior of the relinearization (key switch) operation,
	which is typically performed after each homomorphic multiplication. A smaller decomposition
	bit count makes relinearization slower, but improves the noise growth behavior on multiplication.
	Conversely, a larger decomposition bit count makes homomorphic multiplication faster at the cost
	of increased noise growth.
	*/
	parms.decomposition_bit_count() = 32;

	/*
	We use a constant standard deviation for the error distribution. Using a larger standard
	deviation will result in larger noise growth, but in theory should make the system more secure.
	*/
	parms.noise_standard_deviation() = ChooserEvaluator::default_noise_standard_deviation();

	/*
	For the bound on the error distribution we can also use a constant default value
	which is in fact 5 * ChooserEvaluator::default_noise_standard_deviation()
	*/
	parms.noise_max_deviation() = ChooserEvaluator::default_noise_max_deviation();

	cout << "Encryption parameters specify " << parms.poly_modulus().significant_coeff_count() << " coefficients with "
		<< parms.coeff_modulus().significant_bit_count() << " bits per coefficient" << endl;

	// Encode two integers as polynomials.
	const int value1 = 127;
	const int value2 = -127;
	BalancedEncoder encoder(parms.plain_modulus());
	BigPoly encoded1 = encoder.encode(value1);
	BigPoly encoded2 = encoder.encode(value2);
	cout << "Encoded " << value1 << " as polynomial " << encoded1.to_string() << endl;
	cout << "Encoded " << value2 << " as polynomial " << encoded2.to_string() << endl;

	// Generate keys.
	cout << "Generating keys..." << endl;
	KeyGenerator generator(parms);
	generator.generate();
	cout << "... key generation complete" << endl;
	BigPoly public_key = generator.public_key();
	BigPoly secret_key = generator.secret_key();
	EvaluationKeys evaluation_keys = generator.evaluation_keys();
	//cout << "Public Key = " << public_key.to_string() << endl;
	//cout << "Secret Key = " << secret_key.to_string() << endl;

	// Encrypt values.
	cout << "Encrypting values..." << endl;
	Encryptor encryptor(parms, public_key);
	BigPoly encrypted1 = encryptor.encrypt(encoded1);
	BigPoly encrypted2 = encryptor.encrypt(encoded2);

	// Perform arithmetic on encrypted values.
	cout << "Performing encrypted arithmetic..." << endl;
	Evaluator evaluator(parms, evaluation_keys);
	cout << "... Performing negation..." << endl;
	BigPoly encryptednegated1 = evaluator.negate(encrypted1);
	cout << "... Performing addition..." << endl;
	BigPoly encryptedsum = evaluator.add(encrypted1, encrypted2);
	cout << "... Performing subtraction..." << endl;
	BigPoly encrypteddiff = evaluator.sub(encrypted1, encrypted2);
	cout << "... Performing multiplication..." << endl;
	BigPoly encryptedproduct = evaluator.multiply(encrypted1, encrypted2);

	// Decrypt results.
	cout << "Decrypting results..." << endl;
	Decryptor decryptor(parms, secret_key);
	BigPoly decrypted1 = decryptor.decrypt(encrypted1);
	BigPoly decrypted2 = decryptor.decrypt(encrypted2);
	BigPoly decryptednegated1 = decryptor.decrypt(encryptednegated1);
	BigPoly decryptedsum = decryptor.decrypt(encryptedsum);
	BigPoly decrypteddiff = decryptor.decrypt(encrypteddiff);
	BigPoly decryptedproduct = decryptor.decrypt(encryptedproduct);

	// Decode results.
	int decoded1 = encoder.decode_int32(decrypted1);
	int decoded2 = encoder.decode_int32(decrypted2);
	int decodednegated1 = encoder.decode_int32(decryptednegated1);
	int decodedsum = encoder.decode_int32(decryptedsum);
	int decodeddiff = encoder.decode_int32(decrypteddiff);
	int decodedproduct = encoder.decode_int32(decryptedproduct);
	// delete evaluation_keys;

	// Display results.
	cout << value1 << " after encryption/decryption = " << decoded1 << endl;
	cout << value2 << " after encryption/decryption = " << decoded2 << endl;
	cout << "encrypted negate of " << value1 << " = " << decodednegated1 << endl;
	cout << "encrypted addition of " << value1 << " and " << value2 << " = " << decodedsum << endl;
	cout << "encrypted subtraction of " << value1 << " and " << value2 << " = " << decodeddiff << endl;
	cout << "encrypted multiplication of " << value1 << " and " << value2 << " = " << decodedproduct << endl;

	// How did the noise grow in these operations?
	int max_noise_bit_count = inherent_noise_max(parms).significant_bit_count();
	cout << "Noise in encryption of " << value1 << ": " << inherent_noise(encrypted1, parms, secret_key).significant_bit_count()
		<< "/" << max_noise_bit_count << " bits" << endl;
	cout << "Noise in encryption of " << value2 << ": " << inherent_noise(encrypted2, parms, secret_key).significant_bit_count()
		<< "/" << max_noise_bit_count << " bits" << endl;
	cout << "Noise in the sum: " << inherent_noise(encryptedsum, parms, secret_key).significant_bit_count()
		<< "/" << max_noise_bit_count << " bits" << endl;
	cout << "Noise in the product: " << inherent_noise(encryptedproduct, parms, secret_key).significant_bit_count()
		<< "/" << max_noise_bit_count << " bits" << endl;
}

void example_weighted_average()
{
    print_example_banner("Example: Weighted Average");

    // In this example we demonstrate computing a weighted average of 10 rational numbers.
    
    // The 10 rational numbers we use are:
    const vector<double> rational_numbers { 3.1, 4.159, 2.65, 3.5897, 9.3, 2.3, 8.46, 2.64, 3.383, 2.795 };

    // The 10 weights are:
    const vector<double> coefficients { 0.1, 0.05, 0.05, 0.2, 0.05, 0.3, 0.1, 0.025, 0.075, 0.05 };

    // Create encryption parameters
    EncryptionParameters parms;

    parms.poly_modulus() = "1x^1024 + 1";
    parms.coeff_modulus() = ChooserEvaluator::default_parameter_options().at(1024);
    parms.plain_modulus() = 1 << 8;

    /*
    Since we are not doing any encrypted*encrypted multiplication in this example,
    the decomposition bit count has no practical significance. We set it to the largest
    possible value to make key generation as fast as possible. However, such a large 
    decomposition bit count can not be used to perform any encrypted*encrypted multiplication.
    */
    parms.decomposition_bit_count() = parms.coeff_modulus().bit_count();

    // Set to standard values
    parms.noise_standard_deviation() = ChooserEvaluator::default_noise_standard_deviation();
    parms.noise_max_deviation() = ChooserEvaluator::default_noise_max_deviation();

    cout << "Encryption parameters specify " << parms.poly_modulus().significant_coeff_count() << " coefficients with "
        << parms.coeff_modulus().significant_bit_count() << " bits per coefficient" << endl;

    // Generate keys.
    cout << "Generating keys..." << endl;
    KeyGenerator generator(parms);
    generator.generate();
    cout << "... key generation complete" << endl;
    BigPoly public_key = generator.public_key();
    BigPoly secret_key = generator.secret_key();
    EvaluationKeys evaluation_keys = generator.evaluation_keys();

    /*
    We will need a fractional encoder for dealing with the rational numbers.
    Here we reserve 128 coefficients of the polynomial for the integral part (low-degree terms)
    and 64 coefficients for the fractional part (high-degree terms).
    */
    BalancedFractionalEncoder encoder(parms.plain_modulus(), parms.poly_modulus(), 128, 64);

    // Create the rest of the tools
    Encryptor encryptor(parms, public_key);
    Evaluator evaluator(parms, evaluation_keys);
    Decryptor decryptor(parms, secret_key);

    // First we encrypt the rational numbers
    cout << "Encrypting ... ";
    vector<BigPoly> encrypted_rationals;
    for (int i = 0; i < 10; ++i)
    {
        BigPoly encoded_number = encoder.encode(rational_numbers[i]);
        encrypted_rationals.push_back(encryptor.encrypt(encoded_number));
        cout << to_string(rational_numbers[i]).substr(0,6) + ((i < 9) ? ", " : ".\n");
    }

    // Next we encode the coefficients. There is no reason to encrypt these since they are not private data.
    cout << "Encoding ... ";
    vector<BigPoly> encoded_coefficients;
    for (int i = 0; i < 10; ++i)
    {
        encoded_coefficients.push_back(encoder.encode(coefficients[i]));
        cout << to_string(coefficients[i]).substr(0,6) + ((i < 9) ? ", " : ".\n");
    }
    
    // We also need to encode 0.1. We will multiply the result by this to perform division by 10.
    BigPoly div_by_ten = encoder.encode(0.1);

    // Now compute all the products of the encrypted rational numbers with the plaintext coefficients
    cout << "Computing products ... ";
    vector<BigPoly> encrypted_products;
    for (int i = 0; i < 10; ++i)
    {
        /*
        We use Evaluator::multiply_plain(...) instead of Evaluator::multiply(...) (which would 
        require also the coefficient to be encrypted). This has much better noise growth
        behavior than multiplying two encrypted numbers does.
        */
        BigPoly enc_plain_product = evaluator.multiply_plain(encrypted_rationals[i], encoded_coefficients[i]);
        encrypted_products.push_back(enc_plain_product);
    }
    cout << "done." << endl;

    // Now we add together these products. The most convenient way to do that is
    // to use the function Evaluator::add_many(...).
    cout << "Add up all 10 ciphertexts ... ";
    BigPoly encrypted_dot_product = evaluator.add_many(encrypted_products);
    cout << " done." << endl;

    // Finally we divide by 10 to obtain the result.
    cout << "Divide by 10 ... ";
    BigPoly encrypted_result = evaluator.multiply_plain(encrypted_dot_product, div_by_ten);
    cout << "done." << endl;

    // Decrypt
    cout << "Decrypting ... ";
    BigPoly plain_result = decryptor.decrypt(encrypted_result);
    cout << "done." << endl;

    // Print the answer
    double result = encoder.decode(plain_result);
    cout << "Weighted average: " << result << endl;

    // How much noise did we end up with?
    cout << "Noise in the result: " << inherent_noise(encrypted_result, parms, secret_key).significant_bit_count()
        << "/" << inherent_noise_max(parms).significant_bit_count() << " bits" << endl;
}

void example_parameter_selection()
{
    print_example_banner("Example: Automatic Parameter Selection");

    /*
    Here we demonstrate the automatic parameter selection tool. Suppose we want to find parameters
    that are optimized in a way that allows us to evaluate the polynomial 42x^3-27x+1. We need to know
    the size of the input data, so let's assume that x is an integer with base-3 representation of length
    at most 10.
    */
    cout << "Finding optimized parameters for computing 42x^3-27x+1 ... ";

    ChooserEncoder chooser_encoder;
    ChooserEvaluator chooser_evaluator;
    
    /*
    First create a ChooserPoly representing the input data. You can think of this modeling a freshly
    encrypted cipheretext of a plaintext polynomial with length at most 10 coefficients, where the
    coefficients have absolute value at most 1.
    */
    ChooserPoly cinput(10, 1);

    // Compute the first term
    ChooserPoly ccubed_input = chooser_evaluator.exponentiate(cinput, 3);
    ChooserPoly cterm1 = chooser_evaluator.multiply_plain(ccubed_input, chooser_encoder.encode(42));

    // Compute the second term
    ChooserPoly cterm2 = chooser_evaluator.multiply_plain(cinput, chooser_encoder.encode(27));

    // Subtract the first two terms
    ChooserPoly csum12 = chooser_evaluator.sub(cterm1, cterm2);

    // Add the constant term 1
    ChooserPoly cresult = chooser_evaluator.add_plain(csum12, chooser_encoder.encode(1));

    // To find an optimized set of parameters, we use ChooserEvaluator::select_parameters(...).
    EncryptionParameters optimal_parms;
    chooser_evaluator.select_parameters(cresult, optimal_parms);

    cout << "done." << endl;

    // Let's print these to see what was recommended
    cout << "Selected parameters:" << endl;
    cout << "{ poly_modulus: " << optimal_parms.poly_modulus().to_string() << endl;
    cout << "{ coeff_modulus: " << optimal_parms.coeff_modulus().to_string() << endl;
    cout << "{ plain_modulus: " << optimal_parms.plain_modulus().to_dec_string() << endl;
    cout << "{ decomposition_bit_count: " << optimal_parms.decomposition_bit_count() << endl;
    cout << "{ noise_standard_deviation: " << optimal_parms.noise_standard_deviation() << endl;
    cout << "{ noise_max_deviation: " << optimal_parms.noise_max_deviation() << endl;

    // Let's try to actually perform the homomorphic computation using the recommended parameters.
    // Generate keys.
    cout << "Generating keys..." << endl;
    KeyGenerator generator(optimal_parms);
    generator.generate();
    cout << "... key generation complete" << endl;
    BigPoly public_key = generator.public_key();
    BigPoly secret_key = generator.secret_key();
    EvaluationKeys evaluation_keys = generator.evaluation_keys();

    // Create the encoding/encryption tools
    BalancedEncoder encoder(optimal_parms.plain_modulus());
    Encryptor encryptor(optimal_parms, public_key);
    Evaluator evaluator(optimal_parms, evaluation_keys);
    Decryptor decryptor(optimal_parms, secret_key);

    // Now perform the computations on real encrypted data.
    int input_value = 12345;
    BigPoly plain_input = encoder.encode(input_value);
    cout << "Encoded " << input_value << " as polynomial " << plain_input.to_string() << endl;

    cout << "Encrypting ... ";
    BigPoly input = encryptor.encrypt(plain_input);
    cout << "done." << endl;

    // Compute the first term
    cout << "Computing first term ... ";
    BigPoly cubed_input = evaluator.exponentiate(input, 3);
    BigPoly term1 = evaluator.multiply_plain(cubed_input, encoder.encode(42));
    cout << "done." << endl;

    // Compute the second term
    cout << "Computing second term ... ";
    BigPoly term2 = evaluator.multiply_plain(input, encoder.encode(27));
    cout << "done." << endl;

    // Subtract the first two terms
    cout << "Subtracting first two terms ... ";
    BigPoly sum12 = evaluator.sub(term1, term2);
    cout << "done." << endl;

    // Add the constant term 1
    cout << "Adding one ... ";
    BigPoly result = evaluator.add_plain(sum12, encoder.encode(1));
    cout << "done." << endl;

    // Decrypt and decode
    cout << "Decrypting ... ";
    BigPoly plain_result = decryptor.decrypt(result);
    cout << "done." << endl;
    
    // Finally print the result
    cout << "Polynomial 42x^3-27x+1 evaluated at x=12345: " << encoder.decode_int64(plain_result) << endl;

    // How much noise did we end up with?
    cout << "Noise in the result: " << inherent_noise(result, optimal_parms, secret_key).significant_bit_count()
        << "/" << inherent_noise_max(optimal_parms).significant_bit_count() << " bits" << endl;
};

void example_batching()
{
    print_example_banner("Example: Batching using CRT");

    // Create encryption parameters
    EncryptionParameters parms;

    /*
    For PolyCRTBuilder we need to use a plain modulus congruent to 1 modulo 2*degree(poly_modulus).
    We could use the following parameters:

    parms.poly_modulus() = "1x^4096 + 1";
    parms.coeff_modulus() = ChooserEvaluator::default_parameter_options().at(4096);
    parms.plain_modulus() = 1 073 153;

    However, the primes suggested by ChooserEvaluator::default_parameter_options() are highly
    non-optimal for PolyCRTBuilder. The problem is that the noise in a freshly encrypted ciphertext
    will contain an additive term of the size (coeff_modulus % plain_modulus)*(largest coeff of plaintext).
    In the case of PolyCRTBuilder, the message polynomials typically have very large coefficients
    (of the size plain_modulus) and for a prime plain_modulus the remainder coeff_modulus % plain_modulus
    is typically also of the size of plain_modulus. Thus we get a term of size plain_modulus^2 to
    the noise of a freshly encrypted ciphertext! This is very bad, as normally the initial noise
    is close to size plain_modulus.

    Thus, for improved performance when using PolyCRTBuilder, we recommend the user to use their own
    custom coeff_modulus. The prime should be of the form 2^A - D, where D is as small as possible.
    The plain_modulus should be simultaneously chosen to be a prime so that coeff_modulus % plain_modulus == 1,
    and that it is congruent to 1 modulo 2*degree(poly_modulus). Finally, coeff_modulus should be bounded
    by the following strict upper bounds to ensure security:
    /------------------------------------\
    | poly_modulus | coeff_modulus bound |
    | -------------|---------------------|
    | 1x^1024 + 1  | 48 bits             |
    | 1x^2048 + 1  | 96 bits             |
    | 1x^4096 + 1  | 192 bits            |
    | 1x^8192 + 1  | 384 bits            |
    | 1x^16384 + 1 | 768 bits            |
    \------------------------------------/

    However, one issue with using such primes is that they are never NTT primes, i.e. not congruent 
    to 1 modulo 2*degree(poly_modulus), and hence might not allow for certain optimizations to be 
    used in polynomial arithmetic. Another issue is that the search-to-decision reduction of RLWE 
    does not apply to non-NTT primes, but this is not known to result in any concrete reduction
    in the security level.

    In this example we use the prime 2^190 - 42385533 as our coefficient modulus. The user should 
    try switching between this and ChooserEvaluator::default_parameter_options().at(4096) to see 
    the significant difference in the noise level at the end of the computation.
    */
    parms.poly_modulus() = "1x^4096 + 1";
    parms.coeff_modulus() = "3FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFD793F83";
    //parms.coeff_modulus() = ChooserEvaluator::default_parameter_options().at(4096);
    parms.plain_modulus() = 1073153;

    parms.decomposition_bit_count() = 32;
    parms.noise_standard_deviation() = ChooserEvaluator::default_noise_standard_deviation();
    parms.noise_max_deviation() = ChooserEvaluator::default_noise_max_deviation();

    cout << "Encryption parameters specify " << parms.poly_modulus().significant_coeff_count() << " coefficients with "
        << parms.coeff_modulus().significant_bit_count() << " bits per coefficient" << endl;

    // Create the PolyCRTBuilder
    PolyCRTBuilder crtbuilder(parms.plain_modulus(), parms.poly_modulus());
    size_t slot_count = crtbuilder.get_slot_count();

    // Create a vector of values that are to be stored in the slots. We initialize all values to 0 at this point.
    vector<BigUInt> values(slot_count, BigUInt(parms.plain_modulus().bit_count(), static_cast<uint64_t>(0)));

    // Set the first few entries of the values vector to be non-zero
    values[0] = 122;
    values[1] = 255;
    values[2] = 125;
    values[3] = 7;
    values[4] = 255;
    values[5] = 13;

    // Now compose these into one polynomial using PolyCRTBuilder
    cout << "Plaintext slot contents (slot, value): ";
    for (size_t i = 0; i < 6; ++i)
    {
        string to_write = "(" + to_string(i) + ", " + values[i].to_dec_string() + ")";
        to_write += (i != 5) ? ", " : "\n";
        cout << to_write;
    }
    BigPoly plain_composed_poly = crtbuilder.compose(values);

    // Let's do some homomorphic operations now. First we need all the encryption tools.
    // Generate keys.
    cout << "Generating keys..." << endl;
    KeyGenerator generator(parms);
    generator.generate();
    cout << "... key generation complete" << endl;
    BigPoly public_key = generator.public_key();
    BigPoly secret_key = generator.secret_key();
    EvaluationKeys evaluation_keys = generator.evaluation_keys();

    // Create the encryption tools
    Encryptor encryptor(parms, public_key);
    Evaluator evaluator(parms, evaluation_keys);
    Decryptor decryptor(parms, secret_key);

	Decryptor decryptor_no_relin(parms, secret_key, 2);

    // Encrypt plain_composed_poly
    cout << "Encrypting ... ";
    BigPoly encrypted_composed_poly = encryptor.encrypt(plain_composed_poly);
    cout << "done." << endl;

    // Let's square the encrypted_composed_poly
    cout << "Squaring the encrypted polynomial ... ";
	MyTimer timer;
	timer.start_timer();
    // BigPoly encrypted_square = evaluator.exponentiate(encrypted_composed_poly, 2);
	BigPoly encrypted_square = evaluator.multiply_norelin(encrypted_composed_poly, encrypted_composed_poly);
	cout << "Timp exponentiere la patrat = " << timer.stop_timer() << endl;
    cout << "done." << endl;
    cout << "Decrypting the squared polynomial ... ";
    BigPoly plain_square = decryptor_no_relin.decrypt(encrypted_square);
    cout << "done." << endl;
    
    // Print the squared slots
    cout << "Squared slot contents (slot, value): ";
    for (size_t i = 0; i < 6; ++i)
    {
        string to_write = "(" + to_string(i) + ", " + crtbuilder.get_slot(plain_square, i).to_dec_string() + ")";
        to_write += (i != 5) ? ", " : "\n";
        cout << to_write;
    }

    // Now let's try to multiply the squares with the plaintext coefficients (3, 1, 4, 1, 5, 9, 0, 0, ..., 0).
    // First create the coefficient vector
    vector<BigUInt> plain_coeff_vector(slot_count, BigUInt(parms.plain_modulus().bit_count(), static_cast<uint64_t>(0)));
    plain_coeff_vector[0] = 3;
    plain_coeff_vector[1] = 1;
    plain_coeff_vector[2] = 4;
    plain_coeff_vector[3] = 1;
    plain_coeff_vector[4] = 5;
    plain_coeff_vector[5] = 9;

    // Use PolyCRTBuilder to compose plain_coeff_vector into a polynomial
    BigPoly plain_coeff_poly = crtbuilder.compose(plain_coeff_vector);

    // Print the coefficient vector
    cout << "Coefficient slot contents (slot, value): ";
    for (size_t i = 0; i < 6; ++i)
    {
        string to_write = "(" + to_string(i) + ", " + crtbuilder.get_slot(plain_coeff_poly, i).to_dec_string() + ")";
        to_write += (i != 5) ? ", " : "\n";
        cout << to_write;
    }

    // Now use multiply_plain to multiply each encrypted slot with the corresponding coefficient
    cout << "Multiplying squared slots with the coefficients ... ";
    BigPoly encrypted_scaled_square = evaluator.multiply_plain(encrypted_square, plain_coeff_poly);
    cout << " done." << endl;
    
    cout << "Decrypting the scaled squared polynomial ... ";
    BigPoly plain_scaled_square = decryptor.decrypt(encrypted_scaled_square);
    cout << "done." << endl;

    // Print the scaled squared slots
    cout << "Scaled squared slot contents (slot, value): ";
    for (size_t i = 0; i < 6; ++i)
    {
        string to_write = "(" + to_string(i) + ", " + crtbuilder.get_slot(plain_scaled_square, i).to_dec_string() + ")";
        to_write += (i != 5) ? ", " : "\n";
        cout << to_write;
    }

    // How much noise did we end up with?
    cout << "Noise in the result: " << inherent_noise(encrypted_scaled_square, parms, secret_key).significant_bit_count()
        << "/" << inherent_noise_max(parms).significant_bit_count() << " bits" << endl;
}

