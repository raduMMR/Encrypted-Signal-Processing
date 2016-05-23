#include "HE_Integer.h"

/*************************************************************************************/
HE_Integer::HE_Integer(const EncryptionParameters& params, int t_bits)
	: m_parms(params), m_t_bits(t_bits) {}

HE_Integer::~HE_Integer()
{
	//@todo
}

/*************************************************************************************/
void HE_Integer::encryptBit(BigPoly& ct, int bit, const BigPoly& public_key)
{
	Encryptor encryptor(m_parms, public_key);
	std::uint64_t value[] = { bit };
	BigPoly ptxt(1, 1, value);
	ct = encryptor.encrypt(ptxt);
}

void HE_Integer::encryptIntValue(vector<BigPoly>& vct, int val, const BigPoly& public_key)
{
	vct.clear();
	int bit = 0;
	do
	{
		bit = val % 2;
		BigPoly ct;
		encryptBit(ct, bit, public_key);
		vct.push_back(ct);
		val /= 2;
	} while (val != 0);

	for (int i = vct.size(); i < m_t_bits; i++)
	{
		BigPoly ct;
		encryptBit(ct, 0, public_key);
		vct.push_back(ct);
	}
}

int	HE_Integer::decryptBit(const BigPoly &ct, const BigPoly& secret_key)
{
	Decryptor decryptor(m_parms, secret_key);
	BigPoly decrypted = decryptor.decrypt(ct);
	int bit;
	bit = atoi(decrypted.to_string().c_str());
	return bit;
}

int HE_Integer::decryptIntValue(const vector<BigPoly>& vct, const BigPoly& secret_key)
{
	int value = 0;
	double power_of_two = 1;

	for (int i = 0; i < vct.size(); i++)
	{
		value += power_of_two*decryptBit(vct[i], secret_key);
		power_of_two *= 2;
	}

	return value;
}

void HE_Integer::encryptIntVector(vector<vector<BigPoly> >& vvct, const vector<int>& values
	, const BigPoly& public_key)
{
	vvct.clear();
	vvct.reserve(values.size());
	for (int i = 0; i < values.size(); i++)
	{
		vector<BigPoly> encrypted_integer;
		encryptIntValue(encrypted_integer, values[i], public_key);
		vvct.push_back(encrypted_integer);
	}
}


