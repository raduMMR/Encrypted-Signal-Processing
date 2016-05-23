#pragma once

#include <iostream>
#include <vector>
#include "seal.h"

using namespace std;
using namespace seal;

/* integer bitwise HE */
class HE_Integer
{
	int	m_t_bits;	// working length
	const EncryptionParameters& m_parms;
public:
	HE_Integer(const EncryptionParameters& parms, int t_bits);
	~HE_Integer();

	void	encryptBit(BigPoly& ct, int bit, const BigPoly& public_key);
	void	encryptIntValue(vector<BigPoly>& vct, int val, const BigPoly& public_key);

	int		decryptBit(const BigPoly &ct, const BigPoly& secret_key);
	int 	decryptIntValue(const vector<BigPoly>& vct, const BigPoly& secret_key);

	void	encryptIntVector(vector<vector<BigPoly> >& vvct, 
						const vector<int>& values, const BigPoly& public_key);
};

