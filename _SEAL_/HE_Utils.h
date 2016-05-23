#pragma once

#include <iostream>
#include <vector>
#include "seal.h"

using namespace std;
using namespace seal;

/*************************************************************************************/
/* these are used to evaluate the "comparison" */
BigPoly compute_z(int i, int j, vector<BigPoly>& ct_x, vector<BigPoly>& ct_y);
BigPoly compute_t(int i, int j, vector<BigPoly>& ct_x, vector<BigPoly>& ct_y);
BigPoly compute_s(int i, int j, vector<BigPoly>& ct_x, vector<BigPoly>& ct_y);


/*************************************************************************************/
/* comparisons evaluations */
void evaluate_X_gt_Y(vector<BigPoly>& ct_x, vector<BigPoly>& ct_y, int t_bits);
// void evaluate_X_ge_Y(vector<BigPoly>& ct_x, vector<BigPoly>& ct_y, int t_bits);
// void evaluate_X_eq_Y(vector<BigPoly>& ct_x, vector<BigPoly>& ct_y, int t_bits);

/*************************************************************************************/
/* these are used to evaluate the "maximum" */
vector<BigPoly> select(BigPoly& c, vector<BigPoly>& a, vector<BigPoly>& b);
// vector<BigPoly> getmax(vector<vector<BigPoly> >& vvct);                   // consumes a large number of levels 
vector<BigPoly> getmax(vector<vector<BigPoly> >& vvct, int start, int n); // A tree approach with a small nr. of levels consumed
int gmax(vector<int>& v, int start, int n);