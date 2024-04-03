/*
Auxiliary functions for integer encoding
*/

#pragma once

#include <iostream>
#include <cmath>
#include <vector>
#include <helib/helib.h>
#include <helib/Ctxt.h>
#include <helib/polyEval.h>

using namespace std;
using namespace helib;
using namespace NTL;

inline int intlog(unsigned long base, unsigned long input)
{
	int res = max(static_cast<int>(floor(log2(input)/log2(base))),0);
	if(power_long(base, res+1) == (long int)input)
		return res+1;
	return res;
}

void digit_decomp(vector<long>& decomp, unsigned long input, unsigned long base, int nslots);

// Simple evaluation sum f_i * X^i, assuming that babyStep has enough powers
void simplePolyEval(Ctxt& ret, const NTL::ZZX& poly, DynamicCtxtPowers& babyStep);

// The recursive procedure in the Paterson-Stockmeyer
// polynomial-evaluation algorithm from SIAM J. on Computing, 1973.
// This procedure assumes that poly is monic, deg(poly)=k*(2t-1)+delta
// with t=2^e, and that babyStep contains >= k+delta powers
void PatersonStockmeyer(Ctxt& ret, const NTL::ZZX& poly, long k, long t, long delta, DynamicCtxtPowers& babyStep, DynamicCtxtPowers& giantStep);

// This procedure assumes that k*(2^e +1) > deg(poly) > k*(2^e -1),
// and that babyStep contains >= k + (deg(poly) mod k) powers
void degPowerOfTwo(Ctxt& ret, const NTL::ZZX& poly, long k,
        DynamicCtxtPowers& babyStep, DynamicCtxtPowers& giantStep);

void recursivePolyEval(Ctxt& ret, const NTL::ZZX& poly, long k,
      DynamicCtxtPowers& babyStep, DynamicCtxtPowers& giantStep);


// From Geeks for Geeks
// Function for extended Euclidean Algorithm
int gcdExtended(int a, int b, int* x, int* y);
 
// Function to find modulo inverse of a
int modInverse(int A, int M);

int get_inverse(int nom, int dom, int p);

template<typename T, typename Allocator> void print_vector(const vector<T, Allocator>& vect, int num_entries = 10);

struct Params
{
  const long m, p, r, qbits;
  const std::vector<long> gens;
  const std::vector<long> ords;
  const std::vector<long> mvec;
  Params(long _m,
         long _p,
         long _r,
         long _qbits,
         const std::vector<long>& _gens = {},
         const std::vector<long>& _ords = {},
         const std::vector<long>& _mvec = {}) :
      m(_m), p(_p), r(_r), qbits(_qbits), gens(_gens), ords(_ords), mvec(_mvec)
  {}
  Params(const Params& other) :
      Params(other.m,
             other.p,
             other.r,
             other.qbits,
             other.gens,
             other.ords,
             other.mvec)
  {}
  bool operator!=(Params& other) const { return !(*this == other); }
  bool operator==(Params& other) const
  {
    return m == other.m && p == other.p && r == other.r &&
           qbits == other.qbits && gens == other.gens && ords == other.ords &&
           mvec == other.mvec;
  }
};

struct ContextAndKeys
{
  const Params params;

  helib::Context context;
  helib::SecKey secretKey;
  const helib::PubKey publicKey;
  const helib::EncryptedArray& ea;

  ContextAndKeys(const Params& _params) :
      params(_params),
      context(helib::ContextBuilder<helib::BGV>()
                  .m(params.m)
                  .p(params.p)
                  .r(params.r)
                  .bits(params.qbits)
                  .gens(params.gens)
                  .ords(params.ords)
                  .bootstrappable(!params.mvec.empty())
                  .mvec(params.mvec)
                  .build()),
      secretKey(context),
      publicKey((secretKey.GenSecKey(),
                 helib::addSome1DMatrices(secretKey),
                 secretKey)),
      ea(context.getEA())
  {
  }
};

struct Meta
{
  std::unique_ptr<ContextAndKeys> data;
  Meta& operator()(const Params& params)
  {
    data = std::make_unique<ContextAndKeys>(params);
    return *this;
  }
};