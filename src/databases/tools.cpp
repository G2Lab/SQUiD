#include "tools.hpp"


void digit_decomp(vector<long> &decomp, unsigned long input, unsigned long base, int nslots)
{
  decomp.clear();
  decomp.resize(nslots, 0);
  int power = intlog(base, input) + 1;
  if (power > nslots)
  {
    cout << "Input character is too big to be converted" << endl;
    exit(1);
  }
  unsigned long rest = input;
  unsigned long coeff;

  int i = 0;
  while (i < power)
  {
    coeff = rest % base;
    decomp[i] = coeff;
    rest = (rest - coeff) / base;
    i++;
  }
}

void int_to_slot(ZZX& poly, unsigned long input, unsigned long enc_base, long m_slotDeg)
{ 
    vector<long> decomp;

    //decomposition of a digit
    digit_decomp(decomp, input, enc_base, m_slotDeg);
    poly = ZZX(INIT_MONO, 0, 0);
    for (int iCoef = 0; iCoef < m_slotDeg; iCoef++)
    {
        poly+=ZZX(INIT_MONO, iCoef, decomp[iCoef]);
    }
}

unsigned long digit_compose(const vector<long> &decomp, unsigned long base) {
    unsigned long result = 0;
    unsigned long multiplier = 1;

    for (size_t i = 0; i < decomp.size(); i++) {
        result += decomp[i] * multiplier;
        multiplier *= base;
    }

    return result;
}

unsigned long slot_to_int(const ZZX& poly, unsigned long enc_base, long m_slotDeg) {
    vector<long> decomp(m_slotDeg, 0);

    for (int iCoef = 0; iCoef < m_slotDeg; iCoef++) {
        decomp[iCoef] = NTL::to_long(coeff(poly, iCoef));
    }

    return digit_compose(decomp, enc_base);
}



// Simple evaluation sum f_i * X^i, assuming that babyStep has enough powers
void simplePolyEval(Ctxt &ret, const NTL::ZZX &poly, DynamicCtxtPowers &babyStep)
{
  ret.clear();
  if (deg(poly) < 0)
    return; // the zero polynomial always returns zero

  // OLD: assert(deg(poly)<=babyStep.size()); // ensure that we have enough powers
  helib::assertTrue(deg(poly) <= babyStep.size(), "BabyStep has not enough powers (required more than deg(poly))");

  NTL::ZZ coef;
  NTL::ZZ p = NTL::to_ZZ(babyStep[0].getPtxtSpace());
  for (long i = 1; i <= deg(poly); i++)
  {
    rem(coef, coeff(poly, i), p);
    if (coef > p / 2)
      coef -= p;

    Ctxt tmp = babyStep.getPower(i); // X^i
    tmp.multByConstant(coef);        // f_i X^i
    ret += tmp;
  }
  // Add the free term
  rem(coef, ConstTerm(poly), p);
  if (coef > p / 2)
    coef -= p;
  ret.addConstant(coef);
  //  if (verbose) checkPolyEval(ret, babyStep[0], poly);
}

// The recursive procedure in the Paterson-Stockmeyer
// polynomial-evaluation algorithm from SIAM J. on Computing, 1973.
// This procedure assumes that poly is monic, deg(poly)=k*(2t-1)+delta
// with t=2^e, and that babyStep contains >= k+delta powers
void PatersonStockmeyer(Ctxt &ret, const NTL::ZZX &poly, long k, long t, long delta,
                        DynamicCtxtPowers &babyStep, DynamicCtxtPowers &giantStep)
{
  if (deg(poly) <= babyStep.size())
  { // Edge condition, use simple eval
    simplePolyEval(ret, poly, babyStep);
    return;
  }
  NTL::ZZX r = trunc(poly, k * t);      // degree <= k*2^e-1
  NTL::ZZX q = RightShift(poly, k * t); // degree == k(2^e-1) +delta

  const NTL::ZZ p = NTL::to_ZZ(babyStep[0].getPtxtSpace());
  const NTL::ZZ &coef = coeff(r, deg(q));
  SetCoeff(r, deg(q), coef - 1); // r' = r - X^{deg(q)}

  NTL::ZZX c, s;
  DivRem(c, s, r, q); // r' = c*q + s
  // deg(s)<deg(q), and if c!= 0 then deg(c)<k-delta

  helib::assertTrue(deg(s) < deg(q), "Degree of s is not less than degree of q");
  helib::assertTrue(IsZero(c) || deg(c) < k - delta, "Nonzero c has not degree smaller than k - delta");
  SetCoeff(s, deg(q)); // s' = s + X^{deg(q)}, deg(s)==deg(q)

  // reduce the coefficients modulo p
  for (long i = 0; i <= deg(c); i++)
    rem(c[i], c[i], p);
  c.normalize();
  for (long i = 0; i <= deg(s); i++)
    rem(s[i], s[i], p);
  s.normalize();

  // Evaluate recursively poly = (c+X^{kt})*q + s'
  PatersonStockmeyer(ret, q, k, t / 2, delta, babyStep, giantStep);

  Ctxt tmp(ret.getPubKey(), ret.getPtxtSpace());
  simplePolyEval(tmp, c, babyStep);
  tmp += giantStep.getPower(t);
  ret.multiplyBy(tmp);

  PatersonStockmeyer(tmp, s, k, t / 2, delta, babyStep, giantStep);
  ret += tmp;
}

// This procedure assumes that k*(2^e +1) > deg(poly) > k*(2^e -1),
// and that babyStep contains >= k + (deg(poly) mod k) powers
void degPowerOfTwo(Ctxt &ret, const NTL::ZZX &poly, long k,
                   DynamicCtxtPowers &babyStep, DynamicCtxtPowers &giantStep)
{
  if (deg(poly) <= babyStep.size())
  { // Edge condition, use simple eval
    simplePolyEval(ret, poly, babyStep);
    return;
  }
  long n = divc(deg(poly), k);                // We assume n=2^e or n=2^e -1
  n = 1L << NTL::NextPowerOfTwo(n);           // round up to n=2^e
  NTL::ZZX r = trunc(poly, (n - 1) * k);      // degree <= k(2^e-1)-1
  NTL::ZZX q = RightShift(poly, (n - 1) * k); // 0 < degree < 2k
  SetCoeff(r, (n - 1) * k);                   // monic, degree == k(2^e-1)
  q -= 1;

  PatersonStockmeyer(ret, r, k, n / 2, 0, babyStep, giantStep);

  Ctxt tmp(ret.getPubKey(), ret.getPtxtSpace());
  simplePolyEval(tmp, q, babyStep); // evaluate q

  // multiply by X^{k(n-1)} with minimum depth
  for (long i = 1; i < n; i *= 2)
  {
    tmp.multiplyBy(giantStep.getPower(i));
  }
  ret += tmp;
}

void recursivePolyEval(Ctxt &ret, const NTL::ZZX &poly, long k,
                       DynamicCtxtPowers &babyStep, DynamicCtxtPowers &giantStep)
{
  if (deg(poly) <= babyStep.size())
  { // Edge condition, use simple eval
    simplePolyEval(ret, poly, babyStep);
    return;
  }

  long delta = deg(poly) % k;              // deg(poly) mod k
  long n = divc(deg(poly), k);             // ceil( deg(poly)/k )
  long t = 1L << (NTL::NextPowerOfTwo(n)); // t >= n, so t*k >= deg(poly)

  // Special case for deg(poly) = k * 2^e +delta
  if (n == t)
  {
    degPowerOfTwo(ret, poly, k, babyStep, giantStep);
    return;
  }

  // When deg(poly) = k*(2^e -1) we use the Paterson-Stockmeyer recursion
  if (n == t - 1 && delta == 0)
  {
    PatersonStockmeyer(ret, poly, k, t / 2, delta, babyStep, giantStep);
    return;
  }

  t = t / 2;

  // In any other case we have kt < deg(poly) < k(2t-1). We then set
  // u = deg(poly) - k*(t-1) and poly = q*X^u + r with deg(r)<u
  // and recurse on poly = (q-1)*X^u + (X^u+r)

  long u = deg(poly) - k * (t - 1);
  NTL::ZZX r = trunc(poly, u);      // degree <= u-1
  NTL::ZZX q = RightShift(poly, u); // degree == k*(t-1)
  q -= 1;
  SetCoeff(r, u); // degree == u

  PatersonStockmeyer(ret, q, k, t / 2, 0, babyStep, giantStep);

  Ctxt tmp = giantStep.getPower(u / k);
  if (delta != 0)
  { // if u is not divisible by k then compute it
    tmp.multiplyBy(babyStep.getPower(delta));
  }
  ret.multiplyBy(tmp);

  recursivePolyEval(tmp, r, k, babyStep, giantStep);
  ret += tmp;
}

// Function to find modulo inverse of a
int modInverse(int A, int M)
{
  int x, y;
  int g = gcdExtended(A, M, &x, &y);
  if (g != 1)
    throw invalid_argument("No inverse");
  else
  {

    // m is added to handle negative x
    int res = (x % M + M) % M;
    return res;
  }
}

// Function for extended Euclidean Algorithm
int gcdExtended(int a, int b, int *x, int *y)
{
  // Base Case
  if (a == 0)
  {
    *x = 0, *y = 1;
    return b;
  }

  // To store results of recursive call
  int x1, y1;
  int gcd = gcdExtended(b % a, a, &x1, &y1);

  // Update x and y using results of recursive
  // call
  *x = y1 - (b / a) * x1;
  *y = x1;

  return gcd;
}

int get_inverse(int nom, int dom, int p)
{
  int bottom_inverse = (modInverse(dom, (int)p)) % (int)p;
  if (bottom_inverse < 0)
  {
    bottom_inverse += (int)p;
  }
  int result = (nom * bottom_inverse) % p;
  if (result < 0)
  {
    result += p;
  }
  return result;
}

template <typename T, typename Allocator>
void print_vector(const vector<T, Allocator> &vect, int num_entries)
{
  cout << vect[0];
  for (int i = 1; i < min((int)vect.size(), num_entries); i++)
  {
    cout << ", " << vect[i];
  }
  cout << endl;
}

Ctxt multiplyMany(vector<Ctxt> &v)
{
  uint32_t num_entries = v.size();
  uint32_t depth = ceil(log2(num_entries));

  for (uint32_t d = 0; d < depth; d++)
  {
    uint32_t jump_factor = pow(2, d);
    uint32_t skip_factor = 2 * jump_factor;

    for (uint32_t i = 0; i < num_entries; i += skip_factor)
    {
      v[i].multiplyBy(v[i + jump_factor]);
    }
  }
  return v[0];
}

Ctxt addMany(vector<Ctxt> &v)
{
  uint32_t num_entries = v.size();
  uint32_t depth = ceil(log2(num_entries));

  for (uint32_t d = 0; d < depth; d++)
  {
    uint32_t jump_factor = pow(2, d);
    uint32_t skip_factor = 2 * jump_factor;

    for (uint32_t i = 0; i < num_entries; i += skip_factor)
    {
      v[i] += v[i + jump_factor];
    }
  }
  return v[0];
}

Ctxt addManySafe(vector<Ctxt> &v, const PubKey &pk)
{
  Ctxt result(pk);

  for (uint32_t i = 0; i < v.size(); i++)
  {
    result += v[i];
  }
  return result;
}

void addOneMod2(Ctxt &a)
{
    //   0 -> 1
    //   1 -> 0
    // f(x)-> -x+1
    a.negate();
    a.addConstant(NTL::ZZX(1));
}