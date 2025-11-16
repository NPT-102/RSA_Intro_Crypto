#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <cassert>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <chrono>

std::string readFile(const std::string& filename) {
  std::string line;
  std::ifstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Failed to open file: " << filename << std::endl;
    return "";
  }

  std::getline(file, line);

  // Remove trailing whitespace (space, tab, etc.)
  size_t end = line.find_last_not_of(" \t\r\n");
  if (end != std::string::npos) {
    line.erase(end + 1);
  } else {
    line.clear(); // all whitespace
  }

  std::cout << "File read successfully" << std::endl;
  return line;
}

struct BigInt {
  std::vector<uint64_t> limbs;

  BigInt() {}
  BigInt(uint64_t val) {
    if (val) limbs.push_back(val);
  }

  static BigInt one() { return BigInt(1); }
  static BigInt zero() { return BigInt(0); }

  bool is_zero() const {
    for (auto x : limbs) if (x) return false;
    return true;
  }

  bool is_even() const { return limbs.empty() || !(limbs[0] & 1); }
  BigInt operator-(const BigInt& rhs) const {
    BigInt result;
    result.limbs.resize(limbs.size());
    uint64_t borrow = 0;
    for (size_t i = 0; i < limbs.size(); ++i) {
      uint64_t a = limbs[i];
      uint64_t b = i < rhs.limbs.size() ? rhs.limbs[i] : 0;
      uint64_t sub = a - b - borrow;
      borrow = (a < b + borrow);
      result.limbs[i] = sub;
    }
    while (!result.limbs.empty() && result.limbs.back() == 0)
      result.limbs.pop_back();
    return result;
  }

  BigInt operator%(const BigInt& mod) const {
    BigInt r = *this;
    while (r >= mod)
      r = r - mod;
    return r;
  }

  BigInt operator*(const BigInt& rhs) const {
    BigInt result;
    result.limbs.resize(limbs.size() + rhs.limbs.size(), 0);
    for (size_t i = 0; i < limbs.size(); ++i) {
      uint64_t carry = 0;
      for (size_t j = 0; j < rhs.limbs.size(); ++j) {
        __uint128_t mul = (__uint128_t)limbs[i] * rhs.limbs[j] + result.limbs[i + j] + carry;
        result.limbs[i + j] = (uint64_t)mul;
        carry = (uint64_t)(mul >> 64);
      }
      result.limbs[i + rhs.limbs.size()] += carry;
    }
    while (!result.limbs.empty() && result.limbs.back() == 0)
      result.limbs.pop_back();
    return result;
  }

  BigInt operator>>(int shift) const {
    BigInt result = *this;
    for (int i = 0; i < shift; ++i) {
      if (result.limbs.empty()) break;
      uint64_t carry = 0;
      for (size_t j = result.limbs.size(); j-- > 0; ) {
        uint64_t new_carry = result.limbs[j] & 1;
        result.limbs[j] = (result.limbs[j] >> 1) | (carry << 63);
        carry = new_carry;
      }
      while (!result.limbs.empty() && result.limbs.back() == 0)
        result.limbs.pop_back();
    }
    return result;
  }

  bool operator==(const BigInt& rhs) const {
    return limbs == rhs.limbs;
  }

  bool operator<(const BigInt& rhs) const {
    if (limbs.size() != rhs.limbs.size())
      return limbs.size() < rhs.limbs.size();
    for (int i = limbs.size() - 1; i >= 0; --i) {
      if (limbs[i] != rhs.limbs[i])
        return limbs[i] < rhs.limbs[i];
    }
    return false;
  }

  bool operator>=(const BigInt& rhs) const {
    if (limbs.size() != rhs.limbs.size())
      return limbs.size() > rhs.limbs.size();
    for (int i = limbs.size() - 1; i >= 0; --i) {
      if (limbs[i] != rhs.limbs[i])
        return limbs[i] > rhs.limbs[i];
    }
    return true;
  }

  void print() const {
    for (int i = limbs.size() - 1; i >= 0; --i)
      std::cout << std::hex << limbs[i] << " ";
    std::cout << std::dec << "\n";
  }
  void normalize() {
    while (!limbs.empty() && limbs.back() == 0)
      limbs.pop_back();
  }
};

// ADDITION
BigInt add_bigint(const BigInt &a, const BigInt &b) {
    BigInt r;
    size_t k = std::max(a.limbs.size(), b.limbs.size());
    r.limbs.assign(k, 0);
    unsigned __int128 carry = 0;
    for (size_t i = 0; i < k; ++i) {
        unsigned __int128 av = (i < a.limbs.size()) ? a.limbs[i] : 0;
        unsigned __int128 bv = (i < b.limbs.size()) ? b.limbs[i] : 0;
        unsigned __int128 sum = av + bv + carry;
        r.limbs[i] = (uint64_t)sum;
        carry = sum >> 64;
    }
    if (carry) r.limbs.push_back((uint64_t)carry);
    r.normalize();
    return r;
}

BigInt fromHex(std::string hex) {
  BigInt n;
  n.limbs.clear();
  int len = hex.size();
  int bitsPerLimb = 16; // 16 hex digits = 64 bits
  for (int i = len; i > 0; i -= bitsPerLimb) {
    int start = std::max(0, i - bitsPerLimb);
    std::string chunk = hex.substr(start, i - start);
    uint64_t limb = std::stoull(chunk, nullptr, 16);
    n.limbs.push_back(limb);
  }
  return n;
}

std::string toHex(const BigInt &n) {
  if (n.limbs.empty()) return "0";
  std::stringstream ss;
  ss << std::hex << std::uppercase << std::setfill('0');
  bool started = false;

  for (int i = n.limbs.size() - 1; i >= 0; --i) {
    if (!started) {
      ss << std::hex << std::uppercase << n.limbs[i];
      started = true;
    }
    else {
      ss << std::setw(16) << n.limbs[i];
    }
  }
  return ss.str();
}

BigInt mod_pow(BigInt base, BigInt exp, const BigInt& mod) {
  BigInt result = BigInt::one();
  base = base % mod;
  while (!exp.is_zero()) {
    if (!exp.is_even())
      result = (result * base) % mod;
    exp = exp >> 1;
    base = (base * base) % mod;
  }
  return result;
}

BigInt mod_pow_montgomery(BigInt base, BigInt exp, const BigInt &mod);

bool miller_rabin(const BigInt& n, int k) {
  if (n.is_even() || n == BigInt::one()) return false;

  BigInt d = n - BigInt::one();
  int r = 0;
  while (d.is_even()) {
    d = d >> 1;
    r++;
  }

  for (int i = 0; i < k; ++i) {
    BigInt a = BigInt(2 + i); 
    BigInt x = mod_pow_montgomery(a, d, n);
    if (x == BigInt::one() || x == n - BigInt::one())
      continue;
    bool continue_outer = false;
    for (int j = 0; j < r - 1; ++j) {
      x = mod_pow_montgomery(x, BigInt(2), n);
      if (x == n - BigInt::one()) {
        continue_outer = true;
        break;
      }
    }
    if (continue_outer) continue;
    return false;
  }
  return true;
}

BigInt compute_R2_mod(const BigInt &n) {
    size_t k = n.limbs.size();
    BigInt r = BigInt::one();        
    size_t total_bits = 64 * k * 2;  // R^2 = 2^(64*2*k)
    for (size_t i = 0; i < total_bits; ++i) {
        // r = (r + r) mod n
        r = add_bigint(r, r);
        if (r >= n) r = r - n;
    }
    return r; // r == R^2 mod n
}

uint64_t montgomery_inv64(uint64_t n0) {
  uint64_t x = 1;
  for (int i = 0; i < 6; ++i)
    x*= 2 - n0 * x;
  return ~x + 1;
}

struct MontgomeryContext {
  BigInt n;
  uint64_t n_inv;
  size_t k;
};

MontgomeryContext montgomery_prepare(const BigInt &mod) {
  MontgomeryContext ctx;
  ctx.n = mod;
  ctx.k = mod.limbs.size();
  ctx.n_inv = montgomery_inv64(mod.limbs[0]);
  return ctx;
}

BigInt montgomery_mul(const BigInt &a, const BigInt &b, const MontgomeryContext &ctx) {
  size_t k = ctx.k;
  const auto &n = ctx.n.limbs;

  BigInt aa = a, bb = b;
  while (aa >= ctx.n) aa = aa - ctx.n;
  while (bb >= ctx.n) bb = bb - ctx.n;
  if (aa.limbs.size() < k) aa.limbs.resize(k, 0);
  if (bb.limbs.size() < k) bb.limbs.resize(k, 0);


  BigInt t;
  t.limbs.assign(k * 2 + 3, 0);

  // multiply a * b

  for (size_t i = 0; i < k; ++i) {
    __uint128_t carry = 0;
    for (size_t j = 0; j < k; ++j) {
      __uint128_t cur = (__uint128_t)aa.limbs[i] * bb.limbs[j] + t.limbs[i + j] + carry;
      t.limbs[i + j] = (uint64_t)cur;
      carry = cur >> 64;
    }
    t.limbs[i + k] += (uint64_t)carry;
  }

  //montgomery reducction
  for (size_t i = 0; i < k; ++i){
    uint64_t m = (uint64_t) ((__uint128_t)t.limbs[i] * ctx.n_inv);
    __uint128_t carry = 0;
    for (size_t j = 0; j < k; ++j) {
      __uint128_t cur = (__uint128_t)m * n[j] + t.limbs[i + j] + carry;
      t.limbs[i + j] = (uint64_t)cur;
      carry = cur >> 64;
    }
    size_t pos = i + k;
    __uint128_t cur = (__uint128_t)t.limbs[pos] + carry;
    t.limbs[pos] = (uint64_t)cur;
    carry = cur >> 64;
    if (carry && pos + 1 < t.limbs.size())
      t.limbs[pos + 1] += (uint64_t)carry;
  }

  BigInt res;
  res.limbs.assign(t.limbs.begin() + k, t.limbs.begin() + 2 * k);
  res.normalize();

  if (res >= ctx.n)
    res = res - ctx.n;

  return res;
}


BigInt mod_pow_montgomery(BigInt base, BigInt exp, const BigInt &mod) {
  MontgomeryContext ctx = montgomery_prepare(mod);
  BigInt R2 = compute_R2_mod(mod);
  BigInt baseM = montgomery_mul(base % mod, R2, ctx);
  BigInt resultM = montgomery_mul(BigInt::one(), R2, ctx);

  while (!exp.is_zero()) {
    if (!exp.is_even())
      resultM = montgomery_mul(resultM, baseM, ctx);
    exp = exp >> 1;
    baseM = montgomery_mul(baseM, baseM, ctx);
  }

  //convert back
  BigInt res = montgomery_mul(resultM, BigInt::one(), ctx);
  return res;
}

void makeOutputFile(std::string filename, bool result) {
  std::ofstream fo(filename);
  if (!fo.is_open()) {
    std::cerr << "Failed to open file: " << filename << std::endl;
    return;
  }

  fo << (result ? 0 : 1) << '\n';
  fo.close();
}


int main() {
  std::string hexN = readFile("./project_01_01/test_19.inp");
  BigInt n = fromHex(hexN);

  auto start_time = std::chrono::high_resolution_clock::now();
  bool is_prime = miller_rabin(n, 20);
  auto end_time = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::milli> duration = end_time - start_time;


  makeOutputFile("output.txt", is_prime);
  std::cout << "Number is " << (is_prime ? "prime" : "composite") << std::endl;
  std::cout << "Time: " << duration.count() << "ms\n";
  return 0;
}
