#include <bits/stdc++.h>
using namespace std;

struct BigInt {
    vector<uint32_t> v; // dynamic limbs, little-endian (v[0] = least significant)

    BigInt(int limbs = 0) { v.resize(limbs, 0); }

    int size() const { return (int)v.size(); }
};

// Count bits in hex string
int countBitsFromHex(const string &hexStr) {
    string s = hexStr;
    if (s.size() >= 2 && s[0]=='0' && (s[1]=='x'||s[1]=='X')) s = s.substr(2);
    size_t firstNonZero = s.find_first_not_of('0');
    if (firstNonZero == string::npos) return 0;
    s = s.substr(firstNonZero);
    int bits = s.size() * 4;

    char firstChar = s[0];
    int val = (isdigit(firstChar)) ? firstChar - '0' :
              (isupper(firstChar)) ? firstChar - 'A' + 10 :
              firstChar - 'a' + 10;

    int leadingZeros = 0;
    for (int i = 3; i >= 0; --i) {
        if ((val >> i) & 1) break;
        leadingZeros++;
    }
    bits -= leadingZeros;
    return bits;
}

// Convert hex string (little-endian) -> BigInt
BigInt fromLittleEndianHex(const string &s_in) {
    string s = s_in;
    if (s.size() >= 2 && s[0]=='0' && (s[1]=='x'||s[1]=='X')) s = s.substr(2);
    reverse(s.begin(), s.end());
    if (s.size()%2) s = "0"+s;
    int limbs = (s.size() + 7)/8;
    BigInt r(limbs);
    for (int i = (int)s.size(); i > 0; i -= 8) {
        int len = min(8,i);
        string part = s.substr(i-len, len);
        r.v[(s.size()-i)/8] = strtoul(part.c_str(), nullptr, 16);
    }
    return r;
}

// Convert BigInt -> little-endian hex string
string toLittleEndianHex(const BigInt &a) {
    stringstream ss;
    bool leading = false;
    for (int i = (int)a.v.size()-1; i>=0; --i) {
        if (!leading && a.v[i]==0) continue;
        ss << hex << uppercase << setw(8) << setfill('0') << a.v[i];
        leading = true;
    }
    if (!leading) return "0";
    string s = ss.str();
    while (s.size()>1 && s[0]=='0') s.erase(s.begin());
    reverse(s.begin(), s.end());
    return s;
}

// Compare a >= b
bool geq(const BigInt &a, const BigInt &b) {
    if (a.size() != b.size()) return a.size() > b.size();
    for (int i = (int)a.size()-1; i>=0; --i) {
        if (a.v[i] > b.v[i]) return true;
        if (a.v[i] < b.v[i]) return false;
    }
    return true;
}

// Subtract a -= b (assumes a >= b)
void sub_mod(BigInt &a, const BigInt &b) {
    uint64_t borrow = 0;
    for (size_t i=0; i<a.size(); ++i) {
        uint64_t diff = (uint64_t)a.v[i] - (i<b.size()?b.v[i]:0) - borrow;
        a.v[i] = (uint32_t)diff;
        borrow = (diff > a.v[i]) ? 1 : 0;
    }
}

// Montgomery context
struct MontgomeryCtx {
    BigInt N;
    uint32_t n_inv;
    BigInt R2;
};

// Compute -N^-1 mod 2^32
uint32_t montgomery_inv32(uint32_t n) {
    uint32_t x = n;
    for (int i =0; i<5; ++i) x *= 2 - n*x;
    return ~x +1;
}

// Montgomery multiplication
BigInt mont_mul(const BigInt &a, const BigInt &b, const MontgomeryCtx &ctx) {
    int limbs = (int)a.size();
    vector<uint64_t> T(limbs*2+2,0);

    for (int i=0;i<limbs;++i) {
        uint64_t carry=0;
        for (int j=0;j<limbs;++j) {
            uint64_t prod = (uint64_t)a.v[i]*b.v[j]+T[i+j]+carry;
            T[i+j] = (uint32_t)prod;
            carry = prod >>32;
        }
        T[i+limbs] += carry;
    }

    for (int i=0;i<limbs;++i) {
        uint32_t m = (uint32_t)(T[i]*ctx.n_inv);
        uint64_t carry=0;
        for (int j=0;j<limbs;++j){
            uint64_t prod = (uint64_t)m*ctx.N.v[j]+T[i+j]+carry;
            T[i+j] = (uint32_t)prod;
            carry = prod >>32;
        }
        T[i+limbs] += carry;
    }

    BigInt res(limbs);
    for (int i=0;i<limbs;++i) res.v[i]=(uint32_t)T[i+limbs];

    if (geq(res, ctx.N)) sub_mod(res, ctx.N);
    return res;
}

// Convert to/from Montgomery
BigInt to_mont(const BigInt &a, const MontgomeryCtx &ctx) { return mont_mul(a, ctx.R2, ctx); }
BigInt from_mont(const BigInt &a, const MontgomeryCtx &ctx) {
    BigInt one(a.size());
    one.v[0] = 1;
    return mont_mul(a, one, ctx);
}

// Binary exponentiation
BigInt mont_pow(const BigInt &x, const BigInt &k, const MontgomeryCtx &ctx) {
    BigInt x1 = to_mont(x, ctx);
    BigInt result(x.size());
    result.v[0]=1;
    result = to_mont(result, ctx);

    int bits = (int)k.size()*32;
    for (int i=bits-1;i>=0;--i){
        result = mont_mul(result,result,ctx);
        if ((k.v[i/32]>>(i%32))&1) result = mont_mul(result,x1,ctx);
    }
    return from_mont(result, ctx);
}

// Initialize Montgomery context
MontgomeryCtx mont_init(const BigInt &N) {
    MontgomeryCtx ctx;
    ctx.N = N;
    ctx.n_inv = montgomery_inv32(N.v[0]);

    BigInt R2(N.size());
    R2.v[0]=1;
    int totalBits = N.size()*32*2;
    for (int i=0;i<totalBits;++i){
        uint64_t carry=0;
        for (size_t j=0;j<R2.size();++j){
            uint64_t tmp = ((uint64_t)R2.v[j]<<1)|carry;
            R2.v[j]=(uint32_t)tmp;
            carry=tmp>>32;
        }
        if (geq(R2,N)) sub_mod(R2,N);
    }
    ctx.R2 = R2;
    return ctx;
}

// Demo main
int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    string N_hex, k_hex, x_hex;
    ifstream infile("input.inp");
    infile >> N_hex >> k_hex >> x_hex;
    infile.close();

    int bitsN = countBitsFromHex(N_hex);
    int limbs = (bitsN+31)/32;
    cout << "N has " << bitsN << " bits, using " << limbs << " limbs" << endl;

    BigInt N = fromLittleEndianHex(N_hex);
    BigInt k = fromLittleEndianHex(k_hex);
    BigInt x = fromLittleEndianHex(x_hex);

    MontgomeryCtx ctx = mont_init(N);
    BigInt y = mont_pow(x,k,ctx);

    cout << toLittleEndianHex(y) << "\n";
}
