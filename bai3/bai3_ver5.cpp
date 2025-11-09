#include <bits/stdc++.h>
#include <fstream>
using namespace std;

const int LIMBS = 32; // 1024-bit
struct BigInt
{
    uint32_t v[LIMBS] = {0};
};

// ----- Montgomery Core -----
struct MontgomeryCtx
{
    BigInt N;       // modulus
    uint32_t n_inv; // -N^{-1} mod 2^32
    BigInt R2;      // R^2 mod N
};

// ----- Utils -----
uint32_t montgomery_inv32(uint32_t n)
{
    uint32_t x = n;
    for (int i = 0; i < 5; ++i)
        x *= 2 - n * x;
    return ~x + 1;
}

bool geq(const BigInt &a, const BigInt &b)
{
    for (int i = LIMBS - 1; i >= 0; --i)
    {
        if (a.v[i] > b.v[i])
            return true;
        if (a.v[i] < b.v[i])
            return false;
    }
    return true;
}

void sub_mod(BigInt &a, const BigInt &b)
{
    uint64_t borrow = 0;
    for (int i = 0; i < LIMBS; ++i)
    {
        uint64_t diff = (uint64_t)a.v[i] - b.v[i] - borrow;
        a.v[i] = (uint32_t)diff;
        borrow = (diff > a.v[i] ? 1 : 0);
    }
}

// ----- Debug -----
void printBigInt(const string &name, const BigInt &a)
{
    cout << name << ": ";
    for (int i = LIMBS - 1; i >= 0; --i)
        cout << hex << setw(8) << setfill('0') << a.v[i];
    cout << endl;
}

// ----- Montgomery Multiplication -----
void mont_mul(BigInt &res, const BigInt &a, const BigInt &b, const MontgomeryCtx &ctx)
{
    uint64_t T[LIMBS * 2 + 2] = {0};

    // Step 1: a * b
    for (int i = 0; i < LIMBS; ++i)
    {
        uint64_t carry = 0;
        for (int j = 0; j < LIMBS; ++j)
        {
            uint64_t prod = (uint64_t)a.v[i] * b.v[j] + T[i + j] + carry;
            T[i + j] = (uint32_t)prod;
            carry = prod >> 32;
        }
        T[i + LIMBS] += carry;
    }

    // Step 2: Montgomery reduction
    for (int i = 0; i < LIMBS; ++i)
    {
        uint32_t m = (uint32_t)(T[i] * ctx.n_inv);
        uint64_t carry = 0;
        for (int j = 0; j < LIMBS; ++j)
        {
            uint64_t prod = (uint64_t)m * ctx.N.v[j] + T[i + j] + carry;
            T[i + j] = (uint32_t)prod;
            carry = prod >> 32;
        }
        T[i + LIMBS] += carry;
    }

    // Step 3: copy phần cao
    for (int i = 0; i < LIMBS; ++i)
        res.v[i] = (uint32_t)T[i + LIMBS];

    // Step 4: reduce nếu >= N
    if (geq(res, ctx.N))
        sub_mod(res, ctx.N);
}

// ----- Montgomery conversion -----
void to_mont(BigInt &r, const BigInt &a, const MontgomeryCtx &ctx)
{
    mont_mul(r, a, ctx.R2, ctx);
}

void from_mont(BigInt &r, const BigInt &a, const MontgomeryCtx &ctx)
{
    BigInt one = {};
    one.v[0] = 1;
    mont_mul(r, a, one, ctx);
}

// ----- Binary exponentiation -----
void mont_pow_simple(BigInt &res, const BigInt &x, const BigInt &k, const MontgomeryCtx &ctx)
{
    BigInt x1;
    to_mont(x1, x, ctx);

    BigInt result;
    BigInt one = {};
    one.v[0] = 1;
    to_mont(result, one, ctx);

    for (int i = LIMBS * 32 - 1; i >= 0; --i)
    {
        mont_mul(result, result, result, ctx);
        if ((k.v[i / 32] >> (i % 32)) & 1)
            mont_mul(result, result, x1, ctx);
    }

    from_mont(res, result, ctx);
}

// ----- Montgomery init -----
void mont_init(MontgomeryCtx &ctx)
{
    ctx.n_inv = montgomery_inv32(ctx.N.v[0]);

    // Tính R^2 mod N bằng exponentiation an toàn
    BigInt R = {};
    R.v[LIMBS - 1] = 1; // R = 2^(32*LIMBS)
    BigInt one = {};
    one.v[0] = 1;
    BigInt R2 = one;

    // R2 = R^2 mod N = 2^(2*32*LIMBS) mod N
    // dùng shift an toàn
    for (int i = 0; i < LIMBS * 32 * 2; ++i)
    {
        uint64_t carry = 0;
        for (int j = 0; j < LIMBS; ++j)
        {
            uint64_t tmp = ((uint64_t)R2.v[j] << 1) | carry;
            R2.v[j] = (uint32_t)tmp;
            carry = tmp >> 32;
        }
        if (geq(R2, ctx.N))
            sub_mod(R2, ctx.N);
    }
    ctx.R2 = R2;
}

// ----- Hex <-> BigInt conversion -----
BigInt fromLittleEndianHex(const string &s_in)
{
    string s = s_in;
    if (s.size() >= 2 && s[0] == '0' && (s[1] == 'x' || s[1] == 'X'))
        s = s.substr(2);
    if (s.empty())
        s = "0";

    reverse(s.begin(), s.end());
    if (s.size() % 2)
        s = "0" + s;

    BigInt r;
    for (int i = (int)s.size(); i > 0; i -= 8)
    {
        int len = min(8, i);
        string part = s.substr(i - len, len);
        uint32_t w = strtoul(part.c_str(), nullptr, 16);
        int index = (s.size() - i) / 8;
        if (index < LIMBS)
            r.v[index] = w;
    }
    return r;
}

string toLittleEndianHex(const BigInt &a)
{
    stringstream ss;
    bool leading = false;
    for (int i = LIMBS - 1; i >= 0; --i)
    {
        if (!leading && a.v[i] == 0)
            continue;
        ss << hex << uppercase << setw(8) << setfill('0') << a.v[i];
        leading = true;
    }
    if (!leading)
        return "0";
    string s = ss.str();
    while (s.size() > 1 && s[0] == '0')
        s.erase(s.begin());
    reverse(s.begin(), s.end());
    return s;
}

// ----- Main -----
int main()
{
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    // Test files from test_00.inp to test_09.inp
    for (int test_num = 0; test_num <= 9; ++test_num)
    {
        char test_num_str[4];
        sprintf(test_num_str, "%02d", test_num);
        
        string test_input_file = string("project_01_03/test_") + test_num_str + ".inp";
        string test_output_file = string("project_01_03/test_") + test_num_str + ".out";

        // Read input
        ifstream infile(test_input_file);
        if (!infile.is_open())
        {
            cout << "false" << endl;
            continue;
        }

        string N_hex, k_hex, x_hex;
        infile >> N_hex >> k_hex >> x_hex;
        infile.close();

        // Parse input
        BigInt N = fromLittleEndianHex(N_hex);
        BigInt k = fromLittleEndianHex(k_hex);
        BigInt x = fromLittleEndianHex(x_hex);

        // Initialize Montgomery context
        MontgomeryCtx ctx;
        ctx.N = N;
        mont_init(ctx);

        // Compute result
        BigInt y;
        mont_pow_simple(y, x, k, ctx);
        string result = toLittleEndianHex(y);

        // Read expected output
        ifstream outfile(test_output_file);
        if (!outfile.is_open())
        {
            cout << "false" << endl;
            continue;
        }

        string expected;
        outfile >> expected;
        outfile.close();

        // Compare and print result
        bool matched = (result == expected);
        cout << test_num_str << ": " << (matched ? "true" : "false") << endl;
    }

    return 0;
}
