#include <bits/stdc++.h>
using namespace std;

const int LIMBS = 32; // 1024-bit = 32 * 32
struct BigInt {
    uint32_t v[LIMBS] = {0};
};

// ---- Montgomery Core ----
struct MontgomeryCtx {
    BigInt N;       // modulus
    uint32_t n_inv; // -N^{-1} mod 2^32
    BigInt R2;      // R^2 mod N
};

// tính n_inv = -N^{-1} mod 2^32
uint32_t montgomery_inv32(uint32_t n) {
    uint32_t x = n;
    for (int i = 0; i < 5; ++i)
        x *= 2 - n * x;
    return ~x + 1;
}

// so sánh >=
bool geq(const BigInt &a, const BigInt &b) {
    for (int i = LIMBS - 1; i >= 0; --i) {
        if (a.v[i] > b.v[i]) return true;
        if (a.v[i] < b.v[i]) return false;
    }
    return true;
}

// trừ mod
void sub_mod(BigInt &a, const BigInt &b) {
    uint64_t borrow = 0;
    for (int i = 0; i < LIMBS; ++i) {
        uint64_t diff = (uint64_t)a.v[i] - b.v[i] - borrow;
        a.v[i] = (uint32_t)diff;
        borrow = (diff >> 63);
    }
}

// MontMul (a * b * R^-1 mod N)
void mont_mul(BigInt &res, const BigInt &a, const BigInt &b, const MontgomeryCtx &ctx) {
    uint64_t T[LIMBS * 2] = {0};

    // t = a * b
    for (int i = 0; i < LIMBS; ++i) {
        uint64_t carry = 0;
        for (int j = 0; j < LIMBS; ++j) {
            uint64_t prod = (uint64_t)a.v[i] * b.v[j] + T[i + j] + carry;
            T[i + j] = (uint32_t)prod;
            carry = prod >> 32;
        }
        T[i + LIMBS] += carry;
    }

    // Montgomery reduction
    for (int i = 0; i < LIMBS; ++i) {
        uint32_t m = (uint32_t)(T[i] * ctx.n_inv);
        uint64_t carry = 0;
        for (int j = 0; j < LIMBS; ++j) {
            uint64_t prod = (uint64_t)m * ctx.N.v[j] + T[i + j] + carry;
            T[i + j] = (uint32_t)prod;
            carry = prod >> 32;
        }
        uint64_t sum = T[i + LIMBS] + carry;
        T[i + LIMBS] = (uint32_t)sum;
        // phần cao hơn LIMBS+1 bỏ vì /R
    }

    // copy phần cao (T >> 32*LIMBS)
    for (int i = 0; i < LIMBS; ++i)
        res.v[i] = (uint32_t)T[i + LIMBS];

    // nếu >= N thì trừ N
    if (geq(res, ctx.N)) sub_mod(res, ctx.N);
}

// chuyển vào dạng Montgomery: a * R mod N
void to_mont(BigInt &r, const BigInt &a, const MontgomeryCtx &ctx) {
    mont_mul(r, a, ctx.R2, ctx);
}

// chuyển về dạng thường: a * R^-1 mod N
void from_mont(BigInt &r, const BigInt &a, const MontgomeryCtx &ctx) {
    BigInt one = {};
    one.v[0] = 1;
    mont_mul(r, a, one, ctx);
}

// ---- Sliding Window Exponentiation ----
void mont_pow(BigInt &res, const BigInt &x, const BigInt &k, const MontgomeryCtx &ctx, int w) {
    BigInt x1;
    to_mont(x1, x, ctx);

    int table_size = 1 << w;
    vector<BigInt> table(table_size);
    table[1] = x1;
    for (int i = 3; i < table_size; i += 2)
        mont_mul(table[i], table[i - 2], x1, ctx);

    BigInt result;
    BigInt one = {};
    one.v[0] = 1;
    to_mont(result, one, ctx);

    for (int i = LIMBS * 32 - 1; i >= 0;) {
        if (((k.v[i / 32] >> (i % 32)) & 1) == 0) {
            mont_mul(result, result, result, ctx);
            --i;
        } else {
            int l = max(0, i - w + 1);
            while (((k.v[l / 32] >> (l % 32)) & 1) == 0) ++l;
            int win = 0;
            for (int j = l; j <= i; ++j)
                win = (win << 1) | ((k.v[j / 32] >> (j % 32)) & 1);

            for (int j = 0; j < i - l + 1; ++j)
                mont_mul(result, result, result, ctx);

            mont_mul(result, result, table[win], ctx);
            i = l - 1;
        }
    }

    from_mont(res, result, ctx);
}


void mont_init(MontgomeryCtx &ctx) {
    // n_inv = -N^{-1} mod 2^32
    ctx.n_inv = montgomery_inv32(ctx.N.v[0]);

    // Tính R^2 mod N
    BigInt R2 = {};
    BigInt one = {};
    one.v[0] = 1;

    // R2 = (1 << (64 * LIMBS)) mod N
    // => dùng repeated squaring hoặc shift mod
    BigInt R = {};
    R.v[0] = 1;
    for (int i = 0; i < LIMBS * 32 * 2; i++) {
        // R = (R << 1) mod N
        uint64_t carry = 0;
        for (int j = 0; j < LIMBS; ++j) {
            uint64_t v = ((uint64_t)R.v[j] << 1) | carry;
            R.v[j] = (uint32_t)v;
            carry = v >> 32;
        }
        if (geq(R, ctx.N)) sub_mod(R, ctx.N);
    }
    ctx.R2 = R;
}


// Chuyển little-endian hex string -> BigInt
static BigInt fromLittleEndianHex(const string &s_in) {
    string s = s_in;
    if (s.size() >= 2 && s[0]=='0' && (s[1]=='x'||s[1]=='X'))
        s = s.substr(2);
    if (s.empty()) s = "0";

    // đảo toàn bộ chuỗi theo ký tự
    reverse(s.begin(), s.end());

    // pad để đủ số hex chẵn
    if (s.size() % 2) s = "0" + s;

    BigInt r;
    r.v[LIMBS - 1] = 0; // khởi tạo

    // chia mỗi 8 hex làm 1 word
    for (int i = (int)s.size(); i > 0; i -= 8) {
        int len = min(8, i);
        string part = s.substr(i - len, len);
        uint32_t w = strtoul(part.c_str(), nullptr, 16);
        int index = (s.size() - i) / 8;
        if (index < LIMBS)
            r.v[index] = w;
    }

    return r;
}

// Big-endian hex string -> little-endian hex string
string toBigEndianHex(const string &littleHex)
{
    string s = littleHex;
    if (s.size() <= 2)
        return s; // ✅ không đảo nếu chỉ 1 byte
    if (s.size() % 2)
        s = "0" + s;
    string res;
    for (int i = s.size(); i > 0; i -= 2)
        res += s.substr(i - 2, 2);
    return res;
}

static BigInt fromHex(const string &s_in) {
    string s = s_in;
    if (s.size() >= 2 && s[0]=='0' && (s[1]=='x'||s[1]=='X'))
        s = s.substr(2);
    if (s.empty()) s = "0";

    if (s.size() % 8 != 0)
        s = string(8 - s.size() % 8, '0') + s;

    BigInt r;
    int idx = 0;
    for (int i = (int)s.size(); i > 0; i -= 8) {
        if (idx >= LIMBS) break;
        string part = s.substr(i - 8, 8);
        r.v[idx++] = strtoul(part.c_str(), nullptr, 16);
    }
    return r;
}


string toHex(const BigInt &a) {
    stringstream ss;
    bool leading = false;
    for (int i = LIMBS - 1; i >= 0; --i) {
        if (!leading && a.v[i] == 0) continue;
        ss << hex << uppercase << setw(8) << setfill('0') << a.v[i];
        leading = true;
    }
    if (!leading) return "0";
    string s = ss.str();
    while (s.size() > 1 && s[0] == '0') s.erase(s.begin());
    return s;
}

string fromLittleEndianToDec(const string &bigHex) {
    BigInt a = fromLittleEndianHex(bigHex);
    if (a.v[0] == 0 && std::all_of(a.v, a.v + LIMBS, [](uint32_t x){ return x == 0; }))
        return "0";

    string res;
    BigInt zero = {};
    BigInt ten = {};
    ten.v[0] = 10;

    BigInt temp = a;
    while (std::any_of(temp.v, temp.v + LIMBS, [](uint32_t x){ return x != 0; })) {
        // Chia temp cho 10
        uint64_t rem = 0;
        for (int i = LIMBS - 1; i >= 0; --i) {
            uint64_t cur = (rem << 32) | temp.v[i];
            temp.v[i] = (uint32_t)(cur / 10);
            rem = cur % 10;
        }
        // Ghi lại phần dư
        res.push_back('0' + (char)rem);
    }
    reverse(res.begin(), res.end());
    return res;
}


int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    string N_hex, k_hex, x_hex;
    cin >> N_hex >> k_hex >> x_hex;

    cout << "N:" << fromLittleEndianToDec(N_hex) << "\n";
    cout << "k:" << fromLittleEndianToDec(k_hex) << "\n";
    cout << "x:" << fromLittleEndianToDec(x_hex) << "\n";

    BigInt N = fromHex(N_hex);
    BigInt k = fromHex(k_hex);
    BigInt x = fromHex(x_hex);

    cout << "N: " << toHex(N) << "\n";
    cout << "k: " << toHex(k) << "\n";
    cout << "x: " << toHex(x) << "\n";

    MontgomeryCtx ctx;
    ctx.N = N;
    mont_init(ctx);

    BigInt y;
    mont_pow(y, x, k, ctx, 4); // cửa sổ kích thước 4

    cout << toHex(y) << "\n";
}