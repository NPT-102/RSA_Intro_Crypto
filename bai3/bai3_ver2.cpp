#include <bits/stdc++.h>
using namespace std;

uint64_t modMul(uint64_t a, uint64_t b, uint64_t mod) {
    __uint128_t res = (__uint128_t)a * b;
    return (uint64_t)(res % mod);
}

uint64_t modPow(uint64_t base, uint64_t exp, uint64_t mod) {
    uint64_t result = 1;
    base %= mod;
    while (exp) {
        if (exp & 1) result = modMul(result, base, mod);
        base = modMul(base, base, mod);
        exp >>= 1;
    }
    return result;
}

// Chuyển little-endian hex string -> uint64_t
uint64_t fromLittleEndianHex(const string &s) {
    string hex = s;
    reverse(hex.begin(), hex.end()); // đảo ký tự
    return strtoull(hex.c_str(), nullptr, 16);
}

// uint64_t -> little-endian hex string
string toHexLittleEndian(uint64_t x) {
    stringstream ss;
    ss << hex << uppercase << x;
    string hexstr = ss.str();
    reverse(hexstr.begin(), hexstr.end());
    return hexstr;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    string N_hex, k_hex, x_hex;
    cin >> N_hex >> k_hex >> x_hex;

    uint64_t N = fromLittleEndianHex(N_hex);
    uint64_t k = fromLittleEndianHex(k_hex);
    uint64_t x = fromLittleEndianHex(x_hex);

    uint64_t y = modPow(x, k, N);

    cout << toHexLittleEndian(y) << "\n";
}
