#pragma comment(linker, "/STACK:5000000")
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cstdint>
#include <cctype>
#include <algorithm>
using namespace std;

struct BigInt {
    static const uint64_t BASE = (1ULL << 32);
    vector<uint32_t> a; // little-endian words

    BigInt() { a = {0}; }
    BigInt(uint64_t v) { *this = v; }

    BigInt& operator=(uint64_t v) {
        a.clear();
        if (v == 0) a.push_back(0);
        while (v) { a.push_back(uint32_t(v & 0xFFFFFFFF)); v >>= 32; }
        return *this;
    }

    bool isZero() const { return a.size() == 1 && a[0] == 0; }

    void trim() {
        while (a.size() > 1 && a.back() == 0) a.pop_back();
    }

    static int compare(const BigInt &A, const BigInt &B) {
        if (A.a.size() != B.a.size()) return (A.a.size() < B.a.size()) ? -1 : 1;
        for (int i = (int)A.a.size() - 1; i >= 0; --i)
            if (A.a[i] != B.a[i]) return (A.a[i] < B.a[i]) ? -1 : 1;
        return 0;
    }

    static BigInt add(const BigInt &A, const BigInt &B) {
        BigInt C;
        size_t n = max(A.a.size(), B.a.size());
        C.a.assign(n, 0);
        uint64_t carry = 0;
        for (size_t i = 0; i < n || carry; ++i) {
            if (i == C.a.size()) C.a.push_back(0);
            uint64_t sum = carry;
            if (i < A.a.size()) sum += A.a[i];
            if (i < B.a.size()) sum += B.a[i];
            C.a[i] = uint32_t(sum & 0xFFFFFFFF);
            carry = sum >> 32;
        }
        C.trim();
        return C;
    }

    static BigInt subtract(const BigInt &A, const BigInt &B) {
        BigInt C = A;
        uint64_t borrow = 0;
        for (size_t i = 0; i < B.a.size() || borrow; ++i) {
            uint64_t sub = borrow + (i < B.a.size() ? B.a[i] : 0);
            if (C.a[i] < sub) {
                C.a[i] = uint32_t((uint64_t(C.a[i]) + BASE) - sub);
                borrow = 1;
            } else {
                C.a[i] -= uint32_t(sub);
                borrow = 0;
            }
        }
        C.trim();
        return C;
    }

    static BigInt multiply(const BigInt &A, const BigInt &B) {
        cerr << "[multiply] A=" << A.toHex() << " B=" << B.toHex() << endl;
        BigInt C;
        C.a.assign(A.a.size() + B.a.size(), 0);
        for (size_t i = 0; i < A.a.size(); ++i) {
            uint64_t carry = 0;
            for (size_t j = 0; j < B.a.size() || carry; ++j) {
                if (i + j >= C.a.size()) C.a.push_back(0);
                uint64_t cur = C.a[i + j] + carry;
                if (j < B.a.size()) cur += uint64_t(A.a[i]) * uint64_t(B.a[j]);
                C.a[i + j] = uint32_t(cur & 0xFFFFFFFF);
                carry = cur >> 32;
            }
        }
        C.trim();
        cerr << "[multiply] => " << C.toHex() << endl;
        return C;
    }

    // ---- Simplified modBig for debugging ----
    // Hàm lấy modulo đúng
// Hàm chia lấy modulo: R = A % B
static BigInt modBig(const BigInt &A, const BigInt &B) {
    if (compare(A, B) < 0) return A; // A < B

    BigInt R = A;
    BigInt D = B;

    int n = (int)R.a.size() - (int)D.a.size();
    if (n < 0) return R;

    // Shift D left cho bằng độ dài R
    D.a.insert(D.a.begin(), n, 0);

    for (; n >= 0; --n) {
        // Lặp trừ số lớn
        while (compare(R, D) >= 0) {
            R = subtract(R, D);
        }
        if (!D.a.empty()) D.a.erase(D.a.begin()); // shift right 32-bit word
    }
    return R;
}



    static void shiftRight1(BigInt &A) {
        uint32_t carry = 0;
        for (int i = (int)A.a.size() - 1; i >= 0; --i) {
            uint64_t cur = ((uint64_t)carry << 32) | A.a[i];
            A.a[i] = uint32_t(cur >> 1);
            carry = uint32_t(cur & 1);
        }
        A.trim();
        if (A.a.empty()) A.a.push_back(0);
    }

    static BigInt power(BigInt base, BigInt exp, const BigInt &mod) {
        cerr << "[power] base=" << base.toHex() << " exp=" << exp.toHex() << " mod=" << mod.toHex() << endl;
        BigInt result(1);
        base = modBig(base, mod);
        int iter = 0;
        while (!exp.isZero()) {
            cerr << "  [power] iter=" << iter << " exp=" << exp.toHex() << endl;
            if (exp.a[0] & 1) {
                BigInt mul = multiply(result, base);
                result = modBig(mul, mod);
                cerr << "    [mulmod] result=" << result.toHex() << endl;
            }
            shiftRight1(exp);
            BigInt sq = multiply(base, base);
            base = modBig(sq, mod);
            ++iter;
        }
        cerr << "[power] done => " << result.toHex() << endl;
        return result;
    }

    static BigInt fromHexBigEndian(string s) {
        for (auto &c : s) c = toupper((unsigned char)c);
        if (s.size() % 8) s = string(8 - s.size() % 8, '0') + s;
        BigInt r; r.a.clear();
        for (int i = s.size() - 8; i >= 0; i -= 8) {
            string part = s.substr(i, 8);
            uint32_t w = strtoul(part.c_str(), nullptr, 16);
            r.a.push_back(w);
        }
        if (r.a.empty()) r.a.push_back(0);
        r.trim();
        return r;
    }

    // parse little-endian hex
    // input: little-endian hex string, ví dụ "FA6" hoặc "0FA6"
// output: BigInt chính xác
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
    r.a.clear();

    // chia mỗi 8 hex làm 1 word
    for (int i = (int)s.size(); i > 0; i -= 8) {
        int len = min(8, i);
        string part = s.substr(i - len, len);
        uint32_t w = strtoul(part.c_str(), nullptr, 16);
        r.a.push_back(w);
    }

    if (r.a.empty()) r.a.push_back(0);
    r.trim();
    return r;
}






    string toHex() const {
    if (isZero()) return "0";
    string s;

    // Duyệt từ word cao nhất về thấp nhất
    for (int i = (int)a.size() - 1; i >= 0; --i) {
        char buf[9];
        sprintf(buf, "%08X", a[i]); // chuyển 32-bit word thành 8 hex
        s += buf;
    }

    // loại bỏ leading zeros
    while (s.size() > 1 && s[0] == '0') s.erase(0, 1);
    return s;
}


   string toHexLittleEndian() const {
    string hex = toHex();         // Lấy chuỗi hex chuẩn
    reverse(hex.begin(), hex.end()); // Đảo ngược toàn bộ chuỗi ký tự
    return hex;
}

};

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    ifstream in("single.inp");
    if (!in.is_open()) {
        cerr << "❌ Không mở được single.inp\n";
        return 1;
    }

    string A_hex, B_hex, C_hex;
    in >> A_hex >> B_hex >> C_hex;
    in.close();

    BigInt A = BigInt::fromLittleEndianHex(A_hex);
    BigInt B = BigInt::fromLittleEndianHex(B_hex);
    BigInt C = BigInt::fromLittleEndianHex(C_hex);

    cerr << "DBG A(HEX)=" << A.toHex() << " B(HEX)=" << B.toHex() << " C(HEX)=" << C.toHex() << "\n";

    BigInt R = BigInt::power(C,B,A);
    cout << R.toHexLittleEndian() << "\n";
    return 0;
}
