#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <random>
#include <algorithm>

using namespace std;

struct BigInt
{
    static const uint64_t BASE = (1ULL << 32);
    vector<uint32_t> a; // little-endian: a[0] là word thấp nhất

    BigInt() { a.push_back(0); }
    BigInt(uint64_t v) { *this = v; }

    BigInt &operator=(uint64_t v)
    {
        a.clear();
        if (v == 0)
        {
            a.push_back(0);
            return *this;
        }
        while (v)
        {
            a.push_back(v % BASE);
            v /= BASE;
        }
        return *this;
    }

    bool isZero() const
    {
        return a.size() == 1 && a[0] == 0;
    }

    void trim()
    {
        while (a.size() > 1 && a.back() == 0)
        {
            a.pop_back();
        }
    }

    static int compare(const BigInt &A, const BigInt &B)
    {
        if (A.a.size() != B.a.size())
        {
            return (A.a.size() < B.a.size()) ? -1 : 1;
        }
        for (int i = (int)A.a.size() - 1; i >= 0; i--)
        {
            if (A.a[i] != B.a[i])
            {
                return (A.a[i] < B.a[i]) ? -1 : 1;
            }
        }
        return 0;
    }

    static BigInt add(const BigInt &A, const BigInt &B)
    {
        BigInt C;
        size_t n = max(A.a.size(), B.a.size());
        C.a.resize(n);
        uint64_t carry = 0;
        for (size_t i = 0; i < n || carry; i++)
        {
            uint64_t sum = carry;
            if (i < A.a.size())
                sum += A.a[i];
            if (i < B.a.size())
                sum += B.a[i];
            if (i >= C.a.size())
                C.a.push_back(0);
            C.a[i] = sum % BASE;
            carry = sum / BASE;
        }
        return C;
    }

    static BigInt subtract(const BigInt &A, const BigInt &B)
    {
        BigInt C = A;
        uint64_t borrow = 0;
        for (size_t i = 0; i < B.a.size() || borrow; i++)
        {
            uint64_t sub = borrow;
            if (i < B.a.size())
                sub += B.a[i];
            if (C.a[i] < sub)
            {
                C.a[i] += BASE;
                borrow = 1;
            }
            else
            {
                borrow = 0;
            }
            C.a[i] -= sub;
        }
        C.trim();
        return C;
    }

    static BigInt multiply(const BigInt &A, const BigInt &B)
    {
        BigInt C;
        C.a.resize(A.a.size() + B.a.size());
        for (size_t i = 0; i < A.a.size(); i++)
        {
            uint64_t carry = 0;
            for (size_t j = 0; j < B.a.size() || carry; j++)
            {
                uint64_t prod = C.a[i + j] + carry;
                if (j < B.a.size())
                {
                    prod += uint64_t(A.a[i]) * uint64_t(B.a[j]);
                }
                C.a[i + j] = prod % BASE;
                carry = prod / BASE;
            }
        }
        C.trim();
        return C;
    }

    static BigInt divide(const BigInt &A, const BigInt &B, BigInt &R)
    {
        if (compare(A, B) < 0)
        {
            R = A;
            return BigInt(0);
        }
        if (B.isZero())
            throw runtime_error("Division by zero");

        // ✅ Fast path: if B is small (≤ 32-bit), use divInt
        if (B.a.size() == 1)
        {
            uint32_t divisor = B.a[0];
            BigInt Q = A;
            uint32_t rem = divInt(Q, divisor);
            R = BigInt(rem);
            return Q;
        }

        // ✅ General case: long division with FAST binary search (no multiply per iteration)
        BigInt Q;
        R = BigInt(0);
        Q.a.resize(A.a.size());
        for (int i = (int)A.a.size() - 1; i >= 0; i--)
        {
            R = multiply(R, BigInt(BASE));
            if (R.a.empty())
                R.a.push_back(0);
            R.a[0] = A.a[i];
            R.trim();

            // Binary search for x such that B*x <= R < B*(x+1)
            uint32_t left = 0, right = BASE - 1;
            uint32_t x = 0;
            while (left <= right)
            {
                uint32_t mid = left + (right - left) / 2;
                BigInt t = multiply(B, BigInt(mid));
                int cmp = compare(t, R);
                if (cmp <= 0)
                {
                    x = mid;
                    left = mid + 1;
                }
                else
                {
                    right = mid - 1;
                }
            }
            Q.a[i] = x;
            R = subtract(R, multiply(B, BigInt(x)));
        }
        Q.trim();
        return Q;
    } // static BigInt mod(const BigInt &A, const BigInt &B) {
    //     BigInt R;
    //     R.a.clear();
    //     BigInt cur;
    //     cur.a.clear();
    //     for(int i = (int)A.a.size()*32 - 1; i >= 0; i--) {
    //         uint32_t bit = (A.a[i / 32] >> (i % 32)) & 1;
    //         BigInt tmp;
    //         tmp.a.assign(cur.a.size(), 0);
    //         uint64_t carry = 0 ;
    //         for (size_t j = 0; j < cur.a.size(); j++) {
    //             uint64_t val = ((uint64_t)cur.a[j] << 1) | carry;
    //             cur.a[j] = uint32_t(val);
    //             carry = val >> 32;
    //         }
    //         if (carry) cur.a.push_back(uint32_t(carry));
    //         if (bit) {
    //             if (cur.a.empty()) cur.a.push_back(0);
    //             cur.a[0] |= 1;
    //         }
    //         cur.trim();
    //         if (compare(cur, B) >= 0) {
    //             cur = subtract(cur, B);
    //         }
    //     }
    //     cur.trim();
    //     return cur;
    // }

    static BigInt modBig(const BigInt &A, const BigInt &B)
    {
        if (compare(A, B) < 0)
            return A;
        if (B.isZero())
            throw runtime_error("Modulo by zero");

        // Use divide to get remainder efficiently
        BigInt tempA = A; // make a copy since divide may modify first arg
        BigInt remainder;
        divide(tempA, B, remainder);
        return remainder;
    }

    static BigInt power(BigInt base, BigInt exp, const BigInt &mod)
    {
        BigInt result(1);
        base = modBig(base, mod);
        int iter = 0;

        cerr << "Starting power() exp_words=" << exp.a.size() << "\n";

        while (!exp.isZero())
        {
            cout << "DBG: power loop iter=" << iter << " exp(HEX)=" << exp.toHex() << "\n";
            if (++iter % 1000 == 0)
                cerr << "power: iter=" << iter << " exp_words=" << exp.a.size() << "\n";

            if (exp.a[0] & 1)
                result = modBig(multiply(result, base), mod);

            shiftRight1(exp); // ✅ thay cho chia exp / 2
            base = modBig(multiply(base, base), mod);
        }

        cerr << "Finished power(), total iters=" << iter << "\n";
        return result;
    }

    static void shiftRight1(BigInt &A)
    {
        uint32_t carry = 0;
        for (int i = (int)A.a.size() - 1; i >= 0; i--)
        {
            uint64_t cur = ((uint64_t)carry << 32) | A.a[i];
            A.a[i] = (uint32_t)(cur >> 1);
            carry = (uint32_t)(cur & 1);
        }
        A.trim();
        if (A.a.empty())
            A.a.push_back(0);
    }

    static BigInt fromHex(string s)
    {
        BigInt r(0);
        uint64_t val = 0;
        int cnt = 0;
        for (int i = (int)s.size() - 1; i >= 0; i--)
        {
            char c = s[i];
            uint32_t v = (c >= '0' && c <= '9') ? c - '0' : (toupper(c) - 'A' + 10);
            val |= (uint64_t(v) << (4 * cnt++));
            if (cnt == 8)
            {
                r.a.push_back(uint32_t(val));
                val = 0;
                cnt = 0;
            }
        }
        if (cnt)
            r.a.push_back(uint32_t(val));
        r.trim();
        return r;
    }

    string toHex() const
    {
        if (isZero())
            return "0";
        string s;
        for (size_t i = 0; i < a.size(); i++)
        {
            for (int j = 0; j < 8; j++)
            {
                uint32_t nibble = (a[i] >> (4 * j)) & 0xF;
                char c = (nibble < 10) ? ('0' + nibble) : ('A' + nibble - 10);
                s.push_back(c);
            }
        }
        while (s.size() > 1 && s.back() == '0')
        {
            s.pop_back();
        }
        reverse(s.begin(), s.end());
        return s;
    }

    static uint32_t divInt(BigInt &A, uint32_t B)
    {
        uint64_t rem = 0;
        for (int i = (int)A.a.size() - 1; i >= 0; i--)
        {
            uint64_t cur = (uint64_t)rem * BASE + A.a[i];
            A.a[i] = uint32_t(cur / B);
            rem = cur % B;
        }
        A.trim();
        return uint32_t(rem);
    }
};

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

string toLittleEndianHex(const string &bigHex)
{
    string s = bigHex;
    if (s.size() % 2)
        s = "0" + s;
    string res;
    for (size_t i = 0; i < s.size(); i += 2)
        res = s.substr(i, 2) + res;
    return res;
}

int main()
{
    string folder = "project_01_03";
    int total = 9;
    int passCount = 0;

    for (int i = 0; i < total; i++)
    {
        // Tạo tên file input/output
        string num = (i < 10 ? "0" + to_string(i) : to_string(i));
        string inFile = folder + "/test_" + num + ".inp";
        string outFile = folder + "/test_" + num + ".out";
        string myOutFile = folder + "/test_" + num + "_my.out";

        ifstream fin(inFile.c_str());
        if (!fin.is_open())
        {
            cout << "❌ Không mở được file " << inFile << endl;
            continue;
        }

        // Đọc dữ liệu
        string A_hex, B_hex, C_hex;
        fin >> A_hex >> B_hex >> C_hex;
        fin.close();

        cout << A_hex << endl;
        cout << B_hex << endl;
        cout << C_hex << endl;

        // Chuyển từ little → big endian
        BigInt A, B, C;
        try
        {
            cerr << "DBG: before fromHex A\n";
            A = BigInt::fromHex(toBigEndianHex(A_hex));
            cerr << "DBG: after fromHex A\n";
            cerr << "DBG: before fromHex B\n";
            B = BigInt::fromHex(toBigEndianHex(B_hex));
            cerr << "DBG: after fromHex B\n";
            cerr << "DBG: before fromHex C\n";
            C = BigInt::fromHex(toBigEndianHex(C_hex));
            cerr << "DBG: after fromHex C\n";
        }
        catch (const exception &e)
        {
            cerr << "Exception during fromHex(): " << e.what() << "\n";
            return 1;
        }
        // Tính kết quả
        cerr << "DBG A(HEX)=" << A.toHex() << " B(HEX)=" << B.toHex() << " C(HEX)=" << C.toHex() << "\n";
        BigInt result = BigInt::power(A, B, C);
        string resultHex = toLittleEndianHex(result.toHex());

        // Ghi ra file kết quả
        ofstream fout(myOutFile.c_str());
        fout << resultHex << endl;
        fout.close();

        // So sánh với file output chuẩn
        ifstream refFile(outFile.c_str());
        string expected;
        refFile >> expected;
        refFile.close();

        if (resultHex == expected)
        {
            cout << "Test " << i << ": Pass ✅" << endl;
            passCount++;
        }
        else
        {
            cout << "Test " << i << ": Fail ❌" << endl;
            cout << "Expected: " << expected << endl;
            cout << "Got     : " << resultHex << endl;
        }
    }

    cout << "\nTổng kết: " << passCount << "/" << total << " test case pass." << endl;
    return 0;
}
