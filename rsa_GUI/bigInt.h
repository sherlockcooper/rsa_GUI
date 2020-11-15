#pragma once
//
// Created by zjl on 11/8/2020.
//

#ifndef RSA_BIGINT_H
#define RSA_BIGINT_H

#define MAXK 17
#define MAXN 5000

#include <vector>
#include <list>
#include <string>
#include <algorithm>
#include <iostream>
#include <cstring>
#include <random>
#include <iomanip>
#include <chrono>
#include <cmath>
#include <sstream>
#include <thread>
#include <mutex>

using namespace std;

class bigInt;

typedef int64_t dtype;
typedef vector<dtype> container;
typedef pair<bigInt, bigInt> div_pair;
typedef pair<bigInt, int> mod2_pair;
class bigInt {
public:

    static int base;
    static int basebit;
    static long long modulo;
    static long long** alpha;
    static vector<bigInt> little_primes; // contain primes less than 10000 for speed up

    static bigInt zero;
    static bigInt one;
    static bigInt two;
    static bigInt step;

    static vector<bigInt> smallPrimes;

    bigInt();
    bigInt(const bigInt& a);
    bigInt(const bigInt& a, bool isNegative);// only request equal abs value, customized negativity
    bigInt(long long n);
    bigInt(string decS); // decimal string to big int

    container data;
    bool bNegative; // true for negative, false for positive, 0 is false

    bool isZero() const;
    bool isOne() const;
    bool isTwo() const;
    int absComp(bigInt& b) const;
    bool absLess(const bigInt& b) const;
    bool absEqual(bigInt& b);

    int decimalBits();
    friend ostream& operator << (ostream& out, const bigInt& a);
    static string toString(bigInt& a);

    bigInt operator << (const int shift);
    bigInt operator >> (const int shift);
    bigInt operator -() const;

    bigInt operator +(const bigInt& b) const;
    bigInt operator -(const bigInt& b) const;
    bigInt operator *(const bigInt& b);
    div_pair div(const bigInt& b) const; // first ele quotient, second remain
    static void div_2(const bigInt& dividend, bigInt& quo, dtype& rem); // speed up when divisor is 2
    bigInt operator /(const bigInt& b);
    bigInt operator %(const bigInt& b) const;

    void shrink();
    void initializeAlpha();

    static void extendedGcd(long long a, long long b, long long& gcd, long long& x, long long& y);
    static void gcdForInverse(bigInt& a, bigInt& b, bigInt& gcd, bigInt& x, bigInt& y);
    int base2hundred(const container& baseTen, long long* baseHundred);
    void hundred2base(long long* baseHundred, container& baseTen, int high);
    void convolution(long long* a, long long* b, int bitLength, long long* res);
    void NTT_N2R(long long* a, int nBit);
    void NTT_R2N(long long* a, int nBit);

    static void initializePrimes();
    static bigInt fastExponent(bigInt& a, bigInt& e, bigInt& n);
    static bigInt randomOdd(int bitNumInTen, random_device& generator);
    static bigInt createPrime(int bitNum, int iterCount);
    //    static void findPrimeThread(bigInt& start, mutex& primeMutex, bool &finished, bigInt& res, int iterCount);

private:

};


#endif //RSA_BIGINT_H
