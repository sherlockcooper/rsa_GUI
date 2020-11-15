//
// Created by zjl on 11/8/2020.
//

#include "bigInt.h"

//int bigInt::base = 10;
int bigInt::base = 100000000; // 10^8
int bigInt::basebit = 8;
long long bigInt::modulo = 1325137921;
long long** bigInt::alpha = NULL;
bigInt bigInt::zero(0);
bigInt bigInt::one(1);
bigInt bigInt::two(2);
bigInt bigInt::step(4);
vector<bigInt> bigInt::little_primes(0);


bigInt::bigInt() :bNegative(false) {
    //    data.push_back(0);

}

bigInt::bigInt(const bigInt& a) : data(a.data), bNegative(a.bNegative) {
    shrink();
}

bigInt::bigInt(const bigInt& a, bool isNegative) : data(a.data), bNegative(isNegative) {
    shrink();
}

bigInt::bigInt(string decS) {
    int start = 0;
    bNegative = false;
    if (decS[start] == '-')
    {
        bNegative = true;
        start++;
    }
    int rem_bit = (decS.size() - start) % basebit;
    dtype bit_num = 0;

    for (int i = 0; i < rem_bit; i++)
    {
        bit_num = bit_num * 10 + (decS[start] - '0');
        start++;
    }
    data.push_back(bit_num);
    int tot_bit = (decS.size() - start) / basebit;
    for (int i = 0; i < tot_bit; i++)
    {
        bit_num = 0;
        for (int j = 0; j < basebit; j++)
        {
            bit_num = bit_num * 10 + (decS[start] - '0');
            start++;
        }
        //        data.push_back(decS[start]-'0');
        data.push_back(bit_num);
    }
    shrink();
}

bigInt::bigInt(long long n) {
    long tmp = n;
    vector<int> decomp;
    bNegative = false;
    if (tmp < 0)
    {
        bNegative = true;
        int bit = -(tmp % base);
        decomp.push_back(bit);
        tmp = -(tmp / base);
    }
    dtype bit_num = 0, bit = 0;
    do {
        bit = tmp % base;
        decomp.push_back(bit);
        tmp /= base;
    } while (tmp);
    data = container(decomp.rbegin(), decomp.rend());
    shrink();
}

void bigInt::shrink() {
    int cnt = 0;
    for (auto it = data.begin(); it != data.end(); it++)
    {
        if ((*it) == 0)
            cnt++;
        else
            break;
    }

    if (cnt == data.size())
    {
        cnt--;
        // zero
        if (!data[cnt])
            bNegative = false;
    }
    //    for(int i=0; i<cnt; i++)
    //    {
    //        data.pop_back();
    //    }
    data = container(data.begin() + cnt, data.end());
}

bigInt bigInt::operator+(const bigInt& b) const {
    if (b.isZero())
        return (*this);
    if (isZero())
        return b;

    if (bNegative == b.bNegative)
    {
        bigInt res;
        res.bNegative = bNegative;
        const bigInt* longer = data.size() >= b.data.size() ? this : &b, * shorter = data.size() < b.data.size() ? this : &b;
        int min_len = shorter->data.size(), max_len = longer->data.size(), carry = 0;
        dtype l_bit, s_bit, bit_sum;
        res.data = container(max_len + 1, 0);
        for (int i = 0; i < min_len; i++)
        {
            l_bit = longer->data[max_len - i - 1];
            s_bit = shorter->data[min_len - i - 1];
            bit_sum = l_bit + s_bit + carry;
            res.data[max_len - i] = bit_sum % base;
            carry = bit_sum / base;
        }
        for (int i = min_len; i < max_len; i++)
        {
            l_bit = longer->data[max_len - i - 1];
            bit_sum = l_bit + carry;
            res.data[max_len - i] = bit_sum % base;
            carry = bit_sum / base;
        }
        res.data[0] = carry;
        res.shrink();
        return res;
    }
    else
    {
        bigInt opposite = -b;
        return (*this) - opposite;
    }
}

bigInt bigInt::operator-(const bigInt& b) const {
    if (b.isZero())
        return (*this);
    if (isZero())
        return -b;

    if (bNegative != b.bNegative)
    {
        bigInt opposite = -b;
        return (*this) + opposite;
    }
    else
    {
        const bigInt* larger, * smaller;
        bigInt res;
        if (absLess(b))
        {
            res.bNegative = (!bNegative);
            larger = &b;
            smaller = this;
        }
        else
        {
            res.bNegative = bNegative;
            larger = this;
            smaller = &b;
        }
        int min_len = smaller->data.size(), max_len = larger->data.size(), carry = 0, bit_minus;
        res.data = container(max_len, 0);
        for (int i = 0; i < min_len; i++)
        {
            bit_minus = larger->data[max_len - 1 - i] - smaller->data[min_len - 1 - i] - carry;
            if (bit_minus < 0)
            {
                bit_minus += base;
                carry = 1;
            }
            else
            {
                carry = 0;
            }
            res.data[max_len - 1 - i] = bit_minus;
        }
        for (int i = min_len; i < max_len; i++)
        {
            bit_minus = larger->data[max_len - 1 - i] - carry;
            if (bit_minus < 0)
            {
                bit_minus += base;
                carry = 1;
            }
            else
            {
                carry = 0;
            }
            res.data[max_len - 1 - i] = bit_minus;
        }
        res.shrink();
        return res;
    }
}

bool bigInt::absLess(const bigInt& b) const {
    if (data.size() < b.data.size())
        return true;
    if (data.size() > b.data.size())
        return false;
    int len = data.size();
    for (int i = 0; i < len; i++)
    {
        if (data[i] > b.data[i])
            return false;
        if (data[i] < b.data[i])
            return true;
    }
    return false;
}

bigInt bigInt::operator*(const bigInt& b) {
    bigInt res;
    initializeAlpha();
    if (bNegative == b.bNegative)
        res.bNegative = false;
    else
        res.bNegative = true;
    // hundred-based array, n[0] is LBS
    long long n[2][MAXN] = { 0 };
    int l1 = base2hundred(data, n[0]);
    int l2 = base2hundred(b.data, n[1]);
    int k = 0;
    l1 = l1 < l2 ? l2 : l1;
    while ((1 << k) < l1)
    {
        ++k;
    }
    long long res_data[MAXN] = { 0 };
    convolution(n[0], n[1], k + 1, res_data);
    int high = (1 << (k + 1)) - 1;

    for (int i = 0; i < high; ++i) {
        res_data[i + 1] += res_data[i] / 100;
        res_data[i] %= 100;
    }
    hundred2base(res_data, res.data, high);
    res.shrink();
    return res;
}

div_pair bigInt::div(const bigInt& b) const {
    if (b.isOne())
    {
        return div_pair(bigInt(*this), bigInt(0));
    }
    bigInt res, rem;
    if (absLess(b) || isZero())
    {
        res.data.push_back(0);
        // need to consider sign of remain
        rem = bigInt(*this);
    }
    else
    {
        rem = bigInt(*this);
        if (rem.bNegative)
            rem.bNegative = false;
        bigInt to_add;
        // start stands for pos of highest bit in tmp
        int s_len = b.data.size(), start = 0, pos;
        while (rem.data.size() - start >= s_len)
        {
            if (rem.data[start] > b.data[0])
            {
                dtype bit_att = rem.data[start] / (b.data[0] + 1); // attempt a result
                for (int i = s_len - 1; i >= 0; i--)
                {
                    pos = start + i;
                    rem.data[pos] -= bit_att * b.data[i];
                    if (rem.data[pos] < 0)
                    {
                        dtype t_val = -rem.data[pos] / base + 1;
                        rem.data[pos] += t_val * base;
                        rem.data[pos - 1] -= t_val;
                    }
                }
                to_add.data = container(rem.data.size() - start - s_len + 1, 0);
                to_add.data[0] = bit_att;
                res = res + to_add;
            }
            else
            {
                if (rem.data.size() - start == s_len)
                {
                    int cmp = 0;
                    // compare remain part of tmp with b
                    for (int i = 0; i < s_len; i++)
                    {
                        if (rem.data[start + i] < b.data[i])
                        {
                            cmp = -1;
                            break;
                        }
                        else if (rem.data[start + i] > b.data[i])
                        {
                            cmp = 1;
                            break;
                        }
                    }
                    if (cmp >= 0)
                    {
                        for (int i = s_len - 1; i >= 0; i--)
                        {
                            pos = start + i;
                            rem.data[pos] -= b.data[i];
                            if (rem.data[pos] < 0)
                            {
                                rem.data[pos] += base;
                                rem.data[pos - 1] -= 1;
                            }
                        }
                        res = res + bigInt(1);
                        while (start < rem.data.size() && (!rem.data[start]))
                            start++;
                    }
                    break;
                }
                rem.data[start + 1] += rem.data[start] * base;
                rem.data[start] = 0;
                start++;
            }
        }
    }
    if (bNegative != b.bNegative)
    {
        rem.bNegative = !b.bNegative;
        rem = b + rem;
    }
    res.bNegative = (bNegative != b.bNegative);
    res.shrink();
    rem.shrink();
    return div_pair(res, rem);
}

bigInt bigInt::operator/(const bigInt& b) {
    return div(b).first;
}

bigInt bigInt::operator%(const bigInt& b) const {
    return div(b).second;
}

string bigInt::toString(bigInt& a) {
    ostringstream ss;
    ss << a;
    return ss.str();
}

ostream& operator<<(ostream& out, const bigInt& a) {
    if (a.bNegative)
        out << "-";
    out << a.data[0];
    for (int i = 1; i < a.data.size(); i++)
    {
        out << setw(bigInt::basebit) << setfill('0') << a.data[i];
    }
    return out;
}

void bigInt::initializeAlpha() {
    if (alpha != NULL)
        return;
    alpha = new long long* [2];
    for (int i = 0; i < 2; i++)
        alpha[i] = new long long[MAXK + 1];
    alpha[0][MAXK] = 1101238606;
    long long gcd, y;
    extendedGcd(alpha[0][MAXK], modulo, gcd, alpha[1][MAXK], y);
    for (int i = MAXK - 1; i >= 0; i--)
    {
        alpha[0][i] = alpha[0][i + 1] * alpha[0][i + 1] % modulo;
        alpha[1][i] = alpha[1][i + 1] * alpha[1][i + 1] % modulo;
    }
}

void bigInt::extendedGcd(long long a, long long b, long long& gcd, long long& x, long long& y) {
    if (!b)
    {
        x = 1;
        y = 0;
        gcd = a;
    }
    else {
        extendedGcd(b, a % b, gcd, y, x);
        y -= a / b * x;
    }
}


bigInt bigInt::operator-() const {
    return bigInt((*this), !bNegative);
}

int bigInt::base2hundred(const container& baseContainer, long long* baseHundred) {
    int len = baseContainer.size(), tot = 0, cnt, ret_val;
    dtype top_ele = baseContainer[0];
    while (top_ele)
    {
        top_ele /= 100;
        tot++;
    }
    cnt = tot;
    tot += (len - 1) * (basebit / 2);
    ret_val = tot;
    for (int i = 0; i < baseContainer.size(); i++)
    {
        top_ele = baseContainer[i];
        for (int j = cnt; j > 0; j--)
        {
            baseHundred[tot - j] = top_ele % 100;
            top_ele /= 100;
        }
        tot -= cnt;
        cnt = 4;
    }
    return ret_val;
}

void bigInt::hundred2base(long long int* baseHundred, container& baseContainer, int high) {
    int limit = basebit / 2, len = (high + 1) / limit + 1, bit_hundred;
    baseContainer = container(len, 0);
    dtype ele;
    for (int i = 0; i < len; i++)
    {
        bit_hundred = (4 * i + 3) > (high + 1) ? high + 1 : 4 * i + 3;
        ele = 0;
        for (; bit_hundred >= 4 * i; bit_hundred--)
        {
            ele = ele * 100 + baseHundred[bit_hundred];
        }
        baseContainer[len - 1 - i] = ele;
    }
}

void bigInt::convolution(long long int* a, long long int* b, int bitLength, long long int* res) {
    int& k = bitLength;

    int nLen = 1 << k;

    NTT_N2R(a, k);
    NTT_N2R(b, k);

    for (int i = 0; i < nLen; ++i)
        res[i] = a[i] * b[i] % modulo;

    NTT_R2N(res, k);
}

void bigInt::NTT_N2R(long long int* a, int nBit) {
    int nLen = 1 << nBit;

    // butterfly operation
    for (int layer = nBit; layer > 0; --layer)
    {
        int group = 1 << layer, brother = group >> 1;
        long long kernel = alpha[0][layer];

        for (int k = 0; k < nLen; k += group) {

            long long w = 1;

            for (int j = 0; j < brother; ++j) {

                int cur = k + j, next = cur + brother;
                long long u = a[cur], v = a[next];

                a[cur] = u + v;
                a[cur] = a[cur] < modulo ? a[cur] : a[cur] - modulo;
                // no overflow with current modulo
                a[next] = (u + modulo - v) * w % modulo;
                w = w * kernel % modulo;
            }
        }
    }
}

void bigInt::NTT_R2N(long long int* a, int nBit) {
    int nLen = 1 << nBit;

    // butterfly operation
    for (int layer = 1; layer <= nBit; ++layer)
    {
        int group = 1 << layer, brother = group >> 1;
        long long kernel = alpha[1][layer], & half = alpha[1][1];

        for (int k = 0; k < nLen; k += group) {

            long long w = 1;

            for (int j = 0; j < brother; ++j) {

                int cur = k + j, next = cur + brother;
                long long u = a[cur], t = w * a[next] % modulo;

                a[cur] = u + t;
                a[cur] = a[cur] < modulo ? a[cur] : a[cur] - modulo;

                a[next] = u - t;
                a[next] = a[next] < 0 ? a[next] + modulo : a[next];
                w = w * kernel % modulo;
            }
        }
    }

    long long gcd = 0, invN = 0, invM = 0;

    extendedGcd(nLen, modulo, gcd, invN, invM);

    invN = (invN % modulo + modulo) % modulo;

    for (int i = 0; i < nLen; ++i)
        a[i] = a[i] * invN % modulo;
}

bigInt bigInt::operator<<(const int shift) {
    bigInt res(*this);
    for (int i = 0; i < shift; i++)
    {
        res.data.push_back(0);
    }
    return res;
}

bigInt bigInt::operator>>(const int shift) {
    bigInt res;
    if (shift >= data.size())
        res.data.push_back(0);
    else
    {
        for (int i = 0; i < data.size() - shift; i++)
        {
            res.data.push_back(data[i]);
        }
    }
    return res;
}

bool bigInt::absEqual(bigInt& b) {
    shrink();
    b.shrink();
    if (data.size() != b.data.size())
        return false;
    for (int i = 0; i < data.size(); i++)
    {
        if (data[i] != b.data[i])
            return false;
    }
    return true;
}

bool bigInt::isZero() const {
    return data.size() == 1 && data[0] == 0;
}


bool bigInt::isOne() const {
    return (!bNegative) && data.size() == 1 && data[0] == 1;
}

bool bigInt::isTwo() const {
    return (!bNegative) && data.size() == 1 && data[0] == 2;
}

int bigInt::absComp(bigInt& b) const {
    if (data.size() < b.data.size())
        return -1;
    if (data.size() > b.data.size())
        return 1;
    for (int i = 0; i < data.size(); i++)
    {
        if (data[i] < b.data[i])
            return -1;
        if (data[i] > b.data[i])
            return 1;
    }
    return 0;
}

bigInt bigInt::randomOdd(int bitNumInTen, random_device& generator) {
    bigInt res;
    int len = bitNumInTen / basebit + 1, rem_bit = bitNumInTen % basebit;
    res.data = container(len, 0);
    dtype bit_range = 0, bit_ele;
    for (int i = 0; i < rem_bit; i++)
    {
        bit_range = bit_range * 10 + 10;
    }
    uniform_int_distribution<dtype> msb_distribution(bit_range / 10 + 1, bit_range - 1), oth_distribution(0, base - 1);
    bit_ele = msb_distribution(generator);
    res.data[0] = bit_ele;
    for (int i = 1; i < len; i++)
    {
        bit_ele = oth_distribution(generator);
        res.data[i] = bit_ele;
    }
    if (!(res.data.back() & 1))
    {
        res.data[len - 1] = (res.data[len - 1] + 1) % base;
    }
    return res;
}

void findPrimeThread(bigInt start, mutex* primeMutex, bool* finished, bigInt* res, int iterCount) {
    random_device generator;
    bigInt possible = start, psb_minus1(possible - bigInt::one), q(psb_minus1), quo, checknum;
    dtype rem, tmp_num;
    uniform_int_distribution<dtype> checknum_distribution(3, INT64_MAX);
    int k = 0; //possible = q*2^k+1
    bool passed, early_check;
    while (true)
    {
        possible = possible + bigInt::step;
        psb_minus1 = possible - bigInt::one;
        q = psb_minus1;
        //        primeMutex->lock();
        if (*finished)
            break;
        //        primeMutex->unlock();
                // quick early check
        early_check = true;
        for (int i = 0; i < bigInt::little_primes.size(); i++)
        {
            if ((possible % bigInt::little_primes[i]).isZero())
            {
                early_check = false;
                break;;
            }
        }
        if (!early_check)
            continue;

        bigInt::div_2(q, quo, rem);
        while (!rem)
        {
            k++;
            q = quo;
            bigInt::div_2(q, quo, rem);
        }

        for (int i = 0; i < iterCount; i++)
        {
            passed = false; // pass this iter
            checknum = bigInt(checknum_distribution(generator));
            checknum = bigInt::fastExponent(checknum, q, possible);
            if (checknum.isOne() || checknum.absEqual(psb_minus1))
            {
                passed = true;
                continue;
            }
            for (int j = 1; j < k; j++)
            {
                checknum = bigInt::fastExponent(checknum, bigInt::two, possible);
                if (checknum.isOne())
                {
                    passed = false;
                    break;
                }
                if (checknum.absEqual(psb_minus1))
                {
                    passed = true;
                    break;
                }
            }
            if (!passed)
            {
                break;
            }
        }
        if (passed) // found the prime
        {
            break;
        }

    }
    primeMutex->lock();
    if (!(*finished))
    {
        *res = possible;
        *finished = true;
    }
    primeMutex->unlock();

}

void simple(int i, bigInt r, bool* f, mutex* m, bigInt a)
{
    cout << "thread" << endl;
}

bigInt bigInt::createPrime(int bitNum, int iterCount) {

    random_device num_generator;
    mutex prime_mutex;
    initializePrimes();
    bool finished = false;
    int thread_num = 6;
    vector<thread> workers(thread_num);
    bigInt some_place_to_start = randomOdd(bitNum, num_generator), res(0);
    for (int i = 0; i < thread_num; i++)
    {
        workers.emplace_back(findPrimeThread, some_place_to_start + (i * 256), &prime_mutex, &finished, &res, iterCount);
        //        thread t(findPrimeThread, some_place_to_start+(i*256), &prime_mutex, &finished, &res, iterCount);
        //        thread t(simple);
        //        workers.push_back(t);
        //        cout << i << endl;
        // workers.emplace_back(simple, i, res, &finished, &prime_mutex, some_place_to_start + (i * 256));
        //        workers.push_back(t);
    }
    for (auto& w : workers)
    {
        if (w.joinable())
        {
            w.join();
        }
    }

    return res;
    // original part
//    initializePrimes();
//
//    random_device generator;
//
//    bigInt possible = randomOdd(bitNum, generator), psb_minus1(possible-one),q(psb_minus1), quo, checknum;
//    dtype rem, tmp_num;
//    uniform_int_distribution<dtype> checknum_distribution(3, INT64_MAX);
//    int k=0; //possible = q*2^k+1
//    bool passed, early_check;
//    while (true)
//    {
//        possible = possible+step;
//        psb_minus1=possible-one;
//        q=psb_minus1;
//
//        // quick early check
//        early_check = true;
//        for(int i=0; i<little_primes.size(); i++)
//        {
//            if((possible%little_primes[i]).isZero())
//            {
//                early_check = false;
//                break;;
//            }
//        }
//        if(!early_check)
//            continue;
//
//        div_2(q, quo, rem);
//        while (!rem)
//        {
//            k++;
//            q = quo;
//            div_2(q, quo, rem);
//        }
//
//        for(int i=0; i<iterCount; i++)
//        {
//            passed = false; // pass this iter
//            checknum = bigInt(checknum_distribution(generator));
//            checknum = fastExponent(checknum, q, possible);
//            if(checknum.isOne() || checknum.absEqual(psb_minus1))
//            {
//                passed = true;
//                continue;
//            }
//            for(int j=1; j<k; j++)
//            {
//                checknum = fastExponent(checknum, two, possible);
//                if (checknum.isOne())
//                {
//                    passed = false;
//                    break;
//                }
//                if(checknum.absEqual(psb_minus1))
//                {
//                    passed = true;
//                    break;
//                }
//            }
//            if(!passed)
//            {
//                break;
//            }
//        }
//        if(passed)
//            break;
//
//    }
//
//    return possible;
}


bigInt bigInt::fastExponent(bigInt& a, bigInt& e, bigInt& n) {
    bigInt res(a % n), quo, tmp(e);

    if (res.isZero() || res.isOne() || e.isOne())
        return res;
    dtype rem;
    container trail(0);
    bigInt* p_quo = &quo, * p_diviend = &tmp, * p_tmp;
    //    div_2(tmp, quo, rem);
    //    while (!quo.isZero())
    //    {
    //        trail.push_back(rem);
    //        tmp = quo;
    //        div_2(tmp, quo, rem);
    //    }
    div_2(*p_diviend, *p_quo, rem);
    while (!p_quo->isZero())
    {
        trail.push_back(rem);
        p_tmp = p_diviend;
        p_diviend = p_quo;
        p_quo = p_tmp;
        div_2(*p_diviend, *p_quo, rem);
    }
    bigInt init(res);
    while (!trail.empty())
    {
        rem = trail.back();
        trail.pop_back();
        if (rem)
        {
            res = (((res * res) % n) * init) % n;
            //            res = ((res*res)*init)%n;
        }
        else
        {
            res = (res * res) % n;
        }
    }
    return res;
}


void bigInt::div_2(const bigInt& dividend, bigInt& quo, dtype& rem) {
    dtype carry = 0, bit_num;
    //    rem = dividend.data.back()%2;
    rem = dividend.data.back() & 1;
    quo.data = container(dividend.data.size(), 0);
    for (int i = 0; i < dividend.data.size(); i++)
    {
        bit_num = dividend.data[i];
        if (carry)
        {
            bit_num += base;
        }
        quo.data[i] = bit_num >> 1;
        carry = bit_num & 1;
    }
    quo.shrink();
}

void bigInt::initializePrimes() {
    if (little_primes.size() != 0)
        return;

    vector<int> primes = {
            3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 113,
            193, 241, 257, 337, 353, 401, 433, 449, 577, 593, 641, 673, 769, 881, 929, 977, 1009, 1153, 1201, 1217, 1249,
            1297, 1361, 1409, 1489, 1553, 1601, 1697, 1777, 1873, 1889, 2017, 2081, 2113, 2129, 2161, 2273, 2417, 2593,
            2609, 2657, 2689, 2753, 2801, 2833, 2897, 3041, 3089, 3121, 3137, 3169, 3217, 3313, 3329, 3361, 3457, 3617,
            3697, 3761, 3793, 3889, 4001, 4049, 4129, 4177, 4241, 4273, 4289, 4337, 4481, 4513, 4561, 4657, 4673, 4721,
            4801, 4817, 4993, 5009, 5153, 5233, 5281, 5297, 5393, 5441, 5521, 5569, 5857, 5953, 6113, 6257, 6337, 6353,
            6449, 6481, 6529, 6577, 6673, 6689, 6737, 6833, 6961, 6977, 7057, 7121, 7297, 7393, 7457, 7489, 7537, 7649,
            7681, 7793, 7841, 7873, 7937, 8017, 8081, 8161, 8209, 8273, 8353, 8369, 8513, 8609, 8641, 8689, 8737, 8753,
            8849, 8929, 9041, 9137, 9281, 9377, 9473, 9521, 9601, 9649, 9697, 9857
    };

    for (int i = 0; i < primes.size(); i++)
    {
        little_primes.push_back(bigInt(primes[i]));
    }

}

void bigInt::gcdForInverse(bigInt& a, bigInt& b, bigInt& gcd, bigInt& x, bigInt& y) {
    if (b.isZero())
    {
        x = one;
        y = zero;
        gcd = a;
    }
    else
    {
        bigInt tmp = a % b;
        gcdForInverse(b, tmp, gcd, y, x);
        y = y - (a / b * x);
    }
}

int bigInt::decimalBits() {
    int res = basebit * (data.size() - 1);
    dtype top_ele = data.front();
    while (top_ele != 0)
    {
        res++;
        top_ele /= 10;
    }
    return res;
}










