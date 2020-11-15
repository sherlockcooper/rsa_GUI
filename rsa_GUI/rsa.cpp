//
// Created by zjl on 11/14/2020.
//

#include "rsa.h"

dtype rsa::charEncodeBase = 128;

rsa::rsa() {
    binaryKeyLen = 6;
    iterNum = 8;
    decimalKeyLen = 2;
    initialized = false;
}

rsa::rsa(int requireLen, int iters) :binaryKeyLen(requireLen), iterNum(iters) {
    decimalKeyLen = int(binaryKeyLen * log10(2)) + 2;
    initialized = false;
}

void rsa::setKeyLen(int keyLen)
{
    binaryKeyLen = keyLen;
    decimalKeyLen = int(binaryKeyLen * log10(2)) + 2;
    initialized = false;
}

void rsa::keygen() {

    //    p = bigInt("89296525811403234983453936780142199152015169371315612071397189690674320351218159801470021036070352895676325128237877");
    //    q = bigInt("24374145121443022467352243614741001614881238071714954329114137605116423951751024209501620142174574221970470948982423");
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = chrono::system_clock::now();
    p = bigInt::createPrime((decimalKeyLen) / 2, iterNum);
    q = bigInt::createPrime((decimalKeyLen) / 2, iterNum);
    end = chrono::system_clock::now();
    // cout << "p: " << p << endl;
    // cout << "q: " << q << endl;
    n = p * q;
    // cout << "n: " << n << endl;
    decimalKeyLen = n.decimalBits();
    phi = (p - 1) * (q - 1);
    // cout << "phi: " << phi << endl;
    vector<dtype> candidates = { 65537, 65539, 65543, 65551, 65557, 65563, 65579, 65581, 65587, 65599 };
    for (int i = 0; i < candidates.size(); i++)
    {
        if (!(phi % candidates[i]).isZero())
        {
            e = bigInt(candidates[i]);
            break;
        }
    }
    // cout << "e: " << e << endl;
    //    e=bigInt(65537);
    bigInt gcd, x;
    bigInt::gcdForInverse(e, phi, gcd, d, x);
    d = d % phi;
    
    timeConsumed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    // cout << "private key d is: " << d << endl;
    maxSegmentLength = int((decimalKeyLen - 1) * log(10) / log(charEncodeBase));
    // cout << "max segment length: " << maxSegmentLength << endl;

    initialized = true;
}

void rsa::encryptStr(const string& content, string& res) {
    vector<bigInt> cypher_num(0), segments(0);
    string cypher_segment, zeros_fill;
    bigInt ele(0);
    for (int i = 0; i < content.size(); i++)
    {
        ele = ele * charEncodeBase + (dtype)content[i];
        if (i % maxSegmentLength == maxSegmentLength - 1 || i == content.size() - 1)
        {
            segments.push_back(ele);
            ele = bigInt::zero;
        }
    }
    for (int i = 0; i < segments.size(); i++)
    {
        ele = bigInt::fastExponent(segments[i], e, n);
        //        cypher_num.push_back(ele);
        cypher_segment = bigInt::toString(ele);
        //        cout << cypher_segment << endl;
        zeros_fill = string(decimalKeyLen - cypher_segment.size(), '0');
        res += zeros_fill + cypher_segment;

        //        if(i!=segments.size()-1)
        //        {
        //            res+="\n";
        //        }
    }
}

void rsa::decryptStr(const string& cypher, string& content) {
    vector<bigInt> cypher_segments;
    bigInt ele, single_val;
    string content_segment;
    int ele_num = cypher.size() / decimalKeyLen;
    for (int i = 0; i < ele_num; i++)
    {
        cypher_segments.push_back(bigInt(cypher.substr(i * decimalKeyLen, decimalKeyLen)));
    }
    for (int i = 0; i < ele_num; i++)
    {
        content_segment = "";
        ele = bigInt::fastExponent(cypher_segments[i], d, n);
        while (!ele.isZero())
        {
            single_val = ele % charEncodeBase;
            content_segment.push_back((char)single_val.data[0]);
            ele = ele / charEncodeBase;
        }
        content_segment = string(content_segment.rbegin(), content_segment.rend());
        content += content_segment;
    }

}





