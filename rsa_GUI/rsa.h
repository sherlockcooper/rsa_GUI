#pragma once
#pragma once
//
// Created by zjl on 11/14/2020.
//

#ifndef RSA_RSA_H
#define RSA_RSA_H

#include "bigInt.h"

class rsa
{
public:
    static dtype charEncodeBase;

    int binaryKeyLen; // length of key in base 2, key value <= 2^BinaryKeyLen-1
    int decimalKeyLen; // length of key in base 10
    int iterNum;
    int maxSegmentLength; // max length of string segment

    bool initialized;

    bigInt p;
    bigInt q;
    bigInt n; //n=p*q
    bigInt phi; //Euler value
    bigInt e; // public key
    bigInt d; // private key

    chrono::milliseconds timeConsumed;

    rsa();
    rsa(int requireLen=768, int iters = 10);

    void setKeyLen(int keyLen);

    void keygen();
    void encryptStr(const string& content, string& cypher);
    void decryptStr(const string& cypher, string& content);


};

#endif //RSA_RSA_H
