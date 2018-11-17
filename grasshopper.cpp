#include <algorithm>
#include <functional>
#include <cstdio>
#include <string.h>
#include <fstream>
#include <iostream>
#include "grasshopper.hpp"
#include "tables.hpp"

Grasshopper::Grasshopper()
{
    //GenerateMulTable();
    GenerateCoefTable();
    GenerateLTable();
}

void Grasshopper::Encrypt(std::vector<Block>& data, const Key& key)
{
    KeyPair current_key;
    std::copy(key.begin(), key.begin() + block_size, current_key.first.begin());
    std::copy(key.begin() + block_size, key.end(), current_key.second.begin());
    for (auto& block: data)
    {
        EncryptBlock(block, current_key);
        ApplyX(current_key.first , block);
        ApplyX(current_key.second, block);
    }
}

void Grasshopper::Decrypt(std::vector<Block>& data, const Key& key)
{
    KeyPair current_key;
    std::copy(key.begin(), key.begin() + block_size, current_key.first.begin());
    std::copy(key.begin() + block_size, key.end(), current_key.second.begin());

    KeyPair tmp_key;
    tmp_key = current_key;

    for (auto& block: data)
    {
        ApplyX(tmp_key.first , block);
        ApplyX(tmp_key.second, block);

        DecryptBlock(block, current_key);

        current_key = tmp_key;
    }
}

void Grasshopper::EncryptBlock(Block& data, const KeyPair& key)
{
#ifdef VERBOSE
    printf("before encryption:\n");
    DumpBlock(data);
#endif

    const Keys& keys = GenerateKeys(key);
#ifdef VERBOSE
    printf("keys:\n");
    for (const auto& key: keys)
        DumpBlock(key);
#endif

    for (unsigned i = 0; i < num_rounds - 1; ++i)
    {
        ApplyXSL(data, keys[i]);
    }
    ApplyX(data, keys[num_rounds - 1]);

#ifdef VERBOSE
    printf("after encryption:\n");
    DumpBlock(data);
#endif
}

void Grasshopper::DecryptBlock(Block& data, const KeyPair& key)
{
#ifdef VERBOSE
    printf("before decryption:\n");
    DumpBlock(data);
#endif

    const Keys& keys = GenerateKeys(key);
#ifdef VERBOSE
    printf("keys:\n");
    for (const auto& key: keys)
        DumpBlock(key);
#endif

    for (unsigned i = 0; i < num_rounds - 1; ++i)
        ApplyInvXLS(data, keys[num_rounds - 1 - i]);
    ApplyX(data, keys[0]);

#ifdef VERBOSE
    printf("after decryption:\n");
    DumpBlock(data);
#endif
}

Grasshopper::Keys Grasshopper::GenerateKeys(const KeyPair& key)
{
    Keys keys{{key.first, key.second}};
    for (unsigned i = 1; i < num_rounds / 2; ++i)
    {
        keys[2 * i    ] = keys[2 * i - 2];
        keys[2 * i + 1] = keys[2 * i - 1];
        for (unsigned j = 0; j < 8; ++j)
            ApplyF(keys[2 * i], keys[2 * i + 1], coef_table[i - 1][j]);
    }
    return keys;
}

uint8_t Grasshopper::PolyMul (uint8_t left, uint8_t right)
{
    // p(x) = x^8 + x^7 + x^6 + x + 1 => 0b111000011 => 0xC3 (without MSB)
    uint8_t res = 0;
    while (left && right)
    {
        if (right & 1)
            res ^= left;
        left = (left << 1) ^ (left & 0x80 ? 0xC3 : 0x00);
        right >>= 1;
    }
    return res;
}


void Grasshopper::GenerateMulTable()
{
    for (unsigned i = 0; i < 256; ++i)
        for (unsigned j = 0; j < 16; ++j)
        {
            mul_table[i][j] = PolyMul(i, lin[j]);
        }
}

void Grasshopper::GenerateCoefTable()
{
    for (unsigned i = 0; i < num_rounds / 2 - 1; ++i)
        for (unsigned j = 0; j < 8; ++j)
        {
            Block& C = coef_table[i][j];
            C.fill(0);
            C.back() = i * 8 + j + 1;
            ApplyL(C);
        }
}

void Grasshopper::ApplyF(Block& data1, Block& data0, const Block& key)
{
    Block tmp = data1;
    ApplyXSL(tmp, key);
    ApplyX(tmp, data0);
    std::tie(data1, data0) = std::make_pair(tmp, data1);
}

void Grasshopper::ApplyXSL(Block& data, const Block& key)
{
    ApplyX(data, key);
    //ApplyS(data);
    //ApplyL(data);
    ApplyLS (data);
}

void Grasshopper::ApplyInvXLS(Block& data, const Block& key)
{
    ApplyX(data, key);
    ApplyInvL(data);
    ApplyInvS(data);
}

// Xor transform
void Grasshopper::ApplyX(Block& data, const Block& key)
{
    std::transform(data.begin(), data.end(), key.begin(), data.begin(), std::bit_xor<uint8_t>());
}

// Nonlinear transform
void Grasshopper::ApplyS(Block& data)
{
    std::transform(data.begin(), data.end(), data.begin(), [](uint8_t idx) {return S[idx];});
}

void Grasshopper::ApplyInvS(Block& data)
{
    std::transform(data.begin(), data.end(), data.begin(), [](uint8_t idx) {return invS[idx];});
}

// Linear transform
void Grasshopper::ApplyL(Block& data)
{
    for (unsigned i = 0; i < block_size; ++i)
    {
        uint8_t tmp = 0;
        for (size_t j = 0; j < block_size; ++j)
            tmp ^= mul_table[data[j]][j];
        std::copy_backward(data.begin(), data.end() - 1, data.end());
        data[0] = tmp;
    }
}

void Grasshopper::ApplyInvL(Block& data)
{
    for (unsigned i = 0; i < block_size; ++i)
    {
        uint8_t elem = data[0];
        std::copy(data.begin() + 1, data.end(), data.begin());
        data[block_size - 1] = elem;

        uint8_t tmp = 0;
        for (size_t j = 0; j < block_size; ++j)
            tmp ^= mul_table[data[j]][j];
        data[15] = tmp;
    }
}

void Grasshopper::GenerateLTable()
{
    //initializing the martix
    for (size_t i = 0; i < block_size; i++)
    {
        for (size_t j = 0; j < block_size; j++)
        {
            if (i - 1 == j)
                l_table[i][j] = 1;
            else if (i == 0)
                l_table[i][j] = lin[j];
        }
    }
    
    //powering the matrix
    //NB! in this loop one gets a matrix to degree of 2^p
    uint8_t tmp[16][16] = {};
    for (size_t p = 0; p < 4; p++)
    {
        for (size_t i = 0; i < block_size; i++)
            for (size_t j = 0; j < block_size; j++)
                for (size_t k = 0; k < block_size; k++) 
                    tmp[i][j] ^= PolyMul(l_table[i][k], l_table[k][j]);
        memcpy (l_table, tmp, 256 * sizeof(uint8_t)); 
        memset (tmp, 0, 256 * sizeof (uint8_t));
    }

    for (size_t i = 0; i < block_size; i++)
        for (size_t j = 0; j < 256; j++)
            for (size_t k = 0; k < block_size; k++)
                ls_table[i][j][k] ^= PolyMul(S[j], l_table[k][i]); 


     /*
     //Incorrect brackets order!
     FILE* mtrx;
     mtrx = fopen ("matrix.hpp", "w");
     fprintf (mtrx, "//LS autogenerated precomputed table created from grasshopper.cpp\n");
     fprintf (mtrx, "uint8_t LStable[16][256][16] =\n{");
     for (unsigned int i = 0; i < block_size; i++)
     {
         fprintf(mtrx, "\n\t//Row #%d\n\t{\n", i);
         for (unsigned int j = 0; j < 256; j++)
         {
             fprintf(mtrx, "\t\t{");
             for (unsigned int k = 0; k < block_size; k++)
             {
                 fprintf(mtrx, "0x%02x", ls_table[i][j][k]);
                 if (k < block_size - 1)
                     fprintf(mtrx, ", ");
             }
             fprintf(mtrx, "}, //byte 0x%02x\n", j);
         }
         fprintf(mtrx, "\t},");
     }
     fprintf(mtrx, "};\n");
     fclose(mtrx);
     */
}

//Implementing linear transform via matrices, merging with non-linear
void Grasshopper::ApplyLS(Block& data)
{
    uint8_t tmp[16] = {};
    for (size_t i = 0; i < block_size; i++)
        for (size_t j = 0; j < block_size; j++)
            tmp[j] ^= ls_table[i][data[i]][j];
    for (size_t i = 0; i < block_size; i++)
    {
        data[i] = tmp[i]; 
    }
}

void DumpBlock(const Grasshopper::Block& block)
{
    for (auto x: block)
        printf("%02x", x);
    printf("\n");
}
