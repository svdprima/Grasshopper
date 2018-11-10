#include <algorithm>
#include <functional>
#include <stdio.h>
#include "grasshopper.hpp"
#include "tables.hpp"

void Grasshopper::Encrypt(Block& data, const Key& key)
{
    EncryptBlock(data, key);
}

void Grasshopper::Decrypt(Block& data, const Key& key)
{
    // TODO
}

void Grasshopper::EncryptBlock(Block& data, const Key& key)
{
    const Keys& keys = GenerateKeys(key);
    for (const auto& key: keys)
    {
        for (auto x: key)
            printf("%02x", x);
        printf("\n");
    }
    for (unsigned i = 0; i < num_rounds - 1; ++i)
        ApplyXSL(data, keys[i]);
    ApplyX(data, keys[num_rounds - 1]);
}

void Grasshopper::DecryptBlock(Block& block, const Key& key)
{
    // TODO
}

Grasshopper::Keys Grasshopper::GenerateKeys(const Key& key)
{
    Keys keys;
    std::copy(key.begin(), key.begin() + block_size, keys[0].begin());
    std::copy(key.begin() + block_size, key.end(), keys[1].begin());
    for (unsigned i = 1; i < num_rounds / 2; ++i)
    {
        keys[2 * i    ] = keys[2 * i - 2];
        keys[2 * i + 1] = keys[2 * i - 1];
        for (unsigned j = 0; j < 8; ++j)
        {
            Block C{{}};
            C.back() = (i - 1) * 8 + j + 1;
            ApplyL(C);
            ApplyF(keys[2 * i], keys[2 * i + 1], C);
        }
    }
    return keys;
}

uint8_t Grasshopper::PolynomMul(uint8_t left, uint8_t right)
{
    // p(x) = x^8 + x^7 + x^6 + x + 1 => 0b111000011 => 0xC3 (without MSB)
    uint8_t product = 0;
    while (left && right)
    {
        if (right & 1)
            product ^= left;
        left = (left << 1) ^ (left & 0x80 ? 0xC3 : 0x00);
        right >>= 1;
    }
    return product;
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
    ApplyS(data);
    ApplyL(data);
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

// Linear transform
void Grasshopper::ApplyL(Block& data)
{
    for (unsigned i = 0; i < block_size; ++i)
    {
        uint8_t tmp = 0;
        for (size_t i = 0; i < block_size; ++i)
            tmp ^= PolynomMul(data[i], lin[i]);
        std::copy_backward(data.begin(), data.end() - 1, data.end());
        data[0] = tmp;
    }
}
