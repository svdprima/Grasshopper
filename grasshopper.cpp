#include <algorithm>
#include <functional>
#include <cstdio>
#include "grasshopper.hpp"
#include "tables.hpp"

Grasshopper::Grasshopper()
{
    GenerateMulTable();
    GenerateCoefTable();
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
    // TODO
}

void Grasshopper::EncryptBlock(Block& data, const KeyPair& key)
{
#ifdef VERBOSE
    printf("before:\n");
    DumpBlock(data);
#endif

    const Keys& keys = GenerateKeys(key);
#ifdef VERBOSE
    printf("keys:\n");
    for (const auto& key: keys)
        DumpBlock(key);
#endif

    for (unsigned i = 0; i < num_rounds - 1; ++i)
        ApplyXSL(data, keys[i]);
    ApplyX(data, keys[num_rounds - 1]);

#ifdef VERBOSE
    printf("after:\n");
    DumpBlock(data);
#endif
}

void Grasshopper::DecryptBlock(Block& data, const KeyPair& key)
{
    // TODO
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

void Grasshopper::GenerateMulTable()
{
    for (unsigned i = 0; i < 256; ++i)
        for (unsigned j = 0; j < 16; ++j)
        {
            // p(x) = x^8 + x^7 + x^6 + x + 1 => 0b111000011 => 0xC3 (without MSB)
            uint8_t left = i, right = lin[j],res = 0;
            while (left && right)
            {
                if (right & 1)
                    res ^= left;
                left = (left << 1) ^ (left & 0x80 ? 0xC3 : 0x00);
                right >>= 1;
            }
            mul_table[i][j] = res;
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
        for (size_t j = 0; j < block_size; ++j)
            tmp ^= mul_table[data[j]][j];
        std::copy_backward(data.begin(), data.end() - 1, data.end());
        data[0] = tmp;
    }
}

void DumpBlock(const Grasshopper::Block& block)
{
    for (auto x: block)
        printf("%02x", x);
    printf("\n");
}
