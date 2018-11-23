#include <algorithm>
#include <functional>
#include <cstdio>
#include <xmmintrin.h>
#include <immintrin.h>
#include "grasshopper.hpp"
#include "tables.hpp"

const __m128i& CastBlock(const Grasshopper::Block& block)
{
    return *reinterpret_cast<const __m128i*>(&block);
}

__m128i& CastBlock(Grasshopper::Block& block)
{
    return *reinterpret_cast<__m128i*>(&block);
}

Grasshopper::Grasshopper()
{
    GenerateCoefTable();
    GenerateEncTable();
    GenerateDecTable();
}

void Grasshopper::Encrypt(Data& data, const Key& key, Mode mode)
{
    KeyPair current_key;
    std::copy(key.begin(), key.begin() + block_size, current_key.first.begin());
    std::copy(key.begin() + block_size, key.end(), current_key.second.begin());

    if (mode == Mode::ECB)
    {
        const Keys& keys = GenerateKeys(current_key);
        for (auto& block: data)
            EncryptBlock(block, keys);
    }
    else
    {
        for (auto& block: data)
        {
            const Keys& keys = GenerateKeys(current_key);
            EncryptBlock(block, keys);
            ApplyX(current_key.first , block);
            ApplyX(current_key.second, block);
        }
    }
}

void Grasshopper::Decrypt(Data& data, const Key& key, Mode mode)
{
    KeyPair current_key;
    std::copy(key.begin(), key.begin() + block_size, current_key.first.begin());
    std::copy(key.begin() + block_size, key.end(), current_key.second.begin());

    KeyPair tmp_key(current_key);
    if (mode == Mode::ECB)
    {
        const Keys& keys = GenerateKeys(current_key);
        for (auto& block: data)
            DecryptBlock(block, keys);
    }
    else
    {
        for (auto& block: data)
        {
            ApplyX(tmp_key.first , block);
            ApplyX(tmp_key.second, block);
            const Keys& keys = GenerateKeys(current_key);
            DecryptBlock(block, keys);
            current_key = tmp_key;
        }
    }
}

void Grasshopper::EncryptBlock(Block& data, const Keys& keys)
{

#ifdef VERBOSE
    printf("keys:\n");
    for (const auto& key: keys)
        DumpBlock(key);
    printf("before encryption:\n");
    DumpBlock(data);
#endif

    for (unsigned i = 0; i < num_rounds - 1; ++i)
        ApplyXSL(data, keys[i]);
    ApplyX(data, keys[num_rounds - 1]);

#ifdef VERBOSE
    printf("after encryption:\n");
    DumpBlock(data);
#endif
}

void Grasshopper::DecryptBlock(Block& data, const Keys& keys)
{

#ifdef VERBOSE
    printf("keys:\n");
    for (const auto& key: keys)
        DumpBlock(key);
    printf("before decryption:\n");
    DumpBlock(data);
#endif

    for (unsigned i = 0; i < num_rounds - 1; ++i)
        ApplyInvXLS(data, keys[num_rounds - 1 - i]);
    ApplyX(data, keys[0]);

#ifdef VERBOSE
    printf("after decryption:\n");
    DumpBlock(data);
#endif
}

void Grasshopper::ApplyXSL(Block& data, const Block& key)
{
    ApplyX(data, key);
    __m256i vec1 = _mm256_setzero_si256();
#pragma clang loop unroll(full)
    for (size_t i = 0; i < block_size; i += 2)
    {
        __m256i vec2 = _mm256_castps128_ps256(CastBlock(enc_ls_table[i][data[i]]));
        vec2 = _mm256_insertf128_ps(vec2, CastBlock(enc_ls_table[i + 1][data[i + 1]]), 1);
        vec1 = _mm256_xor_si256(vec1, vec2);
    }
    CastBlock(data) = _mm_xor_si128(_mm256_extracti128_si256(vec1, 0),
                                    _mm256_extracti128_si256(vec1, 1));
}

void Grasshopper::ApplyInvXLS(Block& data, const Block& key)
{
    ApplyX(data, key);
    Block tmp{};
    __m256i vec1 = _mm256_setzero_si256();
#pragma clang loop unroll(full)
    for (size_t i = 0; i < block_size; i += 2)
    {
        __m256i vec2 = _mm256_castps128_ps256(CastBlock(dec_ls_table[i][data[i]]));
        vec2 = _mm256_insertf128_ps(vec2, CastBlock(dec_ls_table[i + 1][data[i + 1]]), 1);
        vec1 = _mm256_xor_si256(vec1, vec2);
    }
    CastBlock(tmp) = _mm_xor_si128(_mm256_extracti128_si256(vec1, 0),
                                   _mm256_extracti128_si256(vec1, 1));
    std::transform(tmp.begin(), tmp.end(), data.begin(),
                   [](uint8_t idx) { return Table::invS[idx]; });
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

void Grasshopper::ApplyF(Block& data1, Block& data0, const Block& key)
{
    Block tmp = data1;
    ApplyXSL(tmp, key);
    ApplyX(tmp, data0);
    std::tie(data1, data0) = std::make_pair(tmp, data1);
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

void Grasshopper::GenerateEncTable()
{
    //initializing the martix
    Matrix l_matrix;
    for (size_t i = 0; i < block_size; ++i)
        for (size_t j = 0; j < block_size; ++j)
            if (i == 0)
                l_matrix[i][j] = Table::lin[j];
            else if (i == j + 1)
                l_matrix[i][j] = 1;
            else
                l_matrix[i][j] = 0;
    // power l_matrix to degree of 2^4
    for (unsigned i = 0; i < 4; ++i)
        l_matrix = SqrMatrix(l_matrix);

    for (size_t i = 0; i < block_size; ++i)
        for (size_t j = 0; j < 256; ++j)
            for (size_t k = 0; k < block_size; ++k)
                enc_ls_table[i][j][k] = PolyMul(Table::S[j], l_matrix[k][i]);
}

void Grasshopper::GenerateDecTable()
{
    //initializing the martix
    Matrix l_matrix;
    for (size_t i = 0; i < block_size; ++i)
        for (size_t j = 0; j < block_size; ++j)
            if (i == block_size - 1)
                l_matrix[i][j] = Table::lin[(j + 15) % block_size];
            else if (i + 1 == j)
                l_matrix[i][j] = 1;
            else
                l_matrix[i][j] = 0;
    // power l_matrix to degree of 2^4
    for (unsigned i = 0; i < 4; ++i)
        l_matrix = SqrMatrix(l_matrix);

    for (size_t i = 0; i < block_size; ++i)
        for (size_t j = 0; j < 256; ++j)
            for (size_t k = 0; k < block_size; ++k)
                dec_ls_table[i][j][k] = PolyMul(j, l_matrix[k][i]);
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

Grasshopper::Matrix Grasshopper::SqrMatrix(const Matrix& mat)
{
    Matrix res{};
    for (size_t i = 0; i < block_size; ++i)
        for (size_t j = 0; j < block_size; ++j)
            for (size_t k = 0; k < block_size; ++k)
                res[i][j] ^= PolyMul(mat[i][k], mat[k][j]);
    return res;
}

void Grasshopper::ApplyX(Block& data, const Block& key)
{
    __m128i& data_128 = CastBlock(data);
    const __m128i& key_128 = CastBlock(key);
    data_128 = _mm_xor_si128(data_128, key_128);
}

void Grasshopper::ApplyL(Block& data)
{
    for (unsigned i = 0; i < block_size; ++i)
    {
        uint8_t tmp = 0;
        for (size_t j = 0; j < block_size; ++j)
            tmp ^= PolyMul(data[j], Table::lin[j]);
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
