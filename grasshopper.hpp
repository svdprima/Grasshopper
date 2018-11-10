#ifndef GRASSHOPPER_H
#define GRASSHOPPER_H

#include <array>
#include <cstdint>

class Grasshopper
{
public:
    static constexpr std::size_t block_size = 16;
    static constexpr unsigned num_rounds = 10;

    using Block = std::array<uint8_t, block_size>;
    using Key = std::array<uint8_t, block_size * 2>;

    void Encrypt(Block& data, const Key& key);
    void Decrypt(Block& data, const Key& key);
private:
    using Keys = std::array<Block, num_rounds>;

    void EncryptBlock(Block& data, const Key& key);
    void DecryptBlock(Block& data, const Key& key);

    Keys GenerateKeys(const Key& key);
    uint8_t PolynomMul(uint8_t left, uint8_t right);

    void ApplyF(Block& data0, Block& data1, const Block& key);
    void ApplyXSL(Block& data, const Block& key);
    void ApplyX(Block& data, const Block& key);
    void ApplyS(Block& data);
    void ApplyL(Block& data);
};

#endif
