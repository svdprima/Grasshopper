#ifndef GRASSHOPPER_H
#define GRASSHOPPER_H

#include <array>
#include <vector>
#include <cstdint>

class Grasshopper
{
public:
    static constexpr std::size_t block_size = 16;
    static constexpr unsigned num_rounds = 10;

    using Block = std::array<uint8_t, block_size>;
    using Key = std::array<uint8_t, block_size * 2>;

    Grasshopper(const Key& key);

    void Encrypt(std::vector<Block>& data);
    void Decrypt(std::vector<Block>& data);
private:
    std::array<Block, num_rounds> keys;
    uint8_t mul_table[256][256];

    void EncryptBlock(Block& data);
    void DecryptBlock(Block& data);

    void GenerateKeys(const Key& key);
    void GenerateMulTable();

    void ApplyF(Block& data0, Block& data1, const Block& key);
    void ApplyXSL(Block& data, const Block& key);
    void ApplyX(Block& data, const Block& key);
    void ApplyS(Block& data);
    void ApplyL(Block& data);
};

void DumpBlock(const Grasshopper::Block& block);

#endif
