#ifndef GRASSHOPPER_H
#define GRASSHOPPER_H

#include <array>
#include <vector>
#include <utility>
#include <cstdint>

class Grasshopper
{
public:
    static constexpr std::size_t block_size = 16;
    static constexpr unsigned num_rounds = 10;

    using Block = std::array<uint8_t, block_size>;
    using Key = std::array<uint8_t, block_size * 2>;

    Grasshopper();

    void Encrypt(std::vector<Block>& data, const Key& key);
    void Decrypt(std::vector<Block>& data, const Key& key);
private:
    using Keys = std::array<Block, num_rounds>;
    using KeyPair = std::pair<Block, Block>;
    uint8_t mul_table[256][256];

    void EncryptBlock(Block& data, const KeyPair& key);
    void DecryptBlock(Block& data, const KeyPair& key);

    Keys GenerateKeys(const KeyPair& key);
    void GenerateMulTable();

    void ApplyF(Block& data0, Block& data1, const Block& key);
    void ApplyXSL(Block& data, const Block& key);
    void ApplyX(Block& data, const Block& key);
    void ApplyS(Block& data);
    void ApplyL(Block& data);
};

void DumpBlock(const Grasshopper::Block& block);

#endif
