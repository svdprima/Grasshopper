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
    uint8_t mul_table[256][16];
    uint8_t l_table[16][16] = {};
    uint8_t ls_table[16][256][16] = {};
    Block coef_table[num_rounds / 2 - 1][8];

    void EncryptBlock(Block& data, const KeyPair& key);
    void DecryptBlock(Block& data, const KeyPair& key);

    Keys GenerateKeys(const KeyPair& key);
    void GenerateMulTable();
    void GenerateCoefTable();
    void GenerateLTable();
    uint8_t PolyMul(uint8_t left, uint8_t right);

    void ApplyF(Block& data0, Block& data1, const Block& key);
    void ApplyXSL(Block& data, const Block& key);
    void ApplyX(Block& data, const Block& key);
    void ApplyS(Block& data);
    void ApplyL(Block& data);

    void ApplyInvXLS(Block& data, const Block& key);
    void ApplyInvS(Block& data);
    void ApplyInvL(Block& data);

    void ApplyLS(Block& data);
};

void DumpBlock(const Grasshopper::Block& block);

#endif
