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
    using Matrix = std::array<Block, block_size>;
    using Key = std::array<uint8_t, block_size * 2>;

    Grasshopper();

    void Encrypt(std::vector<Block>& data, const Key& key);
    void Decrypt(std::vector<Block>& data, const Key& key);
private:
    using Keys = std::array<Block, num_rounds>;
    using KeyPair = std::pair<Block, Block>;
    Block coef_table[num_rounds / 2 - 1][8];
    Block enc_ls_table[block_size][256];
    Block dec_ls_table[block_size][256];

    void EncryptBlock(Block& data, const KeyPair& key);
    void DecryptBlock(Block& data, const KeyPair& key);

    void ApplyXSL(Block& data, const Block& key);
    void ApplyInvXLS(Block& data, const Block& key);

    Keys GenerateKeys(const KeyPair& key);
    void ApplyF(Block& data0, Block& data1, const Block& key);

    void GenerateCoefTable();
    void GenerateEncTable();
    void GenerateDecTable();

    uint8_t PolyMul(uint8_t left, uint8_t right);
    Matrix SqrMatrix(const Matrix& mat);
    void ApplyX(Block& data, const Block& key);
    void ApplyL(Block& data);
};

void DumpBlock(const Grasshopper::Block& block);

#endif
