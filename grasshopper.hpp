#ifndef GRASSHOPPER_H
#define GRASSHOPPER_H

#include <array>
#include <vector>
#include <utility>
#include <cstdlib>
#include <cstdint>

template<typename T>
class AlignedAllocator
{
public:
    typedef T value_type;
    typedef std::size_t size_type;
    typedef std::ptrdiff_t difference_type;

    T* allocate(std::size_t n)
    {
        if (!n)
            return nullptr;
        if (n > max_size)
            throw std::length_error("AlignedAllocator::allocate() - integer overflow");
        void *p = std::aligned_alloc(sizeof(T), n * sizeof(T));
        if (!p)
            throw std::bad_alloc();
        return static_cast<T*>(p);
    }

    void deallocate(T* p, std::size_t n __attribute__((unused)))
    {
        std::free(static_cast<void*>(p));
    }

    bool operator ==(const AlignedAllocator& rhs __attribute__((unused))) const
    {
        return true;
    }

    bool operator !=(const AlignedAllocator& rhs) const
    {
        return !(*this == rhs);
    }

private:
    static constexpr std::size_t max_size =
        (static_cast<std::size_t>(0) - static_cast<std::size_t>(1)) / sizeof(T);
};

class Grasshopper
{
public:
    static constexpr std::size_t block_size = 16;
    static constexpr unsigned num_rounds = 10;

    using Block = struct alignas(block_size) : std::array<uint8_t, block_size> {};
    using Key = std::array<uint8_t, block_size * 2>;
    using Data = std::vector<Block, AlignedAllocator<Block>>;

    Grasshopper();

    void Encrypt(Data& data, const Key& key);
    void Decrypt(Data& data, const Key& key);
private:
    using Matrix = std::array<Block, block_size>;
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
