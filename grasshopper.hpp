#ifndef GRASSHOPPER_H
#define GRASSHOPPER_H

#include <array>
#include <vector>
#include <utility>
#include <cstdlib>
#include <cstdint>
#include <immintrin.h>

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

enum class Mode
{
    ECB,
    CBC,
    CFB,
    OFB
};

class Grasshopper
{
public:
    static constexpr std::size_t block_size = 16;
    static constexpr unsigned num_rounds = 10;

    struct alignas(block_size) Block : std::array<uint8_t, block_size>
    {
        Block() = default;

        Block(__m128i val)
        {
            *reinterpret_cast<__m128i*>(this) = val;
        }

        Block(const std::array<uint8_t, block_size>& arr)
            : std::array<uint8_t, block_size>(arr)
        {
        }

        operator __m128i() const
        {
            return *reinterpret_cast<const __m128i*>(this);
        }

        void Dump()
        {
            for (auto x: *this)
                printf("%02x", x);
            printf("\n");
        }
    };

    struct alignas(block_size * 2) DoubleBlock : std::array<uint8_t, block_size * 2>
    {
        DoubleBlock() = default;

        DoubleBlock(__m256i val)
        {
            *reinterpret_cast<__m256i*>(this) = val;
        }

        DoubleBlock(const std::array<uint8_t, block_size * 2>& arr)
            : std::array<uint8_t, block_size * 2>(arr)
        {
        }

        operator __m256i() const
        {
            return *reinterpret_cast<const __m256i*>(this);
        }

        void Dump()
        {
            for (auto x: *this)
                printf("%02x", x);
            printf("\n");
        }
    };

    using Key = std::array<uint8_t, block_size * 2>;
    using Data = std::vector<Block, AlignedAllocator<Block>>;

    Grasshopper();

    void Encrypt(Data& data, const Key& key, Mode mode);
    void Decrypt(Data& data, const Key& key, Mode mode);
private:
    using Matrix = std::array<Block, block_size>;
    using Keys = std::array<Block, num_rounds>;
    using KeyPair = std::pair<Block, Block>;
    using LookupTable = Block[block_size][256];
    Block coef_table[num_rounds / 2 - 1][8];
    uint8_t mul_table[256][256];
    LookupTable enc_ls_table;
    LookupTable dec_ls_table;

    template<typename BlockType>
    void EncryptBlock(BlockType& data, const Keys& keys);
    template<typename BlockType>
    void DecryptBlock(BlockType& data, const Keys& keys);

    template<typename BlockType>
    void ApplyXSL(BlockType& data, const Block& key);
    template<typename BlockType>
    void ApplyInvXLS(BlockType& data, const Block& key);
    void ApplyLS(Block& data, const LookupTable& lut);
    void ApplyLS(DoubleBlock& data, const LookupTable& lut);

    Keys GenerateKeys(const KeyPair& key);
    void ApplyF(Block& data0, Block& data1, const Block& key);

    void GenerateMulTable();
    void GenerateCoefTable();
    void GenerateEncTable();
    void GenerateDecTable();

    uint8_t PolyMul(uint8_t left, uint8_t right);
    Matrix SqrMatrix(const Matrix& mat);
    void ApplyX(Block& data, const Block& key);
    void ApplyX(DoubleBlock& data, const Block& key);
    void ApplyL(Block& data);
    void ApplyInvL(Block& data);
};

#endif
