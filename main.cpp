#include <chrono>
#include <cstdio>
#include <cstdlib>
#include "grasshopper.hpp"

Grasshopper::Key key {{0x88, 0x99, 0xaa, 0xbb, 0xcc, 0xdd, 0xee, 0xff, 0x00, 0x11, 0x22, 0x33, 0x44, 0x55, 0x66, 0x77,
                       0xfe, 0xdc, 0xba, 0x98, 0x76, 0x54, 0x32, 0x10, 0x01, 0x23, 0x45, 0x67, 0x89, 0xab, 0xcd, 0xef}};

class Timer
{
public:
    void Start()
    {
        t1 = std::chrono::high_resolution_clock::now();
    }
    void Finish()
    {
        t2 = std::chrono::high_resolution_clock::now();
    }
    uint64_t GetMilliseconds()
    {
        return std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    }

private:
    std::chrono::high_resolution_clock::time_point t1, t2;
};

int main()
{
    FILE *f;
    // read data
    if (!(f = fopen("input.txt", "r")))
    {
        printf("Can not open input.txt\n");
        exit(EXIT_FAILURE);
    }
    fseek(f, 0, SEEK_END);
    size_t file_size = ftell(f);
    fseek(f, 0, SEEK_SET);
    size_t num_blocks = (file_size / Grasshopper::block_size) +
                        (file_size % Grasshopper::block_size ? 1 : 0);
    printf("Number of blocks: %lu\n", num_blocks);
    std::vector<Grasshopper::Block> buf(num_blocks);
    fread(buf.data(), Grasshopper::block_size, num_blocks, f);
    fclose(f);
    Grasshopper G;
    Timer timer;
    // encrypt
    timer.Start();
    G.Encrypt(buf, key);
    timer.Finish();
    printf("Encrypt time: %lu ms\n", timer.GetMilliseconds());
    if (!(f = fopen("encrypted.txt", "w")))
    {
        printf("Can not open output.txt\n");
        exit(EXIT_FAILURE);
    }
    fwrite(buf.data(), Grasshopper::block_size, num_blocks, f);
    fclose(f);
    // decrypt
    timer.Start();
    G.Decrypt(buf, key);
    timer.Finish();
    printf("Decrypt time: %lu ms\n", timer.GetMilliseconds());
    if (!(f = fopen("decrypted.txt", "w")))
    {
        printf("Can not open output.txt\n");
        exit(EXIT_FAILURE);
    }
    fwrite(buf.data(), Grasshopper::block_size, num_blocks, f);
    fclose(f);

    return 0;
}
