#include <chrono>
#include "grasshopper.hpp"

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
    Grasshopper::Block buf {{0x11, 0x22, 0x33, 0x44, 0x55, 0x66, 0x77, 0x00, 0xff, 0xee, 0xdd, 0xcc, 0xbb, 0xaa, 0x99, 0x88}};
    Grasshopper::Key key {{0x88, 0x99, 0xaa, 0xbb, 0xcc, 0xdd, 0xee, 0xff, 0x00, 0x11, 0x22, 0x33, 0x44, 0x55, 0x66, 0x77,
                           0xfe, 0xdc, 0xba, 0x98, 0x76, 0x54, 0x32, 0x10, 0x01, 0x23, 0x45, 0x67, 0x89, 0xab, 0xcd, 0xef}};
    Timer timer;
    timer.Start();
    Grasshopper().Encrypt(buf, key);
    timer.Finish();
    printf("Result:\n");
    for (auto c: buf)
        printf("%02x ", c);
    printf("\nTime: %lu ms\n", timer.GetMilliseconds());
    return 0;
}
