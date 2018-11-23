#include <chrono>
#include <iostream>
#include <string>
#include "encryptor.hpp"

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

int main(int argc, char** argv)
{
    Mode mode = Mode::ECB;
    std::string file("input.txt");
    for (int i = 1; i < argc; ++i)
    {
        std::string tmp(argv[i]);
        if (tmp == "-ECB")
            mode = Mode::ECB;
        else if (tmp == "-CBC")
            mode = Mode::CBC;
        else
            file = tmp;
    }

    Encryptor E(mode);

    // encryption
    E.ReadText(file);

    Timer timer;
    // encrypt
    timer.Start();

    E.Encrypt();

    timer.Finish();
    std::cout << "Encrypt time: " << timer.GetMilliseconds() << " ms\n";

    E.SaveText("encrypted.txt");

    // decryption
    E.ReadText("encrypted.txt", true);

    timer.Start();

    E.Decrypt();

    timer.Finish();
    std::cout << "Decrypt time: " << timer.GetMilliseconds() << " ms\n";

    E.SaveText("decrypted.txt", true);

    return 0;
}
