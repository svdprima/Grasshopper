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

unsigned ArgParser (int argc, char** argvec)
{
    std::string tmp;
    for (int i = 1; i < argc; i++)
    {
        tmp = std::string(argvec[i]);
        if (tmp == "-ECB\0")
            return 0;
        if (tmp == "-CBC\0")
            return 1;
    } 
    return 0;
}

int main(int argc, char** argv)
{
    unsigned mode = ArgParser (argc, argv);
    Encryptor E (mode);

    // encryption
    E.ReadText("input.txt");

    Timer timer;
    // encrypt
    timer.Start();

    E.Encrypt(mode);

    timer.Finish();
    std::cout << "Encrypt time: " << timer.GetMilliseconds() << " ms\n";

    E.SaveText("encrypted.txt");

    // decryption
    E.ReadText("encrypted.txt", true);

    timer.Start();

    E.Decrypt(mode);

    timer.Finish();
    std::cout << "Decrypt time: " << timer.GetMilliseconds() << " ms\n";

    E.SaveText("decrypted.txt", true);

    return 0;
}
